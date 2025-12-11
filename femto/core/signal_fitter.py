from ROOT import TFile, TCanvas, TH1F, \
                 RooRealVar, RooCrystalBall, RooFit, RooHistPdf, RooDataHist, RooWorkspace
from torchic import AxisSpec

import sys
sys.path.append('/home/galucia/Lithium4/femto')
from core.fitter import Fitter
from core.utils import write_params_to_text

class SignalFitter(Fitter):

    def __init__(self, name, xvar_spec: AxisSpec, outfile:TFile, workspace:RooWorkspace = None):
        super().__init__(name, xvar_spec, outfile, workspace)
        self._signal_pars = {}
        self._signal_pdf = None

        self._outdir = self._outfile.mkdir('signal')

    @property
    def signal_pars(self):
        return self._signal_pars
    
    @property
    def signal_pdf(self):
        return self._signal_pdf
    
    def _init_signal_from_mc(self, h_signal:TH1F, name:str='signal_pdf'):
        
        h_signal.SetBinContent(h_signal.FindBin(0.), 0.)    # Cure the unphysical first bin behaving like an underflow bin

        xvar = self._roo_workspace.obj(self._xvar_name)
        self._signal_datahist = RooDataHist('signal_dh', 'signal_dh', [xvar], Import=h_signal)
        self._signal_pdf = RooHistPdf(name, name, [xvar], self._signal_datahist)

        frame = xvar.frame()
        self._signal_datahist.plotOn(frame)
        self._signal_pdf.plotOn(frame)

        self._outdir.cd()
        canvas = TCanvas('signal_pdf')
        frame.Draw()
        canvas.Write()

    def _init_signal_crystal_ball(self, name:str='signal_pdf'):
        
        xvar = self._roo_workspace.obj(self._xvar_name)
        self._signal_pars = {
            'mean': RooRealVar('sig_mean', '#mu', 0.081, 0.06, 0.1, 'GeV/c'),
            'sigma': RooRealVar('sig_sigma', '#sigma', 0.01, 0.001, 0.1, 'GeV/c'),
            'aL': RooRealVar('sig_aL', 'a_{L}', 1.3, 0.1, 10.),
            'nL': RooRealVar('sig_nL', 'n_{L}', 2.7, 0.1, 10.),
            'aR': RooRealVar('sig_aR', 'a_{R}', 1.1, 0.1, 10.),
            'nR': RooRealVar('sig_nR', 'n_{R}', 5.7, 0.1, 10.),
        }
        self._signal_pdf = RooCrystalBall(name, name, xvar, *self._signal_pars.values())

    def init_signal(self, mode:str, *args, name:str='signal_pdf'):

        if mode == 'from_mc':           self._init_signal_from_mc(*args, name=name)
        elif mode == 'crystal_ball':    self._init_signal_crystal_ball(name=name)
        else:                           raise ValueError('Only supported modes are "from_mc" and "crystal_ball"')
        

    def fit_signal(self, range_limits:tuple, h_signal:TH1F):

        xvar = self._roo_workspace.obj(self._xvar_name)
        datahist = RooDataHist('signal_dh', 'signal_dh', [xvar], Import=h_signal)

        self._signal_pdf.fitTo(datahist, RooFit.Save(), 
                              RooFit.Range(range_limits[0], range_limits[1]), SumW2Error=True, Extended=True)
        frame = xvar.frame()
        datahist.plotOn(frame)
        self._signal_pdf.plotOn(frame)

        text = write_params_to_text(self._signal_pars.values(), coordinates=(0.7, 0.3, 0.85, 0.5))
        text.AddText(f'#chi^{{2}} / ndf = {frame.chiSquare():.2f}')
        frame.addObject(text)

        self._outdir.cd()
        canvas = TCanvas('fit_signal')
        frame.Draw()
        canvas.Write()

    def save_to_workspace(self):

        getattr(self._roo_workspace, 'import')(self._signal_pdf)
        for param in self._signal_pars.values():
            getattr(self._roo_workspace, 'import')(param)

    