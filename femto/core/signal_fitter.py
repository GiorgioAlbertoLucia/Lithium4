import numpy as np
from ROOT import TFile, TCanvas, TH1F, TTree, \
                 RooRealVar, RooCrystalBall, RooFit, RooHistPdf, RooDataHist, RooWorkspace, RooDataSet, RooKeysPdf
from torchic import AxisSpec

import sys
sys.path.append('/home/galucia/Lithium4/femto')
from core.fitter import Fitter
from core.utils import write_params_to_text

class SignalFitter(Fitter):

    def __init__(self, name, xvar_spec: AxisSpec, outfile:TFile = None, workspace:RooWorkspace = None):
        super().__init__(name, xvar_spec, outfile, workspace)
        self._signal_pars = {}
        self._signal_pdf = None
        self._title = ''

        self._outdir = self._outfile.mkdir('signal') if self._outfile else None

    @property
    def signal_pars(self):
        return self._signal_pars
    
    @property
    def signal_pdf(self):
        return self._signal_pdf
    
    @property
    def title(self):
        return self._title
    
    @title.setter
    def title(self, title:str):
        self._title = title
        if self._signal_pdf is not None:
            self._signal_pdf.SetTitle(self._title)

    
    def _init_signal_from_mc(self, h_signal:TH1F, name:str='signal_pdf'):
        
        h_signal.SetBinContent(h_signal.FindBin(0.), 0.)    # Cure the unphysical first bin behaving like an underflow bin

        xvar = self._roo_workspace.obj(self._xvar_name)
        self._signal_datahist = RooDataHist('signal_dh', 'signal_dh', [xvar], Import=h_signal)
        self._signal_pdf = RooHistPdf(name, name, [xvar], self._signal_datahist)

        frame = xvar.frame()
        self._signal_datahist.plotOn(frame)
        self._signal_pdf.plotOn(frame)

        if self._outdir:
            self._outdir.cd()
            canvas = TCanvas('signal_pdf')
            frame.Draw()
            canvas.Write()

    def _init_signal_from_kde(self, h_signal:TH1F, name:str='signal_pdf', xmin:float=0., xmax:float=0.42, rho:float=2.):
    
        xvar = self._roo_workspace.obj(self._xvar_name)
        old_range = (xvar.getMin(), xvar.getMax())
        xvar.setRange(xmin, xmax)

        x_data, weights = [], []
        for ibin in range(1, h_signal.GetNbinsX()+1):
            x_val = h_signal.GetBinCenter(ibin)
            if not (xmin <= x_val <= xmax):
                continue
            y_val = h_signal.GetBinContent(ibin)
            if y_val <= 0:
                continue
            
            x_data.append(x_val)
            weights.append(y_val)
        
        tree = TTree('tree', 'tree')
        x = np.zeros(1, dtype=np.float64)
        w = np.zeros(1, dtype=np.float64)
        tree.Branch('kstar', x, 'kstar/D')
        tree.Branch('weight', w, 'weight/D')
        
        for x_val, w_val in zip(x_data, weights):
            x[0] = x_val
            w[0] = w_val
            tree.Fill()
        
        weight_var = RooRealVar('weight', 'weight', 0, 1e6)
        self._signal_dataset = RooDataSet(h_signal.GetName()+'_roodata', h_signal.GetName()+'_roodata', 
                                          [xvar, weight_var],
                                          RooFit.Import(tree),
                                          RooFit.WeightVar(weight_var)) #, '', 'weight')
        getattr(self._roo_workspace, 'import')(self._signal_dataset)
        
        self._signal_pdf = RooKeysPdf(name, name, xvar, self._signal_dataset, RooKeysPdf.NoMirror, rho)
        
        frame = xvar.frame(len(x_data))
        self._signal_dataset.plotOn(frame, MarkerStyle=20, MarkerSize=0.8, LineColor=1)
        self._signal_pdf.plotOn(frame, LineColor=2, LineWidth=2)
        canvas = TCanvas(f'cKeysPdf_{h_signal.GetName()}', f'cKeysPdf_{h_signal.GetName()}', 800, 600)
        frame.Draw()

        if self._outdir:
            self._outdir.cd()
            h_signal.Write(f'{h_signal.GetName()}_original')
            self._signal_pdf.Write(f'{h_signal.GetName()}_keyspdf')
            canvas.Write()

        xvar.setRange(*old_range)

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

    def init_signal(self, mode:str, *args, **kwargs):

        if mode == 'from_mc':           self._init_signal_from_mc(*args, **kwargs)
        elif mode == 'from_kde':        self._init_signal_from_kde(*args, **kwargs)
        elif mode == 'crystal_ball':    self._init_signal_crystal_ball(**kwargs)
        else:                           raise ValueError('Only supported modes are "from_mc" and "crystal_ball"')
        

    def fit_signal(self, range_limits:tuple, h_signal:TH1F):

        xvar = self._roo_workspace.obj(self._xvar_name)
        datahist = RooDataHist('signal_dh', 'signal_dh', [xvar], Import=h_signal)

        self._signal_pdf.fitTo(datahist, RooFit.Save(), 
                              RooFit.Range(range_limits[0], range_limits[1]), SumW2Error=True, Extended=True,
                              PrintLevel=-1, Verbose=False)
        frame = xvar.frame()
        datahist.plotOn(frame)
        self._signal_pdf.plotOn(frame)

        text = write_params_to_text(self._signal_pars.values(), coordinates=(0.7, 0.3, 0.85, 0.5))
        text.AddText(f'#chi^{{2}} / ndf = {frame.chiSquare():.2f}')
        frame.addObject(text)

        if self._outdir:
            self._outdir.cd()
            canvas = TCanvas('fit_signal')
            frame.Draw()
            canvas.Write()

    def save_to_workspace(self):

        getattr(self._roo_workspace, 'import')(self._signal_pdf)
        for param in self._signal_pars.values():
            getattr(self._roo_workspace, 'import')(param)

    