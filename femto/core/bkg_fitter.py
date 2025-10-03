from ROOT import TFile, \
                 RooRealVar, RooCrystalBall, RooFit
from torchic import AxisSpec
from fitter import Fitter
from utils import write_params_to_text

class BkgFitter(Fitter):

    def __init__(self, name, xvar_spec: AxisSpec, outfile:TFile):
        super().self.__init__(name, xvar_spec, outfile)
        self._bkg_pars = None
        self._singal_pdf = None

        self._outdir = self._outfile.mkdir('bkg')

    @property
    def bkg_pars(self):
        return self._bkg_pars
    
    @property
    def bkg_pdf(self):
        return self._bkg_pdf

    def init_bkg(self, name:str='bkg_pdf'):
        
        self._bkg_pars = {
            'mean': RooRealVar('sig_mean', '#mu', 0.081, 0.06, 0.1, 'GeV/c'),
            'sigma': RooRealVar('sig_sigma', '#sigma', 0.01, 0.001, 0.1, 'GeV/c'),
            'aL': RooRealVar('sig_aL', 'a_{L}', 1.3, 0.1, 10.),
            'nL': RooRealVar('sig_nL', 'n_{L}', 2.7, 0.1, 10.),
            'aR': RooRealVar('sig_aR', 'a_{R}', 1.1, 0.1, 10.),
            'nR': RooRealVar('sig_nR', 'n_{R}', 5.7, 0.1, 10.),
        }
        self._bkg_pdf = RooCrystalBall(name, name, *self._bkg_pars.values())

    def fit_bkg(self, range_limits:tuple):

        xvar = self._roo_workspace.obj(self._xvar_name)
        datahist = self._roo_workspace.get(self._roo_data_hist_name)
        self._bkg_pdf.fitTo(datahist, RooFit.Save(), 
                              RooFit.Range(range_limits[0], range_limits[1]), SumW2Error=True, Extended=True)
        frame = xvar.frame()
        datahist.plotOn(frame)
        self._bkg_pdf.plotOn(frame)

        text = write_params_to_text(self.sig_params.values(), coordinates=(0.7, 0.3, 0.9, 0.5))
        frame.addObject(text)

        self._outdir.cd()
        frame.Write(f'fit_bkg')

    def save_to_workspace(self):

        getattr(self._roo_workspace, 'import')(self._bkg_pdf)
        for param in self._bkg_pars.values():
            getattr(self._roo_workspace, 'import')(param)

    