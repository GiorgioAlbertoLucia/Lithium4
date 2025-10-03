
from ROOT import TH1F, TFile, \
                 RooDataHist, RooWorkspace, RooRealVar
from abc import ABC
from torchic import AxisSpec

class Fitter(ABC):

    def __init__(self, name, xvar_spec: AxisSpec, outfile:TFile): 
        
        self._hist_data = None
        self._outfile = outfile
        self._roo_workspace = RooWorkspace(name if name else 'roows')

        self._xvar_name = xvar_spec.name
        xvar = RooRealVar(xvar_spec.name, xvar_spec.title, xvar_spec.xmin, xvar_spec.xmax)
        getattr(self._roo_workspace, 'import')(xvar)
        self._roo_data_hist_name = None

    @property
    def roo_workspace(self):
        return self._roo_workspace
    
    @property
    def xvar_name(self):
        return self._xvar_name

    def load_data(self, hist_data: TH1F, name:str='datahist'):
        
        xvar = self._roo_workspace.obj(self.xvar_name)
        self._hist_data = hist_data
        roo_data_hist = RooDataHist(name, name, [xvar], Import=hist_data)
        getattr(self._roo_workspace, 'import')(roo_data_hist)
        self._roo_data_hist_name = name

    def save_to_workspace(self):
        raise NotImplementedError()
