import numpy as np
from ROOT import TFile, TCanvas, TH1F, TTree, \
                 RooRealVar, RooFit, RooHistPdf, RooDataHist, RooWorkspace, RooDataSet, RooKeysPdf, RooExtendPdf
from torchic import AxisSpec


import sys
sys.path.append('/home/galucia/Lithium4/femto')
from core.fitter import Fitter
from core.utils import write_params_to_text

class BkgFitter(Fitter):

    def __init__(self, name, xvar_spec: AxisSpec, outfile:TFile = None, workspace:RooWorkspace = None):
        super().__init__(name, xvar_spec, outfile, workspace)
        self._bkg_pars = {}
        self._bkg_pdf = None
        self._title = ''

        self._bkg_datahist = None
        self._bkg_normalisation = None

        self._outdir = self._outfile.mkdir('bkg') if self._outfile is not None else None

    @property
    def bkg_pars(self):
        return self._bkg_pars
    
    @property
    def bkg_normalisation(self):
        return self._bkg_normalisation
    
    @property
    def bkg_pdf(self):
        return self._bkg_pdf
    
    @property
    def title(self):
        return self._title
    
    @title.setter
    def title(self, title:str):
        self._title = title
        if self._bkg_pdf is not None:
            self._bkg_pdf.SetTitle(self._title)

    def _init_bkg_from_mc(self, h_bkg, name:str='bkg_pdf', extended:bool=False, **kwargs):
        
        xvar = self._roo_workspace.obj(self._xvar_name)
        
        for ibin in range(1, h_bkg.GetNbinsX()+1):
            bin_center = h_bkg.GetBinCenter(ibin)
            if not (xvar.getMin() <= bin_center <= xvar.getMax()):
                h_bkg.SetBinContent(ibin, 0)
        self._bkg_datahist = RooDataHist('bkg_dh', 'bkg_dh', [xvar], Import=h_bkg)
        
        if extended:
            bkg_pdf = RooHistPdf('tmp', 'tmp', [xvar], self._bkg_datahist)
            self._bkg_normalisation = RooRealVar('bkg_normalisation', '#it{N}_{bkg}', 1., 0, 1e4)
            self._bkg_pdf = RooExtendPdf(name, name, bkg_pdf, self._bkg_normalisation)
        else:
            self._bkg_pdf = RooHistPdf(name, name, [xvar], self._bkg_datahist)

        frame = xvar.frame()
        self._bkg_datahist.plotOn(frame)
        self._bkg_pdf.plotOn(frame)

        if self._outdir:
            self._outdir.cd()
            canvas = TCanvas('bkg_pdf')
            frame.Draw()
            canvas.Write()

            del canvas
    
    def _init_bkg_from_kde(self, h_bkg:TH1F, name:str='bkg_pdf', xmin:float=0.01, xmax:float=0.42, rho:float=0.05, **kwargs):
    
        xvar = self._roo_workspace.obj(self._xvar_name)
        old_range = (xvar.getMin(), xvar.getMax())
        xvar.setRange(xmin, xmax)

        x_data, weights = [], []
        for ibin in range(1, h_bkg.GetNbinsX()+1):
            x_val = h_bkg.GetBinCenter(ibin)
            if not (xmin <= x_val <= xmax):
                continue
            y_val = h_bkg.GetBinContent(ibin)
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
        self._bkg_dataset = RooDataSet(h_bkg.GetName()+'_roodata', h_bkg.GetName()+'_roodata', 
                            [xvar, weight_var],
                            RooFit.Import(tree),
                            RooFit.WeightVar(weight_var)) #, '', 'weight')
        getattr(self._roo_workspace, 'import')(self._bkg_dataset)
        
        self._bkg_pdf = RooKeysPdf(name, name, xvar, self._bkg_dataset, RooKeysPdf.NoMirror, rho)
        
        frame = xvar.frame(len(x_data))
        self._bkg_dataset.plotOn(frame, MarkerStyle=20, MarkerSize=0.8, LineColor=1)
        self._bkg_pdf.plotOn(frame, LineColor=2, LineWidth=2)
        canvas = TCanvas(f'cKeysPdf_{h_bkg.GetName()}', f'cKeysPdf_{h_bkg.GetName()}', 800, 600)
        frame.Draw()

        if self._outdir:
            self._outdir.cd()
            h_bkg.Write(f'{h_bkg.GetName()}_original')
            self._bkg_pdf.Write(f'{h_bkg.GetName()}_keyspdf')
            canvas.Write()

        xvar.setRange(*old_range)
        
        del tree
        
    def init_bkg(self, mode:str, *args, **kwargs):

        if mode == 'from_mc':   self._init_bkg_from_mc(*args, **kwargs)
        elif mode == 'from_kde':  self._init_bkg_from_kde(*args, **kwargs)
        else:                   raise ValueError('Only supported modes are "from_mc"')
        
        
    def fit_bkg(self, hist:TH1F, range_limits:tuple=None, range_name:str=None, use_chi2_method:bool=True, save_normalisation_value=True):
        
        xvar = self._roo_workspace.obj(self._xvar_name)
        old_limits = (xvar.getMin(), xvar.getMax())
        range_limits = range_limits if range_limits is not None else old_limits

        if range_limits is not None:
            xvar.setRange(range_limits[0], range_limits[1])
            if range_name is not None and range_limits is not None:
                xvar.setRange(range_name, range_limits[0], range_limits[1])

        fit_options = [RooFit.Save()] #, RooFit.Extended(True)]
        if range_name is not None:
            #fit_options.append(RooFit.NormRange(range_name))
            fit_options.append(RooFit.Range(range_name))
        #if range_limits is not None:
        #    fit_options.append(RooFit.Range(range_limits[0], range_limits[1]))

        datahist = RooDataHist(hist.GetName()+'_datahist', hist.GetName()+'_datahist', [xvar], Import=hist)
        if use_chi2_method:
            self._bkg_pdf.chi2FitTo(datahist, *fit_options, PrintLevel=-1, Verbose=False)
        else:
            self._bkg_pdf.fitTo(datahist, *fit_options, PrintLevel=-1, Verbose=False)

        xvar.setRange("full", range_limits[0], range_limits[1])
        integral_full = self._bkg_pdf.createIntegral([xvar], Range="full").getVal()
        if save_normalisation_value:
            self._bkg_normalisation.setVal(integral_full)

        frame = xvar.frame()
        datahist.plotOn(frame)
        self._bkg_pdf.plotOn(frame)

        text = write_params_to_text(self._bkg_pars.values(), coordinates=(0.7, 0.3, 0.85, 0.5))
        text.AddText(f'{self._bkg_normalisation.GetTitle()} = {self._bkg_normalisation.getVal():.2f} #pm {self._bkg_normalisation.getError():.2f}')
        text.AddText(f'#chi^{{2}} / ndf = {frame.chiSquare():.2f}')
        frame.addObject(text)

        if self._outdir:
            self._outdir.cd()
            canvas = TCanvas('fit_bkg')
            frame.Draw()
            canvas.Write()

            del canvas

        xvar.setRange(old_limits[0], old_limits[1])

        del datahist

    def save_to_workspace(self):

        getattr(self._roo_workspace, 'import')(self._bkg_pdf)
        for param in self._bkg_pars.values():
            getattr(self._roo_workspace, 'import')(param)
        if self._bkg_datahist is not None:
            getattr(self._roo_workspace, 'import')(self._bkg_datahist)
        if self._bkg_normalisation is not None:
            getattr(self._roo_workspace, 'import')(self._bkg_normalisation)

    def cleanup(self):
        self._bkg_pdf = None
        self._bkg_pars = {}
        self._bkg_datahist = None
        self._bkg_normalisation = None
        if hasattr(self, '_bkg_dataset'):
            self._bkg_dataset = None
        super().cleanup()