from ROOT import TFile, TCanvas, TH1F, \
                 RooRealVar, RooCrystalBall, RooFit, RooHistPdf, RooDataHist, RooWorkspace, RooExtendPdf
from torchic import AxisSpec


import sys
sys.path.append('/home/galucia/Lithium4/femto')
from core.fitter import Fitter
from core.utils import write_params_to_text

class BkgFitter(Fitter):

    def __init__(self, name, xvar_spec: AxisSpec, outfile:TFile, workspace:RooWorkspace = None):
        super().__init__(name, xvar_spec, outfile, workspace)
        self._bkg_pars = {}
        self._bkg_pdf = None

        self._bkg_datahist = None
        self._bkg_normalisation = None

        self._outdir = self._outfile.mkdir('bkg')

    @property
    def bkg_pars(self):
        return self._bkg_pars
    
    @property
    def bkg_normalisation(self):
        return self._bkg_normalisation
    
    @property
    def bkg_pdf(self):
        return self._bkg_pdf

    def _init_bkg_from_mc(self, h_bkg, name:str='bkg_pdf', extended:bool=False):
        
        xvar = self._roo_workspace.obj(self._xvar_name)
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

        self._outdir.cd()
        canvas = TCanvas('bkg_pdf')
        frame.Draw()
        canvas.Write()

        del canvas
        
    def init_bkg(self, mode:str, *args, **kwargs):

        if mode == 'from_mc':   self._init_bkg_from_mc(*args, **kwargs)
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
            self._bkg_pdf.chi2FitTo(datahist, *fit_options)
        else:
            self._bkg_pdf.fitTo(datahist, *fit_options)

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

        self._outdir.cd()
        canvas = TCanvas('fit_bkg')
        frame.Draw()
        canvas.Write()

        xvar.setRange(old_limits[0], old_limits[1])

        del datahist, canvas

    def save_to_workspace(self):

        getattr(self._roo_workspace, 'import')(self._bkg_pdf)
        for param in self._bkg_pars.values():
            getattr(self._roo_workspace, 'import')(param)
        if self._bkg_datahist is not None:
            getattr(self._roo_workspace, 'import')(self._bkg_datahist)
        if self._bkg_normalisation is not None:
            getattr(self._roo_workspace, 'import')(self._bkg_normalisation)

    