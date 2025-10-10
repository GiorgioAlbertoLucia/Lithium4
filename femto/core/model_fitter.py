from ROOT import TFile, TCanvas, TH1F, \
                 RooRealVar, RooCrystalBall, RooFit, RooHistPdf, RooDataHist, RooWorkspace, RooAddPdf, RooArgList, RooPlot
from torchic import AxisSpec

import sys
sys.path.append('/home/galucia/Lithium4/femto')
from core.fitter import Fitter
from core.utils import write_params_to_text

class ModelFitter(Fitter):

    def __init__(self, name, xvar_spec: AxisSpec, outfile:TFile, signal_func_names:list = None, bkg_func_names:list = None, workspace:RooWorkspace = None, extended:bool=False):
        
        super().__init__(name, xvar_spec, outfile, workspace)

        self._model_pdf = None
        self.fractions = {}

        self._outdir = outfile.mkdir('model')

        if signal_func_names is None or bkg_func_names is None:
            raise ValueError('At least one signal and one bkg name should be provided')
        
        self._signal_pdfs = {}
        self._bkg_pdfs = {}
        self._model = None

        self._init_model(name, signal_func_names, bkg_func_names, extended)
        

    def _init_model(self, name:str, signal_func_names:list, bkg_func_names:list, extended:bool):

        pdf_list, fraction_list = RooArgList(), RooArgList()
        for signal_name in signal_func_names:
            self._signal_pdfs[signal_name] = self._roo_workspace.obj(signal_name)
            self.fractions[signal_name] = RooRealVar(signal_name+'_frac', signal_name+'_frac', 0.5, 0., 1.)
            if extended:    self.fractions[signal_name].setRange(0., 1e4)
            pdf_list.add(self._signal_pdfs[signal_name])
            fraction_list.add(self.fractions[signal_name])

        for ibkg_name, bkg_name in enumerate(bkg_func_names):
            self._bkg_pdfs[bkg_name] = self._roo_workspace.obj(bkg_name)
            pdf_list.add(self._bkg_pdfs[bkg_name])
            if ibkg_name == len(bkg_func_names)-1 and not extended:
                continue
            self.fractions[bkg_name] = RooRealVar(bkg_name+'_frac', bkg_name+'_frac', 0.5, 0., 1.)
            if extended:    self.fractions[bkg_name].setRange(0., 1e4)
            fraction_list.add(self.fractions[bkg_name])

        self._model_pdf = RooAddPdf(name, name, pdf_list, fraction_list)

    
    @property
    def model_pdf(self):
        return self._model_pdf
    
    def prefit_background(self, h_data:TH1F,range_limits:tuple=None, range_name:str=None, use_chi2_method:bool=True, save_normalisation_value=True):
        
        xvar = self._roo_workspace.obj(self._xvar_name)
        old_limits = (xvar.getMin(), xvar.getMax())
        range_limits = range_limits if range_limits is not None else old_limits

        if range_limits is not None:
            xvar.setRange(range_limits[0], range_limits[1])
            if range_name is not None and range_limits is not None:
                xvar.setRange(range_name, range_limits[0], range_limits[1])

        fit_options = [RooFit.Save(), RooFit.Extended(True)]
        if range_name is not None:
            fit_options.append(RooFit.Range(range_name))

        for signal_name in self._signal_pdfs.keys():
            if signal_name in self.fractions.keys():
                self.fractions[signal_name].setVal(0)
                self.fractions[signal_name].setConstant(True)

        datahist = RooDataHist(h_data.GetName()+'_datahist', h_data.GetName()+'_datahist', [xvar], Import=h_data)
        if use_chi2_method:
            self._model_pdf.chi2FitTo(datahist, *fit_options)
        else:
            self._model_pdf.fitTo(datahist, *fit_options)

        frame = xvar.frame()
        self.plot_model(frame, datahist, 'prefit_bkg')

        xvar.setRange('full', old_limits[0], old_limits[1])
        xvar.setRange('prefit', range_limits[0], range_limits[1])
        I_full = self._model_pdf.createIntegral([xvar], Range='full').getVal()
        I_side = self._model_pdf.createIntegral([xvar], Range='prefit').getVal()
        bkg_full_correction = (I_full / I_side)

        for signal_name in self._signal_pdfs.keys():
            if signal_name in self.fractions.keys():
                self.fractions[signal_name].setConstant(False)
        for bkg_name in self._bkg_pdfs.keys():
            if bkg_name in self.fractions.keys():
                norm = self.fractions[bkg_name].getVal()
                self.fractions[bkg_name].setVal(norm * bkg_full_correction)
                self.fractions[bkg_name].setConstant(True)


        xvar.setRange(old_limits[0], old_limits[1])

        del datahist

    def fit_model(self, h_data:TH1F, use_chi2_fit_method:bool=True, norm_range:str=None):

        xvar = self._roo_workspace.obj(self._xvar_name)
        self._roo_data_hist = RooDataHist(h_data.GetName()+'_datahist', h_data.GetName()+'_datahist', [xvar], Import=h_data)

        fit_options = [RooFit.Save(), RooFit.SumW2Error(True), RooFit.Extended(True)]
        #if norm_range is not None:
            #fit_options.append(RooFit.NormRange(norm_range))

        if use_chi2_fit_method:
            self._model_pdf.chi2FitTo(self._roo_data_hist, *fit_options)
        else:
            self._model_pdf.fitTo(self._roo_data_hist, *fit_options)

        for name, fraction in self.fractions.items():
            print(f'{name=}: {fraction.getVal()=}, {fraction=}')
        
        frame = xvar.frame()
        self.plot_model(frame, self._roo_data_hist, 'fit_signal')

    def plot_model(self, frame:RooPlot, roodatahist:RooDataHist, canvas_name:str):

        line_color = 2
        roodatahist.plotOn(frame)
        self._model_pdf.plotOn(frame, LineColor=line_color)
        line_color += 1
        for signal in self._signal_pdfs.values():
            self._model_pdf.plotOn(frame, Components={signal}, LineColor=line_color)
            line_color = 1
        for bkg in self._bkg_pdfs.values():
            self._model_pdf.plotOn(frame, Components={bkg}, LineColor=line_color)
            line_color += 1

        text = write_params_to_text(self.fractions.values(), coordinates=(0.7, 0.3, 0.85, 0.5))
        text.AddText(f'#chi^{{2}} / ndf = {frame.chiSquare():.2f}')
        frame.addObject(text)

        self._outdir.cd()
        canvas = TCanvas(canvas_name)
        frame.Draw()
        canvas.Write()

    def save_to_workspace(self):

        getattr(self._roo_workspace, 'import')(self._model_pdf)

    