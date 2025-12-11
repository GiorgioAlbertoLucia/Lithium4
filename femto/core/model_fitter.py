import numpy as np
from ROOT import TFile, TCanvas, TH1F, TLegend, \
                 RooRealVar, RooCrystalBall, RooFit, RooHistPdf, RooDataHist, RooWorkspace, RooAddPdf, RooArgList, RooPlot, RooAbsReal
from torchic import AxisSpec
from torchic.utils.root import set_root_object

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

        self.reference_kstar_value_for_normalisation = 0.31 # GeV/c - arbitrary value in the plateau region to perform the normalisation
        self.bkg_normalisations_at_reference_kstar = None

        self._init_model(name, signal_func_names, bkg_func_names, extended)
        

    def _init_model(self, name:str, signal_func_names:list, bkg_func_names:list, extended:bool):

        pdf_list, fraction_list = RooArgList(), RooArgList()
        for signal_name in signal_func_names:
            self._signal_pdfs[signal_name] = self._roo_workspace.obj(signal_name)
            self.fractions[signal_name] = RooRealVar(signal_name+'_frac', f'#it{{f}}_{{{signal_name}}}', 0.5, 0., 1.)
            if extended:
                self.fractions[signal_name].setRange(0., 1e4)
                self.fractions[signal_name].setVal(0.2)
            pdf_list.add(self._signal_pdfs[signal_name])
            fraction_list.add(self.fractions[signal_name])

        for ibkg_name, bkg_name in enumerate(bkg_func_names):
            self._bkg_pdfs[bkg_name] = self._roo_workspace.obj(bkg_name)
            pdf_list.add(self._bkg_pdfs[bkg_name])
            if ibkg_name == len(bkg_func_names)-1 and not extended:
                continue
            self.fractions[bkg_name] = RooRealVar(bkg_name+'_frac', f'#it{{f}}_{{{bkg_name}}}', 0.5, 0., 1.)
            if extended:    self.fractions[bkg_name].setRange(0., 1e4)
            fraction_list.add(self.fractions[bkg_name])

        self._model_pdf = RooAddPdf(name, name, pdf_list, fraction_list)

    @property
    def model_pdf(self):
        return self._model_pdf
    
    def prefit_background(self, h_data:TH1F, range_limits:tuple=None, range_name:str=None,
                          use_chi2_method:bool=True, save_normalisation_value=True):
        
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
            #fit_options.append(RooFit.SumCoefRange(range_name))

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

        for bkg_name in self._bkg_pdfs.keys():
            if bkg_name in self.fractions.keys():

                bkg_pdf = self._bkg_pdfs[bkg_name]
                I_full = bkg_pdf.createIntegral([xvar], NormSet=[xvar], Range='full').getVal()
                I_side = bkg_pdf.createIntegral([xvar], NormSet=[xvar], Range='prefit').getVal()
                bkg_full_correction = (I_full / I_side)

                norm = self.fractions[bkg_name].getVal()
                self.fractions[bkg_name].setVal(norm * bkg_full_correction)
                self.fractions[bkg_name].setConstant(True)

        for signal_name in self._signal_pdfs.keys():
            if signal_name in self.fractions.keys():
                self.fractions[signal_name].setConstant(False)

        bkg_curve = frame.findObject(self._model_pdf.GetName())
        self.bkg_normalisations_at_reference_kstar = bkg_curve.interpolate(self.reference_kstar_value_for_normalisation)

        xvar.setRange(old_limits[0], old_limits[1])

        del datahist

    def fit_model(self, h_data:TH1F, use_chi2_fit_method:bool=True, norm_range:str=None):

        xvar = self._roo_workspace.obj(self._xvar_name)
        self._roo_data_hist = RooDataHist(h_data.GetName()+'_datahist', h_data.GetName()+'_datahist', [xvar], Import=h_data)

        fit_options = [RooFit.Save(), RooFit.SumW2Error(True), RooFit.Extended(True)]
        if norm_range is not None:
            fit_options.append(RooFit.NormRange(norm_range))
            #fit_options.append(RooFit.SumCoefRange(norm_range))

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
        self._model_pdf.plotOn(frame, Name=self._model_pdf.GetName(),  Title=self._model_pdf.GetTitle(),  Normalization=(1.0, RooAbsReal.RelativeExpected), LineColor=line_color)
        line_color += 1
        for signal in self._signal_pdfs.values():
            self._model_pdf.plotOn(frame, Name=signal.GetName(),  Title= signal.GetTitle(), Normalization=(1.0, RooAbsReal.RelativeExpected), Components={signal}, LineColor=line_color)
            line_color = 1
        for bkg in self._bkg_pdfs.values():
            self._model_pdf.plotOn(frame, Name=bkg.GetName(),  Title= bkg.GetTitle(), Normalization=(1.0, RooAbsReal.RelativeExpected), Components={bkg}, LineColor=line_color)
            line_color += 1

        text = write_params_to_text(self.fractions.values(), coordinates=(0.5, 0.2, 0.75, 0.4), size=0.04)
        text.AddText(f'#chi^{{2}} / ndf = {frame.chiSquare():.2f}')
        #frame.addObject(text)

        legend = TLegend(0.5, 0.42, 0.75, 0.64)
        legend.SetBorderSize(0)
        legend.SetTextSize(0.04)
        legend.AddEntry(frame.findObject(self._model_pdf.GetName()), 'Total', 'l')
        for signal in self._signal_pdfs.values():
            legend.AddEntry(frame.findObject(signal.GetName()), signal.GetTitle(), 'l')
        for bkg in self._bkg_pdfs.values():
            legend.AddEntry(frame.findObject(bkg.GetName()), bkg.GetTitle(), 'l')

        self._outdir.cd()
        canvas = TCanvas(canvas_name)
        frame.Draw()
        text.Draw('same')
        legend.Draw()
        canvas.Write()

    def save_to_workspace(self):

        getattr(self._roo_workspace, 'import')(self._model_pdf)

    def compute_chi2(self, h_data:TH1F):
        '''
            Compute the chi2 for data against a background histogram
        '''

        chi2, ndf, nsigma = 0, 0, 0
        h_chi2 = TH1F('chi2', ';#it{k}* (GeV/#it{c}); #chi^{2}', h_data.GetNbinsX(), h_data.GetBinLowEdge(1), h_data.GetBinLowEdge(h_data.GetNbinsX()+1))
        h_nsigma = TH1F('nsigma', ';#it{k}* (GeV/#it{c}); n#sigma', h_data.GetNbinsX(), h_data.GetBinLowEdge(1), h_data.GetBinLowEdge(h_data.GetNbinsX()+1))
        h_bkg_check = TH1F('bkg_check', ';#it{k}* (GeV/#it{c}); C(k*)', h_data.GetNbinsX(), h_data.GetBinLowEdge(1), h_data.GetBinLowEdge(h_data.GetNbinsX()+1))
        h_chi2_ndf = TH1F('chi2_ndf', ';#it{k}* (GeV/#it{c}); #chi^{2} / NDF', h_data.GetNbinsX(), h_data.GetBinLowEdge(1), h_data.GetBinLowEdge(h_data.GetNbinsX()+1))
        xvar = self._roo_workspace.obj(self._xvar_name)

        stored_signal_fractions = {signal_name: self.fractions[signal_name].getVal() for signal_name in self._signal_pdfs.keys()}
        for signal_name in self._signal_pdfs.keys():
            self.fractions[signal_name].setVal(0.)

        for ibin in range(1, h_data.GetNbinsX()+1):
            
            kstar_value = h_data.GetBinCenter(ibin)
            data_value = h_data.GetBinContent(ibin)
            data_error = h_data.GetBinError(ibin)

            xvar.setVal(self.reference_kstar_value_for_normalisation)
            bkg_value_at_reference_kstar = self._model_pdf.getVal(xvar)
            correction = self.bkg_normalisations_at_reference_kstar / bkg_value_at_reference_kstar

            xvar.setVal(kstar_value)
            bkg_value = self._model_pdf.getVal() * correction
            bkg_error = 0

            difference = data_value - bkg_value
            uncertainty = np.sqrt(data_error*data_error +  bkg_error*bkg_error)
            nsigma = difference/uncertainty
            chi2 += difference*difference/(uncertainty*uncertainty)
            ndf += 1

            h_chi2.SetBinContent(ibin, chi2)
            h_chi2_ndf.SetBinContent(ibin, chi2/ndf)
            h_nsigma.SetBinContent(ibin, nsigma)
            
            h_bkg_check.SetBinContent(ibin, bkg_value)
            h_bkg_check.SetBinError(ibin, bkg_error)

        for signal_name in self._signal_pdfs.keys():
            self.fractions[signal_name].setVal(stored_signal_fractions[signal_name])

        canvas = TCanvas('data_bkg_comparison', '')
        set_root_object(h_data, marker_style=20, marker_color=797, line_color=797, title='Data')
        set_root_object(h_bkg_check, marker_style=20, marker_color=420, line_color=420, title='Background')
        h_data.Draw('e1')
        h_bkg_check.Draw('e1 same')
        legend = canvas.BuildLegend(0.5, 0.3, 0.8, 0.5)
        legend.SetBorderSize(0)

        self._outdir.cd()
        h_chi2.Write()
        h_nsigma.Write()
        h_chi2_ndf.Write()
        h_bkg_check.Write()
        canvas.Write()


