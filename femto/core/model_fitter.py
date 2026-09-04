import numpy as np
from ROOT import (TFile, TCanvas, TH1F, TLegend, TPaveText, 
                 RooRealVar, RooFit, RooDataHist, RooWorkspace, 
                 RooAddPdf, RooArgList, RooPlot, RooAbsReal, TGraphErrors)
from torchic import AxisSpec
from torchic.utils.root import set_root_object
from torchic.utils.colors import get_color

import sys
sys.path.append('/home/galucia/Lithium4/femto')
from core.fitter import Fitter
from core.utils import write_params_to_text

class ModelFitter(Fitter):

    def __init__(self, name, xvar_spec: AxisSpec, outfile:TFile = None, signal_func_names:list = None, bkg_func_names:list = None, 
                 workspace:RooWorkspace = None, extended:bool=False, title:str=None):
        
        super().__init__(name, xvar_spec, outfile, workspace)

        self._model_pdf = None
        self.fractions = {}

        self._outdir = outfile.mkdir('model') if outfile is not None else None

        if signal_func_names is None or bkg_func_names is None:
            raise ValueError('At least one signal and one bkg name should be provided')
        
        self._signal_pdfs = {}
        self._bkg_pdfs = {}
        self._model = None
        self.title = title
        self._fit_result = None

        self._data_label = None

        self.REFERENCE_KSTAR_VALUE_FOR_BKG_NORMALIZATION = 0.31 # GeV/c - arbitrary value in the plateau region to perform the normalisation
        self.bkg_normalisations_at_reference_kstar = None

        self.REFERENCE_KSTAR_VALUE_FOR_SIGNAL_NORMALIZATION = 0.07 # GeV/c - arbitrary value in the region where the signal is expected to be dominant to perform the normalisation
        self.signal_normalisations_at_reference_kstar = None
        self.KSTAR_MAX_SIGNIFICANCE = 0.23 # GeV/c - arbitrary value in the region where the signal is expected to be dominant to perform the significance calculation
        self.KSTAR_MAX_SIGNIFICANCE_GS = 0.17 # GeV/c - arbitrary value in the region where the signal is expected to be dominant to perform the significance calculation (assuming a ground state only for li4)

        self._init_model(name, signal_func_names, bkg_func_names, extended)
        

    def _init_model(self, name:str, signal_func_names:list, bkg_func_names:list, extended:bool):

        pdf_list, fraction_list = RooArgList(), RooArgList()
        
        for signal_name in signal_func_names:
            self._signal_pdfs[signal_name] = self._roo_workspace.obj(signal_name)
            title = self._signal_pdfs[signal_name].GetTitle()
            if ';' in title:
                title = title.split(';')[0]
            self.fractions[signal_name] = RooRealVar(signal_name+'_frac', f'#it{{f}}_{{{title}}}', 0.5, 0., 1.)
            
            if extended:
                self.fractions[signal_name].setRange(0., 1e4)
                self.fractions[signal_name].setVal(0.2)
                self.fractions[signal_name].SetTitle(f'#it{{N}}_{{{title}}}')
            pdf_list.add(self._signal_pdfs[signal_name])
            fraction_list.add(self.fractions[signal_name])

        for ibkg_name, bkg_name in enumerate(bkg_func_names):
            self._bkg_pdfs[bkg_name] = self._roo_workspace.obj(bkg_name)
            pdf_list.add(self._bkg_pdfs[bkg_name])
            title = self._bkg_pdfs[bkg_name].GetTitle()
            if ';' in title:
                title = title.split(';')[0]

            if ibkg_name == len(bkg_func_names)-1 and not extended:
                continue
            self.fractions[bkg_name] = RooRealVar(bkg_name+'_frac', f'#it{{f}}_{{{title}}}', 0.5, 0., 1.)
            
            if extended:    
                self.fractions[bkg_name].setRange(0., 1e4)
                self.fractions[bkg_name].SetTitle(f'#it{{N}}_{{{title}}}')
            fraction_list.add(self.fractions[bkg_name])

        self._model_pdf = RooAddPdf(name, name, pdf_list, fraction_list)
        if self.title:
            self._model_pdf.SetTitle(self.title)

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
            self._model_pdf.chi2FitTo(datahist, *fit_options, PrintLevel=-1, Verbose=False)
        else:
            self._model_pdf.fitTo(datahist, *fit_options, PrintLevel=-1, Verbose=False)

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
        self.bkg_normalisations_at_reference_kstar = bkg_curve.interpolate(self.REFERENCE_KSTAR_VALUE_FOR_BKG_NORMALIZATION)

        xvar.setRange(old_limits[0], old_limits[1])

        del datahist

    def fit_model(self, h_data:TH1F, signal_name:str, use_chi2_fit_method:bool=True, norm_range:str=None,
                  data_label:str=None):

        if data_label is not None:
            self._data_label = data_label
        xvar = self._roo_workspace.obj(self._xvar_name)
        self._roo_data_hist = RooDataHist(h_data.GetName()+'_datahist', h_data.GetName()+'_datahist', [xvar], Import=h_data)

        fit_options = [RooFit.Save(), RooFit.SumW2Error(True), RooFit.Extended(True)]
        if norm_range is not None:
            fit_options.append(RooFit.NormRange(norm_range))
            #fit_options.append(RooFit.SumCoefRange(norm_range))

        if use_chi2_fit_method:
            self._fit_result = self._model_pdf.chi2FitTo(self._roo_data_hist, *fit_options, PrintLevel=-1, Verbose=False)
        else:
            self._fit_result = self._model_pdf.fitTo(self._roo_data_hist, *fit_options, PrintLevel=-1, Verbose=False)

        #for name, fraction in self.fractions.items():
        #    print(f'{name=}: {fraction.getVal()=}, {fraction=}')
        
        frame = xvar.frame(xvar.getMin(), xvar.getMax())
        self.plot_model(frame, self._roo_data_hist, 'fit_signal')

        sig_curve = frame.findObject(self._signal_pdfs[signal_name].GetName())
        self.signal_normalisations_at_reference_kstar = sig_curve.interpolate(self.REFERENCE_KSTAR_VALUE_FOR_SIGNAL_NORMALIZATION)
        
        if self._outdir:
            graph = TGraphErrors(1)
            graph.SetPoint(0, 0.5, self.fractions[signal_name].getVal())
            graph.SetPointError(0, 0., self.fractions[signal_name].getError())
            graph.SetMarkerStyle(20)
            self._outdir.cd()
            graph.Write('signal_fraction')
            self._fit_result.Write('fit_result')

    def plot_model(self, frame:RooPlot, roodatahist:RooDataHist, canvas_name:str):

        roodatahist.plotOn(frame, MarkerStyle=20, LineColor=797, MarkerColor=797, MarkerSize=1.7) #, FillColorAlpha=(797, 0.3))
        
        line_color = get_color(0)
        for signal in self._signal_pdfs.values():
            self._model_pdf.plotOn(frame, Name=signal.GetName(),  Title= signal.GetTitle(), Normalization=(1.0, RooAbsReal.RelativeExpected), Components={signal}, LineColor=line_color)
            line_color += 1
        
        line_color = get_color(3)
        for bkg in self._bkg_pdfs.values():
            self._model_pdf.plotOn(frame, Name=bkg.GetName(),  Title= bkg.GetTitle(), Normalization=(1.0, RooAbsReal.RelativeExpected), Components={bkg}, LineColor=line_color)
            line_color += 1
        
        line_color = get_color(1)
        self._model_pdf.plotOn(frame, Name=self._model_pdf.GetName(),  Title=self._model_pdf.GetTitle(),  
                               Normalization=(1.0, RooAbsReal.RelativeExpected), LineColor=line_color)
                               
        #text = write_params_to_text(self.fractions.values(), coordinates=(0.5, 0.2, 0.75, 0.4), size=0.04)
        #text = TPaveText(0.48, 0.48, 0.81, 0.63, 'NDC')
        #text.SetFillColor(0)
        #text.SetBorderSize(0)
        #text.SetTextSize(0.044)
        #text.AddText(f'#bf{{{self.fractions["signal_pdf"].GetTitle()} = {self.fractions["signal_pdf"].getVal():.2f}}}')
        #text.AddText(f'#bf{{{self.fractions["bkg_pdf"].GetTitle()} = {self.fractions["bkg_pdf"].getVal():.2f} (fixed)}}')
        #text.AddText(f'#bf{{#chi^{{2}} / ndf = {frame.chiSquare():.2f}}}')
        #frame.addObject(text)

        legend = TLegend(0.52, 0.18, 0.81, 0.33)
        legend.SetBorderSize(0)
        legend.SetTextSize(0.045)
        h_dummy = TH1F('h_dummy', ';#it{k}* (GeV/#it{c}); C(k*)', 1, 0, 1)
        if self._data_label is not None:
            set_root_object(h_dummy, marker_color=1, line_color=1, 
                            marker_style=20, marker_size=1.7, fill_color_alpha=(797, 0.3))
            legend.AddEntry(h_dummy, self._data_label, 'lep')
        
        legend.AddEntry(frame.findObject(self._model_pdf.GetName()), self._model_pdf.GetTitle(), 'l')
        for signal in self._signal_pdfs.values():
            legend.AddEntry(frame.findObject(signal.GetName()), signal.GetTitle(), 'l')
        for bkg in self._bkg_pdfs.values():
            legend.AddEntry(frame.findObject(bkg.GetName()), bkg.GetTitle(), 'l')
        
        canvas = TCanvas(canvas_name)
        frame.Draw()
        #text.Draw('same')
        legend.Draw()

        if self._outdir:
            self._outdir.cd()
            canvas.Write()
        
        h_dummy.SetDirectory(0)
        del h_dummy

    def save_to_workspace(self):

        getattr(self._roo_workspace, 'import')(self._model_pdf)

    def get_bkg_value_at_reference_kstar(self):
        
        #stored_signal_fractions = {signal_name: self.fractions[signal_name].getVal() for signal_name in self._signal_pdfs.keys()}
        #xvar = self._roo_workspace.obj(self._xvar_name)
        #
        #for signal_name in self._signal_pdfs.keys():
        #    self.fractions[signal_name].setVal(0.)
        #
        #xvar.setVal(self.REFERENCE_KSTAR_VALUE_FOR_BKG_NORMALIZATION)
        #bkg_value_at_reference_kstar = self._model_pdf.getVal(xvar)
        #
        #for signal_name in self._signal_pdfs.keys():
        #    self.fractions[signal_name].setVal(stored_signal_fractions[signal_name])
        #    
        #return bkg_value_at_reference_kstar
        return self.bkg_normalisations_at_reference_kstar
        
    
    def compute_chi2(self, h_data:TH1F, h_systematics:TH1F=None, suffix:str='') -> float:
        '''
            Compute the chi2 for data against a background histogram
        '''

        chi2, chi2_model, ndf, nsigma, nsigma_syst = 0, 0, 0, 0, 0
        h_chi2 = TH1F(f'chi2{suffix}', ';#it{k}* (GeV/#it{c}); #chi^{2}', h_data.GetNbinsX(), h_data.GetBinLowEdge(1), h_data.GetBinLowEdge(h_data.GetNbinsX()+1))
        h_chi2_model = TH1F(f'chi2_model{suffix}', ';#it{k}* (GeV/#it{c}); #chi^{2} (model)', h_data.GetNbinsX(), h_data.GetBinLowEdge(1), h_data.GetBinLowEdge(h_data.GetNbinsX()+1))
        h_chi2_ndf = TH1F(f'chi2_ndf{suffix}', ';#it{k}* (GeV/#it{c}); #chi^{2} / NDF', h_data.GetNbinsX(), h_data.GetBinLowEdge(1), h_data.GetBinLowEdge(h_data.GetNbinsX()+1))
        
        h_nsigma = TH1F(f'nsigma{suffix}', ';#it{k}* (GeV/#it{c}); n#sigma', h_data.GetNbinsX(), h_data.GetBinLowEdge(1), h_data.GetBinLowEdge(h_data.GetNbinsX()+1))
        h_nsigma_syst = TH1F(f'nsigma_syst{suffix}', ';#it{k}* (GeV/#it{c}); n#sigma (systematics)', h_data.GetNbinsX(), h_data.GetBinLowEdge(1), h_data.GetBinLowEdge(h_data.GetNbinsX()+1))
        h_bkg_check = TH1F(f'bkg_check{suffix}', ';#it{k}* (GeV/#it{c}); C(k*)', h_data.GetNbinsX(), h_data.GetBinLowEdge(1), h_data.GetBinLowEdge(h_data.GetNbinsX()+1))
        
        h_nsigma_model = TH1F(f'nsigma_model{suffix}', ';#it{k}* (GeV/#it{c}); n#sigma (model)', h_data.GetNbinsX(), h_data.GetBinLowEdge(1), h_data.GetBinLowEdge(h_data.GetNbinsX()+1))
        h_nsigma_model_syst = TH1F(f'nsigma_model_syst{suffix}', ';#it{k}* (GeV/#it{c}); n#sigma (model + systematics)', h_data.GetNbinsX(), h_data.GetBinLowEdge(1), h_data.GetBinLowEdge(h_data.GetNbinsX()+1))
        h_model_check = TH1F(f'model_check{suffix}', ';#it{k}* (GeV/#it{c}); C(k*)', h_data.GetNbinsX(), h_data.GetBinLowEdge(1), h_data.GetBinLowEdge(h_data.GetNbinsX()+1))
        
        xvar = self._roo_workspace.obj(self._xvar_name)
        stored_signal_fractions = {signal_name: self.fractions[signal_name].getVal() for signal_name in self._signal_pdfs.keys()}

        for ibin in range(1, h_data.GetNbinsX()+1):
            
            kstar_value = h_data.GetBinCenter(ibin)
            if kstar_value < xvar.getMin() or kstar_value > xvar.getMax():
                continue
            data_value = h_data.GetBinContent(ibin)
            data_error = h_data.GetBinError(ibin)
            syst_error = h_systematics.GetBinError(ibin) if h_systematics is not None else 0

            xvar.setVal(self.REFERENCE_KSTAR_VALUE_FOR_BKG_NORMALIZATION)
            bkg_value_at_reference_kstar = self._model_pdf.getVal(xvar)
            correction = self.bkg_normalisations_at_reference_kstar / bkg_value_at_reference_kstar

            xvar.setVal(kstar_value)
            model_value = self._model_pdf.getVal() * correction
            model_error = 0

            difference = data_value - model_value
            uncertainty = np.sqrt(data_error*data_error +  model_error*model_error)
            nsigma = difference / uncertainty if uncertainty > 0 else 0
            uncertainty_syst = np.sqrt(syst_error*syst_error + model_error*model_error + data_error*data_error)
            nsigma_syst = difference / uncertainty_syst if uncertainty_syst > 0 else 0
            chi2_model += nsigma * nsigma

            h_chi2_model.SetBinContent(ibin, chi2_model)
            h_nsigma_model.SetBinContent(ibin, nsigma)
            h_nsigma_model_syst.SetBinContent(ibin, nsigma_syst)
            
            h_model_check.SetBinContent(ibin, model_value)
            h_model_check.SetBinError(ibin, model_error)
        
        
        ### Set signal to zero and compute the background-only chi2
        for signal_name in self._signal_pdfs.keys():
            self.fractions[signal_name].setVal(0.)
        #with open('debug_model_fit.txt', 'w') as debug_file:
            #debug_file.write('#kstar\tdata_value\tdata_error\tbkg_value\tbkg_error\tdifference\tuncertainty\tnsigma\n')

        for ibin in range(1, h_data.GetNbinsX()+1):
            
            kstar_value = h_data.GetBinCenter(ibin)
            if kstar_value < xvar.getMin() or kstar_value > xvar.getMax():
                continue
            data_value = h_data.GetBinContent(ibin)
            data_error = h_data.GetBinError(ibin)
            syst_error = h_systematics.GetBinError(ibin) if h_systematics is not None else 0

            xvar.setVal(self.REFERENCE_KSTAR_VALUE_FOR_BKG_NORMALIZATION)
            bkg_value_at_reference_kstar = self._model_pdf.getVal(xvar)
            correction = self.bkg_normalisations_at_reference_kstar / bkg_value_at_reference_kstar

            xvar.setVal(kstar_value)
            bkg_value = self._model_pdf.getVal() * correction
            bkg_error = 0

            difference = data_value - bkg_value
            uncertainty = np.sqrt(data_error*data_error +  bkg_error*bkg_error)
            nsigma = difference / uncertainty if uncertainty > 0 else 0
            uncertainty_syst = np.sqrt(syst_error*syst_error + bkg_error*bkg_error + data_error*data_error)
            nsigma_syst = difference / uncertainty_syst if uncertainty_syst > 0 else 0
            
            chi2 += nsigma * nsigma
            ndf += 1

            #debug_file.write(f'{kstar_value:.4f}\t{data_value:.4f}\t{data_error:.4f}\t{bkg_value:.4f}\t{bkg_error:.4f}\t{difference:.4f}\t{uncertainty:.4f}\t{nsigma:.4f}\n')

            h_chi2.SetBinContent(ibin, chi2)
            h_chi2_ndf.SetBinContent(ibin, chi2/ndf)
            h_nsigma.SetBinContent(ibin, nsigma)
            h_nsigma_syst.SetBinContent(ibin, nsigma_syst)
            
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

        if self._outdir:
            self._outdir.cd()
            for obj in [h_data, h_bkg_check, h_chi2, h_nsigma, h_nsigma_syst, h_chi2_model,
                        h_chi2_ndf, h_bkg_check, h_nsigma_model, h_nsigma_model_syst, h_model_check, canvas]:
                obj.Write()
        for obj in [h_bkg_check, h_chi2, h_nsigma, h_nsigma_syst, h_chi2_model,
                    h_chi2_ndf, h_bkg_check, h_nsigma_model, h_nsigma_model_syst, h_model_check]:
            obj.SetDirectory(0)
            del obj
        del canvas
        
        return chi2_model
        
    def compute_raw_yield(self, h_same_event:TH1F, h_mixed_event:TH1F, signal_pdf_name:str, bkg_pdf_name:str):

        nsig_stored, nbkg_stored = self.fractions[signal_pdf_name].getVal(), self.fractions[bkg_pdf_name].getVal()
        h_signal_correlation = h_mixed_event.Clone('h_signal_correlation')
        h_background_correlation = h_mixed_event.Clone('h_background_correlation')
        h_same_event_signal = h_mixed_event.Clone('h_same_event_signal')
        h_same_event_total = h_mixed_event.Clone('h_same_event_total')
        h_signal_bin_counting = h_mixed_event.Clone('h_signal_binning')
        nbins = h_mixed_event.GetNbinsX()
        bkg_fraction_error = self.fractions[bkg_pdf_name].getError()

        xvar = self._roo_workspace.obj(self._xvar_name)
        xvar.setVal(self.REFERENCE_KSTAR_VALUE_FOR_BKG_NORMALIZATION)
        self.fractions[signal_pdf_name].setVal(0.)
        self.fractions[bkg_pdf_name].setVal(nbkg_stored)
        bkg_value_at_reference_kstar = self._model_pdf.getVal(xvar)
        bkg_correction = self.bkg_normalisations_at_reference_kstar / bkg_value_at_reference_kstar if bkg_value_at_reference_kstar > 0 else 1.0

        xvar.setVal(self.REFERENCE_KSTAR_VALUE_FOR_SIGNAL_NORMALIZATION)
        self.fractions[signal_pdf_name].setVal(nsig_stored)
        self.fractions[bkg_pdf_name].setVal(0.)
        signal_value_at_reference_kstar = self._model_pdf.getVal(xvar)
        signal_correction = self.signal_normalisations_at_reference_kstar / signal_value_at_reference_kstar if signal_value_at_reference_kstar > 0 else 1.0

        for ibin in range(1, nbins+1):
            kstar_value = h_mixed_event.GetBinCenter(ibin)
            me_value = h_mixed_event.GetBinContent(ibin)
            se_value = h_same_event.GetBinContent(ibin)
            xvar.setVal(kstar_value)

            self.fractions[signal_pdf_name].setVal(0.)
            self.fractions[bkg_pdf_name].setVal(nbkg_stored)
            bkg_value = self._model_pdf.getVal() * bkg_correction
            bkg_error = bkg_fraction_error * bkg_value

            self.fractions[signal_pdf_name].setVal(nsig_stored)
            self.fractions[bkg_pdf_name].setVal(0.)
            signal_value = self._model_pdf.getVal() * signal_correction
            total_value = self._hist_data.GetBinContent(ibin) if self._hist_data is not None else 0.
            signal_subtraction_value = (total_value - bkg_value) * me_value
            signal_subtraction_value_error = (np.sqrt(se_value + bkg_error**2 * me_value**2 + me_value*bkg_value**2)
                                              if se_value + bkg_error**2 * me_value**2 + me_value*bkg_value**2 > 0 else 0)

            h_signal_correlation.SetBinContent(ibin, signal_value)
            h_background_correlation.SetBinContent(ibin, bkg_value)
            h_same_event_signal.SetBinContent(ibin, signal_value * me_value)
            h_same_event_total.SetBinContent(ibin, (signal_value + bkg_value) * me_value)
            h_signal_bin_counting.SetBinContent(ibin, signal_subtraction_value)
            h_signal_bin_counting.SetBinError(ibin, signal_subtraction_value_error)

        last_bin = h_same_event_signal.FindBin(xvar.getMax())
        yield_value = h_same_event_signal.Integral(1, last_bin)
        LAST_BIN_SIGNIFICANCE = h_same_event_total.FindBin(self.KSTAR_MAX_SIGNIFICANCE)
        signal_value = h_same_event_signal.Integral(1, LAST_BIN_SIGNIFICANCE)
        total_value = h_same_event_total.Integral(1, LAST_BIN_SIGNIFICANCE)
        significance = signal_value / np.sqrt(total_value) if total_value > 0 else 0
        
        LAST_BIN_SIGNIFICANCE_GS = h_same_event_total.FindBin(self.KSTAR_MAX_SIGNIFICANCE_GS)
        signal_value_gs = h_same_event_signal.Integral(1, LAST_BIN_SIGNIFICANCE_GS)
        total_value_gs = h_same_event_total.Integral(1, LAST_BIN_SIGNIFICANCE_GS)
        significance_gs = signal_value_gs / np.sqrt(total_value_gs) if total_value_gs > 0 else 0
        
        signal_value_bin_counting, signal_error_bin_counting = 0, 0
        for ibin in range(1, h_signal_bin_counting.FindBin(0.15)+1): # LAST_BIN_SIGNIFICANCE+1):
            signal_value_bin_counting += h_signal_bin_counting.GetBinContent(ibin)
            signal_error_bin_counting += h_signal_bin_counting.GetBinError(ibin)**2
        signal_error_bin_counting = np.sqrt(signal_error_bin_counting)

        canvas = TCanvas('yield_extraction', '')
        set_root_object(h_same_event_signal, marker_style=20, marker_color=797, line_color=797, 
                        title='Same-event signal; #it{k}* (GeV/#it{c}); C(k*)')

        text = TPaveText(0.5, 0.5, 0.8, 0.7, 'ndc')
        set_root_object(text, fill_style=0, border_size=0, text_size=0.04)
        text.AddText(f'Raw yield = {yield_value:.2f}')
        text.AddText(f'Significance = {significance:.2f}')
        text.AddText(f'Significance (GS) = {significance_gs:.2f}')

        h_same_event_signal.Draw('hist')
        text.Draw('same')
        
        canvas_bin_counting = TCanvas('yield_extraction_bin_counting', '')
        set_root_object(h_signal_bin_counting, marker_style=20, marker_color=797, line_color=797, 
                        title='Same-event signal (bin counting); #it{k}* (GeV/#it{c}); C(k*)')
        h_signal_bin_counting.Draw('hist')
        
        text_bin_counting = TPaveText(0.5, 0.5, 0.8, 0.7, 'ndc')
        set_root_object(text_bin_counting, fill_style=0, border_size=0, text_size=0.04)
        text_bin_counting.AddText(f'Raw yield = {signal_value_bin_counting:.2f} #pm {signal_error_bin_counting:.2f}')
        text_bin_counting.Draw('same')
        
        
        if self._outdir:
            self._outdir.cd()
            h_same_event_signal.Write()
            h_same_event_total.Write()
            canvas.Write()
            
            h_signal_bin_counting.Write()
            canvas_bin_counting.Write()
        
        for h in (h_signal_correlation, h_background_correlation, h_same_event_signal, h_same_event_total):
            del h

        return yield_value
    
    def cleanup(self, keep_histograms:bool=False):
        self._model_pdf = None
        self._signal_pdfs = {}
        self._bkg_pdfs = {}
        self.fractions = {}
        if hasattr(self, '_roo_data_hist'):
            self._roo_data_hist = None
        super().cleanup(keep_histograms=keep_histograms)
        