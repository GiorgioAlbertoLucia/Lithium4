'''
    Class to perform a fit to the correlation function
'''

import os
import numpy as np

from torchic.utils.terminal_colors import TerminalColors as tc
from ROOT import RooRealVar, RooDataHist, RooHistPdf, RooCrystalBall, RooFit, RooAddPdf, RooStats, RooWorkspace, \
    RooChi2Var, RooFormulaVar, RooAbsData, RooGenericPdf
from ROOT import TF1, TH1F, TCanvas, TDirectory, TGraph, gInterpreter, TPaveText, TPad, gStyle

from femto.utils import write_params_to_text

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
INCLUDE_DIR = os.path.join(CURRENT_DIR, 'crystal_ball.hh')
gInterpreter.ProcessLine(f'#include "{INCLUDE_DIR}"')
from ROOT import CrystalBall

def get_alice_watermark(xmin, ymin, xmax, ymax):
    '''
        Get the Alice watermark
    '''
    watermark = TPaveText(xmin, ymin, xmax, ymax, 'NDC')
    watermark.SetBorderSize(0)
    watermark.SetFillColor(0)
    watermark.AddText('This work')
    watermark.AddText('#bf{ALICE, Run 3}')
    watermark.AddText('#bf{Pb#topbarPb, #sqrt{s_{NN}} = 5.36 TeV}')
    return watermark

def get_centrality_watermark(xmin, ymin, xmax, ymax, centrality: str):
    '''
        Get the Alice watermark
    '''
    watermark = TPaveText(xmin, ymin, xmax, ymax, 'NDC')
    watermark.SetBorderSize(0)
    watermark.SetFillColor(0)

    if centrality == 'integrated':
        watermark.AddText('#bf{Centrality 0-50%}')
    elif centrality == '0-10':
        watermark.AddText('#bf{Centrality 0-10%}')
    elif centrality == '10-30':
        watermark.AddText('#bf{Centrality 10-30%}')
    elif centrality == '30-50':
        watermark.AddText('#bf{Centrality 30-50%}')

    return watermark

class FitWorkflowCF:

    def __init__(self,  outfile: TDirectory, pdf_outpath: str):

        gStyle.SetOptStat(0)
        
        self.outfile = outfile
        self.outdir = None
        self.pdf_outpath = pdf_outpath
        self.pdf_canvas = TCanvas('pdf_canvas', 'pdf_canvas', 800, 800)
        self.pdf_canvas.Print(f'{self.pdf_outpath}(')

        self.h_bkg = None  # Histogram with the background template
        self.h_CF = None
        self.CF_datahist = None
        self.CF_datahist_significance = None
        self.pull = None

        self.LOWER_LIMIT_FIT = 0.01  # 0.02 # Lower limit for the fit range in GeV/c
        self.UPPER_LIMIT_FIT = 0.4  # Upper limit for the fit range in GeV/c

        # Roofit
        self.kstar = RooRealVar('kstar', '#it{k}* (GeV/#it{c})', self.LOWER_LIMIT_FIT, self.UPPER_LIMIT_FIT)
        self.kstar_significance = RooRealVar('kstar_significance', 'kstar_significance', self.LOWER_LIMIT_FIT, self.UPPER_LIMIT_FIT)
        
        self.sig_counts = RooRealVar('sig_counts', 'sig_counts', 10, 0., 1e3)
        self.signal_pdf = None
        self.sig_params = {}
        self.bkg_counts = RooRealVar('bkg_counts', 'bkg_counts', 0., 1e6)
        self.bkg_datahist = None
        self.bkg_pdf = None
        self.model = None

        self.sign = None  # 'Matter' or 'Antimatter'

        self.roows = None # used for the asymptotic calculator

    def make_directory(self, name: str):
        '''
            Create a directory in the output file
        '''
        self.outdir = self.outfile.mkdir(name)

    def load_hist_data(self, h_CF: TH1F):
        '''
            Load the histogram with the data (signal + bkg)
        '''

        self.h_CF = h_CF.Clone('h_data')
        self.CF_datahist = RooDataHist('datahist', 'datahist', [self.kstar], Import=h_CF)
        self.CF_datahist_significance = RooDataHist('datahist_significance', 'datahist_significance', [self.kstar_significance], Import=h_CF)
        self.bkg_counts.setVal(h_CF.Integral())

    def draw_on_canvas(self, obj):
        self.pdf_canvas.Clear()
        self.pdf_canvas.cd()
        obj.Draw()
        self.pdf_canvas.Print(self.pdf_outpath)

    def draw_cf_canvas(self, frame, logy:bool) -> TCanvas:
        '''
            Create a canvas for the correlation function
        '''

        if self.pull is None:
            self.pull_from_bkg_alt()

        canvas = TCanvas('cf_canvas', 'cf_canvas', 1000, 1000)
        pad1 = TPad('pad1', 'pad1', 0, 0.3, 1, 1)
        pad1.SetBottomMargin(0.005)
        pad1.Draw()
        pad1.cd()
        frame.Draw()
        if logy:
            pad1.SetLogy()

        pad2 = TPad('pad2', 'pad2', 0, 0, 1, 0.3)
        pad2.SetTopMargin(0.0)
        canvas.cd()
        pad2.SetBottomMargin(0.25)
        pad2.Draw()
        pad2.cd()
        self.pull.GetXaxis().SetLabelSize(.08)
        self.pull.GetYaxis().SetLabelSize(.08)
        self.pull.GetXaxis().SetTitleSize(.08)
        self.pull.GetYaxis().SetTitleSize(.1)
        self.pull.GetYaxis().SetTitleOffset(0.5)
        self.pull.Draw('E1 same')
        zero_line = TF1('zero_line', '0', 0, self.kstar.getMax())
        zero_line.SetLineColor(1)
        zero_line.SetLineStyle(2)
        zero_line.Draw('same')

        canvas.SetTopMargin(0.05)
        canvas.SetBottomMargin(0.15)
        canvas.Print(self.pdf_outpath)

        del canvas

    def prepare_signal_fit(self, h_signal: TH1F):
        '''
            Fit the signal peak with a Crystal Ball function
        '''

        self.sig_params = {
            'sig_mean': RooRealVar('sig_mean', '#mu', 0.081, 0.06, 0.1, 'GeV/c'),
            'sig_sigma': RooRealVar('sig_sigma', '#sigma', 0.01, 0.001, 0.1, 'GeV/c'),
            'sig_aL': RooRealVar('sig_aL', 'a_{L}', 1.3, 0.1, 10.),
            'sig_nL': RooRealVar('sig_nL', 'n_{L}', 2.7, 0.1, 10.),
            'sig_aR': RooRealVar('sig_aR', 'a_{R}', 1.1, 0.1, 10.),
            'sig_nR': RooRealVar('sig_nR', 'n_{R}', 5.7, 0.1, 10.),
        }
        self.signal_pdf = RooCrystalBall('signal_pdf', 'signal_pdf', self.kstar, self.sig_params['sig_mean'], self.sig_params['sig_sigma'], self.sig_params['sig_aL'], 
                                self.sig_params['sig_nL'], self.sig_params['sig_aR'], self.sig_params['sig_nR'])

        datahist = RooDataHist('signal_dh', 'signal_dh', [self.kstar], Import=h_signal)
        
        print(f'\n{tc.BOLD+tc.GREEN}Fit to the signal peak {tc.RESET}')
        self.signal_pdf.fitTo(datahist, RooFit.Save(), RooFit.Range(0.01, 0.25), SumW2Error=True, Extended=True)
        #self.signal_pdf.chi2FitTo(datahist, RooFit.Save(), RooFit.Range(0., 0.15))
        frame = self.kstar.frame()
        datahist.plotOn(frame)
        self.signal_pdf.plotOn(frame)

        text = write_params_to_text(self.sig_params.values(), coordinates=(0.7, 0.3, 0.9, 0.5))
        frame.addObject(text)

        self.outdir.cd()
        h_signal.Write(f'signal_hist')
        frame.Write(f'signal_mc')

    def prepare_bkg_template(self, h_bkg: TH1F):
        '''
            Get the Coulomb template for the p-He3 correlation function
        '''

        self.h_bkg = h_bkg.Clone('h_bkg')
        self.bkg_datahist = RooDataHist('bkg_dh', 'bkg_dh', [self.kstar], Import=h_bkg)
        self.bkg_pdf = RooHistPdf('bkg_pdf', 'bkg_pdf', [self.kstar], self.bkg_datahist)

        frame = self.kstar.frame()
        self.bkg_datahist.plotOn(frame)
        self.bkg_pdf.plotOn(frame)

        self.outdir.cd()
        frame.Write(f'bkg_pdf')

    def pull_from_bkg(self, h_bkg: TH1F):
        '''
            Pull the Coulomb template from the histogram
        '''

        self.pull = self.h_CF.Clone('pull')
        
        for ibin in range(1, self.h_CF.GetNbinsX() + 1):
            ikstar = self.h_CF.GetBinCenter(ibin)
            imeasured = self.h_CF.GetBinContent(ibin)
            ierror = self.h_CF.GetBinError(ibin)
            iexpected = h_bkg.GetBinContent(h_bkg.FindBin(ikstar))
            ierr_bkg = h_bkg.GetBinError(h_bkg.FindBin(ikstar))
            ipull = 0. if ierror < 1.e-9 else (imeasured - iexpected) / np.sqrt(ierror**2 + ierr_bkg**2)
            self.pull.SetBinContent(ibin, ipull)
            self.pull.SetBinError(ibin, 1.)

        self.pull.SetTitle(';#it{k}* (GeV/#it{c}); n_{#sigma}')
        self.pull.SetMarkerStyle(20)
        self.pull.SetMarkerSize(1.5)
        self.pull.SetMarkerColor(418)

        self.outdir.cd()
        self.pull.Write(f'pull')

    def pull_from_bkg_alt(self):
        '''
            Pull the Coulomb template from the histogram
        '''

        _sig_counts = self.sig_counts.getVal()
        self.sig_counts.setVal(0.)

        _h_CF = self.h_CF.Clone('h_CF_clone')
        bin_width = _h_CF.GetBinWidth(1)
        _h_CF.GetXaxis().SetRangeUser(self.LOWER_LIMIT_FIT, self.UPPER_LIMIT_FIT)
        self.pull = _h_CF.Clone('pull')

        chi2_value = 0
        ndf = 0
        chi2_far_from_signal = 0
        ndf_far_from_signal = 0
        running_chi2 = []

        for ibin in range(3, _h_CF.GetNbinsX() + 1):
            ikstar = _h_CF.GetBinCenter(ibin)
            imeasured = _h_CF.GetBinContent(ibin)
            ierror = _h_CF.GetBinError(ibin)

            self.kstar.setVal(ikstar)
            iexpected = self.model.getVal(self.kstar) * self.bkg_counts.getVal() * bin_width
            ierr_bkg = 0
            ipull = 0. if ierror < 1.e-9 else (imeasured - iexpected) / np.sqrt(ierror**2 + ierr_bkg**2)
            self.pull.SetBinContent(ibin, ipull)
            self.pull.SetBinError(ibin, 1.)

            chi2_value_bin = (imeasured - iexpected)**2 / (ierror**2 + ierr_bkg**2) if (ierror**2 + ierr_bkg**2) > 0 else 0
            print(f'{(tc.BOLD+tc.RED) if chi2_value_bin > 2 else ""}kstar: {ikstar:.4f}, measured: {imeasured:.4f}, expected: {iexpected:.4f}, error: {ierror:.4f}, pull: {ipull:.4f}, chi2: {chi2_value_bin:.4f} {tc.RESET}')

            if len(running_chi2) < 1:
                running_chi2.append(chi2_value_bin)
            else:
                running_chi2.append(running_chi2[-1] + chi2_value_bin)
            
            if ikstar < 0.15 and ikstar > 0.01:
                chi2_value += chi2_value_bin
                ndf += 1
            else:
                chi2_far_from_signal += chi2_value_bin
                ndf_far_from_signal += 1
            

        print (f'{tc.BOLD+tc.RED}Chi2 value / NDF (0, 0.15) GeV/c: {chi2_value:.4f} / {ndf} {tc.RESET}')
        print (f'{tc.BOLD+tc.RED}Chi2 value / NDF (0.15, 0.4) GeV/c: {chi2_far_from_signal:.4f} / {ndf_far_from_signal} {tc.RESET}')
        print (f'{tc.BOLD+tc.RED}Running chi2: {[f"{x:.3f}" for x in running_chi2]} {tc.RESET}')
        self.pull.GetXaxis().SetRangeUser(self.LOWER_LIMIT_FIT, self.UPPER_LIMIT_FIT)
        self.pull.SetTitle(';#it{k}* (GeV/#it{c}); n_{#sigma}')
        self.pull.SetMarkerStyle(20)
        self.pull.SetMarkerSize(1.5)
        self.pull.SetMarkerColor(418)

        self.outdir.cd()
        self.pull.Write(f'pull')
        self.h_CF.Write(f'h_CF_clone')

        self.sig_counts.setVal(_sig_counts)

    def chi2_running_window(self, window_size: int= 5):
        '''
            Compute the chi2 in a running window of size `window_size` bins.
        '''

        _sig_counts = self.sig_counts.getVal()
        self.sig_counts.setVal(0.)

        _h_CF = self.h_CF.Clone('h_CF_clone')
        bin_width = _h_CF.GetBinWidth(1)
        _h_CF.GetXaxis().SetRangeUser(self.LOWER_LIMIT_FIT, self.UPPER_LIMIT_FIT)
        
        h_running_chi2 = _h_CF.Clone('pull')
        h_running_chi2.Reset()

        print(f'{tc.BOLD+tc.GREEN}\n\n>>>>>>>>>>>>>>Running chi2 with window size {window_size} bins{tc.RESET}')    
        for ibin in range(window_size + 1, _h_CF.GetNbinsX() - window_size + 1):

            ikstar_low = _h_CF.GetBinCenter(ibin - window_size)
            ikstar_high = _h_CF.GetBinCenter(ibin + window_size)

            ichi2 = 0.

            for jbin in range(ibin - window_size, ibin + window_size + 1):
                jkstar = _h_CF.GetBinCenter(jbin)
                jmeasured = _h_CF.GetBinContent(jbin)
                jerror = _h_CF.GetBinError(jbin)

                self.kstar.setVal(jkstar)
                jexpected = self.model.getVal(self.kstar) * self.bkg_counts.getVal() * bin_width
                jerr_bkg = 0
                jpull = 0. if jerror < 1.e-9 else (jmeasured - jexpected) / np.sqrt(jerror**2 + jerr_bkg**2)
                ichi2 += jpull**2

            print(f'k*: ({ikstar_low:.3f}, {ikstar_high:.3f}), chi2: {ichi2}')
            h_running_chi2.SetBinContent(ibin, ichi2)
        print()
        
        h_running_chi2.SetTitle(f'Window width = {window_size*2} bins;#it{{k}}* (GeV/#it{{c}}); #chi^{{2}}')
        h_running_chi2.SetMarkerStyle(20)
        h_running_chi2.SetMarkerSize(1.5)
        h_running_chi2.SetMarkerColor(418)

        self.outdir.cd()
        h_running_chi2.Write(f'running_chi2')

        self.sig_counts.setVal(_sig_counts)



    def prefit_CF(self) -> RooRealVar:
        '''
            Prefit the p-He3 correlation function in a region away from the signal.
            This is used to fix the background parameters.
        '''

        #self.model.fitTo(self.CF_datahist, RooFit.Save(), RooFit.Range(0.2, self.UPPER_LIMIT_FIT))
        self.model.chi2FitTo(self.CF_datahist, RooFit.Save(), RooFit.Range(0.2, self.UPPER_LIMIT_FIT), Extended=True)

        frame = self.kstar.frame()
        self.CF_datahist.plotOn(frame)
        self.model.plotOn(frame, Components={'bkg_pdf'}, LineStyle=2, LineColor=3)

        text = write_params_to_text(self.model.getParameters(self.CF_datahist), coordinates=(0.7, 0.3, 0.9, 0.5))
        frame.addObject(text)

        self.outdir.cd()
        frame.Write(f'prefit_CF')

        n_bkg_counts = self.bkg_counts.getVal()
        correction = self.h_bkg.Integral(self.h_bkg.FindBin(self.LOWER_LIMIT_FIT), self.h_bkg.FindBin(self.UPPER_LIMIT_FIT)) / \
                        self.h_bkg.Integral(self.h_bkg.FindBin(0.2), self.h_bkg.FindBin(self.UPPER_LIMIT_FIT))
        #self.bkg_counts.setVal(n_bkg_counts * correction)

    def fit_CF(self, logy: bool = False, cent: str = 'integrated', use_chi2: bool = False):
        '''
            Fit the p-He3 correlation function
        '''

        print(tc.GREEN + f'Fitting the correlation function' + tc.RESET)
        
        mean_0 = self.sig_params['sig_mean'].getVal()
        sigma_0 = self.sig_params['sig_sigma'].getVal()
        for param in self.sig_params.values():
            print(f'{param=}')
            if 'mean' in param.GetName():
                #param.setRange(mean_0 - 3 * sigma_0, mean_0 + 3 * sigma_0)
                param.setRange(param.getVal() - 3 * param.getError(), param.getVal() + 3 * param.getError())
                #param.setRange(0.6, 0.82) 
                #param.setRange(0.06, param.getVal() + 5 * param.getError())
            elif 'sigma' in param.GetName():
                #param.setRange(0.5 * param.getVal(), 3 * param.getVal())
                #param.setRange(param.getVal(), 4 * param.getVal())
                param.setRange(param.getVal() - 3 * param.getError(), param.getVal() + 3 * param.getError())
            else:
                param.setConstant(True)
            print(f'{param.GetName()} = {param.getVal()}')

        self.model = RooAddPdf('model', 'model', [self.signal_pdf, self.bkg_pdf], [self.sig_counts, self.bkg_counts])

        self.sig_counts.setVal(0.)
        self.sig_counts.setConstant(True)
        self.prefit_CF()
        self.bkg_counts.setConstant(True)
        self.sig_counts.setConstant(False)

        fit_result = None
        if use_chi2:
            fit_result = self.model.chi2FitTo(self.CF_datahist, RooFit.Save(), RooFit.Range(self.LOWER_LIMIT_FIT, self.UPPER_LIMIT_FIT), Extended=True)
        else:
            fit_result = self.model.fitTo(self.CF_datahist, RooFit.Save(), RooFit.Range(self.LOWER_LIMIT_FIT, self.UPPER_LIMIT_FIT), Extended=True, SumW2Error=True)

        print(f'{tc.BOLD+tc.RED}Fit result: {fit_result.status()} {tc.RESET}')
        print(f'{self.bkg_counts.getVal()=}, {self.bkg_counts.getError()=}')

        covariance_matrix = fit_result.covarianceMatrix()
        print(f'{tc.GREEN + "Fit result:" + tc.RESET}')
        fit_result.Print()
        print(f'{tc.GREEN + "Covariance matrix:" + tc.RESET}')
        covariance_matrix.Print()

        sig_counts, mean, sigma = fit_result.floatParsFinal()
        self.sig_counts.setVal(sig_counts.getVal())
        self.sig_counts.setError(sig_counts.getError())

        frame = self.kstar.frame(Title=';#it{k}* (GeV/#it{c});C(#it{k}*)')
        self.CF_datahist.plotOn(frame, MarkerColor=1, MarkerStyle=20, MarkerSize=1.5)
        self.model.plotOn(frame, LineColor=2)
        self.model.plotOn(frame, Components={'signal_pdf'}, LineStyle=2, LineColor=4)
        self.model.plotOn(frame, Components={'bkg_pdf'}, LineStyle=2, LineColor=3)

        print(f'{tc.BOLD+tc.RED}bkg{tc.RESET}{self.bkg_counts.getVal()=}, {self.bkg_counts.getError()=}')

        watermark = get_alice_watermark(0.5, 0.3, 0.8, 0.45)
        centrality_watermark = get_centrality_watermark(0.5, 0.2, 0.8, 0.3, cent)
        sign_text = '#bf{#bar{p}-^{3}#bar{He}}' if self.sign == 'Anti' else '#bf{p-^{3}He}'
        centrality_watermark.AddText(f'{sign_text}')
        text = write_params_to_text({key: self.sig_params[key] for key in ['sig_mean', 'sig_sigma']}.values(),
                                     coordinates=(0.5, 0.03, 0.8, 0.19))
        text.AddText(f'#it{{n}}_{{sig}} = ({self.sig_counts.getVal():.3f} #pm {self.sig_counts.getError():.3f})')
        #text.AddText(f'bkg counts = ({self.bkg_counts.getVal():.4f} #pm {self.bkg_counts.getError():.4f})')
        text.AddText(f'#chi^{{2}} / NDF = {frame.chiSquare():.2f}')
        frame.addObject(watermark)
        frame.addObject(centrality_watermark)
        frame.addObject(text)
        
        print(f'{tc.BOLD+tc.RED}bkg{tc.RESET}{self.bkg_counts.getVal()=}, {self.bkg_counts.getError()=}')

        self.draw_cf_canvas(frame, logy)
        self.outdir.cd()
        frame.Write(f'fit_CF')

        print(f'{tc.BOLD+tc.RED}bkg{tc.RESET}{self.bkg_counts.getVal()=}, {self.bkg_counts.getError()=}')

    def compute_chi2_significance(self):

        self.model.chi2FitTo(self.CF_datahist, RooFit.Save(), RooFit.Range(0.01, self.UPPER_LIMIT_FIT), Extended=True)
        frame = self.kstar.frame(Title=';#it{k}* (GeV/#it{c});C(#it{k}*)')
        self.CF_datahist.plotOn(frame, MarkerColor=1, MarkerStyle=20, MarkerSize=1.5)
        self.model.plotOn(frame, LineColor=2)
        chi2_model = frame.chiSquare()

        nsig_value = self.sig_counts.getVal()
        self.sig_counts.setVal(0.)
        self.sig_counts.setConstant(True)
        
        self.model.chi2FitTo(self.CF_datahist, RooFit.Save(), RooFit.Range(0.01, self.UPPER_LIMIT_FIT), Extended=True)  # fit without signal
        frame_no_signal = self.kstar.frame(Title=';#it{k}* (GeV/#it{c});C(#it{k}*)')
        self.CF_datahist.plotOn(frame_no_signal, MarkerColor=1, MarkerStyle=20, MarkerSize=1.5)
        self.model.plotOn(frame_no_signal, LineColor=2)
        chi2_bkg = frame_no_signal.chiSquare()
        
        self.sig_counts.setVal(nsig_value)
        self.sig_counts.setConstant(False)

        print(f'{tc.BOLD+tc.RED}Chi2 model: {chi2_model:.4f} {tc.RESET}')
        print(f'{tc.BOLD+tc.RED}Chi2 bkg: {chi2_bkg:.4f} {tc.RESET}')

    def integral_CF(self):
        '''
            Compute the integral of the k*^2 C(k*)
            This is used to compute the yield of Li4
        '''

        crystal_ball = TF1('cb', CrystalBall, self.LOWER_LIMIT_FIT, self.UPPER_LIMIT_FIT, 6)
        for iparam, param in enumerate(self.sig_params.values()):
                crystal_ball.SetParameter(iparam, param.getVal())
        def integrand_func(x, par):
            norm = par[0]
            return norm * x[0] * x[0] * crystal_ball.Eval(x[0])
        integrand = TF1('CF_integral', integrand_func, self.LOWER_LIMIT_FIT, 0.25, 1)
        integrand.SetParameter(0, self.sig_counts.getVal())

        integral = integrand.Integral(self.LOWER_LIMIT_FIT, self.UPPER_LIMIT_FIT)
        canvas = TCanvas('integral_canvas', 'integral_canvas', 800, 600)
        canvas.DrawFrame(0, 0, self.UPPER_LIMIT_FIT, 0.007, ';#it{k}* (GeV/#it{c}); #it{k}*^{{2}} C(#it{k}*)')
        integrand.Draw('same')
        text = TPaveText(0.5, 0.6, 0.8, 0.8, 'NDC')
        text.SetFillColor(0)
        text.SetBorderSize(0)
        text.AddText(f'#int_{{0}}^{{{self.UPPER_LIMIT_FIT} GeV/#it{{c}}}} #it{{k}}*^{{2}} C_{{R}}(#it{{k}}*) d#it{{k}}* = {integral:.8f}')
        text.Draw('same')

        self.outdir.cd()
        crystal_ball.Write('crystal_ball')
        integrand.Write('integrand')
        canvas.Write('integral')

    def extract_signal_counts(self, h_mixed_event_centrality: list, normalisation=1., signal_normalisation=None):
        '''
            norm is the normalisation for the total mixed event histogram. This is not necessary if the histograms are already normalised.
        '''
        h_mixed_event = h_mixed_event_centrality[0].Clone('h_mixed_event')
        for ihist in h_mixed_event_centrality[1:]:
            h_mixed_event.Add(ihist)

        h_mixed_event.Scale(normalisation)
        
        #_bkg_counts = self.bkg_counts.getVal()
        #self.bkg_counts.setVal(0)
        #norm = self.model.createIntegral(self.kstar, NormSet={self.kstar}).getVal()

        norm_cf = signal_normalisation if signal_normalisation is not None else self.sig_counts.getVal()
        CF_signal = TF1('cb', CrystalBall, 0.0, self.UPPER_LIMIT_FIT, 6)
        for iparam, param in enumerate(self.sig_params.values()):
                CF_signal.SetParameter(iparam, param.getVal())
                if iparam == 0 or iparam == 1:
                    CF_signal.SetParError(iparam, param.getError())

        def normalised_func(x, par):
            norm = par[0]
            return norm * CF_signal.Eval(x[0])
        normalised_CF_signal = TF1('CF_integral', normalised_func, self.LOWER_LIMIT_FIT, self.UPPER_LIMIT_FIT, 1)
        normalised_CF_signal.SetParameter(0, norm_cf)

        flag = '_upper_limit' if signal_normalisation is not None else ''
        integral_error_signal = norm_cf * CF_signal.IntegralError(self.LOWER_LIMIT_FIT, self.UPPER_LIMIT_FIT)
        print(tc.GREEN + f'Errors from the signal pdf: {flag}' + tc.RESET)
        print(f'{tc.BOLD+tc.RED}Signal pdf error = {integral_error_signal:.4f} {tc.RESET}')

        print(tc.GREEN + f'Extracting signal counts from the mixed event histogram: {flag}' + tc.RESET)
        print(f'{tc.BOLD+tc.RED}Signal counts = {self.sig_counts.getVal():.4f} +/- {self.sig_counts.getError():.4f} {tc.RESET}')
        print(f'{tc.BOLD+tc.RED}Signal normalisation = {norm_cf} {tc.RESET}')

        #self.bkg_counts.setVal(_bkg_counts)

        h_same_event = h_mixed_event.Clone('h_same_event')
        for ibin in range(1, h_same_event.GetNbinsX() + 1):
            kstar = h_same_event.GetBinCenter(ibin)
            bin_content = h_mixed_event.GetBinContent(ibin)
            bin_error = h_mixed_event.GetBinError(ibin)
            signal_weigth = normalised_CF_signal.Eval(kstar)
            h_same_event.SetBinContent(ibin, bin_content * signal_weigth)
            h_same_event.SetBinError(ibin, bin_error * signal_weigth)

        self.outdir.cd()
        h_mixed_event.Write(f'h_mixed_event_norm{normalisation}{flag}')
        normalised_CF_signal.Write(f'CF_signal_norm{normalisation}{flag}')
        h_same_event.Write(f'h_same_event_norm{normalisation}{flag}')

        print(tc.GREEN + f'Extracted signal counts from the mixed event histogram' + tc.RESET)
        print(f'{tc.BOLD+tc.RED}Raw yield {flag}: {h_same_event.Integral():4f} {tc.RESET}')

    def setup_workspace(self):
        '''
            Prepare the model for significance estimation
        '''
        
        if self.roows:
            del self.roows
        self.roows = RooWorkspace('ws')

        getattr(self.roows, 'import')(self.kstar)
        getattr(self.roows, 'import')(self.CF_datahist)

        model_config = RooStats.ModelConfig('model_config', self.roows)
        model_config.SetPdf(self.model)
        model_config.SetParametersOfInterest({self.sig_counts})
        model_config.SetObservables({self.kstar})
        getattr(self.roows, 'import')(model_config)
    
    def get_bkg_config(self):

        '''
            Get the bkg config
        '''

        self.roows.Print()
        model_config = self.roows.obj('model_config')
        poi = model_config.GetParametersOfInterest().first()
        model_config.SetSnapshot({poi})

        bkg_config = model_config.Clone()
        bkg_config.SetName('bkg_config')

        old_poi_val = poi.getVal()
        poi.setVal(0)
        bkg_config.SetSnapshot(poi)
        poi.setVal(old_poi_val)

        bkg_config.Print()

        return bkg_config

    def evaluate_pvalue_and_significance(self, cent: str):

        self.setup_workspace()
        model_config = self.roows.obj('model_config')
        model_config.Print()
        datahist = self.roows.data('datahist')
        bkg_config = self.get_bkg_config()
        self.roows.var('sig_mean').setConstant(True)
        self.roows.var('sig_sigma').setConstant(True)

        print(f'Computing p-value and significance for centrality {cent}')
        mean = self.sig_params['sig_mean'].getVal()
        sigma = self.sig_params['sig_sigma'].getVal()
        #low = mean - 3 * sigma
        low = self.kstar.getMin()
        #high = mean + 3 * sigma
        high = self.kstar.getMax()
        self.kstar.setRange('roi', low, high)
        data_roi = self.CF_datahist.reduce(f'kstar >= {low} && kstar <= {high}')

        print(f'{tc.BOLD+tc.RED}mean/sigma: {mean}/{sigma} {tc.RESET}')
        print(f'{tc.BOLD+tc.RED}kstar range: {low} - {high} {tc.RESET}')

        # chi2 as pseudo-NLL
        chi2 = RooChi2Var("chi2", "chi2", self.model, data_roi, True, RooAbsData.Expected)
        pseudo_pdf = RooGenericPdf("pseudo_pdf", "exp(-0.5 * chi2)", [chi2])
        model_config.SetPdf(pseudo_pdf)

        asymp_calc = RooStats.AsymptoticCalculator(data_roi, model_config, bkg_config)
        #asymp_calc.SetOneSided(True)
        asymp_calc.SetPrintLevel(0)
        asymp_calc_result = asymp_calc.GetHypoTest()
        p_value = asymp_calc_result.NullPValue()
        p_value_err = asymp_calc_result.NullPValueError()
        significance = asymp_calc_result.Significance()
        significance_error = asymp_calc_result.SignificanceError()

        _sig_counts = self.sig_counts.getVal()
        self.sig_counts.setConstant(False)
        test_inverter = RooStats.HypoTestInverter(asymp_calc)
        #test_inverter.SetConfidenceLevel(0.95)
        test_inverter.SetConfidenceLevel(0.9)
        test_inverter.UseCLs(True)

        result = test_inverter.GetInterval()
        upper_limit = result.UpperLimit()
        upper_limit_error = result.UpperLimitEstimatedError()
        
        self.sig_counts.setVal(_sig_counts)
        self.sig_counts.setConstant(True)
        
        
        print(tc.BOLD+tc.WHITE+f'\n-----------------------------------------'+tc.RESET)
        print(tc.RED+f'Fit significance'+tc.RESET)
        print(tc.GREEN+f'centrality = {cent}'+tc.RESET)
        print(f'{tc.BOLD+tc.RED}p-value = {p_value:.10f} +/- {p_value_err:.10f} {tc.RESET}')
        print(f'{tc.BOLD+tc.RED}significance = {significance:.10f} +/- {significance_error:.10f} {tc.RESET}')
        print(f'{tc.BOLD+tc.RED}95% confidence interval  upper limit on signal counts: [{upper_limit:.10f} +/- {upper_limit_error:.10f}] {tc.RESET}')
        print(tc.BOLD+tc.WHITE+f'\n-----------------------------------------'+tc.RESET)
        input()

    def evaluate_confidence_levels(self, cent: str) -> float:
        
        self.setup_workspace()
        model_config = self.roows['model_config']
        datahist = self.roows['datahist']

        # ---------------------------------------------
        # create and use the ProfileLikelihoodCalculator
        # to find and plot the 95% confidence interval
        # on the parameter of interest as specified
        # in the model config

        CONFIDENCE_LEVEL = 0.95
        nScanPoints = 50
        poiXMin = 1
        poiXMax = 0

        pl = RooStats.ProfileLikelihoodCalculator(datahist, model_config)
        pl.SetConfidenceLevel(CONFIDENCE_LEVEL)
        interval = pl.GetInterval()

        # print out the interval on the first Parameter of Interest
        firstPOI = model_config.GetParametersOfInterest().first()
        limit_lower, limit_upper = interval.LowerLimit(firstPOI), interval.UpperLimit(firstPOI)
        #limit_upper_error = interval.UpperLimitEstimatedError(firstPOI)
        print(f"\n>>>> RESULT : {CONFIDENCE_LEVEL * 100}% interval on {firstPOI.GetName()} is : [{limit_lower}, {limit_upper}]\n")
        #print(f'{tc.BOLD+tc.RED}95% confidence interval  upper limit on signal counts: [{limit_upper:.10f} +/- {limit_upper_error:.10f}] {tc.RESET}')

        # make a plot

        print(
            "making a plot of the profile likelihood function ....(if it is taking a lot of time use less points or the "
            "TF1 drawing option)\n"
        )
        plot = RooStats.LikelihoodIntervalPlot(interval)
        plot.SetNPoints(nScanPoints)  # do not use too many points, it could become very slow for some models
        if poiXMin < poiXMax:
            plot.SetRange(poiXMin, poiXMax)
        opt = "tf1"

        self.outdir.cd()

        canvas = TCanvas('confidence_interval', 'confidence_interval', 800, 600)
        plot.SetTitle(';#it{n}_{sig}; -2 log #lambda')
        plot.Draw(opt)
        canvas.Write(f'confidence_interval_{cent}')

        return limit_upper 

    def evaluate_pvalue_and_significance_counts(self, cent: str):
        '''
            Use the expression valid for counts: nsig/ sqrt(nsig + nbkg) = nsig / sqrt(ntot)
        '''

        self.setup_workspace()
        model_config = self.roows.obj('model_config')
        model_config.Print()
        datahist = self.roows.data('datahist')
        bkg_config = self.get_bkg_config()
        self.roows.var('sig_mean').setConstant(True)
        self.roows.var('sig_sigma').setConstant(True)

        print(f'Computing p-value and significance for centrality {cent}')
        mean = self.sig_params['sig_mean'].getVal()
        sigma = self.sig_params['sig_sigma'].getVal()
        low = mean - 4 * sigma
        high = mean + 4 * sigma
        self.kstar.setRange('roi', low, high)

        ntot = self.model.createIntegral(self.kstar, NormSet={self.kstar}, Range='roi')
        nsig =  self.signal_pdf.createIntegral(self.kstar, NormSet={self.kstar}, Range='roi')

        ntot_check = self.model.createIntegral(self.kstar, NormSet={self.kstar})
        norm_tot = self.sig_counts.getVal()+self.bkg_counts.getVal()
        nsig_check = self.signal_pdf.createIntegral(self.kstar, NormSet={self.kstar})
        norm_sig = self.sig_counts.getVal()

        print(f'{tc.BOLD+tc.RED} {nsig_check.getVal()=} {tc.RESET}')
        print(f'{tc.BOLD+tc.RED} {self.sig_counts.getVal()=} {tc.RESET}')
        print(f'{tc.BOLD+tc.RED} {ntot_check.getVal()=} {tc.RESET}')
        print(f'{tc.BOLD+tc.RED} {self.sig_counts.getVal()+self.bkg_counts.getVal()=} {tc.RESET}')

        significance = nsig.getVal() * norm_sig / np.sqrt(ntot.getVal() * norm_tot)
        #significance_error = significance * np.sqrt((nsig.getError()/nsig.getVal())**2 + (ntot.getError()/ntot.getVal())**2)

        print(tc.BOLD+tc.WHITE+f'\n-----------------------------------------'+tc.RESET)
        print(tc.RED+f'Fit significance (counts)' + tc.RESET)
        print(tc.GREEN+f'centrality = {cent}'+tc.RESET)
        print(f'nsig = {nsig.getVal() * norm_sig:.5f}') # +/- {nsig.getError():.10f}')
        print(f'ntot = {ntot.getVal() * norm_tot:.5f}') # +/- {ntot.getError():.10f}')
        print(f'significance = {significance:.5f}') # +/- {significance_error:.10f}')
        print(tc.BOLD+tc.WHITE+f'\n-----------------------------------------'+tc.RESET)
        input()
    
    def scan_pvalue_and_significance(self):
        '''
            Perform a scan in kstar and compute the significance
        '''        

        self.setup_workspace()
        model_config = self.roows.obj('model_config')
        datahist = self.roows.data('datahist')
        bkg_config = self.get_bkg_config()
        self.roows.var('sig_mean').setConstant(True)
        self.roows.var('sig_sigma').setConstant(True)
        kstars, p0_values = [], []

        kstar_array = np.linspace(self.kstar.getMin(), self.kstar.getMax(), 40)
        for ikstar in kstar_array:

            self.roows.var('sig_mean').setVal(ikstar)
            self.roows.var('sig_mean').setConstant(True)

            asymp_calc_scan = RooStats.AsymptoticCalculator(datahist, model_config, bkg_config)
            asymp_calc_scan.SetPrintLevel(-1)
            asymp_calc_scan.SetOneSidedDiscovery(True)
            asym_calc_result_scan = asymp_calc_scan.GetHypoTest()
            kstars.append(ikstar)
            p0_values.append(asym_calc_result_scan.NullPValue())

        g_local_pvalue = TGraph(len(kstars), np.array(kstars), np.array(p0_values))
        g_local_pvalue.SetName('p0_values')
        g_local_pvalue.SetTitle('; #it{k}* (GeV/#it{c}); Local p-value')
        
        g_local_pvalue.SetMarkerStyle(20)
        g_local_pvalue.SetMarkerColor(418)

        self.outdir.cd()
        g_local_pvalue.Write()

    def clear_data(self):
        '''
            Clear the data
        '''

        gStyle.SetOptStat(0)
        
        self.outdir = None

        self.h_CF = None
        self.CF_datahist = None
        self.CF_datahist_significance = None
        self.pull = None

        # Roofit
        self.kstar = RooRealVar('kstar', '#it{k}* (GeV/#it{c})', self.LOWER_LIMIT_FIT, self.UPPER_LIMIT_FIT)
        self.kstar_significance = RooRealVar('kstar_significance', 'kstar_significance', self.LOWER_LIMIT_FIT, self.UPPER_LIMIT_FIT)
        
        self.sig_counts = RooRealVar('sig_counts', 'sig_counts', 0.05, 0., 1e3)
        self.signal_pdf = None
        self.sig_params = {}
        self.bkg_counts = RooRealVar('bkg_counts', 'bkg_counts', 0., 1e6)
        self.bkg_datahist = None
        self.bkg_pdf = None
        self.model = None

        self.roows = None # used for the asymptotic calculator

    def close_canvas(self):
        '''
            Close the canvas
        '''

        self.pdf_canvas.Clear()
        self.pdf_canvas.Print(f'{self.pdf_outpath})')
