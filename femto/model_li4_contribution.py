import numpy as np
from ROOT import TFile, TCanvas, TH1F, TLegend, TH2F, TTree, \
                kOrange, kBlue, kGreen, \
                RooDataSet, RooRealVar, RooDataHist, RooFit, RooFFTConvPdf, RooKeysPdf

from torchic.roopdf import RooSillPdf, RooSillGeneralizedPdf
from torchic.core.histogram import load_hist

from core.roopdf_utils import init_roopdf

class Sampler:

    def __init__(self, mass, intrinsic_width, experimental_width, eth, outfile:TFile, var = None):
        self._mass = mass
        self._intrinsic_width = intrinsic_width
        self._experimental_width = experimental_width
        self._eth = eth
        self._l = 2

        self._var = RooRealVar('invmass', 'm (p+^{3}He)', 3.737, 3.865, 'GeV/#it{c}^{2}')  if var is None else var
        self._sampled_roo_dataset = None
        self._pdf = None
        self._pdf_pars = {}

        self._shape_pdf = None      # pdf to get the shape from the MC (intrinsic resolution is missing)
        self._shape_pdf_pars = {}

        self._outfile = outfile
    
    @property
    def sampled_roo_dataset(self):
        return self._sampled_roo_dataset
    
    def init_shape_from_mc(self, h_mc_signal:TH1F, shape_hypothesis: str):

        mean = RooRealVar('experimental_mean', '#mu', 3.751, 3.746, 3.756, 'GeV/c^{2}')
        sigma = RooRealVar('experimental_sigma', '#sigma', 0.01, 0.0001, 0.1, 'GeV/c^{2}')
        self._shape_pdf, self._shape_pdf_pars = init_roopdf(shape_hypothesis, self._var, mean=mean, sigma=sigma, name='experimental')

        mc_data_hist = RooDataHist('mc_data_hist', 'mc_data_hist', [self._var], Import=h_mc_signal)
        self._shape_pdf.fitTo(mc_data_hist, RooFit.Save(), Range=(3.750, 3.753))

        frame = self._var.frame()
        mc_data_hist.plotOn(frame)
        self._shape_pdf.plotOn(frame, LineColor=kBlue+2)
        self._shape_pdf.paramOn(frame)
        self._var.setRange(3.737, 3.865)
        canvas = TCanvas('c_experimental_points')
        frame.Draw()
        self._outfile.cd()
        canvas.Write()

    def sample(self, shape: str, n_samples: int = 1_000_000):
        if shape == 'sill':                             self._sample_sill(n_samples)
        if shape == 'gaus_conv':                        self._sample_gaus_conv(n_samples)
        if shape == 'sill_gaus_conv':                   self._sample_sill_gaus_conv(n_samples)
        if shape == 'sill_gaus_conv_numpy':             self._sample_sill_gaus_conv_numpy(n_samples)
        if shape == 'sill_generalized_gaus_conv_numpy': self._sample_sill_generalized_gaus_conv_numpy(n_samples)

    def _sample_sill(self, n_samples: int):
        self._pdf_pars = {
            'mass': RooRealVar('mass', 'mass', self._mass),
            'gamma': RooRealVar('gamma', 'gamma', self._intrinsic_width),
            'eth': RooRealVar('eth', 'eth', self._eth)
        }
        self._pdf = RooSillPdf('sill', 'sill', self._var, *self._pdf_pars.values())
        self._sampled_roo_dataset = self._pdf.generate([self._var], NumEvents=n_samples)

    def _sample_gaus_conv(self, n_samples: int):
        if self._shape_pdf is None:
            raise ValueError('You must first initialise the shape with "init_shape_from_mc".')

        mean = RooRealVar('intrinsic_mean', '#mu', 3.751, 3.746, 3.756, 'GeV/c')
        mean.setVal(self._shape_pdf_pars['mean'].getVal())
        sigma = RooRealVar('intrinsic_sigma', '#sigma', 0.003, 0.0001, 0.01, 'GeV/c')  # Add range
        pdf_intrinsic_resolution, pdf_pars_intrinsic_resolution = init_roopdf('gaus', self._var, mean=mean, sigma=sigma, name='intrinsic')

        frame = self._var.frame()
        self._shape_pdf.plotOn(frame, LineColor=kOrange-3)
        pdf_intrinsic_resolution.plotOn(frame, LineColor=kGreen+1)
        canvas = TCanvas('c_gaus_conv')
        frame.Draw()
        self._outfile.cd()
        canvas.Write()

        self._var.setRange(3.73, 3.87)  # Wider range for FFT
        self._var.setBins(4000, "cache")  # More bins for numerical stability

        mean.setConstant(True)
        sigma.setConstant(True)

        self._pdf = RooFFTConvPdf('convolution', 'convolution', self._var, self._shape_pdf, pdf_intrinsic_resolution)

        frame_conv = self._var.frame()
        self._pdf.plotOn(frame_conv, LineColor=kBlue+1)
        canvas_conv = TCanvas('c_convolution_pdf')
        frame_conv.Draw()
        self._outfile.cd()
        canvas_conv.Write()

        self._sampled_roo_dataset = self._pdf.generate([self._var], NumEvents=n_samples)
    
    def _sample_sill_gaus_conv(self, n_samples: int):
        if self._shape_pdf is None:
            raise ValueError('You must first initialise the shape with "init_shape_from_mc".')

        pdf_pars_intrinsic_resolution = {
            'mass': RooRealVar('intrinsic_mass', 'mass', self._mass),
            'gamma': RooRealVar('intrinsic_gamma', 'gamma', self._intrinsic_width),
            'eth': RooRealVar('intrinsic_eth', 'eth', self._eth)
        }
        pdf_intrinsic_resolution = RooSillPdf('sill', 'sill', self._var, *pdf_pars_intrinsic_resolution.values())

        frame = self._var.frame()
        self._shape_pdf.plotOn(frame, LineColor=kOrange-3)
        pdf_intrinsic_resolution.plotOn(frame, LineColor=kGreen+1)
        canvas = TCanvas('c_gaus_conv')
        frame.Draw()
        self._outfile.cd()
        canvas.Write()

        self._var.setRange(3.73, 3.87)  # Wider range for FFT
        self._var.setBins(4000, "cache")  # More bins for numerical stability

        pdf_pars_intrinsic_resolution['mean'].setConstant(True)
        pdf_pars_intrinsic_resolution['gamma'].setConstant(True)

        self._pdf = RooFFTConvPdf('convolution', 'convolution', self._var, self._shape_pdf, pdf_intrinsic_resolution)

        frame_conv = self._var.frame()
        self._pdf.plotOn(frame_conv, LineColor=kBlue+1)
        canvas_conv = TCanvas('c_convolution_pdf')
        frame_conv.Draw()
        self._outfile.cd()
        canvas_conv.Write()

        self._sampled_roo_dataset = self._pdf.generate([self._var], NumEvents=n_samples)

    def _sample_sill_gaus_conv_numpy(self, n_samples: int):

        if self._shape_pdf is None:
            raise ValueError('You must first initialise the shape with "init_shape_from_mc".')

        frame = self._var.frame()
        self._shape_pdf.plotOn(frame, LineColor=kBlue+2)
        self._var.setRange(3.737, 3.865)
        canvas = TCanvas('c_experimental_shape')
        frame.Draw()
        self._outfile.cd()
        canvas.Write()
        del frame

        sill_params = {
            'mass': RooRealVar('mass', 'mass', self._mass),
            'gamma': RooRealVar('gamma', 'gamma', self._intrinsic_width),
            'eth': RooRealVar('eth', 'eth', self._eth)
        }
        sill_pdf = RooSillPdf('sill', 'sill', self._var, *sill_params.values())

        frame = self._var.frame()
        sill_pdf.plotOn(frame, LineColor=kOrange-3)
        canvas = TCanvas('c_intrinsic_shape')
        frame.Draw()
        self._outfile.cd()
        canvas.Write()

        sill_dataset = sill_pdf.generate([self._var], NumEvents=n_samples)
        experimental_sigma = self._shape_pdf_pars['sigma'].getVal()
        self._sampled_roo_dataset = RooDataSet('smeared_data', 'smeared_data', [self._var])

        for ientry in range(sill_dataset.numEntries()):
            
            invmass_true = sill_dataset.get(ientry)['invmass'].getVal()
            invmass_smeared = np.random.normal(invmass_true, experimental_sigma)
            
            if self._var.getMin() <= invmass_smeared <= self._var.getMax():
                self._var.setVal(invmass_smeared)
                self._sampled_roo_dataset.add([self._var])


        self._pdf = sill_pdf
        self._pdf_pars = sill_params

    def _sample_sill_generalized_gaus_conv_numpy(self, n_samples: int):

        if self._shape_pdf is None:
            raise ValueError('You must first initialise the shape with "init_shape_from_mc".')

        frame = self._var.frame()
        self._shape_pdf.plotOn(frame, LineColor=kBlue+2)
        self._var.setRange(3.737, 3.865)
        canvas = TCanvas('c_experimental_shape')
        frame.Draw()
        self._outfile.cd()
        canvas.Write()
        del frame

        sill_params = {
            'mass': RooRealVar('mass', 'mass', self._mass),
            'gamma': RooRealVar('gamma', 'gamma', self._intrinsic_width),
            'eth': RooRealVar('eth', 'eth', self._eth),
            'l': RooRealVar('l', 'l', self._l)
        }
        sill_pdf = RooSillGeneralizedPdf('sill', 'sill', self._var, *sill_params.values())

        frame = self._var.frame()
        sill_pdf.plotOn(frame, LineColor=kOrange-3)
        canvas = TCanvas('c_intrinsic_shape')
        frame.Draw()
        self._outfile.cd()
        canvas.Write()

        sill_dataset = sill_pdf.generate([self._var], NumEvents=n_samples)
        experimental_sigma = self._shape_pdf_pars['sigma'].getVal()
        self._sampled_roo_dataset = RooDataSet('smeared_data', 'smeared_data', [self._var])

        for ientry in range(sill_dataset.numEntries()):
            
            invmass_true = sill_dataset.get(ientry)['invmass'].getVal()
            invmass_smeared = np.random.normal(invmass_true, experimental_sigma)
            
            if self._var.getMin() <= invmass_smeared <= self._var.getMax():
                self._var.setVal(invmass_smeared)
                self._sampled_roo_dataset.add([self._var])


        self._pdf = sill_pdf
        self._pdf_pars = sill_params

    def save_sampling(self):
        
        frame = self._var.frame()
        self._sampled_roo_dataset.plotOn(frame, MarkerColor=kGreen+2, LineColor=kGreen+2, MarkerStyle=20, Name='Convoluted')
        self._pdf.plotOn(frame, LineColor=kOrange-3, Name='Theoretical')#, Normalization=self._sampled_roo_dataset.numEntries())

        legend = TLegend(0.6, 0.42, 0.85, 0.64)
        legend.SetBorderSize(0)
        legend.SetTextSize(0.04)
        legend.AddEntry(frame.findObject('Theoretical'), 'Theoretical', 'l')
        legend.AddEntry(frame.findObject('Convoluted'), 'Convoluted', 'l')
        
        canvas = TCanvas('cSill', '', 800, 600)
        frame.Draw()
        legend.Draw('same')
        self._outfile.cd()
        canvas.Write()


def convert_invmass_to_kstar(invmass: float):
    
    M_PR = 0.938272 # GeV/c^2
    M_HE = 2.80839 # GeV/c^2
    lambda_func = lambda x, y, z: x**2 + y**2 + z**2 -2*x*y - 2*x*z - 2*y*z
    numerator = np.sqrt(lambda_func(invmass**2, M_PR**2, M_HE**2)) if lambda_func(invmass**2, M_PR**2, M_HE**2) >= 0 else 0
    denominator = 2 * invmass
    return numerator / denominator    

def sampling(outfile: TFile, n_samples: int = 1_000_000) -> RooDataSet:

    invmass = RooRealVar('invmass', 'm (p+^{3}He)', 3.743, 3.755, 'GeV/#it{c}^{2}')
    sill_params = {
        'mass': RooRealVar('mass', 'mass', M_LI4),
        'gamma': RooRealVar('mass', 'mass', W_LI4),
        'eth': RooRealVar('eth', 'eth', E_TH)
    }
    sill_pdf = RooSillPdf('sill', 'sill', invmass, *sill_params.values())
    roo_dataset = sill_pdf.generate([invmass], NumEvents=n_samples)

    frame = invmass.frame()
    sill_pdf.plotOn(frame, LineColor=kOrange-3)

    canvas = TCanvas('cSill', '', 800, 600)
    frame.Draw()
    outfile.cd()
    canvas.Write()

    return roo_dataset

def draw_kstar_profile(roo_dataset: RooDataSet, outfile: TFile) -> TH1F:

    h_kstar = TH1F('hKstar', ';#it{k}* (GeV/#it{c});Counts', 400, 0, 0.4)
    h2_kstar_vs_invmass = TH2F('h2KstarVsInvmass', ';m (GeV/#it{c}^{2});#it{k}* (GeV/#it{c})', 800, 3.737, 3.865, 400, 0, 0.4)
    for ientry in range(roo_dataset.numEntries()):
        invmass = roo_dataset.get(ientry)['invmass'].getVal()
        kstar = convert_invmass_to_kstar(invmass)
        h_kstar.Fill(kstar)
        h2_kstar_vs_invmass.Fill(invmass, kstar)

    h_kstar.SetBinContent(1, 0.)    # in the conversion, this bin behaves as an underflow bin

    outfile.cd()
    h_kstar.Write()
    h2_kstar_vs_invmass.Write()

    return h_kstar

def draw_kstar_profile_from_hist(roo_dataset: RooDataSet, outfile: TFile) -> TH1F:

    roo_datahist = roo_dataset.binnedClone()

    h_kstar = TH1F('hKstarHist', ';#it{k}* (GeV/#it{c});Counts', 400, 0, 0.4)
    for ientry in range(1, roo_datahist.getBinnings()[0].numBins()):
        invmass = roo_datahist.get(ientry)['invmass'].getVal()
        weight = roo_datahist.weight(ientry)
        kstar = convert_invmass_to_kstar(invmass)
        h_kstar.Fill(kstar, weight)
    for ibin in range(1, h_kstar.GetNbinsX()+1):
        bin_value = h_kstar.GetBinContent(ibin)
        h_kstar.SetBinError(ibin, np.sqrt(bin_value))

    h_cumualative = h_kstar.GetCumulative()

    outfile.cd()
    h_kstar.Write()
    h_cumualative.Write()

    return h_kstar

def match_bin_width_correlation_function(h_source, h_target):
    """
    Adjust the bin widths of the target histogram to match those of the source histogram.
    
    Args:
        h_source: Histogram with desired bin widths
        h_target: Histogram to be adjusted
    """
    
    h_target_matched = h_source.Clone(f'{h_target.GetName()}_matched')
    for ibin in range(1, h_target_matched.GetNbinsX()+1):
        
        if h_target_matched.GetBinCenter(ibin) > 0.4: 
            continue

        kstar_low = h_target_matched.GetBinLowEdge(ibin)
        kstar_high = h_target_matched.GetBinLowEdge(ibin+1)
        bin_sum, bin_count = 0., 0.

        for jbin in range(1, h_target.GetNbinsX()+1):
            bin_center = h_target.GetBinCenter(jbin)
            if kstar_low <= bin_center < kstar_high:
                bin_sum += h_target.GetBinContent(jbin)
                bin_count += 1
        if bin_count > 0:
            h_target_matched.SetBinContent(ibin, bin_sum / bin_count)
    return h_target_matched

def smoothen_histogram(hist, outfile, xmin:float=0.0, xmax:float=0.42, rho:float=2.):

    kstar = RooRealVar('kstar', 'k*', xmin, xmax)
    
    x_data = []
    weights = []
    
    for ibin in range(1, hist.GetNbinsX()+1):
        kstar_val = hist.GetBinCenter(ibin)
        if not (xmin <= kstar_val <= xmax):
            continue
        y_val = hist.GetBinContent(ibin)
        if y_val <= 0:
            continue
        
        x_data.append(kstar_val)
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
    dataset = RooDataSet(hist.GetName()+'_roodata', hist.GetName()+'_roodata', 
                         tree, [kstar, weight_var], '', 'weight')
    
    keys_pdf = RooKeysPdf(f'{hist.GetName()}_keys', f'{hist.GetName()}_keys', 
                          kstar, dataset, RooKeysPdf.NoMirror, rho)
    

    frame = kstar.frame(len(x_data))
    dataset.plotOn(frame, MarkerStyle=20, MarkerSize=0.8, LineColor=1)
    keys_pdf.plotOn(frame, LineColor=2, LineWidth=2)
    canvas = TCanvas(f'cKeysPdf_{hist.GetName()}', f'cKeysPdf_{hist.GetName()}', 800, 600)
    frame.Draw()

    outfile.cd()
    hist.Write(f'{hist.GetName()}_original')
    keys_pdf.Write(f'{hist.GetName()}_keyspdf')
    canvas.Write()

    frame_cumulative = frame.emptyClone(f'{hist.GetName()}_cumulative')
    cumulative_pdf = keys_pdf.createCdf([kstar], Name=f'{hist.GetName()}_cumulative_pdf')
    cumulative_pdf.plotOn(frame_cumulative, LineColor=2, LineWidth=2)
    canvas_cumulative = TCanvas(f'cKeysPdfCumulative_{hist.GetName()}', f'cKeysPdfCumulative_{hist.GetName()}', 800, 600)
    frame_cumulative.Draw()
    outfile.cd()
    canvas_cumulative.Write()
    
    return keys_pdf

def draw_ck_profile_from_hist(h_signal:TH1F, h_mixed_event:TH1F, outfile:TFile, 
                              h_correlation_function_reference:TH1F=None) -> None:

    h_ck = h_mixed_event.Clone('hCkHist')
    h_ck.SetTitle(';#it{k}* (GeV/#it{c});#it{C}_{^{4}Li}(#it{k}*)')

    for ibin in range(1, h_mixed_event.GetNbinsX()+1):

        bin_center = h_ck.GetBinCenter(ibin)
        value = h_signal.GetBinContent(h_signal.FindBin(bin_center)) / h_mixed_event.GetBinContent(ibin) \
                    if bin_center < 0.4 else 0.
        error = h_signal.GetBinError(h_signal.FindBin(bin_center)) / h_mixed_event.GetBinContent(ibin) \
                    if bin_center < 0.4 else 0.
        h_ck.SetBinContent(ibin, value)
        h_ck.SetBinError(ibin, error)
    
    if h_correlation_function_reference is not None:
        h_ck_matched = match_bin_width_correlation_function(h_correlation_function_reference, h_ck)    
        h_ck_matched.SetName('hCkHist_Matched') 

    smoothen_histogram(h_ck, outfile, xmin=0.0, xmax=0.42, rho=2.)

    outfile.cd()
    h_mixed_event.Write('hMixedEvent')
    h_ck.Write()
    if h_correlation_function_reference is not None:
        h_ck_matched.Write()
    

if __name__ == '__main__':

    M_LI4 = 3.7513 # GeV/c^2
    W_LI4 = 0.005 # GeV/c^2 (li4 width)
    E_TH = 3.7466 # GeV/c^2

    #outfile = TFile.Open('models/li4_contribution.root', 'recreate')
    outfile = TFile.Open('models/li4_contribution_generalized_sill.root', 'recreate')

    sampler = Sampler(M_LI4, W_LI4, 0, E_TH, outfile)
    
    #sampler.sample('sill')
    h_mc_signal = load_hist('/home/galucia/antiLithium4/root_dataframe/output/mc.root', 'InvariantMassAntimatter/hInvariantMassAntimatter')
    #sampler.init_shape_from_mc(h_mc_signal, 'crystal_ball')
    sampler.init_shape_from_mc(h_mc_signal, 'gaus')
    #sampler.sample('sill_gaus_conv_numpy', 10_000_000)
    sampler.sample('sill_generalized_gaus_conv_numpy', 10_000_000)
    
    sampler.save_sampling()
    h_kstar = draw_kstar_profile(sampler.sampled_roo_dataset, outfile)
    draw_kstar_profile_from_hist(sampler.sampled_roo_dataset, outfile)

    #h_mixed_event = load_hist('/home/galucia/Lithium4/preparation/checks/correlation_hadronpid_pass1_pass4_refined_dca.root', 'Correlation/Default/hNormalisedMixedEvent050')
    h_mixed_event = load_hist('/home/galucia/Lithium4/preparation/checks/mixed_event_hadronpid_pass1_pass4_refined_dca_finer_binning.root', 'kstar/hKstar050FinerBinning')
    h_correlation_reference = load_hist('/home/galucia/Lithium4/preparation/checks/correlation_hadronpid_pass1_pass4_refined_dca.root', 'Correlation/Default/hCorrelation050')
    draw_ck_profile_from_hist(h_kstar, h_mixed_event, outfile, h_correlation_reference)
    
    outfile.Close()