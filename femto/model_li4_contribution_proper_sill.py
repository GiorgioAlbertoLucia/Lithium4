import numpy as np
from tqdm import tqdm
from dataclasses import dataclass

from ROOT import TFile, TCanvas, TH1F, TLegend, TH2F, TTree, TLatex, gStyle, \
                kOrange, kBlue, kGreen, \
                RooDataSet, RooRealVar, RooDataHist, RooFit, RooFFTConvPdf, RooKeysPdf

from torchic.roopdf import RooSillPdf, RooSillGeneralizedKstarPdf
from torchic.core.histogram import load_hist
from torchic.utils.root import set_root_object, init_legend
from torchic.utils.colors import get_color

from core.roopdf_utils import init_roopdf

def convert_invmass_to_kstar(invmass: float):
    
    M_PR = 0.938272 # GeV/c^2
    M_HE = 2.80839 # GeV/c^2
    lambda_func = lambda x, y, z: x**2 + y**2 + z**2 -2*x*y - 2*x*z - 2*y*z
    numerator = np.sqrt(lambda_func(invmass**2, M_PR**2, M_HE**2)) if lambda_func(invmass**2, M_PR**2, M_HE**2) >= 0 else 0
    denominator = 2 * invmass
    return numerator / denominator  


class Sampler:

    def __init__(self, mass, intrinsic_width, experimental_width, 
                 mass_daughter_1, mass_daughter_2, l,
                 outfile:TFile, var_invmass = None, var_kstar = None):
        self._mass = mass
        self._intrinsic_width = intrinsic_width
        self._experimental_width = experimental_width
        self._mass_daughter_1 = mass_daughter_1
        self._mass_daughter_2 = mass_daughter_2
        self._l = l

        self._var_invmass = RooRealVar('invmass', 'm (p+^{3}He)', 3.737, 3.865, 'GeV/#it{c}^{2}')  if var_invmass is None else var_invmass
        self._var_kstar = RooRealVar('kstar', '#it{k}* (p+^{3}He)', 0., 0.4, 'GeV/#it{c}')  if var_kstar is None else var_kstar
        self._sampled_roo_dataset = None
        self._pdf = None
        self._pdf_pars = {}
        
        self._shape_pdf_invmass = None      # pdf to get the shape from the MC (intrinsic resolution is missing)
        self._shape_pdf_invmass_pars = {}
        self._shape_pdf = None      # pdf to get the shape from the MC (intrinsic resolution is missing)
        self._shape_pdf_pars = {}

        self._outfile = outfile
    
    @property
    def sampled_roo_dataset(self):
        return self._sampled_roo_dataset
    
    def init_shape_from_mc(self, h_mc_signal:TH1F, shape_hypothesis: str):

        mean = RooRealVar('experimental_mean', '#mu', 3.751, 3.746, 3.756, 'GeV/c^{2}')
        sigma = RooRealVar('experimental_sigma', '#sigma', 0.01, 0.0001, 0.1, 'GeV/c^{2}')
        self._shape_pdf_invmass, self._shape_pdf_invmass_pars = init_roopdf(shape_hypothesis, self._var_invmass, mean=mean, sigma=sigma, name='experimental')

        mc_data_hist = RooDataHist('mc_data_hist', 'mc_data_hist', [self._var_invmass], Import=h_mc_signal)
        self._shape_pdf_invmass.fitTo(mc_data_hist, RooFit.Save(), Range=(3.750, 3.753))

        frame = self._var_invmass.frame()
        mc_data_hist.plotOn(frame)
        self._shape_pdf_invmass.plotOn(frame, LineColor=kBlue+2)
        self._shape_pdf_invmass.paramOn(frame)
        self._var_invmass.setRange(3.737, 3.865)
        canvas = TCanvas('c_experimental_points')
        frame.Draw()
        self._outfile.cd()
        canvas.Write()

    def init_shape_from_mc_in_kstar(self, h_mc_signal:TH1F, shape_hypothesis: str, n_samples: int = 1_000_000):

        h_mc_signal_kstar = TH1F('h_mc_signal_kstar', ';#it{k}* (GeV/#it{c});Counts', 400, 0, 0.4)
        for isample in tqdm(range(n_samples)):
            invmass = h_mc_signal.GetRandom()
            kstar = convert_invmass_to_kstar(invmass)
            h_mc_signal_kstar.Fill(kstar)
        
        mean = RooRealVar('experimental_mean_kstar', '#mu', 0.081, 0.001, 0.1, 'GeV/c')
        sigma = RooRealVar('experimental_sigma_kstar', '#sigma', 0.01, 0.0001, 0.1, 'GeV/c')
        self._shape_pdf, self._shape_pdf_pars = init_roopdf(shape_hypothesis, self._var_kstar, mean=mean, sigma=sigma, name='experimental_kstar')
        mc_data_hist_kstar = RooDataHist('mc_data_hist_kstar', 'mc_data_hist_kstar', [self._var_kstar], Import=h_mc_signal_kstar)
        self._shape_pdf.fitTo(mc_data_hist_kstar, RooFit.Save(), Range=(0.06, 0.1))

        frame = self._var_kstar.frame()
        mc_data_hist_kstar.plotOn(frame)
        self._shape_pdf.plotOn(frame, LineColor=kBlue+2)
    
        canvas = TCanvas('c_experimental_points_kstar')
        frame.Draw()
        self._outfile.cd()
        canvas.Write()

    def sample(self, shape: str, n_samples: int = 1_000_000):
        if shape == 'sill_gaus_conv_numpy': self._sample_sill_gaus_conv_numpy(n_samples)

    def _sample_sill_gaus_conv_numpy(self, n_samples: int):

        if self._shape_pdf is None:
            raise ValueError('You must first initialise the shape with "init_shape_from_mc".')

        frame = self._var_kstar.frame()
        self._shape_pdf.plotOn(frame, LineColor=kBlue+2)
        canvas = TCanvas('c_experimental_shape')
        frame.Draw()
        self._outfile.cd()
        canvas.Write()
        del frame

        pdf_pars_intrinsic_resolution = {
            'mass': RooRealVar('mass', 'mass', self._mass),
            'gamma': RooRealVar('gamma', 'gamma', self._intrinsic_width),
            'mass_daughter_1': RooRealVar('mass_daughter_1', 'mass_daughter_1', self._mass_daughter_1),
            'mass_daughter_2': RooRealVar('mass_daughter_2', 'mass_daughter_2', self._mass_daughter_2),
            'l': RooRealVar('l', 'l', self._l),
        }
        pdf_intrinsic_resolution = RooSillGeneralizedKstarPdf('sill', 'sill', self._var_kstar, *pdf_pars_intrinsic_resolution.values())

        frame = self._var_kstar.frame()
        pdf_intrinsic_resolution.plotOn(frame, LineColor=kOrange-3)
        canvas = TCanvas('c_intrinsic_shape')
        frame.Draw()
        self._outfile.cd()
        canvas.Write()

        sill_dataset = pdf_intrinsic_resolution.generate([self._var_kstar], NumEvents=n_samples)
        experimental_sigma = self._shape_pdf_pars['sigma'].getVal()
        self._sampled_roo_dataset = RooDataSet('smeared_data', 'smeared_data', [self._var_kstar])

        for ientry in tqdm(range(sill_dataset.numEntries())):
            
            invmass_true = sill_dataset.get(ientry)['kstar'].getVal()
            kstar_smeared = np.random.normal(invmass_true, experimental_sigma)
            
            if self._var_kstar.getMin() <= kstar_smeared <= self._var_kstar.getMax():
                self._var_kstar.setVal(kstar_smeared)
                self._sampled_roo_dataset.add([self._var_kstar])


        self._pdf = pdf_intrinsic_resolution
        self._pdf_pars = pdf_pars_intrinsic_resolution

    def save_sampling(self):
        
        frame = self._var_kstar.frame()
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


def draw_profile(roo_dataset: RooDataSet, outfile: TFile, name:str) -> TH1F:

    h_kstar = TH1F('hKstar', ';#it{k}* (GeV/#it{c});Counts', 400, 0, 0.4)
    for ientry in range(roo_dataset.numEntries()):
        kstar = roo_dataset.get(ientry)['kstar'].getVal()
        h_kstar.Fill(kstar)

    h_kstar.SetBinContent(1, 0.)    # in the conversion, this bin behaves as an underflow bin

    outfile.cd()
    h_kstar.Write(f'hKstar_{name}')

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

    if outfile is not None:
        smoothen_histogram(h_ck, outfile, xmin=0.0, xmax=0.42, rho=2.)

        outfile.cd()
        h_mixed_event.Write('hMixedEvent')
        h_ck.Write()

    h_ck_integral = h_ck.Integral()
    h_ck_normalized = h_ck.Clone('hCkHist_Normalized')
    h_ck_normalized.Scale(1./h_ck_integral if h_ck_integral > 0 else 1.)

    gStyle.SetPadTickX(1)
    gStyle.SetPadTickY(1)
    
    canvas = TCanvas('cCk', 'cCk', 700, 800)
    canvas.SetLeftMargin(0.2)
    canvas.SetRightMargin(0.05)
    canvas.SetBottomMargin(0.15)

    hframe = canvas.DrawFrame(0, 0., 0.4, h_ck_normalized.GetMaximum()*1.1, ';#it{k}* (GeV/#it{c});#it{C}_{^{4}Li}(#it{k}*)')
    hframe.GetXaxis().SetTitleSize(0.05)
    hframe.GetYaxis().SetTitleSize(0.05)
    hframe.GetXaxis().SetLabelSize(0.045)
    hframe.GetYaxis().SetLabelSize(0.045)
    #hframe.GetYaxis().SetTitleOffset(0.9)
    
    
    set_root_object(h_ck_normalized, line_color=4, line_width=2)
    h_ck_normalized.Draw('hist same')

    latex = TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.045)
    latex.SetTextFont(42)
    latex.DrawLatex(0.44, 0.6, f'ALICE Simulation')
    latex.DrawLatex(0.44, 0.54, f'Pb#minusPb #sqrt{{#it{{s}}_{{NN}}}} = 5.36 TeV')
    latex.DrawLatex(0.44, 0.48, '^{4}#bar{Li} #rightarrow #bar{p} + ^{3}#bar{He}')

    canvas.SaveAs('models/ck_li4.pdf')


    if h_correlation_function_reference is not None and outfile is not None:
        h_ck_matched.Write()
    

if __name__ == '__main__':

    MASS_PROTON = 0.938272 # GeV/c^2
    MASS_HE3 = 2.80839 # GeV/c^2
    
    @dataclass
    class nucleus:
        mass: float
        width: float
        l: int
        normalization: float
    
    ### --------------------------------------------------
    ### Thermal model yields in PbPb at 5.36 TeV
    ### T = 155 MeV, muB = 0.0007 GeV, 0-5% central
    ### PDG 1000030040,   density = 7.8425e-10 fm^-3
    ### PDG 100003004001, density = 4.6964e-10 fm^-3
    ### PDG 100003004002, density = 1.5490e-10 fm^-3
    ### PDG 100003004003, density = 4.6255e-10 fm^-3
    
    li4_states = {'g.s.': nucleus(mass=3.75073, width=0.006, l=2, normalization=7.8425),
                  'ex. 1': nucleus(mass=3.75105, width=0.00735, l=1, normalization=4.6964),
                  'ex. 2': nucleus(mass=3.75281, width=0.00935, l=0, normalization=1.5490),
                  'ex. 3': nucleus(mass=3.75358, width=0.01351, l=1, normalization=4.6255)}

    outfile = TFile.Open('models/li4_contribution_proper_sill.root', 'recreate')
    h_kstars = []

    for state_name, state in li4_states.items():
        outdir = outfile.mkdir(state_name)
        print(f'Sampling {state_name} state...')
        mass, width, l = state.mass, state.width, state.l
        sampler = Sampler(mass=mass, intrinsic_width=width, experimental_width=0., 
                          mass_daughter_1=MASS_PROTON, mass_daughter_2=MASS_HE3, l=l,
                          outfile=outdir)

        h_mc_signal = load_hist('/home/galucia/antiLithium4/root_dataframe/output/mc.root', 'InvariantMassAntimatter/hInvariantMassAntimatter')
        sampler.init_shape_from_mc_in_kstar(h_mc_signal, 'crystal_ball')
        sampler.sample('sill_gaus_conv_numpy', 10_000_000)

        sampler.save_sampling()
        h_kstar = draw_profile(sampler.sampled_roo_dataset, outfile, state_name)
        h_kstars.append(h_kstar)

    GROUND_STATE_NAME = 'g.s.'
    state_names = list(li4_states.keys())

    h_kstar = None
    h_kstar_ground_state = None
    for ihist, (state_name, hist, states) in enumerate(zip(state_names, h_kstars, li4_states.values())):
        hist.Scale(states.normalization)
        if state_name == GROUND_STATE_NAME:
            h_kstar_ground_state = hist.Clone('hKstar_GroundStateOnly')
        if ihist == 0:
            h_kstar = hist.Clone('hKstar')
        else:
            h_kstar.Add(hist)
    outfile.cd()
    h_kstar.Write('hKstar')
    h_kstar_ground_state.Write('hKstar_GroundStateOnly')

    h_mixed_event = load_hist('/home/galucia/Lithium4/preparation/checks/mixed_event_hadronpid_pass1_pass4_refined_dca_finer_binning.root', 'kstar/hKstar050FinerBinning')
    h_correlation_reference = load_hist('/home/galucia/Lithium4/preparation/checks/correlation_hadronpid_pass1_pass4_refined_dca.root', 'Correlation/Default/hCorrelation050')
    h_kstar = load_hist('models/li4_contribution_proper_sill.root', 'hKstar')
    draw_ck_profile_from_hist(h_kstar, h_mixed_event, outfile, h_correlation_reference)
    
    canvas = TCanvas('cAllStates', 'cAllStates', 800, 600)
    legend = init_legend(0.6, 0.4, 0.8, 0.6)
    for ihist, hist in enumerate(h_kstars):
        hist.SetLineColor(get_color(ihist))
        hist.Draw('hist same')
        legend.AddEntry(hist, f'{list(li4_states.keys())[ihist]}', 'l')
    legend.Draw('same')
    outfile.cd()
    canvas.Write()
    
    if outfile is not None:
        outfile.Close()