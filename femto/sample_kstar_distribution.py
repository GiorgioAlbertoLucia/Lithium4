import numpy as np
from ROOT import TFile, TCanvas, TH1F, \
                kOrange, kBlue, kGreen, \
                RooDataSet, RooRealVar, RooDataHist, RooFit, RooFFTConvPdf

from torchic.roopdf import RooSillPdf
from torchic.core.histogram import load_hist

from core.roopdf_utils import init_roopdf

class Sampler:

    def __init__(self, mass, intrinsic_width, experimental_width, eth, outfile:TFile, var = None):
        self._mass = mass
        self._intrinsic_width = intrinsic_width
        self._experimental_width = experimental_width
        self._eth = eth

        self._var = RooRealVar('invmass', 'm (p+^{3}He)', 3.743, 3.755, 'GeV/#it{c}^{2}')  if var is None else var
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
        self._shape_pdf.fitTo(mc_data_hist, RooFit.Save())

    def sample(self, shape: str, n_samples: int = 1_000_000):
        if shape == 'sill':         self._sample_sill(n_samples)
        if shape == 'gaus_conv':    self._sample_gaus_conv(n_samples)

    def _sample_sill(self, n_samples: int):
        self._pdf_pars = {
            'mass': RooRealVar('mass', 'mass', self._mass),
            'gamma': RooRealVar('mass', 'mass', self._intrinsic_width),
            'eth': RooRealVar('eth', 'eth', self._eth)
        }
        self._pdf = RooSillPdf('sill', 'sill', self._var, *self._pdf_pars.values())
        self._sampled_roo_dataset = self._pdf.generate([self._var], NumEvents=n_samples)

    def _sample_gaus_conv(self, n_samples: int):
        if self._shape_pdf is None:
            raise ValueError('You must first initialise the shape with "init_shape_from_mc".')
        
        mean = RooRealVar('intrinsic_mean', '#mu', 3.751, 'GeV/c')
        mean.setVal(self._shape_pdf_pars['mean'].getVal())
        sigma = RooRealVar('intrinsic_sigma', '#sigma', 0.003, 'GeV/c')
        pdf_intrinsic_resolution, pdf_pars_intrinsic_resolution = init_roopdf('gaus', self._var, mean=mean, sigma=sigma, name='intrinsic')
        mean.setConstant(True)
        sigma.setConstant(True)

        frame = self._var.frame()
        self._shape_pdf.plotOn(frame, LineColor=kOrange-3)
        pdf_intrinsic_resolution.plotOn(frame, LineColor=kGreen+1)
        canvas = TCanvas('c_gaus_conv')
        frame.Draw()
        self._outfile.cd()
        canvas.Write()

        self._var.setRange(3.7, 3.8)
        self._var.setBins(1000, "cache")
        self._pdf = RooFFTConvPdf('convolution', 'convolution', self._var, self._shape_pdf, pdf_intrinsic_resolution, ipOrder=2)
        self._sampled_roo_dataset = self._pdf.generate([self._var], NumEvents=n_samples)
    

    def save_sampling(self):
        
        frame = self._var.frame()
        self._pdf.plotOn(frame, LineColor=kOrange-3)

        canvas = TCanvas('cSill', '', 800, 600)
        frame.Draw()
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

    
    invmass = RooRealVar('invmass', 'm (p+^{3}He)', 3.73, 3.81, 'GeV/#it{c}^{2}')
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

def draw_kstar_profile(roo_dataset: RooDataSet, outfile: TFile) -> None:

    h_kstar = TH1F('hKstar', ';#it{k}* (GeV/#it{c});Counts', 400, 0, 0.4)
    for ientry in range(roo_dataset.numEntries()):
        invmass = roo_dataset.get(ientry)['invmass'].getVal()
        kstar = convert_invmass_to_kstar(invmass)
        h_kstar.Fill(kstar)

    outfile.cd()
    h_kstar.Write()

def draw_kstar_profile_from_hist(roo_dataset: RooDataSet, outfile: TFile) -> None:

    roo_datahist = roo_dataset.binnedClone()
    array_sizes = roo_datahist.getBinnings()
    for size in array_sizes:
        size.Print()

    h_kstar = TH1F('hKstarHist', ';#it{k}* (GeV/#it{c});Counts', 400, 0, 0.4)
    for ientry in range(1, roo_datahist.getBinnings()[0].numBins()):
        invmass = roo_datahist.get(ientry)['invmass'].getVal()
        weight = roo_datahist.weight(ientry)
        kstar = convert_invmass_to_kstar(invmass)
        h_kstar.Fill(kstar, weight)

    outfile.cd()
    h_kstar.Write()

    

if __name__ == '__main__':

    M_LI4 = 3.7513 # GeV/c^2
    W_LI4 = 0.005 # GeV/c^2 (li4 width)
    E_TH = 3.7466 # GeV/c^2

    outfile = TFile.Open('output/sampling.root', 'recreate')

    sampler = Sampler(M_LI4, W_LI4, 0, E_TH, outfile)
    
    #sampler.sample('sill')
    h_mc_signal = load_hist('/home/galucia/antiLithium4/root_dataframe/output/mc.root', 'InvariantMassAntimatter/hInvariantMassAntimatter')
    sampler.init_shape_from_mc(h_mc_signal, 'crystal_ball')
    sampler.sample('gaus_conv')
    
    sampler.save_sampling()
    draw_kstar_profile(sampler.sampled_roo_dataset, outfile)
    draw_kstar_profile_from_hist(sampler.sampled_roo_dataset, outfile)
    
    outfile.Close()