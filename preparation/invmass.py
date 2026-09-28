from math import sqrt

from ROOT import TFile, TH1F, TDirectory, TCanvas, TF1, TBox, TLegend, TTree, TPaveText
from torchic.core.histogram import load_hist
from torchic.utils.root import set_root_object, init_legend
from torchic.utils.colors import get_color
from torchic.physics.particles import PARTICLES

# Invariant mass ranges (GeV/c^2)
# Proton + He3 -> alpha (or other): adjust these to your signal region
NORM_LOW_INVMASS  = 3.78   # right sideband start
NORM_HIGH_INVMASS = 3.88   # right sideband end

NORM_LOW_KSTAR = 0.2
NORM_HIGH_KSTAR = 0.4
OUTSIDE_OF_SIGNAL_KSTAR = 0.25
MAX_KSTAR = 0.4

PEAK_LOW_INVMASS  = 3.74   # signal region start (excluded from normalisation)
PEAK_HIGH_INVMASS = 3.78   # signal region end   (excluded from normalisation)


def invariant_mass_to_kstar(invariant_mass: float,
                            mass_he: float = PARTICLES['He'].mass,
                            mass_proton: float = PARTICLES['Pr'].mass) -> float:
    """Return the two-body relative momentum corresponding to an invariant mass."""
    mass_squared = invariant_mass * invariant_mass
    kallen = (mass_squared - (mass_he + mass_proton) ** 2) \
             * (mass_squared - (mass_he - mass_proton) ** 2)
    if invariant_mass <= 0. or kallen <= 0.:
        return 0.
    return sqrt(kallen) / (2. * invariant_mass)

def kstar_to_inv_mass(kstar: float,
                      mass_he: float = PARTICLES['He'].mass,
                      mass_proton: float = PARTICLES['Pr'].mass) -> float:
    """Return the invariant mass corresponding to a two-body relative momentum."""
    return sqrt((mass_he + mass_proton) ** 2 + 2. * (mass_he + mass_proton) * kstar ** 2 / (mass_he * mass_proton))


def _evaluate_kstar_model(model, kstar: float) -> float:
    if model.InheritsFrom('TH1'):
        return model.Interpolate(kstar)
    if model.InheritsFrom('TF1'):
        return model.Eval(kstar)
    raise TypeError('kstar_model must inherit from TH1 or TF1')


def scale_mixed_by_kstar_model(h_mixed: TH1F, kstar_model,
                               mass_he: float = PARTICLES['He'].mass,
                               mass_proton: float = PARTICLES['Pr'].mass) -> TH1F:
    """Scale an invariant-mass mixed-event histogram with a model defined in k*."""
    h_scaled = h_mixed.Clone(f'{h_mixed.GetName()}KstarModelScaled')
    h_scaled.SetDirectory(0)
    for ibin in range(1, h_scaled.GetNbinsX() + 1):
        kstar = invariant_mass_to_kstar(h_scaled.GetBinCenter(ibin), mass_he, mass_proton)
        scale = _evaluate_kstar_model(kstar_model, kstar)
        h_scaled.SetBinContent(ibin, h_mixed.GetBinContent(ibin) * scale)
        h_scaled.SetBinError(ibin, h_mixed.GetBinError(ibin) * abs(scale))
    return h_scaled


def normalise_and_subtract_with_kstar_model(infile_sames, infile_mixeds, outdir: TDirectory,
                                            mode: str, centrality: str, kstar_model,
                                            rebin: int = 1, suffix: str = '',
                                            mass_he: float = PARTICLES['He'].mass,
                                            mass_proton: float = PARTICLES['Pr'].mass):
    """Apply a k*-space model to mixed events, then normalise and subtract."""
    h_same, h_mixed, _, _ = normalise_and_subtract(
        infile_sames, infile_mixeds, outdir, mode, centrality, rebin, suffix)

    low_bin = h_same.FindBin(NORM_LOW_INVMASS)
    high_bin = h_same.FindBin(NORM_HIGH_INVMASS)
    normalization_factor = h_same.Integral(low_bin, high_bin) / h_mixed.Integral(low_bin, high_bin)
    h_normalised_mixed = h_mixed.Clone(f'hNormalisedMixedEventInvMass{centrality}{suffix}KstarModel')
    h_normalised_mixed.Scale(normalization_factor)
    
    h_mixed_scaled = scale_mixed_by_kstar_model(h_normalised_mixed, kstar_model, mass_he, mass_proton)

    h_signal = h_same.Clone(f'hSignalInvMass{centrality}{suffix}KstarModel')
    h_signal.Add(h_mixed_scaled, -1.)

    if WRITE_HISTS:
        outdir.cd()
        for hist in [h_mixed_scaled, h_normalised_mixed, h_signal]:
            hist.SetTitle(';#it{M} (GeV/#it{c}^{2});')
            hist.Write()
    return h_same, h_mixed, h_mixed_scaled, h_signal

def normalise_and_subtract_no_centrality(infile_sames, infile_mixeds, outdir: TDirectory,
                                         mode: str, rebin: int = 1, suffix: str = ''):

    h_same = infile_sames[0].Get(f'invmass{mode}/hInvariantMass{mode}').Clone(f'hSameEventInvMass{suffix}')
    h_same.SetDirectory(0)
    for infile_same in infile_sames[1:]:
        h_same.Add(infile_same.Get(f'invmass{mode}/hInvariantMass{mode}'))
    if rebin > 1:
        h_same.Rebin(rebin)

    h_mixed = infile_mixeds[0].Get(f'invmass{mode}/hInvariantMass{mode}').Clone(f'hMixedEventInvMass{suffix}')
    h_mixed.SetDirectory(0)
    for infile_mixed in infile_mixeds[1:]:
        h_mixed.Add(infile_mixed.Get(f'invmass{mode}/hInvariantMass{mode}'))
    if rebin > 1:
        h_mixed.Rebin(rebin)
    h_normalised_mixed = h_mixed.Clone(f'hNormalisedMixedEventInvMass{suffix}')

    # Normalise using the right sideband only
    low_bin  = h_same.FindBin(NORM_LOW_INVMASS)
    high_bin = h_same.FindBin(NORM_HIGH_INVMASS)
    normalization_factor = h_same.Integral(low_bin, high_bin) / h_normalised_mixed.Integral(low_bin, high_bin)
    h_normalised_mixed.Scale(normalization_factor)

    h_signal = h_same.Clone(f'hSignalInvMass{suffix}')
    h_signal.Add(h_normalised_mixed, -1.)

    if suffix != '':
        return h_same, h_mixed, h_normalised_mixed, h_signal

    if WRITE_HISTS:
        outdir.cd()
        for hist in [h_same, h_mixed, h_normalised_mixed, h_signal]:
            hist.SetTitle(';#it{M} (GeV/#it{c}^{2});')
            hist.Write()
    return h_same, h_mixed, h_normalised_mixed, h_signal

def normalise_and_subtract(infile_sames, infile_mixeds, outdir: TDirectory,
                           mode: str, centrality: str, rebin: int = 1, suffix: str = ''):

    print(f'Same event histogram: invmass{mode}/hInvariantMass{centrality}{suffix}{mode}')
    print(f'Mixed event histogram: invmass{mode}/hInvariantMass{centrality}{suffix}{mode}')
    
    h_same = infile_sames[0].Get(f'invmass{mode}/hInvariantMass{centrality}{suffix}{mode}').Clone(f'hSameEventInvMass{centrality}{suffix}')
    h_same.SetDirectory(0)
    for infile_same in infile_sames[1:]:
        h_same.Add(infile_same.Get(f'invmass{mode}/hInvariantMass{centrality}{suffix}{mode}'))
    if rebin > 1:
        h_same.Rebin(rebin)

    h_mixed = infile_mixeds[0].Get(f'invmass{mode}/hInvariantMass{centrality}{suffix}{mode}').Clone(f'hMixedEventInvMass{centrality}{suffix}')
    h_mixed.SetDirectory(0)
    for infile_mixed in infile_mixeds[1:]:
        h_mixed.Add(infile_mixed.Get(f'invmass{mode}/hInvariantMass{centrality}{suffix}{mode}'))
    if rebin > 1:
        h_mixed.Rebin(rebin)
    h_normalised_mixed = h_mixed.Clone(f'hNormalisedMixedEventInvMass{centrality}{suffix}')

    # Normalise using the right sideband only
    low_bin  = h_same.FindBin(NORM_LOW_INVMASS)
    high_bin = h_same.FindBin(NORM_HIGH_INVMASS)
    normalization_factor = h_same.Integral(low_bin, high_bin) / h_normalised_mixed.Integral(low_bin, high_bin)
    h_normalised_mixed.Scale(normalization_factor)

    h_signal = h_same.Clone(f'hSignalInvMass{centrality}{suffix}')
    h_signal.Add(h_normalised_mixed, -1.)

    if suffix != '':
        return h_same, h_mixed, h_normalised_mixed, h_signal

    if WRITE_HISTS:
        outdir.cd()
        for hist in [h_same, h_mixed, h_normalised_mixed, h_signal]:
            hist.SetTitle(';#it{M} (GeV/#it{c}^{2});')
            hist.Write()

    return h_same, h_mixed, h_normalised_mixed, h_signal


def plot_invariant_mass(h_same: TH1F, h_normalised_mixed: TH1F, h_signal: TH1F,
                        outdir: TDirectory, centrality: str, suffix: str = ''):

    canvas = TCanvas(f'cInvMass{centrality}{suffix}', '', 800, 600)
    hframe = canvas.DrawFrame(MIN_INVMASS, h_normalised_mixed.GetMinimum(), MAX_INVMASS, max(h_same.GetMaximum(), h_normalised_mixed.GetMaximum()) * 1.2,
                              ';#it{m} (p + ^{3}He + c.c.) (GeV/#it{c}^{2});Counts')
    set_root_object(h_same, marker_style=20, marker_color=get_color(0), line_color=get_color(0))
    set_root_object(h_normalised_mixed, marker_style=24, marker_color=get_color(1), line_color=get_color(1))
    box_normalisation_region = TBox(NORM_LOW_INVMASS, hframe.GetMinimum(), NORM_HIGH_INVMASS, hframe.GetMaximum())
    set_root_object(box_normalisation_region, fill_style=1001, fill_color_alpha=(get_color(6), 0.3))

    h_same.Draw('hist pe1 same')
    h_normalised_mixed.Draw('hist pe1 same')
    box_normalisation_region.Draw('same')

    outdir.cd()
    canvas.Write()

    canvas_signal = TCanvas(f'cSignalInvMass{centrality}{suffix}', '', 800, 600)
    hframe_signal = canvas_signal.DrawFrame(MIN_INVMASS, h_signal.GetMinimum(), MAX_INVMASS, h_signal.GetMaximum() * 1.2,
                                           ';#it{m} (p + ^{3}He + c.c.) (GeV/#it{c}^{2});Counts')
    set_root_object(h_signal, marker_style=20, marker_color=get_color(0), line_color=get_color(0))
    h_signal.Draw('hist pe1 same')
    canvas_signal.Write()
    
    canvas_signal_fit = TCanvas(f'cSignalFitInvMass{centrality}{suffix}', '', 800, 600)
    hframe_signal_fit = canvas_signal_fit.DrawFrame(MIN_INVMASS, h_signal.GetMinimum(), MAX_INVMASS, h_signal.GetMaximum() * 1.2,
                                           ';#it{m} (p + ^{3}He + c.c.) (GeV/#it{c}^{2});Counts')
    set_root_object(h_signal, marker_style=20, marker_color=get_color(0), line_color=get_color(0))
    h_signal.Draw('hist pe1 same')
    OUTSIDE_OF_SIGNAL_INVMASS = kstar_to_inv_mass(OUTSIDE_OF_SIGNAL_KSTAR)  
    pol0_bkg = TF1(f'pol0_bkg_{centrality}{suffix}', '[0]', OUTSIDE_OF_SIGNAL_INVMASS, MAX_INVMASS)
    pol0_bkg.SetParameter(0, 0)
    h_signal.Fit(pol0_bkg, 'RMSQ', '', OUTSIDE_OF_SIGNAL_INVMASS, MAX_INVMASS)
    set_root_object(pol0_bkg, line_color=get_color(1), line_style=2, line_width=2)
    
    text = TPaveText(0.6, 0.7, 0.88, 0.88, 'NDC')
    set_root_object(text, fill_style=0, border_size=0)
    text.SetBorderSize(0)
    text.AddText(f'Centrality: {centrality}%')
    text.AddText(f'pol0 = {pol0_bkg.GetParameter(0):.2f} #pm {pol0_bkg.GetParError(0):.2f}')
    text.AddText(f'#chi^{{2}} / ndf = {pol0_bkg.GetChisquare():.2f} / {pol0_bkg.GetNDF():.2f}')
    text.Draw('same')
    pol0_bkg.Draw('same')
    outdir.cd()
    canvas_signal_fit.Write()

def direct_computation_signal_centrality_integrated(h_sames, h_mixeds, outdir: TDirectory,
                                                    mode: str, centrality: str = '050', suffix: str = '',
                                                    kstar_model=None):
    """
    Sum same-event and mixed-event histograms across centrality classes,
    then normalise and subtract to obtain the signal.
    """
    h_same  = h_sames[0].Clone(f'hSameEventInvMassDirectComputation{centrality}{suffix}')
    h_mixed = h_mixeds[0].Clone(f'hMixedEventInvMassDirectComputation{centrality}{suffix}')

    for h_same_cent, h_mixed_cent in zip(h_sames[1:], h_mixeds[1:]):
        h_same.Add(h_same_cent)
        h_mixed.Add(h_mixed_cent)

    for hist in [h_same, h_mixed]:
        hist.Sumw2()

    h_normalised_mixed = h_mixed.Clone(f'hNormalisedMixedEventInvMassDirectComputation{centrality}{suffix}')

    low_bin  = h_same.FindBin(NORM_LOW_INVMASS)
    high_bin = h_same.FindBin(NORM_HIGH_INVMASS)
    normalization_factor = h_same.Integral(low_bin, high_bin) / h_normalised_mixed.Integral(low_bin, high_bin)
    h_normalised_mixed.Scale(normalization_factor)
    
    if kstar_model is not None:
        h_normalised_mixed_scaled = scale_mixed_by_kstar_model(h_normalised_mixed, kstar_model)
        h_signal = h_same.Clone(f'hSignalInvMassDirectComputation{centrality}{suffix}KstarModel')
        h_signal.Add(h_normalised_mixed_scaled, -1.)
    else:
        h_normalised_mixed_scaled = h_normalised_mixed
        h_signal = h_same.Clone(f'hSignalInvMassDirectComputation{centrality}{suffix}')
        h_signal.Add(h_normalised_mixed, -1.)

    if WRITE_HISTS:
        outdir.cd()
        for hist in [h_same, h_mixed, h_normalised_mixed, h_signal]:
            hist.SetTitle(';#it{M} (GeV/#it{c}^{2});')
            if hist == h_signal:
                hist.SetTitle(';#it{M} (GeV/#it{c}^{2}); Counts')
            hist.Write()
    return h_same, h_normalised_mixed_scaled, h_signal


if __name__ == '__main__':
    
    WRITE_HISTS = False
    NORM_LOW_INVMASS  = kstar_to_inv_mass(NORM_LOW_KSTAR)
    NORM_HIGH_INVMASS = kstar_to_inv_mass(NORM_HIGH_KSTAR)
    MIN_INVMASS = PARTICLES['He'].mass + PARTICLES['Pr'].mass
    MAX_INVMASS = kstar_to_inv_mass(MAX_KSTAR)

    infile_same_paths = [
        'output/PbPb/LHC23_PbPb_pass5_hadronpid_same_invmass.root',
        'output/PbPb/LHC24ar_pass3_hadronpid_same_invmass.root',
        'output/PbPb/LHC25_PbPb_pass1_hadronpid_same_invmass.root',
    ]
    infile_mixed_paths = [
        'output/PbPb/LHC23_PbPb_pass5_hadronpid_event_mixing_invmass.root',
        'output/PbPb/LHC24ar_pass3_hadronpid_event_mixing_invmass.root',
        'output/PbPb/LHC25_PbPb_pass1_hadronpid_event_mixing_invmass.root',
    ]
    outfile_path = (
        #'output/PbPb/invmass_LHC23_PbPb_pass5_hadronpid.root'
        #'output/PbPb/invmass_LHC24ar_pass3_hadronpid.root' 
        #'output/PbPb/invmass_LHC25_PbPb_pass1_hadronpid.root'
        'output/PbPb/invmass_PbPb_hadronpid.root'
        #'output/PbPb/invmass_LHC23_PbPb_pass5_LHC24ar_pass3_hadronpid.root'
        )
    infile_model_path = '/home/galucia/Lithium4/femto/models/lambda_models_LL_10.root'

    infile_sames  = [TFile.Open(p) for p in infile_same_paths]
    infile_mixeds = [TFile.Open(p) for p in infile_mixed_paths]
    outfile = TFile.Open(outfile_path, 'RECREATE')

    for mode in ['Matter', 'Antimatter', '']:

        outdir = outfile.mkdir(f'InvMass{mode}')

        h_sames, h_mixeds, h_normalised_mixeds, h_signals = [], [], [], []

        # Centrality-integrated, no centrality suffix
        #normalise_and_subtract_no_centrality(infile_sames, infile_mixeds, outdir, mode, rebin=2)

        for suffix in ['']:

            outdir_suffix = outdir.mkdir(f'{suffix}' if suffix != '' else 'Default')

            #for centrality in ['010', '1030', '3050', '5080']:
            for centrality in ['010', '1030', '3050', '5080']:

                

                sign = mode if mode != '' else 'Both'
                h_kstar_model = load_hist(infile_model_path, f'{sign}/{centrality}/hLambdaSigmaCorrectedCk_Smeared')
                if h_kstar_model is not None:
                    h_same, h_mixed, h_normalised_mixed, h_signal = normalise_and_subtract_with_kstar_model(
                        infile_sames, infile_mixeds, outdir_suffix,
                        mode, centrality, #rebin=2, 
                        suffix=suffix, kstar_model=h_kstar_model)
                else:
                    h_same, h_mixed, h_normalised_mixed, h_signal = normalise_and_subtract(
                        infile_sames, infile_mixeds, outdir_suffix,
                        mode, centrality, #rebin=2, 
                        suffix=suffix)

                plot_invariant_mass(h_same, h_normalised_mixed, h_signal,
                                    outdir_suffix, centrality, suffix)

                h_sames.append(h_same)
                h_mixeds.append(h_mixed)
                h_normalised_mixeds.append(h_normalised_mixed)
                h_signals.append(h_signal)

            # Centrality-integrated signal via direct sum
            direct_computation_signal_centrality_integrated(
                h_sames[:3], h_mixeds[:3], outdir_suffix, mode, '050', suffix)
            
            h_kstar_model = load_hist(infile_model_path, f'{sign}/1050/hLambdaSigmaCorrectedCk_Smeared')
            h_same, h_normalised_mixed, h_signal = direct_computation_signal_centrality_integrated(
                h_sames[1:3], h_mixeds[1:3], outdir_suffix, mode, '1050', suffix,
                kstar_model=h_kstar_model)
            plot_invariant_mass(h_same, h_normalised_mixed, h_signal,
                outdir_suffix, '1050', suffix)
            
            direct_computation_signal_centrality_integrated(
                h_sames, h_mixeds, outdir_suffix, mode, '080', suffix)
            direct_computation_signal_centrality_integrated(
                h_sames[1:], h_mixeds[1:], outdir_suffix, mode, '1080', suffix)

            h_sames.clear()
            h_mixeds.clear()
            h_normalised_mixeds.clear()
            h_signals.clear()

    outfile.Close()