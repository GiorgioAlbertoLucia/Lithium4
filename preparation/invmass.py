from ROOT import TFile, TH1F, TDirectory, TCanvas, TF1
from torchic.utils.root import set_root_object

# Invariant mass ranges (GeV/c^2)
# Proton + He3 -> alpha (or other): adjust these to your signal region
NORM_LOW_INVMASS  = 3.78   # right sideband start
NORM_HIGH_INVMASS = 3.88   # right sideband end

PEAK_LOW_INVMASS  = 3.74   # signal region start (excluded from normalisation)
PEAK_HIGH_INVMASS = 3.78   # signal region end   (excluded from normalisation)


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

    outdir.cd()
    for hist in [h_same, h_mixed, h_normalised_mixed, h_signal]:
        hist.SetTitle(';#it{M} (GeV/#it{c}^{2});')
        hist.Write()

    return h_same, h_mixed, h_normalised_mixed, h_signal


def plot_invariant_mass(h_same: TH1F, h_normalised_mixed: TH1F, h_signal: TH1F,
                        outdir: TDirectory, centrality: str, suffix: str = ''):

    canvas = TCanvas(f'cInvMass{centrality}{suffix}', '', 800, 600)
    set_root_object(h_same, marker_style=20)
    set_root_object(h_normalised_mixed, marker_style=24)
    set_root_object(h_signal, marker_style=20)

    h_same.Draw('hist pe1')
    h_normalised_mixed.Draw('hist pe1 same')

    outdir.cd()
    canvas.Write()

    canvas_signal = TCanvas(f'cSignalInvMass{centrality}{suffix}', '', 800, 600)
    h_signal.Draw('hist pe1')
    canvas_signal.Write()


def direct_computation_signal_centrality_integrated(h_sames, h_mixeds, outdir: TDirectory,
                                                    mode: str, centrality: str = '050', suffix: str = ''):
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

    h_signal = h_same.Clone(f'hSignalInvMassDirectComputation{centrality}{suffix}')
    h_signal.Add(h_normalised_mixed, -1.)

    outdir.cd()
    for hist in [h_same, h_mixed, h_normalised_mixed, h_signal]:
        hist.SetTitle(';#it{M} (GeV/#it{c}^{2});')
        if hist == h_signal:
            hist.SetTitle(';#it{M} (GeV/#it{c}^{2}); Counts')
        hist.Write()

    return h_same, h_mixed, h_signal


if __name__ == '__main__':

    infile_same_paths = [
        'output/PbPb/LHC23_PbPb_pass5_hadronpid_same.root',
        'output/PbPb/LHC24ar_pass3_hadronpid_same.root',
        'output/PbPb/LHC25_PbPb_pass1_hadronpid_same.root',
    ]
    infile_mixed_paths = [
        'output/PbPb/LHC23_PbPb_pass5_hadronpid_event_mixing.root',
        'output/PbPb/LHC24ar_pass3_hadronpid_event_mixing.root',
        'output/PbPb/LHC25_PbPb_pass1_hadronpid_event_mixing.root',
    ]
    outfile_path = (
        #'output/PbPb/invmass_LHC23_PbPb_pass5_hadronpid.root'
        #'output/PbPb/invmass_LHC24ar_pass3_hadronpid.root' 
        #'output/PbPb/invmass_LHC25_PbPb_pass1_hadronpid.root'
        'output/PbPb/invmass_PbPb_hadronpid.root'
        #'output/PbPb/invmass_LHC23_PbPb_pass5_LHC24ar_pass3_hadronpid.root'
        )

    infile_sames  = [TFile.Open(p) for p in infile_same_paths]
    infile_mixeds = [TFile.Open(p) for p in infile_mixed_paths]
    outfile = TFile.Open(outfile_path, 'RECREATE')

    for mode in ['Matter', 'Antimatter']:

        outdir = outfile.mkdir(f'InvMass{mode}')

        h_sames, h_mixeds, h_normalised_mixeds, h_signals = [], [], [], []

        # Centrality-integrated, no centrality suffix
        #normalise_and_subtract_no_centrality(infile_sames, infile_mixeds, outdir, mode, rebin=2)

        for suffix in ['']:

            outdir_suffix = outdir.mkdir(f'{suffix}' if suffix != '' else 'Default')

            for centrality in ['010', '1030', '3050', '5080']:

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
            direct_computation_signal_centrality_integrated(
                h_sames[1:3], h_mixeds[1:3], outdir_suffix, mode, '1050', suffix)
            direct_computation_signal_centrality_integrated(
                h_sames, h_mixeds, outdir_suffix, mode, '080', suffix)
            direct_computation_signal_centrality_integrated(
                h_sames[1:], h_mixeds[1:], outdir_suffix, mode, '1080', suffix)

            h_sames.clear()
            h_mixeds.clear()
            h_normalised_mixeds.clear()
            h_signals.clear()

    outfile.Close()