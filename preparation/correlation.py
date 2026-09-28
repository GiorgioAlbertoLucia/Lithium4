from array import array
from ROOT import TFile, TH1F, TDirectory, TCanvas, TDirectory, TF1
from torchic.utils.root import set_root_object
from torchic.utils.terminal_colors import TerminalColors as tc

NORM_LOW_KSTAR = 0.2 # 0.25
NORM_HIGH_KSTAR = 0.4 # 0.8
NBINS_KSTAR = 40 # 20

def get_adaptive_bin_edges(x_max: float, mid_edge: float = 0.15,
                            width_low: float = 0.02, width_high: float = 0.01) -> array:
    """
    Build variable-width bin edges: `width_low`-wide bins from 0 up to `mid_edge`,
    then `width_high`-wide bins from `mid_edge` up to `x_max`.
    Returns an array('d', ...) suitable for TH1::Rebin(nbins, name, edges).
    """
    edges = []
    x = 0.
    while x < mid_edge - 1e-9:
        edges.append(round(x, 10))
        x += width_low

    x = mid_edge
    while x < x_max - 1e-9:
        edges.append(round(x, 10))
        x += width_high
    edges.append(round(x_max, 10))

    return array('d', edges)

def normalise_histograms_and_compute_correlation_no_centrality(infile_sames, infile_mixeds, outdir:TDirectory,
                                                 mode:str, rebin:int=1, suffix:str='', adaptive_binning:bool=False):

    h_same = infile_sames[0].Get(f'QA/hKstar{mode}').Clone(f'hSameEvent')
    h_same.SetDirectory(0)  # Detach from file to avoid issues with deletion
    for infile_same in infile_sames[1:]:
        h_same.Add(infile_same.Get(f'QA/hKstar{mode}'))
    if adaptive_binning:
        bin_edges = get_adaptive_bin_edges(h_same.GetXaxis().GetXmax())
        h_same = h_same.Rebin(len(bin_edges) - 1, 'hSameEvent', bin_edges)
        h_same.SetDirectory(0)
    elif rebin > 1:
        h_same.Rebin(rebin)

    print(f'Getting: {tc.CYAN+tc.UNDERLINE}QA/hKstar{mode}{tc.RESET} from {len(infile_mixeds)} mixed event files')
    h_mixed = infile_mixeds[0].Get(f'QA/hKstar{mode}').Clone(f'hMixedEvent')
    h_mixed.SetDirectory(0)  # Detach from file to avoid issues with deletion
    for infile_mixed in infile_mixeds[1:]:
        h_mixed.Add(infile_mixed.Get(f'QA/hKstar{mode}'))
    if adaptive_binning:
        bin_edges = get_adaptive_bin_edges(h_mixed.GetXaxis().GetXmax())
        h_mixed = h_mixed.Rebin(len(bin_edges) - 1, 'hMixedEvent', bin_edges)
        h_mixed.SetDirectory(0)
    elif rebin > 1:
        h_mixed.Rebin(rebin)
    h_normalised_mixed = h_mixed.Clone(f'hNormalisedMixedEvent')

    low_bin = h_same.FindBin(NORM_LOW_KSTAR)
    high_bin = h_same.FindBin(NORM_HIGH_KSTAR)
    normalization_factor = h_same.Integral(low_bin, high_bin) / h_normalised_mixed.Integral(low_bin, high_bin) 
    h_normalised_mixed.Scale(normalization_factor)

    h_corr = h_same.Clone(f'hCorrelation')
    h_corr.Divide(h_normalised_mixed)

    if suffix != '':
        return h_same, h_mixed, h_normalised_mixed, h_corr
        
    outdir.cd()
    for hist in [h_same, h_mixed, h_normalised_mixed, h_corr]:
        hist.SetTitle(';#it{k}* (GeV/#it{c});')
        hist.Write()

    return h_same, h_mixed, h_normalised_mixed, h_corr

def normalise_histograms_and_compute_correlation(infile_sames, infile_mixeds, outdir:TDirectory,
                                                 mode:str, centrality:str, rebin:int=1, suffix:str='',
                                                 adaptive_binning:bool=False):

    h_same = infile_sames[0].Get(f'kstar{mode}/hKstar{centrality}{suffix}{mode}').Clone(f'hSameEvent{centrality}')
    h_same.SetDirectory(0)  # Detach from file to avoid issues with deletion
    for infile_same in infile_sames[1:]:
        h_same.Add(infile_same.Get(f'kstar{mode}/hKstar{centrality}{suffix}{mode}'))
    if adaptive_binning:
        bin_edges = get_adaptive_bin_edges(h_same.GetXaxis().GetXmax())
        h_same = h_same.Rebin(len(bin_edges) - 1, f'hSameEvent{centrality}', bin_edges)
        h_same.SetDirectory(0)
    elif rebin > 1:
        h_same.Rebin(rebin)

    h_mixed = infile_mixeds[0].Get(f'kstar{mode}/hKstar{centrality}{suffix}{mode}').Clone(f'hMixedEvent{centrality}')
    h_mixed.SetDirectory(0)  # Detach from file to avoid issues with deletion
    for infile_mixed in infile_mixeds[1:]:
        h_mixed.Add(infile_mixed.Get(f'kstar{mode}/hKstar{centrality}{suffix}{mode}'))
    if adaptive_binning:
        bin_edges = get_adaptive_bin_edges(h_mixed.GetXaxis().GetXmax())
        h_mixed = h_mixed.Rebin(len(bin_edges) - 1, f'hMixedEvent{centrality}', bin_edges)
        h_mixed.SetDirectory(0)
    elif rebin > 1:
        h_mixed.Rebin(rebin)
    h_normalised_mixed = h_mixed.Clone(f'hNormalisedMixedEvent{centrality}')

    low_bin = h_same.FindBin(NORM_LOW_KSTAR)
    high_bin = h_same.FindBin(NORM_HIGH_KSTAR)
    normalization_factor = h_same.Integral(low_bin, high_bin) / h_normalised_mixed.Integral(low_bin, high_bin)
    h_normalised_mixed.Scale(normalization_factor)

    h_corr = h_same.Clone(f'hCorrelation{centrality}')
    h_corr.Divide(h_normalised_mixed)

    #if suffix != '':
    #    return h_same, h_mixed, h_normalised_mixed, h_corr
        
    outdir.cd()
    for hist in [h_same, h_mixed, h_normalised_mixed, h_corr]:
        hist.SetTitle(';#it{k}* (GeV/#it{c});')
        hist.Write()

    return h_same, h_mixed, h_normalised_mixed, h_corr

def fit_normalisation_region(h_corr:TH1F, outdir:TDirectory, mode:str, centrality:str, suffix:str=''):

    constant = TF1(f'func{centrality}', 'pol0', NORM_LOW_KSTAR, NORM_HIGH_KSTAR)
    h_corr.Fit(constant, 'RMS+')

    canvas = TCanvas(f'cCorrelation{centrality}')
    set_root_object(h_corr, marker_style=20)
    h_corr.Draw('hist pe1')
    constant.Draw('same')
    outdir.cd()
    canvas.Write()


def correlation_function_centrality_integrated(h_sames, h_normalised_mixeds, outdir:TDirectory,
                                               mode:str, centrality:str = '050', suffix:str=''):
    """
    Compute the correlation function for the centrality integrated histograms.
    """
    h_same = h_sames[0].Clone(f'hSameEvent{centrality}')
    h_mixed = h_normalised_mixeds[0].Clone(f'hNormalisedMixedEvent{centrality}')

    for h_same_centrality, h_mixed_centrality in zip(h_sames[1:], h_normalised_mixeds[1:]):
        h_same.Add(h_same_centrality)
        h_mixed.Add(h_mixed_centrality)

    h_corr = h_same.Clone(f'hCorrelation{centrality}')
    h_corr.Divide(h_mixed)

    outdir.cd()
    for hist in [h_same, h_mixed, h_normalised_mixed, h_corr]:
        hist.SetTitle(';#it{k}* (GeV/#it{c});')
        if hist == h_corr:  hist.SetTitle(';#it{k}* (GeV/#it{c}); C(#it{k}*)')
        hist.Write()

    return h_same, h_mixed, h_corr

def direct_computation_correlation_function_centrality_integrated(h_sames, h_mixeds, outdir:TDirectory,
                                               mode:str, centrality:str = '050', suffix:str=''):
    """
    Compute the correlation function for the centrality integrated histograms.
    """
    h_same = h_sames[0].Clone(f'hSameEventDirectComputation{centrality}')
    h_mixed = h_mixeds[0].Clone(f'hMixedEventDirectComputation{centrality}')

    for h_same_centrality, h_mixed_centrality in zip(h_sames[1:], h_mixeds[1:]):
        h_same.Add(h_same_centrality)
        h_mixed.Add(h_mixed_centrality)
    
    for hist in [h_same, h_mixed]:
        hist.Sumw2()  # Ensure proper error calculation

    h_normalised_mixed = h_mixed.Clone(f'hNormalisedMixedEventDirectComputation{centrality}')

    low_bin = h_same.FindBin(NORM_LOW_KSTAR)
    high_bin = h_same.FindBin(NORM_HIGH_KSTAR)
    normalization_factor = h_same.Integral(low_bin, high_bin) / h_normalised_mixed.Integral(low_bin, high_bin)
    h_normalised_mixed.Scale(normalization_factor)

    h_corr = h_same.Clone(f'hCorrelationDirectComputation{centrality}')
    h_corr.Divide(h_normalised_mixed)

    outdir.cd()
    for hist in [h_same, h_mixed, h_normalised_mixed, h_corr]:
        hist.SetTitle(';#it{k}* (GeV/#it{c});')
        if hist == h_corr:  hist.SetTitle(';#it{k}* (GeV/#it{c}); C(#it{k}*)')
        hist.Write()

    return h_same, h_mixed, h_corr

def wighted_sum_correlation_function_centrality_integrated(h_correlations, outdir:TDirectory,
                                               mode:str, centrality:str = '050', suffix:str=''):
    """
    Compute the correlation function for the centrality integrated histograms.
    """

    h_correlation_weighted = h_correlations[0].Clone(f'hCorrelationWeightedSum{centrality}')

    for ibin in range(1, h_correlation_weighted.GetNbinsX()+1):

        value = 0.
        error = 0.
        total_error = 0.
        num = 0

        for h_correlation in h_correlations:
            ivalue = h_correlation.GetBinContent(ibin)
            ierror = h_correlation.GetBinError(ibin)  

            value += ivalue / ierror if ierror > 0 else 0.
            total_error += 1. / ierror if ierror > 0 else 0.
            num += 1
        
        if total_error > 0:
            value /= total_error
            error = (num)**0.5 / total_error
            h_correlation_weighted.SetBinContent(ibin, value)
            h_correlation_weighted.SetBinError(ibin, error)
    

    outdir.cd()
    h_correlation_weighted.SetTitle(';#it{k}* (GeV/#it{c}); C(#it{k}*)')
    h_correlation_weighted.Write()

    return h_correlation_weighted
        

if __name__ == '__main__':
    
    #infile_same_path = 'output/same_event.root'
    #infile_mixed_path = 'output/mixed_event.root'
    #outfile_path = 'output/correlation.root'

    infile_same_paths = [
        #'output/PbPb/LHC23_PbPb_pass5_hadronpid_same.root',
        #'output/PbPb/LHC24ar_pass3_hadronpid_same.root',
        #'output/PbPb/LHC25_PbPb_pass1_hadronpid_same.root',
        #'checks/PbPb/LHC25_PbPb_pass1_hadronpid_same_skip_he3_pt_bin.root',
        'output/pp/LHC24_pass1_skimmed_hadronpid_same.root',
        'output/pp/LHC25_pass1_skimmed_hadronpid_same.root',
        ]
    infile_mixed_paths = [
        #'output/PbPb/LHC23_PbPb_pass5_hadronpid_event_mixing.root',
        #'output/PbPb/LHC24ar_pass3_hadronpid_event_mixing.root',
        #'output/PbPb/LHC25_PbPb_pass1_hadronpid_event_mixing.root',
        #'checks/PbPb/LHC25_PbPb_pass1_hadronpid_event_mixing_skip_he3_pt_bin.root',
        #'checks/PbPb/LHC25_PbPb_pass1_hadronpid_event_mixing_depth5.root',
        'output/pp/LHC24_pass1_skimmed_hadronpid_event_mixing.root',
        'output/pp/LHC25_pass1_skimmed_hadronpid_event_mixing.root',
        ]
    outfile_path = (
        #'output/PbPb/correlation_LHC23_PbPb_pass5_hadronpid.root'
        #'output/PbPb/correlation_LHC24ar_pass3_hadronpid.root' 
        #'output/PbPb/correlation_LHC25_PbPb_pass1_hadronpid.root'
        #'output/PbPb/correlation_PbPb_hadronpid.root'
        #'checks/PbPb/correlation_LHC25_PbPb_pass1_hadronpid_depth5.root'
        #'checks/PbPb/correlation_LHC25_PbPb_pass1_hadronpid_skip_he3_pt_bin.root'
        #'output/PbPb/correlation_LHC23_PbPb_pass5_LHC24ar_pass3_hadronpid.root'
        'output/pp/correlation_pp_hadronpid.root'
        )
    
    ADAPTIVE_BINNING = False

    infile_sames = [TFile.Open(infile_same_path) for infile_same_path in infile_same_paths]
    infile_mixeds = [TFile.Open(infile_mixed_path) for infile_mixed_path in infile_mixed_paths]
    outfile = TFile.Open(outfile_path, 'RECREATE')   

    for mode in ['', 'Matter', 'Antimatter']:
    
        outdir = outfile.mkdir(f'Correlation{mode}')

        h_sames, h_mixeds, h_normalised_mixeds, h_corrs = [], [], [], []

        normalise_histograms_and_compute_correlation_no_centrality(infile_sames, infile_mixeds, outdir, mode, rebin=2, adaptive_binning=ADAPTIVE_BINNING)
        
        if True:
            continue
            
        for suffix in ['', 
                       #'PLiGreaterThan3',
                       #'SharedUnder50', 'SharedUnder40', 'SharedUnder30', 'PLi2', 'PLi3', 'PLi4', 'PLi5',
                       #'PLiUnder2', 
                       #'PLiUnder3', 'PLiUnder4', 'PLiUnder5', 'PtHadronUnder1p0', 'PtHadronUnder1p1', 'PtHadronUnder1p2'
                       ]:
            
            outdir_suffix = outdir.mkdir(f'{suffix}' if suffix != '' else 'Default')
            
            
            for centrality in ['010', '1030', '3050', '5080']:

                h_same, h_mixed, h_normalised_mixed, h_corr = normalise_histograms_and_compute_correlation(
                                                                    infile_sames, infile_mixeds, outdir_suffix, 
                                                                    mode, centrality, rebin=2, suffix=suffix,
                                                                    adaptive_binning=ADAPTIVE_BINNING)
                #if suffix == '':
                fit_normalisation_region(h_corr, outdir_suffix, mode, centrality, suffix)
                h_sames.append(h_same)
                h_mixeds.append(h_mixed)
                h_normalised_mixeds.append(h_normalised_mixed)
                h_corrs.append(h_corr)

            __, __, h_corr = correlation_function_centrality_integrated(h_sames[:3], h_normalised_mixeds[:3], outdir_suffix, mode, '050', suffix)
            fit_normalisation_region(h_corr, outdir_suffix, mode, '050', suffix)

            __, __, h_corr = direct_computation_correlation_function_centrality_integrated(h_sames[:3], h_mixeds[:3], outdir_suffix, mode, '050', suffix)
            __, __, h_corr = direct_computation_correlation_function_centrality_integrated(h_sames[1:3], h_mixeds[1:3], outdir_suffix, mode, '1050', suffix)
            __, __, h_corr = direct_computation_correlation_function_centrality_integrated(h_sames, h_mixeds, outdir_suffix, mode, '080', suffix)
            __, __, h_corr = direct_computation_correlation_function_centrality_integrated(h_sames[1:], h_mixeds[1:], outdir_suffix, mode, '1080', suffix)

            h_corr = wighted_sum_correlation_function_centrality_integrated(h_corrs, outdir_suffix, mode, '050', suffix)
            
            h_sames.clear()
            h_normalised_mixeds.clear()
            h_mixeds.clear()
            h_corrs.clear()

    outfile.Close()
