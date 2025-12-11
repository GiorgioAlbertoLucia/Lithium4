
from ROOT import TFile, TH1F, TDirectory, TCanvas, TDirectory, TF1
from torchic.utils.root import set_root_object

NORM_LOW_KSTAR = 0.2 # 0.25
NORM_HIGH_KSTAR = 0.4 # 0.8
NBINS_KSTAR = 40 # 20

def normalise_histograms_and_compute_correlation(infile_same, infile_mixed, outdir:TDirectory,
                                                 mode:str, centrality:str, rebin:int=1, suffix:str=''):

    h_same = infile_same.Get(f'kstar{mode}/hKstar{centrality}{suffix}{mode}')
    h_same.SetDirectory(0)  # Detach from file to avoid issues with deletion
    h_same.SetName(f'hSameEvent{centrality}{suffix}')
    if rebin > 1:
        h_same.Rebin(rebin)

    h_mixed = infile_mixed.Get(f'kstar{mode}/hKstar{centrality}{suffix}{mode}')
    h_mixed.SetDirectory(0)  # Detach from file to avoid issues with deletion
    h_mixed.SetName(f'hMixedEvent{centrality}{suffix}')
    if rebin > 1:
        h_mixed.Rebin(rebin)


    low_bin = h_same.FindBin(NORM_LOW_KSTAR)
    high_bin = h_same.FindBin(NORM_HIGH_KSTAR)
    normalization_factor = h_same.Integral(low_bin, high_bin) / h_mixed.Integral(low_bin, high_bin)
    h_mixed.Scale(normalization_factor)

    h_corr = h_same.Clone(f'hCorrelation{centrality}{suffix}')
    h_corr.Divide(h_mixed)

    if suffix != '':
        return h_same, h_mixed, h_corr
        
    outdir.cd()
    for hist in [h_same, h_mixed, h_corr]:
        hist.SetTitle(';#it{k}* (GeV/#it{c});')
        hist.Write()

    return h_same, h_mixed, h_corr

def fit_normalisation_region(h_corr:TH1F, outdir:TDirectory, mode:str, centrality:str, suffix:str=''):

    constant = TF1(f'func{centrality}', 'pol0', NORM_LOW_KSTAR, NORM_HIGH_KSTAR)
    h_corr.Fit(constant, 'RMS+')

    canvas = TCanvas(f'cCorrelation{centrality}{suffix}')
    set_root_object(h_corr, marker_style=20)
    h_corr.Draw('hist pe1')
    constant.Draw('same')
    outdir.cd()
    canvas.Write()


def correlation_function_centrality_integrated(h_sames, h_mixeds, outdir:TDirectory,
                                               mode:str, centrality:str = '050', suffix:str=''):
    """
    Compute the correlation function for the centrality integrated histograms.
    """
    h_same = h_sames[0].Clone(f'hSameEvent{centrality}{suffix}')
    h_mixed = h_mixeds[0].Clone(f'hMixedEvent{centrality}{suffix}')

    for h_same_centrality, h_mixed_centrality in zip(h_sames[1:], h_mixeds[1:]):
        h_same.Add(h_same_centrality)
        h_mixed.Add(h_mixed_centrality)

    h_corr = h_same.Clone(f'hCorrelation{centrality}{suffix}')
    h_corr.Divide(h_mixed)

    outdir.cd()
    for hist in [h_same, h_mixed, h_corr]:
        hist.SetTitle(';#it{k}* (GeV/#it{c});')
        hist.Write()

    return h_same, h_mixed, h_corr
        

if __name__ == '__main__':
    
    #infile_same_path = 'output/same_event.root'
    #infile_mixed_path = 'output/mixed_event.root'
    #outfile_path = 'output/correlation.root'

    infile_same_path = 'checks/same_event_hadronpid_pass1_pass4_nohe3pcut_offlinetpc.root'
    infile_mixed_path = 'checks/mixed_event_hadronpid_pass1_pass4_nohe3pcut_offlinetpc.root'
    outfile_path = 'checks/correlation_hadronpid_pass1_pass4_nohe3pcut_offlinetpc.root'

    infile_same = TFile.Open(infile_same_path)
    infile_mixed = TFile.Open(infile_mixed_path)
    outfile = TFile.Open(outfile_path, 'RECREATE')   

    for mode in ['', 'Matter', 'Antimatter']:
    
        outdir = outfile.mkdir(f'Correlation{mode}')

        h_sames, h_mixeds = [], []

        for suffix in ['', 'SharedUnder50', 'SharedUnder40', 'SharedUnder30', 'PLi2', 'PLi3', 'PLi4', 'PLi5',
                       'PLiUnder2', 'PLiUnder3', 'PLiUnder4', 'PLiUnder5', 'PtHadronUnder1p0', 'PtHadronUnder1p1', 'PtHadronUnder1p2']:
            for centrality in ['010', '1030', '3050']:
                h_same, h_mixed, h_corr = normalise_histograms_and_compute_correlation(infile_same, infile_mixed,
                                                                                       outdir, mode, centrality, 
                                                                                       rebin=2, suffix=suffix)
                if suffix == '':
                    fit_normalisation_region(h_corr, outdir, mode, centrality, suffix)
                h_sames.append(h_same)
                h_mixeds.append(h_mixed)

            __, __, h_corr = correlation_function_centrality_integrated(h_sames, h_mixeds, outdir, mode, '050', suffix)
            fit_normalisation_region(h_corr, outdir, mode, '050', suffix)
            
            h_sames.clear()
            h_mixeds.clear()

    outfile.Close()
