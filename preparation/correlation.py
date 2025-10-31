
from ROOT import TFile, TH1F, TDirectory, TCanvas, TDirectory, TF1
from torchic.utils.root import set_root_object

NORM_LOW_KSTAR = 0.2 # 0.25
NORM_HIGH_KSTAR = 0.4 # 0.8
NBINS_KSTAR = 40 # 20

def normalise_histograms_and_compute_correlation(outdir:TDirectory, mode:str, centrality:str, rebin:int=1):

    h_same = file_same.Get(f'kstar{mode}/hKstar{centrality}{mode}')
    h_same.SetDirectory(0)  # Detach from file to avoid issues with deletion
    h_same.SetName(f'hSameEvent{centrality}')
    if rebin > 1:
        h_same.Rebin(rebin)

    h_mixed = file_mixed.Get(f'kstar{mode}/hKstar{centrality}{mode}')
    h_mixed.SetDirectory(0)  # Detach from file to avoid issues with deletion
    h_mixed.SetName(f'hMixedEvent{centrality}')
    if rebin > 1:
        h_mixed.Rebin(rebin)


    low_bin = h_same.FindBin(NORM_LOW_KSTAR)
    high_bin = h_same.FindBin(NORM_HIGH_KSTAR)
    normalization_factor = h_same.Integral(low_bin, high_bin) / h_mixed.Integral(low_bin, high_bin)
    h_mixed.Scale(normalization_factor)

    h_corr = h_same.Clone(f'hCorrelation{centrality}')
    h_corr.Divide(h_mixed)

    outdir.cd()
    for hist in [h_same, h_mixed, h_corr]:
        hist.SetTitle(';#it{k}* (GeV/#it{c});')
        hist.Write()

    return h_same, h_mixed, h_corr

def fit_normalisation_region(h_corr:TH1F, outdir:TDirectory, mode:str, centrality:str):

    constant = TF1(f'func{centrality}', 'pol0', NORM_LOW_KSTAR, NORM_HIGH_KSTAR)
    h_corr.Fit(constant, 'RMS+')

    canvas = TCanvas(f'cCorrelation{centrality}')
    set_root_object(h_corr, marker_style=20)
    h_corr.Draw('hist pe1')
    constant.Draw('same')
    outdir.cd()
    canvas.Write()


def correlation_function_centrality_integrated(h_sames, h_mixeds, outdir:TDirectory, mode:str, centrality:str = '050'):
    """
    Compute the correlation function for the centrality integrated histograms.
    """
    h_same = h_sames[0].Clone(f'hSameEvent{centrality}')
    h_mixed = h_mixeds[0].Clone(f'hMixedEvent{centrality}')

    for h_same_centrality, h_mixed_centrality in zip(h_sames[1:], h_mixeds[1:]):
        h_same.Add(h_same_centrality)
        h_mixed.Add(h_mixed_centrality)

    h_corr = h_same.Clone(f'hCorrelation{centrality}')
    h_corr.Divide(h_mixed)

    outdir.cd()
    for hist in [h_same, h_mixed, h_corr]:
        hist.SetTitle(';#it{k}* (GeV/#it{c});')
        hist.Write()

    return h_same, h_mixed, h_corr
        

if __name__ == '__main__':
    
    #infile_same = 'output/same_event.root'
    #infile_mixed = 'output/mixed_event.root'
    #_outfile = 'output/correlation.root'

    infile_same = 'checks/same_event_hadronpid_pass1_pass4.root'
    infile_mixed = 'checks/mixed_event_hadronpid_pass1_pass4.root'
    _outfile = 'checks/correlation_hadronpid_pass1_pass4.root'

    file_same = TFile.Open(infile_same)
    file_mixed = TFile.Open(infile_mixed)
    outfile = TFile.Open(_outfile, 'RECREATE')   

    for mode in ['', 'Matter', 'Antimatter']:
    
        outdir = outfile.mkdir(f'Correlation{mode}')

        h_sames, h_mixeds = [], []

        for centrality in ['010', '1030', '3050']:
            h_same, h_mixed, h_corr = normalise_histograms_and_compute_correlation(outdir, mode, centrality, rebin=2)
            fit_normalisation_region(h_corr, outdir, mode, centrality)
            h_sames.append(h_same)
            h_mixeds.append(h_mixed)
        
        __, __, h_corr = correlation_function_centrality_integrated(h_sames, h_mixeds, outdir, mode)
        fit_normalisation_region(h_corr, outdir, mode, '050')

    outfile.Close()
