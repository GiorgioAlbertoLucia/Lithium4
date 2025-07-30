
from ROOT import TFile, TH1F, TDirectory

NORM_LOW_KSTAR = 0.2 # 0.25
NORM_HIGH_KSTAR = 0.4 # 0.75
NBINS_KSTAR = 40 # 20

def normalise_histograms_and_compute_correlation(outdir:TDirectory, mode:str, centrality:str):

    h_same = file_same.Get(f'kstar{mode}/hKstar{centrality}{mode}')
    h_same.SetDirectory(0)  # Detach from file to avoid issues with deletion
    h_same.SetName(f'hSameEvent{centrality}')
    #h_same.Rebin()

    h_mixed = file_mixed.Get(f'kstar{mode}/hKstar{centrality}{mode}')
    h_mixed.SetDirectory(0)  # Detach from file to avoid issues with deletion
    h_mixed.SetName(f'hMixedEvent{centrality}')
    #h_mixed.Rebin()


    low_bin = h_same.FindBin(NORM_LOW_KSTAR)
    high_bin = h_same.FindBin(NORM_HIGH_KSTAR)
    normalization_factor = h_same.Integral(low_bin, high_bin) / h_mixed.Integral(low_bin, high_bin)
    h_mixed.Scale(normalization_factor)

    h_corr = h_same.Clone(f'hCorrelation{centrality}')
    h_corr.Divide(h_mixed)

    outdir.cd()
    h_same.Write()
    h_mixed.Write()
    h_corr.Write()

    return h_same, h_mixed

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
    h_same.Write()
    h_mixed.Write()
    h_corr.Write()
        

if __name__ == '__main__':
    
    infile_same = '/home/galucia/antiLithium4/root_dataframe/output/same_event.root'
    infile_mixed = '/home/galucia/antiLithium4/root_dataframe/output/mixed_event.root'
    outfile = TFile.Open('/home/galucia/antiLithium4/root_dataframe/output/correlation.root', 'RECREATE')   

    file_same = TFile.Open(infile_same)
    file_mixed = TFile.Open(infile_mixed)

    for mode in ['', 'Matter', 'Antimatter']:
    
        outdir = outfile.mkdir(f'Correlation{mode}')

        h_sames, h_mixeds = [], []

        for centrality in ['010', '1030', '3050']:
            h_same, h_mixed = normalise_histograms_and_compute_correlation(outdir, mode, centrality)
            h_sames.append(h_same)
            h_mixeds.append(h_mixed)
        
        correlation_function_centrality_integrated(h_sames, h_mixeds, outdir, mode)

    outfile.Close()
