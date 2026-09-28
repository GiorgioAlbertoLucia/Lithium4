from torchic.core.histogram import load_hist
from torchic.core.graph import load_graph
from torchic.utils.root import init_legend, set_root_object, set_alice_global_style
from torchic.utils.colors import get_color

from ROOT import TCanvas, TFile

def rebin_to_match(h1, h2):
    """
    Make the binning of two histograms match. The new values are the average of the original values in the new bins. The new binning is taken from h2.
    """
    original_bin_edges = [h1.GetBinLowEdge(i) for i in range(1, h1.GetNbinsX() + 2)]
    new_bin_edges = [h2.GetBinLowEdge(i) for i in range(1, h2.GetNbinsX() + 2)]
    original_values = [h1.GetBinContent(i) for i in range(1, h1.GetNbinsX() + 1)]
    new_values = []
    h1_rebinned = h2.Clone(f"{h1.GetName()}_rebinned")
    
    for i in range(len(new_bin_edges) - 1):
        new_bin_low = new_bin_edges[i]
        new_bin_high = new_bin_edges[i + 1]
        original_bins_in_range = [j for j in range(len(original_bin_edges) - 1) if original_bin_edges[j] >= new_bin_low and original_bin_edges[j + 1] <= new_bin_high]
        
        if original_bins_in_range:
            avg_value = sum(original_values[j] for j in original_bins_in_range) / len(original_bins_in_range)
            h1_rebinned.SetBinContent(i + 1, avg_value)
        else:
            h1_rebinned.SetBinContent(i + 1, 0)
    
    return h1_rebinned

if __name__ == "__main__":
    
    
    set_alice_global_style()
    outfile = TFile.Open('output/model_comparison.root', 'recreate')
    
    # compare nominals
    radii = {'LL': [4.33, 3.46], #fm
             'SW': [6.12, 4.90], #fm
             }
    
    LL_PATH = '/home/galucia/PhaseShiftAnalysis/numerical_lednicky/pHe/output/pHe3_LL_bands.root'
    SW_PATH = '/home/galucia/phemto/output/pHe3_square_well_bands.root'
    COULOMB_PATH = '/home/galucia/phemto/output/pHe3_coulomb.root'
    
    for iradius, radius in enumerate(radii['SW']):
        g_LL = load_graph(LL_PATH, f"Rs{radii['LL'][iradius]:.2f}/g{radii['LL'][iradius]:.2f}_band")
        set_root_object(g_LL, title=f'#it{{R}}_{{s}} = {radii["SW"][iradius]:.2f} fm;'
                                     f'#it{{k}}* (MeV/#it{{c}}); #it{{C}}(#it{{k}}*)', 
                        line_color=get_color(0), line_width=2,
                        fill_color_alpha=(get_color(0), 0.5), fill_style=1001)
        
        g_SW = load_graph(SW_PATH, f"pHe3_square_well_band_r={radii['SW'][iradius]:.2f}_fm")
        set_root_object(g_SW, line_color=get_color(1), line_width=2,
                        fill_color_alpha=(get_color(1), 0.5), fill_style=1001)
        
        h_coulomb = load_hist(COULOMB_PATH, f"r={radii['SW'][iradius]:.2f}_fm/hcats_CF")
        set_root_object(h_coulomb, line_color=get_color(2), line_width=2, line_style=1)
        
        legend = init_legend(0.4, 0.2, 0.8, 0.4)
        legend.AddEntry(g_LL, f"Lednicky-Lyuboshits", "l")
        legend.AddEntry(g_SW, f"Square-Well", "l")
        legend.AddEntry(h_coulomb, f"Coulomb", "l")
        
        canvas = TCanvas(f"r={radii['SW'][iradius]:.2f}_fm", "canvas", 800, 600)
        h_frame = canvas.DrawFrame(0, 0, 400, 1.2, f'#it{{R}}_{{s}} = {radii["SW"][iradius]:.2f} fm;'
                                                   f'#it{{k}}* (MeV/#it{{c}}); #it{{C}}(#it{{k}}*)')
        g_LL.Draw('l3 same')
        g_SW.Draw('l3 same')
        h_coulomb.Draw('hist same')
        legend.Draw('same')
        canvas.SaveAs(f"output/model_comparison_r={radii['SW'][iradius]:.2f}_fm.pdf")
        
        outfile.cd()
        canvas.Write()
    
    outfile.Close()