from ROOT import TFile, TCanvas, TPad, TPaveText

from torchic.core.histogram import load_hist
from torchic.utils.root import set_root_object, init_legend, set_alice_global_style, set_alice_frame_style
from torchic.utils.colors import get_color

def get_hist_maximum(hist, x_range):
    x_min, x_max = x_range
    bin_min = hist.GetXaxis().FindBin(x_min)
    bin_max = hist.GetXaxis().FindBin(x_max)
    max_value = 0
    for bin in range(bin_min, bin_max + 1):
        value = hist.GetBinContent(bin)
        if value > max_value:
            max_value = value
    return max_value

def get_hist_minimum(hist, x_range):
    x_min, x_max = x_range
    bin_min = hist.GetXaxis().FindBin(x_min)
    bin_max = hist.GetXaxis().FindBin(x_max)
    min_value = float('inf')
    for bin in range(bin_min, bin_max + 1):
        value = hist.GetBinContent(bin)
        if value < min_value:
            min_value = value
    return min_value

def get_hist_maximum_with_error(hist, x_range):
    x_min, x_max = x_range
    bin_min = hist.GetXaxis().FindBin(x_min)
    bin_max = hist.GetXaxis().FindBin(x_max)
    max_value = 0
    for bin in range(bin_min, bin_max + 1):
        value = hist.GetBinContent(bin) + hist.GetBinError(bin)
        if value > max_value:
            max_value = value
    return max_value

def get_hist_minimum_with_error(hist, x_range):
    x_min, x_max = x_range
    bin_min = hist.GetXaxis().FindBin(x_min)
    bin_max = hist.GetXaxis().FindBin(x_max)
    min_value = float('inf')
    for bin in range(bin_min, bin_max + 1):
        value = hist.GetBinContent(bin) - hist.GetBinError(bin)
        if value < min_value:
            min_value = value
    return min_value

def get_kolmogorov_score(hist1, hist2):
    # Ensure both histograms have the same number of bins and range
    if hist1.GetNbinsX() != hist2.GetNbinsX():
        raise ValueError("Histograms must have the same number of bins.")
    
    hist1_normalized = hist1.Clone("hist1_normalized")
    hist2_normalized = hist2.Clone("hist2_normalized")
    hist1_normalized.Scale(1.0 / hist1_normalized.Integral() if hist1_normalized.Integral() != 0 else 1.0)
    hist2_normalized.Scale(1.0 / hist2_normalized.Integral() if hist2_normalized.Integral() != 0 else 1.0)

    kolmogorov_score = hist1_normalized.KolmogorovTest(hist2_normalized)
    return kolmogorov_score

if __name__ == "__main__":
    
    outfile = TFile("compare_event_selections.root", "RECREATE")
    set_alice_global_style()
    
    hist_names = ['he3-hadron-femto/QA/hVtxZ',
                  'he3-hadron-femto/QA/hCentralityFT0C',
                  'he3-hadron-femto/QA/hhe3HadInvMass',
                  'he3-hadron-femto/QA/hhe3HadKstar',
                  'he3-hadron-femto/QA/He3/hHe3Pt',
                  'he3-hadron-femto/QA/He3/hDCAxyHe3',
                  'he3-hadron-femto/QA/He3/hDCAzHe3',
                  'he3-hadron-femto/QA/Had/hHadronPt',
                  'he3-hadron-femto/QA/Had/hDCAxyHad',
                  'he3-hadron-femto/QA/Had/hDCAzHad',]
    x_ranges = [(-11, 11),
                (0, 100),
                (3.74, 3.79),
                (0.033, 0.8),
                (-6, 6),
                (-0.5, 0.5),
                (-1, 1),
                (-3, 3),
                (-0.5, 0.5),
                (-1, 1),]
    rebin = [1, 1, 1, 6, 1, 1, 1, 1, 1, 1]
    
    for year in [
                 '23_PbPb_pass5', 
                 '24ar_pass3', 
                 '25_PbPb_pass1'
                ]:
        
        outdir = outfile.mkdir(year)
        
        for ihist, hist_name in enumerate(hist_names):
            x_range = x_ranges[ihist]
            rebin_factor = rebin[ihist]
            canvas = TCanvas(f"canvas_{year}_{hist_name}", f"canvas_{year}_{hist_name}")
            
            h_default = load_hist(f"/data/galucia/lithium/event_selections/LHC{year}_default_AnalysisResults.root", hist_name)
            set_root_object(h_default, line_color=get_color(0), rebin=rebin_factor)
            h_nosamebunchpileup = load_hist(f"/data/galucia/lithium/event_selections/LHC{year}_nosamebunchpileup_AnalysisResults.root", hist_name)
            set_root_object(h_nosamebunchpileup, line_color=get_color(1), rebin=rebin_factor)
            k_score_nosamebunchpileup = get_kolmogorov_score(h_default, h_nosamebunchpileup)
            h_isgoodzvtxft0vspv = load_hist(f"/data/galucia/lithium/event_selections/LHC{year}_isgoodzvtxft0vspv_AnalysisResults.root", hist_name)
            set_root_object(h_isgoodzvtxft0vspv, line_color=get_color(2), rebin=rebin_factor)
            k_score_isgoodzvtxft0vspv = get_kolmogorov_score(h_default, h_isgoodzvtxft0vspv)
            h_nocollisionrofstandard = load_hist(f"/data/galucia/lithium/event_selections/LHC{year}_nocollisionrofstandard_AnalysisResults.root", hist_name)
            set_root_object(h_nocollisionrofstandard, line_color=get_color(3), rebin=rebin_factor)
            k_score_nocollisionrofstandard = get_kolmogorov_score(h_default, h_nocollisionrofstandard)

            h_ratio_nosamebunchpileup = h_nosamebunchpileup.Clone(f"h_ratio_nosamebunchpileup_{year}_{hist_name}")
            h_ratio_nosamebunchpileup.Divide(h_default)
            set_root_object(h_ratio_nosamebunchpileup, line_color=get_color(1))
            h_ratio_isgoodzvtxft0vspv = h_isgoodzvtxft0vspv.Clone(f"h_ratio_isgoodzvtxft0vspv_{year}_{hist_name}")
            h_ratio_isgoodzvtxft0vspv.Divide(h_default)
            set_root_object(h_ratio_isgoodzvtxft0vspv, line_color=get_color(2))
            h_ratio_nocollisionrofstandard = h_nocollisionrofstandard.Clone(f"h_ratio_nocollisionrofstandard_{year}_{hist_name}")
            h_ratio_nocollisionrofstandard.Divide(h_default)
            set_root_object(h_ratio_nocollisionrofstandard, line_color=get_color(3))

            y_max = max(get_hist_maximum(h_default, x_range), get_hist_maximum(h_nosamebunchpileup, x_range), get_hist_maximum(h_isgoodzvtxft0vspv, x_range), get_hist_maximum(h_nocollisionrofstandard, x_range))
            y_max_ratio = max(get_hist_maximum_with_error(h_ratio_nosamebunchpileup, x_range), get_hist_maximum_with_error(h_ratio_isgoodzvtxft0vspv, x_range), get_hist_maximum_with_error(h_ratio_nocollisionrofstandard, x_range))
            y_min_ratio = min(get_hist_minimum_with_error(h_ratio_nosamebunchpileup, x_range), get_hist_minimum_with_error(h_ratio_isgoodzvtxft0vspv, x_range), get_hist_minimum_with_error(h_ratio_nocollisionrofstandard, x_range))

            kolmogorov_text = TPaveText(0.15, 0.75, 0.5, 0.9, "NDC")
            kolmogorov_text.SetFillStyle(0)
            kolmogorov_text.SetBorderSize(0)
            #kolmogorov_text.SetTextSize(0.04)
            #kolmogorov_text.SetTextFont(42)
            kolmogorov_text.AddText(f"Kolmogorov Test")
            kolmogorov_text.AddText(f"No Same Bunch Pileup: {k_score_nosamebunchpileup:.4f}")
            kolmogorov_text.AddText(f"Is Good Z Vtx FT0 vs PV: {k_score_isgoodzvtxft0vspv:.4f}")
            kolmogorov_text.AddText(f"No Collision ROF Standard: {k_score_nocollisionrofstandard:.4f}")
            
            upper_pad = TPad("upper_pad", "upper_pad", 0, 0.3, 1, 1)
            upper_pad.SetBottomMargin(0.0)
            if 'DCA' in hist_name:
                upper_pad.SetLogy()
            
            lower_pad = TPad("lower_pad", "lower_pad", 0, 0, 1, 0.3)
            lower_pad.SetTopMargin(0.0)
            lower_pad.SetBottomMargin(0.3)
            
            canvas.cd()
            upper_pad.Draw()
            lower_pad.Draw()
            
            upper_pad.cd()
            y_min_upper_pad = 0.1 if 'DCA' in hist_name else 0
            upper_pad.DrawFrame(x_range[0], y_min_upper_pad, x_range[1], y_max * 1.2,
                                f';{h_default.GetXaxis().GetTitle()};Entries')
            h_default.Draw("hist same")
            h_nosamebunchpileup.Draw("hist same")
            h_isgoodzvtxft0vspv.Draw("hist same")
            h_nocollisionrofstandard.Draw("hist same")
            legend = init_legend(0.6, 0.7, 0.89, 0.89)
            legend.AddEntry(h_default, "Default", "l")
            legend.AddEntry(h_nosamebunchpileup, "No Same Bunch Pileup", "l")
            legend.AddEntry(h_isgoodzvtxft0vspv, "Is Good Z Vtx FT0 vs PV", "l")
            legend.AddEntry(h_nocollisionrofstandard, "No Collision ROF Standard", "l")
            legend.Draw()
            #kolmogorov_text.Draw("")
            
            lower_pad.cd()
            hframe_lower = lower_pad.DrawFrame(x_range[0], y_min_ratio * 0.8 if y_min_ratio > 0 else y_min_ratio - 0.3 * abs(y_min_ratio), 
                                               x_range[1], y_max_ratio * 1.2,
                                f';{h_default.GetXaxis().GetTitle()};Ratio to default')
            hframe_lower.GetYaxis().SetLabelSize(0.08)
            hframe_lower.GetXaxis().SetLabelSize(0.08)
            hframe_lower.GetYaxis().SetTitleSize(0.08)
            hframe_lower.GetXaxis().SetTitleSize(0.08)
            hframe_lower.GetYaxis().SetTitleOffset(0.5)
            h_ratio_nosamebunchpileup.Draw("hist same")
            h_ratio_isgoodzvtxft0vspv.Draw("hist same")
            h_ratio_nocollisionrofstandard.Draw("hist same")
            
            outdir.cd()
            canvas.Write()
            canvas.SaveAs(f"figures/compare_event_selections_{year}_{hist_name.replace('/', '_')}.pdf")
        