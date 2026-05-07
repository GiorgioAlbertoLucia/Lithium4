
from ROOT import TCanvas, TLegend, TFile, TPaveText

from torchic.core.histogram import load_hist
from torchic.utils.root import set_root_object
from torchic.utils.colors import get_color

def perform_kolmogorov_smirnov_test(hist1, hist2, outfile, canvas_name, particle):

    kolmogorov_variable = hist1.KolmogorovTest(hist2)

    canvas = TCanvas(canvas_name, "", 800, 600)
    hframe = canvas.DrawFrame(0.0, 0.0, 10.0, max(hist1.GetMaximum(), hist2.GetMaximum())*1.2, 
                              f"Kolmogorov-Smirnov test: {kolmogorov_variable:.4f} ;#it{{p}}_{{T}} (GeV/#it{{c}});Counts")
    if particle == 'Pr':
        set_root_object(hist1, line_color=797, line_width=2, title="p")
        set_root_object(hist2, line_color=632, line_width=2, title="#bar{p}")
    else:
        set_root_object(hist1, line_color=797, line_width=2, title="^{3}He")
        set_root_object(hist2, line_color=632, line_width=2, title="^{3}#bar{He}")
    
    hist1.Draw("same hist e0")
    hist2.Draw("same hist e0")

    legend = TLegend(0.6, 0.7, 0.88, 0.88)
    legend.SetBorderSize(0)
    if particle == 'Pr':
        legend.AddEntry(hist1, "p", "l")
        legend.AddEntry(hist2, "#bar{p}", "l")
    else:
        legend.AddEntry(hist1, "^{3}He", "l")
        legend.AddEntry(hist2, "^{3}#bar{He}", "l")
    legend.Draw()

    text = TPaveText(0.2, 0.75, 0.5, 0.9)
    text.SetFillColor(0)
    text.SetBorderSize(0)
    text.SetTextAlign(12)
    text.SetTextSize(0.04)
    text.AddText(f"KS test: {kolmogorov_variable:.4f}")
    text.Draw('same')

    outfile.cd()
    canvas.Write()

def projections(infile, outdir, particle):

    particle_suffix = 'Had' if particle == 'Pr' else 'He3'
    hist_name = f"QA/hPt{particle_suffix}Kstar"
    print(f"[INFO] Loading histogram: {hist_name} from file: {infile}")
    h2 = load_hist(infile, hist_name)
    kstar_min, kstarmax, kstar_step = 0.1, 0.18, 0.02

    outdir.cd()
    h2.Write()

    dummy_canvas = TCanvas("dummy_canvas", "", 800, 600)
    dummy_canvas.Print(f"output/hPt{particle_suffix}Kstar.pdf[")  # Open PDF for writing

    for kstar_low in [0.1, 0.12, 0.14, 0.16, 0.18]:
        kstar_high = kstar_low + kstar_step
        kstar_low_bin = h2.GetYaxis().FindBin(kstar_low)
        kstar_high_bin = h2.GetYaxis().FindBin(kstar_high) - 1
        h_pt_proj = h2.ProjectionX(
            f"hPt{particle}Kstar_{kstar_low:.2f}_{kstar_high:.2f}",
            kstar_low_bin, kstar_high_bin
        )
        h_pt_proj.Rebin(10)

        h_pt_proj_matter = h_pt_proj.Clone(f"hPt{particle}Kstar_Matter_{kstar_low:.2f}_{kstar_high:.2f}")
        h_pt_proj_antimatter = h_pt_proj.Clone(f"hPt{particle}Kstar_Antimatter_{kstar_low:.2f}_{kstar_high:.2f}")

        for ibin in range(1, h_pt_proj.GetNbinsX() + 1):
            pt = h_pt_proj.GetXaxis().GetBinCenter(ibin)
            if pt < 0.:
                continue
            
            negative_pt_bin = h_pt_proj_antimatter.FindBin(-pt)
            h_pt_proj_antimatter.SetBinContent(ibin, h_pt_proj.GetBinContent(negative_pt_bin))
            h_pt_proj_antimatter.SetBinError(ibin, h_pt_proj.GetBinError(negative_pt_bin))
        
        zero_bin = h_pt_proj_matter.FindBin(0.)
        h_pt_proj_matter.Scale(1./h_pt_proj_matter.Integral(zero_bin, h_pt_proj_matter.GetNbinsX()))
        h_pt_proj_antimatter.Scale(1./h_pt_proj_antimatter.Integral(zero_bin, h_pt_proj_antimatter.GetNbinsX()))

        perform_kolmogorov_smirnov_test(h_pt_proj_matter, h_pt_proj_antimatter, outdir, 
                                        f"c_KS_{particle_suffix}_{kstar_low:.2f}_{kstar_high:.2f}", 
                                        particle)

        canvas = TCanvas(f"c_{particle_suffix}_{kstar_low:.2f}_{kstar_high:.2f}", "", 800, 600)
        hframe = canvas.DrawFrame(0.0, 0.0, 10.0 if particle == 'He' else 4., max(h_pt_proj_matter.GetMaximum(), h_pt_proj_antimatter.GetMaximum())*1.2,
                                  f"{particle}, {kstar_low:.2f} < #it{{k}}* < {kstar_high:.2f} GeV/#it{{c}} ;#it{{p}}_{{T}} (GeV/#it{{c}});Counts")
        set_root_object(h_pt_proj_matter, line_color=get_color(0), line_width=2, title="Matter")
        set_root_object(h_pt_proj_antimatter, line_color=get_color(1), line_width=2, title="Antimatter")
        h_pt_proj_matter.Draw("same hist e0")
        h_pt_proj_antimatter.Draw("same hist e0")

        legend = TLegend(0.6, 0.7, 0.88, 0.88)
        legend.SetBorderSize(0)
        legend.AddEntry(h_pt_proj_matter, "Matter", "l")
        legend.AddEntry(h_pt_proj_antimatter, "Antimatter", "l")
        legend.Draw()

        outdir.cd()
        canvas.Write()
        canvas.Print(f"output/hPt{particle_suffix}Kstar.pdf")  # Save page to PDF

    dummy_canvas.Print(f"output/hPt{particle_suffix}Kstar.pdf]")  # Close PDF

    canvas = TCanvas(f"c_{particle_suffix}_matter", "", 800, 600)
    hs_pt_proj_matter = []

    for kstar_low in [0.1, 0.12, 0.14, 0.16]:
        kstar_high = kstar_low + kstar_step
        kstar_low_bin = h2.GetYaxis().FindBin(kstar_low)
        kstar_high_bin = h2.GetYaxis().FindBin(kstar_high) - 1
        h_pt_proj_matter = h2.ProjectionX(
            f"hPt{particle}KstarMatter_{kstar_low:.2f}_{kstar_high:.2f}", kstar_low_bin, kstar_high_bin)
        h_pt_proj_matter.SetTitle(f"{kstar_low:.2f} < #it{{k}}* < {kstar_high:.2f} GeV/#it{{c}}")
        hs_pt_proj_matter.append(h_pt_proj_matter)
    
    legend = TLegend(0.6, 0.7, 0.88, 0.88)
    legend.SetBorderSize(0)

    maximum = max(h_pt_proj_matter.GetMaximum() for h_pt_proj_matter in hs_pt_proj_matter)
    hframe = canvas.DrawFrame(0.0, 0.0, 10.0 if particle == 'He' else 4., maximum*1.2, f"{particle}, matter ;#it{{p}}_{{T}} (GeV/#it{{c}});Counts")
    for item_index, h_pt_proj_matter in enumerate(hs_pt_proj_matter):
        set_root_object(h_pt_proj_matter, line_width=2, line_color=get_color(item_index))
        h_pt_proj_matter.Draw("same hist e0")
        legend.AddEntry(h_pt_proj_matter, h_pt_proj_matter.GetTitle(), "l")
    legend.Draw()

    outdir.cd()
    canvas.Write()

    canvas = TCanvas(f"c_{particle_suffix}_antimatter", "", 800, 600)
    hs_pt_proj_antimatter = []

    for kstar_low in [0.1, 0.12, 0.14, 0.16, 0.18]:
        kstar_high = kstar_low + kstar_step
        kstar_low_bin = h2.GetYaxis().FindBin(kstar_low)
        kstar_high_bin = h2.GetYaxis().FindBin(kstar_high) - 1
        h_pt_proj = h2.ProjectionX(
            f"hPt{particle}Kstar_{kstar_low:.2f}_{kstar_high:.2f}", kstar_low_bin, kstar_high_bin)
        h_pt_proj_antimatter = h_pt_proj.Clone(f"hPt{particle}Kstar_Antimatter_{kstar_low:.2f}_{kstar_high:.2f}")
        h_pt_proj_antimatter.SetTitle(f"{kstar_low:.2f} < #it{{k}}* < {kstar_high:.2f} GeV/#it{{c}}")
        
        for ibin in range(1, h_pt_proj_antimatter.GetNbinsX() + 1):
            pt = h_pt_proj_antimatter.GetXaxis().GetBinCenter(ibin)
            if pt < 0.:
                continue

            negative_pt_bin = h_pt_proj.FindBin(-pt)
            h_pt_proj_antimatter.SetBinContent(ibin, h_pt_proj.GetBinContent(negative_pt_bin))
            h_pt_proj_antimatter.SetBinError(ibin, h_pt_proj.GetBinError(negative_pt_bin))
            
            negative_pt_bin = h_pt_proj_antimatter.FindBin(-pt)
            h_pt_proj_antimatter.SetBinContent(ibin, h_pt_proj_antimatter.GetBinContent(negative_pt_bin))
            h_pt_proj_antimatter.SetBinError(ibin, h_pt_proj_antimatter.GetBinError(negative_pt_bin))
        
        hs_pt_proj_antimatter.append(h_pt_proj_antimatter)
    
    legend = TLegend(0.6, 0.7, 0.88, 0.88)
    legend.SetBorderSize(0)

    maximum = max(h_pt_proj_antimatter.GetMaximum() for h_pt_proj_antimatter in hs_pt_proj_antimatter)
    hframe = canvas.DrawFrame(0.0, 0.0, 10.0 if particle == 'He' else 4., maximum*1.2, f"{particle}, antimatter ;#it{{p}}_{{T}} (GeV/#it{{c}});Counts")
    for item_index, h_pt_proj_antimatter in enumerate(hs_pt_proj_antimatter):
        set_root_object(h_pt_proj_antimatter, line_width=2, line_color=get_color(item_index))
        h_pt_proj_antimatter.Draw("same hist e0")
        legend.AddEntry(h_pt_proj_antimatter, h_pt_proj_antimatter.GetTitle(), "l")
    legend.Draw()

    outdir.cd()
    canvas.Write()

    hs_ratios = []
    for h_pt_proj_matter, h_pt_proj_antimatter in zip(hs_pt_proj_matter, hs_pt_proj_antimatter):

        h_ratio = h_pt_proj_antimatter.Clone(f"hPt{particle}Kstar_Ratio_{h_pt_proj_matter.GetTitle().replace(' ', '_')}")
        h_ratio.SetTitle(h_pt_proj_matter.GetTitle())

        #for ibin in range(1, h_ratio.GetNbinsX() + 1):
        #    pt = h_ratio.GetXaxis().GetBinCenter(ibin)
        #    if pt < 0.:
        #        continue
        #    negative_pt_bin = h_ratio.FindBin(-pt)
        #    h_ratio.SetBinContent(ibin, h_ratio.GetBinContent(negative_pt_bin))
        #    h_ratio.SetBinError(ibin, h_ratio.GetBinError(negative_pt_bin))
        
        h_ratio.Divide(h_pt_proj_matter)

        hs_ratios.append(h_ratio)
    
    canvas = TCanvas(f"c_{particle_suffix}_ratio", "", 800, 600)
    legend = TLegend(0.6, 0.7, 0.88, 0.88)
    legend.SetBorderSize(0)
    hframe = canvas.DrawFrame(0.0, 0.0, 10.0 if particle == 'He' else 4., 2., f"{particle}, antimatter/matter ratio ;#it{{p}}_{{T}} (GeV/#it{{c}});Counts")
    for item_index, h_ratio in enumerate(hs_ratios):
        set_root_object(h_ratio, line_width=2, line_color=get_color(item_index), 
                        marker_style=20+item_index, marker_color=get_color(item_index))
        canvas.cd()
        h_ratio.Draw("same p e0")
        legend.AddEntry(h_ratio, h_ratio.GetTitle(), "l")

        outdir.cd()
        h_ratio.Write()
        
    legend.Draw()
    outdir.cd()
    canvas.Write()

def compare_pr_he_projections(outdir_event):

    for sign in ['Matter', 'Antimatter']:
        for kstar_low in [0.1, 0.12, 0.14, 0.16]:
            canvas = TCanvas(f"c_compare_pr_he_{sign}_{kstar_low}_{kstar_low+0.02}", "", 800, 600)

            h2_pt_vs_kstar_pr = outdir_event.Get(f"Pr/hPtHadKstar")
            h2_pt_vs_kstar_he = outdir_event.Get(f"He/hPtHe3Kstar")

            kstar_low_bin_pr = h2_pt_vs_kstar_pr.GetYaxis().FindBin(kstar_low)
            kstar_high_bin_pr = h2_pt_vs_kstar_pr.GetYaxis().FindBin(kstar_low + 0.02) - 1
            h_pt_proj_pr = h2_pt_vs_kstar_pr.ProjectionX(f"hPtPrKstar_{kstar_low:.2f}_{kstar_low+0.02:.2f}", kstar_low_bin_pr, kstar_high_bin_pr)
            
            kstar_low_bin_he = h2_pt_vs_kstar_he.GetYaxis().FindBin(kstar_low)
            kstar_high_bin_he = h2_pt_vs_kstar_he.GetYaxis().FindBin(kstar_low + 0.02) - 1
            h_pt_proj_he = h2_pt_vs_kstar_he.ProjectionX(f"hPtHeKstar_{kstar_low:.2f}_{kstar_low+0.02:.2f}", kstar_low_bin_he, kstar_high_bin_he)

            if sign == 'Antimatter':
                for ibin in range(1, h_pt_proj_pr.GetNbinsX() + 1):
                    pt = h_pt_proj_pr.GetXaxis().GetBinCenter(ibin)
                    if pt < 0.:
                        continue
                    
                    negative_pt_bin = h_pt_proj_pr.FindBin(-pt)
                    h_pt_proj_pr.SetBinContent(ibin, h_pt_proj_pr.GetBinContent(negative_pt_bin))
                    h_pt_proj_pr.SetBinError(ibin, h_pt_proj_pr.GetBinError(negative_pt_bin))
                    h_pt_proj_he.SetBinContent(ibin, h_pt_proj_he.GetBinContent(negative_pt_bin))
                    h_pt_proj_he.SetBinError(ibin, h_pt_proj_he.GetBinError(negative_pt_bin))

            set_root_object(h_pt_proj_pr, line_color=797, line_width=2, title="Protons")
            set_root_object(h_pt_proj_he, line_color=418, line_width=2, title="Helium-3")

            maximum = max(h_pt_proj_pr.GetMaximum(), h_pt_proj_he.GetMaximum())
            hframe = canvas.DrawFrame(0.0, 0.0, 10., maximum*1.2, f"{sign} comparison for {kstar_low} < #it{{k}}* < {kstar_low+0.02} GeV/#it{{c}} ;#it{{p}}_{{T}} (GeV/#it{{c}});Counts")
            h_pt_proj_pr.Draw("same hist e0")
            h_pt_proj_he.Draw("same hist e0")

            legend = TLegend(0.6, 0.7, 0.88, 0.88)
            legend.SetBorderSize(0)
            if sign == 'Antimatter':
                legend.AddEntry(h_pt_proj_pr, "#bar{p}", "l")
                legend.AddEntry(h_pt_proj_he, "^{3}#bar{He}", "l")
            if sign == 'Matter':
                legend.AddEntry(h_pt_proj_pr, "p", "l")
                legend.AddEntry(h_pt_proj_he, "^{3}He", "l")
            legend.Draw()

            outdir_event.cd()
            canvas.Write()

    

if __name__ == '__main__':
    
    infiles = {'SameEvent': '/home/galucia/Lithium4/preparation/checks/same_event_hadronpid_pass1_pass4_refined_dca.root',
               'MixedEvent': '/home/galucia/Lithium4/preparation/checks/mixed_event_hadronpid_pass1_pass4_refined_dca.root'
               }
    outfile = TFile.Open('output/pt_vs_kstar_projections.root', "recreate")

    for event_type, infile in infiles.items():
        
        outdir_event = outfile.mkdir(event_type)
        
        for particle in ['Pr', 'He']:
            outdir = outdir_event.mkdir(particle)
            projections(infile, outdir, particle)

        compare_pr_he_projections(outdir_event)
