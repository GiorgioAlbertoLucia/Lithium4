import ROOT

def draw_protons():
    # Open files
    f0 = ROOT.TFile("/home/galucia/Lithium4/preparation/output/single_track_efficiency.root")
    f1 = ROOT.TFile("/home/galucia/Lithium4/preparation/output/single_track_efficiency_no_tof.root")
    f2 = ROOT.TFile("/home/galucia/Lithium4/preparation/output/single_track_efficiency_tof_hit.root")

    # Get histograms
    hist       = f0.Get("hEfficiencyHadAntimatter")
    hist_no_tof  = f1.Get("hEfficiencyHadAntimatter")
    hist_tof_hit = f2.Get("hEfficiencyHadAntimatter")

    hist.SetName("2sigma")

    # Styling
    hist.SetLineColor(ROOT.kOrange - 3)
    hist_no_tof.SetLineColor(ROOT.kAzure + 3)
    hist_tof_hit.SetLineColor(ROOT.kRed + 2)

    for h in [hist, hist_no_tof, hist_tof_hit]:
        h.SetLineWidth(2)

    # Axis title (set on the first drawn histogram)
    hist_no_tof.SetTitle("#bar{p};#it{p}_{T}^{p} (GeV/#it{c}); #varepsilon #times #it{A}")
    hist_no_tof.GetXaxis().SetRangeUser(0, 4)

    # Canvas & draw
    ROOT.gStyle.SetOptStat(0)
    canvas = ROOT.TCanvas()

    hist_no_tof.Draw("hist")
    hist_tof_hit.Draw("hist same")
    hist.Draw("hist same")

    # Legend
    legend = ROOT.TLegend(0.15, 0.68, 0.35, 0.88)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.AddEntry(hist_no_tof,  "No TOF",                "l")
    legend.AddEntry(hist_tof_hit, "TOF signal required",   "l")
    legend.AddEntry(hist,         "|#it{n}#sigma_{TOF}| < 2","l")
    legend.Draw()

    canvas.SaveAs("efficiency_comparison_single_track_protons.pdf")

    f0.Close()
    f1.Close()
    f2.Close()

def draw_li4():
    # Open files
    f0 = ROOT.TFile("/home/galucia/Lithium4/preparation/output/efficiency.root")
    f1 = ROOT.TFile("/home/galucia/Lithium4/preparation/output/efficiency_tof.root")
    f2 = ROOT.TFile("/home/galucia/Lithium4/preparation/output/efficiency_tof_hit.root")

    # Get histograms
    hist       = f0.Get("hEfficiencyAntimatter")
    hist_no_tof  = f1.Get("hEfficiencyAntimatter")
    hist_tof_hit = f2.Get("hEfficiencyAntimatter")

    hist.SetName("2sigma")

    # Styling
    hist.SetLineColor(ROOT.kOrange - 3)
    hist_no_tof.SetLineColor(ROOT.kAzure + 3)
    hist_tof_hit.SetLineColor(ROOT.kRed + 2)

    for h in [hist, hist_no_tof, hist_tof_hit]:
        h.SetLineWidth(2)

    # Axis title (set on the first drawn histogram)
    hist_no_tof.SetTitle("^{4}#bar{Li};#it{p}_{T}^{^{4}Li} (GeV/#it{c}); #varepsilon #times #it{A}")

    # Canvas & draw
    ROOT.gStyle.SetOptStat(0)
    canvas = ROOT.TCanvas()

    hist_no_tof.Draw("hist")
    hist_tof_hit.Draw("hist same")
    hist.Draw("hist same")

    # Legend
    legend = ROOT.TLegend(0.15, 0.68, 0.35, 0.88)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.AddEntry(hist_no_tof,  "No TOF",                "l")
    legend.AddEntry(hist_tof_hit, "TOF signal required",   "l")
    legend.AddEntry(hist,         "|#it{n}#sigma_{TOF}| < 2","l")
    legend.Draw()

    canvas.SaveAs("efficiency_comparison_li4.pdf")

if __name__ == "__main__":

    draw_protons()
    draw_li4()
