import ROOT
from ROOT import TCanvas, TGraphErrors, TH1F, TLegend, TLatex, TArrow, gStyle, \
                    kBlue, kRed, kGreen

from torchic.utils.root import set_root_object

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)


# --- Data ---
# from https://arxiv.org/abs/2311.11758
measurements = {
    "he4":  dict(x=0, y=0.83e-6, ex=0.0, ey_sta=0.22e-6, ey_sys=0.12e-6),
    "ahe4": dict(x=1, y=1.30e-6, ex=0.0, ey_sta=0.28e-6, ey_sys=0.18e-6),
}
upper_limits = {
    "li4":  dict(x=2, ex=0.25, ul=0.67e-6),
    "ali4": dict(x=3, ex=0.25, ul=1.01e-6),
}
shm_predictions = {
    "he4":  dict(x=0, y=0.85e-6, ex=0.25, ey_sta=0.),
    "ahe4": dict(x=1, y=0.82e-6, ex=0.25, ey_sta=0.),
    "li4":  dict(x=2, y=3.50e-6, ex=0.25, ey_sta=0.),
    "ali4": dict(x=3, y=3.38e-6, ex=0.25, ey_sta=0.),
}

# --- Helper functions ---
def make_stat_graph(name, d, color, marker_style, marker_size=1.5, line_style=1, line_width=1):
    g = TGraphErrors(1)
    set_root_object(g, name=name, marker_style=marker_style, marker_size=marker_size, marker_color=color, 
                    line_color=color, line_style=line_style, line_width=line_width)
    g.SetPoint(0, d["x"], d["y"])
    g.SetPointError(0, d["ex"], d["ey_sta"])
    return g

def make_sys_graph(name, d, color, box_width=0.25, alpha=0.3):
    g = TGraphErrors(1)
    set_root_object(g, name=name, marker_color=color, line_color=color, fill_color_alpha=(color, alpha), fill_style=1001)
    g.SetPoint(0, d["x"], d["y"])
    g.SetPointError(0, box_width, d["ey_sys"])
    return g

def make_upper_limit(name, d, color, marker_size=1.8, arrow_fraction=0.4):
    g = TGraphErrors(1)
    set_root_object(g, name=name, marker_size=marker_size, marker_color=color, line_color=color)
    g.SetPoint(0, d["x"], d["ul"])
    g.SetPointError(0, d["ex"], 0.)
    arrow = TArrow(d["x"], d["ul"], d["x"], d["ul"] * arrow_fraction, 0.03, "|>")
    set_root_object(arrow, line_color=color, fill_color=color, line_width=2)
    return g, arrow

def draw_alice_label(latex, x, y, lines):
    latex.SetTextFont(62)
    latex.DrawLatex(x, y, lines[0])
    latex.SetTextFont(42)
    for i, line in enumerate(lines[1:], start=1):
        latex.DrawLatex(x, y - i * 0.05, line)

def draw_four_points():

    canvas = TCanvas("c", "c", 800, 700)
    canvas.SetLogy()
    canvas.SetLeftMargin(0.15)
    canvas.SetBottomMargin(0.13)

    style = {
        "he4":  dict(color=kBlue+1, marker_style=20),
        "ahe4": dict(color=kBlue+1,  marker_style=20),
    }
    graphs_sta = {k: make_stat_graph(f"g_{k}_sta", measurements[k], **style[k]) for k in measurements}
    graphs_sys = {k: make_sys_graph( f"g_{k}_sys", measurements[k], style[k]["color"])  for k in measurements}
    graphs_shm = {k: make_stat_graph(f"g_{k}_shm", shm_predictions[k], color=kRed+1, marker_style=1, line_style=2, line_width=2) for k in shm_predictions}

    ul_color = 797
    graphs_ul, arrows = {}, {}
    for k, d in upper_limits.items():
        graphs_ul[k], arrows[k] = make_upper_limit(f"g_{k}", d, ul_color)

    frame = TH1F("frame", ";; #frac{1}{#it{N}_{events}} #frac{d#it{N}}{d#it{y}}", 4, -0.5, 3.5)
    frame.GetXaxis().SetLabelSize(0.05)
    frame.GetXaxis().SetBinLabel(1, "{}^{4}He")
    frame.GetXaxis().SetBinLabel(2, "{}^{4}#bar{He}")
    frame.GetXaxis().SetBinLabel(3, "{}^{4}Li")
    frame.GetXaxis().SetBinLabel(4, "{}^{4}#bar{Li}")

    frame.GetYaxis().SetRangeUser(1e-8, 1e-3)
    frame.GetYaxis().SetTitleOffset(1.6)
    frame.GetXaxis().SetNdivisions(5)
    frame.Draw()

    for g in graphs_sys.values(): g.Draw("E2 same")
    for g in graphs_sta.values(): g.Draw("P same")
    for g in graphs_ul.values():  g.Draw("P same")
    for g in graphs_shm.values(): g.Draw("L same")
    for a in arrows.values():     a.Draw()

    leg = TLegend(0.5, 0.65, 0.83, 0.88)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.025)
    leg.AddEntry(graphs_sta["he4"], "ALICE Pb#minusPb, #sqrt{#it{s}_{NN}} = 5.02 TeV", "pe")
    leg.AddEntry(graphs_ul["li4"], "#splitline{Pb#minusPb, #sqrt{#it{s}_{NN}} = 5.36 TeV}{95% confidence level}", "l")
    leg.AddEntry(graphs_shm["li4"], "Statistical Hadronisation Model", "l")
    leg.Draw()

    latex = TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.038)
    latex.SetTextFont(42)
    latex.DrawLatex(0.22, 0.83, "ALICE Preliminary")
    latex.DrawLatex(0.22, 0.78, "Pb#minusPb, 0#minus10%")

    canvas.SaveAs("final_plots/yields_matter_antimatter.pdf")

def draw_two_points():
    canvas = TCanvas("c", "c", 800, 700)
    #canvas.SetLogy()
    canvas.SetLeftMargin(0.15)
    canvas.SetBottomMargin(0.13)

    style = {
        "he4":  dict(color=kBlue+1, marker_style=20),
        "ahe4": dict(color=kBlue+1,  marker_style=20),
    }
    measurements_two_points = {k: v for k, v in measurements.items() if k[0] == 'a'}
    for measurement in measurements_two_points.values():
        measurement["x"] = (measurement["x"] - 1)/ 2
    graphs_sta = {k: make_stat_graph(f"g_{k}_sta", measurements[k], **style[k]) for k in measurements_two_points}
    graphs_sys = {k: make_sys_graph( f"g_{k}_sys", measurements[k], style[k]["color"])  for k in measurements_two_points}

    shm_predictions_two_points = {k: v for k, v in shm_predictions.items() if k[0] == 'a'}
    for shm_prediction in shm_predictions_two_points.values():
        shm_prediction["x"] = (shm_prediction["x"] - 1)/ 2
    graphs_shm = {k: make_stat_graph(f"g_{k}_shm", shm_predictions_two_points[k], color=kRed+1, marker_style=1, line_style=2, line_width=2) for k in shm_predictions_two_points}

    ul_color = 1
    graphs_ul, arrows = {}, {}
    for k, d in upper_limits.items():
        if k[0] != 'a':
            continue
        d["x"] = (d["x"] - 1)/ 2
        graphs_ul[k], arrows[k] = make_upper_limit(f"g_{k}", d, ul_color)

    frame = TH1F("frame", ";; #frac{1}{#it{N}_{events}} #frac{d#it{N}}{d#it{y}}", 2, -0.5, 1.5)
    frame.GetXaxis().SetLabelSize(0.06)
    frame.GetXaxis().SetBinLabel(1, "{}^{4}#bar{He}")
    frame.GetXaxis().SetBinLabel(2, "{}^{4}#bar{Li}")

    frame.GetYaxis().SetRangeUser(5e-8, 0.8e-5)
    frame.GetYaxis().SetTitleOffset(1.6)
    frame.GetXaxis().SetNdivisions(5)
    frame.Draw()

    for g in graphs_sys.values(): g.Draw("E2 same")
    for g in graphs_sta.values(): g.Draw("P same")
    for g in graphs_ul.values():  g.Draw("P same")
    for g in graphs_shm.values(): g.Draw("L same")
    for a in arrows.values():     a.Draw()

    leg = TLegend(0.18, 0.5, 0.9, 0.71)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.03)
    leg.SetMargin(0.1)
    leg.SetNColumns(2)
    #leg.AddEntry(graphs_sta["ahe4"], "ALICE Pb#minusPb, #sqrt{#it{s}_{NN}} = 5.02 TeV", "pe")
    #leg.AddEntry(graphs_ul["ali4"], "#splitline{Pb#minusPb, #sqrt{#it{s}_{NN}} = 5.36 TeV}{95% confidence level}", "l")
    leg.AddEntry(graphs_sta["ahe4"], "#splitline{#sqrt{#it{s}_{NN}} = 5.02 TeV}{#it{PLB} 858 (2024) 138943}", "pe")
    leg.AddEntry(graphs_shm["ali4"], "#splitline{#splitline{Thermal-FIST (GCE SHM)}{#it{T} = 156.4 MeV, #it{V} = 4233 fm^{3}}}{Nuclear excitation particle list}", "l")
    leg.AddEntry(graphs_ul["ali4"], "#splitline{#sqrt{#it{s}_{NN}} = 5.36 TeV}{95% confidence level}", "l")
    leg.Draw()

    latex = TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.038)
    latex.SetTextFont(42)
    latex.DrawLatex(0.22, 0.83, "ALICE Preliminary")
    latex.DrawLatex(0.22, 0.78, "Pb#minusPb, FT0C Centrality: 0#minus10%")

    canvas.SaveAs("final_plots/yields_antimatter.pdf")

if __name__ == '__main__':

    gStyle.SetPadTickX(1)
    gStyle.SetPadTickY(1)

    draw_four_points()
    draw_two_points()