
import numpy as np
from numpy import rint
from uncertainties import ufloat
from ROOT import TFile, TCanvas, TLegend, TGraphErrors, TF1, TH1F, TArrow, TPaveText, TGaxis, TLine
from torchic.utils.root import set_root_object, init_legend, set_alice_global_style
from torchic.utils.colors import get_color

EFFICIENCY = {
    'Matter': 0.089,
    'Antimatter': 0.064,
}
NCH_RUN3 = {
    '0-10%':  ufloat(1858, 34),
    '10-30%': ufloat(1051, 21),
    '30-50%': ufloat(455, 12),
    '50-80%': ufloat(123, 5),
}
    
weight_matter = 316374 # entries in the same event histogram for matter in 0-50%
weight_antimatter = 216224 # entries in the same event histogram for antimatter in 0-50%
EFFICIENCY['Both'] = (EFFICIENCY['Matter'] * weight_matter + EFFICIENCY['Antimatter'] * weight_antimatter) / (weight_matter + weight_antimatter)
N_EVENTS = {
    # centrality: (2023 + 2024 + 2025)
    '0-10%': 404485552 + 697710144 + 1099723512,
    '10-30%': 805844504 + 1388453096 + 2196980688,
    '30-50%': 804118540 + 1384085800 + 2193293968,
    '50-80%': 1202232388 + 2064590536 + 3279040944,
    '0-50%': 2014448596 + 3470249040 + 5489998168,
    '10-50%': 1609963044 + 2772538896 + 4390274656,
}
YIELD = {
    'Antimatter': {
        '0-10%': ufloat(113.2, 60.31),
        '10-30%': ufloat(75.93, 59.09),
        '30-50%': ufloat(134.9, 34.19),
        '50-80%': ufloat(2.49, 13.26),
        '0-50%': ufloat(332.4, 90.87),
        '10-50%': ufloat(224.9, 67.45),
    },
    'Matter': {
        '0-10%': ufloat(80.82, 70.28),
        '10-30%': ufloat(336.8, 73.1),
        '30-50%': ufloat(77.34, 38.96),
        '50-80%': ufloat(49.19, 18.4),
        '0-50%': ufloat(473.3, 106.),
        '10-50%': ufloat(406.6, 81.66),
    },
    'Both': {
        '0-10%': ufloat(199.30, 92.07) / 2.,
        '10-30%': ufloat(420.43, 91.97) / 2.,
        '30-50%': ufloat(231.99, 50.76) / 2.,
        '50-80%': ufloat(58.76, 22.19) / 2.,
        '0-50%': ufloat(819.6, 149.) / 2.,
        '10-50%': ufloat(650.18, 105.8) / 2.,
    }
}
YIELD_SYST = {
    'Both': {
        '0-10%': ufloat(199.30, np.sqrt(35.2**2 + ((288-135)/2)**2)) / 2.,
        '10-30%': ufloat(420.43, np.sqrt(28.6**2 + ((540-336)/2)**2)) / 2.,
        '30-50%': ufloat(231.99, np.sqrt(17.2**2 + ((300-187)/2)**2)) / 2.,
        '50-80%': ufloat(58.76, np.sqrt(6.0**2 + ((75-42)/2)**2)) / 2.,
        #'0-50%': ufloat(819.6, np.sqrt(35.2**2 + ((288-135)/2)**2)) / 2.,
        '10-50%': ufloat(650.18, np.sqrt(34.7**2 + ((831-526)/2)**2)) / 2.,
    }
}
UPPER_LIMIT = {
    'Both': {
        '0-10%': ufloat(447, 0.), # 95% CL upper limit with systematics
    }
}

def make_upper_limit(name, d, color, marker_size=1.8, arrow_fraction=0.4):
        g = TGraphErrors(1)
        set_root_object(g, name=name, marker_size=marker_size, marker_color=color, line_color=color)
        g.SetPoint(0, d["x"], d["ul"])
        g.SetPointError(0, d["ex"], 0.)
        arrow = TArrow(d["x"], d["ul"], d["x"], d["ul"] * arrow_fraction, 0.03, "|>")
        set_root_object(arrow, line_color=color, fill_color=color, line_width=2)
        return g, arrow

def correct_yield(sign: str, centrality: str, raw_yield: ufloat = None) -> ufloat:
    """
    Corrects the raw yield based on the centrality and sign.

    Parameters:
    - raw_yield: The raw yield to be corrected.
    - centrality: The centrality class (e.g., '0-10%', '10-20%', etc.).
    - sign: The charge sign ('positive' or 'negative').

    Returns:
    - The corrected yield.
    """
    
    if raw_yield is None:
        raw_yield = YIELD[sign][centrality]
    efficiency = EFFICIENCY[sign]
    n_events = N_EVENTS[centrality]
    corrected_yield = raw_yield / (efficiency * n_events * 2.) # Factor 2 for yield per unit rapidity
    return corrected_yield
            
def draw_yields(outfile: TFile) -> None:
    
    titles = {
        'Matter': '^{4}Li',
        'Antimatter': '^{4}#bar{Li}',
        'Both': '(^{4}Li + ^{4}#bar{Li})/2',
    }
    
    n_centralities = len(N_EVENTS.keys())
    for isign, sign in enumerate(['Matter', 'Antimatter', 'Both']):
        graph = TGraphErrors(n_centralities)
        set_root_object(graph, name=f"g_{sign}", title=titles[sign]+'; #LT d#it{N}_{ch}/ d#it{#eta} #GT^{|#it{#eta}| < 0.5}; #frac{1}{N_{events}} #frac{d#it{N}}{d#it{y}}', 
                        marker_style=20, marker_color=get_color(isign), line_color=get_color(isign), marker_size=1.4)
        graph.SetName(f"g_{sign}")
        
        for icent, (cent, mult) in enumerate(NCH_RUN3.items()):
            nucleus_yield = correct_yield(sign, cent)
            graph.SetPoint(icent,      mult.n, nucleus_yield.n)
            graph.SetPointError(icent, mult.s, nucleus_yield.s)
        
        canvas = TCanvas(f"canvas_{sign}", "yield", 800, 600)
        canvas.SetLogx()
        #canvas.SetLogy()
        hframe = canvas.DrawFrame(100, -2e-7 if sign == 'Antimatter' else 1e-8, 1.2*max(mult.n for mult in NCH_RUN3.values()), 1e-6, 
                                  titles[sign]+'; #LT d#it{N}_{ch}/ d#it{#eta} #GT^{|#it{#eta}| < 0.5}; #frac{1}{N_{events}} #frac{d#it{N}}{d#it{y}}; #frac{1}{N_{events}} #frac{d#it{N}}{d#it{y}}')
        graph.Draw("P SAME")
        outfile.cd()
        canvas.Write()
        
    g_ratio = TGraphErrors(n_centralities)
    set_root_object(g_ratio, name="g_ratio", title="^{4}#bar{Li} / ^{4}Li; #LT d#it{N}_{ch}/ d#it{#eta} #GT^{|#it{#eta}| < 0.5}; Ratio",
                    marker_style=20, marker_color=get_color(2), line_color=get_color(2), marker_size=1.4)
                                                                                     
    g_ratio.SetName("g_ratio")
    g_ratio.SetTitle("^{4}#bar{Li} / ^{4}Li")
    for i, (cent, mult) in enumerate(NCH_RUN3.items()):
        anti_yield = correct_yield('Antimatter', cent)
        matter_yield = correct_yield('Matter', cent)
        if matter_yield.n != 0:
            ratio     = anti_yield.n / matter_yield.n
            ratio_err = ratio * ((anti_yield.s / anti_yield.n)**2 + (matter_yield.s / matter_yield.n)**2)**0.5 if anti_yield.n != 0 else anti_yield.s / matter_yield.n
        else:
            ratio, ratio_err = 0., 0.
        g_ratio.SetPoint(i,      mult.n, ratio)
        g_ratio.SetPointError(i, mult.s, ratio_err)
    canvas = TCanvas("canvas_ratio", "ratio", 800, 600)
    pol0 = TF1("pol0", "pol0", 0, 1.2*max(mult.n for mult in NCH_RUN3.values()))
    g_ratio.Fit("pol0", "S")
    hframe = canvas.DrawFrame(100, -0.1, 1.2*max(mult.n for mult in NCH_RUN3.values()), 4, 
                              "^{4}#bar{Li}/^{4}Li; #LT d#it{N}_{ch}/ d#it{#eta} #GT^{|#it{#eta}| < 0.5}; ^{4}#bar{Li} / ^{4}Li")
    legend_ratio = init_legend(0.1, 0.6, 0.5, 0.8, fill_style=0, border_size=0)
    legend_ratio.AddEntry(g_ratio, "^{4}#bar{Li} / ^{4}Li", "P")
    legend_ratio.AddEntry(pol0, f"pol0: {pol0.GetParameter(0):.2f} #pm {pol0.GetParError(0):.2f}", "L")
    canvas.SetLogx()
    g_ratio.Draw("P SAME")
    legend_ratio.Draw()
    outfile.cd()
    canvas.Write()
    
    ## The plot we decided for
    
    #nucleus_yield_010 = correct_yield('Both', '0-10%')
    #nucleus_upper_limit_010 = ufloat(nucleus_yield_010.n + 2*nucleus_yield_010.s, 0.) # 95% CL upper limit
    
    nucleus_upper_limit_010 = correct_yield('Both', '0-10%', raw_yield=UPPER_LIMIT['Both']['0-10%'])
    graph_yields_010, arrow_010 = make_upper_limit(f"g_Both_010", {"x": 0.5, "ul": nucleus_upper_limit_010.n, "ex": 0.3}, get_color(2),
                                                   arrow_fraction=0.4)
    set_root_object(graph_yields_010, name=f"g_Both_010", title=titles['Both']+'; #LT d#it{N}_{ch}/ d#it{#eta} #GT^{|#it{#eta}| < 0.5}; #frac{1}{N_{events}} #frac{d#it{N}}{d#it{y}}', 
                    marker_style=0, marker_color=get_color(1), line_color=get_color(1), marker_size=1.4)
    set_root_object(arrow_010, line_color=get_color(1), line_width=2, fill_color=get_color(1))
    
    graph_yields_1050 = TGraphErrors(1)
    set_root_object(graph_yields_1050, name=f"g_Both_1050", title=titles['Both']+'; #LT d#it{N}_{ch}/ d#it{#eta} #GT^{|#it{#eta}| < 0.5}; #frac{1}{N_{events}} #frac{d#it{N}}{d#it{y}}', 
                    marker_style=20, marker_color=get_color(2), line_color=get_color(2), marker_size=1.4, line_width=2)
    graph_yields_1050.SetName(f"g_Both_1050")
    nucleus_yield_1050 = correct_yield('Both', '10-50%')
    graph_yields_1050.SetPoint(1, 1.5, nucleus_yield_1050.n)
    graph_yields_1050.SetPointError(1, 0, nucleus_yield_1050.s)
    
    graph_yields_1050_syst = TGraphErrors(1)
    set_root_object(graph_yields_1050_syst, name=f"g_Both_1050_syst", title=titles['Both']+'; #LT d#it{N}_{ch}/ d#it{#eta} #GT^{|#it{#eta}| < 0.5}; #frac{1}{N_{events}} #frac{d#it{N}}{d#it{y}}', 
                    marker_style=20, marker_color=get_color(2), line_color=get_color(2), marker_size=1.4,
                    fill_color_alpha=(get_color(2), 0.3))
    graph_yields_1050_syst.SetName(f"g_Both_1050_syst")
    nucleus_yield_1050_syst = correct_yield('Both', '10-50%', raw_yield=YIELD_SYST['Both']['10-50%'])
    graph_yields_1050_syst.SetPoint(1, 1.5, nucleus_yield_1050_syst.n)
    graph_yields_1050_syst.SetPointError(1, 0.3, nucleus_yield_1050_syst.s)
    
    graph_yields_alpha_010 = TGraphErrors(1)
    set_root_object(graph_yields_alpha_010, name=f"g_alpha_010", 
                    marker_style=20, marker_color=get_color(3), line_color=get_color(3),
                    marker_size=1.5, line_width=2)
    graph_yields_alpha_010.SetName(f"g_alpha_010")
    graph_yields_alpha_010.SetPoint(1, .5, 1.0e-6)
    graph_yields_alpha_010.SetPointError(1, 0, 0.19e-6)
    
    graph_yields_alpha_010_syst = TGraphErrors(1)
    set_root_object(graph_yields_alpha_010_syst, name=f"g_alpha_010",
                    marker_style=20, marker_color=get_color(3), line_color=get_color(3), marker_size=1.4,
                    fill_color_alpha=(get_color(3), 0.3))
    graph_yields_alpha_010_syst.SetName(f"g_alpha_010")
    graph_yields_alpha_010_syst.SetPoint(1, .5, 1.0e-6)
    graph_yields_alpha_010_syst.SetPointError(1, 0.3, 0.1e-6)
    
    graph_yields_li4_fist = TGraphErrors(1)
    set_root_object(graph_yields_li4_fist, name=f"g_li4_fist",
                    marker_style=20, marker_color=get_color(0), line_color=get_color(0), marker_size=0,
                    line_style=2, line_width=2)
    graph_yields_li4_fist.SetName(f"g_li4_fist")
    graph_yields_li4_fist.SetPoint(1, .5, 8.2004e-06)
    graph_yields_li4_fist.SetPointError(1, 0.3, 0)
    
    graph_yields_he4_fist = TGraphErrors(1)
    set_root_object(graph_yields_he4_fist, name=f"g_he4_fist",
                    marker_style=20, marker_color=get_color(4), line_color=get_color(4), marker_size=0,
                    line_style=2, line_width=2)
    graph_yields_he4_fist.SetName(f"g_he4_fist")
    graph_yields_he4_fist.SetPoint(1, .5, 7.9e-07)
    graph_yields_he4_fist.SetPointError(1, 0.3, 0)
    
    graph_yields_he4_fist_1050 = TGraphErrors(1)
    set_root_object(graph_yields_he4_fist_1050, name=f"g_he4_fist_1050",
                    marker_style=20, marker_color=get_color(7), line_color=get_color(7), marker_size=0,
                    line_style=2, line_width=2)
    graph_yields_he4_fist_1050.SetName(f"g_he4_fist_1050")
    graph_yields_he4_fist_1050.SetPoint(1, 1.5, 4.9659e-07)
    graph_yields_he4_fist_1050.SetPointError(1, 0.3, 0)
    
    canvas_yields_centralities = TCanvas(f"canvas_Both_centralities", "yield", 800, 600)
    canvas_yields_centralities.SetLeftMargin(0.15)
    canvas_yields_centralities.SetBottomMargin(0.15)
    canvas_yields_centralities.SetLogy()
    hframe_yields_centralities = TH1F('hfame', '; FT0C Centrality (%); #frac{1}{N_{events}} #frac{d#it{N}}{d#it{y}}; #frac{1}{N_{events}} #frac{d#it{N}}{d#it{y}}', 
                                      4, 0, 4)
    hframe_yields_centralities.GetXaxis().SetBinLabel(1, '0-10%')
    hframe_yields_centralities.GetXaxis().SetBinLabel(2, '10-50%')
    
    hframe_yields_centralities.GetXaxis().SetTitleSize(0.045)
    hframe_yields_centralities.GetXaxis().SetLabelSize(0.08)
    hframe_yields_centralities.GetXaxis().SetTitleOffset(1.3)
    hframe_yields_centralities.GetYaxis().SetTitleSize(0.045)
    hframe_yields_centralities.GetYaxis().SetLabelSize(0.045)
    hframe_yields_centralities.GetYaxis().SetTitleOffset(1.5)
    
    # for log scale
    hframe_yields_centralities.SetMaximum(4e-5)
    hframe_yields_centralities.SetMinimum(1.3e-7)
    
    #hframe_yields_centralities.SetMaximum(1e-5)
    #hframe_yields_centralities.SetMinimum(0.2e-8)
    
    text = TPaveText(0.2, 0.78, 0.5, 0.86, "NDC")
    text.SetBorderSize(0)
    text.SetFillStyle(0)
    text.SetTextSize(0.04)
    text.SetTextFont(42)
    text.AddText("ALICE Pb#minusPb")
    #text.AddText("")
    
    leg = TLegend(0.55, 0.58, 0.86, 0.86)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.03)
    leg.SetMargin(0.1)
    leg.SetNColumns(1)
    leg.SetColumnSeparation(0.05)
    leg.SetMargin(0.1)
    #leg.AddEntry(graphs_sta["ahe4"], "ALICE Pb#minusPb, #sqrt{#it{s}_{NN}} = 5.02 TeV", "pe")
    #leg.AddEntry(graphs_ul["ali4"], "#splitline{Pb#minusPb, #sqrt{#it{s}_{NN}} = 5.36 TeV}{95% confidence level}", "l")
    leg.AddEntry(graph_yields_010, "#splitline{(^{4}Li + ^{4}#bar{Li})/2    #sqrt{#it{s}_{NN}} = 5.36 TeV}{95% confidence level}", "l")
    leg.AddEntry(graph_yields_1050, "(^{4}Li + ^{4}#bar{Li})/2    #sqrt{#it{s}_{NN}} = 5.36 TeV", "pe")
    leg.AddEntry(graph_yields_alpha_010, "#splitline{(^{4}He + ^{4}#bar{He})/2    #sqrt{#it{s}_{NN}} = 5.02 TeV}{#it{PLB} 858 (2024) 138943}", "pe")
    leg.Draw()
    
    leg_fist = TLegend(0.55, 0.2, 0.86, 0.56)
    leg_fist.SetHeader("#splitline{   Thermal-FIST (GCE SHM)}{   Nuclear excitation particle list}")
    leg_fist.SetBorderSize(0)
    leg_fist.SetFillStyle(0)
    leg_fist.SetTextSize(0.03)
    leg_fist.SetMargin(0.1)
    leg_fist.SetNColumns(1)
    leg_fist.SetColumnSeparation(0.05)
    leg_fist.SetMargin(0.1)
    leg_fist.AddEntry(graph_yields_he4_fist, "#splitline{(^{4}He + ^{4}#bar{He})/2}{#it{T} = 156.4 MeV, #it{V} = 4233 fm^{3}}", "l")
    leg_fist.AddEntry(graph_yields_he4_fist_1050, "#splitline{(^{4}He + ^{4}#bar{He})/2}{#it{T} = 158.8 MeV, #it{V} = 1804 fm^{3}}", "l")
    leg_fist.AddEntry(graph_yields_li4_fist, "#splitline{(^{4}Li + ^{4}#bar{Li})/2}{#it{T} = 156.6 MeV, #it{V} = 4459 fm^{3}}", "l")
    leg_fist.Draw()

    nucleus_yield_010_with_syst = correct_yield('Both', '0-10%', raw_yield=ufloat(YIELD['Both']['0-10%'].n, 39))
    print(f"Corrected yield for Both in 0-10%: {nucleus_yield_010_with_syst:.2e} #pm {nucleus_yield_010_with_syst.s:.2e} (stat + syst)")
    # print(f"Corrected yield for Both in 0-10%: {nucleus_yield_010:.2e}")
    print(f"Corrected yield for Both in 10-50%: {nucleus_yield_1050:.2e} #pm {nucleus_yield_1050.s:.2e} (stat) #pm {nucleus_yield_1050_syst.s:.2e} (syst)")
    print(f"95% CL upper limit for Both in 0-10%: {nucleus_upper_limit_010:.2e}")
    # print(f"Corrected yield for (^{{4}}He + ^{{4}}#bar{{He}})/2 in 0-10%: {correct_yield('Both', '0-10%'):.2e}")
    # print(f"Thermal-FIST prediction for ^{{4}}Li: 8.2004e-6")
    print(f"Upper limit to the ^{{4}}Li/^{{4}}He ratio in 0-10%: {nucleus_upper_limit_010.n / correct_yield('Both', '0-10%').n:.2e}")
    
    
    hframe_yields_centralities.Draw()
    graph_yields_010.Draw("P SAME")
    arrow_010.Draw()
    graph_yields_1050.Draw("P SAME")
    graph_yields_1050_syst.Draw("E2 SAME")
    graph_yields_alpha_010.Draw("P SAME")
    graph_yields_alpha_010_syst.Draw("E2 SAME")
    graph_yields_li4_fist.Draw("P SAME")
    graph_yields_he4_fist.Draw("P SAME")
    graph_yields_he4_fist_1050.Draw("P SAME")
    leg.Draw()
    leg_fist.Draw()
    text.Draw('same')
    outfile.cd()
    canvas_yields_centralities.Write()
    canvas_yields_centralities.SaveAs("figures/corrected_yields_centralities.pdf")

def print_corrected_yields():
    
    print('----------------------------------------------------------------')
    for sign in ['Matter', 'Antimatter', 'Both']:
        for centrality in ['0-10%', '10-30%', '30-50%', '50-80%', '0-50%', '10-50%']:
            corrected_yield = correct_yield(sign, centrality, raw_yield=YIELD[sign][centrality])
            corrected_yield_syst = correct_yield(sign, centrality, raw_yield=YIELD_SYST[sign][centrality]) if sign in YIELD_SYST and centrality in YIELD_SYST[sign] else None
            if corrected_yield_syst is not None:
                print(f"Corrected yield for {sign} in {centrality}: {corrected_yield.n:.2e} ± {corrected_yield.s:.2e} (stat) ± {corrected_yield_syst.s:.2e} (syst)")
            else:
                print(f"Corrected yield for {sign} in {centrality}: {corrected_yield.n:.2e} ± {corrected_yield.s:.2e} (stat)")   
    print('----------------------------------------------------------------')
    
if __name__ == "__main__":
    
    set_alice_global_style()

    print_corrected_yields()
            
    outfile = TFile.Open("output/corrected_yields.root", "RECREATE")
    draw_yields(outfile)
    outfile.Close()