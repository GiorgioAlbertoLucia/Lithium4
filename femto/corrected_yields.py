
import numpy as np
from numpy import rint
from uncertainties import ufloat
from ROOT import TFile, TCanvas, TLegend, TGraphErrors, TF1, TH1F, TArrow, TPaveText, TGaxis, TLine
from torchic.utils.root import set_root_object, init_legend, set_alice_global_style
from torchic.utils.colors import get_color

# the efficiency should be corrected with the signal loss
EFFICIENCY_DICT = {
    '0-10%':
    {
        'Matter':           {'2023': 0.049, '2024': 0.088, '2025': 0.095},
        'Antimatter':       {'2023': 0.036, '2024': 0.060, '2025': 0.065},
        # entries in the same event histogram for matter and antimatter in 0-50% centrality
        'WeightMatter':     {'2023': 897538, '2024': 1401744, '2025': 2584925},
        'WeightAntimatter': {'2023': 618707, '2024': 975959, '2025': 1788444},
    },
    '10-50%':
    {
        'Matter':           {'2023': 0.070, '2024': 0.093, '2025': 0.110},
        'Antimatter':       {'2023': 0.051, '2024': 0.066, '2025': 0.076},
        # entries in the same event histogram for matter and antimatter in 0-50% centrality
        'WeightMatter':     {'2023': 969950, '2024': 1488593, '2025': 2747092},
        'WeightAntimatter': {'2023': 689387, '2024': 1058768, '2025': 1910888},
    },
}
EFFICIENCY_DICT['10-30%'] = EFFICIENCY_DICT['10-50%'] # proxy
EFFICIENCY_DICT['30-50%'] = EFFICIENCY_DICT['10-50%'] # proxy
EFFICIENCY_DICT['50-80%'] = EFFICIENCY_DICT['10-50%'] # proxy
EFFICIENCY_DICT['0-50%'] = EFFICIENCY_DICT['10-50%'] # proxy

SIGNAL_LOSS = {
    '0-10%':    {'2023': 0.8081, '2024': 0.7395, '2025': 0.7381},
    '10-30%':   {'2023': 0.8067, '2024': 0.7411, '2025': 0.7422},
    '30-50%':   {'2023': 0.8067, '2024': 0.7405, '2025': 0.7425},
    '50-80%':   {'2023': 0.8065, '2024': 0.7411, '2025': 0.7423},
    '0-50%':    {'2023': 0.8069, '2024': 0.7405, '2025': 0.7425}, # proxy
    '10-50%':   {'2023': 0.8069, '2024': 0.7405, '2025': 0.7425}
}


EFFICIENCY = {
    centrality: {
        'Matter': np.average(
    np.asarray(list(EFFICIENCY_DICT[centrality]['Matter'].values())) *
            np.asarray(list(SIGNAL_LOSS[centrality].values())),
            weights=list(EFFICIENCY_DICT[centrality]['WeightMatter'].values())
        ),

        'Antimatter': np.average(
            np.asarray(list(EFFICIENCY_DICT[centrality]['Antimatter'].values())) *
            np.asarray(list(SIGNAL_LOSS[centrality].values())),
            weights=list(EFFICIENCY_DICT[centrality]['WeightAntimatter'].values())
        ),

        'Both': np.average(
            np.concatenate([
                np.asarray(list(EFFICIENCY_DICT[centrality]['Matter'].values())),
                np.asarray(list(EFFICIENCY_DICT[centrality]['Antimatter'].values()))
            ]) *
            np.concatenate([
                np.asarray(list(SIGNAL_LOSS[centrality].values())),
                np.asarray(list(SIGNAL_LOSS[centrality].values()))
            ]),
            weights=(
                list(EFFICIENCY_DICT[centrality]['WeightMatter'].values()) +
                list(EFFICIENCY_DICT[centrality]['WeightAntimatter'].values())
            )
        )
    } for centrality in EFFICIENCY_DICT.keys()
}

for centrality in EFFICIENCY:
    EFFICIENCY[centrality]['Both_GS'] = EFFICIENCY[centrality]['Both']

print("EFFICIENCY_DICT:", EFFICIENCY_DICT, '\n')
print("EFFICIENCY:", EFFICIENCY)

NCH_RUN3 = {
    '0-10%':  ufloat(1858, 34),
    '10-30%': ufloat(1051, 21),
    '30-50%': ufloat(455, 12),
    '50-80%': ufloat(123, 5),
}


N_EVENTS = {
    # centrality: (2023 + 2024 + 2025)
    '0-10%':    {'2023': 404_485_552, '2024': 697710144, '2025': 1099723512},
    '10-30%':   {'2023': 805_844_504, '2024': 1_388_453_096, '2025': 2_196_980_688},
    '30-50%':   {'2023': 804_118_540, '2024': 1_384_085_800, '2025': 2_193_293_968},
    '50-80%':   {'2023': 1_202_232_388, '2024': 2_064_590_536, '2025': 3_279_040_944},
    '0-50%':    {'2023': 2_014_448_596, '2024': 3_470_249_040, '2025': 5_489_998_168},
    '10-50%':   {'2023': 1_609_963_044, '2024': 2_772_538_896, '2025': 4_390_274_656},
}
# event_efficiency = event_loss * event_splitting
EVENT_EFFICIENCY = {
    '0-10%':    {'2023': 0.7854, '2024': 0.7837, '2025': 0.7847},
    '10-30%':   {'2023': 0.7740, '2024': 0.7649, '2025': 0.7692},
    '30-50%':   {'2023': 0.7707, '2024': 0.7593, '2025': 0.7714},
    '50-80%':   {'2023': 0.7759, '2024': 0.7549, '2025': 0.7647},
    '0-50%':    {'2023': 0.7775, '2024': 0.7463, '2025': 0.7589},
    '10-50%':   {'2023': 0.7775, '2024': 0.7463, '2025': 0.7589}
}
N_EVENTS_CORRECTED = {centrality: {year: N_EVENTS[centrality][year] / EVENT_EFFICIENCY[centrality][year] for year in N_EVENTS[centrality]} for centrality in EVENT_EFFICIENCY}
for centrality in N_EVENTS_CORRECTED:
    N_EVENTS_CORRECTED[centrality]['total'] = sum(N_EVENTS_CORRECTED[centrality].values())

YIELD = {
    'Antimatter': {
        '0-10%': ufloat(113.2, 60.31),
        '10-30%': ufloat(75.93, 59.09),
        '30-50%': ufloat(134.9, 34.19),
        '50-80%': ufloat(2.49, 13.26),
        '0-50%': ufloat(332.4, 84.04),
        # PF
        #'10-50%': ufloat(224.9, 67.45),
        
        # new
        '10-50%': ufloat(332.82, 67.45),
    },
    'Matter': {
        '0-10%': ufloat(80.82, 70.28),
        '10-30%': ufloat(336.8, 73.1),
        '30-50%': ufloat(77.34, 38.96),
        '50-80%': ufloat(49.19, 18.4),
        '0-50%': ufloat(473.3, 106.),
        
        # PF
        #'10-50%': ufloat(406.6, 81.66),
        
        # new
        '10-50%': ufloat(547.91, 101.80),
    },
    'Both': {
        # PF
        ### '0-10%': ufloat(184, 115) / 2.,
        ### '10-30%': ufloat(398, 118) / 2.,
        ### '30-50%': ufloat(221, 63) / 2.,
        ### '50-80%': ufloat(58.76, 22.19) / 2.,
        ### '0-50%': ufloat(819.6, 149.) / 2.,
        ### '10-50%': ufloat(616, 133) / 2.,
        
        # new
        '0-10%': ufloat(312.25, 114.15) / 2.,
        '10-50%': ufloat(878.64, 133.60) / 2.,
    },
    'Both_GS': {
      '0-10%': ufloat(184.23, 90.) / 2.,
      '10-50%': ufloat(493.24, 104.) / 2.,
    }
}
YIELD_SYST = {
    'Antimatter': {
            # value, sqrt( c.v.^2 + lambdaR^2 + width_err^2)
            '10-50%': ufloat(332.82, np.sqrt(20.28**2 + 25.70**2 + 1**2)),
        },
    'Matter': {
            # value, sqrt( c.v.^2 + lambdaR^2 + width_err^2)
            '10-50%': ufloat(547.91, np.sqrt(24.00**2 + 35.87**2 + 0.5**2)),
        },
    
    'Both': {
        
        # PF
        #'0-10%': ufloat(199.30, np.sqrt(35.2**2 + ((288-135)/2)**2)) / 2.,
        # value, sqrt( c.v.^2 + lambdaR^2 + width_err^2 + strong_int^2)
        ### '0-10%': ufloat(184, np.sqrt(35**2 + 44**2 + 3**2)) / 2.,
        ### '10-30%': ufloat(398, np.sqrt(28**2 + 60**2 + 5**2)) / 2.,
        ### '30-50%': ufloat(221, np.sqrt(17**2 + 34**2 + 1**2)) / 2.,
        ### '50-80%': ufloat(58.76, np.sqrt(6.0**2 + ((75-42)/2)**2)) / 2., ## older values
        ### #'0-50%': ufloat(819.6, np.sqrt(35.2**2 + ((288-135)/2)**2)) / 2.,
        ### #'10-50%': ufloat(650.18, np.sqrt(34.7**2 + ((831-526)/2)**2)) / 2.,
        ### '10-50%': ufloat(616, np.sqrt(34**2 + 91**2 + 7**2)) / 2.,
        
        # new 
        # value, sqrt( c.v.^2 + lambdaR^2 + width_err^2 + strong_int^2)
        '0-10%': ufloat(312.25, np.sqrt(34.80**2 + 34.13**2 + 2**2 + 0**2)) / 2.,
        '10-50%': ufloat(878.64, np.sqrt(34.48**2 + 61.50**2 + 1**2 + 0**2)) / 2.,
    },
    'Both_GS': {
      '0-10%': ufloat(184.23, np.sqrt(27.2**2 + (30)**2)) / 2.,
      '10-50%': ufloat(493.24, np.sqrt(26.5**2 + (64)**2)) / 2.,
    }
}
UPPER_LIMIT = {
    'Both': {
        '0-10%': ufloat(446, 0.) / 2., # 95% CL upper limit with systematics
    }
}
TITLES = {
    'Matter': '^{4}Li',
    'Antimatter': '^{4}#bar{Li}',
    'Both': '(^{4}Li + ^{4}#bar{Li})/2',
}

THERMAL_FIST_PREDICTIONS = {
    'Li4': {
        '0-10%': ufloat(8.2004e-06, 0.),
        '10-50%': ufloat((4.9697e-06 + 4.8010e-06)/2, 0.),
    },
    'Li4_gs': {
        '0-10%': ufloat((3.4970e-06 + 3.3763e-06)/2, 0.),
        '10-50%': ufloat((2.0825e-06 + 2.0118e-06)/2, 0.),
    },
    'He4': {
        '0-10%': ufloat(7.9e-07, 0.), # Run 2
        '10-50%': ufloat(4.9659e-07, 0.), # Run 3
    }
}

HE4_PAPER = {
    '0-10%': ufloat(1.0e-6, 0.19e-6), # Run 2
}
HE4_PAPER_SYST = {
    '0-10%': ufloat(1.0e-6, 0.1e-6), # Run 2
}
H4L_PAPER = {
    '0-10%': ufloat(0.78e-6, 0.19e-6), # Run 2
}
H4L_PAPER_SYST = {
    '0-10%': ufloat(0.78e-6, 0.17e-6), # Run 2
}
HE4L_PAPER = {
    '0-10%': ufloat(1.08e-6, 0.34e-6), # Run 2
}
HE4L_PAPER_SYST = {
    '0-10%': ufloat(1.08e-6, 0.20e-6), # Run 2
}
 
# ============================================================================
# Helper functions
# ============================================================================
 
def make_upper_limit(name, d, color, marker_size=1.8, arrow_fraction=0.4):
        g = TGraphErrors(1)
        set_root_object(g, name=name, marker_size=marker_size, marker_color=color, line_color=color)
        g.SetPoint(0, d["x"], d["ul"])
        g.SetPointError(0, d["ex"], 0.)
        arrow = TArrow(d["x"], d["ul"], d["x"], d["ul"] * arrow_fraction, 0.03, "|>")
        set_root_object(arrow, line_color=color, fill_color=color, line_width=2)
        return g, arrow
 
def make_point_graph(name, x, y, ex, ey, **style) -> TGraphErrors:
    """
    Builds a single-point TGraphErrors, e.g. for a measured value or a
    model prediction shown at one x position. `style` is forwarded to
    set_root_object (title, marker_style, marker_color, line_color, ...).
    """
    g = TGraphErrors(1)
    set_root_object(g, name=name, **style)
    g.SetName(name)
    g.SetPoint(1, x, y)
    g.SetPointError(1, ex, ey)
    return g
 
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
    efficiency = EFFICIENCY[centrality][sign]
    n_events = N_EVENTS_CORRECTED[centrality]['total']
    corrected_yield = raw_yield / (efficiency * n_events * 2.) # Factor 2 for yield per unit rapidity
    return corrected_yield
 
# ============================================================================
# Plotting functions
# ============================================================================
 
def draw_yield_vs_multiplicity(outfile: TFile) -> None:
    """
    Draws, for each sign (Matter, Antimatter, Both), the corrected yield
    as a function of <dNch/deta> across the NCH_RUN3 centrality classes.
    """
    n_centralities = len(N_EVENTS.keys())
    for isign, sign in enumerate(['Matter', 'Antimatter', 'Both']):
        graph = TGraphErrors(n_centralities)
        set_root_object(graph, name=f"g_{sign}", title=TITLES[sign]+'; #LT d#it{N}_{ch}/ d#it{#eta} #GT^{|#it{#eta}| < 0.5}; #frac{1}{N_{events}} #frac{d#it{N}}{d#it{y}}', 
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
                                  TITLES[sign]+'; #LT d#it{N}_{ch}/ d#it{#eta} #GT^{|#it{#eta}| < 0.5}; #frac{1}{N_{events}} #frac{d#it{N}}{d#it{y}}; #frac{1}{N_{events}} #frac{d#it{N}}{d#it{y}}')
        graph.Draw("P SAME")
        outfile.cd()
        canvas.Write()
 
def draw_ratio_vs_multiplicity(outfile: TFile) -> None:
    """
    Draws the Antimatter/Matter corrected-yield ratio as a function of
    <dNch/deta>, together with a pol0 fit.
    """
    n_centralities = len(N_EVENTS.keys())
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
 
def draw_yields_vs_centrality(outfile: TFile) -> None:
    """
    Draws the final summary plot: (Li4 + Li4bar)/2 yield (upper limit in
    0-10%, measured point in 10-50%) compared to the He4 measurement and
    to Thermal-FIST predictions.
    """
    
    ## The plot we decided for
    
    #nucleus_yield_010 = correct_yield('Both', '0-10%')
    #nucleus_upper_limit_010 = ufloat(nucleus_yield_010.n + 2*nucleus_yield_010.s, 0.) # 95% CL upper limit
    
    both_title = TITLES['Both']+'; #LT d#it{N}_{ch}/ d#it{#eta} #GT^{|#it{#eta}| < 0.5}; #frac{1}{N_{events}} #frac{d#it{N}}{d#it{y}}'
 
    nucleus_upper_limit_010 = correct_yield('Both', '0-10%', raw_yield=UPPER_LIMIT['Both']['0-10%'])
    graph_yields_010, arrow_010 = make_upper_limit(f"g_Both_010", {"x": 0.5, "ul": nucleus_upper_limit_010.n, "ex": 0.3}, get_color(2),
                                                   arrow_fraction=0.4)
    # style overridden on purpose: make_upper_limit sets an initial style,
    # then this call restyles it (color 2 -> color 1) and adds the title
    set_root_object(graph_yields_010, name=f"g_Both_010", title=both_title,
                    marker_style=0, marker_color=get_color(1), line_color=get_color(1), marker_size=1.4)
    set_root_object(arrow_010, line_color=get_color(1), line_width=2, fill_color=get_color(1))
    
    nucleus_yield_1050 = correct_yield('Both', '10-50%')
    nucleus_yield_1050_syst = correct_yield('Both', '10-50%', raw_yield=YIELD_SYST['Both']['10-50%'])
 
    graph_yields_1050 = make_point_graph(
        "g_Both_1050", x=1.5, y=nucleus_yield_1050.n, ex=0, ey=nucleus_yield_1050.s,
        title=both_title, marker_style=20, marker_color=get_color(1), line_color=get_color(1),
        marker_size=1.4, line_width=2)
    graph_yields_1050_syst = make_point_graph(
        "g_Both_1050_syst", x=1.5, y=nucleus_yield_1050_syst.n, ex=0.3, ey=nucleus_yield_1050_syst.s,
        title=both_title, marker_style=20, marker_color=get_color(1), line_color=get_color(1),
        marker_size=1.4, fill_color_alpha=(get_color(1), 0.3))
 
    graph_yields_alpha_010 = make_point_graph(
        "g_alpha_010", x=.5, y=HE4_PAPER['0-10%'].n, ex=0, ey=HE4_PAPER['0-10%'].s,
        marker_style=33, marker_color=get_color(3), line_color=get_color(3), marker_size=2.3, line_width=2)
    graph_yields_alpha_010_syst = make_point_graph(
        "g_alpha_010", x=.5, y=HE4_PAPER_SYST['0-10%'].n, ex=0.3, ey=HE4_PAPER_SYST['0-10%'].s,
        marker_style=33, marker_color=get_color(3), line_color=get_color(3), marker_size=2.3,
        fill_color_alpha=(get_color(3), 0.3))
    
    graph_yields_h4l_010 = make_point_graph(
        "g_h4l_010", x=.5, y=H4L_PAPER['0-10%'].n, ex=0, ey=H4L_PAPER['0-10%'].s,
        marker_style=34, marker_color=get_color(4), line_color=get_color(4), marker_size=2.3, line_width=2)
    graph_yields_h4l_010_syst = make_point_graph(
        "g_h4l_010_syst", x=.5, y=H4L_PAPER_SYST['0-10%'].n, ex=0.3, ey=H4L_PAPER_SYST['0-10%'].s,
        marker_style=34, marker_color=get_color(4), line_color=get_color(4), marker_size=2.3,
        fill_color_alpha=(get_color(4), 0.3))
    
    graph_yields_he4l_010 = make_point_graph(
        "g_he4l_010", x=.5, y=HE4L_PAPER['0-10%'].n, ex=0, ey=HE4L_PAPER_SYST['0-10%'].s,
        marker_style=35, marker_color=get_color(8), line_color=get_color(8), marker_size=2.3, line_width=2)
    graph_yields_he4l_010_syst = make_point_graph(
        "g_he4l_010_syst", x=.5, y=HE4L_PAPER_SYST['0-10%'].n, ex=0.3, ey=HE4L_PAPER_SYST['0-10%'].s,
        marker_style=35, marker_color=get_color(8), line_color=get_color(8), marker_size=2.3,
        fill_color_alpha=(get_color(8), 0.3))
 
    graph_yields_li4_fist = make_point_graph(
        "g_li4_fist", x=.5, y=THERMAL_FIST_PREDICTIONS['Li4']['0-10%'].n, ex=0.3, ey=0,
        marker_style=20, marker_color=get_color(2), line_color=get_color(2), marker_size=0,
        line_style=2, line_width=2)
    graph_yields_li4_fist_gs = make_point_graph(
        "g_li4_fist_gs", x=.5, y=THERMAL_FIST_PREDICTIONS['Li4_gs']['0-10%'].n, ex=0.3, ey=0,
        marker_style=20, marker_color=get_color(0), line_color=get_color(0), marker_size=0,
        line_style=9, line_width=2)
    
    graph_yields_li4_fist_1050 = make_point_graph(
        "g_li4_fist_1050", x=1.5, y=THERMAL_FIST_PREDICTIONS['Li4']['10-50%'].n, ex=0.3, ey=0,
        marker_style=20, marker_color=get_color(2), line_color=get_color(2), marker_size=0,
        line_style=2, line_width=2)
    graph_yields_li4_fist_gs_1050 = make_point_graph(
        "g_li4_fist_gs_1050", x=1.5, y=THERMAL_FIST_PREDICTIONS['Li4_gs']['10-50%'].n, ex=0.3, ey=0,
        marker_style=20, marker_color=get_color(0), line_color=get_color(0), marker_size=0,
        line_style=9, line_width=2)
 
    graph_yields_he4_fist = make_point_graph(
        "g_he4_fist", x=.5, y=THERMAL_FIST_PREDICTIONS['He4']['0-10%'].n, ex=0.3, ey=0,
        marker_style=20, marker_color=get_color(3), line_color=get_color(3), marker_size=0,
        line_style=5, line_width=2)
 
    graph_yields_he4_fist_1050 = make_point_graph(
        "g_he4_fist_1050", x=1.5, y=THERMAL_FIST_PREDICTIONS['He4']['10-50%'].n, ex=0.3, ey=0,
        marker_style=20, marker_color=get_color(3), line_color=get_color(3), marker_size=0,
        line_style=5, line_width=2)
    
    canvas_yields_centralities = TCanvas(f"canvas_Both_centralities", "yield", 900, 650)
    canvas_yields_centralities.SetLeftMargin(0.15)
    canvas_yields_centralities.SetBottomMargin(0.15)
    canvas_yields_centralities.SetLogy()
    hframe_yields_centralities = TH1F('hfame', '; FT0C Centrality (%); #frac{1}{#it{N}_{events}} #frac{d#it{N}}{d#it{y}}; #frac{1}{N_{events}} #frac{d#it{N}}{d#it{y}}', 
                                      4, 0, 4)
                                      #3, np.array([0, 1, 2, 4], dtype=np.float32))
    hframe_yields_centralities.GetXaxis().SetBinLabel(1, '0-10%')
    hframe_yields_centralities.GetXaxis().SetBinLabel(2, '10-50%')
    
    hframe_yields_centralities.GetXaxis().SetTitleSize(0.045)
    hframe_yields_centralities.GetXaxis().SetLabelSize(0.08)
    hframe_yields_centralities.GetXaxis().SetTitleOffset(1.3)
    hframe_yields_centralities.GetYaxis().SetTitleSize(0.045)
    hframe_yields_centralities.GetYaxis().SetLabelSize(0.045)
    hframe_yields_centralities.GetYaxis().SetTitleOffset(1.5)
    
    # for log scale
    hframe_yields_centralities.SetMaximum(1.0e-5)
    hframe_yields_centralities.SetMinimum(1.e-7)
    
    #hframe_yields_centralities.SetMaximum(1e-5)
    #hframe_yields_centralities.SetMinimum(0.2e-8)
    
    text = TPaveText(0.55, 0.82, 0.86, 0.88, "NDC")
    text.SetBorderSize(0)
    text.SetFillStyle(0)
    text.SetTextSize(0.04)
    text.SetTextFont(42)
    text.AddText("ALICE Pb#minusPb")
    #text.AddText("")
    
    leg = init_legend(0.52, 0.38, 0.86, 0.82, fill_style=0, border_size=0, text_size=0.03, 
                      margin=0.1, column_separation=0.05)
    #leg.AddEntry(graphs_sta["ahe4"], "ALICE Pb#minusPb, #sqrt{#it{s}_{NN}} = 5.02 TeV", "pe")
    #leg.AddEntry(graphs_ul["ali4"], "#splitline{Pb#minusPb, #sqrt{#it{s}_{NN}} = 5.36 TeV}{95% confidence level}", "l")
    leg.AddEntry(graph_yields_010, "#splitline{(^{4}Li + ^{4}#bar{Li})/2    #sqrt{#it{s}_{NN}} = 5.36 TeV}{95% confidence level}", "l")
    leg.AddEntry(graph_yields_1050, "(^{4}Li + ^{4}#bar{Li})/2    #sqrt{#it{s}_{NN}} = 5.36 TeV", "pe")
    leg.AddEntry(graph_yields_alpha_010, "#splitline{(^{4}He + ^{4}#bar{He})/2    #sqrt{#it{s}_{NN}} = 5.02 TeV}{#it{PLB} 858 (2024) 138943}", "pe")
    leg.AddEntry(graph_yields_h4l_010, "#splitline{(^{4}_{#Lambda}H + ^{4}_{#bar{#Lambda}}#bar{H})/2    #sqrt{#it{s}_{NN}} = 5.02 TeV}{#it{PRL} 134 (2025) 162301}", "pe")
    leg.AddEntry(graph_yields_he4l_010, "#splitline{(^{4}_{#Lambda}He + ^{4}_{#bar{#Lambda}}#bar{He})/2    #sqrt{#it{s}_{NN}} = 5.02 TeV}{#it{PRL} 134 (2025) 162301}", "pe")
    leg.Draw()
    
    leg_fist = init_legend(0.57, 0.16, 0.86, 0.36, fill_style=0, border_size=0, text_size=0.03, 
                           margin=0.1, column_separation=0.05)
    leg_fist.SetHeader("#splitline{   Thermal-FIST (GCE SHM)}{   Nuclear excitation particle list}")
    leg_fist.AddEntry(graph_yields_li4_fist, "(^{4}Li + ^{4}#bar{Li})/2", "l") #           Run 3", "l") #}{#it{T} = 156.6 MeV, #it{V} = 4459 fm^{3}}", "l")
    leg_fist.AddEntry(graph_yields_li4_fist_gs, "(^{4}Li + ^{4}#bar{Li})/2 g.s.", "l") #,   Run 3", "l") #}{#it{T} = 156.6 MeV, #it{V} = 4459 fm^{3}}", "l")
    leg_fist.AddEntry(graph_yields_he4_fist, "(^{4}He + ^{4}#bar{He})/2", "l") #        Run 2", "l") #}{#it{T} = 156.4 MeV, #it{V} = 4233 fm^{3}}", "l")
    #leg_fist.AddEntry(graph_yields_he4_fist_1050, "(^{4}He + ^{4}#bar{He})/2        Run 3", "l") #}{#it{T} = 158.8 MeV, #it{V} = 1804 fm^{3}}", "l")
    leg_fist.Draw()
 
    nucleus_yield_010_with_syst = correct_yield('Both', '0-10%', raw_yield=ufloat(YIELD['Both']['0-10%'].n, 39))
    print(f"Corrected yield for Both in 0-10%: {nucleus_yield_010_with_syst:.2e} #pm {nucleus_yield_010_with_syst.s:.2e} (stat + syst)")
    # print(f"Corrected yield for Both in 0-10%: {nucleus_yield_010:.2e}")
    print(f"Corrected yield for Both in 10-50%: {nucleus_yield_1050:.2e} #pm {nucleus_yield_1050.s:.2e} (stat) #pm {nucleus_yield_1050_syst.s:.2e} (syst)")
    print(f"95% CL upper limit for Both in 0-10%: {nucleus_upper_limit_010:.2e}")
    # print(f"Corrected yield for (^{{4}}He + ^{{4}}#bar{{He}})/2 in 0-10%: {correct_yield('Both', '0-10%'):.2e}")
    # print(f"Thermal-FIST prediction for ^{{4}}Li: 8.2004e-6")
    print(f"Upper limit to the ^{{4}}Li/^{{4}}He ratio in 0-10%: {(nucleus_upper_limit_010.n / 1.0e-6):.2e}")
    
    ratio_1050 = nucleus_yield_1050 / THERMAL_FIST_PREDICTIONS['He4']['10-50%']
    ratio_1050_syst = nucleus_yield_1050_syst / THERMAL_FIST_PREDICTIONS['He4']['10-50%']
    print(f"^{{4}}Li/^{{4}}He ratio in 10-50%: {ratio_1050.n:.2e} #pm {ratio_1050.s:.2e} (stat) #pm {ratio_1050_syst.s:.2e} (syst)")
    
    hframe_yields_centralities.Draw()
    
    arrow_010.Draw()
    for graph in [graph_yields_1050_syst, graph_yields_alpha_010_syst, graph_yields_h4l_010_syst, graph_yields_he4l_010_syst]:
        graph.Draw("E2 SAME")
    for graph in [graph_yields_010, graph_yields_1050, graph_yields_alpha_010, graph_yields_h4l_010, graph_yields_he4l_010,
                  graph_yields_li4_fist, graph_yields_he4_fist, graph_yields_li4_fist_gs,
                  graph_yields_li4_fist_1050, graph_yields_li4_fist_gs_1050, graph_yields_he4_fist_1050]:
        graph.Draw("P SAME")
    
    leg.Draw()
    leg_fist.Draw()
    text.Draw('same')
    outfile.cd()
    canvas_yields_centralities.Write()
    canvas_yields_centralities.SaveAs("figures/corrected_yields_centralities.pdf")
    
def compute_antimatter_matter_ratio() -> None:
    
    yield_antimatter_1050 = correct_yield('Antimatter', '10-50%')
    yield_matter_1050 = correct_yield('Matter', '10-50%')
    
    raw_yield_antimatter_1050_syst = ufloat(332.82, 29.5) # errors: lambda
    raw_yield_matter_1050_syst = ufloat(547.91, 40.5)
    yield_antimatter_1050_syst = correct_yield('Antimatter', '10-50%', raw_yield=raw_yield_antimatter_1050_syst)
    yield_matter_1050_syst = correct_yield('Matter', '10-50%', raw_yield=raw_yield_matter_1050_syst)
    
    ratio_1050 = yield_antimatter_1050 / yield_matter_1050
    ratio_1050_syst = yield_antimatter_1050_syst / yield_matter_1050_syst
    
    syst_unc_li4_width = 0.0015 # uncertainty to the uncorrected ratio
    syst_unc_R = 0. # uncertainty to the uncorrected ratio
    syst_unc_cut_var = 0.047 # uncertainty to the uncorrected ratio
    syst_unc_total = np.sqrt(syst_unc_li4_width**2 + syst_unc_R**2 + syst_unc_cut_var**2)
    corrected_syst_unc = syst_unc_total * EFFICIENCY['10-50%']['Matter'] / EFFICIENCY['10-50%']['Antimatter']
    ratio_1050_syst_unc = np.sqrt(ratio_1050_syst.s**2 + corrected_syst_unc**2)
    print(f"\nAntimatter/Matter ratio in 10-50%: {ratio_1050.n:.2e} #pm {ratio_1050.s:.2e} (stat) #pm {ratio_1050_syst_unc:.2e} (syst)")
 
def draw_yields(outfile: TFile) -> None:
    """
    Orchestrator: produces all yield-related plots.
    """
    ### draw_yield_vs_multiplicity(outfile)
    ### draw_ratio_vs_multiplicity(outfile)
    draw_yields_vs_centrality(outfile)
    compute_antimatter_matter_ratio()
    
def print_corrected_yields():
    
    print('----------------------------------------------------------------')
    for sign in ['Matter', 'Antimatter', 'Both', 'Both_GS']:
        for centrality in ['0-10%', '10-30%', '30-50%', '50-80%', '0-50%', '10-50%']:
            if centrality not in YIELD[sign]:
                continue
            corrected_yield = correct_yield(sign, centrality, raw_yield=YIELD[sign][centrality])
            corrected_yield_syst = correct_yield(sign, centrality, raw_yield=YIELD_SYST[sign][centrality]) if sign in YIELD_SYST and centrality in YIELD_SYST[sign] else None
            if corrected_yield_syst is not None:
                print(f"Corrected yield for {sign} in {centrality}: {corrected_yield.n:.2e} ± {corrected_yield.s:.2e} (stat) ± {corrected_yield_syst.s:.2e} (syst)")
            else:
                print(f"Corrected yield for {sign} in {centrality}: {corrected_yield.n:.2e} ± {corrected_yield.s:.2e} (stat)")   
    print('----------------------------------------------------------------')
 
# ============================================================================
# Main
# ============================================================================
 
if __name__ == "__main__":
    
    set_alice_global_style()
 
    print_corrected_yields()
            
    outfile = TFile.Open("output/corrected_yields.root", "RECREATE")
    draw_yields(outfile)
    outfile.Close()