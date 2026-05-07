
import numpy as np
from uncertainties import ufloat
from ROOT import TF1, TCanvas, TGraphErrors, TLegend, gStyle, TGraphAsymmErrors, TGraphErrors, TPaveText
from ROOT import kRed, kBlack, kBlue, kGreen, kMagenta, kCyan

import sys
sys.path.append('..')
from utils.particles import ParticleMasses
from torchic.utils.root import set_root_object

def pt_to_mt(pt:float, particle:str):
    return ufloat(np.sqrt( pt**2 + ParticleMasses[particle]**2 ), 0)

PARAMETERISATION = {
    '0-10%': {
        'a': ufloat(3.59, 0.05),
        'b': ufloat(2.74, 0.04),
        'c': ufloat(-1.91, 0.09),
    },
    '10-30%': {
        'a': ufloat(3.38, 0.02),
        'b': ufloat(2.21, 0.08),
        'c': ufloat(-3.16, 0.14),
    },
    '30-50%': {
        'a': ufloat(2.50, 0.03),
        'b': ufloat(1.28, 0.03),
        'c': ufloat(-2.13, 0.14),
    },
}

MEASURED_POINTS_PT = {
    'Pr': {
        '0-10%': ufloat(0.9584, 0.3755),
        '10-30%': ufloat(0.9274, 0.3569),
        '30-50%': ufloat(0.8532, 0.3152),
    },
    'He': {
        '0-10%': ufloat(2.899, 1.8083),
        '10-30%': ufloat(2.801, 1.032),
        '30-50%': ufloat(2.578, 0.9818),
    },
}

MEASURED_POINTS = {
    'Pr': {
        centrality: pt_to_mt(MEASURED_POINTS_PT['Pr'][centrality].n, 'Pr')
        for centrality in MEASURED_POINTS_PT['Pr']
    },
    'He': {
        centrality: pt_to_mt(MEASURED_POINTS_PT['He'][centrality].n, 'He')
        for centrality in MEASURED_POINTS_PT['He']
    },
}


def draw_source_radius():

    funcs = {}
    funcs_plus_error = {}
    bands = {}
    points = { 'Pr': {}, 'He': {} }

    XMIN, XMAX = 0.9, 7.5
    NPOINTS = 100
    X = np.linspace(XMIN, XMAX, NPOINTS)

    for centrality in ['0-10%', '10-30%', '30-50%']:

        funcs[centrality] = TF1(f'func_{centrality}', f'[0] + [1]*x^[2]', XMIN, XMAX)
        funcs[centrality].SetParameters(
            PARAMETERISATION[centrality]['a'].n,
            PARAMETERISATION[centrality]['b'].n,
            PARAMETERISATION[centrality]['c'].n
        )
        funcs[centrality].SetParErrors(
            np.array([
                PARAMETERISATION[centrality]['a'].s,
                PARAMETERISATION[centrality]['b'].s,
                PARAMETERISATION[centrality]['c'].s
            ]))
        
        funcs_plus_error[centrality] = TF1(f'func_plus_error{centrality}', f'[0] + [1]*x^[2]', XMIN, XMAX)
        funcs_plus_error[centrality].SetParameters(
            PARAMETERISATION[centrality]['a'].n + PARAMETERISATION[centrality]['a'].s,
            PARAMETERISATION[centrality]['b'].n + PARAMETERISATION[centrality]['b'].s,
            PARAMETERISATION[centrality]['c'].n + PARAMETERISATION[centrality]['c'].s
        )
        
        for particle in ['Pr', 'He']:
            points[particle][centrality] = TGraphErrors(1)
            points[particle][centrality].SetPoint(0, MEASURED_POINTS[particle][centrality].n, funcs[centrality].Eval(MEASURED_POINTS[particle][centrality].n))
        print(f'Measured {centrality}: {MEASURED_POINTS[particle][centrality].n}, Rp: {points['Pr'][centrality].GetY()[0]}, Rhe3: {points['He'][centrality].GetY()[0]}')

        bands[centrality] = TGraphAsymmErrors(NPOINTS)
        for ix, x in enumerate(X):
            #dfa = 1
            #dfb = x**PARAMETERISATION[centrality]['c'].n if x > 0 else 0
            #dfc = np.log(x) * PARAMETERISATION[centrality]['b'].n * x**(PARAMETERISATION[centrality]['c'].n - 1) if x > 0 else 0
            #err = np.sqrt((dfa * PARAMETERISATION[centrality]['a'].s)**2 + 
            #              (dfb * PARAMETERISATION[centrality]['b'].s)**2 + 
            #              (dfc * PARAMETERISATION[centrality]['c'].s)**2)

            err = funcs_plus_error[centrality].Eval(x) - funcs[centrality].Eval(x)

            bands[centrality].SetPoint(ix, x, funcs[centrality].Eval(x))
            bands[centrality].SetPointError(ix, 0., 0., err, err)



    canvas = TCanvas('canvas', 'Source Radius Parameterisation', 800, 600)
    gStyle.SetOptStat(0)
    hframe = canvas.DrawFrame(XMIN, 2.2, XMAX, 7.5, ';#LT#it{m}_{T}#GT (GeV/#it{c}^{2});#it{R} (fm)')

    watermark = TPaveText(0.15, 0.76, 0.45, 0.88, 'NDC')
    watermark.SetFillColor(0)
    watermark.SetBorderSize(0)
    watermark.AddText('This work')
    watermark.AddText('#bf{ALICE Run 3}')
    watermark.AddText('#bf{Pb-Pb  #it{#sqrt{s_{NN}}} = 5.36 TeV}')
    
    legend = TLegend(0.5, 0.68, 0.88, 0.88)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetNColumns(3)
    legend.SetTextSize(0.03)

    results_panel = TPaveText(0.38, 0.48, 0.88, 0.64, 'NDC')
    results_panel.SetFillColor(0)
    results_panel.SetBorderSize(0)
    results_panel.SetTextSize(0.038)
    results_panel.AddText('       #bf{#bf{#it{R}_{p} (fm)}  #bf{#it{R}_{^{3}He} (fm)}}')

    for color, centrality in zip([kRed, kBlue, kGreen], ['0-10%', '10-30%', '30-50%']):

        bands[centrality].SetFillColor(color)
        bands[centrality].SetFillStyle(3001)
        bands[centrality].Draw('3 same')

        set_root_object(points['Pr'][centrality], marker_color=color, marker_size=2, marker_style=20)
        points['Pr'][centrality].Draw('p same')
        set_root_object(points['He'][centrality], marker_color=color, marker_size=3, marker_style=33)
        points['He'][centrality].Draw('p same')
        
        legend.AddEntry(bands[centrality], centrality, 'f')
        legend.AddEntry(points['Pr'][centrality], f'#it{{R}}_{{p}}', 'p')
        legend.AddEntry(points['He'][centrality], f'#it{{R}}_{{^{{3}}He}}', 'p')

        results_panel.AddText(f'#bf{{{centrality}:}}   #bf{{{points['Pr'][centrality].GetY()[0]:.2f}}}    #bf{{{points['He'][centrality].GetY()[0]:.2f}}}')

        print(f'Centrality {centrality}: Rp = {points["Pr"][centrality].GetY()[0]:.2f} fm, Rhe3 = {points["He"][centrality].GetY()[0]:.2f} fm, Rsource = {np.sqrt((points["Pr"][centrality].GetY()[0])**2 + (points["He"][centrality].GetY()[0])**2):.2f} fm')

    watermark.Draw()
    results_panel.Draw()
    legend.Draw('same')
    canvas.SaveAs('figures/source_radius_parameterisation.pdf')

if __name__ == '__main__':

    draw_source_radius()
    print("Source radius parameterisation plot saved as 'source_radius_parameterisation.pdf'.")
