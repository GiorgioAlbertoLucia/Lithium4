
import numpy as np
from uncertainties import ufloat
from ROOT import TF1, TCanvas, TGraphErrors, TLegend, gStyle, TGraphAsymmErrors, TGraphErrors, TPaveText
from ROOT import kRed, kBlack, kBlue, kGreen, kMagenta, kCyan

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

MEASURED_POINTS = {
    '0-10%': ufloat(1.416, 0.301),
    '10-30%': ufloat(1.388, 0.299),
    '30-50%': ufloat(1.333, 0.282),
}

def draw_source_radius():

    funcs = {}
    bands = {}
    points = {}

    XMIN, XMAX = 0.9, 3.5
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
        
        points[centrality] = TGraphErrors(1)
        points[centrality].SetPoint(0, MEASURED_POINTS[centrality].n, funcs[centrality].Eval(MEASURED_POINTS[centrality].n))
        print(f'Measured {centrality}: {MEASURED_POINTS[centrality].n}, Rsource: {points[centrality].GetY()[0]}')

        bands[centrality] = TGraphAsymmErrors(NPOINTS)
        for ix, x in enumerate(X):
            dfa = 1
            dfb = x**PARAMETERISATION[centrality]['c'].n if x > 0 else 0
            dfc = np.log(x) * PARAMETERISATION[centrality]['b'].n * x**(PARAMETERISATION[centrality]['c'].n - 1) if x > 0 else 0
            err = np.sqrt((dfa * PARAMETERISATION[centrality]['a'].s)**2 + 
                          (dfb * PARAMETERISATION[centrality]['b'].s)**2 + 
                          (dfc * PARAMETERISATION[centrality]['c'].s)**2)
            bands[centrality].SetPoint(ix, x, funcs[centrality].Eval(x))
            bands[centrality].SetPointError(ix, 0., 0., err, err)



    canvas = TCanvas('canvas', 'Source Radius Parameterisation', 800, 600)
    gStyle.SetOptStat(0)
    hframe = canvas.DrawFrame(XMIN, 2.2, XMAX, 7.5, ';#LT#it{m}_{T}#GT (GeV/#it{c}^{2});#it{R}_{source} (fm)')

    watermark = TPaveText(0.15, 0.76, 0.45, 0.88, 'NDC')
    watermark.SetFillColor(0)
    watermark.SetBorderSize(0)
    watermark.AddText('This work')
    watermark.AddText('#bf{ALICE Run 3}')
    watermark.AddText('#bf{Pb-Pb  #it{#sqrt{s_{NN}}} = 5.36 TeV}')
    
    legend = TLegend(0.5, 0.6, 0.88, 0.88)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetNColumns(2)
    legend.SetTextSize(0.03)

    results_panel = TPaveText(0.5, 0.44, 0.88, 0.6, 'NDC')
    results_panel.SetFillColor(0)
    results_panel.SetBorderSize(0)
    results_panel.AddText('         #bf{#LT#it{m}_{T}#GT (GeV/#it{c}^{2})}  #bf{#it{R}_{source} (fm)}')

    for color, centrality in zip([kRed, kBlue, kGreen], ['0-10%', '10-30%', '30-50%']):

        bands[centrality].SetFillColor(color)
        bands[centrality].SetFillStyle(3001)
        bands[centrality].Draw('3 same')

        points[centrality].SetMarkerColor(color+3)
        points[centrality].SetMarkerSize(2)
        points[centrality].SetMarkerStyle(33)
        points[centrality].Draw('p same')

        legend.AddEntry(bands[centrality], centrality, 'f')
        legend.AddEntry(points[centrality], f'Measured {centrality}', 'p')

        results_panel.AddText(f'#bf{{{centrality}:}}    #bf{{{MEASURED_POINTS[centrality].n:.2f}}}     #bf{{{points[centrality].GetY()[0]:.2f}}}')

    watermark.Draw()
    results_panel.Draw()
    legend.Draw('same')
    canvas.SaveAs('figures/source_radius_parameterisation.pdf')

if __name__ == '__main__':

    draw_source_radius()
    print("Source radius parameterisation plot saved as 'source_radius_parameterisation.pdf'.")
