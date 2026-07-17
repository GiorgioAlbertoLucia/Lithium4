from itertools import product

import numpy as np
from uncertainties import ufloat
from ROOT import TF1, TCanvas, TGraphErrors, TLegend, gStyle, TGraphAsymmErrors, TGraphErrors, TPaveText, TFile
from ROOT import kRed, kBlack, kBlue, kGreen, kMagenta, kCyan

import sys
sys.path.append('..')
from utils.particles import ParticleMasses
from torchic.utils.root import set_root_object, init_legend
from torchic.utils.colors import get_color

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

NCH_RUN2 = {
    '0-10%':  ufloat(1764, 24),
    '10-30%': ufloat(983, 19),
    '30-50%': ufloat(415, 10),
}

NCH_RUN3 = {
    '0-10%':  ufloat(1858, 34),
    '10-30%': ufloat(1051, 21),
    '30-50%': ufloat(455, 12),
    '50-80%': ufloat(123, 5),
}

# Reweighted for events with at least a 3He in the final state
NCH_RUN3_PHE3 = {
    '0-10%':  ufloat(1871, 35),
    '10-30%': ufloat(1084, 22),
    '30-50%': ufloat(480, 13),
    '50-80%': ufloat(169, 7),
    '10-50%': ufloat(890, 16),
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

MEASURED_POINTS_MT = {
    'Pr': {
        '0-10%':  ufloat(1.74, 0.),
        '10-30%': ufloat(1.74, 0.),
        '30-50%': ufloat(1.69, 0.),
        '50-80%': ufloat(1.59, 0.),
        '10-50%': ufloat(1.72, 0.),
    },
    'He': {
        '0-10%':  ufloat(4.97, 0.),
        '10-30%': ufloat(4.83, 0.),
        '30-50%': ufloat(4.60, 0.),
        '50-80%': ufloat(4.07, 0.),
        '10-50%': ufloat(4.60, 0.),
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
    
    outfile = TFile('source_radius_parameterisation.root', 'RECREATE')
    
    for centrality_to_draw in ['0-10%', '10-30%', '30-50%', '50-80%', '10-50%']:
        for icentrality, centrality in enumerate(['0-10%', '10-30%', '30-50%']):
            
            a, b, c = PARAMETERISATION[centrality]['a'].n, PARAMETERISATION[centrality]['b'].n, PARAMETERISATION[centrality]['c'].n
            a_err, b_err, c_err = PARAMETERISATION[centrality]['a'].s, PARAMETERISATION[centrality]['b'].s, PARAMETERISATION[centrality]['c'].s
            combos = np.array(list(product([a-a_err, a+a_err], [b-b_err, b+b_err], [c-c_err, c+c_err])))

            pr_point = TGraphAsymmErrors(1)
            he_point = TGraphAsymmErrors(1)
            
            funcs[centrality] = TF1(f'func_{centrality}', f'[0] + [1]*x^[2]', XMIN, XMAX)
            funcs[centrality].SetParameters(a, b, c)
            
            x_pr = MEASURED_POINTS_MT['Pr'][centrality_to_draw].n
            y_pr = funcs[centrality].Eval(x_pr)
            max_y_pr = np.max(combos[:,0] + combos[:,1] * x_pr**combos[:,2])
            min_y_pr = np.min(combos[:,0] + combos[:,1] * x_pr**combos[:,2])
        
            pr_point.SetPoint(0, x_pr, y_pr)
            pr_point.SetPointError(0, 0, 0, y_pr - min_y_pr, max_y_pr - y_pr)
            points['Pr'][centrality] = pr_point
            
            x_he = MEASURED_POINTS_MT['He'][centrality_to_draw].n
            y_he = funcs[centrality].Eval(x_he)
            max_y_he = np.max(combos[:,0] + combos[:,1] * x_he**combos[:,2])
            min_y_he = np.min(combos[:,0] + combos[:,1] * x_he**combos[:,2])
        
            he_point.SetPoint(0, x_he, y_he)
            he_point.SetPointError(0, 0, 0, y_he - min_y_he, max_y_he - y_he)
            points['He'][centrality] = he_point
            
            bands[centrality] = TGraphAsymmErrors(NPOINTS)
            for i, x in enumerate(X):
                y = funcs[centrality].Eval(x)
                max_y = np.max(combos[:,0] + combos[:,1] * x**combos[:,2])
                min_y = np.min(combos[:,0] + combos[:,1] * x**combos[:,2])
                bands[centrality].SetPoint(i, x, y)
                bands[centrality].SetPointError(i, 0, 0, y - min_y, max_y - y)
        

        canvas = TCanvas(f'canvas_{centrality_to_draw}', 'Source Radius Parameterisation', 800, 600)
        gStyle.SetOptStat(0)
        hframe = canvas.DrawFrame(XMIN, 2.2, XMAX, 7.5, f'#LT#it{{m}}_{{T}}#GT for {centrality_to_draw};#LT#it{{m}}_{{T}}#GT (GeV/#it{{c}}^{{2}});#it{{R}} (fm)')

        watermark = TPaveText(0.15, 0.76, 0.45, 0.88, 'NDC')
        watermark.SetFillColor(0)
        watermark.SetBorderSize(0)
        watermark.AddText('This work')
        watermark.AddText('#bf{ALICE Run 3}')
        watermark.AddText('#bf{Pb-Pb  #it{#sqrt{s_{NN}}} = 5.36 TeV}')
        
        legend = init_legend(0.5, 0.68, 0.88, 0.88, text_size=0.03, n_columns=3)
        
        results_panel = TPaveText(0.38, 0.48, 0.88, 0.64, 'NDC')
        results_panel.SetFillColor(0)
        results_panel.SetBorderSize(0)
        results_panel.SetTextSize(0.038)
        results_panel.AddText('       #bf{#bf{#it{R}_{p} (fm)}  #bf{#it{R}_{^{3}He} (fm)}}')

        for icentrality, centrality in enumerate(['0-10%', '10-30%', '30-50%']):

            #bands[centrality].SetFillColor(color)
            #bands[centrality].SetFillStyle(3001)
            #bands[centrality].Draw('3 same')
            
            set_root_object(funcs[centrality], line_color=get_color(icentrality), line_width=2)
            funcs[centrality].Draw('same')

            set_root_object(points['Pr'][centrality], marker_color=get_color(icentrality), marker_size=1.5, marker_style=20)
            points['Pr'][centrality].Draw('p same')
            set_root_object(points['He'][centrality], marker_color=get_color(icentrality), marker_size=2, marker_style=33)
            points['He'][centrality].Draw('p same')
            
            set_root_object(bands[centrality], fill_color=get_color(icentrality), line_color=0, fill_style=3001, fill_color_alpha=(get_color(icentrality), 0.3))
            bands[centrality].Draw('3 same')
            
            #legend.AddEntry(bands[centrality], centrality, 'f')
            legend.AddEntry(points['Pr'][centrality], f'#it{{R}}_{{p}}', 'p')
            legend.AddEntry(points['He'][centrality], f'#it{{R}}_{{^{{3}}He}}', 'p')

            results_panel.AddText(f'#bf{{{centrality}:}}   #bf{{{points['Pr'][centrality].GetY()[0]:.2f}}}    #bf{{{points['He'][centrality].GetY()[0]:.2f}}}')
            print(f'Centrality {centrality}: Rp = {points["Pr"][centrality].GetY()[0]:.2f} fm, Rhe3 = {points["He"][centrality].GetY()[0]:.2f} fm, Rsource = {np.sqrt((points["Pr"][centrality].GetY()[0])**2 + (points["He"][centrality].GetY()[0])**2):.2f} fm')

        watermark.Draw()
        results_panel.Draw()
        legend.Draw('same')
        canvas.SaveAs(f'source_radius_parameterisation_{centrality_to_draw}.pdf')
        outfile.cd()
        canvas.Write(f'canvas_{centrality_to_draw}')
    

        Nch_pr_graph = TGraphAsymmErrors(3)
        Nch_he_graph = TGraphAsymmErrors(3)
        for icentrality, centrality in enumerate(['0-10%', '10-30%', '30-50%']):
            
            nch = NCH_RUN2[centrality]
            
            Nch_pr_graph.SetPoint(icentrality, nch.n**(1/3), points['Pr'][centrality].GetY()[0])
            Nch_pr_graph.SetPointError(icentrality, nch.s * (1/3) * nch.n**(-2/3), nch.s * (1/3) * nch.n**(-2/3), 
                                    points['Pr'][centrality].GetErrorYlow(0), points['Pr'][centrality].GetErrorYhigh(0))
            
            Nch_he_graph.SetPoint(icentrality, nch.n**(1/3), points['He'][centrality].GetY()[0])
            Nch_he_graph.SetPointError(icentrality, nch.s * (1/3) * nch.n**(-2/3), nch.s * (1/3) * nch.n**(-2/3), 
                                    points['He'][centrality].GetErrorYlow(0), points['He'][centrality].GetErrorYhigh(0))
                
        canvas_nch = TCanvas(f'canvas_nch_{centrality_to_draw}', 'Source Radius Parameterisation Nch', 800, 600)
        gStyle.SetOptStat(0)
        hframe_nch = canvas_nch.DrawFrame(3, 1.2, 15, 6.5, f'{centrality_to_draw};#LTd#it{{N}}_{{ch}}/d#it{{#eta}} #GT^{{1/3}}_{{|#it{{#eta}}|<0.5}};#it{{R}} (fm)')
        legend_nch = init_legend(0.54, 0.18, 0.88, 0.34, text_size=0.03, n_columns=2)
        NCH_CUBEROOT_RANGE = np.linspace(3, 15, 100)
        
        funcs_nch = {}
        run3_points = {'Pr': TGraphAsymmErrors(4), 'He': TGraphAsymmErrors(4) }
        bands = {'Pr': TGraphAsymmErrors(NCH_CUBEROOT_RANGE.size), 'He': TGraphAsymmErrors(NCH_CUBEROOT_RANGE.size)}
        label = {'Pr': '#it{R}_{p}', 'He': '#it{R}_{^{3}He}'}
        
        text_nch = TPaveText(0.15, 0.66, 0.45, 0.88, 'NDC')
        text_nch.SetFillColor(0)
        text_nch.SetBorderSize(0)
        text_nch.SetTextSize(0.038)
        text_nch.SetTextFont(42)
        text_nch.AddText('ALICE')
        text_nch.AddText('Pb-Pb  #it{#sqrt{s_{NN}}} = 5.36 TeV')
        
        source_radius, upper_source_radius, lower_source_radius = 0, 0, 0
        for iparticle, particle in enumerate(['Pr', 'He']):
            
            graph = Nch_pr_graph if particle == 'Pr' else Nch_he_graph
            funcs_nch[particle] = TF1(f'func_nch_{particle}', f'[0] + [1]*x', 3, 15)
            fit_result = graph.Fit(funcs_nch[particle], 'RMS+')
            covariance_matrix = fit_result.GetCovarianceMatrix()
            covariance = covariance_matrix(0,1)
            
            print('Parameterisation for particle', particle)
            covariance_matrix.Print()
            
            nch = NCH_RUN3_PHE3[centrality_to_draw]
            radius = funcs_nch[particle].Eval(nch.n**(1/3))
            
            #a, b = funcs_nch[particle].GetParameter(0), funcs_nch[particle].GetParameter(1) 
            #a_err, b_err = funcs_nch[particle].GetParError(0), funcs_nch[particle].GetParError(1)
            #combos = np.array(list(product([a-a_err, a+a_err], [b-b_err, b+b_err])))
            #max_y = np.max(combos[:,0] + combos[:,1] * nch.n**(1/3))
            #min_y = np.min(combos[:,0] + combos[:,1] * nch.n**(1/3))
            #
            #run3_points[particle].SetPoint(icentrality, nch.n**(1/3), radius)
            #run3_points[particle].SetPointError(icentrality, nch.s * (1/3) * nch.n**(-2/3), nch.s * (1/3) * nch.n**(-2/3), radius - min_y, max_y - radius)
            #print(f'Centrality {centrality}: Nch = {nch.n:.0f}, Radius = {radius:.2f} + {max_y - radius:.2f} / - {radius - min_y:.2f} fm')
            
            radius_err = np.sqrt( (funcs_nch[particle].GetParError(0))**2 + (nch.n**(1/3) * funcs_nch[particle].GetParError(1))**2 - 2 * nch.n**(1/3) * covariance )
            run3_points[particle].SetPoint(0, nch.n**(1/3), radius)
            run3_points[particle].SetPointError(0, nch.s * (1/3) * nch.n**(-2/3), nch.s * (1/3) * nch.n**(-2/3), radius_err, radius_err)
            print(f'Centrality {centrality_to_draw}: Nch = {nch.n:.0f}, (Nch)^{{1/3}} = {nch.n**(1/3):.2f}, Radius = {radius:.2f} ± {radius_err:.2f} fm')
            particle_label = 'p' if particle == 'Pr' else '^{3}He'
            text_nch.AddText(f'#it{{R}}_{{{particle_label}}} ({centrality_to_draw}) = ({radius:.2f} #pm {radius_err:.2f}) fm')
            source_radius += radius**2
            upper_source_radius += (radius + radius_err)**2
            lower_source_radius += (radius - radius_err)**2
                
            for i, nch_cuberoot in enumerate(NCH_CUBEROOT_RANGE):
                y = funcs_nch[particle].Eval(nch_cuberoot)
                
                #a, b = funcs_nch[particle].GetParameter(0), funcs_nch[particle].GetParameter(1) 
                #a_err, b_err = funcs_nch[particle].GetParError(0), funcs_nch[particle].GetParError(1)
                #combos = np.array(list(product([a-a_err, a+a_err], [b-b_err, b+b_err])))
                #max_y = np.max(combos[:,0] + combos[:,1] * nch)
                #min_y = np.min(combos[:,0] + combos[:,1] * nch)
                #
                #bands[particle].SetPoint(i, nch, y)
                #bands[particle].SetPointError(i, 0, 0, y - min_y, max_y - y)
                
                y_err = np.sqrt( (funcs_nch[particle].GetParError(0))**2 + (nch_cuberoot * funcs_nch[particle].GetParError(1))**2 - 2 * nch_cuberoot * covariance )
                bands[particle].SetPoint(i, nch_cuberoot, y)
                bands[particle].SetPointError(i, 0, 0, y_err, y_err)
            
            set_root_object(graph, marker_color=get_color(iparticle), marker_size=(1.5 if particle == 'Pr' else 2), 
                            marker_style=20 if particle == 'Pr' else 33, line_color=get_color(iparticle))
            set_root_object(funcs_nch[particle], line_color=get_color(iparticle), line_width=2)
            set_root_object(run3_points[particle], marker_color=get_color(iparticle), marker_size=(1.5 if particle == 'Pr' else 2), 
                            marker_style=24 if particle == 'Pr' else 27, line_color=get_color(iparticle))
            set_root_object(bands[particle], fill_color=get_color(iparticle), line_color=0, fill_style=3001, fill_color_alpha=(get_color(iparticle), 0.3))
            
            graph.Draw('p same')
            funcs_nch[particle].Draw('same')
            run3_points[particle].Draw('p same')
            bands[particle].Draw('3 same')
            legend_nch.AddEntry(graph, f'Run 2: {label[particle]}', 'p')
            legend_nch.AddEntry(run3_points[particle], f'Run 3: {label[particle]}', 'p')
        
        source_radius = np.sqrt(source_radius)
        upper_source_radius = np.sqrt(upper_source_radius)
        lower_source_radius = np.sqrt(lower_source_radius)
        upper_source_radius_err = upper_source_radius - source_radius
        lower_source_radius_err = source_radius - lower_source_radius
        print(f'Centrality {centrality_to_draw}: Source radius = {source_radius:.2f} + {upper_source_radius_err:.2f} / - {lower_source_radius_err:.2f} fm')
        #input("Press Enter to continue...")
        
        text_nch.Draw()
        legend_nch.Draw('same')
        canvas_nch.SaveAs(f'source_radius_parameterisation_nch_{centrality_to_draw}.pdf')


if __name__ == '__main__':

    draw_source_radius()
    print("Source radius parameterisation plot saved as 'source_radius_parameterisation.pdf'.")
