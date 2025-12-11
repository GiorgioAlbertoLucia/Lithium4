'''
    Study of the dca distribution for p-He3 pairs from data
'''

import numpy as np
import pandas as pd
import sys
from enum import Enum

from ROOT import TFile, TDirectory, TH2F, TCanvas, TLegend, TPaveText,TGraphErrors, \
                 RooRealVar, RooCrystalBall, RooGaussian, RooDataHist, RooFFTConvPdf, RooAddPdf

from torchic.core.histogram import load_hist
from torchic.core.graph import create_graph
from torchic.utils.terminal_colors import TerminalColors as tc
from torchic.utils.root import set_root_object

def get_dca_one_sigma(pt:float, particle:str):

    if particle == 'He':
        return 0.0032 * np.exp(- np.abs(pt) * 0.5206) + 0.0012
    elif particle == 'Pr':
        return 0.0118 * np.exp(- np.abs(pt) * 0.6889) + 0.0017
    
    return 0.

def build_model(x, suffix:str='', **kwargs):

    pars = {
        'mean': kwargs.get('mean', RooRealVar(f'mean{suffix}', f'#mu^{suffix}', 0., -0.01, 0.01)),
        'sigma': RooRealVar(f'sigma{suffix}', f'#sigma^{suffix}', 0.01, 1.e-4, 0.05),
        'aL': RooRealVar(f'aL{suffix}', f'#alpha_{{L}}^{suffix}', 1.2, 0.1, 10.),
        'nL': RooRealVar(f'nL{suffix}', f'n_{{L}}^{suffix}', 3., 0.5, 30.),
        'aR': RooRealVar(f'aR{suffix}', f'#alpha_{{R}}^{suffix}', 1.2, 0.1, 10.),
        'nR': RooRealVar(f'nR{suffix}', f'n_{{R}}^{suffix}', 3., 0.5, 30.)

    }
    pdf = RooCrystalBall(f'pdf_{suffix}', f'pdf_{suffix}', x, *pars.values())
    return pdf, pars

def _prepare_template(x, pdf, params, h_mc, gaussian_smearing, bin_outdir #, gaussian_smearing_pars, sigma_core
                      ):

    dh = RooDataHist(f'tmp', '', [x], Import=h_mc)
    pdf.fitTo(dh, PrintLevel=1, Save=True)

    for param_name, param in params.items():
        param.setConstant(True)
    
    frame = x.frame()
    dh.plotOn(frame)
    pdf.plotOn(frame, Name=pdf.GetName(), LineColor=2)
    pdf.paramOn(frame, ShowConstants=True)
    
    pdf_convoluted = RooFFTConvPdf(pdf.GetName()+'_convoluted', pdf.GetTitle()+'_convoluted', x, 
                                   pdf, gaussian_smearing, ipOrder=2)
    pdf_convoluted.plotOn(frame, Name=pdf_convoluted.GetName(), LineColor=3, LineStyle='--')

    legend = TLegend(0.11, 0.6, 0.4, 0.89)
    legend.SetBorderSize(0)

    dh.plotOn(frame)
    legend.AddEntry(frame.findObject(pdf.GetName()), pdf.GetTitle(), "l")
    legend.AddEntry(frame.findObject(pdf_convoluted.GetName()), pdf_convoluted.GetTitle(), "l")

    canvas = TCanvas(pdf.GetName(), '')
    frame.Draw()
    legend.Draw('same')
    
    bin_outdir.cd()
    canvas.Write()

    del dh
    return pdf_convoluted

def _fit_slice(h2_data:TH2F, h2_mc:dict, 
               pdfs:dict, pdf_params:dict, gaussian_core, gaussian_core_pars:dict, dca:RooRealVar,
               pt_bin:int, particle:str,
               bin_outdir:TDirectory):
    
    pt = h2_data.GetXaxis().GetBinCenter(pt_bin)
    h_dca = h2_data.ProjectionY(f'h_dca_{pt:.2f}', pt_bin, pt_bin, 'e')
    if h_dca.GetEntries() < 100:
        return None

    pt_low_edge = h2_data.GetXaxis().GetBinLowEdge(pt_bin)
    pt_high_edge = h2_data.GetXaxis().GetBinLowEdge(pt_bin+1)


    dh = RooDataHist(f'dh_{pt:.2f}', 'dh', [dca], Import=h_dca)
    gaussian_core.fitTo(dh, PrintLevel=1, Save=True, Range='core')
    
    gaussian_smearing_pars = {
        'mean': RooRealVar('mean_smearning', f'#mu^{{smearing}}', 0., -0.05, 0.05),
        'sigma': RooRealVar('sigma_smearning', f'#sigma^{{smearing}}', 0.003, 1.e-5, 0.005),
        }
    gaussian_smearing = RooGaussian('gaussian_smearing', 'gaussian_smearing', 
                                    dca, *gaussian_smearing_pars.values())    
    gaussian_smearing_pars['mean'].setVal(gaussian_core_pars['mean'].getVal())
    gaussian_smearing_pars['mean'].setConstant(True)
        

    convoluted_pdfs, normalisations = {}, {}
    yield_mc = {}
    total_yield_mc = 0

    for iflag, flag in enumerate(h2_mc.keys()):

        h_dca_flag = h2_mc[flag].ProjectionY(f'h_mc', pt_bin, pt_bin, 'e')
        convoluted_pdfs[flag] = _prepare_template(dca, pdfs[flag], pdf_params[flag], h_dca_flag, 
                                                  gaussian_smearing, bin_outdir #, gaussian_smearing_pars,
                                                  #gaussian_core_pars['sigma'].getVal()
                                                  )
        yield_mc[flag] = h_dca_flag.GetEntries()
        total_yield_mc += h_dca_flag.GetEntries()
        del h_dca_flag

    for flag in h2_mc.keys():
        normalisations[flag] = RooRealVar(f'frac_{flag}', f'#it{{f}}_{{{flag}}}', yield_mc[flag]/total_yield_mc * h_dca.GetEntries(), 0., h_dca.GetEntries()) # extended
    
    if 'material' in yield_mc.keys() and 'weak_decay'in yield_mc.keys():# and particle == 'Pr':
        normalisations['material'].setVal( yield_mc['material'] / yield_mc['weak_decay'] * normalisations['weak_decay'].getVal() )
        normalisations['material'].setRange(0., 5. * yield_mc['material'] / yield_mc['weak_decay'] * normalisations['weak_decay'].getVal() )

    model = RooAddPdf('total_fit', 'total fit', list(convoluted_pdfs.values()), list(normalisations.values()))
    model.fitTo(dh, PrintLevel=1, Save=True)

    frame = dca.frame(Title=f'{pt_low_edge:.2f} < #it{{p}}_{{T}} < {pt_high_edge:.2f} GeV/#it{{c}}')

    legend = TLegend(0.11, 0.6, 0.4, 0.89)
    legend.SetBorderSize(0)

    dh.plotOn(frame)
    model.plotOn(frame, Name=model.GetName(), LineColor=2)
    model.paramOn(frame)
    legend.AddEntry(frame.findObject(model.GetName()), model.GetTitle(), "l")

    for iflag, flag in enumerate(h2_mc.keys()):
        model.plotOn(frame, Name=convoluted_pdfs[flag].GetName(), Components={convoluted_pdfs[flag]}, LineColor=3+iflag, LineStyle='--')
        legend.AddEntry(frame.findObject(convoluted_pdfs[flag].GetName()), convoluted_pdfs[flag].GetTitle(), "l")
    
    frame.addObject(legend)

    selection_window = (-3*get_dca_one_sigma(pt, particle), 3*get_dca_one_sigma(pt, particle))
    dca.setRange('integral_range', selection_window[0], selection_window[1])
    primaries_integral = convoluted_pdfs['primaries'].createIntegral(dca, dca, 'integral_range').getVal() * normalisations['primaries'].getVal()
    total_normalisation = 0
    for normalisation in normalisations.values():
        total_normalisation += normalisation.getVal() 
    total_integral = model.createIntegral(dca, dca, 'integral_range').getVal() * total_normalisation

    ifit_results = {
                    'pt': np.abs(pt),
                    'primaries_integral': primaries_integral,
                    'primaries_integral_mc': yield_mc['primaries'],
                    'total_integral': total_integral,
                    'primary_fraction': primaries_integral/total_integral,
                    'primary_fraction_mc': yield_mc['primaries'] / total_yield_mc if total_yield_mc > 0. else 0.
                    }

    return frame, ifit_results

def draw_results(h2_nsigma: TH2F, fit_results: pd.DataFrame, out_dir: TFile, sign: str):
    '''
        Draw the results of the fit
    '''
    fit_results['pt_err'] = h2_nsigma.GetXaxis().GetBinWidth(1)/2.
    
    graph_primaries = create_graph(fit_results, 'pt', 'primaries_integral', 'pt_err', 0., name=f'g_primaries_{sign}', title=';#it{p}_{T} (GeV/c); #it{N}_{p}')
    graph_primaries_mc = create_graph(fit_results, 'pt', 'primaries_integral_mc', 'pt_err', 0., name=f'g_primaries_mc_{sign}', title=';#it{p}_{T} (GeV/c); #it{N}_{p}')
    graph_total = create_graph(fit_results, 'pt', 'total_integral', 'pt_err', 0., name=f'g_total_{sign}', title=';#it{p}_{T} (GeV/c); #it{N}_{total}')
    graph_primary_fraction = create_graph(fit_results, 'pt', 'primary_fraction', 'pt_err', 0., name=f'g_primary_fraction_{sign}', title=';#it{p}_{T} (GeV/c); f_{p}')
    graph_primary_fraction_mc = create_graph(fit_results, 'pt', 'primary_fraction_mc', 'pt_err', 0., name=f'g_primary_fraction_mc_{sign}', title='MC;#it{p}_{T} (GeV/c); f_{p}')

    for graph in [graph_primaries, graph_total, graph_primary_fraction, graph_primary_fraction_mc]:
        set_root_object(graph, marker_style=20)

    out_dir.cd()
    graph_primaries.Write(f'g_primaries_integral_{sign}')
    graph_primaries_mc.Write(f'g_primaries_integral_mc_{sign}')
    graph_total.Write(f'g_total_integral_{sign}')
    graph_primary_fraction.Write(f'g_primary_fraction_{sign}')
    graph_primary_fraction_mc.Write(f'g_primary_fraction_mc_{sign}')


def template_fitting_routine(h2_data:TH2F, h2_mc:dict, outdir:TDirectory, particle:str):

    canvas = TCanvas(f'c', '')

    tmp_pt_min, tmp_pt_max = (1.6, 5) if particle == 'He' else (0.4, 4)
    dca_min, dca_max = -0.05, 0.05#-0.15, 0.15
    
    dca = RooRealVar('dca', f'DCA_#it{{xy}}', dca_min, dca_max, 'cm')

    gaussian_core_pars = {
        'mean': RooRealVar(f'mean_smearning', f'#mu^{{core}}', 0., -0.05, 0.05),
        'sigma': RooRealVar(f'sigma_smearning', f'#sigma^{{core}}', 1.e-3, 1.e-5, 0.1),
    }
    gaussian_core = RooGaussian('gaussian_core', 'gaussian_core', dca, *gaussian_core_pars.values())
    dca.setRange('core', -0.004, 0.006)

    for sign in ['matter', 'antimatter']:
        
        fit_results = None

        if sign == 'antimatter':
            pt_min, pt_max = -tmp_pt_max, -tmp_pt_min
            h2_mc.pop('material')
        else:
            pt_min, pt_max = tmp_pt_min, tmp_pt_max

        for pt_bin in range(h2_data.GetXaxis().FindBin(pt_min), h2_data.GetXaxis().FindBin(pt_max)):
        
            pdfs, pdf_params = {}, {}
            for flag in h2_mc.keys():
                if sign == 'antimatter' and flag == 'material':
                    continue
                pdfs[flag], pdf_params[flag] = build_model(dca, flag, mean=gaussian_core_pars['mean'])

            pt = h2_data.GetXaxis().GetBinCenter(pt_bin)
            h_dca = h2_data.ProjectionY(f'h_dca_{pt:.2f}', pt_bin, pt_bin, 'e')
            if h_dca.GetEntries() < 100:
                continue
            
            bin_outdir = outdir.mkdir(f'pt_{pt:.2f}')

            dca_frame, ifit_results = _fit_slice(h2_data, h2_mc, pdfs, pdf_params,
                                                 gaussian_core, gaussian_core_pars,
                                                 dca, pt_bin, particle, bin_outdir)
            
            if fit_results is None:
                fit_results = pd.DataFrame.from_dict([ifit_results])
            else:
                fit_results = pd.concat([fit_results, pd.DataFrame.from_dict([ifit_results])], ignore_index=True)

            #text = write_params_to_text({**signal_pars#, **bkg_pol1_pars
            #                             }.values(), coordinates=[0.4, 0.15, 0.63, 0.4])
            text = TPaveText(0.45, 0.2, 0.55, 0.25, 'ndc')
            text.SetBorderSize(0)
            text.SetFillColor(0)
            text.AddText(f'#chi^{{2}} / NDF = {dca_frame.chiSquare():.2f}')
            dca_frame.addObject(text)
            dca_frame.SetMinimum(1)

            dca_frame.Draw()
            canvas.SetLogy()

            bin_outdir.cd()
            canvas.Write(f'dca_frame_{sign}_{pt:.2f}')
            canvas.Clear()

            pdfs.clear()
            pdf_params.clear()

        draw_results(h2_data, fit_results, outdir, sign)

    del canvas

def matter_antimatter_ratio(particle:str, outdir: TDirectory):

    g_matter = outdir.Get(f'g_primaries_integral_matter')
    g_antimatter = outdir.Get(f'g_primaries_integral_antimatter')

    g_matter_mc = outdir.Get(f'g_primaries_integral_mc_matter')
    g_antimatter_mc = outdir.Get(f'g_primaries_integral_mc_antimatter')

    npoints = g_matter.GetN()
    g_ratio = TGraphErrors(npoints)
    g_ratio_mc = TGraphErrors(npoints)
    
    for ipoint in range(npoints):
        x_value = g_matter.GetPointX(ipoint)

        y_value = g_antimatter.Eval(x_value)/g_matter.Eval(x_value) if g_matter.Eval(x_value) > 0. else 0.
        y_error = y_value * np.sqrt( 1./g_antimatter.Eval(x_value) + 1./g_matter.Eval(x_value) ) if (g_matter.Eval(x_value) > 0. and g_antimatter.Eval(x_value) > 0.) else 0.
        g_ratio.SetPoint(ipoint, x_value, y_value)
        g_ratio.SetPointError(ipoint, g_matter.GetErrorX(ipoint), y_error)

        y_value_mc = g_antimatter_mc.Eval(x_value)/g_matter_mc.Eval(x_value) if g_matter_mc.Eval(x_value) > 0. else 0.
        y_error_mc = y_value_mc * np.sqrt( 1./g_antimatter_mc.Eval(x_value) + 1./g_matter_mc.Eval(x_value) ) if (g_matter_mc.Eval(x_value) > 0. and g_antimatter_mc.Eval(x_value) > 0.) else 0.
        g_ratio_mc.SetPoint(ipoint, x_value, y_value_mc)
        g_ratio_mc.SetPointError(ipoint, g_matter_mc.GetErrorX(ipoint), y_error_mc)

    set_root_object(g_ratio, name='g_ratio', title=';#it{p}_{T} (GeV/#it{c}); Antimatter / Matter', marker_style=20)
    set_root_object(g_ratio_mc, name='g_ratio_mc', title='MC;#it{p}_{T} (GeV/#it{c}); Antimatter / Matter', marker_style=20)
    
    outdir.cd()
    g_ratio.Write()
    g_ratio_mc.Write()


def main():
    
    outfile = TFile.Open('output/primary_fraction_estimation.root', 'RECREATE')
    
    for particle in ['He', 'Pr']:

        outdir = outfile.mkdir(particle)
        
        h2_mc = {}
        mc_flags = {
            'primaries': '_IsPhysicalPrimary', 'material': '_IsSecondaryFromMaterial',
            'weak_decay': '_IsSecondaryFromWeakDecay', 
            #'li4': '_IsFromLi4'
        } if particle == 'He' else {
            'primaries': '_IsPhysicalPrimary', 'material': '_IsSecondaryFromMaterial',
            'weak_decay': '_IsFromLambda0', #'sigmaplus': '_IsFromSigmaPlus', 
            #'li4': '_IsFromLi4'
        }
        
        for flag_name, flag in mc_flags.items():
            h2_dca_vs_pt_mc = load_hist('output/dca_mc_template.root', f'{particle}/h2DCAxyPt{particle}{flag}')
            h2_mc[flag_name] = h2_dca_vs_pt_mc

        h2_data = load_hist('output/dca_data_template.root', f'{particle}/h2DCAxyPt{particle}')

        template_fitting_routine(h2_data, h2_mc, outdir, particle)
        matter_antimatter_ratio(particle, outdir)

    outfile.Close()


if __name__ == '__main__':

    main()
