"""
Fitting routines and result management
"""

import numpy as np
from typing import Dict, Tuple, Optional
from ROOT import (TH2F, TDirectory, RooRealVar, RooAddPdf, RooProduct, 
                  RooDataHist, TCanvas, TLegend, TPaveText, TLine, TH2D, TF1,
                   kRed, kGreen, kBlue, kMagenta, kCyan)

from torchic.core.histogram import load_hist
from torchic.utils.root import set_root_object

from core.primary_fraction_models import build_crystal_ball_model, build_gaussian_model, build_pol0_model
from core.primary_fraction_models import build_smearing_gaussian
from core.primary_fraction_templates import prepare_convoluted_template, prepare_template
from core.primary_fraction_config import FitConfig


def get_dca_one_sigma(pt: float, particle: str, is_mc: bool = False) -> float:
    """Get 1-sigma DCA selection window"""
    # mean, 1 sigma
    if particle == 'Pr':
        return (8.18e-5, (0.0040 + 0.0026 / np.abs(pt)**(1.1714))) if not is_mc \
            else (-3.60e-5, (0.0013 + 0.0021 / np.abs(pt)**(1.4722)))
    elif particle == 'He':
        return (1.09e-4, (0.0011 + 0.0065 / np.abs(pt)**(1.0399))) if not is_mc \
            else (-4.54e-5, (0.0015 + 0.0081 / np.abs(pt)**(1.6000)))
    return 0.


def extract_covariance_matrix(fit_result) -> Optional[np.ndarray]:
    """
    Extract covariance matrix from RooFit result
    
    Args:
        fit_result: RooFitResult object
    
    Returns:
        Covariance matrix as numpy array, or None if extraction fails
    """
    if not fit_result:
        return None
    
    try:
        cov_matrix = fit_result.covarianceMatrix()
        n_params = cov_matrix.GetNrows()
        cov_array = np.zeros((n_params, n_params))
        
        for i in range(n_params):
            for j in range(n_params):
                cov_array[i, j] = cov_matrix[i][j]
        
        return cov_array
    except:
        return None


def fit_slice(h2_data: TH2F, h2_mc: dict, 
              pdfs: dict, pdf_params: dict, 
              gaussian_core, gaussian_core_pars: dict, 
              dca: RooRealVar, pt_bin: int, particle: str,
              bin_outdir: TDirectory, fit_config: FitConfig,
              h_fraction_he3: Optional[TH2F] = None) -> tuple:
    """
    Fit a single pt slice
    
    Returns:
        Tuple of (frame, fit_results_dict) or (None, None) if fit fails
    """
    pt = h2_data.GetXaxis().GetBinCenter(pt_bin)
    h_dca = h2_data.ProjectionY(f'h_dca_{pt:.2f}', pt_bin, pt_bin, 'e')
    gaus_core = TF1('gaus_core', 'gaus', -0.01, 0.01)
    h_dca.Fit(gaus_core, 'RMS+Q', '', -0.01, 0.01)
    gaussian_core_pars['mean'].setVal(gaus_core.GetParameter(1))
    
    if h_dca.GetEntries() < fit_config.min_entries_for_fit:
        return None, None

    pt_low_edge = h2_data.GetXaxis().GetBinLowEdge(pt_bin)
    pt_high_edge = h2_data.GetXaxis().GetBinLowEdge(pt_bin + 1)

    # Fit core Gaussian for smearing
    dh = RooDataHist(f'dh_{pt:.2f}', 'dh', [dca], Import=h_dca)
    #gaussian_core.fitTo(dh, PrintLevel=1, Save=True, Range='core')
    
    # Build smearing Gaussian
    gaussian_smearing, gaussian_smearing_pars = build_smearing_gaussian(
        dca, gaussian_core_pars['mean'].getVal()
    )

    if np.abs(pt) > fit_config.max_pt_material_template and 'material' in pdfs:
        # Remove material template for high pT if configured
        del pdfs['material']
        del pdf_params['material']

    stored_pdfs_material = []
    if fit_config.add_pol0_to_mc_template_for_material and 'material' in pdfs:
        pol0_material, pol0_material_pars = build_pol0_model(dca, '_material')
        frac_pol0_material = RooRealVar('frac_pol0_material', 'Fraction of pol0 in material', 0.9, 0., 1.)
        combined_material_pdf = RooAddPdf('combined_material', 'combined_material', 
                                         [pdfs['material'], pol0_material], 
                                         [frac_pol0_material])
        
        stored_pdfs_material.append(pdfs['material'])
        stored_pdfs_material.append(pol0_material)
        
        pdfs['material'] = combined_material_pdf

        pdf_params['material']['frac_pol0_material'] = frac_pol0_material
        for par_name, par in pol0_material_pars.items():
            pdf_params['material'][f'pol0_{par_name}'] = par

    stored_pdfs_primaries = []
    if fit_config.add_pol0_to_mc_template_for_primaries and 'primaries' in pdfs:
        pol0_primaries, pol0_primaries_pars = build_pol0_model(dca, '_primaries')
        pol0_primaries_pars['c0'].setVal(1.)
        pol0_primaries_pars['c0'].setConstant(True)
        frac_pol0_primaries = RooRealVar('frac_pol0_primaries', 'Fraction of pol0 in primaries', 0.999, 0.99, 1.)
        combined_primaries_pdf = RooAddPdf('combined_primaries', 'combined_primaries', 
                                         [pdfs['primaries'], pol0_primaries], 
                                         [frac_pol0_primaries])
        
        stored_pdfs_primaries.append(pdfs['primaries'])
        stored_pdfs_primaries.append(pol0_primaries)
        
        pdfs['primaries'] = combined_primaries_pdf

        pdf_params['primaries']['frac_pol0_primaries'] = frac_pol0_primaries
        for par_name, par in pol0_primaries_pars.items():
            pdf_params['primaries'][f'pol0_{par_name}'] = par

    # Prepare templates
    convoluted_pdfs, normalisations = {}, {}
    yield_mc = {}
    total_yield_mc = 0
    fix_material_yield = fit_config.fix_material_yield
    very_low_count_rejection = False

    for flag in h2_mc.keys():

        if flag not in pdfs.keys():
            continue

        pt_bin_flag = h2_mc[flag].GetXaxis().FindBin(pt)
        h_dca_flag = h2_mc[flag].ProjectionY('h_mc', pt_bin_flag, pt_bin_flag, 'e')
        bin_min = h_dca_flag.FindBin(dca.getMin())
        bin_max = h_dca_flag.FindBin(dca.getMax())
        
        if h_dca_flag.Integral(bin_min, bin_max) < fit_config.min_entries_for_template:
            very_low_count_rejection = True
            continue

        if fit_config.use_convolutional_smearing:
            
            convoluted_pdfs[flag] = prepare_convoluted_template(
                dca, pdfs[flag], pdf_params[flag],
                h_dca_flag, gaussian_smearing, gaussian_smearing_pars,
                bin_outdir, flag, particle
            )
        else:
            convoluted_pdfs[flag] = prepare_template(
                dca, pdfs[flag], pdf_params[flag],
                h_dca_flag, bin_outdir, flag, particle
            )
        
        pdf_params[flag]['mean'].setVal(gaussian_core_pars['mean'].getVal())
        
        yield_mc[flag] = h_dca_flag.Integral(bin_min, bin_max)
        total_yield_mc += yield_mc[flag]
        del h_dca_flag

    # Handle missing material template
    if 'material' not in convoluted_pdfs.keys() and not very_low_count_rejection and 0 < pt < fit_config.max_pt_material_template:
        if fit_config.material_background_type == 'gaussian':
            pdf_material, params_material = build_gaussian_model(dca, 'material')
        elif fit_config.material_background_type == 'pol0':
            pdf_material, params_material = build_pol0_model(dca, 'material')
        else:  # gaussian fallback
            pdf_material, params_material = build_gaussian_model(dca, 'material')
        
        pdf_params['material'] = params_material
        convoluted_pdfs['material'] = pdf_material
        yield_mc['material'] = yield_mc.get('primaries', 1.0) * 0.05

    # Build normalisations
    fraction_material_to_primary = None
    fraction_weak_decay_to_primary = None

    for flag in yield_mc.keys():
        init_val = yield_mc[flag] / total_yield_mc * h_dca.GetEntries() if total_yield_mc > 0 else 0
        normalisations[flag] = RooRealVar(f'frac_{flag}', f'#it{{f}}_{{{flag}}}', init_val, 0., h_dca.GetEntries())

        if flag == 'weak_decay':
            if particle == 'Pr':
                fraction_weak_decay_to_primary = RooRealVar(
                    'fraction_weak_decay_to_primary', 'fraction_weak_decay_to_primary',
                    yield_mc['weak_decay'] / yield_mc['primaries'] if yield_mc['primaries'] > 0 else 0.1,
                    0., 1.
                )
                if fit_config.fix_weak_decay_yield:
                    fraction_weak_decay_to_primary.setConstant(True)
                    
                normalisations['weak_decay'] = RooProduct(
                    'frac_weak_decay', '#it{f}_{weak_decay}',
                    [normalisations['primaries'], fraction_weak_decay_to_primary]
                )
            elif particle == 'He' and h_fraction_he3 is not None and not fit_config.initialise_from_mc_yields:
                value = h_fraction_he3.GetBinContent(h_fraction_he3.FindBin(pt))
                error = h_fraction_he3.GetBinError(h_fraction_he3.FindBin(pt))
                fraction_weak_decay_to_primary = RooRealVar(
                    'fraction_weak_decay_to_primary', 'fraction_weak_decay_to_primary',
                    value, value - 2 * error, value + 2 * error
                )
                normalisations['weak_decay'] = RooProduct(
                    'frac_weak_decay', '#it{f}_{weak_decay}',
                    [normalisations['primaries'], fraction_weak_decay_to_primary]
                )
            elif particle == 'He' and fit_config.initialise_from_mc_yields:
                fraction_value = yield_mc['weak_decay'] / yield_mc['primaries'] if yield_mc['primaries'] > 0 else 0.1
                fraction_weak_decay_to_primary = RooRealVar(
                    'fraction_weak_decay_to_primary', 'fraction_weak_decay_to_primary',
                    fraction_value, 0.9 * fraction_value, 1.1 * fraction_value
                )
                if fit_config.fix_weak_decay_yield:
                    fraction_weak_decay_to_primary.setConstant(True)
                    
                normalisations['weak_decay'] = RooProduct(
                    'frac_weak_decay', '#it{f}_{weak_decay}',
                    [normalisations['primaries'], fraction_weak_decay_to_primary]
                )

        if flag == 'material':
            if fit_config.fix_material_yield:
                normalisations[flag].setConstant(True)
            
            if particle == 'Pr':
                fraction_material_to_primary = RooRealVar(
                    'fraction_material_to_primary', 'fraction_material_to_primary',
                    yield_mc['material'] / yield_mc['primaries'] if yield_mc['primaries'] > 0 else 0.1
                )
                if fit_config.fix_material_yield:
                    fraction_material_to_primary.setConstant(True)
                    
                normalisations['material'] = RooProduct(
                    'frac_material', '#it{f}_{material}',
                    [normalisations['primaries'], fraction_material_to_primary]
                )
            elif particle == 'He':
                fraction_value = yield_mc['material'] / yield_mc['primaries'] if yield_mc['primaries'] > 0 else 0.1
                fraction_material_to_primary = RooRealVar(
                    'fraction_material_to_primary', 'fraction_material_to_primary',
                    #fraction_value, 0.9 * fraction_value, 1.1 * fraction_value
                    fraction_value, 0.1 * fraction_value, 5. * fraction_value
                )
                if fit_config.fix_material_yield:
                    fraction_material_to_primary.setConstant(True)
                    
                normalisations['material'] = RooProduct(
                    'frac_material', '#it{f}_{material}',
                    [normalisations['primaries'], fraction_material_to_primary]
                )
        

    # Perform fit
    model = RooAddPdf('total_fit', 'total fit', list(convoluted_pdfs.values()), list(normalisations.values()))
    fit_result = model.fitTo(dh, PrintLevel=1, Save=True)

    # Extract covariance if requested
    covariance = None
    if fit_config.store_covariance:
        covariance = extract_covariance_matrix(fit_result)

    # Create frame for visualization
    frame = dca.frame(Title=f'{pt_low_edge:.2f} < #it{{p}}_{{T}} < {pt_high_edge:.2f} GeV/#it{{c}}')

    legend = TLegend(0.15, 0.74, 0.35, 0.88)
    legend.SetBorderSize(0)
    
    colors  = [kRed, kGreen+2, kBlue+2, kMagenta+2, kCyan+2]

    dh.plotOn(frame)
    model.plotOn(frame, Name=model.GetName(), LineColor=colors[0], Precision=1e-6)
    #model.paramOn(frame, Layout=(0.6, 0.88, 0.88))
    legend.AddEntry(frame.findObject(model.GetName()), model.GetTitle(), "l")

    params_text = TPaveText(0.6, 0.74, 0.88, 0.88, 'ndc')
    params_text.SetBorderSize(0)
    params_text.SetFillColor(0)
    params_text.SetTextSize(0.035)


    for iflag, flag in enumerate(convoluted_pdfs.keys()):
        model.plotOn(frame, Name=convoluted_pdfs[flag].GetName(), Components={convoluted_pdfs[flag]}, Precision=1e-6,
                    LineColor=colors[1+iflag], LineStyle='--')
        legend.AddEntry(frame.findObject(convoluted_pdfs[flag].GetName()), convoluted_pdfs[flag].GetTitle(), "l")
        if flag == 'primaries':
            params_text.AddText(f'#it{{N}}_{{{flag}}} = {normalisations[flag].getVal():.0f} #pm {normalisations[flag].getError():.0f}') 
        else:
            params_text.AddText(f'#it{{N}}_{{{flag}}} = {normalisations[flag].getVal():.0f}') 
    
    frame.addObject(legend)
    frame.addObject(params_text)

    # Calculate integrals
    dca.setRange('full_range', dca.getMin(), dca.getMax())
    mean, sigma = get_dca_one_sigma(pt, particle, is_mc=fit_config.is_mc)
    selection_window = (mean -3 * sigma, mean + 3 * sigma)
    dca.setRange('integral_range', selection_window[0], selection_window[1])
    
    primaries_integral = (convoluted_pdfs['primaries'].createIntegral(dca, dca, 'integral_range').getVal() 
                         * normalisations['primaries'].getVal())
    primaries_integral_full = (convoluted_pdfs['primaries'].createIntegral(dca, dca, 'full_range').getVal() 
                         * normalisations['primaries'].getVal())
    primaries_integral_error = (convoluted_pdfs['primaries'].createIntegral(dca, dca, 'integral_range').getVal() 
                             * normalisations['primaries'].getError())
    
    total_normalisation = sum(norm.getVal() for norm in normalisations.values())
    total_integral = model.createIntegral(dca, dca, 'integral_range').getVal() * total_normalisation

    # Calculate fractions
    material_integral, materrial_integral_error = 0., 0.
    weak_decay_integral, weak_decay_integral_error = 0., 0.

    
    if 'material' in convoluted_pdfs:
        material_integral = (convoluted_pdfs['material'].createIntegral(dca, dca, 'integral_range').getVal()
                           * normalisations['material'].getVal())
        materrial_integral_error = (convoluted_pdfs['material'].createIntegral(dca, dca, 'integral_range').getVal()
                             * normalisations['primaries'].getVal() * fraction_material_to_primary.getError() if fraction_material_to_primary else 0.)
    
    if 'weak_decay' in convoluted_pdfs:
        weak_decay_integral = (convoluted_pdfs['weak_decay'].createIntegral(dca, dca, 'integral_range').getVal()
                             * normalisations['weak_decay'].getVal())
        weak_decay_integral_error = (convoluted_pdfs['weak_decay'].createIntegral(dca, dca, 'integral_range').getVal()
                                * normalisations['primaries'].getVal() * fraction_weak_decay_to_primary.getError() if fraction_weak_decay_to_primary else 0.)

    # Add selection lines
    lines = [
        TLine(selection_window[0], 1., selection_window[0], frame.GetMaximum()),
        TLine(selection_window[1], 1., selection_window[1], frame.GetMaximum())
    ]
    for line in lines:
        set_root_object(line, line_style=2, line_color=4)
        frame.addObject(line)

    # Add chi2 text
    text = TPaveText(0.6, 0.68, 0.88, 0.72, 'ndc')
    text.SetBorderSize(0)
    text.SetFillColor(0)
    text.AddText(f'#chi^{{2}} / NDF = {frame.chiSquare():.2f}')
    frame.addObject(text)
    frame.SetMinimum(1)

    # Prepare results
    fit_results = {
        'pt': np.abs(pt),
        'primaries_integral': primaries_integral,
        'primaries_integral_error': primaries_integral_error,
        'primaries_fraction_in_range': primaries_integral / primaries_integral_full,
        'primaries_integral_mc': yield_mc.get('primaries', 0.),
        'material_integral': material_integral,
        'material_integral_error': materrial_integral_error,
        'material_integral_mc': yield_mc.get('material', 0.),
        'weak_decay_integral': weak_decay_integral,
        'weak_decay_integral_error': weak_decay_integral_error,
        'weak_decay_integral_mc': yield_mc.get('weak_decay', 0.),
        'total_integral': total_integral,
        'primary_fraction': primaries_integral / total_integral if total_integral > 0 else 0.,
        'primary_fraction_error': (primaries_integral_error / total_integral
                                    if total_integral > 0. else 0.),
        'primary_fraction_mc': (yield_mc.get('primaries', 0.) / total_yield_mc 
                               if total_yield_mc > 0. else 0.),
        'material_fraction': material_integral / total_integral if total_integral > 0 else 0.,
        'material_fraction_error': (materrial_integral_error / total_integral 
                                    if total_integral > 0. else 0.),
        'material_fraction_mc': (yield_mc.get('material', 0.) / total_yield_mc 
                                if total_yield_mc > 0. else 0.),
        'weak_decay_fraction': weak_decay_integral / total_integral if total_integral > 0 else 0.,
        'weak_decay_fraction_error': (weak_decay_integral_error / total_integral 
                                     if total_integral > 0. else 0.),
        'weak_decay_fraction_mc': (yield_mc.get('weak_decay', 0.) / total_yield_mc 
                                  if total_yield_mc > 0. else 0.),
        'covariance': covariance,
        'mean': gaus_core.GetParameter(1),
        'mean_error': gaus_core.GetParError(1),
        'sigma': gaus_core.GetParameter(2),
        'sigma_error': gaus_core.GetParError(2)
    }

    del dh, model, fit_result, gaus_core, lines, text, params_text, legend, h_dca

    return frame, fit_results