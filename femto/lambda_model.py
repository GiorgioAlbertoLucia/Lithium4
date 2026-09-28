'''
    Create a model that includes the presence of lambda parameters.
'''

import numpy as np

from ROOT import (TFile, TDirectory, TCanvas, TLegend, TTree, TF1, \
                  RooDataSet, RooKeysPdf, RooRealVar, RooFit, RooDataHist,
                  RooMsgService, RooFit)

from torchic.core.histogram import load_hist, HistLoadInfo
from torchic.utils.root import set_root_object, init_legend
from torchic.utils.colors import get_color
from torchic.roopdf.roopdf_utils import init_roopdf
import argparse
from core.config_loader import load_yaml, build_hist_load_info_dict, build_hist_load_info_variations_radius, build_hist_load_info_variations_model_params

LAMBDA_MODIFICATION_FACTOR = None
LAMBDA_VARIATION = None

INPUT_CK_PATH = None
INPUT_SIGMA_CK_PATH = None
INPUT_CK_VARIATIONS = None
INPUT_CK_MODEL_PARAM_VARIATIONS = None
INPUT_SIGMA_CK_VARIATIONS = None

CENTRALITY_BINS = None

INPUT_LAMBDA_PARAMETER_PATH = None
INPUT_RESOLUTION_PATH = None
INPUT_MIXED_EVENT_REFERENCE_PATH = None
OUTPUT_LAMBDA_MODEL_PATH = None


def match_bin_width_correlation_function(h_source, h_target, kstar_threshold:float=0.4):
    """
    Adjust the bin widths of the target histogram to match those of the source histogram.
    
    Args:
        h_source: Histogram with desired bin widths
        h_target: Histogram to be adjusted
    """
    
    h_target_matched = h_source.Clone(f'{h_target.GetName()}_matched')
    for ibin in range(1, h_target_matched.GetNbinsX()+1):
        

        if h_target_matched.GetBinCenter(ibin) > kstar_threshold: 
            continue

        kstar_low = h_target_matched.GetBinLowEdge(ibin)
        kstar_high = h_target_matched.GetBinLowEdge(ibin+1)
        bin_sum, bin_count = 0., 0.

        for jbin in range(1, h_target.GetNbinsX()+1):
            bin_center = h_target.GetBinCenter(jbin)
            if kstar_low <= bin_center < kstar_high:
                bin_sum += h_target.GetBinContent(jbin)
                bin_count += 1
        
        if ibin == h_target_matched.GetNbinsX():
            print(f'Last bin: kstar = {h_target_matched.GetBinCenter(ibin):.3f} GeV/c, content = {h_target_matched.GetBinContent(ibin):.3f}')
            print(f'CK: kstar = {bin_sum} GeV/c, content = {bin_count}')
        
        if bin_count > 0:
            h_target_matched.SetBinContent(ibin, bin_sum / bin_count)
        else:
            h_target_matched.SetBinContent(ibin, 0.)
    
    return h_target_matched

def smoothen_histogram(hist, outfile, n_events:int=100_000, xmin:float=0.02, xmax:float=0.42, rho:float=0.05):

    data, data_weights = np.array([]), np.array([])
    kstar = RooRealVar('kstar', 'k*', xmin, xmax)
    
    # Sample x-values uniformly, weight by the function value
    x_data = []
    weights = []
    
    for ibin in range(1, hist.GetNbinsX()+1):
        kstar_val = hist.GetBinCenter(ibin)
        if not (xmin <= kstar_val <= xmax):
            continue
        y_val = hist.GetBinContent(ibin)
        if y_val <= 0:
            continue
        
        x_data.append(kstar_val)
        weights.append(y_val)
    
    # Create weighted RooDataSet
    tree = TTree('tree', 'tree')
    x = np.zeros(1, dtype=np.float64)
    w = np.zeros(1, dtype=np.float64)
    tree.Branch('kstar', x, 'kstar/D')
    tree.Branch('weight', w, 'weight/D')
    
    for x_val, w_val in zip(x_data, weights):
        x[0] = x_val
        w[0] = w_val
        tree.Fill()
    
    # Create weighted dataset
    weight_var = RooRealVar('weight', 'weight', 0, 1e6)
    dataset = RooDataSet(hist.GetName()+'_roodata', hist.GetName()+'_roodata', 
                         [kstar, weight_var], 
                         RooFit.Import(tree), 
                         RooFit.WeightVar('weight'))
    
    # Create RooKeysPdf with the weighted data
    keys_pdf = RooKeysPdf(f'{hist.GetName()}_keys', f'{hist.GetName()}_keys', 
                          kstar, dataset, RooKeysPdf.NoMirror, rho)
    

    frame = kstar.frame()
    dataset.plotOn(frame, MarkerStyle=20, MarkerSize=0.8, LineColor=1)
    keys_pdf.plotOn(frame, LineColor=2, LineWidth=2)
    canvas = TCanvas(f'cKeysPdf_{hist.GetName()}', f'cKeysPdf_{hist.GetName()}', 800, 600)
    frame.Draw()

    
    outfile.cd()
    hist.Write(f'{hist.GetName()}_original')
    keys_pdf.Write(f'{hist.GetName()}_keyspdf')
    canvas.Write()
    
    return keys_pdf

def compute_lambda_corrected_values(kstars, ck_values, lambda_param_hist, scale: float = 1.0,
                                     sigma_ck_values=None, lambda_sigma_hist=None):
    '''
    Compute the lambda-corrected correlation value at each kstar point.
    If sigma_ck_values/lambda_sigma_hist are provided, includes the Sigma contamination term;
    otherwise applies the plain lambda*ck + (1-lambda) correction.
    '''
    corrected = []
    for i, (kstar, ck) in enumerate(zip(kstars, ck_values)):
        lam = lambda_param_hist.GetBinContent(lambda_param_hist.FindBin(kstar)) * scale
        lam = min(max(lam, 0.), 1.)  # Ensure lambda is between 0 and 1
        if (lam-0) < 1e-12:
            lam = 1. # lambda = 0 happens for no entries in the kstar histogram

        if lambda_sigma_hist is not None:
            lam_s = lambda_sigma_hist.GetBinContent(lambda_sigma_hist.FindBin(kstar))
            sig = sigma_ck_values[i]
        else:
            lam_s, sig = 0., 0.

        corrected.append(lam * ck + lam_s * sig + (1 - lam - lam_s))
    return corrected

def fill_hist_from_weighted_average(h_out, kstars, values):
    '''
    Fill h_out bin-by-bin by averaging all (kstar, value) pairs whose kstar falls
    within that bin's edges. Works regardless of whether `kstars` shares h_out's binning.
    '''
    for ibin in range(1, h_out.GetNbinsX()+1):
        kstar_low = h_out.GetBinLowEdge(ibin)
        kstar_high = h_out.GetBinLowEdge(ibin+1)

        val_sum, count = 0., 0
        for kstar, val in zip(kstars, values):
            if kstar_low <= kstar < kstar_high:
                val_sum += val
                count += 1

        if count > 0:
            h_out.SetBinContent(ibin, val_sum / count)
    return h_out

def apply_lambda_correction(h_out, h_ck, h_sigma_ck, lambda_param_hist, lambda_sigma_hist, scale: float = 1.0):
    kstars = [h_ck.GetBinCenter(jbin) for jbin in range(1, h_ck.GetNbinsX()+1)]
    ck_values = [h_ck.GetBinContent(jbin) for jbin in range(1, h_ck.GetNbinsX()+1)]
    sigma_values = [h_sigma_ck.GetBinContent(h_sigma_ck.FindBin(kstar)) for kstar in kstars]

    corrected = compute_lambda_corrected_values(kstars, ck_values, lambda_param_hist, scale=scale,
                                                 sigma_ck_values=sigma_values, lambda_sigma_hist=lambda_sigma_hist)
    fill_hist_from_weighted_average(h_out, kstars, corrected)
        
def produce_lambda_with_modified_values(sign:str, centrality:str, h_theoretical_Ck, h_theoretical_Sigma_Ck,
                                        outdir:TDirectory, modification_factor:float):
    '''
        Produce a lambda-corrected model by systematically changing the lambda value by x%, 
        where x is the modification_factor. This can be used to understand the sensitivity of the model to changes in lambda.
        The lambda is changed both in lambda + x% and lambda - x% to understand the effect in both directions.
    '''
    
    if not INPUT_LAMBDA_PARAMETER_PATH:
        raise ValueError("INPUT_LAMBDA_PARAMETER_PATH is not set. Please set it before calling this function.")

    centrality_dir = 'centrality_0_10' if centrality == '010' else 'centrality_10_50'
    h_lambda_parameter = load_hist(INPUT_LAMBDA_PARAMETER_PATH, f'{centrality_dir}/{sign}/hLambdaParameters')
    h_lambda_Sigma_parameter = load_hist(INPUT_LAMBDA_PARAMETER_PATH, f'{centrality_dir}/{sign}/hLambdaSigmaParameters')
    
    h_lambda_parameter_higher_lambda = h_lambda_parameter.Clone(f'hLambdaParameter_HigherLambda')
    h_lambda_parameter_lower_lambda = h_lambda_parameter.Clone(f'hLambdaParameter_LowerLambda')

    h_lambda_corrected_Ck_higher_lambda = h_theoretical_Ck.Clone(f'hLambdaCorrectedCk_Higher')
    apply_lambda_correction(h_lambda_corrected_Ck_higher_lambda, h_theoretical_Ck, h_theoretical_Sigma_Ck,
                            h_lambda_parameter, h_lambda_Sigma_parameter, scale=1 + modification_factor)
    h_lambda_corrected_Ck_lower_lambda = h_theoretical_Ck.Clone(f'hLambdaCorrectedCk_Lower')
    apply_lambda_correction(h_lambda_corrected_Ck_lower_lambda, h_theoretical_Ck, h_theoretical_Sigma_Ck,
                            h_lambda_parameter, h_lambda_Sigma_parameter, scale=1 - modification_factor)

    h_lambda_Sigma_corrected_Ck_higher_lambda = h_theoretical_Ck.Clone(f'hLambdaSigmaCorrectedCk_Higher')
    apply_lambda_correction(h_lambda_Sigma_corrected_Ck_higher_lambda, h_theoretical_Ck, h_theoretical_Sigma_Ck,
                            h_lambda_parameter, h_lambda_Sigma_parameter, scale=1 + modification_factor)
    h_lambda_Sigma_corrected_Ck_lower_lambda = h_theoretical_Ck.Clone(f'hLambdaSigmaCorrectedCk_Lower')
    apply_lambda_correction(h_lambda_Sigma_corrected_Ck_lower_lambda, h_theoretical_Ck, h_theoretical_Sigma_Ck,
                            h_lambda_parameter, h_lambda_Sigma_parameter, scale=1 - modification_factor)

    return (h_lambda_parameter_higher_lambda, h_lambda_parameter_lower_lambda,
            h_lambda_corrected_Ck_lower_lambda, h_lambda_corrected_Ck_higher_lambda, 
            h_lambda_Sigma_corrected_Ck_lower_lambda, h_lambda_Sigma_corrected_Ck_higher_lambda)

def precompute_resolution_fits(outfile: TFile) -> dict:
    """
    Fit each kstar slice of the resolution matrix with a RooFit Crystal Ball PDF.
    Returns a dict mapping kstar bin index (1-based) -> callable(x) using the fitted PDF.
    Saves all fits to outfile under 'ResolutionFits/'.
    """
    
    if not INPUT_RESOLUTION_PATH:
        raise ValueError("INPUT_RESOLUTION_PATH is not set. Please set it before calling this function.")
    
    h_resolution = load_hist(INPUT_RESOLUTION_PATH, 'he3-hadron-femto/QA/hKstarRecVsKstarGen')

    outdir_fits = outfile.mkdir('ResolutionFits')
    fits = {}

    x = RooRealVar(f'x', 'x', 0., 0.7)
    n_bins_x = h_resolution.GetNbinsX()
    for ibin in range(1, n_bins_x + 1):
        h_slice = h_resolution.ProjectionY(f'hResSlice_bin{ibin}', ibin, ibin)

        if h_slice.GetEntries() < 10:
            continue

        peak  = h_slice.GetBinCenter(h_slice.GetMaximumBin())
        sigma_est = max(h_slice.GetRMS(), 1e-4)
        
        crystal_ball, pars = init_roopdf('crystal_ball', x, 
                                   mean=RooRealVar(f'mean_{ibin}', 'mean', peak),
                                   sigma=RooRealVar(f'sigma_{ibin}', 'sigma', sigma_est, 0.00001, 0.1),
                                   aL=RooRealVar(f'alpha_{ibin}', 'alpha', 1.5, 0.5, 5.0),
                                   nL=RooRealVar(f'n_{ibin}', 'n', 25.0, 20., 100.0),
                                   aR=RooRealVar(f'alphaR_{ibin}', 'alphaR', 1.5, 0.5, 5.0),
                                   nR=RooRealVar(f'nR_{ibin}', 'nR', 25.0, 20., 100.0),)    

        dataset = RooDataHist(f'ds_{ibin}', f'ds_{ibin}', x, Import=h_slice)
        crystal_ball.fitTo(dataset, RooFit.PrintLevel(-1))
        for par in pars.values():
            par.setConstant(True)

        fits[ibin] = (x, crystal_ball, pars)

        frame = x.frame()
        frame.SetTitle(f'#it{{kstar}} = {h_resolution.GetXaxis().GetBinCenter(ibin):.3f} GeV/#it{{c}}')
        dataset.plotOn(frame, MarkerStyle=20, MarkerSize=0.8, LineColor=1)
        crystal_ball.plotOn(frame, LineColor=2, LineWidth=2)
        crystal_ball.paramOn(frame, Layout=(0.55, 0.9, 0.9))
        canvas = TCanvas(f'cCrystalBallFit_bin{ibin}', f'cCrystalBallFit_bin{ibin}', 800, 600)
        frame.Draw()
        
        outdir_fits.cd()
        canvas.Write(f'cCrystalBallFit_bin{ibin}')

    return fits

def apply_resolution_smearing(h_correlation_function, outdir:TDirectory, resolution_fits:dict):

    if not INPUT_RESOLUTION_PATH:
        raise ValueError("INPUT_RESOLUTION_PATH is not set. Please set it before calling this function.")
    if not INPUT_MIXED_EVENT_REFERENCE_PATH:
        raise ValueError("INPUT_MIXED_EVENT_REFERENCE_PATH is not set. Please set it before calling this function.")

    h_resolution = load_hist(INPUT_RESOLUTION_PATH, 'he3-hadron-femto/QA/hKstarRecVsKstarGen')
    h_mixed_event = load_hist(INPUT_MIXED_EVENT_REFERENCE_PATH, 'QA/hKstar')
    mixed_fit = TF1('mixed_fit', 'pol3', 0.01, 0.4)
    h_mixed_event.Fit(mixed_fit, 'RMS+')

    # no binning-matching needed anymore: the resolution slices are described by
    # analytic Crystal Ball fits, so we can sample h_correlation_function directly
    # at any kstar_gen, independent of its own binning.
    h_smeared_correlation_function = h_correlation_function.Clone(f'{h_correlation_function.GetName()}_smeared')

    for ibin in range(1, h_smeared_correlation_function.GetNbinsX()+1):

        smeared_value, weight, total_weight = 0., 0., 0.
        kstar = h_smeared_correlation_function.GetBinCenter(ibin)
        resolution_bin = h_resolution.GetXaxis().FindBin(kstar)
        h_resolution_slice = h_resolution.ProjectionX(f'hResolutionSlice_kstar_{kstar:.3f}', resolution_bin, resolution_bin)

        for jbin in range(1, h_resolution_slice.GetNbinsX()+1):

            kstar_gen = h_resolution_slice.GetBinCenter(jbin)
            mixed_weight = mixed_fit.Eval(kstar_gen)

            fit_entry = resolution_fits.get(resolution_bin)
            slice_val = 0.
            if fit_entry is not None:
                x_var, crystal_ball_pdf, __ = fit_entry
                x_var.setVal(kstar_gen)
                slice_val = crystal_ball_pdf.getVal()
            else:
                slice_val = h_resolution_slice.GetBinContent(jbin)
            weight = slice_val * mixed_weight if 0.01 < kstar_gen < 0.7 else 0.

            # evaluate the correlation function directly at kstar_gen, regardless of its binning
            correlation_value = h_correlation_function.GetBinContent(h_correlation_function.FindBin(kstar_gen))

            smeared_value += correlation_value * weight
            total_weight += weight

        if total_weight > 0:
            smeared_value /= total_weight
            h_smeared_correlation_function.SetBinContent(ibin, smeared_value)

    return h_smeared_correlation_function

def produce_lambda_models(sign:str, centrality:str, outdir:TDirectory, resolution_fits:dict):
    
    if not INPUT_LAMBDA_PARAMETER_PATH:
            raise ValueError("INPUT_LAMBDA_PARAMETER_PATH is not set. Please set it before calling this function.")
    if not LAMBDA_MODIFICATION_FACTOR:
            raise ValueError("LAMBDA_MODIFICATION_FACTOR is not set. Please set it before calling this function.")
    if not LAMBDA_VARIATION:
            raise ValueError("LAMBDA_VARIATION is not set. Please set it before calling this function.")
        
    # only 010 and 1050 are computed 
    centrality_dir = 'centrality_0_10' if centrality == '010' else 'centrality_10_50'
    h_lambda_parameter = load_hist(INPUT_LAMBDA_PARAMETER_PATH, f'{centrality_dir}/{sign}/hLambdaParameters')
    h_lambda_Sigma_parameter = load_hist(INPUT_LAMBDA_PARAMETER_PATH, f'{centrality_dir}/{sign}/hLambdaSigmaParameters')
    h_theoretical_Ck = load_hist(INPUT_CK_PATH[centrality])
    h_theoretical_Sigma_Ck = load_hist(INPUT_SIGMA_CK_PATH[centrality])
    
    h_lambda_corrected_Ck = h_theoretical_Ck.Clone(f'hLambdaCorrectedCk')
    h_lambda_Sigma_corrected_Ck = h_theoretical_Ck.Clone(f'hLambdaSigmaCorrectedCk')

    kstars = [h_theoretical_Ck.GetBinCenter(jbin) for jbin in range(1, h_theoretical_Ck.GetNbinsX()+1)]
    ck_values = [h_theoretical_Ck.GetBinContent(jbin) for jbin in range(1, h_theoretical_Ck.GetNbinsX()+1)]
    sigma_values = [h_theoretical_Sigma_Ck.GetBinContent(h_theoretical_Sigma_Ck.FindBin(kstar)) for kstar in kstars]

    theoretical_corrected = compute_lambda_corrected_values(kstars, ck_values, h_lambda_parameter)
    theoretical_Sigma_corrected = compute_lambda_corrected_values(kstars, ck_values, h_lambda_parameter,
                                                                    sigma_ck_values=sigma_values, lambda_sigma_hist=h_lambda_Sigma_parameter)

    fill_hist_from_weighted_average(h_lambda_corrected_Ck, kstars, theoretical_corrected)
    fill_hist_from_weighted_average(h_lambda_Sigma_corrected_Ck, kstars, theoretical_Sigma_corrected)

    modification_factor = LAMBDA_VARIATION
    (h_lambda_parameter_higher_lambda, h_lambda_parameter_lower_lambda,
    h_lambda_corrected_Ck_lower_lambda, h_lambda_corrected_Ck_higher_lambda,
    h_lambda_Sigma_corrected_Ck_higher_lambda, h_lambda_Sigma_corrected_Ck_lower_lambda) \
        = produce_lambda_with_modified_values(sign, centrality, outdir=outdir,
                                            h_theoretical_Ck=h_theoretical_Ck, h_theoretical_Sigma_Ck=h_theoretical_Sigma_Ck,
                                             modification_factor=LAMBDA_MODIFICATION_FACTOR)
    
    #h_correlation_reference = load_hist(EXPERIMENTAL_CK_PATH, EXPERIMENTAL_CK_NAME)
    #h_lambda_corrected_Ck_matched = match_bin_width_correlation_function(h_correlation_reference, h_lambda_corrected_Ck)
    #h_lambda_Sigma_corrected_Ck_matched = match_bin_width_correlation_function(h_correlation_reference, h_lambda_Sigma_corrected_Ck)
    
    h_lambda_Sigma_smeared_Ck = apply_resolution_smearing(h_lambda_Sigma_corrected_Ck, outdir, resolution_fits)
    
    smoothen_histogram(h_lambda_corrected_Ck, outdir)
    smoothen_histogram(h_lambda_Sigma_corrected_Ck, outdir)
    smoothen_histogram(h_lambda_corrected_Ck_higher_lambda, outdir)
    
    outdir.cd()

    for ihist, hist in enumerate([h_lambda_parameter, h_theoretical_Ck, 
                            h_lambda_corrected_Ck, h_lambda_corrected_Ck_higher_lambda, h_lambda_corrected_Ck_lower_lambda,
                             h_lambda_Sigma_corrected_Ck, h_lambda_Sigma_corrected_Ck_higher_lambda, h_lambda_Sigma_corrected_Ck_lower_lambda,
                             h_lambda_Sigma_smeared_Ck]):
        set_root_object(hist, title='; #it{k}* (MeV/c); C(#it{k}*)', line_width=2, line_color=get_color(ihist)) 
        
    h_lambda_parameter.Write('hLambdaParameters')
    h_theoretical_Ck.Write('hTheoreticalCk')
    h_lambda_corrected_Ck.Write()
    h_lambda_Sigma_corrected_Ck.Write('hLambdaSigmaCorrectedCk')

    #h_lambda_corrected_Ck_matched.Write('hLambdaCorrectedCk_Matched')
    #h_lambda_Sigma_corrected_Ck_matched.Write('hLambdaSigmaCorrectedCk_Matched')

    h_lambda_parameter_higher_lambda.Write('hLambdaParameter_HigherLambda')
    h_lambda_parameter_lower_lambda.Write('hLambdaParameter_LowerLambda')
    h_lambda_corrected_Ck_higher_lambda.Write('hLambdaCorrectedCk_HigherLambda')
    h_lambda_corrected_Ck_lower_lambda.Write('hLambdaCorrectedCk_LowerLambda')
    h_lambda_Sigma_corrected_Ck_higher_lambda.Write('hLambdaSigmaCorrectedCk_HigherLambda')
    h_lambda_Sigma_corrected_Ck_lower_lambda.Write('hLambdaSigmaCorrectedCk_LowerLambda')

    h_lambda_Sigma_smeared_Ck.Write('hLambdaSigmaCorrectedCk_Smeared')

    canvas = TCanvas(f'cLambdaModel_{sign}', f'cLambdaModel_{sign}', 800, 600)
    hframe = canvas.DrawFrame(0.01, 0., 0.4, 1.08, '; #it{k}* (GeV/c); C(#it{k}*)')
    h_theoretical_Ck.Draw('HIST SAME')
    h_lambda_corrected_Ck.Draw('HIST SAME')
    h_lambda_Sigma_corrected_Ck.Draw('HIST SAME')
    h_lambda_Sigma_smeared_Ck.Draw('HIST SAME')

    legend = TLegend(0.6, 0.2, 0.88, 0.4)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)

    legend.AddEntry(h_theoretical_Ck, 'C_{genuine}(k*)', 'l')
    legend.AddEntry(h_lambda_corrected_Ck, 'C_{full model}(k*)', 'l')
    legend.AddEntry(h_lambda_Sigma_corrected_Ck, 'C_{full model + #Sigma}(k*)', 'l')
    legend.AddEntry(h_lambda_Sigma_smeared_Ck, 'C_{full model + #Sigma + resolution}(k*)', 'l')
    legend.Draw()

    canvas.Write()

    del canvas

    canvas = TCanvas(f'cLambdaModel_{sign}_SystematicChange', f'cLambdaModel_{sign}_SystematicChange', 800, 600)
    hframe = canvas.DrawFrame(0.01, 0., 0.4, 1.08, '; #it{k}* (GeV/c); C(#it{k}*)')
    h_lambda_corrected_Ck_higher_lambda.Draw('HIST SAME')
    h_lambda_corrected_Ck.Draw('HIST SAME')
    h_lambda_corrected_Ck_lower_lambda.Draw('HIST SAME')

    legend = init_legend(0.6, 0.2, 0.88, 0.4, border_size=0, fill_style=0)
    legend.AddEntry(h_lambda_corrected_Ck_higher_lambda, f'#lambda_{{nominal}} + {modification_factor*100.:.0f}%', 'l')
    legend.AddEntry(h_lambda_corrected_Ck, '#lambda_{nominal}', 'l')
    legend.AddEntry(h_lambda_corrected_Ck_lower_lambda, f'#lambda_{{nominal}} - {modification_factor*100.:.0f}%', 'l')
    legend.Draw()

    canvas.Write()

    del canvas

    canvas = TCanvas(f'cLambdaSigmaModel_{sign}_SystematicChange', f'cLambdaSigmaModel_{sign}_SystematicChange', 800, 600)
    hframe = canvas.DrawFrame(0.01, 0., 0.4, 1.08, '; #it{k}* (GeV/c); C(#it{k}*)')
    h_lambda_Sigma_corrected_Ck_higher_lambda.Draw('HIST SAME')
    h_lambda_Sigma_corrected_Ck.Draw('HIST SAME')
    h_lambda_Sigma_corrected_Ck_lower_lambda.Draw('HIST SAME')

    legend = init_legend(0.6, 0.2, 0.88, 0.4, border_size=0, fill_style=0)
    legend.AddEntry(h_lambda_Sigma_corrected_Ck_higher_lambda, f'#lambda_{{nominal}} + {modification_factor*100.:.0f}%', 'l')
    legend.AddEntry(h_lambda_Sigma_corrected_Ck, '#lambda_{nominal}', 'l')
    legend.AddEntry(h_lambda_Sigma_corrected_Ck_lower_lambda, f'#lambda_{{nominal}} - {modification_factor*100.:.0f}%', 'l')
    legend.Draw()

    canvas.Write()

    del canvas

def produce_lambda_models_with_variations(sign: str, centrality: str, outdir: TDirectory, resolution_fits: dict):

    if not INPUT_LAMBDA_PARAMETER_PATH:
        raise ValueError("INPUT_LAMBDA_PARAMETER_PATH is not set. Please set it before calling this function.")
    if not LAMBDA_MODIFICATION_FACTOR:
        raise ValueError("LAMBDA_MODIFICATION_FACTOR is not set. Please set it before calling this function.")
    if not INPUT_CK_VARIATIONS or not INPUT_SIGMA_CK_VARIATIONS:
        raise ValueError("INPUT_CK_VARIATIONS and INPUT_SIGMA_CK_VARIATIONS must be set. Please set them before calling this function.")
    
    # only 010 and 1050 are computed 
    centrality_dir = 'centrality_0_10' if centrality == '010' else 'centrality_10_50'
    h_lambda_parameter = load_hist(INPUT_LAMBDA_PARAMETER_PATH, f'{centrality_dir}/{sign}/hLambdaParameters')
    h_lambda_Sigma_parameter = load_hist(INPUT_LAMBDA_PARAMETER_PATH, f'{centrality_dir}/{sign}/hLambdaSigmaParameters')
    
    hist_radius_variation = {}

    for variation_name in INPUT_CK_VARIATIONS[centrality]:

        variation_dir = outdir.mkdir(variation_name)

        h_theoretical_Ck       = load_hist(INPUT_CK_VARIATIONS[centrality][variation_name])
        h_theoretical_Sigma_Ck = load_hist(INPUT_SIGMA_CK_VARIATIONS[centrality][variation_name])

        h_lambda_Sigma_corrected_Ck = h_theoretical_Ck.Clone('hLambdaSigmaCorrectedCk')
        apply_lambda_correction(h_lambda_Sigma_corrected_Ck, h_theoretical_Ck, h_theoretical_Sigma_Ck,
                                h_lambda_parameter, h_lambda_Sigma_parameter)

        *_, h_higher, h_lower = produce_lambda_with_modified_values(
            sign, centrality, outdir=variation_dir,
            h_theoretical_Ck=h_theoretical_Ck, h_theoretical_Sigma_Ck=h_theoretical_Sigma_Ck,
            modification_factor=LAMBDA_MODIFICATION_FACTOR
        )
        
        hist_radius_variation[variation_name] = h_lambda_Sigma_corrected_Ck

        smeared = {
            f'hLambdaSigmaCorrectedCk_Smeared_{variation_name}':        apply_resolution_smearing(h_lambda_Sigma_corrected_Ck, outdir, resolution_fits),
            f'hLambdaSigmaCorrectedCk_Smeared_{variation_name}_higher': apply_resolution_smearing(h_higher, outdir, resolution_fits),
            f'hLambdaSigmaCorrectedCk_Smeared_{variation_name}_lower':  apply_resolution_smearing(h_lower, outdir, resolution_fits),
        }

        variation_dir.cd()
        for name, hist in smeared.items():
            set_root_object(hist, title='; #it{k}* (MeV/c); C(#it{k}*)', line_width=2)
            hist.Write(name)
            
    smeared_hists = {}
    for variation_name in INPUT_CK_VARIATIONS[centrality]:
        variation_dir = outdir.Get(variation_name)
        hist_name = f'hLambdaSigmaCorrectedCk_Smeared_{variation_name}'
        h = variation_dir.Get(hist_name)
        if h:
            h.SetDirectory(0)
            smeared_hists[variation_name] = h

    if hist_radius_variation:
        colors = {'nominal': get_color(0), 'lower': get_color(1), 'upper': get_color(2)}
        labels = {'nominal': '#it{R}_{s}', 'lower': '#it{R}_{s} - #sigma_{R}', 'upper': '#it{R}_{s} + #sigma_{R}'}

        canvas_summary = TCanvas(f'cRadiusVariations_{sign}_{centrality}',
                                 f'Radius variations {sign} {centrality}', 800, 600)
        hframe = canvas_summary.DrawFrame(0.01, 0., 0.4, 1.08,
                                          f'{sign} {centrality}; #it{{k}}* (GeV/#it{{c}}); C(#it{{k}}*)')
        legend_summary = TLegend(0.55, 0.65, 0.88, 0.88)
        legend_summary.SetBorderSize(0)
        legend_summary.SetFillStyle(0)

        for variation_name, h in hist_radius_variation.items():
            set_root_object(h, line_color=colors.get(variation_name, 1), line_width=2)
            h.Draw('HIST SAME')
            legend_summary.AddEntry(h, labels.get(variation_name, variation_name), 'l')

        legend_summary.Draw()
        outdir.cd()
        canvas_summary.Write(f'cRadiusVariations_{sign}_{centrality}')

def produce_lambda_models_with_model_param_variations(sign: str, centrality: str, outdir: TDirectory, resolution_fits: dict):

    if not INPUT_LAMBDA_PARAMETER_PATH:
        raise ValueError("INPUT_LAMBDA_PARAMETER_PATH is not set. Please set it before calling this function.")
    if not INPUT_CK_MODEL_PARAM_VARIATIONS:
        raise ValueError("INPUT_CK_MODEL_PARAM_VARIATIONS must be set. Please set it before calling this function.")

    centrality_dir = 'centrality_0_10' if centrality == '010' else 'centrality_10_50'
    h_lambda_parameter = load_hist(INPUT_LAMBDA_PARAMETER_PATH, f'{centrality_dir}/{sign}/hLambdaParameters')
    h_lambda_Sigma_parameter = load_hist(INPUT_LAMBDA_PARAMETER_PATH, f'{centrality_dir}/{sign}/hLambdaSigmaParameters')

    # Sigma contamination does not depend on the model parameters, so the nominal
    # Sigma Ck is reused for every variation.
    h_theoretical_Sigma_Ck = load_hist(INPUT_SIGMA_CK_PATH[centrality])

    for variation_name in INPUT_CK_MODEL_PARAM_VARIATIONS[centrality]:

        variation_dir = outdir.mkdir(variation_name)

        h_theoretical_Ck = load_hist(INPUT_CK_MODEL_PARAM_VARIATIONS[centrality][variation_name])

        h_lambda_Sigma_corrected_Ck = h_theoretical_Ck.Clone('hLambdaSigmaCorrectedCk')
        apply_lambda_correction(h_lambda_Sigma_corrected_Ck, h_theoretical_Ck, h_theoretical_Sigma_Ck,
                                h_lambda_parameter, h_lambda_Sigma_parameter)

        h_lambda_Sigma_smeared_Ck = apply_resolution_smearing(h_lambda_Sigma_corrected_Ck, variation_dir, resolution_fits)

        variation_dir.cd()
        set_root_object(h_lambda_Sigma_corrected_Ck, title='; #it{k}* (MeV/c); C(#it{k}*)', line_width=2)
        set_root_object(h_lambda_Sigma_smeared_Ck, title='; #it{k}* (MeV/c); C(#it{k}*)', line_width=2)
        h_lambda_Sigma_corrected_Ck.Write('hLambdaSigmaCorrectedCk')
        h_lambda_Sigma_smeared_Ck.Write(f'hLambdaSigmaCorrectedCk_Smeared_{variation_name}')


if __name__ == '__main__':
    
    RooMsgService.instance().setGlobalKillBelow(5) # 3 = WARNING, 4 = ERROR, 5 = FATAL
    RooFit.PrintLevel(-1)
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', default='configs/lambda_model.yaml',
                        help='Path to YAML config file')
    args, _ = parser.parse_known_args()

    cfg = load_yaml(args.config)

    LAMBDA_MODIFICATION_FACTOR = cfg['lambda_modification_factor']
    LAMBDA_VARIATION = cfg['lambda_variation']

    INPUT_CK_PATH = build_hist_load_info_dict(cfg['ck_input'])
    INPUT_SIGMA_CK_PATH = build_hist_load_info_dict(cfg['sigma_ck_input'])
    INPUT_CK_VARIATIONS = build_hist_load_info_variations_radius(cfg['ck_input']) if 'variations_radius' in cfg['ck_input'] else None
    INPUT_CK_MODEL_PARAM_VARIATIONS = build_hist_load_info_variations_model_params(cfg['ck_input']) if 'variations_model_params' in cfg['ck_input'] else None
    INPUT_SIGMA_CK_VARIATIONS = build_hist_load_info_variations_radius(cfg['sigma_ck_input']) if 'variations_radius' in cfg['sigma_ck_input'] else None

    INPUT_LAMBDA_PARAMETER_PATH = cfg['paths']['lambda_parameters']
    INPUT_RESOLUTION_PATH = cfg['paths']['resolution']
    INPUT_MIXED_EVENT_REFERENCE_PATH = cfg['paths']['mixed_event_reference']
    OUTPUT_LAMBDA_MODEL_PATH = cfg['paths']['output_model']
    
    CENTRALITY_BINS = cfg['centrality_bins']

    outfile = TFile.Open(OUTPUT_LAMBDA_MODEL_PATH, 'recreate')
    resolution_fits = precompute_resolution_fits(outfile)

    for sign in ['Both', 'Matter', 'Antimatter']:
        print(f"\n{'='*60}")
        print(f"Processing {sign}")

        outdir_sign = outfile.mkdir(sign)

        for centrality in CENTRALITY_BINS:
        #for centrality in ['010', '1030', '3050']:
            print(f"\n{'-'*40}")
            print(f"Processing centrality {centrality}")

            outdir = outdir_sign.mkdir(f'{centrality}')
            produce_lambda_models(sign, centrality, outdir, resolution_fits)

            if INPUT_CK_VARIATIONS and INPUT_SIGMA_CK_VARIATIONS:
                produce_lambda_models_with_variations(sign, centrality, outdir, resolution_fits)
            
            if INPUT_CK_MODEL_PARAM_VARIATIONS:
                produce_lambda_models_with_model_param_variations(sign, centrality, outdir, resolution_fits)
            
    outfile.Close()