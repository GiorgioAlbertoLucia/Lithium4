'''
    Create a model that includes the presence of lambda parameters.
'''

import numpy as np

from ROOT import TFile, TDirectory, TCanvas, RooRealVar, RooFit, RooDataHist, TF1

from torchic.core.histogram import load_hist
from torchic.roopdf.roopdf_utils import init_roopdf
from torchic.utils.root import set_root_object, silence_roofit

LAMBDA_MODIFICATION_FACTOR = 0.1  # 10% change in lambda

RADIUS_VALUES = np.arange(3.5, 8.0 + 0.05, 0.1)
CENTRALITIES_TO_SCAN = ['010', '1050']
### USE_LL = False
### CK_SCAN_PATH = '/home/galucia/phemto/output/pHe3_square_well_scan_rescaled.root'
USE_LL = True
CK_SCAN_PATH = '/home/galucia/PhaseShiftAnalysis/numerical_lednicky/pHe/output/pHe3_LL_radius_span_GeV.root'
SIGMA_CK_SCAN_PATH = '/home/galucia/phemto/output/he3Sigma_scan_rescaled.root'

LAMBDA_VARIATION = 0.10 # typical is 10%, we use 15% 
INPUT_LAMBDA_PARAMETER_PATH = '/home/galucia/Lithium4/calibration/output/lambda_parameters.root'

INPUT_RESOLUTION_PATH = '/data/galucia/lithium/MC/AnalysisResults_LHC25g11.root'
INPUT_MIXED_EVENT_REFERENCE_PATH = '/home/galucia/Lithium4/preparation/output/PbPb/LHC25_PbPb_pass1_hadronpid_event_mixing.root'

OUTPUT_LAMBDA_MODEL_PATH = '/home/galucia/Lithium4/femto/models/lambda_models_scan_LL.root'

def get_ck_scan_name(radius: float, use_LL: bool) -> str:
    if not use_LL:
        return f'hcats_CF_r={radius:.2f}_fm'
    else:
        radius_LL = radius / np.sqrt(2)  # Convert to LL radius
        return f'Rs{radius_LL:.2f}/combined/combined_h{radius_LL:.2f}_aRe-11.26_aIm0.00_r1.65_GeV_h{radius_LL:.2f}_aRe-9.06_aIm0.00_r1.36_GeV'

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

def apply_lambda_correction(h_out, h_ck, h_sigma_ck, lambda_param_hist, lambda_sigma_hist, scale: float = 1.0):
    for ibin in range(1, h_out.GetNbinsX() + 1):
        kstar = h_out.GetBinCenter(ibin)
        lam = lambda_param_hist.GetBinContent(lambda_param_hist.FindBin(kstar)) * scale
        lam = min(max(lam, 0.), 1.)  # Ensure lambda is between 0 and 1
        lam_s = lambda_sigma_hist.GetBinContent(lambda_sigma_hist.FindBin(kstar))
        ck = h_ck.GetBinContent(ibin)
        sig = h_sigma_ck.GetBinContent(h_sigma_ck.FindBin(kstar))
        h_out.SetBinContent(ibin, lam * ck + lam_s * sig + (1 - lam - lam_s))

def precompute_resolution_fits(outfile: TFile) -> dict:
    """
    Fit each kstar slice of the resolution matrix with a RooFit Crystal Ball PDF.
    Returns a dict mapping kstar bin index (1-based) -> callable(x) using the fitted PDF.
    Saves all fits to outfile under 'ResolutionFits/'.
    """
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

def apply_resolution_smearing(h_correlation_function, resolution_fits:dict, outdir:TDirectory=None):

    h_resolution = load_hist(INPUT_RESOLUTION_PATH, 'he3-hadron-femto/QA/hKstarRecVsKstarGen')
    h_mixed_event = load_hist(INPUT_MIXED_EVENT_REFERENCE_PATH, 'QA/hKstar')
    mixed_fit = TF1('mixed_fit', 'pol3', 0.01, 0.4)
    h_mixed_event.Fit(mixed_fit, 'RMSQ+')
    
    # match the binning of the resolution histogram to that of the correlation function
    
    h_resolution_reference = h_resolution.ProjectionX('hResolutionReference', 2, 2)
    h_correlation_function_matched = h_correlation_function.Clone(f'{h_correlation_function.GetName()}_matched_for_smearing')
    for ibin in range(1, h_correlation_function_matched.GetNbinsX()+1):
        if h_correlation_function_matched.GetBinCenter(ibin) < 0.01: 
            h_correlation_function_matched.SetBinContent(ibin, 0.)

    h_correlation_function_matched = match_bin_width_correlation_function(h_resolution_reference, h_correlation_function_matched, kstar_threshold=0.7)

    #outdir_check = outdir.mkdir('ResolutionSmearingCheck')
    #outdir_check.cd()
    #h_resolution_reference.Write('hResolutionReference')
    #h_correlation_function_matched.Write('hCorrelationFunction_MatchedForSmearing')
    #h_mixed_event.Write('hMixedEventForSmearing')

    h_smeared_correlation_function = h_correlation_function_matched.Clone(f'{h_correlation_function_matched.GetName()}_smeared')

    for ibin in range(1, h_correlation_function_matched.GetNbinsX()+1):

        smeared_value, weight, total_weight = 0., 0., 0.
        kstar = h_correlation_function_matched.GetBinCenter(ibin)
        resolution_bin = h_resolution.GetXaxis().FindBin(kstar)
        h_resolution_slice = h_resolution.ProjectionX(f'hResolutionSlice_kstar_{kstar:.3f}', resolution_bin, resolution_bin)

        #outdir_check.cd()
        #h_resolution_slice.Write()
        
        for jbin in range(1, h_resolution_slice.GetNbinsX()+1):
            
            kstar_gen = h_resolution_slice.GetBinCenter(jbin)
            #mixed_weight = h_mixed_event.GetBinContent(h_mixed_event.FindBin(kstar_gen))
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
            
            #weight = h_resolution_slice.GetBinContent(jbin) * mixed_weight if 0.01 < kstar_gen < 0.7 else 0.  
            
            # skip the region where the corrected correlation function is not defined
            correlation_value = h_correlation_function_matched.GetBinContent(jbin)
            
            smeared_value += correlation_value * weight
            total_weight += weight

        if total_weight > 0:
            smeared_value /= total_weight
            h_smeared_correlation_function.SetBinContent(ibin, smeared_value)


    return h_smeared_correlation_function

def produce_lambda_with_modified_values(sign:str, centrality:str, h_theoretical_ck, h_theoretical_sigma_ck,
                                        modification_factor:float):
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

    h_lambda_corrected_Ck_higher_lambda = h_theoretical_ck.Clone(f'hLambdaCorrectedCk_Higher')
    apply_lambda_correction(h_lambda_corrected_Ck_higher_lambda, h_theoretical_ck, h_theoretical_sigma_ck,
                            h_lambda_parameter, h_lambda_Sigma_parameter, scale=1 + modification_factor)
    h_lambda_corrected_Ck_lower_lambda = h_theoretical_ck.Clone(f'hLambdaCorrectedCk_Lower')
    apply_lambda_correction(h_lambda_corrected_Ck_lower_lambda, h_theoretical_ck, h_theoretical_sigma_ck,
                            h_lambda_parameter, h_lambda_Sigma_parameter, scale=1 - modification_factor)

    h_lambda_Sigma_corrected_Ck_higher_lambda = h_theoretical_ck.Clone(f'hLambdaSigmaCorrectedCk_Higher')
    apply_lambda_correction(h_lambda_Sigma_corrected_Ck_higher_lambda, h_theoretical_ck, h_theoretical_sigma_ck,
                            h_lambda_parameter, h_lambda_Sigma_parameter, scale=1 + modification_factor)
    h_lambda_Sigma_corrected_Ck_lower_lambda = h_theoretical_ck.Clone(f'hLambdaSigmaCorrectedCk_Lower')
    apply_lambda_correction(h_lambda_Sigma_corrected_Ck_lower_lambda, h_theoretical_ck, h_theoretical_sigma_ck,
                            h_lambda_parameter, h_lambda_Sigma_parameter, scale=1 - modification_factor)

    return (h_lambda_parameter_higher_lambda, h_lambda_parameter_lower_lambda,
            h_lambda_corrected_Ck_lower_lambda, h_lambda_corrected_Ck_higher_lambda, 
            h_lambda_Sigma_corrected_Ck_lower_lambda, h_lambda_Sigma_corrected_Ck_higher_lambda)


def produce_lambda_models_radius_scan(sign: str, centrality: str, outdir: TDirectory, resolution_fits: dict):

    centrality_dir = 'centrality_0_10' if centrality == '010' else 'centrality_10_50'

    h_lambda_parameter = load_hist(INPUT_LAMBDA_PARAMETER_PATH, f'{centrality_dir}/{sign}/hLambdaParameters')
    h_lambda_sigma_parameter = load_hist(INPUT_LAMBDA_PARAMETER_PATH, f'{centrality_dir}/{sign}/hLambdaSigmaParameters')

    for radius in RADIUS_VALUES:

        radius_str = f'{radius:.2f}'
        print(f'Processing {sign}, centrality {centrality}, R = {radius_str} fm')

        h_theoretical_ck = load_hist(CK_SCAN_PATH, get_ck_scan_name(radius, USE_LL))
        h_theoretical_sigma_ck = load_hist(SIGMA_CK_SCAN_PATH, f'hhe3_Sigma_plus_CF_r={radius_str}_fm')

        # Apply lambda correction and combine the contributions
        h_lambda_sigma_corrected_ck = h_theoretical_ck.Clone(f'hLambdaSigmaCorrectedCk_R_{radius_str}')
        apply_lambda_correction(h_lambda_sigma_corrected_ck, h_theoretical_ck, h_theoretical_sigma_ck, h_lambda_parameter, h_lambda_sigma_parameter)
        h_lambda_sigma_smeared_ck = apply_resolution_smearing(h_lambda_sigma_corrected_ck, resolution_fits)

        # Set plotting style
        set_root_object(h_theoretical_ck, title='; #it{k}* (GeV/#it{c}); C(#it{k}*)', line_width=2)
        set_root_object(h_theoretical_sigma_ck, title='; #it{k}* (GeV/#it{c}); C(#it{k}*)', line_width=2)
        set_root_object(h_lambda_sigma_corrected_ck, title='; #it{k}* (GeV/#it{c}); C(#it{k}*)', line_width=2)
        set_root_object(h_lambda_sigma_smeared_ck, title='; #it{k}* (GeV/#it{c}); C(#it{k}*)', line_width=2)

        # Write everything
        outdir.cd(f'R_{radius_str}')
        h_theoretical_ck.Write(f'hTheoreticalCk_R_{radius_str}')
        h_theoretical_sigma_ck.Write(f'hTheoreticalSigmaCk_R_{radius_str}')
        h_lambda_sigma_corrected_ck.Write(f'hLambdaSigmaCorrectedCk_R_{radius_str}')
        h_lambda_sigma_smeared_ck.Write(f'hLambdaSigmaCorrectedCk_Smeared_R_{radius_str}')

def produce_lambda_models_with_variations(sign: str, centrality: str, outdir: TDirectory, resolution_fits: dict):

    if not INPUT_LAMBDA_PARAMETER_PATH:
        raise ValueError("INPUT_LAMBDA_PARAMETER_PATH is not set. Please set it before calling this function.")
    if not LAMBDA_MODIFICATION_FACTOR:
        raise ValueError("LAMBDA_MODIFICATION_FACTOR is not set. Please set it before calling this function.")
    
    # only 010 and 1050 are computed 
    centrality_dir = 'centrality_0_10' if centrality == '010' else 'centrality_10_50'
    h_lambda_parameter = load_hist(INPUT_LAMBDA_PARAMETER_PATH, f'{centrality_dir}/{sign}/hLambdaParameters')
    h_lambda_sigma_parameter = load_hist(INPUT_LAMBDA_PARAMETER_PATH, f'{centrality_dir}/{sign}/hLambdaSigmaParameters')
    
    for radius in RADIUS_VALUES:

        radius_str = f'{radius:.2f}'
        print(f'Processing {sign}, centrality {centrality}, R = {radius_str} fm')

        outdir.cd(f'R_{radius_str}')
        
        h_theoretical_ck = load_hist(CK_SCAN_PATH, get_ck_scan_name(radius, USE_LL))
        h_theoretical_sigma_ck = load_hist(SIGMA_CK_SCAN_PATH, f'hhe3_Sigma_plus_CF_r={radius_str}_fm')
    
        h_lambda_sigma_corrected_ck = h_theoretical_ck.Clone(f'hLambdaSigmaCorrectedCk_R_{radius_str}')
        apply_lambda_correction(h_lambda_sigma_corrected_ck, h_theoretical_ck, h_theoretical_sigma_ck,
                                h_lambda_parameter, h_lambda_sigma_parameter)

        *_, h_higher, h_lower = produce_lambda_with_modified_values(
            sign, centrality,
            h_theoretical_ck=h_theoretical_ck, h_theoretical_sigma_ck=h_theoretical_sigma_ck,
            modification_factor=LAMBDA_MODIFICATION_FACTOR
        )

        smeared = {
            f'hLambdaSigmaCorrectedCk_Smeared':        apply_resolution_smearing(h_lambda_sigma_corrected_ck, resolution_fits, outdir),
            f'hLambdaSigmaCorrectedCk_Smeared_higher': apply_resolution_smearing(h_higher, resolution_fits, outdir),
            f'hLambdaSigmaCorrectedCk_Smeared_lower':  apply_resolution_smearing(h_lower, resolution_fits, outdir),
        }

        outdir.mkdir(f'R_{radius_str}/lambda_variations')
        outdir.cd(f'R_{radius_str}/lambda_variations')
        for name, hist in smeared.items():
            set_root_object(hist, title='; #it{k}* (MeV/c); C(#it{k}*)', line_width=2)
            hist.Write(name)

if __name__ == '__main__':

    silence_roofit()
    
    outfile = TFile.Open(OUTPUT_LAMBDA_MODEL_PATH, 'recreate')
    resolution_fits = precompute_resolution_fits(outfile)

    for sign in ['Both', 'Matter', 'Antimatter']:
        print(f"\n{'='*60}")
        print(f"Processing {sign}")

        outdir_sign = outfile.mkdir(sign)

        for centrality in CENTRALITIES_TO_SCAN:

            print(f"\n{'-'*40}")
            print(f"Processing centrality {centrality}")

            outdir = outdir_sign.mkdir(centrality)
            for radius in RADIUS_VALUES:
                radius_str = f'{radius:.2f}'
                outdir.mkdir(f'R_{radius_str}')

            produce_lambda_models_radius_scan(sign, centrality, outdir, resolution_fits)
            produce_lambda_models_with_variations(sign, centrality, outdir, resolution_fits)
            
    outfile.Close()