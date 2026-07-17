'''
    Create a model that includes the presence of lambda parameters.
'''

import numpy as np

from ROOT import TFile, TDirectory, TCanvas, TLegend, TTree, TF1, \
                 RooDataSet, RooKeysPdf, RooRealVar, RooFit, RooDataHist

from torchic.core.histogram import load_hist, HistLoadInfo
from torchic.utils.root import set_root_object, init_legend
from torchic.utils.colors import get_color
from torchic.roopdf.roopdf_utils import init_roopdf

LAMBDA_MODIFICATION_FACTOR = 0.1  # 10% change in lambda

INPUT_CK_PATH = {
    '010':  HistLoadInfo('/home/galucia/phemto/output/pHe3_square_well_recomputed_rescaled.root', 'hcats_CF_r=6.12_fm'),
    '1030': HistLoadInfo('/home/galucia/phemto/output/pHe3_square_well_recomputed_rescaled.root', 'hcats_CF_r=5.18_fm'),
    '3050': HistLoadInfo('/home/galucia/phemto/output/pHe3_square_well_recomputed_rescaled.root', 'hcats_CF_r=4.09_fm'),
    '5080': HistLoadInfo('/home/galucia/phemto/output/pHe3_square_well_recomputed_rescaled.root', 'hcats_CF_r=3.05_fm'),
    '050':  HistLoadInfo('/home/galucia/phemto/output/pHe3_square_well_recomputed_rescaled.root', 'hcats_CF_r=5.52_fm'),
    '080':  HistLoadInfo('/home/galucia/phemto/output/pHe3_square_well_recomputed_rescaled.root', 'hcats_CF_r=5.50_fm'),
    '1050': HistLoadInfo('/home/galucia/phemto/output/pHe3_square_well_recomputed_rescaled.root', 'hcats_CF_r=4.90_fm'),
    '1080': HistLoadInfo('/home/galucia/phemto/output/pHe3_square_well_recomputed_rescaled.root', 'hcats_CF_r=4.91_fm'),
}
INPUT_SIGMA_CK_PATH = {
    '010':  HistLoadInfo('/home/galucia/phemto/output/he3Sigma_recomputed_rescaled.root', 'hhe3_Sigma_plus_CF_r=6.12_fm'),
    '1030': HistLoadInfo('/home/galucia/phemto/output/he3Sigma_recomputed_rescaled.root', 'hhe3_Sigma_plus_CF_r=5.18_fm'),
    '3050': HistLoadInfo('/home/galucia/phemto/output/he3Sigma_recomputed_rescaled.root', 'hhe3_Sigma_plus_CF_r=4.09_fm'),
    '5080': HistLoadInfo('/home/galucia/phemto/output/he3Sigma_recomputed_rescaled.root', 'hhe3_Sigma_plus_CF_r=3.05_fm'),
    '050':  HistLoadInfo('/home/galucia/phemto/output/he3Sigma_recomputed_rescaled.root', 'hhe3_Sigma_plus_CF_r=5.52_fm'),
    '080':  HistLoadInfo('/home/galucia/phemto/output/he3Sigma_recomputed_rescaled.root', 'hhe3_Sigma_plus_CF_r=5.50_fm'),
    '1050': HistLoadInfo('/home/galucia/phemto/output/he3Sigma_recomputed_rescaled.root', 'hhe3_Sigma_plus_CF_r=4.90_fm'),
    '1080': HistLoadInfo('/home/galucia/phemto/output/he3Sigma_recomputed_rescaled.root', 'hhe3_Sigma_plus_CF_r=4.91_fm'),
}
INPUT_CK_VARIATIONS = {
    '010': {
        'nominal': HistLoadInfo('/home/galucia/phemto/output/pHe3_square_well_recomputed_rescaled.root', 'hcats_CF_r=6.12_fm'),
        'lower':   HistLoadInfo('/home/galucia/phemto/output/pHe3_square_well_recomputed_rescaled.root', 'hcats_CF_r=5.44_fm'),
        'upper':   HistLoadInfo('/home/galucia/phemto/output/pHe3_square_well_recomputed_rescaled.root', 'hcats_CF_r=6.80_fm'),
    },
    '1030': {
        'nominal': HistLoadInfo('/home/galucia/phemto/output/pHe3_square_well_recomputed_rescaled.root', 'hcats_CF_r=5.18_fm'),
        'lower':   HistLoadInfo('/home/galucia/phemto/output/pHe3_square_well_recomputed_rescaled.root', 'hcats_CF_r=4.57_fm'),
        'upper':   HistLoadInfo('/home/galucia/phemto/output/pHe3_square_well_recomputed_rescaled.root', 'hcats_CF_r=5.80_fm'),
    },
    '3050': {
        'nominal': HistLoadInfo('/home/galucia/phemto/output/pHe3_square_well_recomputed_rescaled.root', 'hcats_CF_r=4.09_fm'),
        'lower':   HistLoadInfo('/home/galucia/phemto/output/pHe3_square_well_recomputed_rescaled.root', 'hcats_CF_r=3.54_fm'),
        'upper':   HistLoadInfo('/home/galucia/phemto/output/pHe3_square_well_recomputed_rescaled.root', 'hcats_CF_r=4.64_fm'),
    },
    '5080': {
        'nominal': HistLoadInfo('/home/galucia/phemto/output/pHe3_square_well_recomputed_rescaled.root', 'hcats_CF_r=3.05_fm'),
        'lower':   HistLoadInfo('/home/galucia/phemto/output/pHe3_square_well_recomputed_rescaled.root', 'hcats_CF_r=2.56_fm'),
        'upper':   HistLoadInfo('/home/galucia/phemto/output/pHe3_square_well_recomputed_rescaled.root', 'hcats_CF_r=3.54_fm'),
    },
    '050': {
        'nominal': HistLoadInfo('/home/galucia/phemto/output/pHe3_square_well_recomputed_rescaled.root', 'hcats_CF_r=5.52_fm'),
        'lower':   HistLoadInfo('/home/galucia/phemto/output/pHe3_square_well_recomputed_rescaled.root', 'hcats_CF_r=5.14_fm'),
        'upper':   HistLoadInfo('/home/galucia/phemto/output/pHe3_square_well_recomputed_rescaled.root', 'hcats_CF_r=5.90_fm'),
    },
    '080': {
        'nominal': HistLoadInfo('/home/galucia/phemto/output/pHe3_square_well_recomputed_rescaled.root', 'hcats_CF_r=5.50_fm'),
        'lower':   HistLoadInfo('/home/galucia/phemto/output/pHe3_square_well_recomputed_rescaled.root', 'hcats_CF_r=5.12_fm'),
        'upper':   HistLoadInfo('/home/galucia/phemto/output/pHe3_square_well_recomputed_rescaled.root', 'hcats_CF_r=5.88_fm'),
    },
    '1050': {
        'nominal': HistLoadInfo('/home/galucia/phemto/output/pHe3_square_well_recomputed_rescaled.root', 'hcats_CF_r=4.90_fm'),
        'lower':   HistLoadInfo('/home/galucia/phemto/output/pHe3_square_well_recomputed_rescaled.root', 'hcats_CF_r=4.30_fm'),
        'upper':   HistLoadInfo('/home/galucia/phemto/output/pHe3_square_well_recomputed_rescaled.root', 'hcats_CF_r=5.50_fm'),
    },
    '1080': {
        'nominal': HistLoadInfo('/home/galucia/phemto/output/pHe3_square_well_recomputed_rescaled.root', 'hcats_CF_r=4.91_fm'),
        'lower':   HistLoadInfo('/home/galucia/phemto/output/pHe3_square_well_recomputed_rescaled.root', 'hcats_CF_r=4.42_fm'),
        'upper':   HistLoadInfo('/home/galucia/phemto/output/pHe3_square_well_recomputed_rescaled.root', 'hcats_CF_r=5.40_fm'),
    },
}
INPUT_SIGMA_CK_VARIATIONS = {
    '010': {
        'nominal': HistLoadInfo('/home/galucia/phemto/output/he3Sigma_recomputed_rescaled.root', 'hhe3_Sigma_plus_CF_r=6.12_fm'),
        'lower':   HistLoadInfo('/home/galucia/phemto/output/he3Sigma_recomputed_rescaled.root', 'hhe3_Sigma_plus_CF_r=5.44_fm'),
        'upper':   HistLoadInfo('/home/galucia/phemto/output/he3Sigma_recomputed_rescaled.root', 'hhe3_Sigma_plus_CF_r=6.80_fm'),
    },
    '1030': {
        'nominal': HistLoadInfo('/home/galucia/phemto/output/he3Sigma_recomputed_rescaled.root', 'hhe3_Sigma_plus_CF_r=5.18_fm'),
        'lower':   HistLoadInfo('/home/galucia/phemto/output/he3Sigma_recomputed_rescaled.root', 'hhe3_Sigma_plus_CF_r=4.57_fm'),
        'upper':   HistLoadInfo('/home/galucia/phemto/output/he3Sigma_recomputed_rescaled.root', 'hhe3_Sigma_plus_CF_r=5.80_fm'),
    },
    '3050': {
        'nominal': HistLoadInfo('/home/galucia/phemto/output/he3Sigma_recomputed_rescaled.root', 'hhe3_Sigma_plus_CF_r=4.09_fm'),
        'lower':   HistLoadInfo('/home/galucia/phemto/output/he3Sigma_recomputed_rescaled.root', 'hhe3_Sigma_plus_CF_r=3.54_fm'),
        'upper':   HistLoadInfo('/home/galucia/phemto/output/he3Sigma_recomputed_rescaled.root', 'hhe3_Sigma_plus_CF_r=4.64_fm'),
    },
    '5080': {
        'nominal': HistLoadInfo('/home/galucia/phemto/output/he3Sigma_recomputed_rescaled.root', 'hhe3_Sigma_plus_CF_r=3.05_fm'),
        'lower':   HistLoadInfo('/home/galucia/phemto/output/he3Sigma_recomputed_rescaled.root', 'hhe3_Sigma_plus_CF_r=2.56_fm'),
        'upper':   HistLoadInfo('/home/galucia/phemto/output/he3Sigma_recomputed_rescaled.root', 'hhe3_Sigma_plus_CF_r=3.54_fm'),
    },
    '050': {
        'nominal': HistLoadInfo('/home/galucia/phemto/output/he3Sigma_recomputed_rescaled.root', 'hhe3_Sigma_plus_CF_r=5.52_fm'),
        'lower':   HistLoadInfo('/home/galucia/phemto/output/he3Sigma_recomputed_rescaled.root', 'hhe3_Sigma_plus_CF_r=5.14_fm'),
        'upper':   HistLoadInfo('/home/galucia/phemto/output/he3Sigma_recomputed_rescaled.root', 'hhe3_Sigma_plus_CF_r=5.90_fm'),
    },
    '080': {
        'nominal': HistLoadInfo('/home/galucia/phemto/output/he3Sigma_recomputed_rescaled.root', 'hhe3_Sigma_plus_CF_r=5.50_fm'),
        'lower':   HistLoadInfo('/home/galucia/phemto/output/he3Sigma_recomputed_rescaled.root', 'hhe3_Sigma_plus_CF_r=5.12_fm'),
        'upper':   HistLoadInfo('/home/galucia/phemto/output/he3Sigma_recomputed_rescaled.root', 'hhe3_Sigma_plus_CF_r=5.88_fm'),
    },
    '1050': {
        'nominal': HistLoadInfo('/home/galucia/phemto/output/he3Sigma_recomputed_rescaled.root', 'hhe3_Sigma_plus_CF_r=4.90_fm'),
        'lower':   HistLoadInfo('/home/galucia/phemto/output/he3Sigma_recomputed_rescaled.root', 'hhe3_Sigma_plus_CF_r=4.30_fm'),
        'upper':   HistLoadInfo('/home/galucia/phemto/output/he3Sigma_recomputed_rescaled.root', 'hhe3_Sigma_plus_CF_r=5.50_fm'),
    },
    '1080': {
        'nominal': HistLoadInfo('/home/galucia/phemto/output/he3Sigma_recomputed_rescaled.root', 'hhe3_Sigma_plus_CF_r=4.91_fm'),
        'lower':   HistLoadInfo('/home/galucia/phemto/output/he3Sigma_recomputed_rescaled.root', 'hhe3_Sigma_plus_CF_r=4.42_fm'),
        'upper':   HistLoadInfo('/home/galucia/phemto/output/he3Sigma_recomputed_rescaled.root', 'hhe3_Sigma_plus_CF_r=5.40_fm'),
    },
}
INPUT_LAMBDA_PARAMETER_PATH = '/home/galucia/Lithium4/calibration/output/LHC25_PbPb_pass1_lambda_parameters.root'
EXPERIMENTAL_CK_PATH = '/home/galucia/Lithium4/preparation/checks/correlation_hadronpid_pass1_pass4_refined_dca.root'
EXPERIMENTAL_CK_NAME = 'Correlation/Default/hCorrelation010'

INPUT_RESOLUTION_PATH = '/data/galucia/lithium/MC/AnalysisResults_LHC25g11.root'
INPUT_MIXED_EVENT_REFERENCE_PATH = '/home/galucia/Lithium4/preparation/output/PbPb/LHC25_PbPb_pass1_hadronpid_event_mixing.root'

OUTPUT_LAMBDA_MODEL_PATH = '/home/galucia/Lithium4/femto/models/LHC25_PbPb_pass1_lambda_models.root'

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

def apply_lambda_correction(h_out, h_ck, h_sigma_ck, lambda_param_hist, lambda_sigma_hist, scale: float = 1.0):
    for ibin in range(1, h_out.GetNbinsX() + 1):
        kstar = h_out.GetBinCenter(ibin)
        lam = lambda_param_hist.GetBinContent(lambda_param_hist.FindBin(kstar)) * scale
        lam = min(max(lam, 0.), 1.)  # Ensure lambda is between 0 and 1
        lam_s = lambda_sigma_hist.GetBinContent(lambda_sigma_hist.FindBin(kstar))
        ck = h_ck.GetBinContent(ibin)
        sig = h_sigma_ck.GetBinContent(h_sigma_ck.FindBin(kstar))
        h_out.SetBinContent(ibin, lam * ck + lam_s * sig + (1 - lam - lam_s))
        
def produce_lambda_with_modified_values(sign:str, centrality:str, h_theoretical_Ck, h_theoretical_Sigma_Ck,
                                        outdir:TDirectory, modification_factor:float):
    '''
        Produce a lambda-corrected model by systematically changing the lambda value by x%, 
        where x is the modification_factor. This can be used to understand the sensitivity of the model to changes in lambda.
        The lambda is changed both in lambda + x% and lambda - x% to understand the effect in both directions.
    '''

    h_lambda_parameter = load_hist(INPUT_LAMBDA_PARAMETER_PATH, f'{sign}/hLambdaParameters')
    h_lambda_Sigma_parameter = load_hist(INPUT_LAMBDA_PARAMETER_PATH, f'{sign}/hLambdaSigmaParameters')
    
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

    h_resolution = load_hist(INPUT_RESOLUTION_PATH, 'he3-hadron-femto/QA/hKstarRecVsKstarGen')
    h_mixed_event = load_hist(INPUT_MIXED_EVENT_REFERENCE_PATH, 'QA/hKstar')
    mixed_fit = TF1('mixed_fit', 'pol3', 0.01, 0.4)
    h_mixed_event.Fit(mixed_fit, 'RMS+')
    
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


def produce_lambda_models(sign:str, centrality:str, outdir:TDirectory, resolution_fits:dict):
    
    h_lambda_parameter = load_hist(INPUT_LAMBDA_PARAMETER_PATH, f'{sign}/hLambdaParameters')
    h_lambda_Sigma_parameter = load_hist(INPUT_LAMBDA_PARAMETER_PATH, f'{sign}/hLambdaSigmaParameters')
    h_theoretical_Ck = load_hist(INPUT_CK_PATH[centrality])
    h_theoretical_Sigma_Ck = load_hist(INPUT_SIGMA_CK_PATH[centrality])
    
    h_lambda_corrected_Ck = h_theoretical_Ck.Clone(f'hLambdaCorrectedCk')
    h_lambda_Sigma_corrected_Ck = h_theoretical_Ck.Clone(f'hLambdaSigmaCorrectedCk')

    for ibin in range(1, h_lambda_corrected_Ck.GetNbinsX()+1):
        
        kstar = h_lambda_corrected_Ck.GetBinCenter(ibin)
        lambda_param = h_lambda_parameter.GetBinContent(h_lambda_parameter.FindBin(kstar))
        original_value = h_theoretical_Ck.GetBinContent(ibin)

        lambda_Sigma_param = h_lambda_Sigma_parameter.GetBinContent(h_lambda_Sigma_parameter.FindBin(kstar))
        Sigma_value = h_theoretical_Sigma_Ck.GetBinContent(h_theoretical_Sigma_Ck.FindBin(kstar))
        
        lambda_corrected_value = lambda_param * original_value + (1.0 - lambda_param)
        h_lambda_corrected_Ck.SetBinContent(ibin, lambda_corrected_value)

        lambda_Sigma_corrected_value = lambda_param * original_value + lambda_Sigma_param * Sigma_value + (1 - lambda_param - lambda_Sigma_param) * 1.0
        h_lambda_Sigma_corrected_Ck.SetBinContent(ibin, lambda_Sigma_corrected_value)

    modification_factor = 0.1  # 10% change in lambda
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

    h_lambda_parameter = load_hist(INPUT_LAMBDA_PARAMETER_PATH, f'{sign}/hLambdaParameters')
    h_lambda_Sigma_parameter = load_hist(INPUT_LAMBDA_PARAMETER_PATH, f'{sign}/hLambdaSigmaParameters')
    
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

if __name__ == '__main__':

    outfile = TFile.Open(OUTPUT_LAMBDA_MODEL_PATH, 'recreate')
    resolution_fits = precompute_resolution_fits(outfile)

    for sign in ['Both', 'Matter', 'Antimatter']:
        print(f"\n{'='*60}")
        print(f"Processing {sign}")

        outdir_sign = outfile.mkdir(sign)

        for centrality in ['010', '1030', '3050', '5080', '050', '080', '1050', '1080']:
        #for centrality in ['010', '1030', '3050']:
            print(f"\n{'-'*40}")
            print(f"Processing centrality {centrality}")

            outdir = outdir_sign.mkdir(f'{centrality}')
            produce_lambda_models(sign, centrality, outdir, resolution_fits)

            produce_lambda_models_with_variations(sign, centrality, outdir, resolution_fits)
            
    outfile.Close()