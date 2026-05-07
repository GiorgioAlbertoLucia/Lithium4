'''
    Create a model that includes the presence of lambda parameters.
'''

import numpy as np

from alive_progress import alive_bar

from ROOT import TFile, TDirectory, TCanvas, TLegend, TKDE, TTree, TF1, \
                 RooDataSet, RooKeysPdf, RooRealVar, RooDataHist, \
                 kRed, kAzure, kGreen, kOrange, kViolet, kYellow, kCyan, kPink, kSpring, kMagenta, kTeal

from torchic.core.histogram import load_hist
from torchic.utils.root import set_root_object
from torchic.utils.colors import get_color

COLORS: list[int] = [
        kRed    + 1,   # bright red
        kAzure  + 1,   # rich azure-blue
        kGreen  + 2,   # vivid green
        kOrange + 1,   # warm orange
        kViolet + 1,   # purple-violet
        #kYellow - 7,   # golden yellow (base is pale; -7 gives a richer tone)
        kCyan   + 1,   # cyan-teal
        kPink   - 4,   # rose-pink
        kSpring + 4,   # yellow-green
        kMagenta+ 2,   # deep magenta
        kTeal   + 2,   # dark teal
        kRed    - 7,   # soft coral (lighter red variant)
    ]


input_Ck_path = {
    '010': '/home/galucia/Lithium4/femto/models/pHe3_square_well_010_GeV.root',
    '050': '/home/galucia/Lithium4/femto/models/pHe3_square_well_050_GeV.root',
    '1050': '/home/galucia/Lithium4/femto/models/pHe3_square_well_1050_GeV.root',
    #'1030': '/home/galucia/Lithium4/femto/models/CATS_cent10_30_converted.root',
    #'3050': '/home/galucia/Lithium4/femto/models/CATS_cent30_50_converted.root',
}
input_Sigma_Ck_path = {
    '010': '/home/galucia/Lithium4/femto/models/SigmaHe3_Coulomb_010_GeV.root',
    '050': '/home/galucia/Lithium4/femto/models/SigmaHe3_Coulomb_050_GeV.root',
    '1050': '/home/galucia/Lithium4/femto/models/SigmaHe3_Coulomb_1050_GeV.root',
    #'1030': '/home/galucia/Lithium4/femto/models/CATS_cent10_30_converted.root',
    #'3050': '/home/galucia/Lithium4/femto/models/CATS_cent30_50_converted.root',
}
input_Ck_path_variations = {
    '010': '/home/galucia/Lithium4/femto/models/pHe3_square_well_variations_ready.root',
}
input_Sigma_Ck_path_variations = {
    '010': '/home/galucia/Lithium4/femto/models/he3_Sigma_plus_Coulomb_variations_ready.root',
}
experimental_Ck_path = '/home/galucia/Lithium4/preparation/checks/correlation_hadronpid_pass1_pass4_refined_dca.root'
experimental_Ck_name = 'Correlation/Default/hCorrelation010'

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
                         tree, [kstar, weight_var], '', 'weight')
    
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

def produce_lambda_with_modified_values(sign:str, centrality:str, h_theoretical_Ck, h_theoretical_Sigma_Ck,
                                        outdir:TDirectory, modification_factor:float):
    '''
        Produce a lambda-corrected model by systematically changing the lambda value by x%, 
        where x is the modification_factor. This can be used to understand the sensitivity of the model to changes in lambda.
        The lambda is changed both in lambda + x% and lambda - x% to understand the effect in both directions.
    '''

    h_lambda_parameter = load_hist('/home/galucia/Lithium4/calibration/output/lambda_parameters_smaller_tolerance_pr.root',
                                        f'{sign}/hLambdaParameters')

    h_lambda_Sigma_parameter = load_hist('/home/galucia/Lithium4/calibration/output/lambda_parameters_smaller_tolerance_pr.root',
                                        f'{sign}/hLambdaSigmaParameters')
    
    h_lambda_parameter_higher_lambda = h_lambda_parameter.Clone(f'hLambdaParameter_HigherLambda')
    h_lambda_parameter_lower_lambda = h_lambda_parameter.Clone(f'hLambdaParameter_LowerLambda')

    h_lambda_corrected_Ck_higher_lambda = h_theoretical_Ck.Clone(f'hLambdaCorrectedCk_Higher')
    h_lambda_corrected_Ck_lower_lambda = h_theoretical_Ck.Clone(f'hLambdaCorrectedCk_Lower')

    h_lambda_Sigma_corrected_Ck_higher_lambda = h_theoretical_Ck.Clone(f'hLambdaSigmaCorrectedCk_Higher')
    h_lambda_Sigma_corrected_Ck_lower_lambda = h_theoretical_Ck.Clone(f'hLambdaSigmaCorrectedCk_Lower')

    for ibin in range(1, h_lambda_corrected_Ck_higher_lambda.GetNbinsX()+1):

        kstar = h_lambda_corrected_Ck_higher_lambda.GetBinCenter(ibin)
        lambda_param = h_lambda_parameter.GetBinContent(h_lambda_parameter.FindBin(kstar))
        lambda_Sigma_param = h_lambda_Sigma_parameter.GetBinContent(h_lambda_Sigma_parameter.FindBin(kstar))
        original_value = h_theoretical_Ck.GetBinContent(ibin)
        original_Sigma_value = h_theoretical_Sigma_Ck.GetBinContent(h_theoretical_Sigma_Ck.FindBin(kstar))
        
        modified_lambda_higher = lambda_param * (1 + modification_factor)
        modified_lambda_lower = lambda_param * (1 - modification_factor)
        
        corrected_value_higher = 1.0 + modified_lambda_higher * (original_value - 1.0)
        corrected_value_lower = 1.0 + modified_lambda_lower * (original_value - 1.0)

        corrected_Sigma_value_higher = modified_lambda_higher * original_value + lambda_Sigma_param * original_Sigma_value + (1 - modified_lambda_higher - lambda_Sigma_param) * 1.0
        corrected_Sigma_value_lower = modified_lambda_lower * original_value + lambda_Sigma_param * original_Sigma_value + (1 - modified_lambda_lower - lambda_Sigma_param) * 1.0

        h_lambda_parameter_higher_lambda.SetBinContent(h_lambda_parameter_higher_lambda.FindBin(kstar), modified_lambda_higher)
        h_lambda_parameter_lower_lambda.SetBinContent(h_lambda_parameter_lower_lambda.FindBin(kstar), modified_lambda_lower)
        
        h_lambda_corrected_Ck_higher_lambda.SetBinContent(ibin, corrected_value_higher)
        h_lambda_corrected_Ck_lower_lambda.SetBinContent(ibin, corrected_value_lower)

        h_lambda_Sigma_corrected_Ck_higher_lambda.SetBinContent(ibin, corrected_Sigma_value_higher)
        h_lambda_Sigma_corrected_Ck_lower_lambda.SetBinContent(ibin, corrected_Sigma_value_lower)

    return (h_lambda_parameter_higher_lambda, h_lambda_corrected_Ck_lower_lambda, h_lambda_corrected_Ck_higher_lambda, h_lambda_corrected_Ck_lower_lambda,
            h_lambda_Sigma_corrected_Ck_higher_lambda, h_lambda_Sigma_corrected_Ck_lower_lambda)

def apply_resolution_smearing(h_correlation_function, outdir:TDirectory):

    h_resolution = load_hist('/data/galucia/lithium_local/MC/AnalysisResults_LHC24i5.root',
                            'he3-hadron-femto/QA/hKstarRecVsKstarGen')
    h_mixed_event = load_hist('/home/galucia/Lithium4/preparation/checks/mixed_event_hadronpid_pass1_pass4_reject_multiples.root',
                            'QA/hKstar')
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
            weight = h_resolution_slice.GetBinContent(jbin) * mixed_weight if 0.01 < kstar_gen < 0.7 else 0.  
            # skip the region where the corrected correlation function is not defined
            correlation_value = h_correlation_function_matched.GetBinContent(jbin)
            
            smeared_value += correlation_value * weight
            total_weight += weight

        if total_weight > 0:
            smeared_value /= total_weight
            h_smeared_correlation_function.SetBinContent(ibin, smeared_value)


    return h_smeared_correlation_function


def produce_lambda_models(sign:str, centrality:str, outdir:TDirectory):
    
    h_lambda_parameter = load_hist('/home/galucia/Lithium4/calibration/output/lambda_parameters_smaller_tolerance_pr.root',
                                        f'{sign}/hLambdaParameters')
    h_lambda_Sigma_parameter = load_hist('/home/galucia/Lithium4/calibration/output/lambda_parameters_smaller_tolerance_pr.root',
                                        f'{sign}/hLambdaSigmaParameters')
    h_theoretical_Ck = load_hist(input_Ck_path[centrality], 'hHe3_p_Coul_CF')
    h_theoretical_Sigma_Ck = load_hist(input_Sigma_Ck_path[centrality], 'hhe3_Sigma_plus_CF')
    
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
    h_lambda_parameter_higher_lambda, h_lambda_parameter_lower_lambda, \
        h_lambda_corrected_Ck_higher_lambda, h_lambda_corrected_Ck_lower_lambda, \
        h_lambda_Sigma_corrected_Ck_higher_lambda, h_lambda_Sigma_corrected_Ck_lower_lambda = \
        produce_lambda_with_modified_values(sign, centrality, outdir=outdir,
                                            h_theoretical_Ck=h_theoretical_Ck, h_theoretical_Sigma_Ck=h_theoretical_Sigma_Ck,
                                             modification_factor=modification_factor)
    
    h_correlation_reference = load_hist(experimental_Ck_path, experimental_Ck_name)
    h_lambda_corrected_Ck_matched = match_bin_width_correlation_function(h_correlation_reference, h_lambda_corrected_Ck)
    h_lambda_Sigma_corrected_Ck_matched = match_bin_width_correlation_function(h_correlation_reference, h_lambda_Sigma_corrected_Ck)
    
    h_lambda_Sigma_smeared_Ck = apply_resolution_smearing(h_lambda_Sigma_corrected_Ck, outdir)
    
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

    h_lambda_corrected_Ck_matched.Write('hLambdaCorrectedCk_Matched')
    h_lambda_Sigma_corrected_Ck_matched.Write('hLambdaSigmaCorrectedCk_Matched')

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

    legend = TLegend(0.6, 0.2, 0.88, 0.4)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)

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

    legend = TLegend(0.6, 0.2, 0.88, 0.4)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)

    legend.AddEntry(h_lambda_Sigma_corrected_Ck_higher_lambda, f'#lambda_{{nominal}} + {modification_factor*100.:.0f}%', 'l')
    legend.AddEntry(h_lambda_Sigma_corrected_Ck, '#lambda_{nominal}', 'l')
    legend.AddEntry(h_lambda_Sigma_corrected_Ck_lower_lambda, f'#lambda_{{nominal}} - {modification_factor*100.:.0f}%', 'l')
    legend.Draw()

    canvas.Write()

    del canvas

def produce_lambda_models_with_variations(sign:str, centrality:str, outdir:TDirectory):

    variation_dict = {'nominal': '7.12',
                      'upper': '7.16',
                      'lower': '7.06'}
    
    h_lambda_parameter = load_hist('/home/galucia/Lithium4/calibration/output/lambda_parameters_smaller_tolerance_pr.root',
                                        f'{sign}/hLambdaParameters')
    h_lambda_Sigma_parameter = load_hist('/home/galucia/Lithium4/calibration/output/lambda_parameters_smaller_tolerance_pr.root',
                                        f'{sign}/hLambdaSigmaParameters')
    
    for variation_name, variation_value in variation_dict.items():  

        variation_dir = outdir.mkdir(variation_name)
        input_variation_dir = f'r={variation_value}_fm'

        h_theoretical_Ck = load_hist(input_Ck_path_variations[centrality], f'{input_variation_dir}/hHe3_p_Coul_CF')
        h_theoretical_Sigma_Ck = load_hist(input_Sigma_Ck_path_variations[centrality], f'{input_variation_dir}/hhe3_Sigma_plus_CF')
    
        h_lambda_corrected_Ck = h_theoretical_Ck.Clone(f'hLambdaCorrectedCk')
        h_lambda_Sigma_corrected_Ck = h_theoretical_Ck.Clone(f'hLambdaSigmaCorrectedCk')

        for ibin in range(1, h_lambda_corrected_Ck.GetNbinsX()+1):
            
            kstar = h_lambda_corrected_Ck.GetBinCenter(ibin)
            lambda_param = h_lambda_parameter.GetBinContent(h_lambda_parameter.FindBin(kstar))
            original_value = h_theoretical_Ck.GetBinContent(ibin)

            lambda_Sigma_param = h_lambda_Sigma_parameter.GetBinContent(h_lambda_Sigma_parameter.FindBin(kstar))
            Sigma_value = h_theoretical_Sigma_Ck.GetBinContent(h_theoretical_Sigma_Ck.FindBin(kstar))

            lambda_Sigma_corrected_value = lambda_param * original_value + lambda_Sigma_param * Sigma_value + (1 - lambda_param - lambda_Sigma_param) * 1.0
            h_lambda_Sigma_corrected_Ck.SetBinContent(ibin, lambda_Sigma_corrected_value)

        modification_factor = 0.1  # 10% change in lambda
        __, __, __, __, h_lambda_Sigma_corrected_Ck_higher_lambda, h_lambda_Sigma_corrected_Ck_lower_lambda = \
            produce_lambda_with_modified_values(sign, centrality, outdir=variation_dir,
                                                h_theoretical_Ck=h_theoretical_Ck, h_theoretical_Sigma_Ck=h_theoretical_Sigma_Ck,
                                                modification_factor=modification_factor)
    
        h_lambda_Sigma_smeared_Ck = apply_resolution_smearing(h_lambda_Sigma_corrected_Ck, outdir)
        h_lambda_Sigma_smeared_Ck_higher_lambda = apply_resolution_smearing(h_lambda_Sigma_corrected_Ck_higher_lambda, outdir)
        h_lambda_Sigma_smeared_Ck_lower_lambda = apply_resolution_smearing(h_lambda_Sigma_corrected_Ck_lower_lambda, outdir)
    
        variation_dir.cd()
        for hist in [h_lambda_Sigma_smeared_Ck, h_lambda_Sigma_smeared_Ck_higher_lambda, h_lambda_Sigma_smeared_Ck_lower_lambda]:
            set_root_object(hist, title='; #it{k}* (MeV/c); C(#it{k}*)', line_width=2) 
            hist.Write(hist.GetName())

if __name__ == '__main__':

    outfile = TFile.Open('models/lambda_models.root', 'recreate')

    for sign in ['Both', 'Matter', 'Antimatter']:
        print(f"\n{'='*60}")
        print(f"Processing {sign}")

        outdir_sign = outfile.mkdir(sign)

        for centrality in ['050', '010', '1050']: #, '1030', '3050']:
            print(f"\n{'-'*40}")
            print(f"Processing centrality {centrality}")

            outdir = outdir_sign.mkdir(f'{centrality}')
            produce_lambda_models(sign, centrality, outdir)

            if sign == 'Both':
                continue

            if centrality == '010': 
                produce_lambda_models_with_variations(sign, centrality, outdir)
            
    outfile.Close()