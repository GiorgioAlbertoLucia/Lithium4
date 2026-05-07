import numpy as np

from ROOT import TFile, TH1F


from torchic.core.histogram import load_hist, build_efficiency

def histogram_division(h_numerator:TH1F, h_denominator:TH1F):
    '''
        The numerator dictates the final binning
    '''

    h_reference_numerator = h_numerator.Clone('reference')

    for ibin in range(1, h_numerator.GetNbinsX()+1):
        
        bin_center = h_numerator.GetBinCenter(ibin)
        print(f'\n{bin_center=}')
        
        nearest_bin_denominator = h_denominator.FindBin(bin_center)
        print(f'{nearest_bin_denominator=}, {h_denominator.GetBinCenter(nearest_bin_denominator)=}')
        center_difference = bin_center - h_denominator.GetBinCenter(nearest_bin_denominator)
        print(f'{center_difference=}')
        second_nearest_bin = nearest_bin_denominator + 1 if center_difference > 0 else nearest_bin_denominator - 1
        print(f'{second_nearest_bin=}, {h_denominator.GetBinCenter(second_nearest_bin)=}')
        weight = np.abs(center_difference) / np.abs(h_denominator.GetBinCenter(nearest_bin_denominator) - h_denominator.GetBinCenter(second_nearest_bin))
        denominator = (1 - weight) * h_denominator.GetBinContent(nearest_bin_denominator) + weight * h_denominator.GetBinContent(second_nearest_bin)
        print(f'{weight=}, {denominator=}')

        numerator = h_reference_numerator.GetBinContent(ibin)
        value = numerator / denominator
        h_numerator.SetBinContent(ibin, value)

def histogram_multiplication(h_numerator:TH1F, h_denominator:TH1F):
    '''
        The numerator dictates the final binning
    '''

    h_reference_numerator = h_numerator.Clone('reference')

    for ibin in range(1, h_numerator.GetNbinsX()+1):
        
        bin_center = h_numerator.GetBinCenter(ibin)
        print(f'\n{bin_center=}')
        
        nearest_bin_denominator = h_denominator.FindBin(bin_center)
        print(f'{nearest_bin_denominator=}, {h_denominator.GetBinCenter(nearest_bin_denominator)=}')
        center_difference = bin_center - h_denominator.GetBinCenter(nearest_bin_denominator)
        print(f'{center_difference=}')
        second_nearest_bin = nearest_bin_denominator + 1 if center_difference > 0 else nearest_bin_denominator - 1
        print(f'{second_nearest_bin=}, {h_denominator.GetBinCenter(second_nearest_bin)=}')
        weight = center_difference / np.abs(h_denominator.GetBinCenter(nearest_bin_denominator) - h_denominator.GetBinCenter(second_nearest_bin))
        denominator = (1 - weight) * h_denominator.GetBinContent(nearest_bin_denominator) + weight * h_denominator.GetBinContent(second_nearest_bin)
        print(f'{weight=}, {denominator=}')

        numerator = h_reference_numerator.GetBinContent(ibin)
        value = numerator * denominator
        h_numerator.SetBinContent(ibin, value)

def interpolate_hist_two_bins(hist:TH1F, x):
    bin_c = hist.FindBin(x)

    # protect boundaries
    if bin_c <= 1:
        return hist.GetBinContent(1), hist.GetBinError(1)
    if bin_c >= hist.GetNbinsX():
        return hist.GetBinContent(hist.GetNbinsX()), hist.GetBinError(hist.GetNbinsX())

    x1 = hist.GetBinCenter(bin_c)
    x0 = hist.GetBinCenter(bin_c - 1)

    y1 = hist.GetBinContent(bin_c)
    y0 = hist.GetBinContent(bin_c - 1)

    e1 = hist.GetBinError(bin_c)
    e0 = hist.GetBinError(bin_c - 1)

    # linear weights
    w1 = (x - x0) / (x1 - x0)
    w0 = 1.0 - w1

    val = w0 * y0 + w1 * y1
    err = np.sqrt((w0 * e0)**2 + (w1 * e1)**2)

    return val, err


if __name__ == '__main__':

    outfile = TFile.Open('output/he3_from_hypertriton.root', 'recreate')

    BR = 0.25
    BR_error = 0.023

    h_pt_he3_gen = load_hist('/home/galucia/Lithium4/phase_space_studies/output/single_track_efficiency.root',
                             'He/hPtGenPrimaryHe')
    h_pt_he3_rec = load_hist('/home/galucia/Lithium4/phase_space_studies/output/single_track_efficiency.root',
                             'He/hPtRecPrimaryHe')
    
    h_efficiency_he3 = build_efficiency(h_pt_he3_gen, h_pt_he3_rec, name='h_efficiency_he3', ytitle='Efficiency')

    h_pt_he3_from_h3l_gen = load_hist('/home/galucia/Lithium4/phase_space_studies/output/single_track_efficiency.root',
                             'He/hPtGenFromHypertritonHe')
    h_pt_he3_from_h3l_rec = load_hist('/home/galucia/Lithium4/phase_space_studies/output/single_track_efficiency.root',
                             'He/hPtRecFromHypertritonHe')
    
    h_efficiency_he3_from_h3l = build_efficiency(h_pt_he3_from_h3l_gen, h_pt_he3_from_h3l_rec,
                                                 name='h_efficiency_he3_from_h3l', ytitle='Efficiency')
    
    # only for visualisation, will not be used for computation
    h_pt_he3_from_h3l_gen_map = load_hist('/home/galucia/Lithium4/phase_space_studies/input/He3FromHypertritonMap.root',
                             'nuclei-from-hypertriton-map/registryMC/he3SecPtGen_from_hypertriton')
    h_pt_h3l_into_he3_gen = load_hist('/home/galucia/Lithium4/phase_space_studies/input/He3FromHypertritonMap.root',
                             'nuclei-from-hypertriton-map/registryMC/hypertritonPtGen')
    h_correction_factor = h_pt_he3_from_h3l_gen_map.Clone('h_correction_factor')
    h_correction_factor.Divide(h_pt_h3l_into_he3_gen)

    spectrum_file = TFile.Open('output/hypertriton_spectra.root')
    f_spectrum_he3 = spectrum_file.Get('fBlastWaveHe3')
    f_spectrum_h3l = spectrum_file.Get('fBlastWaveHypertriton')
    fit_result_he3 = spectrum_file.Get('resultBlastWaveHe3')
    fit_result_h3l = spectrum_file.Get('resultBlastWaveHypertriton')

    pt_min, pt_max, pt_step = -8., 8., 0.5
    nbins = int((pt_max - pt_min) / pt_step)
    pt_array = np.arange(pt_min, pt_max, pt_step)

    ratio_spectrum = TH1F('h_ratio', '^{3}_{#Lambda}H #times #Gamma(^{3}_{#Lambda}H #rightarrow ^{3}He) / ^{3}He; #it{p}_{T} (GeV/#it{c}); ^{3}_{#Lambda}H #times #Gamma(^{3}_{#Lambda}H #rightarrow ^{3}He) / ^{3}He',
                          nbins, pt_min, pt_max)
    fraction_spectrum = TH1F('h_fraction', ' d#it{N}/d#it{p}_{T}(^{3}He #leftarrow ^{3}_{#Lambda}H) / d#it{N}/d#it{p}_{T}(^{3}He)_{primary}; #it{p}_{T} (GeV/#it{c}); d#it{N}/d#it{p}_{T}(^{3}He #leftarrow ^{3}_{#Lambda}H) / d#it{N}/d#it{p}_{T}(^{3}He)_{primary}',
                          nbins, pt_min, pt_max)
    hypertriton_to_he3_map = TH1F('h_map', 'd#it{N}/d#it{p}_{T} (^{3}He #leftarrow ^{3}_{#Lambda}H) / d#it{N}/d#it{p}_{T} (^{3}_{#Lambda}H #rightarrow ^{3}He); #it{p}_{T} (GeV/#it{c}); d#it{N}/d#it{p}_{T} (^{3}He #leftarrow ^{3}_{#Lambda}H) / d#it{N}/d#it{p}_{T} (^{3}_{#Lambda}H #rightarrow ^{3}He)',
                          nbins, pt_min, pt_max)

    # ratio of the spectra
    for ibin in range(1, ratio_spectrum.GetNbinsX()+1):

        pt = ratio_spectrum.GetBinCenter(ibin)
        pt_low = ratio_spectrum.GetBinLowEdge(ibin)
        pt_high = ratio_spectrum.GetBinLowEdge(ibin+1)
        #pt = 0.5 * (pt_array[ibin] + pt_array[ibin-1])

        # swap and use absolute for negative pt
        if np.abs(pt_low) > np.abs(pt_high):
            pt_low, pt_high = np.abs(pt_high), np.abs(pt_low)

        hypertriton_yield = f_spectrum_h3l.Integral(pt_low, pt_high-0.001)
        hypertriton_yield_error = f_spectrum_h3l.IntegralError(pt_low, pt_high-0.001,
                                                               fit_result_h3l.GetParams(), fit_result_h3l.GetCovarianceMatrix().GetMatrixArray())
        he3_yield = f_spectrum_he3.Integral(pt_low, pt_high-0.001)
        he3_yield_error = f_spectrum_he3.IntegralError(pt_low, pt_high-0.001,
                                                               fit_result_he3.GetParams(), fit_result_he3.GetCovarianceMatrix().GetMatrixArray())
        
        value = hypertriton_yield * BR / he3_yield if he3_yield > 1e-7 else 0.
        error = np.sqrt( (hypertriton_yield_error * BR /he3_yield)**2 +
                        (hypertriton_yield * BR_error /he3_yield)**2 +
                        (hypertriton_yield * BR * he3_yield_error/(he3_yield*he3_yield))**2
                        ) if he3_yield > 1e-7 else 0.

        ratio_spectrum.SetBinContent(ibin, value)
        ratio_spectrum.SetBinError(ibin, error)

    # hypertriton to he3 map
    for ibin in range(1, hypertriton_to_he3_map.GetNbinsX()+1):

        pt = ratio_spectrum.GetBinCenter(ibin)
        pt_low = ratio_spectrum.GetBinLowEdge(ibin)
        pt_high = ratio_spectrum.GetBinLowEdge(ibin+1)

        N_he3_from_h3l, N_he3_from_h3l_error = interpolate_hist_two_bins(h_pt_he3_from_h3l_gen_map, pt)
        N_h3l_into_he3, N_h3l_into_he3_error = interpolate_hist_two_bins(h_pt_h3l_into_he3_gen, pt)
        
        if pt < 0.:
            print(f'{pt=}, {pt_low=}, {pt_high=}')
            print(f'{N_he3_from_h3l=}, {N_h3l_into_he3=}')

        value = N_he3_from_h3l / N_h3l_into_he3 if N_h3l_into_he3 > 1e-7 else 0.
        error = np.sqrt( (N_he3_from_h3l_error / N_h3l_into_he3)**2 + 
                         (N_he3_from_h3l * N_h3l_into_he3_error / (N_h3l_into_he3*N_h3l_into_he3))**2 
                        ) if N_h3l_into_he3 > 1e-7 else 0.

        hypertriton_to_he3_map.SetBinContent(ibin, value)
        hypertriton_to_he3_map.SetBinError(ibin, error)
        
        
    for ibin in range(1, fraction_spectrum.GetNbinsX()+1):

        pt = fraction_spectrum.GetBinCenter(ibin)
        #pt = 0.5 * (pt_array[ibin] + pt_array[ibin-1])

        ratio_spectrum_value = ratio_spectrum.GetBinContent(ibin)
        ratio_spectrum_error = ratio_spectrum.GetBinError(ibin)

        h3l_into_he3_generated = h_pt_h3l_into_he3_gen.GetBinContent(h_pt_h3l_into_he3_gen.FindBin(pt))
        h3l_into_he3_error = h_pt_h3l_into_he3_gen.GetBinError(h_pt_h3l_into_he3_gen.FindBin(pt))
        he3_from_h3l_generated = h_pt_he3_from_h3l_gen.GetBinContent(h_pt_he3_from_h3l_gen.FindBin(pt))
        he3_from_h3l_error = h_pt_he3_from_h3l_gen.GetBinError(h_pt_he3_from_h3l_gen.FindBin(pt))

        efficiency = h_efficiency_he3.GetBinContent(h_efficiency_he3.FindBin(pt))
        efficiency_error = h_efficiency_he3.GetBinError(h_efficiency_he3.FindBin(pt))

        efficiency_from_h3l = h_efficiency_he3_from_h3l.GetBinContent(h_efficiency_he3_from_h3l.FindBin(pt))
        efficiency_from_h3l_error = h_efficiency_he3_from_h3l.GetBinError(h_efficiency_he3_from_h3l.FindBin(pt))

        map_value = hypertriton_to_he3_map.GetBinContent(hypertriton_to_he3_map.FindBin(pt))
        map_error = hypertriton_to_he3_map.GetBinError(hypertriton_to_he3_map.FindBin(pt))

        value = ratio_spectrum_value * efficiency_from_h3l / efficiency if efficiency > 1e-7 else 0.
        error = np.sqrt( (ratio_spectrum_error * map_value * efficiency_from_h3l / efficiency)**2 +
                        (ratio_spectrum_value * map_error * efficiency_from_h3l / efficiency)**2 +
                        (ratio_spectrum_value * map_value * efficiency_from_h3l_error / efficiency)**2 +
                        (ratio_spectrum_value * map_value * efficiency_from_h3l * efficiency_error / (efficiency * efficiency))**2
                        ) if efficiency > 1e-7 else 0.

        fraction_spectrum.SetBinContent(ibin, value)
        fraction_spectrum.SetBinError(ibin, error)

    
    outfile.cd()
    
    h_efficiency_he3.Write()
    h_efficiency_he3_from_h3l.Write()

    h_correction_factor.Write()

    hypertriton_to_he3_map.Write()

    h_pt_he3_gen.Write()
    h_pt_he3_from_h3l_gen.Write()
    h_pt_h3l_into_he3_gen.Write()
    h_pt_he3_from_h3l_rec.Write()
    
    ratio_spectrum.Write()
    fraction_spectrum.Write()
    h_efficiency_he3_from_h3l.Write()
