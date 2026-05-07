import numpy as np
import pandas as pd
from ROOT import TFile, RooRealVar, TDirectory, TGraphErrors, TF1, RooDataHist, TCanvas

from torchic import Dataset, AxisSpec
from torchic.roopdf.roopdf_utils import init_roopdf

def parametrise_resolution(h2_res, xmin, xmax, ymin, ymax, output_file: TDirectory):

    fits_dir = output_file.mkdir('fits')
    fit_results = pd.DataFrame()

    x = RooRealVar('x', '#it{p}_{T} (GeV/#it{c})', ymin, ymax)
    #roofitter = Roofitter(x, ['exp_mod_gaus', 'crystal_ball'])
    
    #pars['rlife'].setRange(-1., 0.)
    
    low_bin = h2_res.GetXaxis().FindBin(xmin)
    high_bin = h2_res.GetXaxis().FindBin(xmax)

    for ibin in range(low_bin, high_bin+1):
        
        gausexp, pars = init_roopdf('gaus_exp', x, name='exp_mod_gaus')
        pars['mean'].setRange(ymin, ymax)
        pars['sigma'].setRange(0.0001, 0.5)
        pars['mean'].setVal(0.05)
        pars['sigma'].setVal(0.1)
        pars['mean'].setConstant(False)
        pars['sigma'].setConstant(False)
    
        h1 = h2_res.ProjectionY(f'h1_{ibin}', ibin, ibin)
        data = RooDataHist(f'data_{ibin}', f'data_{ibin}', [x], Import=h1)
        gausexp.fitTo(data, Range=(ymin, ymax))
        
        frame = x.frame()
        data.plotOn(frame)
        gausexp.plotOn(frame)
        gausexp.paramOn(frame, Format='NEU', Layout=(0.6, 0.9, 0.9))
        chi2 = frame.chiSquare()
        print(f'Bin {ibin}, chi2/ndf = {chi2:.2f}')
        canvas = TCanvas(f'c_{ibin}', f'c_{ibin}', 800, 600)
        frame.Draw()
        fits_dir.cd()
        canvas.Write()
        
        bin_fit_results = pd.DataFrame.from_dict({
            'exp_mod_gaus_0_mean': [pars['mean'].getVal()],
            'exp_mod_gaus_0_sigma': [pars['sigma'].getVal()],
            #'integral': [gausexp.expectedEvents(Range=(ymin, ymax))],
        })
        bin_fit_results['bin_center'] = h2_res.GetXaxis().GetBinCenter(ibin)
        bin_fit_results['unnorm_integral'] = h1.Integral() #* bin_fit_results['integral'] #* fractions[0].getVal()

        fit_results = pd.concat([fit_results, bin_fit_results], ignore_index=True)
    
    bin_error = (fit_results['bin_center'][1] - fit_results['bin_center'].values[0])/2.
    fit_results['bin_error'] = bin_error

    fit_results['mean_error'] = fit_results['exp_mod_gaus_0_sigma'] / fit_results['unnorm_integral']
    graph_mean = TGraphErrors(len(fit_results), np.array(fit_results['bin_center'], dtype=np.float32), np.array(fit_results['exp_mod_gaus_0_mean'], dtype=np.float32), np.array(fit_results['bin_error'], dtype=np.float32), np.array(fit_results['mean_error'], dtype=np.float32))
    graph_mean.SetName('gResolutionMean')
    graph_mean.SetTitle(';#it{p}_{T}^{rec} (GeV/#it{c});(#it{p}_{T}^{rec} - #it{p}_{T}^{gen}) / #it{p}_{T}^{rec}')

    xmin_fit = 1.6
    xmax_fit = 2.4
    # fit_res = TF1('fit_res', '[0] / (1 + exp((x - [1])/[2]) )', xmin_fit, xmax_fit)
    # fit_res.SetParameter(0, 0.1)
    # fit_res.SetParameter(1, 1.9)
    # fit_res.SetParameter(2, 0.1)
    # fit_res = TF1('fit_res', '[0] + [1] * x + [2] * x^2', xmin_fit, xmax_fit)
    fit_res = TF1('fit_res', '[0] + [1] * x', xmin_fit, xmax_fit)
    fit_res.SetParameter(0, 0.1)
    fit_res.SetParameter(1, -1.)
    # fit_res.SetParameter(2, 0.1)
    graph_mean.Fit(fit_res, 'RMS+')

    # fit_res_results = [fit_res.GetParameter(i) for i in range(3)]
    fit_res_results = [fit_res.GetParameter(i) for i in range(2)]

    output_file.cd()
    graph_mean.Write()
    fit_res.Write()

    return fit_res_results

def main(input_file_path, output_file_path):

    dataset_mc = Dataset.from_root(input_file_path,
                                tree_name='O2he3hadtablemc',
                                folder_name='DF*')
    dataset_rec = Dataset.from_root(input_file_path,
                                tree_name='O2he3hadtable',
                                folder_name='DF*')
    dataset = dataset_rec.concat(dataset_mc, axis=1)
    print(f'{dataset.columns=}')

    # select he3
    #dataset.query('abs(fPDGcode) == 1000020030', inplace=True)
    #readPidTracking = lambda x: (x >> 12) & 0x1F
    #readPidTracking = np.vectorize(readPidTracking)
    #dataset['fPidTracking'] = readPidTracking(dataset['fFlags'])

    dataset.add_subset('H3', dataset['fPIDtrkHe3'] == 6)
    dataset.add_subset('He4', dataset['fPIDtrkHe3'] == 8)
    #dataset.loc[dataset['fPDGcode'] < 0, 'fPt'] = - dataset['fPt']

    #dataset['fPt'] = 2 * abs(dataset['fPt'])
    dataset['fPtHe3'] = abs(dataset['fPtHe3'])
    dataset['fPtDiff'] = dataset['fPtHe3'] - dataset['fPtMCHe3']
    dataset['fPtRes'] = dataset['fPtDiff'] / dataset['fPtHe3']
    
    dataset['fPHe3'] = dataset['fPtHe3'] * np.cosh(dataset['fEtaHe3'])
    dataset['fPMCHe3'] = dataset['fPtMCHe3'] * np.cosh(dataset['fEtaHe3'])
    dataset['fPDiff'] = dataset['fPHe3'] - dataset['fPMCHe3']
    dataset['fPRes'] = dataset['fPDiff'] / dataset['fPHe3']

    
    axis_spec_pt = AxisSpec(100, 0., 5., 'fPtHe3', ';#it{p}_{T}^{rec} (GeV/#it{c})')
    axis_spec_gpt = AxisSpec(100, 0., 5., 'fPtMCHe3', ';#it{p}_{T}^{gen} (GeV/#it{c})')

    h_pt = dataset.build_th1('fPtHe3', axis_spec_pt)
    h_gpt = dataset.build_th1('fPtMCHe3', axis_spec_gpt)
    h2_res = dataset.build_th2('fPtHe3', 'fPtRes', axis_spec_pt, AxisSpec(100, -0.5, 0.5, 'fPtRes', ';#it{p}_{T}^{rec} (GeV/#it{c});(#it{p}_{T}^{rec} - #it{p}_{T}^{gen}) / #it{p}_{T}^{rec}'))
    h2_resh3 = dataset.build_th2('fPtHe3', 'fPtRes', axis_spec_pt, AxisSpec(100, -0.5, 0.5, 'fPtResH3', ';#it{p}_{T}^{rec} (GeV/#it{c});(#it{p}_{T}^{rec} - #it{p}_{T}^{gen}) / #it{p}_{T}^{rec}'), subset='H3')
    h2_reshe4 = dataset.build_th2('fPtHe3', 'fPtRes', axis_spec_pt, AxisSpec(100, -0.5, 0.5, 'fPtResHe4', ';#it{p}_{T}^{rec} (GeV/#it{c});(#it{p}_{T}^{rec} - #it{p}_{T}^{gen}) / #it{p}_{T}^{rec}'), subset='He4')
    
    axis_spec_p = AxisSpec(100, 0., 10., 'fPHe3', ';#it{p}^{rec} (GeV/#it{c})')
    axis_spec_gp = AxisSpec(100, 0., 10., 'fPMCHe3', ';#it{p}^{gen} (GeV/#it{c})')

    h_p = dataset.build_th1('fPHe3', axis_spec_p)
    h_gp = dataset.build_th1('fPMCHe3', axis_spec_gp)
    h2_res_p = dataset.build_th2('fPHe3', 'fPRes', axis_spec_p, AxisSpec(100, -0.5, 0.5, 'fPRes', ';#it{p}^{rec} (GeV/#it{c});(#it{p}^{rec} - #it{p}^{gen}) / #it{p}^{rec}'))
    h2_res_ph3 = dataset.build_th2('fPHe3', 'fPRes', axis_spec_p, AxisSpec(100, -0.5, 0.5, 'fPResH3', ';#it{p}^{rec} (GeV/#it{c});(#it{p}^{rec} - #it{p}^{gen}) / #it{p}^{rec}'), subset='H3')

    output_file = TFile.Open(output_file_path, 'recreate')
    output_dir_pt = output_file.mkdir('pt_resolution')
    output_dir_p = output_file.mkdir('p_resolution')
    
    output_dir_pt.cd()
    h_pt.Write()
    h_gpt.Write()
    h2_res.Write()
    h2_resh3.Write()
    h2_reshe4.Write()
    
    output_dir_p.cd()
    h_p.Write()
    h_gp.Write()
    h2_res_p.Write()
    h2_res_ph3.Write()

    fit_res_results = parametrise_resolution(h2_resh3, 1.55, 2.2, -0.5, 0.5, output_dir_pt)

    # pol1 correction
    dataset.loc[dataset['fPIDtrkHe3'] == 6, 'fPtHe3'] = dataset['fPtHe3'] - dataset['fPtHe3']*(fit_res_results[0] + fit_res_results[1] * dataset['fPtHe3'])
    # pol2 correction
    # dataset.loc[dataset['fPidTracking'] == 6, 'fPt'] = dataset['fPt'] - dataset['fPt']*(fit_res_results[0] + fit_res_results[1] * dataset['fPt'] + fit_res_results[2] * dataset['fPt']**2)
    
    dataset['fPtDiff'] = dataset['fPtHe3'] - dataset['fPtMCHe3']
    dataset['fPtRes'] = dataset['fPtDiff'] / dataset['fPtHe3']
    h2_res_corrected = dataset.build_th2('fPtHe3', 'fPtRes', axis_spec_pt, AxisSpec(100, -0.5, 0.5, 'fPtResCorrected', ';#it{p}_{T}^{rec} (GeV/#it{c});(#it{p}_{T}^{rec} - #it{p}_{T}^{gen}) / #it{p}_{T}^{rec}'))
    h2_resh3_corrected = dataset.build_th2('fPtHe3', 'fPtRes', axis_spec_pt, AxisSpec(100, -0.5, 0.5, 'fPtResH3Corrected', ';#it{p}_{T}^{rec} (GeV/#it{c});(#it{p}_{T}^{rec} - #it{p}_{T}^{gen}) / #it{p}_{T}^{rec}'), subset='H3')
    
    output_dir_pt.cd()
    h2_res_corrected.Write()
    h2_resh3_corrected.Write()
    
    fit_res_results_p = parametrise_resolution(h2_res_ph3, 1.55, 2.2, -0.5, 0.5, output_dir_p)

    dataset.loc[dataset['fPIDtrkHe3'] == 6, 'fPHe3'] = dataset['fPHe3'] - dataset['fPHe3'] * (fit_res_results_p[0] + fit_res_results_p[1] * dataset['fPHe3'])

    dataset['fPDiff'] = dataset['fPHe3'] - dataset['fPMCHe3']
    dataset['fPRes'] = dataset['fPDiff'] / dataset['fPHe3']
    h2_res_p_corrected = dataset.build_th2('fPHe3', 'fPRes', axis_spec_p, AxisSpec(100, -0.5, 0.5, 'fPResCorrected', ';#it{p}^{rec} (GeV/#it{c});(#it{p}^{rec} - #it{p}^{gen}) / #it{p}^{rec}'))
    h2_res_ph3_corrected = dataset.build_th2('fPHe3', 'fPRes', axis_spec_p, AxisSpec(100, -0.5, 0.5, 'fPResH3Corrected', ';#it{p}^{rec} (GeV/#it{c});(#it{p}^{rec} - #it{p}^{gen}) / #it{p}^{rec}'), subset='H3')

    output_dir_p.cd()
    h2_res_p_corrected.Write()
    h2_res_ph3_corrected.Write()

if __name__ == '__main__':

    #input_file_path = '/data/galucia/lithium_local/MC/LHC26c6.root'
    #output_file_path = 'output/he3_pid_trk_correction.root'
    
    input_file_path = '/home/galucia/Lithium4/task/he3HadronFemto/AO2D_lit_mc.root'
    output_file_path = 'output/he3_pid_trk_correction_LHC25g11.root'
    main(input_file_path, output_file_path)
    