'''
    Code to run calibration of ITS and tof parametrisations
'''

import numpy as np
import pandas as pd
from ROOT import TFile, TCanvas, TF1, TH2F, TDirectory
from ROOT import RooRealVar, RooCrystalBall, RooAddPdf, RooGaussian, RooAbsPdf, RooExponential, RooGenericPdf

from particle import Particle

from torchic import Dataset, AxisSpec
from torchic.roopdf import RooGausExp, RooGausDExp

from torchic.core.graph import create_graph

import sys
sys.path.append('..')
from utils.utils import initialize_means_and_covariances, calibration_fit_slice

PT_MIN = 0.
PT_MAX = 4.0
PT_BKG_FIT = {'matter':  0., #2.29,
              'antimatter': 0., #2.29 # 0.8
              } 

def init_signal_roofit(nsigma_tof: RooRealVar, function: str = 'crystalball'):

    if function == 'crystalball':
        signal_pars = {
            'mean': RooRealVar('mean', '#mu_{sig}', 1., 0.7, 1.2, ''),
            'sigma': RooRealVar('sigma', '#sigma_{sig}', 0.1, 0.0001, 0.3, ''),
            'aL': RooRealVar('aL', '#alpha_{L, sig}', 0.7, 30.),
            'nL': RooRealVar('nL', '#it{n}_{L, sig}', 0.3, 30.),
            'aR': RooRealVar('aR', '#alpha_{R, sig}', 0.7, 30.),
            'nR': RooRealVar('nR', '#it{n}_{R, sig}', 0.3, 30.),
        }
        signal = RooCrystalBall('signal', 'signal', nsigma_tof, signal_pars['mean'], signal_pars['sigma'],
                                signal_pars['aL'], signal_pars['nL'], doubleSided=True) #
                                #signal_pars['aR'], signal_pars['nR'])

        return signal, signal_pars
    
    elif function == 'gausexp':
        signal_pars = {
            'mean': RooRealVar('mean', '#mu_{sig}', 1., 0.7, 1.2, ''),
            'sigma': RooRealVar('sigma', '#sigma_{sig}', 0.1, 0.0001, 0.3, ''),
            'rlife': RooRealVar('rlife', '#tau_{sig}', 0., 10.),
        }
        signal = RooGausExp('signal', 'signal', nsigma_tof, *signal_pars.values())
        return signal, signal_pars
    
    elif function == 'gausdexp':
        signal_pars = {
            'mean': RooRealVar('mean', '#mu_{sig}', 0., -1., 1., ''),
            'sigma': RooRealVar('sigma', '#sigma_{sig}', 1, 0.3, 2., ''),
            'rlife0': RooRealVar('rlife0', '#tau_{0, sig}', -10., 0.),
            'rlife1': RooRealVar('rlife1', '#tau_{1, sig}', 0., 10.),
        }
        signal = RooGausDExp('signal', 'signal', nsigma_tof, *signal_pars.values())
        return signal, signal_pars
    
    elif function == 'gaus':
        signal_pars = {
            'mean': RooRealVar('mean', 'mean', 1., 0.7, 1.2, ''),
            'sigma': RooRealVar('sigma', 'sigma', 0.1, 0.0001, 0.3, ''),
        }
        signal = RooGaussian('signal', 'signal', nsigma_tof, *signal_pars.values())
        return signal, signal_pars
    
    else:
        raise ValueError(f'Unknown function: {function}. Supported functions are "crystalball" and "gausexp".')

def init_background_roofit(nsigma_tof: RooRealVar, function: str = 'gaus'):

    if function == 'gausexp':
        bkg_pars = {
            'mean': RooRealVar('bkg_mean', '#mu_{bkg}', 1., 0.7, 1.2, ''),
            'sigma': RooRealVar('bkg_sigma', '#sigma_{bkg}', 0.1, 0.0001, 0.3, ''),
            'rlife': RooRealVar('bkg_rlife', '#tau_{bkg}', 0., 10.),
        }
        bkg_pdf = RooGausExp('bkg', 'bkg', nsigma_tof, *bkg_pars.values())
        return bkg_pdf, bkg_pars

    elif function == 'gaus':
        bkg_pars = {
            'mean': RooRealVar('bkg_mean', '#mu_{bkg}', 1., 0.7, 1.2, ''),
            'sigma': RooRealVar('bkg_sigma', '#sigma_{bkg}', 0.1, 0.0001, 0.3, ''),
        }
        bkg_pdf = RooGaussian('bkg', 'bkg', nsigma_tof, *bkg_pars.values())
        return bkg_pdf, bkg_pars
    
    elif function == 'exp':
        bkg_pars = {
            'alpha': RooRealVar('bkg_alpha', '#alpha_{bkg}', -0.1, -10, 10),
        }
        bkg_pdf = RooExponential('bkg', 'bkg', nsigma_tof, *bkg_pars.values())
        return bkg_pdf, bkg_pars
    
    elif function == 'exp':
        pars = {
            'alpha': RooRealVar('alpha', 'alpha', 1., 0.01, 100.),
            'offset': RooRealVar('offset', 'offset', -5., -20., 0.)
        }
        pdf = RooGenericPdf(f'bkg', 'bkg', 
                           f'exp(-alpha * ({nsigma_tof.GetName()} - offset))', 
                           [nsigma_tof, *pars.values()])
        return pdf, pars
    
    elif function == 'dexp':
        pars = {
            'alpha_0': RooRealVar('alpha_0', 'alpha_0', 0.05, 0., 0.15),
            'offset_0': RooRealVar('offset_0', 'offset_0', -5., -20., 0.),
            'alpha_1': RooRealVar('alpha_1', 'alpha_1', 1., 0.9, 2.),
            'offset_1': RooRealVar('offset_1', 'offset_1', -5., -20., 0.),
        }
        exp2 = RooGenericPdf(f'bkg_exp_0', 'bkg', 
                           f'exp(-alpha_0 * ({nsigma_tof.GetName()} - offset_0))', 
                           [nsigma_tof, pars['alpha_0'], pars['offset_0']])
        exp1 = RooGenericPdf(f'bkg_exp_1', 'bkg', 
                           f'exp(-alpha_1 * ({nsigma_tof.GetName()} - offset_1))', 
                           [nsigma_tof, pars['alpha_1'], pars['offset_1']])
        
        return [exp1, exp2], pars
    
    elif function == 'pol0':
        pars = {
        }
        pol0 = RooGenericPdf('bkg_pol0', 'bkg', '1', [nsigma_tof])
        
        return pol0, pars
    
    else: 
        raise ValueError(f'Unknown function: {function}. Supported functions are "gausexp" and "gaus".')

def fit_parametrisation(fit_results: pd.DataFrame, sign: str, outfile: TFile):

    g_mean = create_graph(fit_results, 'pt', 'mean', 'pt_err', 'mean_err', 
                            f'g_mean_{sign}', ';#it{p}_{T} (GeV/#it{c}); #mu_{TOF}')
    
    PT_MIN, PT_MAX = 0.4, 4
    f_mean = TF1('f_mean', '[0] + [1]*x + [2]*x^2', PT_MIN, PT_MAX, 5)
    f_mean.SetParameters(0.94, 0.005, 0.01)
    g_mean.Fit(f_mean, 'RMS+')
    fit_params = [f_mean.GetParameter(iparam) for iparam in range(3)]
    
    g_sigma = create_graph(fit_results, 'pt', 'sigma', 'pt_err', 'sigma_err',
                            f'g_sigma_{sign}', ';#it{p}_{T} (GeV/#it{c}); #sigma_{TOF}')
    f_sigma = TF1('f_sigma', '[0]', PT_MIN, PT_MAX)
    f_sigma.SetParameters(1)
    g_sigma.Fit(f_sigma, 'RMS+')

    g_resolution = create_graph(fit_results, 'pt', 'resolution', 'pt_err', 'resolution_err', 
                                f'g_resolution_{sign}', ';#it{p}_{T} (GeV/#it{c});#sigma_{TOF} / #it{m}_{TOF}')

    PT_MIN_RES, PT_MAX_RES = 0.7, 4
    f_resolution = TF1('f_resolution', '[0] + [1]*x', PT_MIN_RES, PT_MAX_RES)
    f_resolution.SetParameters(0., 0.025)
    g_resolution.Fit(f_resolution, 'RMS+')
    resolution_params = [f_resolution.GetParameter(iparam) for iparam in range(2)]
        
    outfile.cd()
    g_mean.Write()
    g_resolution.Write()
    g_sigma.Write()

    return fit_params, resolution_params

def _bin_calibration(h2_tof: TH2F, pt_bin: float, pt_step: float, 
                     signal_pdf: RooAbsPdf, signal_pars: dict, bkg_pdf: RooAbsPdf, bkg_pars: dict, 
                     tof_mass: RooRealVar, fit_results: list, tof_dir: TDirectory, sign: str):

    pt = h2_tof.GetXaxis().GetBinCenter(pt_bin)
    pt_low_edge = h2_tof.GetXaxis().GetBinLowEdge(pt_bin)
    pt_high_edge = h2_tof.GetXaxis().GetBinLowEdge(pt_bin+1)
    
    h_tof = h2_tof.ProjectionY(f'tof_mass_{pt:.2f}', pt_bin, pt_bin, 'e')
    if h_tof.GetEntries() < 30:
        return

    if np.abs(pt) > PT_BKG_FIT[sign]:
        sig_frac = RooRealVar('sig_frac', 'sig_frac', 0.7, 0., 1.)
        bkg_fracs = []
        if isinstance(bkg_pdf, list):
            for idx in range(len(bkg_pdf)-1):
                bkg_frac = RooRealVar(f'bkg_frac_{idx}', f'bkg_frac_{idx}', 0.5, 0., 1.)
                bkg_fracs.append(bkg_frac)
            model = RooAddPdf('model', 'model', [signal_pdf, *bkg_pdf], [sig_frac, *bkg_fracs])
        else:
            model = RooAddPdf('model', 'model', [signal_pdf, bkg_pdf], [sig_frac])
        
        #estimation = initialize_means_and_covariances(h_tof, 2, method='kmeans')
        estimation = initialize_means_and_covariances(h_tof, 1, method='kmeans')
        if estimation is not None:
            means, covariances = estimation
            signal_pars['mean'].setVal(means[0])
            signal_pars['sigma'].setVal(np.sqrt(covariances[0]))
            #signal_pars['mean'].setVal(means[1])
            #signal_pars['sigma'].setVal(np.sqrt(covariances[1]))
            #bkg_pars['mean'].setVal(means[0])
            #bkg_pars['sigma'].setVal(np.sqrt(covariances[0]))

    else:
        model = signal_pdf
        estimation = initialize_means_and_covariances(h_tof, 1, method='kmeans')
        if estimation is not None:
            means, covariances = estimation
            signal_pars['mean'].setVal(means[0])
            signal_pars['sigma'].setVal(np.sqrt(covariances[0]))

    iframe, ifit_result = calibration_fit_slice(model, h_tof, tof_mass, signal_pars, pt_low_edge, pt_high_edge)
    ifit_result['pt'] = np.abs(pt)
    ifit_result['pt_err'] = pt_step / 2
    fit_results.append(ifit_result)

    canvas = TCanvas(f'toffit_{pt:.2f}', f'#it{{p}}_{{T}} = {pt:.2f} (GeV/#it{{c}})', 800, 600)
    iframe.Draw()
    tof_dir.cd()
    canvas.Write()

    del iframe, canvas, h_tof

def TOF_calibration(dataset: Dataset, outfile:TFile):

    #axis_spec_tofmass = AxisSpec(100, 0.7, 1.2, 'tof_mass', ';#it{p}_{T} (GeV/#it{c});#it{m}_{TOF} (GeV/#it{c}^{2})')
    axis_spec_nsigma = AxisSpec(100, -3, 3, 'tof_mass', ';#it{p}_{T} (GeV/#it{c});#it{m}_{TOF} (GeV/#it{c}^{2})')
    axis_spec_pt = AxisSpec(160, -8, 8, 'pt', ';#it{p}_{T} (GeV/#it{c});#it{m}_{TOF} (GeV/#it{c}^{2})')
    #h2_tof = dataset.build_th2('fPtHad', 'fMassTOFHad', axis_spec_pt, axis_spec_tofmass)
    h2_tof = dataset.build_th2('fPtHad', 'fNSigmaTOFHadPr', axis_spec_pt, axis_spec_nsigma)
    
    #tof_mass = RooRealVar('fMassTOF', '#it{m}_{TOF} (GeV/#it{c}^{2})', 0.5, 1.5,)
    nsigma_tof = RooRealVar('fNSigmaTOF', 'n#sigma_{TOF}', -3, 3)
    signal_pdf, signal_pars = init_signal_roofit(nsigma_tof, function='gausdexp')

    fit_params, resolution_params, model = {}, {}, None

    for sign in ['matter', 'antimatter']:

        if sign == 'matter':
            slice_range = [PT_MIN, PT_MAX]
        else:
            slice_range = [-PT_MAX, -PT_MIN]

        fit_results = []

        tof_dir = outfile.mkdir(f'TOF_{sign}')

        pt_bin_min = h2_tof.GetXaxis().FindBin(slice_range[0])
        pt_bin_max = h2_tof.GetXaxis().FindBin(slice_range[1])
        pt_step = h2_tof.GetXaxis().GetBinWidth(1)

        for pt_bin in range(pt_bin_min, pt_bin_max):
            
            bkg_pdf_func = 'pol0' if np.abs(h2_tof.GetXaxis().GetBinCenter(pt_bin)) < 2.8 else 'exp'
            bkg_pdf, bkg_pars = init_background_roofit(nsigma_tof, function=bkg_pdf_func)
            
            _bin_calibration(h2_tof, pt_bin, pt_step,
                             signal_pdf, signal_pars, bkg_pdf, bkg_pars,
                             nsigma_tof, fit_results, tof_dir, sign)
            
        fit_results_df = pd.DataFrame(fit_results)
        
        fit_params[sign], resolution_params[sign] = fit_parametrisation(fit_results_df, sign, tof_dir)

    return fit_params['matter'], resolution_params['matter']

def get_resolution(pt:float, resolution_params:list):
    return resolution_params[0] + resolution_params[1] * pt

def get_expected_tof_mass(pt:float, fit_params:list):
    return fit_params[0] + fit_params[1] * pt + fit_params[2] * pt**2

def visualize_distributions_and_fit(dataset: Dataset, outfile: TFile, fit_params:list, resolution_params:float):

    dataset['fExpTOFSignal'] = get_expected_tof_mass(np.abs(dataset['fPtHad'].values), fit_params)
    dataset['fTOFResolution'] = get_resolution(np.abs(dataset['fPtHad'].values), resolution_params)
    dataset['fNSigmaTOF'] = (dataset['fMassTOFHad'] - dataset['fExpTOFSignal']) / (dataset['fExpTOFSignal'] * dataset['fTOFResolution'])

    axis_spec_pt = AxisSpec(160, -8, 8, 'pt', ';#it{p}_{T} (GeV/#it{c});#it{m}_{TOF} (GeV/#it{c}^{2})')
    axis_spec_tofmass = AxisSpec(100, 0.7, 1.2, 'tof_mass', ';#it{p}_{T} (GeV/#it{c});#it{m}_{TOF} (GeV/#it{c}^{2})')
    axis_spec_nsigmatof = AxisSpec(100, -5, 5, 'nsigma_tof', ';#it{p}_{T} (GeV/#it{c});n#sigma_{TOF}')

    h2_nsigmatof = dataset.build_th2('fPtHad', 'fNSigmaTOF', axis_spec_pt, axis_spec_nsigmatof, title=';#it{p}_{T} (GeV/#it{c});n#sigma_{TOF}')
    h2_exptof = dataset.build_th2('fPtHad', 'fExpTOFSignal', axis_spec_pt, axis_spec_tofmass, title=';#it{p}_{T} (GeV/#it{c});#mathbb{E}(#it{m}_{TOF}) (GeV/#it{c}^{2})')
    h2_tof = dataset.build_th2('fPtHad', 'fMassTOFHad', axis_spec_pt, axis_spec_tofmass, title=';#it{p}_{T} (GeV/#it{c});#it{m}_{TOF} (GeV/#it{c}^{2})')

    f_fit_matter = TF1('f_fit_matter', '[0] * std::exp(- std::abs(x) * [1]) + [2]', PT_MIN, PT_MAX, 5)
    f_fit_matter.SetParameters(*fit_params)
    f_fit_antimatter = TF1('f_fit_antimatter', '[0] * std::exp(- std::abs(x) * [1]) + [2]', -PT_MAX, -PT_MIN, 5)

    outfile.cd()
    
    h2_tof.Write()
    h2_nsigmatof.Write()
    h2_exptof.Write('exp_tof_mass')
    
    canvas = TCanvas('cNSigmatof', 'cNSigmatof', 800, 600)
    h2_tof.Draw('colz')
    f_fit_matter.Draw('same')
    f_fit_antimatter.Draw('same')
    canvas.Write()


if __name__ == '__main__':

    #infile_path = '/data/galucia/lithium_local/MC/LHC25a4_with_primaries_and_mixed.root',
    infile_path = [#'/data/galucia/lithium_local/same/LHC23_PbPb_pass5_same.root',
                   '/data/galucia/lithium_local/same/LHC24ar_pass3_same.root',]
    folder_name = 'DF*'
    tree_name = 'O2he3hadtable'
    dataset = Dataset.from_root(infile_path, tree_name, folder_name, columns=['fNSigmaTOFHadPr', 'fPtHad'])
    print(dataset.columns)

    #['fPtHe3', 'fEtaHe3', 'fPhiHe3', 'fPtHad', 'fEtaHad', 'fPhiHad',
    #   'fDCAxyHe3', 'fDCAzHe3', 'fDCAxyHad', 'fDCAzHad', 'fDCApair',
    #   'fSignalTPCHe3', 'fInnerParamTPCHe3', 'fSignalTPCHad',
    #   'fInnerParamTPCHad', 'fNClsTPCHe3', 'fNSigmaTPCHe3', 'fNSigmaTPCHadPi',
    #   'fNSigmaTPCHadKa', 'fNSigmaTPCHadPr', 'fNSigmaTOFHadPi',
    #   'fNSigmaTOFHadKa', 'fNSigmaTOFHadPr', 'fChi2TPCHe3', 'fChi2TPCHad',
    #   'fMassTOFHe3', 'fMassTOFHad', 'fPIDtrkHe3', 'fPIDtrkHad',
    #   'fItsClusterSizeHe3', 'fItsClusterSizeHad', 'fSharedClustersHe3',
    #   'fSharedClustersHad'],

    #exit(0)

    #dataset.query(f'0.5 < fChi2TPCHad < 4', inplace=True)
    
    outfile = TFile('output/TOF_pr_from_central.root', 'recreate')

    fit_params, resolution_params = TOF_calibration(dataset, outfile)
    #visualize_distributions_and_fit(dataset, outfile, fit_params, resolution_params)

    outfile.Close()
    