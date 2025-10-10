'''
    Code to run calibration of ITS and TPC parametrisations
'''

import numpy as np
import pandas as pd
from ROOT import TFile, TF1, TCanvas
from ROOT import RooRealVar, RooCrystalBall, RooGaussian, RooAddPdf

from torchic import Dataset, AxisSpec
from torchic.roopdf import RooGausExp
from torchic.core.graph import create_graph
from torchic.physics.ITS import average_cluster_size, expected_cluster_size, sigma_its

from particle import Particle

import sys
sys.path.append('..')
from utils.pid_routine import PARTICLE_ID, PDG_CODE
from utils.utils import calibration_fit_slice, initialize_means_and_covariances

CONF = {
    'beta_gamma': {
        'x_min': 0,
        'x_max': 5,

        'Pi': {
            'x_min_fit': 1.5,
            'x_max_fit': 5,
            'x_nbins': 50,
        },
        'Ka': {
            'x_min_fit': 0.7,
            'x_max_fit': 4,
            'x_nbins': 50,
        },
        'Pr': {
            'x_min_fit': 0.5,
            'x_max_fit': 4.5,
            'x_max_bkg': 1.1,
            'x_nbins': 50,
        },
        'De': {
            'x_min_fit': 0.4,
            'x_max_fit': 3.5,
            'x_nbins': 50,
        },
        'He': {
            'x_min_fit': 0.7,
            'x_max_fit': 3,
            'x_max_bkg': 2.6,
            'x_nbins': 50,
        }
    },
    'p': {
        'x_min': 0,
        'x_max': 6,

        'Pi': {
            'x_min_fit': 0,
            'x_max_fit': 1.3,
            'x_nbins': 120
        },
        'Ka': {
            'x_min_fit': 0,
            'x_max_fit': 3,
            'x_nbins': 60
        },
        'Pr': {
            'x_min_fit': 0.,
            'x_max_fit': 4.,
            'x_nbins': 60
        },
        'De': {
            'x_min_fit': 0.,
            'x_max_fit': 6,
            'x_nbins': 60
        },
        'He': {
            'x_min_fit': 1.7,
            'x_max_fit': 5.9,
            'x_max_bkg': 6,
            'x_nbins': 60
        }
    }
}

DATASET_COLUMN_NAMES = {
    'P':'fP',
    'Pt':'fPt',
    'Eta':'fEta',
    'Phi':'fPhi',
    'TofNSigma':'fTofNSigma',
    'Chi2TPC':'fChi2TPC',
    'PIDtracking': 'fPIDinTrk',
    'ItsClusterSize':'fItsClusterSize',
}

SUFFIX = {
    'Pr': 'Had',
    'He': 'He3'
}

def prepare_dataset(dataset: Dataset, particle: str, mode: str = 'truncated'):

    dataset[f'fItsClusterSize{SUFFIX[particle]}'] = np.array(dataset[f'fItsClusterSize{SUFFIX[particle]}'], np.uint64)
    dataset[f'fAvgClSizeCosLam{SUFFIX[particle]}'], dataset[f'fNHitsIts{SUFFIX[particle]}'] = np.zeros(dataset.shape[0], dtype=float), np.zeros(dataset.shape[0], dtype=int)
    if mode == 'truncated':
        dataset[f'fAvgClusterSize{SUFFIX[particle]}'], dataset[f'fNHitsIts{SUFFIX[particle]}'] = average_cluster_size(dataset[f'fItsClusterSize{SUFFIX[particle]}'], do_truncated=True)
    elif mode == 'mean':
        dataset[f'fAvgClusterSize{SUFFIX[particle]}'], dataset[f'fNHitsIts{SUFFIX[particle]}'] = average_cluster_size(dataset[f'fItsClusterSize{SUFFIX[particle]}'], do_truncated=False)

    dataset.query(f'fNHitsIts{SUFFIX[particle]} > 5', inplace=True)
    dataset[f'fAvgClSizeCosLam{SUFFIX[particle]}'] = dataset[f'fAvgClusterSize{SUFFIX[particle]}'] / np.cosh(dataset[f'fEta{SUFFIX[particle]}'])
    dataset[f'fBetaGamma{SUFFIX[particle]}'] = abs(dataset[f'fP{SUFFIX[particle]}']) / (Particle.from_pdgid(PDG_CODE[particle]).mass / 1_000)
    dataset[f'fP{SUFFIX[particle]}'] = np.abs(dataset[f'fP{SUFFIX[particle]}'])

def init_signal_roofit(clsize: RooRealVar, function: str = 'crystalball'):

    if function == 'crystalball':
        signal_pars = {
            'mean': RooRealVar('mean', 'mean', 1., 15, ''),
            'sigma': RooRealVar('sigma', 'sigma', 0.01, 10, ''),
            'aL': RooRealVar('aL', 'aL', 0.7, 30.),
            'nL': RooRealVar('nL', 'nL', 0.3, 30.),
            'aR': RooRealVar('aR', 'aR', 0.7, 30.),
            'nR': RooRealVar('nR', 'nR', 0.3, 30.),
        }
        signal = RooCrystalBall('signal', 'signal', clsize, signal_pars['mean'], signal_pars['sigma'],
                                signal_pars['aL'], signal_pars['nL'], doubleSided=True) #
                                #signal_pars['aR'], signal_pars['nR'])

        return signal, signal_pars
    
    elif function == 'gausexp':
        signal_pars = {
            'mean': RooRealVar('mean', 'mean', 1., 15, ''),
            'sigma': RooRealVar('sigma', 'sigma', 0.01, 10, ''),
            'rlife': RooRealVar('rlife', 'rlife', 2., 0., 10.),
        }
        signal = RooGausExp('signal', 'signal', clsize, *signal_pars.values())
        return signal, signal_pars
    
    elif function == 'gaus':
        signal_pars = {
            'mean': RooRealVar('mean', 'mean', 1., 15, ''),
            'sigma': RooRealVar('sigma', 'sigma', 0.01, 10, ''),
        }
        signal = RooGaussian('signal', 'signal', clsize, *signal_pars.values())
        return signal, signal_pars
    
    else:
        raise ValueError(f'Unknown function: {function}. Supported functions are "crystalball" and "gausexp".')

def init_background_roofit(clsize: RooRealVar, particle: str, function: str = 'gaus'):

    if function == 'gausexp':
        bkg_pars = {
            'mean': RooRealVar('bkg_mean', 'bkg_mean', 0., 3.5, ''),
            'sigma': RooRealVar('bkg_sigma', 'bkg_sigma', 0.1, 0.8, ''),
            'rlife': RooRealVar('bkg_rlife', 'rlife', 2., 0., 10.),
        }
        if particle == 'He':
            bkg_pars['mean'] = RooRealVar('bkg_mean', 'bkg_mean', 0., 3, '')
        bkg = RooGausExp('bkg', 'bkg', clsize, *bkg_pars.values())
        return bkg, bkg_pars
    elif function == 'gaus':
        bkg_pars = {
            'mean': RooRealVar('bkg_mean', 'bkg_mean', 0., 1, ''),
            'sigma': RooRealVar('bkg_sigma', 'bkg_sigma', 0.1, 0.8, ''),
            #'rlife': RooRealVar('rlife', 'rlife', 0., 10.),
        }
        if particle == 'He':
            bkg_pars['mean'] = RooRealVar('bkg_mean', 'bkg_mean', 0., 3, '')
        bkg = RooGaussian('bkg', 'bkg', clsize, bkg_pars['mean'], bkg_pars['sigma'])

        return bkg, bkg_pars
    else:
        raise ValueError(f'Unknown function: {function}. Supported functions are "gausexp" and "gaus".')


X_DICT = {
    'beta_gamma': {
        'axis_name': 'bg',
        'axis_title': '#beta#gamma',
        'var_name': 'fBetaGamma',
    },
    'p': {
        'axis_name': 'p',
        'axis_title': '#it{p} (GeV/c)',
        'var_name': 'fP',
    }
}

def visualize_fit_results(dataset, fit_results_df, particle, bg_min, bg_max, particle_dir, x:str):

    x_dict = X_DICT[x]
    g_mean = create_graph(fit_results_df, 'x', 'mean', 'x_error', 'mean_err', 
                                f'g_mean', f';{x_dict["axis_title"]};#LT ITS Cluster Size #GT #times cos #LT #lambda #GT')
    f_mean = TF1('simil_bethe_bloch_func', '[0]/x^[1] + [2]', bg_min, bg_max)
    f_mean.SetParameters(2.6, 3.6, 2)
    f_mean.SetParLimits(0, 0, 3)
    if particle == 'He':
        f_mean.SetParameters(2.3, 1.7, 4.5)
    g_mean.Fit(f_mean, 'RMS+')
    c_mean = TCanvas('c_mean', 'c_mean', 800, 600)
    g_mean.Draw('ap')
    f_mean.Draw('same')

    g_sigma = create_graph(fit_results_df, 'x', 'sigma', 'x_error', 'sigma_err', 
                                f'g_sigma', f';{x_dict["axis_title"]};#sigma / #mu')

    g_resolution = create_graph(fit_results_df, 'x', 'resolution', 'x_error', 'resolution_err', 
                                f'g_resolution', f';{x_dict["axis_title"]};#sigma / #mu')
    
    #f_resolution = TF1('resolution_fit', '[0]*ROOT::Math::erf((x - [1])/[2])', bg_min, bg_max)
    #f_resolution.SetParameters(0.24, -0.32, 1.53)
    f_resolution = TF1('resolution_fit', '[0]', bg_min, bg_max)
    f_resolution.SetParameters(0.155)
    
    if particle == 'He':    
        f_resolution = TF1('resolution_fit', '[0] + x*[1] + x*x * [2]', bg_min, bg_max)
        f_resolution.SetParameters(0.116, -0.01/1.7, 0)
        #f_resolution.FixParameter(2, 0)
    
    g_resolution.Fit(f_resolution, 'RMS+')
    c_resolution = TCanvas('c_resolution', 'c_resolution', 800, 600)
    g_resolution.Draw('ap')
    f_resolution.Draw('same')
    
    pid_params = (f_mean.GetParameter(0), f_mean.GetParameter(1), f_mean.GetParameter(2),
                  f_resolution.GetParameter(0), f_resolution.GetParameter(1), f_resolution.GetParameter(2))
    dataset[f'fExpClSizeCosLam{SUFFIX[particle]}'] = expected_cluster_size(dataset[f'fBetaGamma{SUFFIX[particle]}'], pid_params)
    if particle == 'He':    
        dataset[f'fSigmaITS{SUFFIX[particle]}'] = dataset[f'fExpClSizeCosLam{SUFFIX[particle]}'] * (pid_params[3] + pid_params[4]*dataset[f'fBetaGamma{SUFFIX[particle]}'])
    else:                   
        #dataset[f'fSigmaITS{SUFFIX[particle]}'] = sigma_its(dataset[f'fBetaGamma{SUFFIX[particle]}'], pid_params)
        dataset[f'fSigmaITS{SUFFIX[particle]}'] = pid_params[3] * dataset[f'fExpClSizeCosLam{SUFFIX[particle]}']
    dataset[f'fNSigmaITS{SUFFIX[particle]}'] = (dataset[f'fAvgClSizeCosLam{SUFFIX[particle]}'] - dataset[f'fExpClSizeCosLam{SUFFIX[particle]}']) / dataset[f'fSigmaITS{SUFFIX[particle]}']

    axis_spec_bg = AxisSpec(50, 0, 5, x_dict['axis_name'], ';;')
    axis_spec_nsigma = AxisSpec(100, -5, 5, 'nsigma', f';{x_dict["axis_title"]};n#sigma_{{ITS}}')
    h2_nsigma = dataset.build_th2(f'fBetaGamma{SUFFIX[particle]}', f'fNSigmaITS{SUFFIX[particle]}', axis_spec_bg, axis_spec_nsigma)

    particle_dir.cd()
    c_mean.Write()
    g_sigma.Write()
    c_resolution.Write()
    h2_nsigma.Write()

def calibration_routine(dataset: Dataset, outfile: TFile, particle: str, x: str = 'beta_gamma'): 

    x_dict = X_DICT[x]

    cfg = CONF[x][particle]

    axis_spec_x = AxisSpec(cfg['x_nbins'], CONF[x]['x_min'], CONF[x]['x_max'], x_dict['axis_name'], ';;')
    axis_spec_clsize = AxisSpec(60, 0, 15, 'cluster_size_cal', f';{x_dict["axis_title"]};#LT ITS Cluster Size #GT #times cos #LT #lambda #GT')

    h2_clsize = dataset.build_th2(x_dict['var_name']+SUFFIX[particle], f'fAvgClSizeCosLam{SUFFIX[particle]}',
                                  axis_spec_x, axis_spec_clsize)

    # RooFit initialization
    clsize = RooRealVar('fClSizeCosLam', '#LT Cluster size #GT #LT cos#lambda #GT', 0., 15.)
    signal, signal_pars = init_signal_roofit(clsize, function='gausexp')
    bkg, bkg_pars = init_background_roofit(clsize, particle, function='gausexp')

    x_min = cfg['x_min_fit']
    x_max = cfg['x_max_fit']

    fit_results_df = None

    x_bin_min = h2_clsize.GetXaxis().FindBin(x_min)
    x_bin_max = h2_clsize.GetXaxis().FindBin(x_max)
    for x_bin in range(x_bin_min, x_bin_max+1):
        
        ix = h2_clsize.GetXaxis().GetBinCenter(x_bin)
        x_error = h2_clsize.GetXaxis().GetBinWidth(x_bin) / 2.
        x_low_edge = h2_clsize.GetXaxis().GetBinLowEdge(x_bin)
        x_high_edge = h2_clsize.GetXaxis().GetBinLowEdge(x_bin+1)
        
        h_clsize = h2_clsize.ProjectionY(f'clsize_{ix:.2f}', x_bin, x_bin, 'e')
        if h_clsize.GetEntries() <= 0:
            print(f'No entries for particle {particle}, {x_dict["axis_name"]} = {ix:.2f}, skipping...')
            continue

        model = None
        if (particle == 'Pr' and ix < cfg.get('x_max_bkg', 0)):
            sig_frac = RooRealVar('sig_frac', 'sig_frac', 0.5, 0., 1.)
            model = RooAddPdf('model', 'signal + bkg', [signal, bkg], [sig_frac])

            if h_clsize.GetEntries() > 30:
                means, covariances = initialize_means_and_covariances(h_clsize, 2, method='kmeans')
                mean_sig, sigma_sig = (means[1], np.sqrt(covariances[1]))
                mean_bkg, sigma_bkg = (means[0], np.sqrt(covariances[0]))

                signal_pars['mean'].setVal(mean_sig)
                signal_pars['sigma'].setVal(sigma_sig)
                bkg_pars['mean'].setVal(mean_bkg)
                bkg_pars['sigma'].setVal(sigma_bkg)

        else:
            model = signal
            if h_clsize.GetEntries() > 30:
                means, sigmas = initialize_means_and_covariances(h_clsize, 1, method='kmeans')
                signal_pars['mean'].setVal(means[0])
                signal_pars['sigma'].setVal(np.sqrt(sigmas[0]))

        frame, fit_results = calibration_fit_slice(model, h_clsize, clsize, signal_pars, x_low_edge, x_high_edge)
        fit_results['x'] = np.abs(ix)
        fit_results['x_error'] = x_error
        if fit_results_df is None:
            fit_results_df = pd.DataFrame.from_dict([fit_results])
        else:
            fit_results_df = pd.concat([fit_results_df, pd.DataFrame.from_dict([fit_results])], ignore_index=True)

        canvas = TCanvas(f'cClSizeCosLam_{ix:.2f}', f'cClSizeCosLam_{ix:.2f}', 800, 600)
        frame.Draw()
        outfile.cd()
        canvas.Write()

    if fit_results_df is None:
        print(f'No fit results for particle {particle}, skipping...')
        return
    visualize_fit_results(dataset, fit_results_df, particle, x_min, x_max, outfile, x)

    outfile.cd()
    h2_clsize.Write()

    del h2_clsize, clsize, signal, signal_pars, bkg, bkg_pars

def main_routine(dataset: Dataset, mode: str, x: str):

    outfile_path = f'output/ITS.root'
    outfile = TFile(outfile_path, 'RECREATE')
    
    particles = ['Pr', 'He']
    dataset.query(f'fPIDtrkHe3 == 7', inplace=True)
    dataset['fPHad'] = dataset['fPtHad'] * np.cosh(dataset['fEtaHad'])
    dataset['fPHe3'] = dataset['fPtHe3'] * np.cosh(dataset['fEtaHe3'])
    
    for particle in particles:
        
        prepare_dataset(dataset, particle, mode=mode)
        particle_dir = outfile.mkdir(particle)
        calibration_routine(dataset, particle_dir, particle, x)

    outfile.Close()


if __name__ == '__main__':

    infile_path = ['/data/galucia/lithium_local/same/LHC23_PbPb_pass5_same.root',
                   '/data/galucia/lithium_local/same/LHC24_PbPb_pass2_same.root']

    folder_name = 'DF*'
    tree_name = 'O2he3hadtable'

    
    dataset = Dataset.from_root(infile_path, tree_name=tree_name, folder_name=folder_name)

    print(f'{dataset.shape=}\n{dataset.columns=}')

    main_routine(dataset, mode='mean', x='beta_gamma')