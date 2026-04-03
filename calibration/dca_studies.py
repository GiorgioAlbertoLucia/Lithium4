'''
    Study of the dca distribution for p-He3 pairs from MC
'''

import numpy as np
import pandas as pd
from ROOT import TFile, TDirectory, RooRealVar, RooGaussian, RooCrystalBall, TCanvas, TF1, RooAddPdf, RooChebychev, RooPolynomial, gStyle
from torchic import Dataset, AxisSpec
from torchic.roopdf import RooGausExp   
from torchic.core.fit import calibration_fit_slice, initialize_means_and_covariances
from torchic.core.graph import create_graph

FLAGS = ['BothPrimaries', 'BothFromLi4', 'BothFromHypertriton', 'MixedPairs']

CONFIG = {
    'He3': {
        'xy': {
            'axis_title': 'DCA_{#it{xy}} (cm)',
            'x_min_fit': -4.,
            'x_max_fit': -1.5,
            'x_max_bkg': 1.5,
            'y_nbins': 50,
            'y_min': -0.05,
            'y_max': 0.05,
        },
        'z': {
            'axis_title': 'DCA_{#it{z}} (cm)',
            'x_min_fit': -4.,
            'x_max_fit': -1.5,
            'x_max_bkg': 1.5,
            'y_min': -0.15,
            'y_max': 0.15,
            'y_nbins': 50
        }
    },
    'Had': {
        'xy': {
            'axis_title': 'DCA_{#it{xy}} (cm)',
            'x_min_fit': -4.,
            'x_max_fit': -0.4,
            'x_max_bkg': 1.5,
            'y_nbins': 50,
            'y_min': -0.05,
            'y_max': 0.05,
        },
        'z': {
            'axis_title': 'DCA_{#it{z}} (cm)',
            'x_min_fit': -4.0,
            'x_max_fit': -0.4,
            'x_max_bkg': 1.5,
            'y_min': -0.15,
            'y_max': 0.15,
            'y_nbins': 50
        }
    }
}

def visualise(dataset:Dataset, outfile:TDirectory, particle:str = 'He3', is_mc:bool = True):

    cfg = CONFIG[particle]
    axis_spec_pt = AxisSpec(100, -5, 5, 'pt', '')
    axis_spec_dcaxy = AxisSpec(cfg['xy']['y_nbins'], cfg['xy']['y_min'], cfg['xy']['y_max'], 'DCAxy', ';#it{p}_{T};DCA_{xy} (cm)')
    axis_spec_dcaz = AxisSpec(cfg['z']['y_nbins'], cfg['z']['y_min'], cfg['z']['y_max'], 'DCAz', ';#it{p}_{T};DCA_{z} (cm)')

    h2_dcaxy_pt = dataset.build_th2(f'fPt{particle}', f'fDCAxy{particle}', axis_spec_pt, axis_spec_dcaxy, name='h2DCAxyPt', title=';#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm)')
    h2_dcaz_pt = dataset.build_th2(f'fPt{particle}', f'fDCAz{particle}', axis_spec_pt, axis_spec_dcaz, name='h2DCAzPt', title=';#it{p}_{T} (GeV/#it{c});DCA_{z} (cm)')

    outfile.cd()
    h2_dcaxy_pt.Write()
    h2_dcaz_pt.Write()

    if is_mc:
        for flag in FLAGS:
            h_dcaxy_flag = dataset.build_th1(f'fDCAxy{particle}', axis_spec_dcaxy, subset=flag, name=f'hDCAxy{particle}{flag}', title=';DCA_{xy} (cm);Counts')
            h_dcaz_flag = dataset.build_th1(f'fDCAz{particle}', axis_spec_dcaz, subset=flag, name=f'hDCAz{particle}{flag}', title= ';DCA_{z} (cm);Counts')
            h_dcaxy_flag.Write()
            h_dcaz_flag.Write()

    return h2_dcaxy_pt, h2_dcaz_pt

def init_signal_roofit(x: RooRealVar, function: str = 'crystalball'):

    if function == 'crystalball':
        signal_pars = {
            'mean': RooRealVar('mean', 'mean', -0.5, 0.5, ''),
            'sigma': RooRealVar('sigma', 'sigma', 0.0001, 0.01, ''),
            'aL': RooRealVar('aL', 'aL', 0.7, 30.),
            'nL': RooRealVar('nL', 'nL', 0.3, 30.),
            'aR': RooRealVar('aR', 'aR', 0.7, 30.),
            'nR': RooRealVar('nR', 'nR', 0.3, 30.),
        }
        signal = RooCrystalBall('signal', 'signal', x, signal_pars['mean'], signal_pars['sigma'],
                                signal_pars['aL'], signal_pars['nL'], doubleSided=True) #
                                #signal_pars['aR'], signal_pars['nR'])

        return signal, signal_pars
    
    elif function == 'gausexp':
        signal_pars = {
            'mean': RooRealVar('mean', 'mean', -0.5, 0.5, ''),
            'sigma': RooRealVar('sigma', 'sigma', 0.0001, 0.1, ''),
            'rlife': RooRealVar('rlife', 'rlife', 2., 0., 10.),
        }
        signal = RooGausExp('signal', 'signal', x, *signal_pars.values())
        return signal, signal_pars
    
    elif function == 'gaus':
        signal_pars = {
            'mean': RooRealVar('mean', 'mean', -0.5, 0.5, ''),
            'sigma': RooRealVar('sigma', 'sigma', 0.0001, 0.1, ''),
        }
        signal = RooGaussian('signal', 'signal', x, *signal_pars.values())
        return signal, signal_pars
    
    else:
        raise ValueError(f'Unknown function: {function}. Supported functions are "crystalball" and "gausexp".')

def init_background_roofit(x: RooRealVar, particle: str, function: str = 'gaus'):

    if function == 'crystalball':
        bkg_pars = {
            'mean': RooRealVar('bkg_mean', 'bkg_mean', -0.5, 0.5, ''),
            'sigma': RooRealVar('bkg_sigma', 'bkg_sigma', 0.0001, 0.01, ''),
            'aL': RooRealVar('bkg_aL', 'bkg_aL', 0.7, 30.),
            'nL': RooRealVar('bkg_nL', 'bkg_nL', 0.3, 30.),
            'aR': RooRealVar('bkg_aR', 'bkg_aR', 0.7, 30.),
            'nR': RooRealVar('bkg_nR', 'bkg_nR', 0.3, 30.),
        }
        bkg = RooCrystalBall('bkg', 'bkg', x, bkg_pars['mean'], bkg_pars['sigma'],
                                bkg_pars['aL'], bkg_pars['nL'], doubleSided=True) #
                                #bkg_pars['aR'], bkg_pars['nR'])

        return bkg, bkg_pars

    elif function == 'gausexp':
        bkg_pars = {
            'mean': RooRealVar('bkg_mean', 'bkg_mean', 0., 1, ''),
            'sigma': RooRealVar('bkg_sigma', 'bkg_sigma', 0.1, 0.8, ''),
            'rlife': RooRealVar('bkg_rlife', 'rlife', 2., 0., 10.),
        }
        if particle == 'He':
            bkg_pars['mean'] = RooRealVar('bkg_mean', 'bkg_mean', 0., 3, '')
        bkg = RooGausExp('bkg', 'bkg', x, *bkg_pars.values())
        return bkg, bkg_pars
    
    elif function == 'gaus':
        bkg_pars = {
            'mean': RooRealVar('bkg_mean', 'bkg_mean', 0., 1, ''),
            'sigma': RooRealVar('bkg_sigma', 'bkg_sigma', 0.1, 0.8, ''),
            #'rlife': RooRealVar('rlife', 'rlife', 0., 10.),
        }
        if particle == 'He':
            bkg_pars['mean'] = RooRealVar('bkg_mean', 'bkg_mean', 0., 3, '')
        bkg = RooGaussian('bkg', 'bkg', x, bkg_pars['mean'], bkg_pars['sigma'])

        return bkg, bkg_pars
    
    elif function == 'cheb1':
        bkg_pars = {
            'p0': RooRealVar('bkg_p0', 'bkg_p0', -0.1, -1, 1, ''),
            'p1': RooRealVar('bkg_p1', 'bkg_p1', 0.001, -1, 1, ''),
            #'rlife': RooRealVar('rlife', 'rlife', 0., 10.),
        }
        bkg = RooChebychev('bkg', 'bkg', x, list(bkg_pars.values()))

        return bkg, bkg_pars
    
    elif function == 'pol0':
        pars = {
            'c0': RooRealVar('c0{bkg}', 'c_{0}^{bkg}', 1.0, 0.0, 10.0)
        }

        pdf = RooPolynomial('pdf_{bkg}', 'pdf_{bkg}', x, [pars['c0']])
        return pdf, pars

    else:
        raise ValueError(f'Unknown function: {function}. Supported functions are "gausexp" and "gaus".')

def visualize_fit_results(fit_results_df, particle, particle_dir, x:str):

    x_antimatter = [x for x in fit_results_df['x'] if x < 0]
    x_min = min(x_antimatter) if len(x_antimatter) > 0 else -1
    x_max = max(x_antimatter) if len(x_antimatter) > 0 else -0.1
    
    g_mean = create_graph(fit_results_df, 'x', 'mean', 'x_error', 'mean_err', 
                                f'g_mean', f'{particle};#it{{p}}_{{T}} (GeV/#it{{c}});{CONFIG[particle][x]["axis_title"]}')
    f_mean = TF1('f_mean', '[0]', x_min, x_max)
    f_mean.SetParameter(0, 0.)
    g_mean.Fit(f_mean, 'RMS+')
    c_mean = TCanvas('c_mean', 'c_mean', 800, 600)
    g_mean.Draw('ap')
    f_mean.Draw('same')

    g_sigma = create_graph(fit_results_df, 'x', 'sigma', 'x_error', 'sigma_err', 
                                f'g_sigma', f'{particle};#it{{p}}_{{T}} (GeV/#it{{c}});#sigma {CONFIG[particle][x]["axis_title"]}')
    #f_sigma = TF1('f_sigma', '[0]*exp(-abs(x)*[1]) + [2]', x_min, x_max)
    f_sigma = TF1('f_sigma', 'sqrt([0]*[0] + ([1]/x)*([1]/x) + ([2]/(x*x))*([2]/(x*x)))', x_min, x_max)
    #f_sigma = TF1('f_sigma', '[0] + [1]/x^[2]', x_min, x_max)

    #f_sigma.SetParameters(0.005, 0.5, 0.0025)
    f_sigma.SetParameters(1.8e-3, 2.2e-3, 8e-4)
    #f_sigma.SetParameters(1.8e-3, 2.2e-3, 1.5)
    #f_sigma.SetParLimits(2, 0., 8e-3)
    g_sigma.Fit(f_sigma, 'RMS+')
    c_sigma = TCanvas('c_sigma', 'c_sigma', 800, 600)
    g_sigma.Draw('ap')
    f_sigma.Draw('same')

    particle_dir.cd()
    c_mean.Write()
    c_sigma.Write()

    del g_mean, c_mean, g_sigma, f_sigma

def _bin_calibration(h2, outfile: TFile, particle: str, x_bin: int, model, signal_pars, dca: RooRealVar,
                     fit_results_df: pd.DataFrame = None, x: str = 'xy'):

    ix = h2.GetXaxis().GetBinCenter(x_bin)
    x_error = h2.GetXaxis().GetBinWidth(x_bin) / 2.
    x_low_edge = h2.GetXaxis().GetBinLowEdge(x_bin)
    x_high_edge = h2.GetXaxis().GetBinLowEdge(x_bin+1)
    
    h_dca = h2.ProjectionY(f'pt_{ix:.2f}', x_bin, x_bin, 'e')
    if h_dca.GetEntries() <= 0:
        print(f'No entries for particle {particle}, pt = {ix:.2f}, skipping...')
        return fit_results_df

    if h_dca.GetEntries() > 30:
        means, sigmas = initialize_means_and_covariances(h_dca, 1, method='kmeans')
        signal_pars['mean'].setVal(means[0])
        signal_pars['sigma'].setVal(np.sqrt(sigmas[0]))

    frame, fit_results = calibration_fit_slice(model, h_dca, dca, signal_pars, x_low_edge, x_high_edge)
    fit_results['x'] = ix
    fit_results['x_error'] = x_error
    if fit_results_df is None:
        fit_results_df = pd.DataFrame.from_dict([fit_results])
    else:
        fit_results_df = pd.concat([fit_results_df, pd.DataFrame.from_dict([fit_results])], ignore_index=True)

    canvas = TCanvas(f'cDCA{x}_pt_{ix:.2f}', f'cDCA{x}_pt_{ix:.2f}', 800, 600)
    canvas.SetLogy()
    frame.SetMinimum(1)
    frame.Draw()
    outfile.cd()
    canvas.Write()

    del h_dca, frame, canvas

    return fit_results_df

def calibration_routine(h2, outfile: TFile, particle: str, x: str = 'xy'): 

    cfg = CONFIG[particle][x]

    dca = RooRealVar('fDCA', cfg["axis_title"], cfg['y_min'], cfg['y_max'])
    signal, signal_pars = init_signal_roofit(dca, function='crystalball')
    bkg, bkg_pars = init_background_roofit(dca, particle, function='pol0')

    x_min = cfg['x_min_fit']
    x_max = cfg['x_max_fit']

    fit_results_df = None

    x_bin_min = h2.GetXaxis().FindBin(x_min)
    x_bin_max = h2.GetXaxis().FindBin(x_max)
    for x_bin in range(x_bin_min, x_bin_max+1):
        
        signal_fraction = RooRealVar('signal_fraction', 'signal_fraction', 0.5, 0., 1.)
        model = RooAddPdf('model', 'signal + bkg', [signal, bkg], [signal_fraction])
        
        # antimatter
        ix = h2.GetXaxis().GetBinCenter(x_bin)
        #fit_results_df = _bin_calibration(h2, outfile, particle, x_bin, signal, signal_pars, dca, fit_results_df, x)
        fit_results_df = _bin_calibration(h2, outfile, particle, x_bin, model, signal_pars, dca, fit_results_df, x)

        # matter
        positive_x_bin = h2.GetXaxis().FindBin(-ix)
        for par in ['mean', 'aL', 'nL']:
            signal_pars[par].setConstant(True)
            #bkg_pars[par].setVal(signal_pars[par].getVal())
            #bkg_pars[par].setConstant(True)
        #bkg_pars['sigma'].setVal(3*signal_pars['sigma'].getVal())
        #bkg_pars['sigma'].setRange(1.2*signal_pars['sigma'].getVal(),1)

        fit_results_df = _bin_calibration(h2, outfile, particle, positive_x_bin, model, signal_pars, dca, fit_results_df, x)
        signal_pars['aL'].setConstant(False)
        signal_pars['nL'].setConstant(False)

    if fit_results_df is None:
        print(f'No fit results for particle {particle}, skipping...')
        return
    visualize_fit_results(fit_results_df, particle, outfile, x)

    outfile.cd()
    h2.Write()

    del h2, dca, signal, signal_pars, bkg, bkg_pars


def _bin_calibration_simple_gaus(h2, outfile: TFile, particle: str, x_bin: int,
                                  fit_results_df: pd.DataFrame = None, x: str = 'xy'):
    """
    Helper function for simple Gaussian calibration of a single bin using TF1.
    """
    ix = h2.GetXaxis().GetBinCenter(x_bin)
    x_error = h2.GetXaxis().GetBinWidth(x_bin) / 2.
    
    h_dca = h2.ProjectionY(f'pt_{ix:.2f}', x_bin, x_bin, 'e')
    if h_dca.GetEntries() <= 0:
        print(f'No entries for particle {particle}, pt = {ix:.2f}, skipping...')
        return fit_results_df
    
    # Create Gaussian TF1
    xmin, xmax = (-0.01, 0.01) if x == 'xy' else (-0.02, 0.02)
    gaus_fit = TF1(f'gaus_fit_pt_{ix:.2f}', 'gaus', xmin, xmax)
    gaus_fit.SetParameters(h_dca.GetMaximum(), 0., 0.005)  # Initial parameters: amplitude, mean, sigma
    
    # Perform fit
    fit_result = h_dca.Fit(gaus_fit, 'RMS+', '', xmin, xmax)
    
    # Extract fit parameters
    mean = gaus_fit.GetParameter(1)
    mean_err = gaus_fit.GetParError(1)
    sigma = gaus_fit.GetParameter(2)
    sigma_err = gaus_fit.GetParError(2)
    
    fit_results = {
        'x': ix,
        'x_error': x_error,
        'mean': mean,
        'mean_err': mean_err,
        'sigma': sigma,
        'sigma_err': sigma_err
    }
    
    if fit_results_df is None:
        fit_results_df = pd.DataFrame.from_dict([fit_results])
    else:
        fit_results_df = pd.concat([fit_results_df, pd.DataFrame.from_dict([fit_results])], ignore_index=True)
    
    # Create and save canvas
    canvas = TCanvas(f'cDCA{x}_pt_{ix:.2f}', f'cDCA{x}_pt_{ix:.2f}', 800, 600)
    canvas.SetLogy()
    h_dca.SetMinimum(1)
    h_dca.Draw()
    gaus_fit.Draw('same')
    outfile.cd()
    canvas.Write()
    
    del h_dca, gaus_fit, canvas
    
    return fit_results_df

def simple_gaussian_calibration(h2, outfile: TFile, particle: str, x: str = 'xy'):
    """
    Perform simple Gaussian fits in the range [-0.015, 0.015] for calibration using TF1.
    """
    cfg = CONFIG[particle][x]
    
    x_min = cfg['x_min_fit']
    x_max = cfg['x_max_fit']
    
    fit_results_df = None
    
    x_bin_min = h2.GetXaxis().FindBin(x_min)
    x_bin_max = h2.GetXaxis().FindBin(x_max)
    
    for x_bin in range(x_bin_min, x_bin_max+1):
        # antimatter
        ix = h2.GetXaxis().GetBinCenter(x_bin)
        fit_results_df = _bin_calibration_simple_gaus(h2, outfile, particle, x_bin, fit_results_df, x)
        
        # matter
        positive_x_bin = h2.GetXaxis().FindBin(-ix)
        fit_results_df = _bin_calibration_simple_gaus(h2, outfile, particle, positive_x_bin, fit_results_df, x)
    
    if fit_results_df is None:
        print(f'No fit results for particle {particle}, skipping...')
        return
    
    visualize_fit_results(fit_results_df, particle, outfile, x)
    
    outfile.cd()
    h2.Write()
    
    del h2


def main_mc():
    input_files = ['/data/galucia/lithium_local/MC/LHC25a4_with_primaries_and_mixed.root',]
    tree_names = ['O2he3hadtable', 'O2he3hadtablemc',# 'O2he3hadmult'
                  ]
    folder_name = 'DF*'
    columns = ['fPtHe3', 'fPtHad', 'fDCAxyHe3', 'fDCAzHe3', 'fDCAxyHad', 'fDCAzHad', 'fFlags']
    datasets = []

    for tree_name in tree_names:
        datasets.append(Dataset.from_root(input_files, tree_name=tree_name, folder_name=folder_name, columns=columns))
    dataset = datasets[0].concat(datasets[1:], axis=1)

    read_bit = lambda x, bit: (x >> bit) & 0b1
    for ibit in range(4):
        dataset[FLAGS[ibit]] = dataset['fFlags'].apply(read_bit, args=(ibit,))
        dataset.add_subset(FLAGS[ibit], dataset[FLAGS[ibit]] == 1)

    print(f'{dataset["fFlags"].unique()}')
    print(f'{dataset.data[[*FLAGS]].describe()}')

    outfile = TFile.Open('output/dca_studies_mc.root', 'RECREATE')
    
    for particle in ['He3', 'Had']:
        outdir = outfile.mkdir(f'{particle}')
        h2_dcaxy, h2_dcaz = visualise(dataset, outdir, particle)
        fitdir = outdir.mkdir('fit_results')
        DCAxy_dir = fitdir.mkdir('DCAxy')
        calibration_routine(h2_dcaxy, DCAxy_dir, particle, x='xy')
        DCAz_dir = fitdir.mkdir('DCAz')
        calibration_routine(h2_dcaz, DCAz_dir, particle, x='z')

    outfile.Close()

def main_data():
    input_files = ['/data/galucia/lithium_local/same/LHC23_PbPb_pass4_hadronpid_same.root',
                   '/data/galucia/lithium_local/same/LHC24ar_pass1_hadronpid_same.root',
                   '/data/galucia/lithium_local/same/LHC24as_pass1_hadronpid_same.root',]
    #input_files = ['/data/galucia/lithium_local/same/LHC23_PbPb_pass5_same.root',
    #               '/data/galucia/lithium_local/same/LHC24_PbPb_pass2_same.root',
    #               #'/data/galucia/lithium_local/same/LHC24as_pass1_same.root',
    #               ]
    tree_names = ['O2he3hadtable'
                  ]
    folder_name = 'DF*'
    columns = ['fPtHe3', 'fPtHad', 'fDCAxyHe3', 'fDCAzHe3', 'fDCAxyHad', 'fDCAzHad']
    datasets = []

    for tree_name in tree_names:
        datasets.append(Dataset.from_root(input_files, tree_name=tree_name, folder_name=folder_name, columns=columns))
    dataset = datasets[0].concat(datasets[1:], axis=1)

    outfile = TFile.Open('output/dca_studies.root', 'RECREATE')
    
    for particle in ['He3', 'Had']:
        outdir = outfile.mkdir(f'{particle}')
        h2_dcaxy, h2_dcaz = visualise(dataset, outdir, particle, is_mc=False)
        fitdir = outdir.mkdir('fit_results')
        DCAxy_dir = fitdir.mkdir('DCAxy')
        calibration_routine(h2_dcaxy, DCAxy_dir, particle, x='xy')

        DCAz_dir = fitdir.mkdir('DCAz')
        calibration_routine(h2_dcaz, DCAz_dir, particle, x='z')

    outfile.Close()

def main_data_simple_gaus():
    """
    Alternative data processing routine using simple Gaussian fits.
    """
    input_files = ['/data/galucia/lithium_local/same/LHC23_PbPb_pass4_hadronpid_same.root',
                   '/data/galucia/lithium_local/same/LHC24ar_pass1_hadronpid_same.root',
                   '/data/galucia/lithium_local/same/LHC24as_pass1_hadronpid_same.root',]
    tree_names = ['O2he3hadtable']
    folder_name = 'DF*'
    columns = ['fPtHe3', 'fPtHad', 'fDCAxyHe3', 'fDCAzHe3', 'fDCAxyHad', 'fDCAzHad']
    datasets = []
    
    for tree_name in tree_names:
        datasets.append(Dataset.from_root(input_files, tree_name=tree_name, folder_name=folder_name, columns=columns))
    dataset = datasets[0].concat(datasets[1:], axis=1)

    gStyle.SetOptFit(1)
    
    outfile = TFile.Open('output/dca_studies_simple_gaus.root', 'RECREATE')
    
    for particle in ['He3', 'Had']:
        outdir = outfile.mkdir(f'{particle}')
        h2_dcaxy, h2_dcaz = visualise(dataset, outdir, particle, is_mc=False)
        fitdir = outdir.mkdir('fit_results')
        DCAxy_dir = fitdir.mkdir('DCAxy')
        simple_gaussian_calibration(h2_dcaxy, DCAxy_dir, particle, x='xy')
        
        DCAz_dir = fitdir.mkdir('DCAz')
        simple_gaussian_calibration(h2_dcaz, DCAz_dir, particle, x='z')
    
    outfile.Close()


if __name__ == '__main__':

    #main_data()
    #main_mc()
    main_data_simple_gaus()
