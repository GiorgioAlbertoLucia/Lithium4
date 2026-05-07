import sys
import yaml
import numpy as np
import ROOT
from ROOT import TFile, TChain, gInterpreter, RDataFrame, RooWorkspace, TH1F, TF1, RooMsgService, RooMinimizer, TH1, gROOT
from tqdm import tqdm

import time
import gc

from torchic.utils.terminal_colors import TerminalColors as tc
from torchic.utils.timeit import timeit
from torchic.core.histogram import load_hist, AxisSpec, HistLoadInfo
from torchic.utils.root import set_root_object

sys.path.append('..')
from utils.particles import ParticleMasses

from femto.core.signal_fitter import SignalFitter
from femto.core.bkg_fitter import BkgFitter
from femto.core.model_fitter import ModelFitter

gInterpreter.ProcessLine(f'#include "../include/Common.h"')
gInterpreter.ProcessLine(f'#include "../include/Variation.h"')
gInterpreter.ProcessLine(f'#include "../include/Systematics.h"')
from ROOT import ComputeAllSystematics

ROOT.EnableImplicitMT(10)
ROOT.gROOT.SetBatch(True)

selections = [
    '((fSignedPtHe3 > 0 && fSignedPtHad > 0) || (fSignedPtHe3 < 0 && fSignedPtHad < 0))',
    #'((fDeltaPhi/0.006)*(fDeltaPhi/0.006) + (fDeltaEta/0.012)*(fDeltaEta/0.012) > 1)',

    #'((fChi2TPCHe3 > 0.5) || (fIs23 == false))', # the cut on chi2tpc should be only applied to 2023
    '(fChi2TPCHe3 > 0.5)',
    '(fChi2TPCHe3 < 4)',
    '(fChi2TPCHad < 4)',
    '(std::abs(fEtaHe3) < 0.9)',
    '(std::abs(fEtaHad) < 0.9)',
    
    '(fNClsTPCHe3 > 110)',
    '((fPtHe3 < 2.5) || (fPIDtrkHe3 == 7))',
    '((fPIDtrkHe3 == 6) || (fPIDtrkHe3 == 7) || (fPIDtrkHe3 == 8))',

    '(std::abs(fNSigmaDCAxyHe3) < 3)',
    '(std::abs(fNSigmaDCAzHe3) < 3)',
    '(std::abs(fNSigmaDCAxyHad) < 3)',
    '(std::abs(fNSigmaDCAzHad) < 3)',
]
#selections = conf['selections']
base_selection = selections[0]
for sel in selections[1:]:
    base_selection += (' && ' + sel)
print(f'Selection: {base_selection}')

def find_good_df_folders(file_name: str, tree_name: str) -> list:
    import ROOT
    good_folders = []
    f = ROOT.TFile.Open(file_name)
    
    for key in f.GetListOfKeys():
        key_name = key.GetName()
        if 'DF_' not in key_name:
            continue
        
        tree = f.Get(f'{key_name}/{tree_name}')
        if not tree or tree.IsZombie():
            print(f'[WARN] Skipping zombie tree: {key_name}')
            continue

        tree.SetBranchStatus('*', 1)  # enable all branches
        n = tree.GetEntries()
        if n == 0:
            good_folders.append(key_name)
            continue

        is_good = True
        for entry in [0, n - 1]:
            result = tree.GetEntry(entry)
            if result < 0:  # GetEntry returns -1 on error
                print(f'[WARN] Skipping bad folder {key_name} at entry {entry}')
                is_good = False
                break
        
        if is_good:
            good_folders.append(key_name)
    
    f.Close()
    return good_folders

def prepare_input_tchain(config:dict):

    input_data = config['input_data']
    tree_names = config['tree_names']
    mode = config['mode']

    file_data_list = input_data if isinstance(input_data, list) else [input_data]
    tree_name = tree_names if isinstance(tree_names, str) else tree_names[0]
    chain_data = TChain('tchain')

    additional_chains = []
    if isinstance(tree_names, list) and len(tree_names) > 1:
        for idx, tname in enumerate(tree_names[1:]):
            additional_chain = TChain(f'tchain_friend_{idx}')
            additional_chains.append(additional_chain)

    for file_name in file_data_list:
        fileData = TFile(file_name)

        if mode == 'DF':
            for key in fileData.GetListOfKeys():
                key_name = key.GetName()
                if 'DF_' in key_name :
                    chain_data.Add(f'{file_name}/{key_name}/{tree_name}')
                    for idx, additional_chain in enumerate(additional_chains):
                        additional_chain.Add(f'{file_name}/{key_name}/{tree_names[idx+1]}')
            
        elif mode == 'tree':
            print(f'Adding {tc.CYAN+tc.UNDERLINE}{file_name}/{tree_name}{tc.RESET} to the chain')
            chain_data.Add(f'{file_name}/{tree_name}')
            for idx, additional_chain in enumerate(additional_chains):
                additional_chain.Add(f'{file_name}/{tree_names[idx+1]}')

    for idx, additional_chain in enumerate(additional_chains):
        chain_data.AddFriend(additional_chain)

    return chain_data, additional_chains

def prepare_rdataframe(chain_data: TChain, base_selection: str, selection: str):
   
    rdf = RDataFrame(chain_data)
    print(tc.GREEN+'\nDataset columns'+tc.RESET)
    print(tc.UNDERLINE+tc.CYAN+f'{rdf.GetColumnNames()}'+tc.RESET)
    
    # TPC
    if 'fNSigmaTPCHadPr' in rdf.GetColumnNames():
        rdf = rdf.Define('fNSigmaTPCHad', 'fNSigmaTPCHadPr')
    
    # TOF
    if 'fNSigmaTOFHadPr' not in rdf.GetColumnNames():
        rdf = rdf.Define('fNSigmaTOFHad', 'ComputeNsigmaTOFPr(std::abs(fPtHad), fMassTOFHad)')
    else:
       rdf = rdf.Define('fNSigmaTOFHad', 'fNSigmaTOFHadPr')
    
      # Recalibration
      # Correct for PID in tracking
      
    rdf = rdf.Define('fSignedPtHad', 'fPtHad') \
      .Define('fSignHe3', 'fPtHe3/std::abs(fPtHe3)') \
      .Redefine('fPtHe3', 'std::abs(fPtHe3)') \
      .Redefine('fPtHad', 'std::abs(fPtHad)') \
      .Redefine('fPtHe3', '(fPIDtrkHe3 == 7) || (fPIDtrkHe3 == 8) || (fPtHe3 > 2.5) ? fPtHe3 : CorrectPidTrkHe(fPtHe3)') \
      .Define('fSignedPtHe3', 'fPtHe3 * fSignHe3') \
      .Define(f'fEHe3', f'std::sqrt((fPtHe3 * std::cosh(fEtaHe3))*(fPtHe3 * std::cosh(fEtaHe3)) + {ParticleMasses["He"]}*{ParticleMasses["He"]})') \
      .Define(f'fEHad', f'std::sqrt((fPtHad * std::cosh(fEtaHad))*(fPtHad * std::cosh(fEtaHad)) + {ParticleMasses["Pr"]}*{ParticleMasses["Pr"]})') \
      .Define('fDeltaEta', 'fEtaHe3 - fEtaHad') \
      .Define('fDeltaPhi', 'fPhiHe3 - fPhiHad') \
      .Redefine('fInnerParamTPCHe3', 'fInnerParamTPCHe3 * 2') \
      .Define('fClusterSizeCosLamHe3', 'ComputeAverageClusterSize(fItsClusterSizeHe3) / cosh(fEtaHe3)') \
      .Define('fClusterSizeCosLamHad', 'ComputeAverageClusterSize(fItsClusterSizeHad) / cosh(fEtaHad)') \
      .Define('fNSigmaITSHe3', 'ComputeNsigmaITSHe(fPtHe3 * std::cosh(fEtaHe3), fClusterSizeCosLamHe3)') \
      .Define('fNSigmaITSHad', 'ComputeNsigmaITSPr(fPtHad * std::cosh(fEtaHad), fClusterSizeCosLamHad)') \
      .Redefine('fNSigmaTPCHe3', 'ComputeNsigmaTPCHe(std::abs(fInnerParamTPCHe3), fSignalTPCHe3)') \
      .Define('fNSigmaDCAxyHe3', 'ComputeNsigmaDCAxyHe(fPtHe3, fDCAxyHe3)') \
      .Define('fNSigmaDCAzHe3', 'ComputeNsigmaDCAzHe(fPtHe3, fDCAzHe3)') \
      .Define('fNSigmaDCAxyHad', 'ComputeNsigmaDCAxyPr(fPtHad, fDCAxyHad)') \
      .Define('fNSigmaDCAzHad', 'ComputeNsigmaDCAzPr(fPtHad, fDCAzHad)') \
      .Filter(base_selection).Filter(selection) \
      .Filter('(fCentralityFT0C < 10)') \
      .Define('fKstar', f'ComputeKstar(fPtHe3, fEtaHe3, fPhiHe3, {ParticleMasses["He"]}, fPtHad, fEtaHad, fPhiHad, {ParticleMasses["Pr"]})')
    
    return rdf

def load_same():
    
    config_same = {
        'input_data':   [ 
                #'/data/galucia/lithium_local/same_merged/LHC23_PbPb_pass5_same.root',
                #'/data/galucia/lithium_local/same_merged/LHC24_PbPb_pass2_same.root',
                '/data/galucia/lithium_local/same/LHC23_PbPb_pass4_hadronpid_same.root',
                '/data/galucia/lithium_local/same/LHC24ar_pass1_hadronpid_same.root',
                '/data/galucia/lithium_local/same/LHC24as_pass1_hadronpid_same.root',
              ],
        'tree_names':    ['O2he3hadtable', 'O2he3hadmult'],
        'mode':         'DF',
    }
    chain_data_same, additional_chains = prepare_input_tchain(config_same)
    rdf_same = prepare_rdataframe(chain_data_same, base_selection, 'true')
    
    return rdf_same, chain_data_same, additional_chains

def load_mixed():
    
    config_mixed = {
        'input_data':   [ 
                '/data/galucia/lithium_local/mixing/LHC23_PbPb_pass4_hadronpid_event_mixing_batch256_reject_multiples.root',
                '/data/galucia/lithium_local/mixing/LHC24ar_pass1_hadronpid_event_mixing_batch256_reject_multiples.root',
                '/data/galucia/lithium_local/mixing/LHC24as_pass1_hadronpid_event_mixing_batch256_reject_multiples.root',
                #'/data/galucia/lithium_local/mixing/LHC23_PbPb_pass4_hadronpid_event_mixing_batch105_refined_dca.root',
                #'/data/galucia/lithium_local/mixing/LHC24ar_pass1_hadronpid_event_mixing_batch105_refined_dca.root',
                #'/data/galucia/lithium_local/mixing/LHC24as_pass1_hadronpid_event_mixing_batch105_refined_dca.root',
                #'/data/galucia/lithium_local/mixing/LHC23_PbPb_pass4_hadronpid_event_mixing_batch42_refined_dca.root',
                #'/data/galucia/lithium_local/mixing/LHC24ar_pass1_hadronpid_event_mixing_batch42_refined_dca.root',
                #'/data/galucia/lithium_local/mixing/LHC24as_pass1_hadronpid_event_mixing_batch42_refined_dca.root',
                #'/data/galucia/lithium_local/mixing/LHC23_PbPb_pass4_hadronpid_event_mixing_batch256_refined_dca.root',
                #'/data/galucia/lithium_local/mixing/LHC24ar_pass1_hadronpid_event_mixing_batch256_refined_dca.root',
                #'/data/galucia/lithium_local/mixing/LHC24as_pass1_hadronpid_event_mixing_batch256_refined_dca.root',
                #'/data/galucia/lithium_local/mixing/LHC23_PbPb_pass4_hadronpid_event_mixing_batch3112_refined_dca.root',
                #'/data/galucia/lithium_local/mixing/LHC24ar_pass1_hadronpid_event_mixing_batch3112_refined_dca.root',
                #'/data/galucia/lithium_local/mixing/LHC24as_pass1_hadronpid_event_mixing_batch3112_refined_dca.root',
              ],
        'tree_names':    'MixedTree',
        'mode':         'tree',
    }
    chain_data_mixed, additional_chains = prepare_input_tchain(config_mixed)
    rdf_mixed = prepare_rdataframe(chain_data_mixed, base_selection, 'true')
    
    return rdf_mixed, chain_data_mixed, additional_chains

def poisson_sampling(hist:TH1F):

    for ibin in range(1, hist.GetNbinsX()+1):
        content = hist.GetBinContent(ibin)

        new_content = np.random.poisson(content)

        hist.SetBinContent(ibin, new_content)
        hist.SetBinError(ibin, np.sqrt(new_content))

    return hist
    

@timeit
def prepare_systematics_histograms(rdf:RDataFrame, n_variations:int, hist_name_suffix:str='', variable:str=''):

    variational_rdf = None
    if variable == 'all':
        variational_rdf = rdf.Define('fNSigmaITSHe3Systematic', 'fNSigmaITSHe3') \
                         .Define('fAbsNSigmaTPCHe3Systematic', 'std::abs(fNSigmaTPCHe3)') \
                         .Define('fAbsNSigmaDCAxyHe3Systematic', 'std::abs(fNSigmaDCAxyHe3)') \
                         .Define('fAbsNSigmaDCAzHe3Systematic', 'std::abs(fNSigmaDCAzHe3)') \
                         .Define('fNClsTPCHe3Systematic', 'static_cast<float>(fNClsTPCHe3)') \
                         .Define('fChi2TPCHe3Systematic', 'fChi2TPCHe3') \
                         .Define('fNSigmaITSHadSystematic', 'fNSigmaITSHad') \
                         .Define('fAbsNSigmaTPCHadSystematic', 'std::abs(fNSigmaTPCHad)') \
                         .Define('fAbsNSigmaTOFHadSystematic', 'std::abs(fNSigmaTOFHad)') \
                         .Define('fAbsNSigmaDCAxyHadSystematic', 'std::abs(fNSigmaDCAxyHad)') \
                         .Define('fAbsNSigmaDCAzHadSystematic', 'std::abs(fNSigmaDCAzHad)') \
                         .Define('fChi2TPCHadSystematic', 'fChi2TPCHad') \
                         .Define('fAbsZVertexSystematic', 'std::abs(fZVertex)') \
                         .Vary(['fNSigmaITSHe3Systematic', 'fAbsNSigmaTPCHe3Systematic',
                                'fAbsNSigmaDCAxyHe3Systematic', 'fAbsNSigmaDCAzHe3Systematic',
                                'fNClsTPCHe3Systematic', 'fChi2TPCHe3Systematic',
                                'fNSigmaITSHadSystematic', 'fAbsNSigmaTPCHadSystematic', 'fAbsNSigmaTOFHadSystematic',
                                'fAbsNSigmaDCAxyHadSystematic', 'fAbsNSigmaDCAzHadSystematic', 'fAbsZVertexSystematic'],
                                f'systematicCuts({n_variations}, fNSigmaITSHe3Systematic, fAbsNSigmaTPCHe3Systematic, \
                                fAbsNSigmaDCAxyHe3Systematic, fAbsNSigmaDCAzHe3Systematic, fNClsTPCHe3Systematic, fChi2TPCHe3Systematic, \
                                fNSigmaITSHadSystematic, fAbsNSigmaTPCHadSystematic, fAbsNSigmaTOFHadSystematic, \
                                fAbsNSigmaDCAxyHadSystematic, fAbsNSigmaDCAzHadSystematic, fChi2TPCHadSystematic, \
                                fAbsZVertexSystematic)',
                                n_variations, 'systematic') \
                         .Filter('(fNSigmaITSHe3Systematic > 0) && \
                                  (fAbsNSigmaTPCHe3Systematic < 0) && \
                                  (fNSigmaITSHadSystematic > 0) && \
                                  (fAbsNSigmaTPCHadSystematic < 0) && \
                                  ((std::abs(fPtHad) < 0.8) || (fAbsNSigmaTOFHadSystematic < 0))')
    else:
        variational_rdf = rdf.Define('fNSigmaITSHe3Systematic', 'fNSigmaITSHe3') \
                            .Define('fAbsNSigmaTPCHe3Systematic', 'std::abs(fNSigmaTPCHe3)') \
                            .Define('fAbsNSigmaDCAxyHe3Systematic', 'std::abs(fNSigmaDCAxyHe3)') \
                            .Define('fAbsNSigmaDCAzHe3Systematic', 'std::abs(fNSigmaDCAzHe3)') \
                            .Define('fNClsTPCHe3Systematic', 'static_cast<float>(fNClsTPCHe3)') \
                            .Define('fChi2TPCHe3Systematic', 'fChi2TPCHe3') \
                            .Define('fNSigmaITSHadSystematic', 'fNSigmaITSHad') \
                            .Define('fAbsNSigmaTPCHadSystematic', 'std::abs(fNSigmaTPCHad)') \
                            .Define('fAbsNSigmaTOFHadSystematic', 'std::abs(fNSigmaTOFHad)') \
                            .Define('fAbsNSigmaDCAxyHadSystematic', 'std::abs(fNSigmaDCAxyHad)') \
                            .Define('fAbsNSigmaDCAzHadSystematic', 'std::abs(fNSigmaDCAzHad)') \
                            .Define('fChi2TPCHadSystematic', 'fChi2TPCHad') \
                            .Define('fAbsZVertexSystematic', 'std::abs(fZVertex)') \
                            .Vary(['fNSigmaITSHe3Systematic', 'fAbsNSigmaTPCHe3Systematic',
                                    'fAbsNSigmaDCAxyHe3Systematic', 'fAbsNSigmaDCAzHe3Systematic',
                                    'fNClsTPCHe3Systematic', 'fChi2TPCHe3Systematic',
                                    'fNSigmaITSHadSystematic', 'fAbsNSigmaTPCHadSystematic', 'fAbsNSigmaTOFHadSystematic',
                                    'fAbsNSigmaDCAxyHadSystematic', 'fAbsNSigmaDCAzHadSystematic', 'fAbsZVertexSystematic'],
                                    f'systematicCutsSingleVariable({n_variations}, "{variable}", fNSigmaITSHe3Systematic, fAbsNSigmaTPCHe3Systematic, \
                                    fAbsNSigmaDCAxyHe3Systematic, fAbsNSigmaDCAzHe3Systematic, fNClsTPCHe3Systematic, fChi2TPCHe3Systematic, \
                                    fNSigmaITSHadSystematic, fAbsNSigmaTPCHadSystematic, fAbsNSigmaTOFHadSystematic, \
                                    fAbsNSigmaDCAxyHadSystematic, fAbsNSigmaDCAzHadSystematic, fChi2TPCHadSystematic, \
                                    fAbsZVertexSystematic)',
                                    n_variations, 'systematic') \
                            .Filter('(fNSigmaITSHe3Systematic > 0) && \
                                  (fAbsNSigmaTPCHe3Systematic < 0) && \
                                  (fNSigmaITSHadSystematic > 0) && \
                                  (fAbsNSigmaTPCHadSystematic < 0) && \
                                  ((std::abs(fPtHad) < 0.8) || (fAbsNSigmaTOFHadSystematic < 0)) ')

    nominal_hist_010 = variational_rdf.Filter('fCentralityFT0C < 10') \
                                      .Histo1D((f'hKstar010{hist_name_suffix}', ';#it{k}* (GeV/#it{c});', 20, 0, 0.4), 'fKstar')
    variational_hists_010 = ROOT.RDF.Experimental.VariationsFor(nominal_hist_010)

    #nominal_hist_1030 = variational_rdf.Filter('fCentralityFT0C >= 10 && fCentralityFT0C < 30') \
    #                                  .Histo1D((f'hKstar1030{hist_name_suffix}', ';#it{k}* (GeV/#it{c});', 20, 0, 0.4), 'fKstar')
    #variational_hists_1030 = ROOT.RDF.Experimental.VariationsFor(nominal_hist_1030)
    #
    #nominal_hist_3050 = variational_rdf.Filter('fCentralityFT0C >= 30 && fCentralityFT0C < 50') \
    #                                  .Histo1D((f'hKstar3050{hist_name_suffix}', ';#it{k}* (GeV/#it{c});', 20, 0, 0.4), 'fKstar')
    #variational_hists_3050 = ROOT.RDF.Experimental.VariationsFor(nominal_hist_3050)

    return variational_hists_010 #, variational_hists_1030, variational_hists_3050

def normalise_histogram(h_same, h_mixed, NORM_LOW_KSTAR, NORM_HIGH_KSTAR):

    low_bin = h_same.FindBin(NORM_LOW_KSTAR)
    high_bin = h_same.FindBin(NORM_HIGH_KSTAR)
    normalization_factor = h_same.Integral(low_bin, high_bin) / h_mixed.Integral(low_bin, high_bin)

    h_mixed_normalised = h_mixed.Clone(f'{h_mixed.GetName()}_Normalised')
    h_mixed_normalised.Scale(normalization_factor)

    return h_mixed_normalised

def correlation_function_centrality_integrated(h_sames, h_mixeds, suffix:str):
    """
    Compute the correlation function for the centrality integrated histograms.
    """

    h_same = h_sames[0].Clone()
    h_mixed = h_mixeds[0].Clone()

    for h_same_centrality, h_mixed_centrality in zip(h_sames[1:], h_mixeds[1:]):
        h_same.Add(h_same_centrality)
        h_mixed.Add(h_mixed_centrality)

    h_corr = h_same.Clone(f'hCorrelation{suffix}')
    h_corr.Divide(h_mixed)

    return h_corr

def fitting_routine(outfile, bkg_input_path:str, h_bkg_name:str,
                    h_data:TH1F, h_mixed_event:TH1F, 
                    sign:str, centrality:str, use_smoothening:bool=False):

    #SIGNAL_HIST_LOAD_INFO = HistLoadInfo('models/li4_contribution.root', 'hCkHist')
    #SIGNAL_HIST_LOAD_INFO = HistLoadInfo('/home/galucia/Lithium4/femto/models/li4_contribution_finer_binning.root', 'hCkHist')
    SIGNAL_HIST_LOAD_INFO = HistLoadInfo('/home/galucia/Lithium4/femto/models/li4_contribution_proper_sill.root', 'hCkHist')
    
    workspace = RooWorkspace('roows')

    h_bkg = load_hist(bkg_input_path, h_bkg_name)
    h_signal = load_hist(SIGNAL_HIST_LOAD_INFO)
    h_correlation_function = h_data

    KSTAR_MIN, KSTAR_MAX = ((0.01, 0.4) if ('010' not in centrality or sign != 'Antimatter' ) 
                            else (0.02, 0.4))
    kstar_spec = AxisSpec(100, KSTAR_MIN, KSTAR_MAX, 'kstar', '#it{k}* (GeV/#it{c})')
    
    signal_fitter = SignalFitter('signal', kstar_spec, outfile, workspace)
    signal_init_mode = 'from_mc' if not use_smoothening else 'from_kde'
    signal_fitter.init_signal(signal_init_mode, h_signal)
    signal_fitter.title = '^{4}Li'
    signal_fitter.save_to_workspace()

    bkg_fitter = BkgFitter('bkg', kstar_spec, outfile, workspace)
    bkg_init_mode = 'from_mc' if not use_smoothening else 'from_kde'
    bkg_fitter.init_bkg(bkg_init_mode, h_bkg, rho=(0.05 if '010' not in centrality else 0.1)) #, extended=True)
    bkg_fitter.title = 'Full model' 
    bkg_fitter.save_to_workspace()

    model_fitter = ModelFitter('model', kstar_spec, outfile, ['signal_pdf'], ['bkg_pdf'], workspace, extended=True)
    #model_fitter.fractions['signal_pdf'].setRange(0., 1.)
    model_fitter.load_data(h_correlation_function, h_correlation_function.GetName())
    model_fitter.prefit_background(h_correlation_function, range_limits=(0.2, 0.4), range_name='bkg_fit_range',
                                   save_normalisation_value=True) #, use_chi2_method=False)
    model_fitter.fractions['signal_pdf'].setVal(0.3)
    model_fitter.fractions['signal_pdf'].setRange(-1e4, 1e4)
    model_fitter.fit_model(h_correlation_function, signal_name='signal_pdf', norm_range='bkg_fit_range')
    model_fitter.save_to_workspace()
    model_fitter.compute_chi2(h_correlation_function)
    xvar = workspace.obj(kstar_spec.name)
    xvar.setRange(KSTAR_MIN, 0.39)
    raw_yield_value = model_fitter.compute_raw_yield(h_mixed_event, 'signal_pdf', 'bkg_pdf')

    del workspace, signal_fitter, bkg_fitter, model_fitter

    return raw_yield_value



def prepare_histograms():

    rdf_same, _chain_same, _chain_additional_same = load_same()
    rdf_mixed, _chain_mixed, _chain_additional_mixed = load_mixed()

    N_ITERATIONS = 100
    outFile = TFile("output/hist_systematics_with_upper_limit_010.root", "RECREATE")

    for sign, condition in {'Matter': 'fSignedPtHe3 > 0', 'Antimatter': 'fSignedPtHe3 < 0'}.items():

        outDir = outFile.mkdir(sign)

        tmp_rdf_same = rdf_same.Filter(condition)
        tmp_rdf_mixed = rdf_mixed.Filter(condition)

        variational_hists_same = prepare_systematics_histograms(tmp_rdf_same, n_variations=1, hist_name_suffix='Same', variable='')
        variational_hists_mixed = prepare_systematics_histograms(tmp_rdf_mixed, n_variations=1, hist_name_suffix='Mixed', variable='')

        NORM_LOW_KSTAR, NORM_HIGH_KSTAR = 0.2, 0.4

        #same_010, same_1030, same_3050 = variational_hists_same[0][f'systematic:{0}'], \
        #                                variational_hists_same[1][f'systematic:{0}'], variational_hists_same[2][f'systematic:{0}']
        #mixed_010, mixed_1030, mixed_3050 = variational_hists_mixed[0][f'systematic:{0}'], \
        #                                variational_hists_mixed[1][f'systematic:{0}'], variational_hists_mixed[2][f'systematic:{0}']

        same_010 = variational_hists_same[f'systematic:{0}']
        mixed_010 = variational_hists_mixed[f'systematic:{0}']

        mixed_normalised_010 = normalise_histogram(same_010, mixed_010, NORM_LOW_KSTAR, NORM_HIGH_KSTAR)
        #mixed_normalised_1030 = normalise_histogram(same_1030, mixed_1030, NORM_LOW_KSTAR, NORM_HIGH_KSTAR)
        #mixed_normalised_3050 = normalise_histogram(same_3050, mixed_3050, NORM_LOW_KSTAR, NORM_HIGH_KSTAR)

        #h_correlation_ref = correlation_function_centrality_integrated([same_010, same_1030, same_3050],
        #                                                                [mixed_normalised_010, mixed_normalised_1030, mixed_normalised_3050],
        #                                                                0)
        outDir.cd()
        #h_correlation_ref.Write('hCorrelationReference')

        gc.collect()

        variational_hists_same = prepare_systematics_histograms(tmp_rdf_same, n_variations=N_ITERATIONS, hist_name_suffix='Same', variable='all')
        variational_hists_mixed = prepare_systematics_histograms(tmp_rdf_mixed, n_variations=N_ITERATIONS, hist_name_suffix='Mixed', variable='all')

        h_sames, h_mixeds, h_mixeds_normalised = [], [], []
        NORM_LOW_KSTAR, NORM_HIGH_KSTAR = 0.2, 0.4

        for iter in range(N_ITERATIONS):

            print(f'Processing iteration {iter+1}/{N_ITERATIONS} for {sign}...')

            #same_010, same_1030, same_3050 = variational_hists_same[0][f'systematic:{iter}'], \
            #                                variational_hists_same[1][f'systematic:{iter}'], variational_hists_same[2][f'systematic:{iter}']
            #mixed_010, mixed_1030, mixed_3050 = variational_hists_mixed[0][f'systematic:{iter}'], \
            #                                variational_hists_mixed[1][f'systematic:{iter}'], variational_hists_mixed[2][f'systematic:{iter}']
            same_010 = variational_hists_same[f'systematic:{iter}']
            mixed_010 = variational_hists_mixed[f'systematic:{iter}']

            mixed_normalised_010 = normalise_histogram(same_010, mixed_010, NORM_LOW_KSTAR, NORM_HIGH_KSTAR)
            #mixed_normalised_1030 = normalise_histogram(same_1030, mixed_1030, NORM_LOW_KSTAR, NORM_HIGH_KSTAR)
            #mixed_normalised_3050 = normalise_histogram(same_3050, mixed_3050, NORM_LOW_KSTAR, NORM_HIGH_KSTAR)

            #h_correlation_iter = correlation_function_centrality_integrated([same_010, same_1030, same_3050],
            #                                                                [mixed_normalised_010, mixed_normalised_1030, mixed_normalised_3050],
            #                                                                iter)

            h_sames.append({'010': same_010})                           #, '1030': same_1030, '3050': same_3050})
            h_mixeds.append({'010': mixed_010,})                        # '1030': mixed_1030, '3050': mixed_3050})
            h_mixeds_normalised.append({'010': mixed_normalised_010})   #, '1030': mixed_normalised_1030, '3050': mixed_normalised_3050})

        for iter, (sames, mixeds, mixeds_normalised) in enumerate(zip(h_sames, h_mixeds, h_mixeds_normalised)):
            outDirIter = outDir.mkdir(f'iter_{iter}')
            outDirIter.cd()
            for centrality in ['010']: #, '1030', '3050']:
                sames[centrality].Write(f'hSame_{centrality}')
                mixeds[centrality].Write(f'hMixed_{centrality}')
                mixeds_normalised[centrality].Write(f'hMixedNormalised_{centrality}')
    
    outFile.Close()

def upper_limit_systematic_routine():

    N_ITERATIONS = 20
    N_UPPER_LIMIT_ITERATIONS = 500
    TH1.AddDirectory(False) # Prevent ROOT from automatically associating histograms to the current directory, which can lead to memory issues

    infile = TFile.Open("output/hist_systematics_with_upper_limit_010.root")
    outFile = TFile.Open("output/systematics_with_upper_limit_010_v3.root", "RECREATE")

    available_bkgs = ['nominal/hLambdaSigmaCorrectedCk_matched_for_smearing_matched_smeared', 
                      'nominal/hLambdaSigmaCorrectedCk_Higher_matched_for_smearing_matched_smeared', 
                      'nominal/hLambdaSigmaCorrectedCk_Lower_matched_for_smearing_matched_smeared', 
                      'upper/hLambdaSigmaCorrectedCk_matched_for_smearing_matched_smeared',
                      'upper/hLambdaSigmaCorrectedCk_Higher_matched_for_smearing_matched_smeared', 
                      'upper/hLambdaSigmaCorrectedCk_Lower_matched_for_smearing_matched_smeared', 
                      'lower/hLambdaSigmaCorrectedCk_matched_for_smearing_matched_smeared', 
                      'lower/hLambdaSigmaCorrectedCk_Higher_matched_for_smearing_matched_smeared',
                      'lower/hLambdaSigmaCorrectedCk_Lower_matched_for_smearing_matched_smeared',]

    #for sign in ['Matter', 'Antimatter']:
    for sign in ['Matter']:

        inDir = infile.Get(sign)
        outDir = outFile.mkdir(sign)

        h_upper_limits_010 = TH1F(f'hUpperLimits_010', ';[#it{N}_{^{4}Li}^{raw}]_{upper limit};', 600, 0, 600)

        for iter in tqdm(range(N_ITERATIONS)):

            outDirIter = outDir.mkdir(f'iter_{iter}')

            h_same_010 = inDir.Get(f'iter_{iter}/hSame_010')
            h_mixed_010 = inDir.Get(f'iter_{iter}/hMixed_010')
            h_mixed_normalised_010 = inDir.Get(f'iter_{iter}/hMixedNormalised_010')
            
            h_correlation_010 = h_same_010.Clone(f'hCorrelation_010_{iter}')
            h_correlation_010.Divide(h_mixed_normalised_010)
            
            outDirNominal = outDirIter.mkdir(f'nominal')

            random_bkg_name = np.random.choice(available_bkgs)

            nominal_raw_yield = fitting_routine(outDirNominal, 
                                bkg_input_path='/home/galucia/Lithium4/femto/models/lambda_models.root', 
                                h_bkg_name=f'{sign}/010/{random_bkg_name}',
                                h_data=h_correlation_010, h_mixed_event=h_mixed_normalised_010, 
                                sign=sign, centrality='010', use_smoothening=True)

            h_raw_yields_010_iter = TH1F(f'hRawYields_010_Iter_{iter}', ';Raw yield;', 300, -200, 400)

            for iter_upper in range(N_UPPER_LIMIT_ITERATIONS):

                h_same_010_iter_upper = h_same_010.Clone(f'hSame_010_Iter_Upper_{iter_upper}')
                h_mixed_normalised_010_iter_upper = h_mixed_normalised_010.Clone(f'hMixed_010_Iter_Upper_{iter_upper}')

                h_same_010_iter_upper = poisson_sampling(h_same_010_iter_upper)
                h_mixed_normalised_010_iter_upper = poisson_sampling(h_mixed_normalised_010_iter_upper)
                h_correlation_010_iter_upper = h_same_010_iter_upper.Clone(f'hCorrelation_010_Iter_Upper_{iter_upper}')
                h_correlation_010_iter_upper.Divide(h_mixed_normalised_010_iter_upper)

                outDir_to_pass = outDirIter.mkdir(f'inner_iter_{iter_upper}') if iter_upper == 0 else None
                raw_yield = fitting_routine(outDir_to_pass, 
                                bkg_input_path='/home/galucia/Lithium4/femto/models/lambda_models.root', 
                                h_bkg_name=f'{sign}/010/hLambdaSigmaCorrectedCk',
                                h_data=h_correlation_010_iter_upper, h_mixed_event=h_mixed_normalised_010_iter_upper, 
                                sign=sign, centrality='010', use_smoothening=True)

                h_raw_yields_010_iter.Fill(raw_yield)

                #outDir.cd()
                #h_same_010_iter_upper.Write()
                #h_mixed_normalised_010_iter_upper.Write()

                del h_same_010_iter_upper, h_mixed_normalised_010_iter_upper, h_correlation_010_iter_upper


            fit_func_010 = TF1(f'fit_func_010_{iter}', 'gaus', -200, 400)
            fit_func_010.SetParameters(h_raw_yields_010_iter.GetMaximum(), h_raw_yields_010_iter.GetMean(), h_raw_yields_010_iter.GetRMS())
            h_raw_yields_010_iter.Fit('gaus', 'RMS+', '', -200, 400)

            outDirIter.cd()
            h_raw_yields_010_iter.Write()
            
            upper_limit_010 = nominal_raw_yield + 1.96 * fit_func_010.GetParameter(2) # 95% confidence level upper limit
            h_upper_limits_010.Fill(upper_limit_010)

            #gROOT.GetListOfCleanups().Clear()
            gc.collect()

            del h_same_010, h_mixed_010, h_mixed_normalised_010, h_correlation_010, h_raw_yields_010_iter, fit_func_010
        
        outDir.cd()
        h_upper_limits_010.Write()

        del h_upper_limits_010

    infile.Close()
    outFile.Close()


if __name__ == "__main__":

    RooMsgService.instance().setGlobalKillBelow(5) # 3 = WARNING, 4 = ERROR, 5 = FATAL
    ROOT.RooFit.PrintLevel(-1)

    #prepare_histograms()

    upper_limit_systematic_routine()

    print("Systematics analysis completed successfully.")