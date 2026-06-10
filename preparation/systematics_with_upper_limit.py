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
from include.load_parameters import load_parametrisation
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

PARAMETRISATIONS = {
    '2023': {
        'kDCAxyResolutionParams_Pr': [8.19e-5, 0.004, 0.0026, 1.1741],
        'kDCAxyResolutionParams_He': [1.09e-4, 0.0011, 0.0065, 1.0399],
        'kDCAzResolutionParams_Pr': [1.18e-4, 0.0020, 0.0025, 1.3460],
        'kDCAzResolutionParams_He': [9.36e-5, 0.0019, 0.0080, 1.416],

        # PbPb 23 pass5
        'kHeTPCParams': [-244.627, 0.100002, 1.14056, 0.365353, 3.26339],
        'kHeTPCResolution': 0.061,
        'kHeTPCResidualType': 'double_gaus',
        'kHeTPCParamsResiduals': [87.1178, 0.438363, 0.0814248, -16.9173, 0.754105, 0.285692],

        'kITSParams_Pr': [1.0228, 1.9634, 2.2081], # Pr from 'Ka' fitting
        'kITSParams_He': [2.5734, 1.3620, 5.1383],  # He 
        'kITSResolutionParams_Pr': [0.1575, 0., 0.],
        'kITSResolutionParams_He': [0.1124, -0.0100, 0.],

        'kHePidTrkParamsPt': [0.3101, -0.1759, 0.0262],
        'kHePidTrkParamsP': [0., 0., 0.],
    }, 
    '2024': {
        'kDCAxyResolutionParams_Pr': [8.19e-5, 0.004, 0.0026, 1.1741],
        'kDCAxyResolutionParams_He': [1.09e-4, 0.0011, 0.0065, 1.0399],
        'kDCAzResolutionParams_Pr': [1.18e-4, 0.0020, 0.0025, 1.3460],
        'kDCAzResolutionParams_He': [9.36e-5, 0.0019, 0.0080, 1.416],

        # 24_pass3
        'kHeTPCParams': [-190.723, 0.100000, 1.45622, 0.944291, 2.96807],
        'kHeTPCResolution': 0.059,
        'kHeTPCResidualType': 'exp',
        'kHeTPCParamsResiduals': [-19.1089, 5.27807, 0.698301, 0, 0, 0],

        'kITSParams_Pr': [1.0228, 1.9634, 2.2081], # Pr from 'Ka' fitting
        'kITSParams_He': [2.5734, 1.3620, 5.1383],  # He 
        'kITSResolutionParams_Pr': [0.1575, 0., 0.],
        'kITSResolutionParams_He': [0.1124, -0.0100, 0.],

        'kHePidTrkParamsPt': [0.3101, -0.1759, 0.0262],
        'kHePidTrkParamsP': [0., 0., 0.],
    }, 
    '2025': {
        'kDCAxyResolutionParams_Pr': [8.19e-5, 0.004, 0.0026, 1.1741],
        'kDCAxyResolutionParams_He': [1.09e-4, 0.0011, 0.0065, 1.0399],
        'kDCAzResolutionParams_Pr': [1.18e-4, 0.0020, 0.0025, 1.3460],
        'kDCAzResolutionParams_He': [9.36e-5, 0.0019, 0.0080, 1.416],

        # 25_pass1
        'kHeTPCParams': [-208.928, 0.100000, 1.29761, 0.756458, 3.04256],
        'kHeTPCResolution': 0.059,
        'kHeTPCResidualType': 'double_gaus',
        'kHeTPCParamsResiduals': [200.000, 0.362588, 0.0235587, -32.5300, 0.400000, 0.247265],

        'kITSParams_Pr': [1.0228, 1.9634, 2.2081], # Pr from 'Ka' fitting
        'kITSParams_He': [2.5734, 1.3620, 5.1383],  # He 
        'kITSResolutionParams_Pr': [0.1575, 0., 0.],
        'kITSResolutionParams_He': [0.1124, -0.0100, 0.],

        'kHePidTrkParamsPt': [0.3101, -0.1759, 0.0262],
        'kHePidTrkParamsP': [0., 0., 0.],
    }, 
        
    }

selections = [
    '((fSignedPtHe3 > 0 && fSignedPtHad > 0) || (fSignedPtHe3 < 0 && fSignedPtHad < 0))',
    #'((fDeltaPhi/0.006)*(fDeltaPhi/0.006) + (fDeltaEta/0.012)*(fDeltaEta/0.012) > 1)',

    '(fChi2TPCHe3 < 4)',
    '(fChi2TPCHad < 4)',
    '(std::abs(fEtaHe3) < 0.9)',
    '(std::abs(fEtaHad) < 0.9)',
    
    '(fNClsTPCHe3 > 110)',
    '((fPtHe3 < 2.5) || (fPIDtrkHe3 == 7))',
    '((fPIDtrkHe3 == 6) || (fPIDtrkHe3 == 7) || (fPIDtrkHe3 == 8))',

    #'(std::abs(fNSigmaDCAxyHe3) < 5)',
    #'(std::abs(fNSigmaDCAzHe3) < 5)',
    #'(std::abs(fNSigmaDCAxyHad) < 5)',
    #'(std::abs(fNSigmaDCAzHad) < 5)',
]
#selections = conf['selections']
base_selection = selections[0]
for sel in selections[1:]:
    base_selection += (' && ' + sel)
print(f'Selection: {base_selection}')

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
    elif 'fNSigmaTPCHad' not in rdf.GetColumnNames():
        raise RuntimeError("Neither fNSigmaTPCHadPr nor fNSigmaTPCHad found in dataset!")
    
    # TOF
    if 'fNSigmaTOFHadPr' in rdf.GetColumnNames() and 'fNSigmaTOFHad' in rdf.GetColumnNames():
        rdf = rdf.Define('fNSigmaTOFHad', 'fNSigmaTOFHadPr')
    elif 'fNSigmaTOFHad' not in rdf.GetColumnNames():
        #if 'fNSigmaTOFHadPr' not in rdf.GetColumnNames():
        #    rdf = rdf.Define('fNSigmaTOFHad', 'ComputeNsigmaTOFPr(std::abs(fPtHad), fMassTOFHad)')
        #else:
        #    rdf = rdf.Define('fNSigmaTOFHad', 'fNSigmaTOFHadPr')
       rdf = rdf.Define('fNSigmaTOFHad', 'ComputeNsigmaTOFPr(std::abs(fPtHad), fMassTOFHad)')
    
      # Recalibration
      # Correct for PID in tracking
      
      # close pair rejection
      #.Define('fDeltaEta', 'fEtaHe3 - fEtaHad') \
      #.Define('fDeltaPhi', 'fPhiHe3 - fPhiHad') \
      
    rdf = rdf.Define('fSignedPtHad', 'fPtHad') \
      .Define('fSignHe3', 'fPtHe3/std::abs(fPtHe3)') \
      .Redefine('fPtHe3', 'std::abs(fPtHe3)') \
      .Redefine('fPtHad', 'std::abs(fPtHad)') \
      .Redefine('fPtHe3', '(fPIDtrkHe3 == 7) || (fPIDtrkHe3 == 8) || (fPtHe3 > 2.5) ? fPtHe3 : CorrectPidTrkHe(fPtHe3)') \
      .Redefine('fInnerParamTPCHe3', 'fInnerParamTPCHe3 * 2') \
      .Redefine('fInnerParamTPCHe3', '(fPIDtrkHe3 == 7) || (fPIDtrkHe3 == 8) || (fPtHe3 > 2.5) ? fInnerParamTPCHe3 : CorrectPidTrkHe(fInnerParamTPCHe3, false)') \
      .Define('fSignedPtHe3', 'fPtHe3 * fSignHe3') \
      .Define(f'fEHe3', f'std::sqrt((fPtHe3 * std::cosh(fEtaHe3))*(fPtHe3 * std::cosh(fEtaHe3)) + {ParticleMasses["He"]}*{ParticleMasses["He"]})') \
      .Define(f'fEHad', f'std::sqrt((fPtHad * std::cosh(fEtaHad))*(fPtHad * std::cosh(fEtaHad)) + {ParticleMasses["Pr"]}*{ParticleMasses["Pr"]})') \
      .Define('fClusterSizeCosLamHe3', 'ComputeAverageClusterSize(fItsClusterSizeHe3) / cosh(fEtaHe3)') \
      .Define('fClusterSizeCosLamHad', 'ComputeAverageClusterSize(fItsClusterSizeHad) / cosh(fEtaHad)') \
      .Define('fExpectedClusterSizeHe3', 'ComputeExpectedClusterSizeCosLambdaHe(fPtHe3 * std::cosh(fEtaHe3))') \
      .Define('fExpectedClusterSizeHad', 'ComputeExpectedClusterSizeCosLambdaPr(fPtHad * std::cosh(fEtaHad))') \
      .Define('fNSigmaITSHe3', 'ComputeNsigmaITSHe(fPtHe3 * std::cosh(fEtaHe3), fClusterSizeCosLamHe3)') \
      .Define('fNSigmaITSHad', 'ComputeNsigmaITSPr(fPtHad * std::cosh(fEtaHad), fClusterSizeCosLamHad)') \
      .Redefine('fNSigmaTPCHe3', 'ComputeNsigmaTPCHe(std::abs(fInnerParamTPCHe3), fSignalTPCHe3, false, true)') \
      .Define('fNSigmaTPCPi', 'ComputeNsigmaTPCPi(std::abs(fInnerParamTPCHad), fSignalTPCHad)') \
      .Define('fNSigmaDCAxyHe3', 'ComputeNsigmaDCAxyHe(fPtHe3, fDCAxyHe3)') \
      .Define('fNSigmaDCAzHe3', 'ComputeNsigmaDCAzHe(fPtHe3, fDCAzHe3)') \
      .Define('fNSigmaDCAxyHad', 'ComputeNsigmaDCAxyPr(fPtHad, fDCAxyHad)') \
      .Define('fNSigmaDCAzHad', 'ComputeNsigmaDCAzPr(fPtHad, fDCAzHad)') \
      .Filter(base_selection).Filter(selection) \
      .Define('fKstar', f'ComputeKstar(fPtHe3, fEtaHe3, fPhiHe3, {ParticleMasses["He"]}, fPtHad, fEtaHad, fPhiHad, {ParticleMasses["Pr"]})')
    
    return rdf

def load_same():
    
    config_same = {
        '2023': {'input_data':'/data/galucia/lithium/same/LHC23_PbPb_pass5_hadronpid_same.root',},
        '2024': {'input_data':'/data/galucia/lithium/same/LHC24ar_pass3_hadronpid_same.root',},
        '2025': {'input_data':'/data/galucia/lithium/same/LHC25_PbPb_pass1_hadronpid_same.root',},
    }
    
    for year in config_same.keys():
        config_same[year].update({
            'tree_names':   ['O2he3hadtable', 'O2he3hadmult'],
            'mode':         'DF',
            'parameterisation': PARAMETRISATIONS[year],
        })
        
        load_parametrisation(config_same[year])
        chain_data_same, __ = prepare_input_tchain(config_same[year])
        rdf_same = prepare_rdataframe(chain_data_same, base_selection, 'true')
        rdf_same.Snapshot(f'rdf_same', f'/data/galucia/lithium/tmp/rdf_same_{year}.root')
        
    return RDataFrame('rdf_same', [f'/data/galucia/lithium/tmp/rdf_same_{year}.root' for year in config_same.keys()])

def load_mixed():
    
    config_mixed = {
            '2023': {'input_data': '/data/galucia/lithium/mixing/LHC23_PbPb_pass5_hadronpid_event_mixing_batch256.root'},
            '2024': {'input_data': '/data/galucia/lithium/mixing/LHC24ar_pass3_hadronpid_event_mixing_batch256.root'},
            '2025': {'input_data': '/data/galucia/lithium/mixing/LHC25_PbPb_pass1_hadronpid_event_mixing_batch256.root'},
    }
    
    for year in config_mixed.keys():
        config_mixed[year].update({
            'tree_names':    'MixedTree',
            'mode':         'tree',
            'parameterisation': PARAMETRISATIONS[year],
        })
        
        load_parametrisation(config_mixed[year])
        chain_data_mixed, __ = prepare_input_tchain(config_mixed[year])
        rdf_mixed = prepare_rdataframe(chain_data_mixed, base_selection, 'true')
        rdf_mixed.Snapshot(f'rdf_mixed', f'/data/galucia/lithium/tmp/rdf_mixed_{year}.root')
    
    return RDataFrame('rdf_mixed', [f'/data/galucia/lithium/tmp/rdf_mixed_{year}.root' for year in config_mixed.keys()])

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
                         .Define('fNsigmaTpcMinHe3Systematic', 'std::abs(fNSigmaTPCHe3)') \
                         .Define('fNsigmaTpcMaxHe3Systematic', 'std::abs(fNSigmaTPCHe3)') \
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
                         .Vary(['fNSigmaITSHe3Systematic', 'fNsigmaTpcMinHe3Systematic', 'fNsigmaTpcMaxHe3Systematic',
                                'fAbsNSigmaDCAxyHe3Systematic', 'fAbsNSigmaDCAzHe3Systematic',
                                'fNClsTPCHe3Systematic', 'fChi2TPCHe3Systematic',
                                'fNSigmaITSHadSystematic', 'fAbsNSigmaTPCHadSystematic', 'fAbsNSigmaTOFHadSystematic',
                                'fAbsNSigmaDCAxyHadSystematic', 'fAbsNSigmaDCAzHadSystematic', 'fAbsZVertexSystematic'],
                                f'systematicCuts({n_variations}, fNSigmaITSHe3Systematic, fNsigmaTpcMinHe3Systematic, fNsigmaTpcMaxHe3Systematic, \
                                fAbsNSigmaDCAxyHe3Systematic, fAbsNSigmaDCAzHe3Systematic, fNClsTPCHe3Systematic, fChi2TPCHe3Systematic, \
                                fNSigmaITSHadSystematic, fAbsNSigmaTPCHadSystematic, fAbsNSigmaTOFHadSystematic, \
                                fAbsNSigmaDCAxyHadSystematic, fAbsNSigmaDCAzHadSystematic, fChi2TPCHadSystematic, \
                                fAbsZVertexSystematic)',
                                n_variations, 'systematic') \
                         .Filter('(fNSigmaITSHe3Systematic > 0) && \
                                  (fNsigmaTpcMinHe3Systematic > 0) && \
                                  (fNsigmaTpcMaxHe3Systematic < 0) && \
                                  (fNSigmaITSHadSystematic > 0) && \
                                  (fAbsNSigmaTPCHadSystematic < 0) && \
                                  ((std::abs(fPtHad) < 0.8) || (fAbsNSigmaTOFHadSystematic < 0)) && \
                                  (fAbsZVertexSystematic < 0) \
                                  (fAbsNSigmaDCAxyHe3Systematic < 0) && \
                                  (fAbsNSigmaDCAzHe3Systematic < 0) && \
                                  (fAbsNSigmaDCAxyHadSystematic < 0) && \
                                  (fAbsNSigmaDCAzHadSystematic < 0) && \
                                  (fChi2TPCHe3Systematic < 0) && \
                                  (fChi2TPCHadSystematic < 0) && \
                                  (fNClsTPCHe3Systematic < 0) ) ')
    else:
        variational_rdf = rdf.Define('fNSigmaITSHe3Systematic', 'fNSigmaITSHe3') \
                            .Define('fNsigmaTpcMinHe3Systematic', 'fNSigmaTPCHe3') \
                            .Define('fNsigmaTpcMaxHe3Systematic', 'fNSigmaTPCHe3') \
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
                            .Vary(['fNSigmaITSHe3Systematic', 'fNsigmaTpcMinHe3Systematic', 'fNsigmaTpcMaxHe3Systematic',
                                    'fAbsNSigmaDCAxyHe3Systematic', 'fAbsNSigmaDCAzHe3Systematic',
                                    'fNClsTPCHe3Systematic', 'fChi2TPCHe3Systematic',
                                    'fNSigmaITSHadSystematic', 'fAbsNSigmaTPCHadSystematic', 'fAbsNSigmaTOFHadSystematic',
                                    'fAbsNSigmaDCAxyHadSystematic', 'fAbsNSigmaDCAzHadSystematic', 'fAbsZVertexSystematic'],
                                    f'systematicCutsSingleVariable({n_variations}, "{variable}", fNSigmaITSHe3Systematic, fNsigmaTpcMinHe3Systematic, fNsigmaTpcMaxHe3Systematic, \
                                    fAbsNSigmaDCAxyHe3Systematic, fAbsNSigmaDCAzHe3Systematic, fNClsTPCHe3Systematic, fChi2TPCHe3Systematic, \
                                    fNSigmaITSHadSystematic, fAbsNSigmaTPCHadSystematic, fAbsNSigmaTOFHadSystematic, \
                                    fAbsNSigmaDCAxyHadSystematic, fAbsNSigmaDCAzHadSystematic, fChi2TPCHadSystematic, \
                                    fAbsZVertexSystematic)',
                                    n_variations, 'systematic') \
                            .Filter('(fNSigmaITSHe3Systematic > 0) && \
                                  (fNsigmaTpcMinHe3Systematic < 0) && \
                                  (fNsigmaTpcMaxHe3Systematic < 0) && \
                                  (fNSigmaITSHadSystematic > 0) && \
                                  (fAbsNSigmaTPCHadSystematic < 0) && \
                                  ((std::abs(fPtHad) < 0.8) || (fAbsNSigmaTOFHadSystematic < 0)) ')
                            
    nominal_hists = {}
    variational_hists = {}
                            
    for centrality, centrality_condition in zip(['010', '1030', '3050', '5080'], 
                                                ['fCentralityFT0C < 10', 
                                                 'fCentralityFT0C >= 10 && fCentralityFT0C < 30', 
                                                 'fCentralityFT0C >= 30 && fCentralityFT0C < 50',
                                                 'fCentralityFT0C >= 50 && fCentralityFT0C < 80']):
        nominal_hists[centrality] = variational_rdf.Filter(centrality_condition) \
                                        .Histo1D((f'hKstar{centrality}{hist_name_suffix}', ';#it{k}* (GeV/#it{c});', 20, 0, 0.4), 'fKstar') 
        variational_hists[centrality] = ROOT.RDF.Experimental.VariationsFor(nominal_hists[centrality])
        
    return variational_hists

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

    SIGNAL_HIST_LOAD_INFO = HistLoadInfo('/home/galucia/Lithium4/femto/models/li4_contribution_proper_sill.root', 'hCkHist')
    
    workspace = RooWorkspace('roows')

    h_bkg = load_hist(bkg_input_path, h_bkg_name)
    h_signal = load_hist(SIGNAL_HIST_LOAD_INFO)
    h_correlation_function = h_data

    KSTAR_MIN, KSTAR_MAX = (0.01, 0.4) if h_correlation_function.GetBinWidth(1) > 1e-5 else (0.02, 0.4)
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

    rdf_same = load_same()
    rdf_mixed = load_mixed()

    N_ITERATIONS = 100
    outFile = TFile("output/hist_systematics_with_upper_limit.root", "RECREATE")

    for sign, condition in {'Matter': 'fSignedPtHe3 > 0', 'Antimatter': 'fSignedPtHe3 < 0', 'Both': 'true'}.items():

        outDir = outFile.mkdir(sign)

        tmp_rdf_same = rdf_same.Filter(condition)
        tmp_rdf_mixed = rdf_mixed.Filter(condition)

        variational_hists_same = prepare_systematics_histograms(tmp_rdf_same, n_variations=1, hist_name_suffix='Same', variable='')
        variational_hists_mixed = prepare_systematics_histograms(tmp_rdf_mixed, n_variations=1, hist_name_suffix='Mixed', variable='')

        NORM_LOW_KSTAR, NORM_HIGH_KSTAR = 0.2, 0.4
        
        for centrality in ['010']: #, '1030', '3050', '5080']:
            
            outDirCentrality = outDir.mkdir(centrality)
             
            check_procedure = False
            if check_procedure:
                h_same = variational_hists_same[centrality][f'systematic:{0}']
                h_mixed = variational_hists_mixed[centrality][f'systematic:{0}']
                h_mixeds_normalised = normalise_histogram(h_same, h_mixed, NORM_LOW_KSTAR, NORM_HIGH_KSTAR)
                
                outDirCentrality.cd()
                gc.collect()

            variational_hists_same = prepare_systematics_histograms(tmp_rdf_same, n_variations=N_ITERATIONS, hist_name_suffix='Same', variable='all')
            variational_hists_mixed = prepare_systematics_histograms(tmp_rdf_mixed, n_variations=N_ITERATIONS, hist_name_suffix='Mixed', variable='all')

            h_sames, h_mixeds, h_mixeds_normalised = [], [], []
            NORM_LOW_KSTAR, NORM_HIGH_KSTAR = 0.2, 0.4

            for iter in range(N_ITERATIONS):

                print(f'Processing iteration {iter+1}/{N_ITERATIONS} for {sign}...')

                same = variational_hists_same[centrality][f'systematic:{iter}']
                mixed = variational_hists_mixed[centrality][f'systematic:{iter}']

                mixed_normalised = normalise_histogram(same, mixed, NORM_LOW_KSTAR, NORM_HIGH_KSTAR)

                h_sames.append(same)                        
                h_mixeds.append(mixed)                      
                h_mixeds_normalised.append(mixed_normalised)

            for iter, (sames, mixeds, mixeds_normalised) in enumerate(zip(h_sames, h_mixeds, h_mixeds_normalised)):
                outDirIter = outDirCentrality.mkdir(f'iter_{iter}')
                outDirIter.cd()
                sames.Write(f'hSame')
                mixeds.Write(f'hMixed')
                mixeds_normalised.Write(f'hMixedNormalised')
    
    outFile.Close()

def upper_limit_systematic_routine():

    N_ITERATIONS = 20
    N_UPPER_LIMIT_ITERATIONS = 500
    TH1.AddDirectory(False)

    infile = TFile.Open("output/hist_systematics_with_upper_limit_010.root")
    outFile = TFile.Open("output/systematics_with_upper_limit_010_v3.root", "RECREATE")

    # (nominal, upper, lower): radius variations
    # (-blank-, higher, lower): lambda variations
    AVAILABLE_BKGS = [
        f'nominal/hLambdaSigmaCorrectedCk_Smeared_nominal',
        f'nominal/hLambdaSigmaCorrectedCk_Smeared_nominal_higher',
        f'nominal/hLambdaSigmaCorrectedCk_Smeared_nominal_lower',
        f'upper/hLambdaSigmaCorrectedCk_Smeared_lower',
        f'upper/hLambdaSigmaCorrectedCk_Smeared_lower_higher',
        f'upper/hLambdaSigmaCorrectedCk_Smeared_lower_lower',
        f'lower/hLambdaSigmaCorrectedCk_Smeared_upper',
        f'lower/hLambdaSigmaCorrectedCk_Smeared_upper_higher',
        f'lower/hLambdaSigmaCorrectedCk_Smeared_upper_lower',
        ]
    BKG_PATH = '/home/galucia/Lithium4/femto/models/LHC25_PbPb_pass1_lambda_models.root'

    for sign in ['Matter']:

        inDir = infile.Get(sign)
        outDir = outFile.mkdir(sign)

        # --- new: centrality dict, same pattern as prepare_histograms ---
        h_upper_limits = {}
        for centrality in ['010']:  # extend as needed: '1030', '3050', '5080'
            

            h_upper_limits[centrality] = TH1F(f'hUpperLimits_{centrality}', ';[#it{N}_{^{4}Li}^{raw}]_{upper limit};', 600, 0, 600)
            inDirCentrality = inDir.Get(centrality)  # assumes subdirs per centrality in infile
            outDirCentrality = outDir.mkdir(centrality)

            for iter in tqdm(range(N_ITERATIONS)):

                outDirIter = outDirCentrality.mkdir(f'{centrality}/iter_{iter}')

                h_same             = inDirCentrality.Get(f'iter_{iter}/hSame')
                h_mixed            = inDirCentrality.Get(f'iter_{iter}/hMixed')
                h_mixed_normalised = inDirCentrality.Get(f'iter_{iter}/hMixedNormalised')

                h_correlation = h_same.Clone(f'hCorrelation_{centrality}_{iter}')
                h_correlation.Divide(h_mixed_normalised)

                outDirNominal = outDirIter.mkdir('nominal')
                random_bkg_name = np.random.choice(AVAILABLE_BKGS)

                nominal_raw_yield = fitting_routine(outDirNominal,
                                bkg_input_path=BKG_PATH,
                                h_bkg_name=f'{sign}/{centrality}/{random_bkg_name}',
                                h_data=h_correlation, h_mixed_event=h_mixed_normalised,
                                sign=sign, centrality=centrality, use_smoothening=True)

                h_raw_yields_iter = TH1F(f'hRawYields_{centrality}_Iter_{iter}', ';Raw yield;', 300, -200, 400)

                for iter_upper in range(N_UPPER_LIMIT_ITERATIONS):

                    h_same_iter_upper             = h_same.Clone(f'hSame_{centrality}_Iter_Upper_{iter_upper}')
                    h_mixed_normalised_iter_upper = h_mixed_normalised.Clone(f'hMixed_{centrality}_Iter_Upper_{iter_upper}')

                    h_same_iter_upper             = poisson_sampling(h_same_iter_upper)
                    h_mixed_normalised_iter_upper = poisson_sampling(h_mixed_normalised_iter_upper)

                    h_correlation_iter_upper = h_same_iter_upper.Clone(f'hCorrelation_{centrality}_Iter_Upper_{iter_upper}')
                    h_correlation_iter_upper.Divide(h_mixed_normalised_iter_upper)

                    outDir_to_pass = outDirIter.mkdir(f'inner_iter_{iter_upper}') if iter_upper == 0 else None
                    raw_yield = fitting_routine(outDir_to_pass,
                                    bkg_input_path=BKG_PATH,
                                    h_bkg_name=f'{sign}/{centrality}/hLambdaSigmaCorrectedCk',
                                    h_data=h_correlation_iter_upper, h_mixed_event=h_mixed_normalised_iter_upper,
                                    sign=sign, centrality=centrality, use_smoothening=True)

                    h_raw_yields_iter.Fill(raw_yield)
                    del h_same_iter_upper, h_mixed_normalised_iter_upper, h_correlation_iter_upper

                fit_func = TF1(f'fit_func_{centrality}_{iter}', 'gaus', -200, 400)
                fit_func.SetParameters(h_raw_yields_iter.GetMaximum(), h_raw_yields_iter.GetMean(), h_raw_yields_iter.GetRMS())
                h_raw_yields_iter.Fit('gaus', 'RMS+', '', -200, 400)

                outDirIter.cd()
                h_raw_yields_iter.Write()

                upper_limit = nominal_raw_yield + 1.96 * fit_func.GetParameter(2)
                h_upper_limits[centrality].Fill(upper_limit)

                gc.collect()
                del h_same, h_mixed, h_mixed_normalised, h_correlation, h_raw_yields_iter, fit_func

            outDir.cd()
            h_upper_limits[centrality].Write()
            del h_upper_limits[centrality]

    infile.Close()
    outFile.Close()




if __name__ == "__main__":

    RooMsgService.instance().setGlobalKillBelow(5) # 3 = WARNING, 4 = ERROR, 5 = FATAL
    ROOT.RooFit.PrintLevel(-1)

    prepare_histograms()

    #upper_limit_systematic_routine()

    print("Systematics analysis completed successfully.")