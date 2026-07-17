import sys
import yaml
import argparse

import ROOT
from ROOT import TFile, TChain, gInterpreter, RDataFrame

from torchic.utils.terminal_colors import TerminalColors as tc

sys.path.append('..')
from utils.particles import ParticleMasses

from include.load_parameters import load_parametrisation

gInterpreter.ProcessLine(f'#include "../include/Common.h"')

ROOT.EnableImplicitMT(30)
ROOT.gROOT.SetBatch(True)


def prepare_selections(selections):

    base_selection = ''
    if config.get('like_sign', True):  
        base_selection = '(fSignedPtHe3 > 0 && fSignedPtHad > 0) || (fSignedPtHe3 < 0 && fSignedPtHad < 0)'
    else:
        base_selection = '(fSignedPtHe3 > 0 && fSignedPtHad < 0) || (fSignedPtHe3 < 0 && fSignedPtHad > 0)'
    
    selection = selections[0]
    for sel in selections[1:]:
        selection += (' && ' + sel)

    return base_selection, selection

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
      .Define('fDeltaEta', 'fEtaHe3 - fEtaHad') \
      .Define('fDeltaPhi', 'fPhiHe3 - fPhiHad') \
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
      .Define('fMtHad', f'std::sqrt((fPtHad)*(fPtHad) + {ParticleMasses["Pr"]}*{ParticleMasses["Pr"]})') \
      .Define('fMtHe3', f'std::sqrt((fPtHe3)*(fPtHe3) + {ParticleMasses["He"]}*{ParticleMasses["He"]})')
    
    return rdf



if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Prepare data for femtoscopic analysis')
    parser.add_argument('--config', type=str, default='config/config_prepare.yml', help='Path to the configuration YAML file')  
    parser.add_argument('--mode', type=str, default='same_event', choices=['same_event', 'mixed_event'], help='Mode of analysis: same_event or mixed_event')
    args = parser.parse_args()

    config_file = args.config
    config = yaml.safe_load(open(config_file, 'r'))
    
    load_parametrisation(config)  # Load parametrisation into Common.h
    
    selections = config.get('selections', [])
    config = config[args.mode]  # Use the specific mode section from the config

    base_selection, selection = prepare_selections(selections)
    chain_data, additional_chain_data = prepare_input_tchain(config)
    rdf = prepare_rdataframe(chain_data, base_selection, selection)

    output_file = ROOT.TFile('output/mt.root', "RECREATE")
    
    for centrality, centrality_selection in {'010': 'fCentralityFT0C >= 0 && fCentralityFT0C < 10', 
                                             '1030': 'fCentralityFT0C >= 10 && fCentralityFT0C < 30', 
                                             '3050': 'fCentralityFT0C >= 30 && fCentralityFT0C < 50',
                                             '5080': 'fCentralityFT0C >= 50 && fCentralityFT0C < 80',
                                             '1050': 'fCentralityFT0C >= 10 && fCentralityFT0C < 50'}.items():
        output_dir = output_file.mkdir(centrality)
        hist_centrality_had = rdf.Filter(centrality_selection).Histo1D(('hMtHad', ';m_{T} (GeV/c^{2}); Counts', 200, 0, 10), 'fMtHad').GetValue()
        hist_centrality_he3 = rdf.Filter(centrality_selection).Histo1D(('hMtHe3', ';m_{T} (GeV/c^{2}); Counts', 200, 0, 10), 'fMtHe3').GetValue()
        
        output_dir.cd()
        hist_centrality_had.Write()
        hist_centrality_he3.Write()
        
    output_file.Close()

  