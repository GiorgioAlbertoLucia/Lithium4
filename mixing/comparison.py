'''
    Run a comparison between event mixing produced locally and the one produced on the grid
'''

import sys
import yaml
import ROOT
import itertools
from ROOT import TFile, TChain, gInterpreter, RDataFrame, TCanvas, TText

from torchic.utils.terminal_colors import TerminalColors as tc
from torchic.utils.root import set_root_object

sys.path.append('..')
from utils.particles import ParticleMasses
from utils.histogram_registry import HistogramRegistry
from utils.histogram_archive import register_qa_histograms, select_qa_histograms

gInterpreter.ProcessLine(f'#include "../include/Common.h"')
from ROOT import ComputeNsigmaTPCHe, ComputeNsigmaITSHe, ComputeNsigmaITSPr, ComputeNsigmaTOFPr, ComputeAverageClusterSize, CorrectPidTrkHe, \
  ComputeKstar, ComputeExpectedClusterSizeCosLambdaHe, ComputeExpectedClusterSizeCosLambdaPr, ComputeNsigmaDCAxyHe, ComputeNsigmaDCAxyPr, ComputeNsigmaDCAzHe, ComputeNsigmaDCAzPr 

ROOT.EnableImplicitMT(30)
ROOT.gROOT.SetBatch(True)

HISTOGRAMS = ['hKstar', 'hPtHe', 'hEtaHe', 'hPhiHe']


def prepare_selections(config):

    base_selection = ''
    if config.get('like_sign', True):  
        base_selection = '(fSignedPtHe3 > 0 && fSignedPtHad > 0) || (fSignedPtHe3 < 0 && fSignedPtHad < 0)'
    else:
        base_selection = '(fSignedPtHe3 > 0 && fSignedPtHad < 0) || (fSignedPtHe3 < 0 && fSignedPtHad > 0)'

    selections = config['selections']
    
    selection = selections[0]
    for sel in selections[1:]:
        selection += (' && ' + sel)

    return base_selection, selection

def prepare_input_tchain(input_data, tree_name, mode):

    file_data_list = input_data if isinstance(input_data, list) else [input_data]
    chain_data = TChain('tchain')

    for file_name in file_data_list:
        fileData = TFile(file_name)

        if mode == 'DF':
            for key in fileData.GetListOfKeys():
                key_name = key.GetName()
                if 'DF_' in key_name :
                    print(f'Adding {tc.CYAN+tc.UNDERLINE}{file_name}/{key_name}/{tree_name}{tc.RESET} to the chain')
                    chain_data.Add(f'{file_name}/{key_name}/{tree_name}')
        elif mode == 'tree':
            print(f'Adding {tc.CYAN+tc.UNDERLINE}{file_name}/{tree_name}{tc.RESET} to the chain')
            chain_data.Add(f'{file_name}/{tree_name}')
        else:
            raise ValueError('Invalid mode name')
    
    return chain_data

def prepare_rdataframe(chain_data: TChain, base_selection: str, selection: str):
   
    rdf = RDataFrame(chain_data) 
    if 'fNSigmaTPCHadPr' in rdf.GetColumnNames():
        rdf = rdf.Define('fNSigmaTPCHad', 'fNSigmaTPCHadPr')
    if 'fNSigmaTOFHadPr' not in rdf.GetColumnNames():
        rdf = rdf.Define('fNSigmaTOFHad', 'ComputeNsigmaTOFPr(std;;abs(fPtHad), fMassTOFHad)')
    else:
       rdf = rdf.Define('fNSigmaTOFHad', 'fNSigmaTOFHadPr')
    
      # Recalibration
      #.Redefine('fNSigmaTPCHe3', 'ComputeNsigmaTPCHe(std::abs(fInnerParamTPCHe3), fSignalTPCHe3)') \ 
      # Correct for PID in tracking
      
    rdf = rdf.Define('fSignedPtHad', 'fPtHad') \
      .Define('fSignHe3', 'fPtHe3/std::abs(fPtHe3)') \
      .Redefine('fPtHe3', 'std::abs(fPtHe3)') \
      .Redefine('fPtHad', 'std::abs(fPtHad)') \
      .Define('fSignedPtHe3', 'fPtHe3 * fSignHe3') \
      .Redefine('fPtHe3', '(fPIDtrkHe3 == 7) || (fPIDtrkHe3 == 8) || (fPtHe3 > 2.5) ? fPtHe3 : CorrectPidTrkHe(fPtHe3)') \
      .Define(f'fEHe3', f'std::sqrt((fPtHe3 * std::cosh(fEtaHe3))*(fPtHe3 * std::cosh(fEtaHe3)) + {ParticleMasses["He"]}*{ParticleMasses["He"]})') \
      .Define(f'fEHad', f'std::sqrt((fPtHad * std::cosh(fEtaHad))*(fPtHad * std::cosh(fEtaHad)) + {ParticleMasses["Pr"]}*{ParticleMasses["Pr"]})') \
      .Define('fDeltaEta', 'fEtaHe3 - fEtaHad') \
      .Define('fDeltaPhi', 'fPhiHe3 - fPhiHad') \
      .Redefine('fInnerParamTPCHe3', 'fInnerParamTPCHe3 * 2') \
      .Define('fClusterSizeCosLamHe3', 'ComputeAverageClusterSize(fItsClusterSizeHe3) / cosh(fEtaHe3)') \
      .Define('fClusterSizeCosLamHad', 'ComputeAverageClusterSize(fItsClusterSizeHad) / cosh(fEtaHad)') \
      .Define('fExpectedClusterSizeHe3', 'ComputeExpectedClusterSizeCosLambdaHe(fPtHe3 * std::cosh(fEtaHe3))') \
      .Define('fExpectedClusterSizeHad', 'ComputeExpectedClusterSizeCosLambdaPr(fPtHad * std::cosh(fEtaHad))') \
      .Define('fNSigmaITSHe3', 'ComputeNsigmaITSHe(fPtHe3 * std::cosh(fEtaHe3), fClusterSizeCosLamHe3)') \
      .Define('fNSigmaITSHad', 'ComputeNsigmaITSPr(fPtHad * std::cosh(fEtaHad), fClusterSizeCosLamHad)') \
      .Define('fNSigmaDCAxyHe3', 'ComputeNsigmaDCAxyHe(fPtHe3, fDCAxyHe3)') \
      .Define('fNSigmaDCAzHe3', 'ComputeNsigmaDCAzHe(fPtHe3, fDCAzHe3)') \
      .Define('fNSigmaDCAxyHad', 'ComputeNsigmaDCAxyPr(fPtHad, fDCAxyHad)') \
      .Define('fNSigmaDCAzHad', 'ComputeNsigmaDCAzPr(fPtHad, fDCAzHad)') \
      .Filter(base_selection).Filter(selection) \
      .Define('fKstar', f'ComputeKstar(fPtHe3, fEtaHe3, fPhiHe3, {ParticleMasses["He"]}, fPtHad, fEtaHad, fPhiHad, {ParticleMasses["Pr"]})')
    
    return rdf

def prepare_registry(rdf, output_file: TFile):
   
    print(f'\n{tc.GREEN}Creating histograms{tc.RESET}')
    histogram_registry = HistogramRegistry()

    select_qa_histograms(HISTOGRAMS)
    register_qa_histograms(histogram_registry)
    histogram_registry.prepare_directories(output_file)
    histogram_registry.draw_histogram(rdf)

    histogram_registry.save_histograms(output_file)

    return histogram_registry



if __name__ == '__main__':

    config_file = '/home/galucia/Lithium4/preparation/config/config_prepare.yml'
    config = yaml.safe_load(open(config_file, 'r'))

    base_selection, selection = prepare_selections(config)

    input_datas = [
        ['/data/galucia/lithium_local/mixing/LHC23_PbPb_pass4_all_event_mixing_batch1995.root',
         '/data/galucia/lithium_local/mixing/LHC24ar_pass1_all_event_mixing_batch1995.root',
         '/data/galucia/lithium_local/mixing/LHC24as_pass1_all_event_mixing_batch1995.root',
        ],
        ['/data/galucia/lithium_local/mixing/LHC23_PbPb_pass4_all_event_mixing_batch105.root',
         '/data/galucia/lithium_local/mixing/LHC24ar_pass1_all_event_mixing_batch105.root',
         '/data/galucia/lithium_local/mixing/LHC24as_pass1_all_event_mixing_batch105.root',
         ],
        ['/data/galucia/lithium_local/mixing/LHC23_PbPb_pass4_all_event_mixing_batch42.root',
         '/data/galucia/lithium_local/mixing/LHC24ar_pass1_all_event_mixing_batch42.root',
         '/data/galucia/lithium_local/mixing/LHC24as_pass1_all_event_mixing_batch42.root',
         ],
         ['/data/galucia/lithium_local/mixing_merged/LHC23_PbPb_pass4_all_event_mixing_grid.root',
          '/data/galucia/lithium_local/mixing_merged/LHC24ar_pass1_all_event_mixing_grid.root',
          '/data/galucia/lithium_local/mixing_merged/LHC24as_pass1_all_event_mixing_grid.root'
          ]
        ]
    tree_names = ['MixedTree', 'MixedTree', 'MixedTree', 'O2he3hadtable']
    modes = ['tree', 'tree', 'tree', 'DF']
    folder_names = ['1995', '105', '42', 'grid']
    
    output_file_path = 'output/comparison.root'
    output_file = ROOT.TFile(output_file_path, "RECREATE")

    histogram_registries = {}

    for input_data, tree_name, mode, folder_name in zip(input_datas, tree_names, modes, folder_names):
       
        outdir = output_file.mkdir(folder_name)

        chain_data = prepare_input_tchain(input_data, tree_name, mode)
        rdf = prepare_rdataframe(chain_data, base_selection, selection)

        histogram_registries[folder_name] = prepare_registry(rdf, outdir)

    for pair in itertools.combinations(folder_names, 2):
        folder1, folder2 = pair
        outdir = output_file.mkdir(f'{folder1}_{folder2}')
        
        for histogram in HISTOGRAMS:

            canvas = TCanvas(f'c_{histogram}')

            set_root_object(histogram_registries[folder1][histogram], line_color=797, line_width=2)
            histogram_registries[folder1][histogram].Scale(1./histogram_registries[folder1][histogram].Integral())
            histogram_registries[folder1][histogram].Draw('hist')
            
            set_root_object(histogram_registries[folder2][histogram], line_color=402, line_width=2)
            histogram_registries[folder2][histogram].Scale(1./histogram_registries[folder2][histogram].Integral())
            histogram_registries[folder2][histogram].Draw('hist same')

            hist2 = histogram_registries[folder2][histogram].GetValue()

            kolmogorov_probability = \
             histogram_registries[folder1][histogram].KolmogorovTest(hist2)
            text = TText(0.25, 0.2, f'KS probability = {kolmogorov_probability:.2f}')
            text.SetNDC()
            text.Draw('same')

            ratio = histogram_registries[folder1][histogram].Clone(f'{histogram}_ratio')
            ratio.Divide(hist2)

            outdir.cd()
            canvas.Write()
            ratio.Write()
    

    output_file.Close()

  
