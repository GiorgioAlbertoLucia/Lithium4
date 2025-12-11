'''
    Study of the dca distribution for p-He3 pairs from data
'''

import sys
from enum import Enum

from ROOT import TFile, TChain, RDataFrame, TDirectory, \
    gInterpreter, gROOT, EnableImplicitMT

from torchic.utils.terminal_colors import TerminalColors as tc

gInterpreter.ProcessLine(f'#include "../include/Common.h"')
from ROOT import ComputeNsigmaTPCHe, ComputeNsigmaITSHe, ComputeNsigmaITSPr, ComputeNsigmaTOFPr, ComputeAverageClusterSize, CorrectPidTrkHe

EnableImplicitMT(5)
gROOT.SetBatch(True)

sys.path.append('..')
from utils.particles import ParticlePDG, ParticleMasses

class Flags(Enum):
    kProton = 0
    kDeuteron = 1
    kTriton = 2
    kHe3 = 3
    kHe4 = 4
    kHasTOF = 5
    kHasTRD = 6
    kIsAmbiguous = 7
    kITSrof = 8
    kIsPhysicalPrimary = 9
    kIsSecondaryFromMaterial = 10
    kIsSecondaryFromWeakDecay = 11

def visualise(rdf: RDataFrame, outfile:TDirectory, particle:str = 'He'):

    particle_suffix = 'He3' if particle == 'He' else 'Had'
    h2_dcaxy_pt = rdf.Histo2D((f'h2DCAxyPt{particle}', ';#it{p}_{T} (GeV/#it{c});DCA_{#it{xy}} (cm)', 
                               100, -5, 5, 150, -0.15, 0.15), f'fSignedPt{particle_suffix}', f'fDCAxy{particle_suffix}')
    h2_dcaz_pt = rdf.Histo2D((f'h2DCAzPt{particle}', ';#it{p}_{T} (GeV/#it{c});DCA_{#it{z}} (cm)', 
                               100, -5, 5, 150, -0.3, 0.3), f'fSignedPt{particle_suffix}', f'fDCAz{particle_suffix}')
    
    histos = []
    histos.append(h2_dcaxy_pt)
    histos.append(h2_dcaz_pt)

    outfile.cd()
    for hist in histos:
        hist.Write()

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

def prepare_input_tchain(config):

    input_data = config['input_data']
    tree_name = config['tree_name']
    mode = config['mode']

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
    
    return chain_data

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
      #.Redefine('fNSigmaTPCHe3', 'ComputeNsigmaTPCHe(std::abs(fInnerParamTPCHe3), fSignalTPCHe3)') \
      # Correct for PID in tracking
      
    rdf = rdf.Define('fSignedPtHad', 'fPtHad') \
      .Define('fSignHe3', 'fPtHe3/std::abs(fPtHe3)') \
      .Redefine('fPtHe3', 'std::abs(fPtHe3)') \
      .Redefine('fPtHad', 'std::abs(fPtHad)') \
      .Redefine('fPtHe3', '(fPIDtrkHe3 == 7) || (fPIDtrkHe3 == 8) || (fPtHe3 > 2.5) ? fPtHe3 : CorrectPidTrkHe(fPtHe3)') \
      .Define('fSignedPtHe3', 'fPtHe3 * fSignHe3') \
      .Define(f'fEHe3', f'std::sqrt((fPtHe3 * std::cosh(fEtaHe3))*(fPtHe3 * std::cosh(fEtaHe3)) + {ParticleMasses["He"]}*{ParticleMasses["He"]})') \
      .Define(f'fEHad', f'std::sqrt((fPtHad * std::cosh(fEtaHad))*(fPtHad * std::cosh(fEtaHad)) + {ParticleMasses["Pr"]}*{ParticleMasses["Pr"]})') \
      .Redefine('fInnerParamTPCHe3', 'fInnerParamTPCHe3 * 2') \
      .Define('fClusterSizeCosLamHe3', 'ComputeAverageClusterSize(fItsClusterSizeHe3) / cosh(fEtaHe3)') \
      .Define('fClusterSizeCosLamHad', 'ComputeAverageClusterSize(fItsClusterSizeHad) / cosh(fEtaHad)') \
      .Define('fNSigmaITSHe3', 'ComputeNsigmaITSHe(fPtHe3 * std::cosh(fEtaHe3), fClusterSizeCosLamHe3)') \
      .Define('fNSigmaITSHad', 'ComputeNsigmaITSPr(fPtHad * std::cosh(fEtaHad), fClusterSizeCosLamHad)') \
      .Filter(base_selection).Filter(selection)
    
    return rdf

def main():
    
    config = {
                'input_data': [ '/data/galucia/lithium_local/same_merged/LHC23_PbPb_pass4_all_same.root',
                                '/data/galucia/lithium_local/same_merged/LHC24ar_pass1_all_same.root',
                                '/data/galucia/lithium_local/same_merged/LHC24as_pass1_all_same.root',
                                ],
                'tree_name': 'O2he3hadtable',
                'mode': 'DF',
             }

    chain_data = prepare_input_tchain(config)
    selections_dict = {'selections': [
                # quality cuts
                 #'((fDeltaPhi/0.006)*(fDeltaPhi/0.006) + (fDeltaEta/0.012)*(fDeltaEta/0.012) > 1)',
                 '(((fChi2TPCHe3 > 0.5) || (fIs23 == false)) && (fChi2TPCHe3 < 4))', # the cut on chi2tpc should be only applied to 2023
                 '((fChi2TPCHad < 4))',
                 '(std::abs(fEtaHe3) < 0.9)',
                 '(std::abs(fEtaHad) < 0.9)',
                 '(fNClsTPCHe3 > 110)',
                 
                 '((fPtHe3 < 2.5) || (fPIDtrkHe3 == 7))',
                 '((fPIDtrkHe3 == 6) || (fPIDtrkHe3 == 7) || (fPIDtrkHe3 == 8))',
                 #'((fPIDtrkHe3 == 7))',
                 '((-2. < fNSigmaTPCHe3) && (fNSigmaTPCHe3 < 2.))',
                 #'(fClusterSizeCosLamHe3 > 4)',
                 '(fNSigmaITSHe3 > -1.5)',
                 #'(fMassTOFHe3 > 2 || fMassTOFHe3 < 0)',

                 '(std::abs(fNSigmaTPCHad) < 2)',
                 #'(std::abs(fNSigmaTPCHadPr) < 2)',
                 '((std::abs(fPtHad) < 0.8) || (std::abs(fNSigmaTOFHad) < 2))',  # TOF hard cut
                 #'((std::abs(fPtHad) < 0.8) || (std::abs(fNSigmaTOFHad) < 2) || (fNSigmaTOFHad < -9.9))', # TOF veto
                 '(fNSigmaITSHad > -2.)',
                ]}
    base_selections, selections = prepare_selections(selections_dict)
    rdf = prepare_rdataframe(chain_data, base_selections, selections)

    outfile = TFile.Open('output/dca_data_template.root', 'RECREATE')
    
    for particle in ['Pr', 'He']:
        outdir = outfile.mkdir(f'{particle}')
        visualise(rdf, outdir, particle)
        
    outfile.Close()


if __name__ == '__main__':

    main()
