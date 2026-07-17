'''
    Study of the dca distribution for p-He3 pairs from data
'''

import sys
from enum import Enum
import numpy as np

from ROOT import TFile, TChain, RDataFrame, TDirectory, \
    gInterpreter, gROOT, EnableImplicitMT

from torchic.utils.terminal_colors import TerminalColors as tc


EnableImplicitMT(1)
gROOT.SetBatch(True)

sys.path.append('..')
from include.load_parameters import load_parametrisation
gInterpreter.ProcessLine(f'#include "../include/Common.h"')
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

# Custom binning for pT and DCA histograms
NBINS_PT = 200
PT_MIN = -10
PT_MAX = 10
PT_BINNING = np.linspace(PT_MIN, PT_MAX, NBINS_PT + 1)

#NBINS_DCA = 26
#DCA_BINNING = np.array([-0.1, -0.07, -0.05, -0.04, -0.03, -0.025, -0.02, -0.015, -0.0125, -0.01, -0.0075, -0.005, -0.0025, 0,
#                        0.0025, 0.005, 0.0075, 0.01, 0.0125, 0.015, 0.02, 0.025, 0.03, 0.04, 0.05, 0.07, 0.1])

#DCA_BINNING = np.concatenate((np.arange(-0.1, -0.05, 0.005), np.arange(-0.05, 0.05, 0.0025), np.arange(0.05, 0.1, 0.005)))
DCA_BINNING = np.concatenate((np.arange(-0.1, -0.05, 0.005), np.arange(-0.05, -0.025, 0.0025), np.arange(-0.025, 0.025, 0.00125), np.arange(0.025, 0.05, 0.0025), np.arange(0.05, 0.1, 0.005)))
NBINS_DCA = len(DCA_BINNING) - 1

    

def filter_duplicates(rdf: RDataFrame, variable: str, tolerance: float):
    """
    Filter out candidates that have duplicate values within tolerance.
    This works by taking a snapshot, then filtering duplicates in a second pass.
    """
    
    # Define C++ code to filter duplicates
    filter_code = f"""
    std::unordered_set<size_t> seen_duplicates_{variable.replace("f", "")};
    std::vector<float> all_values_{variable.replace("f", "")};
    
    bool FilterDuplicates_{variable.replace("f", "")}(float current_val, ULong64_t entry) {{
        // Check against all previously seen values
        for (size_t i = 0; i < all_values_{variable.replace("f", "")}.size(); ++i) {{
            if (std::abs(all_values_{variable.replace("f", "")}.at(i) - current_val) < {tolerance}) {{
                seen_duplicates_{variable.replace("f", "")}.insert(entry);
                return false;  // Found duplicate, exclude this candidate
            }}
        }}
        all_values_{variable.replace("f", "")}.push_back(current_val);
        return true;  // No duplicate found, keep this candidate
    }}
    """
    
    gInterpreter.Declare(filter_code)
    
    # Apply the filter
    rdf = rdf.Filter(f"FilterDuplicates_{variable.replace('f', '')}({variable}, rdfentry_)")
    
    return rdf

def adjust_bin_content_for_variable_bin_width(h2d):
    minimum_bin_width = 0.
    for ibin in range(1, h2d.GetNbinsY() + 1):
        bin_width = h2d.GetYaxis().GetBinUpEdge(ibin) - h2d.GetYaxis().GetBinLowEdge(ibin)
        if minimum_bin_width == 0. or bin_width < minimum_bin_width:
            minimum_bin_width = bin_width
        for jbin in range(1, h2d.GetNbinsX() + 1):
            content = h2d.GetBinContent(jbin, ibin)
            h2d.SetBinContent(jbin, ibin, content / bin_width)
    h2d.Scale(minimum_bin_width)

def visualise(rdf: RDataFrame, outfile:TDirectory):

    h_centrality = rdf.Histo1D((f'hCentrality', ';Centrality (%)', 100, 0, 100), 'fCentralityFT0C')
    h_pt_he3 = rdf.Histo1D((f'hPtHe3', ';p_{T} He3 (GeV/c)', NBINS_PT, PT_BINNING), 'fPtHe3')
    h_pt_had = rdf.Histo1D((f'hPtHad', ';p_{T} Hadron (GeV/c)', NBINS_PT, PT_BINNING), 'fPtHad')

    outfile.cd()
    h_centrality.Write()
    h_pt_he3.Write()
    h_pt_had.Write()

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
      
    rdf = (rdf.Define('fSignedPtHad', 'fPtHad')
      .Define('fSignHe3', 'fPtHe3/std::abs(fPtHe3)')
      .Redefine('fPtHe3', 'std::abs(fPtHe3)')
      .Redefine('fPtHad', 'std::abs(fPtHad)')
      .Redefine('fPtHe3', '(fPIDtrkHe3 == 7) || (fPIDtrkHe3 == 8) ? fPtHe3 : CorrectPidTrkHe(fPtHe3)')
      .Define('fSignedPtHe3', 'fPtHe3 * fSignHe3')
      .Redefine('fInnerParamTPCHe3', 'fInnerParamTPCHe3 * 2')
      .Define('fClusterSizeCosLamHe3', 'ComputeAverageClusterSize(fItsClusterSizeHe3) / cosh(fEtaHe3)')
      .Define('fClusterSizeCosLamHad', 'ComputeAverageClusterSize(fItsClusterSizeHad) / cosh(fEtaHad)')
      .Define('fNSigmaITSHe3', 'ComputeNsigmaITSHe(fPtHe3 * std::cosh(fEtaHe3), fClusterSizeCosLamHe3)')
      .Define('fNSigmaITSHad', 'ComputeNsigmaITSPr(fPtHad * std::cosh(fEtaHad), fClusterSizeCosLamHad)')
      .Redefine('fNSigmaTPCHe3', 'ComputeNsigmaTPCHe(std::abs(fInnerParamTPCHe3), fSignalTPCHe3, false, true)')
      .Filter(base_selection).Filter(selection))
    
    return rdf

CONFIG = {
        '2023': {
            'input_data': '/data/galucia/lithium/same/LHC23_PbPb_pass5_hadronpid_same.root',
            'tree_name': ['O2he3hadtable', 'O2he3hadmult'],
            'mode': 'DF',
            'parameterisation': {
                    'kHeTPCParams': [-244.627, 0.100002, 1.14056, 0.365353, 3.26339],
                    'kHeTPCResolution': 0.061,
                    'kHeTPCResidualType': 'double_gaus',
                    'kHeTPCParamsResiduals': [87.1178, 0.438363, 0.0814248, -16.9173, 0.754105, 0.285692],

                    'kITSParams_Pr': [1.0228, 1.9634, 2.2081], # Pr from Ka fitting
                    'kITSParams_He': [2.5734, 1.3620, 5.1383],  # He 

                    'kITSResolutionParams_Pr': [0.1575, 0., 0.],
                    'kITSResolutionParams_He': [0.1124, -0.0100, 0.],

                    'kHePidTrkParamsPt': [0.3101, -0.1759, 0.0262],
                }
        },
        '2024': {
            'input_data': '/data/galucia/lithium/same/LHC24ar_pass3_hadronpid_same.root',
            'tree_name': ['O2he3hadtable', 'O2he3hadmult'],
            'mode': 'DF',
            'parameterisation': {
                    'kHeTPCParams': [-200.84, 0.3000, 1.6669, 1.1523, 2.8760],
                    'kHeTPCResolution': 0.059,
                    'kHeTPCResidualType': 'exp',
                    'kHeTPCParamsResiduals': [-162.97, 11.118, 0.4094, 0., 0., 0.],

                    'kITSParams_Pr': [1.0228, 1.9634, 2.2081], # Pr from Ka fitting
                    'kITSParams_He': [2.5734, 1.3620, 5.1383],  # He 

                    'kITSResolutionParams_Pr': [0.1575, 0., 0.],
                    'kITSResolutionParams_He': [0.1124, -0.0100, 0.],

                    'kHePidTrkParamsPt': [0.3101, -0.1759, 0.0262],
            }
        },
        '2025': {
            'input_data': '/data/galucia/lithium/same/LHC25_PbPb_pass1_hadronpid_same.root',
            'tree_name': ['O2he3hadtable', 'O2he3hadmult'],
            'mode': 'DF',
            'parameterisation': {
                    'kHeTPCParams': [-208.928, 0.100000, 1.29761, 0.756458, 3.04256],
                    'kHeTPCResolution': 0.059,
                    'kHeTPCResidualType': 'double_gaus',
                    'kHeTPCParamsResiduals': [200.000, 0.362588, 0.0235587, -32.5300, 0.400000, 0.247265],

                    'kITSParams_Pr': [1.0228, 1.9634, 2.2081], # Pr from Ka fitting
                    'kITSParams_He': [2.5734, 1.3620, 5.1383],  # He 

                    'kITSResolutionParams_Pr': [0.1575, 0., 0.],
                    'kITSResolutionParams_He': [0.1124, -0.0100, 0.],

                    'kHePidTrkParamsPt': [0.3101, -0.1759, 0.0262],
            }
        },
    }

def main(prepare_years: bool = True):
    
    base_selection = '((fSignedPtHe3 > 0 && fSignedPtHad > 0) || (fSignedPtHe3 < 0 && fSignedPtHad < 0))'
    if prepare_years:    
        for year in CONFIG.keys():
            
            load_parametrisation(CONFIG[year])
            chain_data_same, __ = prepare_input_tchain(CONFIG[year])
            rdf_same = prepare_rdataframe(chain_data_same, base_selection, 'true')
            rdf_same.Snapshot(f'rdf_same', f'/data/galucia/lithium/tmp/rdf_same_{year}.root')
        
    rdf = RDataFrame('rdf_same', [f'/data/galucia/lithium/tmp/rdf_same_{year}.root' for year in CONFIG.keys()])

    selections_dict = {'selections': [
                # quality cuts
                 #'((fDeltaPhi/0.006)*(fDeltaPhi/0.006) + (fDeltaEta/0.012)*(fDeltaEta/0.012) > 1)',
                 #'(((fChi2TPCHe3 > 0.5) || (fIs23 == false)) && (fChi2TPCHe3 < 4))', # the cut on chi2tpc should be only applied to 2023
                 '((fChi2TPCHad < 4))',
                 '(std::abs(fEtaHe3) < 0.9)',
                 '(std::abs(fEtaHad) < 0.9)',
                 '(fNClsTPCHe3 > 110)',
                 
                 '((fPtHe3 < 2.5) || (fPIDtrkHe3 == 7))',
                 '((fPIDtrkHe3 == 6) || (fPIDtrkHe3 == 7) || (fPIDtrkHe3 == 8))',
                 #'((fPIDtrkHe3 == 7))',
                 '((-2. < fNSigmaTPCHe3) && (fNSigmaTPCHe3 < 3.))',
                 #'(fClusterSizeCosLamHe3 > 4)',
                 '(fNSigmaITSHe3 > -2.)',
                 #'(fMassTOFHe3 > 2 || fMassTOFHe3 < 0)',

                 '(std::abs(fNSigmaTPCHad) < 2)',
                 #'(std::abs(fNSigmaTPCHadPr) < 2)',
                 '((std::abs(fPtHad) < 0.8) || (std::abs(fNSigmaTOFHad) < 2))',  # TOF hard cut
                 #'((std::abs(fPtHad) < 0.8) || (std::abs(fNSigmaTOFHad) < 2) || (fNSigmaTOFHad < -9.9))', # TOF veto
                 '(fNSigmaITSHad > -2.)',
                ]}

    __, selections = prepare_selections(selections_dict)
    rdf = rdf.Filter(selections)

    filter_code = f"""
    std::unordered_set<std::string> seen_candidates;
    std::atomic<int> discarded_count{{0}};
    std::atomic<int> total_count{{0}};

    bool FilterDuplicates(float z_vertex) {{
        total_count++;

        // Round values to tolerance precision and create hash key
        int z_vertex_key = static_cast<int>(z_vertex / 1e-6);

        std::string key = std::to_string(z_vertex_key);

        if (seen_candidates.count(key) > 0) {{
            discarded_count++;
            return false;  // Duplicate found
        }}

        seen_candidates.insert(key);
        return true;
    }}
    """
    
    gInterpreter.Declare(filter_code)
    rdf = rdf.Filter(f"FilterDuplicates(fZVertex)")
    
    outfile = TFile.Open('centrality_pairs.root', 'RECREATE')
    visualise(rdf, outfile)

    total = int(gROOT.ProcessLine(f"total_count.load();"))
    discarded = int(gROOT.ProcessLine(f"discarded_count.load();"))
    print(f"Total candidates: {total}, Discarded: {discarded}, Kept: {total - discarded}")
        
    outfile.Close()


if __name__ == '__main__':

    main(prepare_years=False)
