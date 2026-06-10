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

NBINS_DCA = 22
DCA_BINNING = np.array([-0.05, -0.04, -0.03, -0.025, -0.02, -0.015, -0.0125, -0.01, -0.0075, -0.005, -0.0025, 0,
                        0.0025, 0.005, 0.0075, 0.01, 0.0125, 0.015, 0.02, 0.025, 0.03, 0.04, 0.05])

    

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

def visualise(rdf: RDataFrame, outfile:TDirectory, particle:str = 'He'):

    particle_suffix = 'He3' if particle == 'He' else 'Had'

    h_pt = rdf.Histo1D((f'hPt{particle}', ';#it{p}_{T} (GeV/#it{c})', 200, -10, 10), f'fSignedPt{particle_suffix}')

    h2_dcaxy_pt = rdf.Histo2D((f'h2DCAxyPt{particle}', ';#it{p}_{T} (GeV/#it{c});DCA_{#it{xy}} (cm)',
                               100, -5, 5, 120, -0.15, 0.15), f'fSignedPt{particle_suffix}', f'fDCAxy{particle_suffix}')
                               #NBINS_PT, PT_BINNING, NBINS_DCA, DCA_BINNING), f'fSignedPt{particle_suffix}', f'fDCAxy{particle_suffix}')
    h2_dcaz_pt = rdf.Histo2D((f'h2DCAzPt{particle}', ';#it{p}_{T} (GeV/#it{c});DCA_{#it{z}} (cm)', 
                               100, -5, 5, 120, -0.3, 0.3), f'fSignedPt{particle_suffix}', f'fDCAz{particle_suffix}')
    
    histos = []
    histos.append(h2_dcaxy_pt)
    histos.append(h2_dcaz_pt)
    histos.append(h_pt)

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
      .Redefine('fInnerParamTPCHe3', '(fPIDtrkHe3 == 7) || (fPIDtrkHe3 == 8) || (fInnerParamTPCHe3 > 2.4) ? fInnerParamTPCHe3 : CorrectPidTrkHe(fInnerParamTPCHe3, false)') \
      .Define('fClusterSizeCosLamHe3', 'ComputeAverageClusterSize(fItsClusterSizeHe3) / cosh(fEtaHe3)') \
      .Define('fClusterSizeCosLamHad', 'ComputeAverageClusterSize(fItsClusterSizeHad) / cosh(fEtaHad)') \
      .Define('fNSigmaITSHe3', 'ComputeNsigmaITSHe(fPtHe3 * std::cosh(fEtaHe3), fClusterSizeCosLamHe3)') \
      .Define('fNSigmaITSHad', 'ComputeNsigmaITSPr(fPtHad * std::cosh(fEtaHad), fClusterSizeCosLamHad)') \
      .Redefine('fNSigmaTPCHe3', 'ComputeNsigmaTPCHe(std::abs(fInnerParamTPCHe3), fSignalTPCHe3, false, true)') \
      .Filter(base_selection).Filter(selection)
    
    return rdf

def main():
    
    config = {
                'input_data': [ #'/data/galucia/lithium/same/LHC23_PbPb_pass5_hadronpid_same.root',
                                '/data/galucia/lithium/same/LHC24ar_pass3_hadronpid_same.root',
                                #'/data/galucia/lithium/same/LHC25_PbPb_pass1_hadronpid_same.root',
                                ],
                'tree_name': 'O2he3hadtable',
                'mode': 'DF',
                
                'parameterisation': {

                    # PbPb 23 pass5
                    ### 'kHeTPCParams': [-335.570, 0.7066, 1.5731, 0.5202, 2.8735],
                    ### 'kHeTPCResolution': 0.060,
                    ### 'kHeTPCParamsResiduals': [31.8927, 0.6830, 0.0750, -8.3483, 0.9692, 0.1287],

                    # 24_pass3
                    'kHeTPCParams': [-200.84, 0.3000, 1.6669, 1.1523, 2.8760],
                    'kHeTPCResolution': 0.059,
                    'kHeTPCResidualType': 'exp',
                    'kHeTPCParamsResiduals': [-162.97, 11.118, 0.4094, 0., 0., 0.],

                    # 25_pass1
                    ### 'kHeTPCParams': [-444.875, 1.03555, 1.79299, 1.49951, 2.26467],
                    ### 'kHeTPCResolution': 0.055,
                    ### 'kHeTPCResidualType': 'double_gaus',
                    ### 'kHeTPCParamsResiduals': [417.113, 0.532354, 0.0620908, -43.7363, 0.710629, 0.286426],

                    # PbPb 23 pass5
                    'kITSParams_Pr': [1.0228, 1.9634, 2.2081], # Pr from Ka fitting
                    'kITSParams_He': [2.5734, 1.3620, 5.1383],  # He 

                    'kITSResolutionParams_Pr': [0.1575, 0., 0.],
                    'kITSResolutionParams_He': [0.1124, -0.0100, 0.],

                    'kHePidTrkParamsPt': [0.3101, -0.1759, 0.0262],
                    
                    # from LHC25g11
                    'kHePidTrkParamsP': [0.0514009, -0.245084, 0.120208],
                    ####'kHePidTrkParamsP': [0., 0., 0.],
                }
             }

    chain_data = prepare_input_tchain(config)
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
    
    load_parametrisation(config)  # Load parametrisation into Common.h
    base_selections, selections = prepare_selections(selections_dict)
    rdf = prepare_rdataframe(chain_data, base_selections, selections)

    outfile = TFile.Open('output/dca_data_template.root', 'RECREATE')
    
    for particle_name, particle in zip(['Had', 'He3'], ['Pr', 'He']):

        #filter_code = f"""
        #std::vector<std::tuple<float, float, float>> seen_candidates_{particle_name};
        #
        #bool FilterDuplicates_{particle_name}(float pt, float eta, float phi) {{
        #    for (const auto& cand : seen_candidates_{particle_name}) {{
        #        if (std::abs(std::get<0>(cand) - pt) < 1e-6 &&
        #            std::abs(std::get<1>(cand) - eta) < 1e-7 &&
        #            std::abs(std::get<2>(cand) - phi) < 1e-7) {{
        #            return false;  // All three match, exclude duplicate
        #        }}
        #    }}
        #    seen_candidates_{particle_name}.push_back(std::make_tuple(pt, eta, phi));
        #    return true;
        #}}
        #"""

        filter_code = f"""
        std::unordered_set<std::string> seen_candidates_{particle_name};
        std::atomic<int> discarded_count_{particle_name}{{0}};
        std::atomic<int> total_count_{particle_name}{{0}};

        bool FilterDuplicates_{particle_name}(float pt, float eta, float phi) {{
            total_count_{particle_name}++;

            // Round values to tolerance precision and create hash key
            int pt_key = static_cast<int>(pt / 1e-6);
            int eta_key = static_cast<int>(eta / 1e-6);
            int phi_key = static_cast<int>(phi / 1e-6);

            std::string key = std::to_string(pt_key) + "_" + 
                             std::to_string(eta_key) + "_" + 
                             std::to_string(phi_key);

            if (seen_candidates_{particle_name}.count(key) > 0) {{
                discarded_count_{particle_name}++;
                return false;  // Duplicate found
            }}

            seen_candidates_{particle_name}.insert(key);
            return true;
        }}
        """
        
        gInterpreter.Declare(filter_code)
        
        rdf_particle = rdf.Filter(f"FilterDuplicates_{particle_name}(fPt{particle_name}, fEta{particle_name}, fPhi{particle_name})")


        #rdf_particle = filter_duplicates(rdf, f'fPt{particle_name}', 1e-6)
        #rdf_particle = filter_duplicates(rdf_particle, f'fEta{particle_name}', 1e-7)
        #rdf_particle = filter_duplicates(rdf_particle, f'fPhi{particle_name}', 1e-7)

        outdir = outfile.mkdir(f'{particle}')
        visualise(rdf_particle, outdir, particle)

        total = int(gROOT.ProcessLine(f"total_count_{particle_name}.load();"))
        discarded = int(gROOT.ProcessLine(f"discarded_count_{particle_name}.load();"))
        print(f"{tc.CYAN}{particle_name}{tc.RESET}: Total candidates: {total}, Discarded: {discarded}, Kept: {total - discarded}")
        
    outfile.Close()


if __name__ == '__main__':

    main()
