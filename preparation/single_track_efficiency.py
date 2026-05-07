import sys
import yaml
import ROOT
import numpy as np
from ROOT import TFile, TChain, gInterpreter, RDataFrame

from torchic.utils.terminal_colors import TerminalColors as tc

sys.path.append('..')
from utils.particles import ParticleMasses
from utils.histogram_registry import HistogramRegistry
from utils.histogram_archive import register_qa_histograms, Archive

gInterpreter.ProcessLine(f'#include "../include/Common.h"')
ROOT.EnableImplicitMT(10)
ROOT.gROOT.SetBatch(True)

def declare_duplicate_filter(particle_name, tolerance=1e-6):
    filter_code = f"""
            std::unordered_set<std::string> seen_candidates_{particle_name};
            std::atomic<int> discarded_count_{particle_name}{{0}};
            std::atomic<int> total_count_{particle_name}{{0}};

            bool FilterDuplicates_{particle_name}(float pt, float eta, float phi) {{
                total_count_{particle_name}++;

                // Round values to tolerance precision and create hash key
                int pt_key = static_cast<int>(pt / {tolerance});
                int eta_key = static_cast<int>(eta / {tolerance});
                int phi_key = static_cast<int>(phi / {tolerance});

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
    else:
        rdf = rdf.Define('fNSigmaTPCHadPr', 'fNSigmaTPCHad')
    
    # TOF
    if 'fNSigmaTOFHadPr' not in rdf.GetColumnNames():
        rdf = rdf.Define('fNSigmaTOFHad', 'ComputeNsigmaTOFPr(std::abs(fPtHad), fMassTOFHad)')
    else:
       rdf = rdf.Define('fNSigmaTOFHad', 'fNSigmaTOFHadPr')
    
      # Recalibration
      # Correct for PID in tracking
    
    rdf_gen = rdf.Define('fSignedPtHe3', 'fPtHe3') \
        .Define('fSignedPtHad', 'fPtHad') \
        .Define('fSignedPtHe3MC', 'fPtMCHe3') \
        .Define('fSignedPtHadMC', 'fPtMCHad') \
        .Redefine('fPtMCHe3', 'std::abs(fPtMCHe3)') \
        .Redefine('fPtMCHad', 'std::abs(fPtMCHad)') \
        .Filter(base_selection) \
        .Define('fSignHe3', 'fPtMCHe3/std::abs(fPtMCHe3)') \
      
    rdf_rec = rdf.Define('fSignedPtHad', 'fPtHad') \
      .Define('fSignHe3', 'fPtHe3/std::abs(fPtHe3)') \
      .Redefine('fPtHe3', 'std::abs(fPtHe3)') \
      .Redefine('fPtHad', 'std::abs(fPtHad)') \
      .Redefine('fPtMCHe3', 'std::abs(fPtMCHe3)') \
      .Redefine('fPtMCHad', 'std::abs(fPtMCHad)') \
      .Redefine('fPtHe3', '(fPIDtrkHe3 == 7) || (fPIDtrkHe3 == 8) || (fPtHe3 > 2.5) ? fPtHe3 : CorrectPidTrkHe(fPtHe3)') \
      .Define('fSignedPtHe3', 'fPtHe3 * fSignHe3') \
      .Define(f'fEHe3', f'std::sqrt((fPtHe3 * std::cosh(fEtaHe3))*(fPtHe3 * std::cosh(fEtaHe3)) + {ParticleMasses["He"]}*{ParticleMasses["He"]})') \
      .Define(f'fEHad', f'std::sqrt((fPtHad * std::cosh(fEtaHad))*(fPtHad * std::cosh(fEtaHad)) + {ParticleMasses["Pr"]}*{ParticleMasses["Pr"]})') \
      .Redefine('fInnerParamTPCHe3', 'fInnerParamTPCHe3 * 2') \
      .Redefine('fNSigmaTPCHe3', 'ComputeNsigmaTPCHe(std::abs(fInnerParamTPCHe3), fSignalTPCHe3, true)') \
      .Define('fClusterSizeCosLamHe3', 'ComputeAverageClusterSize(fItsClusterSizeHe3) / cosh(fEtaHe3)') \
      .Define('fClusterSizeCosLamHad', 'ComputeAverageClusterSize(fItsClusterSizeHad) / cosh(fEtaHad)') \
      .Define('fNSigmaITSHe3', 'ComputeNsigmaITSHe(fPtHe3 * std::cosh(fEtaHe3), fClusterSizeCosLamHe3, true)') \
      .Define('fNSigmaITSHad', 'ComputeNsigmaITSPr(fPtHad * std::cosh(fEtaHad), fClusterSizeCosLamHad, true)') \
      .Define('fNSigmaDCAxyHe3', 'ComputeNsigmaDCAxyHe(fPtHe3, fDCAxyHe3, true)') \
      .Define('fNSigmaDCAzHe3', 'ComputeNsigmaDCAzHe(fPtHe3, fDCAzHe3, true)') \
      .Define('fNSigmaDCAxyHad', 'ComputeNsigmaDCAxyPr(fPtHad, fDCAxyHad, true)') \
      .Define('fNSigmaDCAzHad', 'ComputeNsigmaDCAzPr(fPtHad, fDCAzHad, true)') \
      .Filter(base_selection).Filter(selection)
    
    return rdf_rec, rdf_gen

def visualise(rdf, output_file: TFile, mode:str='rec'):
   
    print(f'\n{tc.GREEN}Creating histograms{tc.RESET}')
    histogram_registry = HistogramRegistry()

    if mode == 'rec':
        register_qa_histograms(histogram_registry)

    if mode == 'gen':
        gen_qa_histograms_archive = {key: value for key, value in Archive.QA_HISTOGRAMS.items()}
        for name, entry in gen_qa_histograms_archive.items():
            if 'Pt' in entry.xvar:    entry.xvar += 'MC'
            histogram_registry.register(entry)
    
    print(f'\tQA Histograms created!')
        
    histogram_registry.prepare_directories(output_file)
    histogram_registry.draw_histogram(rdf)
    histogram_registry.save_histograms(output_file)

def produce_efficiency_histogram(hist_rec, hist_gen, name):
    """
    Calculate the efficiency histogram by dividing the reconstructed histogram by the generated histogram.
    """
    efficiency_hist = hist_rec.Clone(name)

    for ibin in range(1, hist_rec.GetNbinsX() + 1):
        if hist_gen.GetBinContent(ibin) > 0:
            efficiency = hist_rec.GetBinContent(ibin) / hist_gen.GetBinContent(ibin)
            if efficiency > 1:
                print(f'Warning: Efficiency greater than 1 in bin {ibin}. Setting efficiency to 1.')
                efficiency = 1
            efficiency_error = np.sqrt( efficiency * (1 - efficiency) / hist_gen.GetBinContent(ibin))
            efficiency_hist.SetBinContent(ibin, efficiency)
            efficiency_hist.SetBinError(ibin, efficiency_error)
        else:
            efficiency_hist.SetBinContent(ibin, 0)
            efficiency_hist.SetBinError(ibin, 0)

    return efficiency_hist

def efficiency_histograms(rdf_rec, rdf_gen, output_file: TFile, particle: str):

    pt_bins, pt_min, pt_max = 100, 0, 10
    hPtMatter = rdf_rec.Filter(f'fSignedPt{particle} > 0').Histo1D((f"hPt{particle}Matter", f";#it{{p}}_{{T}} ({particle}) (GeV/#it{{c}});", pt_bins, pt_min, pt_max), f"fPt{particle}").GetValue()
    hPtAntimatter = rdf_rec.Filter(f'fSignedPt{particle} < 0').Histo1D((f"hPt{particle}Antimatter", f";#it{{p}}_{{T}} ({particle}) (GeV/#it{{c}});", pt_bins, pt_min, pt_max), f"fPt{particle}").GetValue()
    hPtMCMatter = rdf_gen.Filter(f'fSignedPt{particle}MC > 0').Histo1D((f"hPt{particle}MCMatter", f";#it{{p}}_{{T}} (anti{particle}) (GeV/#it{{c}});", pt_bins, pt_min, pt_max), f"fPtMC{particle}").GetValue()
    hPtMCAntimatter = rdf_gen.Filter(f'fSignedPt{particle}MC < 0').Histo1D((f"hPt{particle}MCAntimatter", f";#it{{p}}_{{T}} (anti{particle}) (GeV/#it{{c}});", pt_bins, pt_min, pt_max), f"fPtMC{particle}").GetValue()

    hEfficiencyMatter = produce_efficiency_histogram(hPtMatter, hPtMCMatter, f"hEfficiency{particle}Matter")
    hEfficiencyAntimatter = produce_efficiency_histogram(hPtAntimatter, hPtMCAntimatter, f"hEfficiency{particle}Antimatter")

    output_file.cd()
    for hist in [hPtMatter, hPtAntimatter, hPtMCMatter, hPtMCAntimatter,
                 hEfficiencyMatter, hEfficiencyAntimatter]:
        hist.Write()


if __name__ == '__main__':

    config_file = 'config/config_efficiency.yml'
    config = yaml.safe_load(open(config_file, 'r'))

    base_selection, selection = prepare_selections(config)
    chain_data = prepare_input_tchain(config)
    rdf_rec, rdf_gen = prepare_rdataframe(chain_data, base_selection, selection)

    declare_duplicate_filter('Had')
    rdf_rec_protons = rdf_rec.Filter('FilterDuplicates_Had(fPtMCHad, fEtaMCHad, fPhiMCHad)')
    rdf_gen_protons = rdf_gen.Filter('FilterDuplicates_Had(fPtMCHad, fEtaMCHad, fPhiMCHad)')

    declare_duplicate_filter('He3')
    rdf_rec_he3 = rdf_rec.Filter('FilterDuplicates_He3(fPtMCHe3, fEtaMCHe3, fPhiMCHe3)')
    rdf_gen_he3 = rdf_gen.Filter('FilterDuplicates_He3(fPtMCHe3, fEtaMCHe3, fPhiMCHe3)')

    output_file_path = 'output/single_track_efficiency_tof_hit.root'
    output_file = ROOT.TFile(output_file_path, "RECREATE")
    output_dir_rec = output_file.mkdir('Reconstructed')
    output_dir_gen = output_file.mkdir('Generated')

    #visualise(rdf_rec_protons, output_dir_rec, mode='rec')
    #visualise(rdf_gen_protons, output_dir_gen, mode='gen')

    efficiency_histograms(rdf_rec, rdf_gen, output_file, 'Had')
    efficiency_histograms(rdf_rec, rdf_gen, output_file, 'He3')

    output_file.Close()

  