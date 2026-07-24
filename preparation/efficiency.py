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

from include.load_parameters import load_parametrisation

gInterpreter.ProcessLine(f'#include "../include/Common.h"')
ROOT.EnableImplicitMT(10)
ROOT.gROOT.SetBatch(True)


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
    print(tc.GREEN+'nDataset columns'+tc.RESET)
    print(tc.UNDERLINE+tc.CYAN+f'{rdf.GetColumnNames()}'+tc.RESET)
    
    # TPC
    if 'fNSigmaTPCHadPr' in rdf.GetColumnNames():
        rdf = rdf.Define('fNSigmaTPCHad', 'fNSigmaTPCHadPr')
    else:
        rdf = rdf.Define('fNSigmaTPCHadPr', 'fNSigmaTPCHad')
    
    # TOF
    if 'fNSigmaTOFHad' not in rdf.GetColumnNames() and 'fNSigmaTOFHadPr' not in rdf.GetColumnNames():
        rdf = rdf.Define('fNSigmaTOFHad', 'ComputeNsigmaTOFPr(std::abs(fPtHad), fMassTOFHad)')
    elif 'fNSigmaTOFHad' not in rdf.GetColumnNames() and 'fNSigmaTOFHadPr' in rdf.GetColumnNames():
       rdf = rdf.Define('fNSigmaTOFHad', 'fNSigmaTOFHadPr')
    
      # Recalibration
      # Correct for PID in tracking
    
    rdf_gen = (rdf.Define('fSignedPtHe3', 'fPtHe3') 
        .Define('fSignedPtHad', 'fPtHad') 
        .Define('fSignedPtHe3MC', 'fPtMCHe3') 
        .Define('fSignedPtHadMC', 'fPtMCHad') 
        .Filter(base_selection) 
        .Define('fSignHe3', 'fPtMCHe3/std::abs(fPtMCHe3)') 
        .Define('fPtLiMC', 'std::abs(fSignedPtMC)')
    )
      
      
    rdf_rec = (rdf.Define('fSignedPtHad', 'fPtHad') 
      .Define('fSignHe3', 'fPtHe3/std::abs(fPtHe3)') 
      .Redefine('fPtHe3', 'std::abs(fPtHe3)') 
      .Redefine('fPtHe3', '(fPIDtrkHe3 == 7) || (fPIDtrkHe3 == 8) || (fPtHe3 > 2.5) ? fPtHe3 : CorrectPidTrkHe(fPtHe3)') 
      .Redefine('fPtHad', 'std::abs(fPtHad)') 
      .Define('fSignedPtHe3', 'fPtHe3 * fSignHe3') 
      .Define(f'fEHe3', f'std::sqrt((fPtHe3 * std::cosh(fEtaHe3))*(fPtHe3 * std::cosh(fEtaHe3)) + {ParticleMasses["He"]}*{ParticleMasses["He"]})') 
      .Define(f'fEHad', f'std::sqrt((fPtHad * std::cosh(fEtaHad))*(fPtHad * std::cosh(fEtaHad)) + {ParticleMasses["Pr"]}*{ParticleMasses["Pr"]})') 
      .Define('fDeltaEta', 'fEtaHe3 - fEtaHad') 
      .Define('fDeltaPhi', 'fPhiHe3 - fPhiHad') 
      #.Redefine('fInnerParamTPCHe3', 'fInnerParamTPCHe3 * 2') 
      #.Redefine('fInnerParamTPCHe3', '(fPIDtrkHe3 == 7) || (fPIDtrkHe3 == 8) || (fPtHe3 > 2.5) ? fInnerParamTPCHe3 : CorrectPidTrkHe(fInnerParamTPCHe3, false)') 
      .Redefine('fNSigmaTPCHe3', 'ComputeNsigmaTPCHe(std::abs(fInnerParamTPCHe3), fSignalTPCHe3, true, true)') 
      .Define('fClusterSizeCosLamHe3', 'ComputeAverageClusterSize(fItsClusterSizeHe3) / cosh(fEtaHe3)') 
      .Define('fClusterSizeCosLamHad', 'ComputeAverageClusterSize(fItsClusterSizeHad) / cosh(fEtaHad)') 
      .Define('fExpectedClusterSizeHe3', 'ComputeExpectedClusterSizeCosLambdaHe(fPtHe3 * std::cosh(fEtaHe3))') 
      .Define('fExpectedClusterSizeHad', 'ComputeExpectedClusterSizeCosLambdaPr(fPtHad * std::cosh(fEtaHad))') 
      .Define('fNSigmaITSHe3', 'ComputeNsigmaITSHe(fPtHe3 * std::cosh(fEtaHe3), fClusterSizeCosLamHe3, true)') 
      .Define('fNSigmaITSHad', 'ComputeNsigmaITSPr(fPtHad * std::cosh(fEtaHad), fClusterSizeCosLamHad, true)') 
      .Define('fNSigmaDCAxyHe3', 'ComputeNsigmaDCAxyHe(fPtHe3, fDCAxyHe3, true)') 
      .Define('fNSigmaDCAzHe3', 'ComputeNsigmaDCAzHe(fPtHe3, fDCAzHe3, true)') 
      .Define('fNSigmaDCAxyHad', 'ComputeNsigmaDCAxyPr(fPtHad, fDCAxyHad, true)') 
      .Define('fNSigmaDCAzHad', 'ComputeNsigmaDCAzPr(fPtHad, fDCAzHad, true)') 
      .Filter(base_selection).Filter(selection) 
      .Define('fPxLi', 'fPtHe3 * std::cos(fPhiHe3) + fPtHad * std::cos(fPhiHad)') 
      .Define('fPyLi', 'fPtHe3 * std::sin(fPhiHe3) + fPtHad * std::sin(fPhiHad)') 
      .Define('fPzLi', 'fPtHe3 * std::sinh(fEtaHe3) + fPtHad * std::sinh(fEtaHad)')
      .Define('fPtLi', 'sqrt(fPxLi*fPxLi + fPyLi*fPyLi)')
      .Define('fSignedPtLi', 'fPtLi * fSignHe3') 
      .Define('fELi', 'fEHe3 + fEHad')
      .Define('fPLi', 'sqrt(fPxLi*fPxLi + fPyLi*fPyLi + fPzLi*fPzLi)')
      .Define('fMassInvLi', 'std::sqrt(fELi*fELi - fPLi*fPLi)')
    )
      
    return rdf_rec, rdf_gen

def visualise(rdf, output_file: TFile, mode:str='rec'):
   
    print(f'n{tc.GREEN}Creating histograms{tc.RESET}')
    histogram_registry = HistogramRegistry()

    if mode == 'rec':
        register_qa_histograms(histogram_registry)

    if mode == 'gen':
        gen_qa_histograms_archive = {key: value for key, value in Archive.QA_HISTOGRAMS.items()}
        for name, entry in gen_qa_histograms_archive.items():
            if 'Pt' in entry.xvar:    entry.xvar += 'MC'
            histogram_registry.register(entry)
    
    print(f'tQA Histograms created!')
        
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
            efficiency_error = np.sqrt( efficiency * (1 - efficiency) / hist_gen.GetBinContent(ibin)) if efficiency < 1 else 0
            efficiency_hist.SetBinContent(ibin, efficiency)
            efficiency_hist.SetBinError(ibin, efficiency_error)
        else:
            efficiency_hist.SetBinContent(ibin, 0)
            efficiency_hist.SetBinError(ibin, 0)

    efficiency_hist.GetYaxis().SetTitle('Efficiency')
    return efficiency_hist

def efficiency_histograms(rdf_rec, rdf_gen, output_file: TFile):

    pt_bins, pt_min, pt_max = 20, 0, 10
    hPtLiMatter = rdf_rec.Filter('fSignHe3 > 0').Histo1D(("hPtLiMatter", ";#it{p}_{T} (^{4}Li) (GeV/#it{c});", pt_bins, pt_min, pt_max), "fPtLi").GetValue()
    hPtLiAntimatter = rdf_rec.Filter('fSignHe3 < 0').Histo1D(("hPtLiAntimatter", ";#it{p}_{T} (^{4}Li) (GeV/#it{c});", pt_bins, pt_min, pt_max), "fPtLi").GetValue()
    hPtLiMCMatter = rdf_gen.Filter('fSignHe3 > 0').Histo1D(("hPtLiMCMatter", ";#it{p}_{T} (^{4}#bar{Li}) (GeV/#it{c});", pt_bins, pt_min, pt_max), "fPtLiMC").GetValue()
    hPtLiMCAntimatter = rdf_gen.Filter('fSignHe3 < 0').Histo1D(("hPtLiMCAntimatter", ";#it{p}_{T} (^{4}#bar{Li}) (GeV/#it{c});", pt_bins, pt_min, pt_max), "fPtLiMC").GetValue()
    
    hPtLi = rdf_rec.Histo1D(("hPtLi", ";#it{p}_{T} (^{4}#bar{Li}) (GeV/#it{c});", 2*pt_bins, -pt_max, pt_max), "fSignedPtLi").GetValue()
    hPtLi010 = rdf_rec.Filter('fCentralityFT0C < 10').Histo1D(("hPtLi010", "0-10% FT0C;#it{p}_{T} (^{4}#bar{Li}) (GeV/#it{c});", 2*pt_bins, -pt_max, pt_max), "fSignedPtLi").GetValue()
    hPtLi1050 = rdf_rec.Filter('10 < fCentralityFT0C && fCentralityFT0C < 50').Histo1D(("hPtLi1050", "10-50% FT0C;#it{p}_{T} (^{4}#bar{Li}) (GeV/#it{c});", 2*pt_bins, -pt_max, pt_max), "fSignedPtLi").GetValue()
    hPtLiMC = rdf_gen.Histo1D(("hPtLiMC", ";#it{p}_{T} (^{4}#bar{Li}) (GeV/#it{c});", 2*pt_bins, -pt_max, pt_max), "fSignedPtMC").GetValue()
    hPtLiMC010 = rdf_gen.Filter('fCentralityFT0C < 10').Histo1D(("hPtLiMC010", "0-10% FT0C;#it{p}_{T} (^{4}#bar{Li}) (GeV/#it{c});", 2*pt_bins, -pt_max, pt_max), "fSignedPtMC").GetValue()
    hPtLiMC1050 = rdf_gen.Filter('10 < fCentralityFT0C && fCentralityFT0C < 50').Histo1D(("hPtLiMC1050", "10-50% FT0C;#it{p}_{T} (^{4}#bar{Li}) (GeV/#it{c});", 2*pt_bins, -pt_max, pt_max), "fSignedPtMC").GetValue()

    hEfficiencyMatter = produce_efficiency_histogram(hPtLiMatter, hPtLiMCMatter, "hEfficiencyMatter")
    hEfficiencyAntimatter = produce_efficiency_histogram(hPtLiAntimatter, hPtLiMCAntimatter, "hEfficiencyAntimatter")
    hEfficiency = produce_efficiency_histogram(hPtLi, hPtLiMC, "hEfficiency")
    hEfficiency010 = produce_efficiency_histogram(hPtLi010, hPtLiMC010, "hEfficiency010")
    hEfficiency1050 = produce_efficiency_histogram(hPtLi1050, hPtLiMC1050, "hEfficiency1050")
    

    output_file.mkdir('Efficiency')
    output_file.cd('Efficiency')
    for hist in [hPtLiMatter, hPtLiAntimatter, hPtLiMCMatter, hPtLiMCAntimatter, hPtLi, hPtLi010, hPtLi1050, hPtLiMC, hPtLiMC010, hPtLiMC1050,
                 hEfficiencyMatter, hEfficiencyAntimatter, hEfficiency, hEfficiency010, hEfficiency1050]:
        hist.Write()


if __name__ == '__main__':

    config_file = 'config/config_efficiency.yml'
    config = yaml.safe_load(open(config_file, 'r'))

    load_parametrisation(config)  # Load parametrisation into Common.h

    base_selection, selection = prepare_selections(config)
    chain_data, additional_chain_data = prepare_input_tchain(config)
    rdf_rec, rdf_gen = prepare_rdataframe(chain_data, base_selection, selection)

    output_file_path = config['output_file']
    output_file = ROOT.TFile(output_file_path, "RECREATE")
    output_dir_rec = output_file.mkdir('Reconstructed')
    output_dir_gen = output_file.mkdir('Generated')

    visualise(rdf_rec, output_dir_rec, mode='rec')
    visualise(rdf_gen, output_dir_gen, mode='gen')

    efficiency_histograms(rdf_rec, rdf_gen, output_file)

    output_file.Close()

  