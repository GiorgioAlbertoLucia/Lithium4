import sys
import yaml
import ROOT
from ROOT import TFile, TChain

from torchic.utils.terminal_colors import TerminalColors as tc

sys.path.append('..')
from utils.particles import ParticleMasses
from utils.histogram_registry import HistogramRegistry
from utils.histogram_archive import register_qa_histograms, register_kstar_histograms, register_kstar_matter_histograms, \
    register_kstar_antimatter_histograms

ROOT.gROOT.LoadMacro('../include/Common.h++')
from ROOT import NsigmaTpcHe, NsigmaITSHe, NsigmaITSPr, NsigmaTOFPr, averageClusterSize, CorrectPidTrkHe, \
  Kstar, ExpectedClusterSizeCosLambdaHe, ExpectedClusterSizeCosLambdaPr

ROOT.EnableImplicitMT(30)
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
   
    rdf = ROOT.ROOT.RDataFrame(chain_data) \
      .Define('fSignedPtHad', 'fPtHad') \
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
      .Redefine('fNSigmaTPCHe3', 'NsigmaTpcHe(std::abs(fInnerParamTPCHe3), fSignalTPCHe3)') \
      .Define('fNSigmaTOFHad', 'NsigmaTOFPr(fPtHad * std::cosh(fEtaHad), fMassTOFHad)') \
      .Define('fClusterSizeCosLamHe3', 'averageClusterSize(fItsClusterSizeHe3) / cosh(fEtaHe3)') \
      .Define('fClusterSizeCosLamHad', 'averageClusterSize(fItsClusterSizeHad) / cosh(fEtaHad)') \
      .Define('fExpectedClusterSizeHe3', 'ExpectedClusterSizeCosLambdaHe(fPtHe3 * std::cosh(fEtaHe3))') \
      .Define('fExpectedClusterSizeHad', 'ExpectedClusterSizeCosLambdaPr(fPtHad * std::cosh(fEtaHad))') \
      .Define('fNSigmaITSHe3', 'NsigmaITSHe(fPtHe3 * std::cosh(fEtaHe3), fClusterSizeCosLamHe3)') \
      .Define('fNSigmaITSHad', 'NsigmaITSPr(fPtHad * std::cosh(fEtaHad), fClusterSizeCosLamHad)') \
      .Filter(base_selection).Filter(selection) \
      .Define('fKstar', f'Kstar(fPtHe3, fEtaHe3, fPhiHe3, {ParticleMasses["He"]}, fPtHad, fEtaHad, fPhiHad, {ParticleMasses["Pr"]})') \
      .Define('fPxLi', 'fPtHe3 * std::cos(fPhiHe3) + fPtHad * std::cos(fPhiHad)') \
      .Define('fPyLi', 'fPtHe3 * std::sin(fPhiHe3) + fPtHad * std::sin(fPhiHad)') \
      .Define('fPzLi', 'fPtHe3 * std::sinh(fEtaHe3) + fPtHad * std::sinh(fEtaHad)') \
      .Define('fELi', 'fEHe3 + fEHad') \
      .Define('fPLi', 'sqrt(fPxLi*fPxLi + fPyLi*fPyLi + fPzLi*fPzLi)') \
      .Define('fPtLi', 'sqrt(fPxLi*fPxLi + fPyLi*fPyLi)') \
      .Define('fKt', 'fPtLi / 2') \
      .Define('fMt', f'std::sqrt((fPtLi/4)*(fPtLi/4) + {ParticleMasses["Pr"]}*{ParticleMasses["Pr"]})') \
      .Define('fEtaLi', 'std::acosh(fPLi / fELi)') \
      .Define('fPhiLi', 'std::atan2(fPyLi, fPxLi)') \
      .Define('fSignedPtLi', 'fPtLi * fSignHe3') \
      .Define('fMassInvLi', 'std::sqrt(fELi*fELi - fPLi*fPLi)') \
      .Define('fMassTLi', 'std::sqrt(fELi*fELi - fPtLi*fPtLi)') \
    
    return rdf

def visualise(rdf, output_file: TFile):
   
    print(f'\n{tc.GREEN}Creating histograms{tc.RESET}')
    histogram_registry = HistogramRegistry()

    register_qa_histograms(histogram_registry)
    print(f'\tQA Histograms created!')

    register_kstar_histograms(histogram_registry)
    print(f'\t(anti)matter histograms created!')

    register_kstar_matter_histograms(histogram_registry)
    print(f'\tmatter histograms created!')

    register_kstar_antimatter_histograms(histogram_registry)
    print(f'\tantimatter histograms created!')

    histogram_registry.prepare_directories(output_file)
    histogram_registry.draw_histogram(rdf)
    histogram_registry.save_histograms(output_file)



if __name__ == '__main__':

    config_file = 'config/config_prepare.yml'
    config = yaml.safe_load(open(config_file, 'r'))

    base_selection, selection = prepare_selections(config)
    chain_data = prepare_input_tchain(config)
    rdf = prepare_rdataframe(chain_data, base_selection, selection)

    output_file_path = config['output_file']
    output_file = ROOT.TFile(output_file_path, "RECREATE")

    visualise(rdf, output_file)    

    output_file.Close()

  