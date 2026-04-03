import sys
import yaml
import ROOT
from ROOT import TFile, TChain, gInterpreter, RDataFrame

from torchic.core.graph import load_graph
from torchic.core.histogram import load_hist, hist_to_graph
from torchic.utils.terminal_colors import TerminalColors as tc

sys.path.append('..')
from utils.particles import ParticleMasses

gInterpreter.ProcessLine(f'#include "../include/Common.h"')
gInterpreter.ProcessLine(f'#include "../include/LambdaUtils.h"')

ROOT.EnableImplicitMT(10)
ROOT.gROOT.SetBatch(True)


def prepare_selections(config):

    base_selection = ''
    if config.get('like_sign', True):  
        base_selection = '(fSignedPtHe3 > 0 && fSignedPtHad > 0) || (fSignedPtHe3 < 0 && fSignedPtHad < 0)'
    else:
        base_selection = '(fSignedPtHe3 > 0 && fSignedPtHad < 0) || (fSignedPtHe3 < 0 && fSignedPtHad > 0)'

    purity_primary_selection = '(fPtHad > 0.4) && (fPtHad < 4.0) && (fPtHe3 > 1.6) && (fPtHe3 < 4.0)'

    selections = config['selections']
    
    selection = selections[0]
    for sel in selections[1:]:
        selection += (' && ' + sel)

    return base_selection, selection, purity_primary_selection

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

def prepare_rdataframe(chain_data: TChain, base_selection: str, selection: str, purity_primary_selection: str='') -> RDataFrame:
   
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
      .Filter(base_selection).Filter(selection).Filter(purity_primary_selection) \
      .Define('fKstar', f'ComputeKstar(fPtHe3, fEtaHe3, fPhiHe3, {ParticleMasses["He"]}, fPtHad, fEtaHad, fPhiHad, {ParticleMasses["Pr"]})') \
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

def purity_and_primary_fraction_histograms(rdf:RDataFrame, outfile:TFile):
   
    graphs = {
        'Matter': {},
        'Antimatter': {}
    }

    #input_primary_fraction_file = '/home/galucia/Lithium4/calibration/output/dca/primary_fraction_results_new.root'
    input_primary_fraction_file = '/home/galucia/Lithium4/calibration/output/dca/primary_fraction_results_new_smaller_tolerance_pr.root'

    graphs['Matter']['primary_he3'] = load_graph(input_primary_fraction_file, 'He/DCAxy/g_primary_fraction_matter')
    graphs['Matter']['primary_p'] = load_graph(input_primary_fraction_file, 'Pr/DCAxy/g_primary_fraction_matter')
    graphs['Matter']['secondary_p_from_weak_decay'] = load_graph(input_primary_fraction_file, 'Pr/DCAxy/g_weak_decay_fraction_matter')
    graphs['Matter']['purity_he3'] = load_graph('/home/galucia/Lithium4/calibration/output/purity/purity_23_pass4_he3.root',
                                    'He3/TPC/g_purity_matter')
    graphs['Matter']['purity_p_tpc'] = load_graph('/home/galucia/Lithium4/calibration/output/purity/purity_23_pass4_p.root',
                                    'Had/TPC/g_purity_matter')
    graphs['Matter']['purity_p_tof'] = load_graph('/home/galucia/Lithium4/calibration/output/purity/purity_23_pass4_p.root',
                                    'Had/TPC/g_purity_matter')
    graphs['Antimatter']['primary_he3'] = load_graph(input_primary_fraction_file, 'He/DCAxy/g_primary_fraction_antimatter')
    graphs['Antimatter']['primary_p'] = load_graph(input_primary_fraction_file, 'Pr/DCAxy/g_primary_fraction_antimatter')
    graphs['Antimatter']['secondary_p_from_weak_decay'] = load_graph(input_primary_fraction_file, 'Pr/DCAxy/g_weak_decay_fraction_antimatter')
    graphs['Antimatter']['purity_he3'] = load_graph('/home/galucia/Lithium4/calibration/output/purity/purity_23_pass4_he3.root',
                                    'He3/TPC/g_purity_antimatter')
    graphs['Antimatter']['purity_p_tpc'] = load_graph('/home/galucia/Lithium4/calibration/output/purity/purity_23_pass4_p.root',
                                    'Had/TPC/g_purity_antimatter')
    graphs['Antimatter']['purity_p_tof'] = load_graph('/home/galucia/Lithium4/calibration/output/purity/purity_23_pass4_p.root',
                                    'Had/TPC/g_purity_antimatter')
    
    h_protons_from_sigma = load_hist('/home/galucia/Lithium4/calibration/output/dca/dca_mc_template_check.root',
                                     'Pr/hPtPr_IsFromSigmaPlus')
    h_protons_from_lambda = load_hist('/home/galucia/Lithium4/calibration/output/dca/dca_mc_template_check.root',
                                     'Pr/hPtPr_IsFromLambda0')
    h_sigma_protons_to_weak_decay_protons = h_protons_from_sigma.Clone('h_sigma_to_weak_decay_protons_ratio')
    
    h_weak_decay_protons = h_protons_from_lambda.Clone('h_weak_decay_protons')
    h_weak_decay_protons.Add(h_protons_from_sigma)

    h_sigma_protons_to_weak_decay_protons.Divide(h_weak_decay_protons)
    g_sigma_protons_to_weak_decay_protons = hist_to_graph(h_sigma_protons_to_weak_decay_protons, 'g_sigma_to_weak_decay_protons_ratio')
    
    ROOT.gInterpreter.Declare('''
                              TGraph * graph_purity_he3_matter = nullptr;
                              TGraph * graph_purity_p_tpc_matter = nullptr;
                              TGraph * graph_purity_p_tof_matter = nullptr;
                              TGraph * graph_primary_he3_matter = nullptr;
                              TGraph * graph_primary_p_matter = nullptr;
                              TGraph * graph_secondary_p_from_weak_decay_matter = nullptr;
                              
                              TGraph * graph_purity_he3_antimatter = nullptr;
                              TGraph * graph_purity_p_tpc_antimatter = nullptr;
                              TGraph * graph_purity_p_tof_antimatter = nullptr;
                              TGraph * graph_primary_he3_antimatter = nullptr;
                              TGraph * graph_primary_p_antimatter = nullptr;
                              TGraph * graph_secondary_p_from_weak_decay_antimatter = nullptr;

                              TGraph * graph_sigma_protons_to_weak_decay_protons = nullptr;
                              ''')
    
    ROOT.graph_purity_he3_matter = graphs['Matter']['purity_he3']
    ROOT.graph_purity_p_tpc_matter = graphs['Matter']['purity_p_tpc']
    ROOT.graph_purity_p_tof_matter = graphs['Matter']['purity_p_tof']
    ROOT.graph_primary_he3_matter = graphs['Matter']['primary_he3']
    ROOT.graph_primary_p_matter = graphs['Matter']['primary_p']
    ROOT.graph_secondary_p_from_weak_decay_matter = graphs['Matter']['secondary_p_from_weak_decay']

    ROOT.graph_purity_he3_antimatter = graphs['Matter']['purity_he3']
    ROOT.graph_purity_p_tpc_antimatter = graphs['Matter']['purity_p_tpc']
    ROOT.graph_purity_p_tof_antimatter = graphs['Matter']['purity_p_tof']
    ROOT.graph_primary_he3_antimatter = graphs['Matter']['primary_he3']
    ROOT.graph_primary_p_antimatter = graphs['Matter']['primary_p']
    ROOT.graph_secondary_p_from_weak_decay_antimatter = graphs['Matter']['secondary_p_from_weak_decay']

    ROOT.graph_sigma_protons_to_weak_decay_protons = g_sigma_protons_to_weak_decay_protons

    
    rdf = rdf.Define('fPurityHe3', f'fSignedPtHe3 > 0 ? readFromTGraph(graph_purity_he3_matter, fPtHe3) : \
                     readFromTGraph(graph_purity_he3_antimatter, fPtHe3)') \
             .Define('fPurityHad', f'fSignedPtHad > 0 ? (fPtHad < 0.8 ? readFromTGraph(graph_purity_p_tpc_matter, fPtHad) : readFromTGraph(graph_purity_p_tof_matter, fPtHad)) : \
                     (fPtHad < 0.9 ? readFromTGraph(graph_purity_p_tpc_antimatter, fPtHad) :readFromTGraph(graph_purity_p_tof_antimatter, fPtHad))') \
             .Define('fPrimaryHe3', f'fSignedPtHe3 > 0 ? readFromTGraph(graph_primary_he3_matter, fPtHe3) : \
                     readFromTGraph(graph_primary_he3_antimatter, fPtHe3)') \
             .Define('fPrimaryHad', f'fSignedPtHad > 0 ? readFromTGraph(graph_primary_p_matter, fPtHad) : \
                     readFromTGraph(graph_primary_p_antimatter, fPtHad)') \
             .Define('fSigmaProtonsToWeakDecayProtons', 'readFromTGraph(graph_sigma_protons_to_weak_decay_protons, fPtHad)') \
             .Define('fSecondaryHadFromWeakDecay', f'fSignedPtHad > 0 ? readFromTGraph(graph_secondary_p_from_weak_decay_matter, fPtHad) : \
                     readFromTGraph(graph_secondary_p_from_weak_decay_antimatter, fPtHad)') \
             .Define('fSecondaryHadFromSigma', f'fSigmaProtonsToWeakDecayProtons * fSecondaryHadFromWeakDecay')
    
    for sign in ['Both', 'Matter', 'Antimatter']:
        
        filter_string = 'fSignedPtHe3 > 0' if sign == 'Matter' else 'fSignedPtHe3 < 0' if sign == 'Antimatter' else '1' 
        tmp_rdf = rdf.Filter(filter_string)
        
        outhists = {}
        outhists['kstar'] = tmp_rdf.Histo1D(("hKstar", "^;#it{k}* (GeV/#it{c});", 200, 0., 2.), "fKstar").GetValue()
        outhists['purity_he3'] = tmp_rdf.Histo1D(("hPurityHe3", "^{3}He ;#it{k}* (GeV/#it{c}); Purity (^{3}He)", 200, 0., 2.), "fKstar", "fPurityHe3").GetValue()
        outhists['purity_p'] = tmp_rdf.Histo1D(("hPurityHad", "p ;#it{k}* (GeV/#it{c}); Purity (p)", 200, 0., 2.), "fKstar", "fPurityHad").GetValue()
        outhists['primary_he3'] = tmp_rdf.Histo1D(("hPrimaryHe3", "^{3}He ;#it{k}* (GeV/#it{c}); Primary fraction (^{3}He)", 200, 0., 2.), "fKstar", "fPrimaryHe3").GetValue()
        outhists['primary_p'] = tmp_rdf.Histo1D(("hPrimaryHad", "p ;#it{k}* (GeV/#it{c}); Primary fraction (p)", 200, 0., 2.), "fKstar", "fPrimaryHad").GetValue()
        outhists['secondary_p_from_sigma'] = tmp_rdf.Histo1D(("hSecondaryHadFromSigma", "p ;#it{k}* (GeV/#it{c}); Weak decay fraction (p)", 200, 0., 2.), "fKstar", "fSecondaryHadFromSigma").GetValue()

        outdir = outfile.mkdir(sign)
        outdir.cd()
        for name, hist in outhists.items():
            if name != 'kstar':
                hist.Divide(outhists['kstar'])

            hist.Write()

        outhists['lambda'] = outhists['purity_he3'].Clone('hLambdaParameters')
        outhists['lambda'].Multiply(outhists['purity_p'])
        outhists['lambda'].Multiply(outhists['primary_he3'])
        outhists['lambda'].Multiply(outhists['primary_p'])
        outhists['lambda'].SetTitle('p ;#it{k}* (GeV/#it{c}); #lambda_{p ^{3}He}')

        outhists['lambda_sigma'] = outhists['purity_he3'].Clone('hLambdaSigmaParameters')
        outhists['lambda_sigma'].Multiply(outhists['purity_p'])
        outhists['lambda_sigma'].Multiply(outhists['primary_he3'])
        outhists['lambda_sigma'].Multiply(outhists['secondary_p_from_sigma'])
        outhists['lambda_sigma'].SetTitle('p ;#it{k}* (GeV/#it{c}); #lambda_{#Sigma^{+} ^{3}He}')

        outdir.cd()
        outhists['lambda'].Write()
        outhists['lambda_sigma'].Write()

if __name__ == '__main__':

    proton_from_Sigma_to_proton_from_Lambda_ratio = 0.25

    config_file = '/home/galucia/Lithium4/preparation/config/config_prepare.yml'
    config = yaml.safe_load(open(config_file, 'r'))

    base_selection, selection, purity_primary_selection = prepare_selections(config)
    chain_data = prepare_input_tchain(config)
    rdf = prepare_rdataframe(chain_data, base_selection, selection, purity_primary_selection)

    #output_file_path = 'output/lambda_parameters.root'
    output_file_path = 'output/lambda_parameters_smaller_tolerance_pr.root'
    output_file = ROOT.TFile(output_file_path, "RECREATE")

    purity_and_primary_fraction_histograms(rdf, output_file)    

    output_file.Close()

  