import sys
import argparse
import yaml
import ROOT
from ROOT import TFile, TChain, gInterpreter, RDataFrame

from torchic.core.graph import load_graph
from torchic.core.histogram import load_hist, hist_to_graph
from torchic.utils.terminal_colors import TerminalColors as tc

sys.path.append('..')
from utils.particles import ParticleMasses
from include.load_parameters import load_parametrisation

gInterpreter.ProcessLine(f'#include "../include/Common.h"')
gInterpreter.ProcessLine(f'#include "../include/LambdaUtils.h"')

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
    '((fChi2TPCHad < 4))',
    '(std::abs(fEtaHe3) < 0.9)',
    '(std::abs(fEtaHad) < 0.9)',
    '(fNClsTPCHe3 > 110)',

    '((fPtHe3 < 2.5) || (fPIDtrkHe3 == 7))',
    '((-2. < fNSigmaTPCHe3) && (fNSigmaTPCHe3 < 3.))',
    '(std::abs(fNSigmaDCAxyHe3) < 5)',
    '(std::abs(fNSigmaDCAzHe3) < 5)',
    '(fNSigmaITSHe3 > -2.)',

    '(std::abs(fNSigmaTPCHad) < 2)',
    '((std::abs(fPtHad) < 0.8) || (std::abs(fNSigmaTOFHad) < 2))',  # TOF hard cut
    '(fNSigmaITSHad > -2.)',
    '(std::abs(fNSigmaDCAxyHad) < 5)',
    '(std::abs(fNSigmaDCAzHad) < 5)',
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
      #.Define('fDeltaEta', 'fEtaHe3 - fEtaHad')
      #.Define('fDeltaPhi', 'fPhiHe3 - fPhiHad')
      
    rdf = (rdf.Define('fSignedPtHad', 'fPtHad')
      .Define('fSignHe3', 'fPtHe3/std::abs(fPtHe3)')
      .Redefine('fPtHe3', 'std::abs(fPtHe3)')
      .Redefine('fPtHad', 'std::abs(fPtHad)')
      .Redefine('fPtHe3', '(fPIDtrkHe3 == 7) || (fPIDtrkHe3 == 8) || (fPtHe3 > 2.5) ? fPtHe3 : CorrectPidTrkHe(fPtHe3)')
      .Redefine('fInnerParamTPCHe3', 'fInnerParamTPCHe3 * 2')
      .Define('fSignedPtHe3', 'fPtHe3 * fSignHe3')
      .Define(f'fEHe3', f'std::sqrt((fPtHe3 * std::cosh(fEtaHe3))*(fPtHe3 * std::cosh(fEtaHe3)) + {ParticleMasses["He"]}*{ParticleMasses["He"]})')
      .Define(f'fEHad', f'std::sqrt((fPtHad * std::cosh(fEtaHad))*(fPtHad * std::cosh(fEtaHad)) + {ParticleMasses["Pr"]}*{ParticleMasses["Pr"]})')
      .Define('fClusterSizeCosLamHe3', 'ComputeAverageClusterSize(fItsClusterSizeHe3) / cosh(fEtaHe3)')
      .Define('fClusterSizeCosLamHad', 'ComputeAverageClusterSize(fItsClusterSizeHad) / cosh(fEtaHad)')
      .Define('fExpectedClusterSizeHe3', 'ComputeExpectedClusterSizeCosLambdaHe(fPtHe3 * std::cosh(fEtaHe3))')
      .Define('fExpectedClusterSizeHad', 'ComputeExpectedClusterSizeCosLambdaPr(fPtHad * std::cosh(fEtaHad))')
      .Define('fNSigmaITSHe3', 'ComputeNsigmaITSHe(fPtHe3 * std::cosh(fEtaHe3), fClusterSizeCosLamHe3)')
      .Define('fNSigmaITSHad', 'ComputeNsigmaITSPr(fPtHad * std::cosh(fEtaHad), fClusterSizeCosLamHad)')
      .Redefine('fNSigmaTPCHe3', 'ComputeNsigmaTPCHe(std::abs(fInnerParamTPCHe3), fSignalTPCHe3, false, true)')
      .Define('fNSigmaTPCPi', 'ComputeNsigmaTPCPi(std::abs(fInnerParamTPCHad), fSignalTPCHad)')
      .Define('fNSigmaDCAxyHe3', 'ComputeNsigmaDCAxyHe(fPtHe3, fDCAxyHe3)')
      .Define('fNSigmaDCAzHe3', 'ComputeNsigmaDCAzHe(fPtHe3, fDCAzHe3)')
      .Define('fNSigmaDCAxyHad', 'ComputeNsigmaDCAxyPr(fPtHad, fDCAxyHad)')
      .Define('fNSigmaDCAzHad', 'ComputeNsigmaDCAzPr(fPtHad, fDCAzHad)')
      .Filter(base_selection).Filter(selection)
      .Define('fKstar', f'ComputeKstar(fPtHe3, fEtaHe3, fPhiHe3, {ParticleMasses["He"]}, fPtHad, fEtaHad, fPhiHad, {ParticleMasses["Pr"]})')
    )
    return rdf

def load_same(prepare_years:bool=True):
    
    config_same = {
        '2023': {'input_data':'/data/galucia/lithium/same/LHC23_PbPb_pass5_hadronpid_same.root',},
        '2024': {'input_data':'/data/galucia/lithium/same/LHC24ar_pass3_hadronpid_same.root',},
        '2025': {'input_data':'/data/galucia/lithium/same/LHC25_PbPb_pass1_hadronpid_same.root',},
    }
    
    if prepare_years:    
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

def purity_and_primary_fraction_histograms(rdf:RDataFrame, outfile:TFile, input_file:str):
   
    graphs = {
        'Matter': {},
        'Antimatter': {}
    }

    graphs['Matter']['primary_he3'] = load_graph(input_file, 'primary_fraction/He/g_primary_fraction_matter_averaged')
    graphs['Matter']['primary_p'] = load_graph(input_file, 'primary_fraction/Pr/g_primary_fraction_matter_averaged')
    graphs['Matter']['secondary_p_from_weak_decay'] = load_graph(input_file, 'primary_fraction/Pr/g_weak_decay_fraction_matter_averaged')
    graphs['Matter']['purity_he3'] = load_graph(input_file, 'purity/He3/TPC/g_purity_matter_averaged')
    graphs['Matter']['purity_p_tpc'] = load_graph(input_file, 'purity/Had/TPC/g_purity_matter_averaged')
    graphs['Matter']['purity_p_tof'] = load_graph(input_file, 'purity/Had/TOF/g_purity_matter_averaged')
    graphs['Antimatter']['primary_he3'] = load_graph(input_file, 'primary_fraction/He/g_primary_fraction_antimatter_averaged')
    graphs['Antimatter']['primary_p'] = load_graph(input_file, 'primary_fraction/Pr/g_primary_fraction_antimatter_averaged')
    graphs['Antimatter']['secondary_p_from_weak_decay'] = load_graph(input_file, 'primary_fraction/Pr/g_weak_decay_fraction_antimatter_averaged')
    graphs['Antimatter']['purity_he3'] = load_graph(input_file, 'purity/He3/TPC/g_purity_antimatter_averaged')
    graphs['Antimatter']['purity_p_tpc'] = load_graph(input_file, 'purity/Had/TPC/g_purity_antimatter_averaged')
    graphs['Antimatter']['purity_p_tof'] = load_graph(input_file, 'purity/Had/TOF/g_purity_antimatter_averaged')
    
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

    ROOT.graph_purity_he3_antimatter = graphs['Antimatter']['purity_he3']
    ROOT.graph_purity_p_tpc_antimatter = graphs['Antimatter']['purity_p_tpc']
    ROOT.graph_purity_p_tof_antimatter = graphs['Antimatter']['purity_p_tof']
    ROOT.graph_primary_he3_antimatter = graphs['Antimatter']['primary_he3']
    ROOT.graph_primary_p_antimatter = graphs['Antimatter']['primary_p']
    ROOT.graph_secondary_p_from_weak_decay_antimatter = graphs['Antimatter']['secondary_p_from_weak_decay']

    ROOT.graph_sigma_protons_to_weak_decay_protons = g_sigma_protons_to_weak_decay_protons

    
    rdf = (rdf.Define('fPurityHe3', f'fSignedPtHe3 > 0 ? readFromTGraph(graph_purity_he3_matter, fPtHe3) : \
                      readFromTGraph(graph_purity_he3_antimatter, fPtHe3)')
             .Define('fPurityHad', f'fSignedPtHad > 0 ? (std::abs(fPtHad) < 0.8 ? readFromTGraph(graph_purity_p_tpc_matter, std::abs(fPtHad)) : readFromTGraph(graph_purity_p_tof_matter, std::abs(fPtHad))) : \
                     (std::abs(fPtHad) < 0.8 ? readFromTGraph(graph_purity_p_tpc_antimatter, std::abs(fPtHad)) :readFromTGraph(graph_purity_p_tof_antimatter, std::abs(fPtHad)))')
             .Define('fPrimaryHe3', f'fSignedPtHe3 > 0 ? readFromTGraph(graph_primary_he3_matter, fPtHe3) : \
                     readFromTGraph(graph_primary_he3_antimatter, fPtHe3)')
             .Define('fPrimaryHad', f'fSignedPtHad > 0 ? readFromTGraph(graph_primary_p_matter, fPtHad) : \
                     readFromTGraph(graph_primary_p_antimatter, fPtHad)')
             .Define('fSigmaProtonsToWeakDecayProtons', 'readFromTGraph(graph_sigma_protons_to_weak_decay_protons, fPtHad)')
             .Define('fSecondaryHadFromWeakDecay', f'fSignedPtHad > 0 ? readFromTGraph(graph_secondary_p_from_weak_decay_matter, fPtHad) : \
                     readFromTGraph(graph_secondary_p_from_weak_decay_antimatter, fPtHad)')
             .Define('fSecondaryHadFromSigma', f'fSigmaProtonsToWeakDecayProtons * fSecondaryHadFromWeakDecay')
        )     
    
    all_hists = {}
    for sign in ['Both', 'Matter', 'Antimatter']:
        
        filter_string = 'fSignedPtHe3 > 0' if sign == 'Matter' else 'fSignedPtHe3 < 0' if sign == 'Antimatter' else '1' 
        tmp_rdf = rdf.Filter(filter_string)
        
        all_hists[sign] = {
            'kstar': tmp_rdf.Histo1D(("hKstar", "^; #it{k}* (GeV/#it{c});", 200, 0., 2.), "fKstar"),
            'pt_had': tmp_rdf.Histo1D(("hPtHad", "p ;#it{p}_{T} (GeV/#it{c});", 200, 0., 10.), "fPtHad"),
            'purity_he3': tmp_rdf.Histo1D(("hPurityHe3", "^{3}He ;#it{k}* (GeV/#it{c}); Purity (^{3}He)", 200, 0., 2.), "fKstar", "fPurityHe3"),
            'purity_p': tmp_rdf.Histo1D(("hPurityHad", "p ;#it{k}* (GeV/#it{c}); Purity (p)", 200, 0., 2.), "fKstar", "fPurityHad"),
            'purity_p_pt': tmp_rdf.Histo1D(("hPurityHadPt", "p ;#it{p}_{T} (GeV/#it{c}); Purity (p)", 200, 0., 10.), "fPtHad", "fPurityHad"),   # check
            'primary_he3': tmp_rdf.Histo1D(("hPrimaryHe3", "^{3}He ;#it{k}* (GeV/#it{c}); Primary fraction (^{3}He)", 200, 0., 2.), "fKstar", "fPrimaryHe3"),
            'primary_p': tmp_rdf.Histo1D(("hPrimaryHad", "p ;#it{k}* (GeV/#it{c}); Primary fraction (p)", 200, 0., 2.), "fKstar", "fPrimaryHad"),
            'secondary_p_from_sigma': tmp_rdf.Histo1D(("hSecondaryHadFromSigma", "p ;#it{k}* (GeV/#it{c}); Weak decay fraction (p)", 200, 0., 2.), "fKstar", "fSecondaryHadFromSigma")
        }
    ROOT.RDF.RunGraphs([h for hists in all_hists.values() for h in hists.values()])
    
    for sign in ['Both', 'Matter', 'Antimatter']:
        
        outhists = all_hists[sign]    
        outdir = outfile.mkdir(sign)
        outdir.cd()
        for name, hist in outhists.items():
            if name == 'kstar' or name == 'pt_had':
                pass
            elif 'pt' in name:
                hist.Divide(outhists['pt_had'].GetValue())
            else:
                hist.Divide(outhists['kstar'].GetValue())

            hist.Write()

        outhists['lambda'] = outhists['purity_he3'].Clone('hLambdaParameters')
        outhists['lambda'].Multiply(outhists['purity_p'].GetValue())
        outhists['lambda'].Multiply(outhists['primary_he3'].GetValue())
        outhists['lambda'].Multiply(outhists['primary_p'].GetValue())
        outhists['lambda'].SetTitle('p ;#it{k}* (GeV/#it{c}); #lambda_{p ^{3}He}')

        outhists['lambda_sigma'] = outhists['purity_he3'].Clone('hLambdaSigmaParameters')
        outhists['lambda_sigma'].Multiply(outhists['purity_p'].GetValue())
        outhists['lambda_sigma'].Multiply(outhists['primary_he3'].GetValue())
        outhists['lambda_sigma'].Multiply(outhists['secondary_p_from_sigma'].GetValue())
        outhists['lambda_sigma'].SetTitle('p ;#it{k}* (GeV/#it{c}); #lambda_{#Sigma^{+} ^{3}He}')

        outdir.cd()
        outhists['lambda'].Write()
        outhists['lambda_sigma'].Write()
    
    outfile.cd()
    graphs['Matter']['purity_he3'].Write('graph_purity_he3_matter')
    graphs['Matter']['purity_p_tpc'].Write('graph_purity_p_tpc_matter')
    graphs['Matter']['purity_p_tof'].Write('graph_purity_p_tof_matter')
    graphs['Matter']['primary_he3'].Write('graph_primary_he3_matter')
    graphs['Matter']['primary_p'].Write('graph_primary_p_matter')
    graphs['Matter']['secondary_p_from_weak_decay'].Write('graph_secondary_p_from_weak_decay_matter')
    graphs['Antimatter']['purity_he3'].Write('graph_purity_he3_antimatter')
    graphs['Antimatter']['purity_p_tpc'].Write('graph_purity_p_tpc_antimatter')
    graphs['Antimatter']['purity_p_tof'].Write('graph_purity_p_tof_antimatter')
    graphs['Antimatter']['primary_he3'].Write('graph_primary_he3_antimatter')
    graphs['Antimatter']['primary_p'].Write('graph_primary_p_antimatter')
    graphs['Antimatter']['secondary_p_from_weak_decay'].Write('graph_secondary_p_from_weak_decay_antimatter')


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Compute lambda parameters for He3-p correlation analysis.')
    parser.add_argument('--input', type=str, default='calibration/output/average_purity_and_primary_fraction.root', help='Path to the input configuration YAML file.')
    parser.add_argument('--output', type=str, default='calibration/output/lambda_parameters.root', help='Path to the output ROOT file.')
    args = parser.parse_args()

    proton_from_Sigma_to_proton_from_Lambda_ratio = 0.25
    rdf = load_same(prepare_years=False)
    rdf = rdf.Filter(base_selection)

    output_file = ROOT.TFile(args.output, "RECREATE")

    purity_and_primary_fraction_histograms(rdf, output_file, args.input)

    output_file.Close()

  