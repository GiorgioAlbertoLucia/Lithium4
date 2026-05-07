import sys
import yaml
import numpy as np
import ROOT
from ROOT import TFile, TChain, gInterpreter, RDataFrame, TCanvas, TLegend, gStyle
from alive_progress import alive_bar

import time
import gc

from torchic.utils.terminal_colors import TerminalColors as tc
from torchic.utils.timeit import timeit
from torchic.core.histogram import load_hist
from torchic.utils.root import set_root_object

sys.path.append('..')
from utils.particles import ParticleMasses

gInterpreter.ProcessLine(f'#include "../include/Common.h"')
gInterpreter.ProcessLine(f'#include "../include/Variation.h"')
gInterpreter.ProcessLine(f'#include "../include/Systematics.h"')
from ROOT import ComputeAllSystematics

ROOT.EnableImplicitMT(10)
ROOT.gROOT.SetBatch(True)

selections = [
    '((fSignedPtHe3 > 0 && fSignedPtHad > 0) || (fSignedPtHe3 < 0 && fSignedPtHad < 0))',
    #'((fDeltaPhi/0.006)*(fDeltaPhi/0.006) + (fDeltaEta/0.012)*(fDeltaEta/0.012) > 1)',

    '((fChi2TPCHe3 > 0.5) || (fIs23 == false))', # the cut on chi2tpc should be only applied to 2023
    #'(fChi2TPCHe3 < 4))'
    #'((fChi2TPCHad < 4))',
    '(std::abs(fEtaHe3) < 0.9)',
    '(std::abs(fEtaHad) < 0.9)',
    
    #'(fNClsTPCHe3 > 110)',
    '((fPtHe3 < 2.5) || (fPIDtrkHe3 == 7))',
    '((fPIDtrkHe3 == 6) || (fPIDtrkHe3 == 7) || (fPIDtrkHe3 == 8))',

    #'(std::abs(fNSigmaDCAxyHe3) < 3)',
    #'(std::abs(fNSigmaDCAzHe3) < 3)',
    #'(std::abs(fNSigmaDCAxyHad) < 3)',
    #'(std::abs(fNSigmaDCAzHad) < 3)',
]
#selections = conf['selections']
base_selection = selections[0]
for sel in selections[1:]:
    base_selection += (' && ' + sel)
print(f'Selection: {base_selection}')

def prepare_input_tchain(config:dict):

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
      .Filter(base_selection).Filter(selection) \
      .Filter('(fCentralityFT0C < 50)') \
      .Define('fKstar', f'ComputeKstar(fPtHe3, fEtaHe3, fPhiHe3, {ParticleMasses["He"]}, fPtHad, fEtaHad, fPhiHad, {ParticleMasses["Pr"]})')
    
    return rdf

def load_same():
    
    config_same = {
        'input_data':   [ 
                #'/data/galucia/lithium_local/same_merged/LHC23_PbPb_pass5_same.root',
                #'/data/galucia/lithium_local/same_merged/LHC24_PbPb_pass2_same.root',
                '/data/galucia/lithium_local/same_merged/LHC23_PbPb_pass4_hadronpid_same.root',
                '/data/galucia/lithium_local/same_merged/LHC24ar_pass1_hadronpid_same.root',
                '/data/galucia/lithium_local/same_merged/LHC24as_pass1_hadronpid_same.root',
              ],
        'tree_name':    'O2he3hadtable',
        'mode':         'DF',
    }
    chain_data_same = prepare_input_tchain(config_same)
    rdf_same = prepare_rdataframe(chain_data_same, base_selection, 'true')
    
    return rdf_same, chain_data_same

def load_mixed():
    
    config_mixed = {
        'input_data':   [ 
                '/data/galucia/lithium_local/mixing/LHC23_PbPb_pass4_hadronpid_event_mixing_batch1995_refined_dca.root',
                '/data/galucia/lithium_local/mixing/LHC24ar_pass1_hadronpid_event_mixing_batch1995_refined_dca.root',
                '/data/galucia/lithium_local/mixing/LHC24as_pass1_hadronpid_event_mixing_batch1995_refined_dca.root',
                #'/data/galucia/lithium_local/mixing/LHC23_PbPb_pass4_hadronpid_event_mixing_batch105_refined_dca.root',
                #'/data/galucia/lithium_local/mixing/LHC24ar_pass1_hadronpid_event_mixing_batch105_refined_dca.root',
                #'/data/galucia/lithium_local/mixing/LHC24as_pass1_hadronpid_event_mixing_batch105_refined_dca.root',
                '/data/galucia/lithium_local/mixing/LHC23_PbPb_pass4_hadronpid_event_mixing_batch42_refined_dca.root',
                '/data/galucia/lithium_local/mixing/LHC24ar_pass1_hadronpid_event_mixing_batch42_refined_dca.root',
                '/data/galucia/lithium_local/mixing/LHC24as_pass1_hadronpid_event_mixing_batch42_refined_dca.root',
                '/data/galucia/lithium_local/mixing/LHC23_PbPb_pass4_hadronpid_event_mixing_batch256_refined_dca.root',
                '/data/galucia/lithium_local/mixing/LHC24ar_pass1_hadronpid_event_mixing_batch256_refined_dca.root',
                '/data/galucia/lithium_local/mixing/LHC24as_pass1_hadronpid_event_mixing_batch256_refined_dca.root',
                '/data/galucia/lithium_local/mixing/LHC23_PbPb_pass4_hadronpid_event_mixing_batch3112_refined_dca.root',
                '/data/galucia/lithium_local/mixing/LHC24ar_pass1_hadronpid_event_mixing_batch3112_refined_dca.root',
                '/data/galucia/lithium_local/mixing/LHC24as_pass1_hadronpid_event_mixing_batch3112_refined_dca.root',
              ],
        'tree_name':    'MixedTree',
        'mode':         'tree',
    }
    chain_data_mixed = prepare_input_tchain(config_mixed)
    rdf_mixed = prepare_rdataframe(chain_data_mixed, base_selection, 'true')
    
    return rdf_mixed, chain_data_mixed

@timeit
def prepare_systematics_histograms(rdf:RDataFrame, n_variations:int, hist_name_suffix:str='', variable:str=''):

    variational_rdf = None
    if variable == 'all':
        variational_rdf = rdf.Define('fNSigmaITSHe3Systematic', 'fNSigmaITSHe3') \
                         .Define('fAbsNSigmaTPCHe3Systematic', 'std::abs(fNSigmaTPCHe3)') \
                         .Define('fAbsNSigmaDCAxyHe3Systematic', 'std::abs(fNSigmaDCAxyHe3)') \
                         .Define('fAbsNSigmaDCAzHe3Systematic', 'std::abs(fNSigmaDCAzHe3)') \
                         .Define('fNClsTPCHe3Systematic', 'static_cast<float>(fNClsTPCHe3)') \
                         .Define('fChi2TPCHe3Systematic', 'fChi2TPCHe3') \
                         .Define('fNSigmaITSHadSystematic', 'fNSigmaITSHad') \
                         .Define('fAbsNSigmaTPCHadSystematic', 'std::abs(fNSigmaTPCHad)') \
                         .Define('fAbsNSigmaTOFHadSystematic', 'std::abs(fNSigmaTOFHad)') \
                         .Define('fAbsNSigmaDCAxyHadSystematic', 'std::abs(fNSigmaDCAxyHad)') \
                         .Define('fAbsNSigmaDCAzHadSystematic', 'std::abs(fNSigmaDCAzHad)') \
                         .Define('fChi2TPCHadSystematic', 'fChi2TPCHad') \
                         .Define('fAbsZVertexSystematic', 'std::abs(fZVertex)') \
                         .Vary(['fNSigmaITSHe3Systematic', 'fAbsNSigmaTPCHe3Systematic',
                                'fAbsNSigmaDCAxyHe3Systematic', 'fAbsNSigmaDCAzHe3Systematic',
                                'fNClsTPCHe3Systematic', 'fChi2TPCHe3Systematic',
                                'fNSigmaITSHadSystematic', 'fAbsNSigmaTPCHadSystematic', 'fAbsNSigmaTOFHadSystematic',
                                'fAbsNSigmaDCAxyHadSystematic', 'fAbsNSigmaDCAzHadSystematic', 'fAbsZVertexSystematic'],
                                f'systematicCuts({n_variations}, fNSigmaITSHe3Systematic, fAbsNSigmaTPCHe3Systematic, \
                                fAbsNSigmaDCAxyHe3Systematic, fAbsNSigmaDCAzHe3Systematic, fNClsTPCHe3Systematic, fChi2TPCHe3Systematic, \
                                fNSigmaITSHadSystematic, fAbsNSigmaTPCHadSystematic, fAbsNSigmaTOFHadSystematic, \
                                fAbsNSigmaDCAxyHadSystematic, fAbsNSigmaDCAzHadSystematic, fChi2TPCHadSystematic, \
                                fAbsZVertexSystematic)',
                                n_variations, 'systematic') \
                         .Filter('(fNSigmaITSHe3Systematic > 0) && \
                                  (fAbsNSigmaTPCHe3Systematic < 0) && \
                                  (fAbsNSigmaDCAxyHe3Systematic < 0) && \
                                  (fAbsNSigmaDCAzHe3Systematic < 0) && \
                                  (fNClsTPCHe3Systematic > 0) && \
                                  (fChi2TPCHe3Systematic < 0) && \
                                  (fNSigmaITSHadSystematic > 0) && \
                                  (fAbsNSigmaTPCHadSystematic < 0) && \
                                  ((std::abs(fPtHad) < 0.8) || (fAbsNSigmaTOFHadSystematic < 0)) && \
                                  (fAbsNSigmaDCAxyHadSystematic < 0) && \
                                  (fAbsNSigmaDCAzHadSystematic < 0) && \
                                  (fAbsZVertexSystematic < 0)')
    else:
        variational_rdf = rdf.Define('fNSigmaITSHe3Systematic', 'fNSigmaITSHe3') \
                            .Define('fAbsNSigmaTPCHe3Systematic', 'std::abs(fNSigmaTPCHe3)') \
                            .Define('fAbsNSigmaDCAxyHe3Systematic', 'std::abs(fNSigmaDCAxyHe3)') \
                            .Define('fAbsNSigmaDCAzHe3Systematic', 'std::abs(fNSigmaDCAzHe3)') \
                            .Define('fNClsTPCHe3Systematic', 'static_cast<float>(fNClsTPCHe3)') \
                            .Define('fChi2TPCHe3Systematic', 'fChi2TPCHe3') \
                            .Define('fNSigmaITSHadSystematic', 'fNSigmaITSHad') \
                            .Define('fAbsNSigmaTPCHadSystematic', 'std::abs(fNSigmaTPCHad)') \
                            .Define('fAbsNSigmaTOFHadSystematic', 'std::abs(fNSigmaTOFHad)') \
                            .Define('fAbsNSigmaDCAxyHadSystematic', 'std::abs(fNSigmaDCAxyHad)') \
                            .Define('fAbsNSigmaDCAzHadSystematic', 'std::abs(fNSigmaDCAzHad)') \
                            .Define('fChi2TPCHadSystematic', 'fChi2TPCHad') \
                            .Define('fAbsZVertexSystematic', 'std::abs(fZVertex)') \
                            .Vary(['fNSigmaITSHe3Systematic', 'fAbsNSigmaTPCHe3Systematic',
                                    'fAbsNSigmaDCAxyHe3Systematic', 'fAbsNSigmaDCAzHe3Systematic',
                                    'fNClsTPCHe3Systematic', 'fChi2TPCHe3Systematic',
                                    'fNSigmaITSHadSystematic', 'fAbsNSigmaTPCHadSystematic', 'fAbsNSigmaTOFHadSystematic',
                                    'fAbsNSigmaDCAxyHadSystematic', 'fAbsNSigmaDCAzHadSystematic', 'fAbsZVertexSystematic'],
                                    f'systematicCutsSingleVariable({n_variations}, "{variable}", fNSigmaITSHe3Systematic, fAbsNSigmaTPCHe3Systematic, \
                                    fAbsNSigmaDCAxyHe3Systematic, fAbsNSigmaDCAzHe3Systematic, fNClsTPCHe3Systematic, fChi2TPCHe3Systematic, \
                                    fNSigmaITSHadSystematic, fAbsNSigmaTPCHadSystematic, fAbsNSigmaTOFHadSystematic, \
                                    fAbsNSigmaDCAxyHadSystematic, fAbsNSigmaDCAzHadSystematic, fChi2TPCHadSystematic, \
                                    fAbsZVertexSystematic)',
                                    n_variations, 'systematic') \
                            .Filter('(fNSigmaITSHe3Systematic > 0) && \
                                  (fAbsNSigmaTPCHe3Systematic < 0) && \
                                  (fAbsNSigmaDCAxyHe3Systematic < 0) && \
                                  (fAbsNSigmaDCAzHe3Systematic < 0) && \
                                  (fNClsTPCHe3Systematic > 0) && \
                                  (fChi2TPCHe3Systematic < 0) && \
                                  (fNSigmaITSHadSystematic > 0) && \
                                  (fAbsNSigmaTPCHadSystematic < 0) && \
                                  ((std::abs(fPtHad) < 0.8) || (fAbsNSigmaTOFHadSystematic < 0)) && \
                                  (fAbsNSigmaDCAxyHadSystematic < 0) && \
                                  (fAbsNSigmaDCAzHadSystematic < 0) && \
                                  (fAbsZVertexSystematic < 0)')

    nominal_hist_010 = variational_rdf.Filter('fCentralityFT0C < 10') \
                                      .Histo1D((f'hKstar010{hist_name_suffix}', ';#it{k}* (GeV/#it{c});', 20, 0, 0.4), 'fKstar')
    variational_hists_010 = ROOT.RDF.Experimental.VariationsFor(nominal_hist_010)

    nominal_hist_1030 = variational_rdf.Filter('fCentralityFT0C >= 10 && fCentralityFT0C < 30') \
                                      .Histo1D((f'hKstar1030{hist_name_suffix}', ';#it{k}* (GeV/#it{c});', 20, 0, 0.4), 'fKstar')
    variational_hists_1030 = ROOT.RDF.Experimental.VariationsFor(nominal_hist_1030)
    
    nominal_hist_3050 = variational_rdf.Filter('fCentralityFT0C >= 30 && fCentralityFT0C < 50') \
                                      .Histo1D((f'hKstar3050{hist_name_suffix}', ';#it{k}* (GeV/#it{c});', 20, 0, 0.4), 'fKstar')
    variational_hists_3050 = ROOT.RDF.Experimental.VariationsFor(nominal_hist_3050)

    return variational_hists_010, variational_hists_1030, variational_hists_3050

def normalise_histogram(h_same, h_mixed, NORM_LOW_KSTAR, NORM_HIGH_KSTAR):

    low_bin = h_same.FindBin(NORM_LOW_KSTAR)
    high_bin = h_same.FindBin(NORM_HIGH_KSTAR)
    normalization_factor = h_same.Integral(low_bin, high_bin) / h_mixed.Integral(low_bin, high_bin)
    h_mixed.Scale(normalization_factor)

    return h_mixed

def correlation_function_centrality_integrated(h_sames, h_mixeds, suffix:str):
    """
    Compute the correlation function for the centrality integrated histograms.
    """

    h_same = h_sames[0].Clone()
    h_mixed = h_mixeds[0].Clone()

    for h_same_centrality, h_mixed_centrality in zip(h_sames[1:], h_mixeds[1:]):
        h_same.Add(h_same_centrality)
        h_mixed.Add(h_mixed_centrality)

    h_corr = h_same.Clone(f'hCorrelation{suffix}')
    h_corr.Divide(h_mixed)

    return h_corr

def run_systematics():

    rdf_same, _chain_same = load_same()
    rdf_mixed, _chain_mixed = load_mixed()

    N_ITERATIONS = 200
    outFile = TFile("output/systematics_cpp.root", "RECREATE")

    for sign, condition in {'Matter': 'fSignedPtHe3 > 0', 'Antimatter': 'fSignedPtHe3 < 0'}.items():

        outDir = outFile.mkdir(sign)

        tmp_rdf_same = rdf_same.Filter(condition)
        tmp_rdf_mixed = rdf_mixed.Filter(condition)

        variational_hists_same = prepare_systematics_histograms(tmp_rdf_same, n_variations=N_ITERATIONS, hist_name_suffix='Same')
        variational_hists_mixed = prepare_systematics_histograms(tmp_rdf_mixed, n_variations=N_ITERATIONS, hist_name_suffix='Mixed')

        h_correlations = []
        NORM_LOW_KSTAR, NORM_HIGH_KSTAR = 0.2, 0.4
        #
        #vec_same = ROOT.std.vector[ROOT.RDF.Experimental.RResultMap["TH1D"]]()
        #for elem in variational_hists_same:
        #    vec_same.push_back(elem)
        #   
        #vec_mixed = ROOT.std.vector[ROOT.RDF.Experimental.RResultMap["TH1D"]]()
        #for elem in variational_hists_mixed:
        #    vec_mixed.push_back(elem)
        #
        #h_correlations = ComputeAllSystematics(vec_same, vec_mixed, N_ITERATIONS, NORM_LOW_KSTAR, NORM_HIGH_KSTAR)

        with alive_bar(N_ITERATIONS, title='Running systematics...') as bar:
            for iter in range(N_ITERATIONS):

                same_010, same_1030, same_3050 = variational_hists_same[0][f'systematic:{iter}'], \
                                                variational_hists_same[1][f'systematic:{iter}'], variational_hists_same[2][f'systematic:{iter}']
                mixed_010, mixed_1030, mixed_3050 = variational_hists_mixed[0][f'systematic:{iter}'], \
                                                variational_hists_mixed[1][f'systematic:{iter}'], variational_hists_mixed[2][f'systematic:{iter}']
                
                mixed_normalised_010 = normalise_histogram(same_010, mixed_010, NORM_LOW_KSTAR, NORM_HIGH_KSTAR)
                mixed_normalised_1030 = normalise_histogram(same_1030, mixed_1030, NORM_LOW_KSTAR, NORM_HIGH_KSTAR)
                mixed_normalised_3050 = normalise_histogram(same_3050, mixed_3050, NORM_LOW_KSTAR, NORM_HIGH_KSTAR)

                h_correlation_iter = correlation_function_centrality_integrated([same_010, same_1030, same_3050],
                                                                                [mixed_normalised_010, mixed_normalised_1030, mixed_normalised_3050],
                                                                                iter)
                
                h_correlations.append(h_correlation_iter)
                bar()

        point_positions = {}
        subDir = outDir.mkdir('correlations')
        subDir.cd()
        with alive_bar(h_correlations[0].GetNbinsX(), title='Evaluating systematics from histograms...') as bar:
            for ibin in range(1, h_correlations[0].GetNbinsX()+1):
                for hist in h_correlations:
                    if ibin == 1:
                        hist.Write()
                    point_positions[ibin] = point_positions.get(ibin, []) + [hist.GetBinContent(ibin)]
                bar()

        hist_systematics = h_correlations[0].Clone('hSystematics')
        hist_systematics.Reset()
        for ibin in range(1, hist_systematics.GetNbinsX()+1):
            point_systematics = np.std(point_positions[ibin])
            hist_systematics.SetBinContent(ibin, point_systematics)

        outDir.cd()
        hist_systematics.Write()
    
    outFile.Close()

def run_barlow_test():

    rdf_same, _chain_same = load_same()
    rdf_mixed, _chain_mixed = load_mixed()

    N_ITERATIONS = 20
    outFile = TFile("output/systematics_barlow.root", "RECREATE")

    for sign, condition in {'Matter': 'fSignedPtHe3 > 0', 'Antimatter': 'fSignedPtHe3 < 0'}.items():

        outDir = outFile.mkdir(sign)

        tmp_rdf_same = rdf_same.Filter(condition)
        tmp_rdf_mixed = rdf_mixed.Filter(condition)

        variational_hists_same = prepare_systematics_histograms(tmp_rdf_same, n_variations=1, hist_name_suffix='Same', variable='')
        variational_hists_mixed = prepare_systematics_histograms(tmp_rdf_mixed, n_variations=1, hist_name_suffix='Mixed', variable='')

        NORM_LOW_KSTAR, NORM_HIGH_KSTAR = 0.2, 0.4

        same_010, same_1030, same_3050 = variational_hists_same[0][f'systematic:{0}'], \
                                        variational_hists_same[1][f'systematic:{0}'], variational_hists_same[2][f'systematic:{0}']
        mixed_010, mixed_1030, mixed_3050 = variational_hists_mixed[0][f'systematic:{0}'], \
                                        variational_hists_mixed[1][f'systematic:{0}'], variational_hists_mixed[2][f'systematic:{0}']

        mixed_clone_010 = mixed_010.Clone()
        mixed_clone_1030 = mixed_1030.Clone()
        mixed_clone_3050 = mixed_3050.Clone()

        mixed_normalised_010 = normalise_histogram(same_010, mixed_clone_010, NORM_LOW_KSTAR, NORM_HIGH_KSTAR)
        mixed_normalised_1030 = normalise_histogram(same_1030, mixed_clone_1030, NORM_LOW_KSTAR, NORM_HIGH_KSTAR)
        mixed_normalised_3050 = normalise_histogram(same_3050, mixed_clone_3050, NORM_LOW_KSTAR, NORM_HIGH_KSTAR)

        h_correlation_ref = correlation_function_centrality_integrated([same_010, same_1030, same_3050],
                                                                        [mixed_normalised_010, mixed_normalised_1030, mixed_normalised_3050],
                                                                        0)
        outDir.cd()
        h_correlation_ref.Write('hCorrelationReference')

        for hist in [same_010, same_1030, same_3050, mixed_normalised_010, mixed_normalised_1030, mixed_normalised_3050]:
            del hist
        

        for variable in ['fNSigmaITSHe3', 'fNSigmaTPCHe3', 'fNSigmaITSHad', 'fNSigmaTPCHad', 'fNSigmaTOFHad']:

            gc.collect()

            print(tc.YELLOW+f'\nRunning Barlow test for systematic variable: {variable}'+tc.RESET)
            outDirBarlow = outDir.mkdir(variable)

            variational_hists_same = prepare_systematics_histograms(tmp_rdf_same, n_variations=N_ITERATIONS, hist_name_suffix='Same', variable=variable)
            variational_hists_mixed = prepare_systematics_histograms(tmp_rdf_mixed, n_variations=N_ITERATIONS, hist_name_suffix='Mixed', variable=variable)

            h_correlations = []
            NORM_LOW_KSTAR, NORM_HIGH_KSTAR = 0.2, 0.4

            with alive_bar(N_ITERATIONS, title=f'Running Barlow for {variable}...') as bar:
                for iter in range(N_ITERATIONS):

                    same_010, same_1030, same_3050 = variational_hists_same[0][f'systematic:{iter}'], \
                                                    variational_hists_same[1][f'systematic:{iter}'], variational_hists_same[2][f'systematic:{iter}']
                    mixed_010, mixed_1030, mixed_3050 = variational_hists_mixed[0][f'systematic:{iter}'], \
                                                    variational_hists_mixed[1][f'systematic:{iter}'], variational_hists_mixed[2][f'systematic:{iter}']
                    
                    mixed_clone_010 = mixed_010.Clone()
                    mixed_clone_1030 = mixed_1030.Clone()
                    mixed_clone_3050 = mixed_3050.Clone()

                    mixed_normalised_010 = normalise_histogram(same_010, mixed_clone_010, NORM_LOW_KSTAR, NORM_HIGH_KSTAR)
                    mixed_normalised_1030 = normalise_histogram(same_1030, mixed_clone_1030, NORM_LOW_KSTAR, NORM_HIGH_KSTAR)
                    mixed_normalised_3050 = normalise_histogram(same_3050, mixed_clone_3050, NORM_LOW_KSTAR, NORM_HIGH_KSTAR)

                    h_correlation_iter = correlation_function_centrality_integrated([same_010, same_1030, same_3050],
                                                                                    [mixed_normalised_010, mixed_normalised_1030, mixed_normalised_3050],
                                                                                    iter)

                    h_correlations.append(h_correlation_iter)

                    for hist in [same_010, same_1030, same_3050, mixed_normalised_010, mixed_normalised_1030, mixed_normalised_3050]:
                        del hist

                    bar()

            point_positions = {}
            subDir = outDirBarlow.mkdir('correlations')
            subDir.cd()
            with alive_bar(h_correlations[0].GetNbinsX(), title='Evaluating systematics from histograms...') as bar:
                for ibin in range(1, h_correlations[0].GetNbinsX()+1):
                    for hist in h_correlations:
                        if ibin == 1:
                            hist.Write()
                        point_positions[ibin] = point_positions.get(ibin, []) + [hist.GetBinContent(ibin)]
                    bar()

            hist_mean = h_correlations[0].Clone('hMean')
            hist_std = h_correlations[0].Clone('hStd')
            hist_barlow = h_correlations[0].Clone('hBarlow')
            hist_mean.Reset()
            hist_std.Reset()
            for ibin in range(1, hist_mean.GetNbinsX()+1):
                point_mean = np.mean(point_positions[ibin])
                hist_mean.SetBinContent(ibin, point_mean)
                point_std = np.std(point_positions[ibin])
                hist_std.SetBinContent(ibin, point_std)
                point_barlow = (point_mean - h_correlation_ref.GetBinContent(ibin)) / np.sqrt(np.abs(h_correlation_ref.GetBinError(ibin)*h_correlation_ref.GetBinError(ibin) - point_std*point_std)) if \
                                np.sqrt(np.abs(h_correlation_ref.GetBinError(ibin)*h_correlation_ref.GetBinError(ibin) - point_std*point_std)) > 0 else 999
                hist_barlow.SetBinContent(ibin, point_barlow)

            outDirBarlow.cd()
            hist_mean.Write()
            hist_std.Write()
            hist_barlow.Write()
    
    outFile.Close()

def run_barlow_test_all_variables():

    rdf_same, _chain_same = load_same()
    rdf_mixed, _chain_mixed = load_mixed()

    N_ITERATIONS = 100
    outFile = TFile("output/systematics_barlow_all_variables_extended.root", "RECREATE")

    for sign, condition in {'Matter': 'fSignedPtHe3 > 0', 'Antimatter': 'fSignedPtHe3 < 0'}.items():

        outDir = outFile.mkdir(sign)

        tmp_rdf_same = rdf_same.Filter(condition)
        tmp_rdf_mixed = rdf_mixed.Filter(condition)

        variational_hists_same = prepare_systematics_histograms(tmp_rdf_same, n_variations=1, hist_name_suffix='Same', variable='')
        variational_hists_mixed = prepare_systematics_histograms(tmp_rdf_mixed, n_variations=1, hist_name_suffix='Mixed', variable='')

        NORM_LOW_KSTAR, NORM_HIGH_KSTAR = 0.2, 0.4

        same_010, same_1030, same_3050 = variational_hists_same[0][f'systematic:{0}'], \
                                        variational_hists_same[1][f'systematic:{0}'], variational_hists_same[2][f'systematic:{0}']
        mixed_010, mixed_1030, mixed_3050 = variational_hists_mixed[0][f'systematic:{0}'], \
                                        variational_hists_mixed[1][f'systematic:{0}'], variational_hists_mixed[2][f'systematic:{0}']

        mixed_normalised_010 = normalise_histogram(same_010, mixed_010, NORM_LOW_KSTAR, NORM_HIGH_KSTAR)
        mixed_normalised_1030 = normalise_histogram(same_1030, mixed_1030, NORM_LOW_KSTAR, NORM_HIGH_KSTAR)
        mixed_normalised_3050 = normalise_histogram(same_3050, mixed_3050, NORM_LOW_KSTAR, NORM_HIGH_KSTAR)

        h_correlation_ref = correlation_function_centrality_integrated([same_010, same_1030, same_3050],
                                                                        [mixed_normalised_010, mixed_normalised_1030, mixed_normalised_3050],
                                                                        0)
        outDir.cd()
        h_correlation_ref.Write('hCorrelationReference')

        gc.collect()

        variational_hists_same = prepare_systematics_histograms(tmp_rdf_same, n_variations=N_ITERATIONS, hist_name_suffix='Same', variable='all')
        variational_hists_mixed = prepare_systematics_histograms(tmp_rdf_mixed, n_variations=N_ITERATIONS, hist_name_suffix='Mixed', variable='all')

        h_correlations = []
        NORM_LOW_KSTAR, NORM_HIGH_KSTAR = 0.2, 0.4

        with alive_bar(N_ITERATIONS, title=f'Running Barlow for all PID variables...') as bar:
            for iter in range(N_ITERATIONS):

                same_010, same_1030, same_3050 = variational_hists_same[0][f'systematic:{iter}'], \
                                                variational_hists_same[1][f'systematic:{iter}'], variational_hists_same[2][f'systematic:{iter}']
                mixed_010, mixed_1030, mixed_3050 = variational_hists_mixed[0][f'systematic:{iter}'], \
                                                variational_hists_mixed[1][f'systematic:{iter}'], variational_hists_mixed[2][f'systematic:{iter}']

                mixed_normalised_010 = normalise_histogram(same_010, mixed_010, NORM_LOW_KSTAR, NORM_HIGH_KSTAR)
                mixed_normalised_1030 = normalise_histogram(same_1030, mixed_1030, NORM_LOW_KSTAR, NORM_HIGH_KSTAR)
                mixed_normalised_3050 = normalise_histogram(same_3050, mixed_3050, NORM_LOW_KSTAR, NORM_HIGH_KSTAR)

                h_correlation_iter = correlation_function_centrality_integrated([same_010, same_1030, same_3050],
                                                                                [mixed_normalised_010, mixed_normalised_1030, mixed_normalised_3050],
                                                                                iter)

                h_correlations.append(h_correlation_iter)
                bar()

        point_positions = {}
        subDir = outDir.mkdir('correlations')
        subDir.cd()
        with alive_bar(h_correlations[0].GetNbinsX(), title='Evaluating systematics from histograms...') as bar:
            for ibin in range(1, h_correlations[0].GetNbinsX()+1):
                for hist in h_correlations:
                    if ibin == 1:
                        hist.Write()
                    point_positions[ibin] = point_positions.get(ibin, []) + [hist.GetBinContent(ibin)]
                bar()

        hist_mean = h_correlations[0].Clone('hMean')
        hist_std = h_correlations[0].Clone('hStd')
        hist_barlow = h_correlations[0].Clone('hBarlow')
        hist_mean.Reset()
        hist_std.Reset()
        for ibin in range(1, hist_mean.GetNbinsX()+1):
            point_mean = np.mean(point_positions[ibin])
            hist_mean.SetBinContent(ibin, point_mean)
            point_std = np.std(point_positions[ibin])
            hist_std.SetBinContent(ibin, point_std)
            point_barlow = (point_mean - h_correlation_ref.GetBinContent(ibin)) / np.sqrt(np.abs(h_correlation_ref.GetBinError(ibin)*h_correlation_ref.GetBinError(ibin) - point_std*point_std)) if \
                            np.sqrt(np.abs(h_correlation_ref.GetBinError(ibin)*h_correlation_ref.GetBinError(ibin) - point_std*point_std)) > 0 else 999
            hist_barlow.SetBinContent(ibin, point_barlow)

        outDir.cd()
        hist_mean.Write()
        hist_std.Write()
        hist_barlow.Write()
    
    outFile.Close()

def display_relative_systematics():

    h_systematics = load_hist('/home/galucia/Lithium4/preparation/output/systematics_barlow_all_variables.root',
                              'Antimatter/hStd')
    h_correlation_ref = load_hist('/home/galucia/Lithium4/preparation/output/systematics_barlow_all_variables.root',
                              'Antimatter/hCorrelationReference')
    
    h_relative_systematics = h_systematics.Clone('hRelativeSystematics')
    h_relative_statistics = h_correlation_ref.Clone('hRelativeStatistics')
    for ibin in range(1, h_relative_systematics.GetNbinsX()+1):
        relative_systematic = h_systematics.GetBinContent(ibin) / h_correlation_ref.GetBinContent(ibin) if h_correlation_ref.GetBinContent(ibin) != 0 else 0
        h_relative_systematics.SetBinContent(ibin, relative_systematic)
        relative_statistic = h_correlation_ref.GetBinError(ibin) / h_correlation_ref.GetBinContent(ibin) if h_correlation_ref.GetBinContent(ibin) != 0 else 0
        h_relative_statistics.SetBinContent(ibin, relative_statistic)

    gStyle.SetOptStat(0)
    canvas = TCanvas('cRelativeSystematics', 'cRelativeSystematics', 800, 600)
    set_root_object(h_relative_systematics, line_color=797, line_width=2, title=';#it{k}* (GeV/#it{c});Relative Uncertainty')
    set_root_object(h_relative_statistics, line_color=632, line_width=2, title=';#it{k}* (GeV/#it{c});Relative Uncertainty')
    h_relative_statistics.Draw('hist')
    h_relative_systematics.Draw('hist same')

    legend = TLegend(0.55, 0.65, 0.85, 0.75)
    legend.SetBorderSize(0)
    legend.AddEntry(h_relative_statistics, '#sigma_{stat} / C(#it{k}*)', 'l')
    legend.AddEntry(h_relative_systematics, '#sigma_{syst} / C(#it{k}*)', 'l')
    legend.Draw()

    canvas.SaveAs('/home/galucia/Lithium4/figures/systematics/hRelativeSystematics.pdf')

if __name__ == "__main__":

    #run_systematics()
    
    ## These work
    #run_barlow_test()
    run_barlow_test_all_variables()

    # display
    #display_relative_systematics()



    #run_indivisual_systematics()
    #display_systematics()
    #display_individual_systematics()
    print("Systematics analysis completed successfully.")