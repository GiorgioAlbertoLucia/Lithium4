import sys
import yaml
import argparse
from collections import Counter

import ROOT
from ROOT import TFile, TChain, gInterpreter, RDataFrame

from torchic.utils.terminal_colors import TerminalColors as tc

sys.path.append('..')
from utils.particles import ParticleMasses

from include.load_parameters import load_parametrisation
gInterpreter.ProcessLine(f'#include "../include/Common.h"')

ROOT.EnableImplicitMT(30)
ROOT.gROOT.SetBatch(True)

# kstar binning used for the "multiplicity vs kstar" 2D histograms.
# Adjust to match whatever binning you use in histogram_archive.py.
KSTAR_NBINS = 20
KSTAR_MIN = 0.0
KSTAR_MAX = 0.4  # GeV/c

# multiplicity axis binning
MULT_NBINS = 1000
MULT_MIN = 0.5
MULT_MAX = 1000.5


def prepare_selections(selections):

    base_selection = ''
    if config.get('like_sign', True):
        base_selection = '(fSignedPtHe3 > 0 && fSignedPtHad > 0) || (fSignedPtHe3 < 0 && fSignedPtHad < 0)'
    else:
        base_selection = '(fSignedPtHe3 > 0 && fSignedPtHad < 0) || (fSignedPtHe3 < 0 && fSignedPtHad > 0)'

    selection = selections[0]
    for sel in selections[1:]:
        selection += (' && ' + sel)

    return base_selection, selection


def prepare_input_tchain(config: dict):

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
                if 'DF_' in key_name:
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


def prepare_rdataframe_with_multiplicity(chain_data: TChain, base_selection: str, selection: str):

    rdf = RDataFrame(chain_data)
    # NOTE: unlike the original script, we do NOT filter out duplicate He3 entries here -
    # we need every occurrence to know how many times each He3 appears, and to correlate
    # that count with the kstar of each individual pair.
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
        rdf = rdf.Define('fNSigmaTOFHad', 'ComputeNsigmaTOFPr(std::abs(fPtHad), fMassTOFHad)')

    rdf = (rdf.Define('fSignedPtHad', 'fPtHad')
      .Define('fSignHe3', 'fPtHe3/std::abs(fPtHe3)')
      .Redefine('fPtHe3', 'std::abs(fPtHe3)')
      .Redefine('fPtHad', 'std::abs(fPtHad)')
      .Redefine('fPtHe3', '(fPIDtrkHe3 == 7) || (fPIDtrkHe3 == 8) || (fPtHe3 > 2.5) ? fPtHe3 : CorrectPidTrkHe(fPtHe3)')
      .Redefine('fInnerParamTPCHe3', 'fInnerParamTPCHe3 * 2')
      .Redefine('fInnerParamTPCHe3', '(fPIDtrkHe3 == 7) || (fPIDtrkHe3 == 8) || (fInnerParamTPCHe3 > 2.4) ? fInnerParamTPCHe3 : CorrectPidTrkHe(fInnerParamTPCHe3, false)')
      .Define('fSignedPtHe3', 'fPtHe3 * fSignHe3')
      .Define(f'fEHe3', f'std::sqrt((fPtHe3 * std::cosh(fEtaHe3))*(fPtHe3 * std::cosh(fEtaHe3)) + {ParticleMasses["He"]}*{ParticleMasses["He"]})')
      .Define(f'fEHad', f'std::sqrt((fPtHad * std::cosh(fEtaHad))*(fPtHad * std::cosh(fEtaHad)) + {ParticleMasses["Pr"]}*{ParticleMasses["Pr"]})')
      .Define('fDeltaEta', 'fEtaHe3 - fEtaHad')
      .Define('fDeltaPhi', 'fPhiHe3 - fPhiHad')
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
      # multiplicity of this He3 (how many rows share its (fZVertex, fPtHe3) key), and whether
      # this row is the first occurrence of that He3 in the chain
      .Define('fHe3Multiplicity', 'get_he3_multiplicity(rdfentry_)')
      .Define('fHe3IsFirstOccurrence', 'is_first_occurrence(rdfentry_)')
      )

    return rdf


def fill_multiplicity_histograms(rdf, output_file: TFile):

    print(f'\n{tc.GREEN}Filling He3 multiplicity histograms{tc.RESET}')

    h_mult_matter = ROOT.TH1F('hHe3Multiplicity_Matter', 'He3 multiplicity (matter);N occurrences;counts',
                               MULT_NBINS, MULT_MIN, MULT_MAX)
    h_mult_antimatter = ROOT.TH1F('hHe3Multiplicity_Antimatter', 'He3 multiplicity (antimatter);N occurrences;counts',
                                   MULT_NBINS, MULT_MIN, MULT_MAX)

    h_mult_vs_kstar_matter = ROOT.TH2F('hHe3Multiplicity_vs_Kstar_Matter',
                                        'He3 multiplicity vs k* (matter);k* (GeV/c);N occurrences',
                                        KSTAR_NBINS, KSTAR_MIN, KSTAR_MAX, MULT_NBINS, MULT_MIN, MULT_MAX)
    h_mult_vs_kstar_antimatter = ROOT.TH2F('hHe3Multiplicity_vs_Kstar_Antimatter',
                                            'He3 multiplicity vs k* (antimatter);k* (GeV/c);N occurrences',
                                            KSTAR_NBINS, KSTAR_MIN, KSTAR_MAX, MULT_NBINS, MULT_MIN, MULT_MAX)
    
    if 'fCollisionId' in rdf.GetColumnNames():
        max_coll_id = int(rdf.Max('fCollisionId').GetValue())
        h_coll_id_vs_kstar_matter = (rdf.Filter('fSignHe3 > 0')
                                        .Histo2D(('hCollisionId_vs_Kstar_Matter', ';k* (GeV/c);Collision Id',
                                                KSTAR_NBINS, KSTAR_MIN, KSTAR_MAX, max_coll_id + 1, -0.5, max_coll_id + 0.5), 
                                                'fKstar', 'fCollisionId').GetValue())
        h_coll_id_vs_kstar_antimatter = (rdf.Filter('fSignHe3 < 0')
                                            .Histo2D(('hCollisionId_vs_Kstar_Antimatter', ';k* (GeV/c);Collision Id',
                                                    KSTAR_NBINS, KSTAR_MIN, KSTAR_MAX, max_coll_id + 1, -0.5, max_coll_id + 0.5), 
                                                    'fKstar', 'fCollisionId').GetValue())
    else: 
        h_coll_id_vs_kstar_matter = None
        h_coll_id_vs_kstar_antimatter = None
        
    h_z_vtx_vs_kstar_matter = (rdf.Filter('fSignHe3 > 0')
                                   .Histo2D(('hZvertex_vs_Kstar_Matter', ';k* (GeV/c);#it{z}_{}vtx} (cm)',
                                            KSTAR_NBINS, KSTAR_MIN, KSTAR_MAX, 200, -10, 10), 
                                                   'fKstar', 'fZVertex').GetValue())
    h_z_vtx_vs_kstar_antimatter = (rdf.Filter('fSignHe3 < 0')
                                    .Histo2D(('hZvertex_vs_Kstar_Antimatter', ';k* (GeV/c);#it{z}_{vtx} (cm)',
                                                KSTAR_NBINS, KSTAR_MIN, KSTAR_MAX, 200, -10, 10), 
                                                'fKstar', 'fZVertex').GetValue())

    # overall multiplicity: one entry per unique He3 (first occurrence only)
    rdf_matter_first = rdf.Filter('fSignHe3 > 0 && fHe3IsFirstOccurrence')
    rdf_antimatter_first = rdf.Filter('fSignHe3 < 0 && fHe3IsFirstOccurrence')

    mult_matter_vals = rdf_matter_first.Take['int']('fHe3Multiplicity').GetValue()
    mult_antimatter_vals = rdf_antimatter_first.Take['int']('fHe3Multiplicity').GetValue()

    for v in mult_matter_vals:
        h_mult_matter.Fill(v)
    for v in mult_antimatter_vals:
        h_mult_antimatter.Fill(v)

    # multiplicity vs kstar: one entry per row/pair
    rdf_matter = rdf.Filter('fSignHe3 > 0')
    rdf_antimatter = rdf.Filter('fSignHe3 < 0')

    # RDataFrame doesn't Take pairs directly, so grab the two columns separately and zip them
    kstar_matter_vals = rdf_matter.Take['float']('fKstar').GetValue()
    mult_matter_row_vals = rdf_matter.Take['int']('fHe3Multiplicity').GetValue()
    kstar_antimatter_vals = rdf_antimatter.Take['float']('fKstar').GetValue()
    mult_antimatter_row_vals = rdf_antimatter.Take['int']('fHe3Multiplicity').GetValue()

    for k, m in zip(kstar_matter_vals, mult_matter_row_vals):
        h_mult_vs_kstar_matter.Fill(k, m)
    for k, m in zip(kstar_antimatter_vals, mult_antimatter_row_vals):
        h_mult_vs_kstar_antimatter.Fill(k, m)

    output_file.cd()
    h_mult_matter.Write()
    h_mult_antimatter.Write()
    h_mult_vs_kstar_matter.Write()
    h_mult_vs_kstar_antimatter.Write()
    if h_coll_id_vs_kstar_matter:
        h_coll_id_vs_kstar_matter.Write()
    if h_coll_id_vs_kstar_antimatter:
        h_coll_id_vs_kstar_antimatter.Write()
    h_z_vtx_vs_kstar_matter.Write()
    h_z_vtx_vs_kstar_antimatter.Write()

    print(f'\t{tc.CYAN}hHe3Multiplicity_Matter{tc.RESET}: {h_mult_matter.GetEntries()} unique He3')
    print(f'\t{tc.CYAN}hHe3Multiplicity_Antimatter{tc.RESET}: {h_mult_antimatter.GetEntries()} unique He3')


def fill_local_multiplicity_vs_kstar_histograms(rdf, output_file: TFile):

    print(f'\n{tc.GREEN}Filling He3 multiplicity histograms (multiplicity computed within each k* bin){tc.RESET}')

    h_local_mult_vs_kstar_matter = ROOT.TH2F(
        'hHe3LocalMultiplicity_vs_Kstar_Matter',
        'He3 multiplicity within k* bin (matter);k* (GeV/c);N occurrences in this k* bin',
        KSTAR_NBINS, KSTAR_MIN, KSTAR_MAX, MULT_NBINS, MULT_MIN, MULT_MAX)
    h_local_mult_vs_kstar_antimatter = ROOT.TH2F(
        'hHe3LocalMultiplicity_vs_Kstar_Antimatter',
        'He3 multiplicity within k* bin (antimatter);k* (GeV/c);N occurrences in this k* bin',
        KSTAR_NBINS, KSTAR_MIN, KSTAR_MAX, MULT_NBINS, MULT_MIN, MULT_MAX)

    kstar_vals = rdf.Take['float']('fKstar').GetValue()
    zvtx_vals = rdf.Take['float']('fZVertex').GetValue()
    pt_vals = rdf.Take['float']('fPtHe3').GetValue()
    sign_vals = rdf.Take['float']('fSignHe3').GetValue()

    def kstar_bin_of(k):
        # mirrors the uniform binning of the TH2 above (KSTAR_NBINS, KSTAR_MIN, KSTAR_MAX)
        if k < KSTAR_MIN:
            return -1  # underflow
        if k >= KSTAR_MAX:
            return KSTAR_NBINS  # overflow
        return int((k - KSTAR_MIN) / (KSTAR_MAX - KSTAR_MIN) * KSTAR_NBINS)

    # how many times does each He3 (zvtx, pt) appear within each individual k* bin
    local_counts = Counter()
    for k, zvtx, pt in zip(kstar_vals, zvtx_vals, pt_vals):
        local_counts[(kstar_bin_of(k), zvtx, pt)] += 1

    # fill once per unique (He3, k* bin) combination: if a He3 appears N times in a given
    # k* bin, that bin's multiplicity histogram gets exactly one entry at N
    seen = set()
    for k, zvtx, pt, sign in zip(kstar_vals, zvtx_vals, pt_vals, sign_vals):
        key = (kstar_bin_of(k), zvtx, pt)
        if key in seen:
            continue
        seen.add(key)
        n_local = local_counts[key]
        if sign > 0:
            h_local_mult_vs_kstar_matter.Fill(k, n_local)
        else:
            h_local_mult_vs_kstar_antimatter.Fill(k, n_local)

    output_file.cd()
    h_local_mult_vs_kstar_matter.Write()
    h_local_mult_vs_kstar_antimatter.Write()


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Study He3 duplicate multiplicity for femtoscopic analysis')
    parser.add_argument('--config', type=str, default='config/config_prepare.yml', help='Path to the configuration YAML file')
    parser.add_argument('--mode', type=str, default='same_event', choices=['same_event', 'mixed_event'], help='Mode of analysis: same_event or mixed_event')
    parser.add_argument('--output', type=str, default='output_multiplicity.root', help='Output ROOT file name')
    args = parser.parse_args()

    config_file = args.config
    config = yaml.safe_load(open(config_file, 'r'))

    load_parametrisation(config)  # Load parametrisation into Common.h

    selections = config.get('selections', [])
    config = config[args.mode]  # Use the specific mode section from the config

    base_selection, selection = prepare_selections(selections)
    selection += f' && (fCentralityFT0C > 10 && fCentralityFT0C < 30)'
    chain_data, additional_chain_data = prepare_input_tchain(config)

    tmp_rdf = RDataFrame(chain_data)
    print(tc.GREEN+'\nDataset columns'+tc.RESET)
    print(tc.UNDERLINE+tc.CYAN+f'{tmp_rdf.GetColumnNames()}'+tc.RESET)
    zvtx_vals = tmp_rdf.Take["float"]("fZVertex").GetValue()
    pt_vals = tmp_rdf.Take["float"]("fPtHe3").GetValue()

    # count how many rows share the same (fZVertex, fPtHe3) key -> "same He3" multiplicity,
    # and record the first row index at which each key appears
    key_counts = Counter()
    first_occurrence_entry = {}
    for i, (zvtx, pt) in enumerate(zip(zvtx_vals, pt_vals)):
        key = (zvtx, pt)
        key_counts[key] += 1
        if key not in first_occurrence_entry:
            first_occurrence_entry[key] = i

    ROOT.gInterpreter.Declare("""
    #include <unordered_map>
    std::unordered_map<ULong64_t, int> entry_multiplicity;
    std::unordered_map<ULong64_t, bool> entry_is_first;
    int get_he3_multiplicity(ULong64_t entry) {
        auto it = entry_multiplicity.find(entry);
        return it != entry_multiplicity.end() ? it->second : 0;
    }
    bool is_first_occurrence(ULong64_t entry) {
        auto it = entry_is_first.find(entry);
        return it != entry_is_first.end() ? it->second : false;
    }
    """)

    first_entries = set(first_occurrence_entry.values())
    for i, (zvtx, pt) in enumerate(zip(zvtx_vals, pt_vals)):
        ROOT.entry_multiplicity[i] = key_counts[(zvtx, pt)]
        ROOT.entry_is_first[i] = (i in first_entries)

    rdf = prepare_rdataframe_with_multiplicity(chain_data, base_selection, selection)

    output_file_path = args.output
    output_file = ROOT.TFile(output_file_path, "RECREATE")

    fill_multiplicity_histograms(rdf, output_file)
    fill_local_multiplicity_vs_kstar_histograms(rdf, output_file)

    output_file.Close()