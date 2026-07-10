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

#NBINS_DCA = 60
#DCA_BINNING = np.concatenate((np.arange(-0.1, -0.05, 0.005), np.arange(-0.05, 0.05, 0.0025), np.arange(0.05, 0.1, 0.005)))
DCA_BINNING = np.concatenate((np.arange(-0.1, -0.05, 0.005), np.arange(-0.05, -0.025, 0.0025), np.arange(-0.025, 0.025, 0.00125), np.arange(0.025, 0.05, 0.0025), np.arange(0.05, 0.1, 0.005)))
NBINS_DCA = len(DCA_BINNING) - 1

def prepare_selections(config):
    
    selections = config['selections']
    
    selection = selections[0]
    for sel in selections[1:]:
        selection += (' && ' + sel)

    return  selection

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
    
def visualise(rdf: RDataFrame, outfile:TDirectory, particle:str = 'He'):

    particle_suffix = 'He3' if particle == 'He' else 'Had'

    h_pt = rdf.Histo1D((f'hPt{particle}', ';#it{p}_{T} (GeV/#it{c})', 200, -10, 10), f'signedPt')

    h2_dcaxy_pt = rdf.Histo2D((f'h2DCAxyPt{particle}', ';#it{p}_{T} (GeV/#it{c});DCA_{#it{xy}} (cm)',
                               NBINS_PT, PT_BINNING, NBINS_DCA, DCA_BINNING), f'signedPt', f'fDCAxy')
                               #100, -5, 5, 120, -0.15, 0.15), f'signedPt', f'fDCAxy')
    h2_dcaz_pt = rdf.Histo2D((f'h2DCAzPt{particle}', ';#it{p}_{T} (GeV/#it{c});DCA_{#it{z}} (cm)', 
                               NBINS_PT, PT_BINNING, NBINS_DCA, DCA_BINNING), f'signedPt', f'fDCAz')
                               #100, -5, 5, 120, -0.3, 0.3), f'signedPt', f'fDCAz')
                               
    adjust_bin_content_for_variable_bin_width(h2_dcaxy_pt)
    adjust_bin_content_for_variable_bin_width(h2_dcaz_pt)
    
    
    histos = []
    histos.append(h2_dcaxy_pt)
    histos.append(h2_dcaz_pt)
    histos.append(h_pt)

    outfile.cd()
    for hist in histos:
        hist.Write()

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

def prepare_rdataframe(chain_data: TChain, selection: str):
   
    rdf = RDataFrame(chain_data)
    print(tc.GREEN+'\nDataset columns'+tc.RESET)
    print(tc.UNDERLINE+tc.CYAN+f'{rdf.GetColumnNames()}'+tc.RESET)
    
    rdf = (rdf.Define("ptUncorr", "2 * std::abs(fPt)")
            .Define("matter", "fPt > 0")
            .Define("pidForTracking", "fFlags >> 12")
            .Define("pt", '(pidForTracking == 7) || (pidForTracking == 8) || (ptUncorr > 2.5) ? ptUncorr : CorrectPidTrkHe(ptUncorr)') \
            .Define("signedPt", "pt / fPt > 0 ? pt : -pt")
            .Define("nsigmaHe3", "ComputeNsigmaTPCHe(std::abs(fTPCInnerParam) * 2, fTPCsignal, false, true)")
            .Define("avgClSizeIts", "ComputeAverageClusterSize(fITSclusterSizes)")
            .Define("clSizeItsCosLam", "avgClSizeIts / std::cosh(fEta)")
            .Define("nsigmaItsHe3", "ComputeNsigmaITSHe(std::abs(pt), clSizeItsCosLam)")
            .Define("nITSclsIB", "int(0) + bool(fITSclsMap & 1) + bool(fITSclsMap & 2) + bool(fITSclsMap & 4)")
            .Define("nITScls", "nITSclsIB + bool(fITSclsMap & 8) + bool(fITSclsMap & 16) + bool(fITSclsMap & 32) + bool(fITSclsMap & 64)")
            .Filter(selection)
            )
    
    return rdf

def main():
    
    config = {
                'input_data': [ '/data/galucia/he4/PbPb/LHC23_pass5_centrality.root',
                                ],
                'tree_name': 'O2nucleitable',
                'mode': 'DF',
                
                'parameterisation': {

                    # PbPb 23 pass5
                    'kHeTPCParams': [-335.570, 0.7066, 1.5731, 0.5202, 2.8735],
                    'kHeTPCResolution': 0.060,
                    'kHeTPCParamsResiduals': [31.8927, 0.6830, 0.0750, -8.3483, 0.9692, 0.1287],

                    # 24_pass3
                    ### 'kHeTPCParams': [-200.84, 0.3000, 1.6669, 1.1523, 2.8760],
                    ### 'kHeTPCResolution': 0.059,
                    ### 'kHeTPCResidualType': 'exp',
                    ### 'kHeTPCParamsResiduals': [-162.97, 11.118, 0.4094, 0., 0., 0.],

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
                "fTPCnCls >= 110",
                "nITScls >= 5",
                "std::abs(fEta) < 0.9",
                "nsigmaItsHe3 > -2",
                "((pt < 2.5) || (pidForTracking == 7) || (pidForTracking == 8))",
                "nsigmaHe3 > -2",
                "nsigmaHe3 < 3",
                ]}
    
    load_parametrisation(config)  # Load parametrisation into Common.h
    selections = prepare_selections(selections_dict)
    rdf = prepare_rdataframe(chain_data, selections)

    outfile = TFile.Open('output/dca/dca_data_template_nuclei_spectra.root', 'RECREATE')
    
    for particle_name, particle in zip(['He3'], ['He']):

        outdir = outfile.mkdir(f'{particle}')
        visualise(rdf, outdir, particle)
        
    outfile.Close()


if __name__ == '__main__':

    main()
