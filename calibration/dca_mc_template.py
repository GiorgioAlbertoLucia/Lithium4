'''
    Study of the dca distribution for p-He3 pairs from MC
'''

import sys
from enum import Enum
import numpy as np

from ROOT import TFile, TChain, RDataFrame, TDirectory, \
    gInterpreter, gROOT, EnableImplicitMT

from torchic.utils.terminal_colors import TerminalColors as tc

gInterpreter.ProcessLine(f'#include "../include/Common.h"')
from ROOT import ComputeNsigmaTOFPr, ComputeAverageClusterSize, CorrectPidTrkHe, ReadPidTrkFromFlags, ReadBitFromFlags

EnableImplicitMT(5)
gROOT.SetBatch(True)

sys.path.append('..')
from utils.particles import ParticlePDG

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

def visualise(rdf: RDataFrame, outfile:TDirectory, particle:str = 'He3'):

    h_pt = rdf.Histo1D((f'hPt{particle}', ';#it{p}_{T} (GeV/#it{c})', 200, -10, 10), 'fPt')
    
    h2_dcaxy_pt = rdf.Histo2D((f'h2DCAxyPt{particle}', ';#it{p}_{T} (GeV/#it{c});DCA_{#it{xy}} (cm)', 
                               100, -5, 5, 120, -0.15, 0.15), 'fPt', 'fDCAxy')
                               #NBINS_PT, PT_BINNING, NBINS_DCA, DCA_BINNING), 'fPt', 'fDCAxy')
    h2_dcaz_pt = rdf.Histo2D((f'h2DCAzPt{particle}', ';#it{p}_{T} (GeV/#it{c});DCA_{#it{z}} (cm)', 
                               100, -5, 5, 120, -0.3, 0.3), 'fPt', 'fDCAz')
    
    histos = []
    histos.append(h2_dcaxy_pt)
    histos.append(h2_dcaz_pt)
    histos.append(h_pt)

    for flag in ['IsPhysicalPrimary', 'IsSecondaryFromMaterial', 'IsSecondaryFromWeakDecay', 'IsFromLi4', 'IsFromSigmaPlus', 'IsFromLambda0']:
        h2_dcaxy_pt = rdf.Filter(f'f{flag} == true') \
                         .Histo2D((f'h2DCAxyPt{particle}_{flag}', ';#it{p}_{T} (GeV/#it{c});DCA_{#it{xy}} (cm)', 
                                   200, -10, 10, 120, -0.15, 0.15), 'fPt', 'fDCAxy')
                                   #NBINS_PT, PT_BINNING, NBINS_DCA, DCA_BINNING), 'fPt', 'fDCAxy')
        h2_dcaz_pt = rdf.Filter(f'f{flag} == true') \
                         .Histo2D((f'h2DCAzPt{particle}_{flag}', ';#it{p}_{T} (GeV/#it{c});DCA_{#it{z}} (cm)', 
                                   200, -10, 10, 120, -0.3, 0.3), 'fPt', 'fDCAz')
        
        histos.append(h2_dcaxy_pt)
        histos.append(h2_dcaz_pt)

        if flag in ['IsFromSigmaPlus', 'IsFromLambda0']:
            h_pt = rdf.Filter(f'f{flag} == true') \
                    .Histo1D((f'hPt{particle}_{flag}', ';#it{p}_{T} (GeV/#it{c})', 200, -10, 10), 'fPt')
            histos.append(h_pt)
    
    outfile.cd()

    for hist in histos:
        hist.Write()


def main():
    
    input_files = {'Pr': [#'/data/galucia/lithium_local/MC/nucleiqc/LHC25a4_protons.root'
                            #'/data/galucia/lithium_local/MC/nucleiqc/LHC24e2d_protons.root'
                            '/data/galucia/lithium_local/MC/alimonitor/nucleiQC/protons/train_LHC24g3/AO2D.root'],
                    'Pr_as_He': [#'/data/galucia/lithium_local/MC/nucleiqc/LHC25a4_protons.root'
                            #'/data/galucia/lithium_local/MC/nucleiqc/LHC24e2d_protons.root'
                            '/data/galucia/lithium_local/MC/alimonitor/nucleiQC/protons/train_LHC24g3/AO2D.root'],
                   'De_as_He': ['/data/galucia/lithium_local/MC/alimonitor/nucleiQC/deuterons/train_LHC24g3/AO2D.root'],
                   'He': ['/data/galucia/lithium_local/MC/nucleiqc/LHC25a4_he3.root',
                          '/data/galucia/lithium_local/MC/nucleiqc/LHC24i5_he3.root',]}
    #input_files = {'Pr': '/home/galucia/Lithium4/task/nucleiQC/local_test.root',
    #               'He': '/home/galucia/Lithium4/task/nucleiQC/local_test.root',}
    tree_name = 'O2nucleitablered'

    chain_data = {'Pr': TChain('tchainPr'),
                    'Pr_as_He': TChain('tchainPr_as_He'),
                  'De_as_He': TChain('tchainDe'),
                  'He': TChain('tchainHe')}
    
    fileDatas = {}
    for particle, file_names in input_files.items():
        for file_name in file_names:
            fileDatas[particle] = TFile(file_name)

        for key in fileDatas[particle].GetListOfKeys():
          key_name = key.GetName()
          if 'DF_' in key_name :
              print(f'Adding {tc.CYAN+tc.UNDERLINE}{file_name}/{key_name}/{tree_name}{tc.RESET} to the chain')
              chain_data[particle].Add(f'{file_name}/{key_name}/{tree_name}')
 
    print(f'{Flags.kIsPhysicalPrimary.value=}, {Flags.kIsSecondaryFromMaterial.value=}, {Flags.kIsSecondaryFromWeakDecay.value=}')

    mothers_from_material = [211, 321, 130, 22, 2212, 2112, 1000100200, 1000070150, 1000010020, 1000010030, 1000010040,
                             1000010050, 1000020040, 1000020050, 1000020060, 1000040060, 1000050080,]
    required_mother_from_material = ''
    for imother, mother in enumerate(mothers_from_material):
        required_mother_from_material += f'std::abs(fMotherPDGcode) == {mother}'
        if imother < len(mothers_from_material)-1:
            required_mother_from_material += ' || '
    required_mother_from_material = '(' + required_mother_from_material + ')'

    rdfs = {}
    PDGcode = {'Pr': 2212, 'Pr_as_He': 2212, 'De_as_He': 1000010020, 'He': 1000020030}
    #.Define('fIsSecondaryFromMaterial', f'ReadBitFromFlags(fFlags, {Flags.kIsSecondaryFromMaterial.value}) && std::abs(fMotherPDGcode) != {ParticlePDG["Li4"]} && {required_mother_from_material}') \
    #.Define('fIsSecondaryFromWeakDecay', f'fMcProcess != 0 && fMcProcess != 4 && std::abs(fMotherPDGcode) != {ParticlePDG["Li4"]}') \

    # TMCProcess:
    # kPPrimary = 0
    # kPDecay = 4
    for particle in ['Pr', 'Pr_as_He', 'De_as_He', 'He']:
        rdfs[particle] = RDataFrame(chain_data[particle]) \
                           .Filter(f'std::abs(fPDGcode) == {PDGcode[particle]}') \
                           .Define('fPidForTracking', 'ReadPidTrkFromFlags(fFlags)') \
                           .Define('fIsPhysicalPrimary', f'ReadBitFromFlags(fFlags, {Flags.kIsPhysicalPrimary.value}) && std::abs(fMotherPDGcode) != {ParticlePDG["Li4"]}') \
                           .Define('fIsSecondaryFromWeakDecay', f'ReadBitFromFlags(fFlags, {Flags.kIsSecondaryFromWeakDecay.value}) && std::abs(fMotherPDGcode) != {ParticlePDG["Li4"]}') \
                           .Define('fIsSecondaryFromMaterial', f'ReadBitFromFlags(fFlags, {Flags.kIsSecondaryFromMaterial.value}) && std::abs(fMotherPDGcode) != {ParticlePDG["Li4"]}') \
                           .Define('fIsFromLi4', f'std::abs(fMotherPDGcode) == {ParticlePDG["Li4"]}') \
                           .Define('fIsFromSigmaPlus', f'std::abs(fMotherPDGcode) == {ParticlePDG["Sigma+"]}') \
                           .Define('fIsFromLambda0', f'std::abs(fMotherPDGcode) == {ParticlePDG["Lambda0"]}')

        if particle == 'He' or particle == 'De_as_He' or particle == 'Pr_as_He':
            rdfs[particle] = rdfs[particle].Redefine('fPt', 'fPt * 2')

        if False:   # check the content - mother pdg code and mc process
            for boolean_variable in ['fIsSecondaryFromMaterial']:#, 'fIsPhysicalPrimary', 'fIsSecondaryFromWeakDecay']:

                filtered_rdf = rdfs[particle].Filter(f'{boolean_variable} == true')
                col1ValsRP = filtered_rdf.Take["Int_t"]("fMotherPDGcode")
                #col2ValsRP = filtered_rdf.Take["Int_t"]("fPDGcode")
                #col3ValsRP = filtered_rdf.Take["ULong64_t"]("fMcProcess")
                uniqueCol1Vals = list(set(col1ValsRP.GetValue()))
                #uniqueCol2Vals = list(set(col2ValsRP.GetValue()))
                #uniqueCol3Vals = list(set(col3ValsRP.GetValue()))

                #subdfs = [rdfs[particle].Filter(f'fMotherPDGcode == {uniqueColVal}') for uniqueColVal in uniqueCol1Vals]

                #evtsPerCol1ValsRP = [subdf.Count().GetValue() for subdf in subdfs]

                print(f'\n{boolean_variable}')
                #print(f'unique pdg codes = {uniqueCol2Vals}')
                print(f'unique mother pdg codes = {uniqueCol1Vals}')
                #print(f'unique mc process = {uniqueCol3Vals}')
                #for uniqueVal, nUniqueVal in zip(uniqueCol1Vals, evtsPerCol1ValsRP):
                #    print(f'{particle} fMotherPDGCode: {uniqueVal}, n entries: {nUniqueVal}.')
    
    outfile = TFile.Open('output/dca_mc_template_check.root', 'RECREATE')
    
    for particle in ['Pr', 'Pr_as_He', 'De_as_He', 'He']:
        outdir = outfile.mkdir(f'{particle}')
        visualise(rdfs[particle], outdir, particle)
        
    outfile.Close()


if __name__ == '__main__':

    main()
