'''
    Study of the dca distribution for p-He3 pairs from MC
'''

from enum import Enum

from ROOT import TFile, TChain, RDataFrame, TDirectory, \
    gInterpreter, gROOT, EnableImplicitMT

from torchic.utils.terminal_colors import TerminalColors as tc

gInterpreter.ProcessLine(f'#include "../include/Common.h"')
from ROOT import ComputeNsigmaTOFPr, ComputeAverageClusterSize, CorrectPidTrkHe, ReadPidTrkFromFlags, ReadBitFromFlags

EnableImplicitMT(5)
gROOT.SetBatch(True)

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

def visualise(rdf: RDataFrame, outfile:TDirectory, particle:str = 'He3'):

    
    h2_dcaxy_pt = rdf.Histo2D((f'h2DCAxyPt{particle}', ';#it{p}_{T} (GeV/#it{c});DCA_{#it{xy}} (cm)', 
                               100, -5, 5, 150, -0.15, 0.15), 'fPt', 'fDCAxy')
    h2_dcaz_pt = rdf.Histo2D((f'h2DCAzPt{particle}', ';#it{p}_{T} (GeV/#it{c});DCA_{#it{z}} (cm)', 
                               100, -5, 5, 150, -0.3, 0.3), 'fPt', 'fDCAz')
    
    histos = []
    histos.append(h2_dcaxy_pt)
    histos.append(h2_dcaz_pt)

    for flag in ['IsPhysicalPrimary', 'IsSecondaryFromMaterial', 'IsSecondaryFromWeakDecay']:
        h2_dcaxy_pt = rdf.Filter(f'f{flag} == true') \
                         .Histo2D((f'h2DCAxyPt{particle}_{flag}', ';#it{p}_{T} (GeV/#it{c});DCA_{#it{xy}} (cm)', 
                                   100, -5, 5, 150, -0.15, 0.15), 'fPt', 'fDCAxy')
        h2_dcaz_pt = rdf.Filter(f'f{flag} == true') \
                         .Histo2D((f'h2DCAzPt{particle}_{flag}', ';#it{p}_{T} (GeV/#it{c});DCA_{#it{z}} (cm)', 
                                   100, -5, 5, 150, -0.3, 0.3), 'fPt', 'fDCAz')
        histos.append(h2_dcaxy_pt)
        histos.append(h2_dcaz_pt)
    
    outfile.cd()
    for hist in histos:
        hist.Write()


def main():
    
    #input_files = {'Pr': '/data/galucia/lithium_local/MC/nucleiqc/LHC25a4_protons.root',
    #               'He': '/data/galucia/lithium_local/MC/nucleiqc/LHC25a4_he3.root',}
    input_files = {'Pr': '/home/galucia/Lithium4/task/nucleiQC/mixed_event_test.root',
                   'He': '/home/galucia/Lithium4/task/nucleiQC/mixed_event_test.root',}
    tree_name = 'O2nucleitablered'

    chain_data = {'Pr': TChain('tchainPr'),
                  'He': TChain('tchainHe')}
    
    fileDatas = {}
    for particle, file_name in input_files.items():
        fileDatas[particle] = TFile(file_name)

        for key in fileDatas[particle].GetListOfKeys():
          key_name = key.GetName()
          if 'DF_' in key_name :
              print(f'Adding {tc.CYAN+tc.UNDERLINE}{file_name}/{key_name}/{tree_name}{tc.RESET} to the chain')
              chain_data[particle].Add(f'{file_name}/{key_name}/{tree_name}')
 
    print(f'{Flags.kIsPhysicalPrimary.value=}, {Flags.kIsSecondaryFromMaterial.value=}, {Flags.kIsSecondaryFromWeakDecay.value=}')

    rdfs = {}    
    PDGcode = {'Pr': 2212, 'He': 1000020030}
    for particle in ['Pr', 'He']:
        rdfs[particle] = RDataFrame(chain_data[particle]) \
                           .Filter(f'std::abs(fPDGcode) == {PDGcode[particle]}') \
                           .Define('fPidForTracking', 'ReadPidTrkFromFlags(fFlags)') \
                           .Define('fIsPhysicalPrimary', f'ReadBitFromFlags(fFlags, {Flags.kIsPhysicalPrimary.value})') \
                           .Define('fIsSecondaryFromMaterial', f'ReadBitFromFlags(fFlags, {Flags.kIsSecondaryFromMaterial.value})') \
                           .Define('fIsSecondaryFromWeakDecay', f'ReadBitFromFlags(fFlags, {Flags.kIsSecondaryFromWeakDecay.value})')

        col1ValsRP = rdfs[particle].Take["Int_t"]("fMotherPDGcode")
        col2ValsRP = rdfs[particle].Take["Int_t"]("fPDGcode")
        uniqueCol1Vals = list(set(col1ValsRP.GetValue()))
        uniqueCol2Vals = list(set(col2ValsRP.GetValue()))

        #subdfs = [rdfs[particle].Filter(f'fMotherPDGcode == {uniqueColVal}') for uniqueColVal in uniqueCol1Vals]

        #evtsPerCol1ValsRP = [subdf.Count().GetValue() for subdf in subdfs]

        print()
        print(f'unique pdg codes = {uniqueCol2Vals}')
        print(f'unique mother pdg codes = {uniqueCol1Vals}')
        #for uniqueVal, nUniqueVal in zip(uniqueCol1Vals, evtsPerCol1ValsRP):
        #    print(f'{particle} fMotherPDGCode: {uniqueVal}, n entries: {nUniqueVal}.')
    
    outfile = TFile.Open('output/dca_mc_template_check.root', 'RECREATE')
    
    for particle in ['Pr', 'He']:
        outdir = outfile.mkdir(f'{particle}')
        visualise(rdfs[particle], outdir, particle)
        
    outfile.Close()


if __name__ == '__main__':

    main()
