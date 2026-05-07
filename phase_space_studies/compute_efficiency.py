'''
    Study of the dca distribution for p-He3 pairs from MC
'''

import sys
from enum import Enum

from ROOT import TFile, TChain, RDataFrame, TDirectory, \
    gInterpreter, gROOT, EnableImplicitMT
    #kPPrimary

from torchic.core.histogram import build_efficiency
from torchic.utils.terminal_colors import TerminalColors as tc

gInterpreter.ProcessLine(f'#include "../include/Common.h"')

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

def visualise(rdf: RDataFrame, outfile:TDirectory, particle:str = 'He3'):

    h_pt_rec = rdf.Histo1D((f'hPtRec{particle}', 'Total;#it{p}_{T}^{rec} (GeV/#it{c});',
                                 100, -10, 10), 'fPt').GetValue()
    
    h_pt_gen = rdf.Histo1D((f'hPtGen{particle}', 'Total;#it{p}_{T}^{gen} (GeV/#it{c});',
                                 100, -10, 10), 'fgPt').GetValue()

    h_pt_rec_prim = rdf.Filter('fIsPrimary == true') \
                       .Histo1D((f'hPtRecPrimary{particle}', 'Primary;#it{p}_{T}^{rec} (GeV/#it{c});',
                                 100, -10, 10), 'fPt').GetValue()
    
    h_pt_gen_prim = rdf.Filter('fIsPrimary == true') \
                       .Histo1D((f'hPtGenPrimary{particle}', 'Primary;#it{p}_{T}^{gen} (GeV/#it{c});',
                                 100, -10, 10), 'fgPt').GetValue()
    
    h_pt_rec_hypertriton = rdf.Filter('fIsFromHypertriton == true') \
                              .Histo1D((f'hPtRecFromHypertriton{particle}', 'Secondary from ^{3}_{#Lambda}H;#it{p}_{T}^{rec} (GeV/#it{c});',
                                        100, -10, 10), 'fPt').GetValue()
    
    h_pt_gen_hypertriton = rdf.Filter('fIsFromHypertriton == true') \
                              .Histo1D((f'hPtGenFromHypertriton{particle}', 'Secondary from ^{3}_{#Lambda}H;#it{p}_{T}^{gen} (GeV/#it{c});',
                                        100, -10, 10), 'fgPt').GetValue()

    h_pt_rec_material = rdf.Filter('fIsFromMaterial == true') \
                              .Histo1D((f'hPtRecFromMaterial{particle}', 'Secondary from material;#it{p}_{T}^{rec} (GeV/#it{c});',
                                        100, -10, 10), 'fPt').GetValue()
    
    h_pt_gen_material = rdf.Filter('fIsFromMaterial == true') \
                              .Histo1D((f'hPtGenFromMaterial{particle}', 'Secondary from material;#it{p}_{T}^{gen} (GeV/#it{c});',
                                        100, -10, 10), 'fgPt').GetValue()
    
    h_efficiency = build_efficiency(h_pt_gen, h_pt_rec, name='h_efficiency_he3', ytitle='Efficiency')
    
    h_efficiency_prim = build_efficiency(h_pt_gen_prim, h_pt_rec_prim,
                                         name='h_efficiency_he3_prim', ytitle='Efficiency')

    h_efficiency_he3_from_h3l = build_efficiency(h_pt_gen_hypertriton, h_pt_rec_hypertriton,
                                                 name='h_efficiency_he3_from_h3l', ytitle='Efficiency')

    h_efficiency_material = build_efficiency(h_pt_gen_material, h_pt_rec_material,
                                         name='h_efficiency_he3_material', ytitle='Efficiency')
    
    outfile.cd()
    for hist in [h_pt_rec, h_pt_gen, h_efficiency,
                 h_pt_rec_prim, h_pt_gen_prim, h_efficiency_prim, 
                 h_pt_rec_hypertriton, h_pt_gen_hypertriton, h_efficiency_he3_from_h3l, 
                 h_pt_rec_material, h_pt_gen_material, h_efficiency_material]:
        hist.Write()


def main():
    
    input_files = { 'Pr': [#'/data/galucia/lithium_local/MC/nucleiqc/LHC25a4_protons.root',
                           '/data/galucia/lithium_local/MC/nucleiqc/LHC24e2d_protons.root'],
                   'He': [#'/data/galucia/lithium_local/MC/nucleiqc/LHC25a4_he3.root',
                          '/data/galucia/lithium_local/MC/nucleiqc/LHC24i5_he3.root',]}
    #input_files = {'Pr': '/home/galucia/Lithium4/task/nucleiQC/local_test.root',
    #               'He': '/home/galucia/Lithium4/task/nucleiQC/local_test.root',}
    tree_name = 'O2nucleitablered'

    chain_data = {'Pr': TChain('tchainPr'),
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

    mothers_from_material = [211,321, 130, 22, 2212, 2112, 1000100200, 1000070150, 1000010020, 1000010030, 1000010040,
                             1000010050, 1000020040, 1000020050, 1000020060, 1000040060, 1000050080,]
    required_mother_from_material = ''
    for imother, mother in enumerate(mothers_from_material):
        required_mother_from_material += f'std::abs(fMotherPDGcode) == {mother}'
        if imother < len(mothers_from_material)-1:
            required_mother_from_material += ' || '
    required_mother_from_material = '(' + required_mother_from_material + ')'

    rdfs = {}
    PDGcode = {'Pr': 2212, 'He': 1000020030}

    kPPrimary = 0

    for particle in ['He', 'Pr']:

        print(f'Processing {tc.CYAN+tc.UNDERLINE}{particle}{tc.RESET}...')

        isHe3 = 'true ' if (particle == 'He') else 'false'
        rdfs[particle] = RDataFrame(chain_data[particle]) \
                           .Filter(f'std::abs(fPDGcode) == {PDGcode[particle]}') \
                           .Define('fPidForTracking', 'ReadPidTrkFromFlags(fFlags)') \
                           .Redefine('fPt', f'{isHe3} ? 2 * fPt : fPt') \
                           .Redefine('fPt', f'{isHe3} ? ((fPidForTracking == 7) || (fPidForTracking == 8) || (fPt > 2.5) ? fPt : CorrectPidTrkHe(fPt)) : fPt') \
                           .Define('fIsPrimary', f'fMcProcess == {kPPrimary}') \
                           .Define('fIsFromHypertriton', f'std::abs(fMotherPDGcode) == {ParticlePDG["Hypertriton"]}') \
                           .Define('fIsFromMaterial', f'ReadBitFromFlags(fFlags, {Flags.kIsSecondaryFromMaterial.value}) && std::abs(fMotherPDGcode) != {ParticlePDG["Li4"]}') # && {required_mother_from_material}') \

        print(f'{rdfs[particle].GetColumnNames()=}')

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
    
    outfile = TFile.Open('output/single_track_efficiency.root', 'RECREATE')
    
    for particle in ['He', 'Pr']:
        outdir = outfile.mkdir(f'{particle}')
        visualise(rdfs[particle], outdir, particle)
        
    outfile.Close()


if __name__ == '__main__':

    main()
