import numpy as np
from ROOT import TF1, TH1F, TFile, TCanvas, gInterpreter
from torchic.physics.simulations import RunTwoBodyDecaySimulation
from torchic.core.histogram import load_hist

from particle import Particle

import sys
sys.path.append('..')
from utils.particles import ParticleMasses

gInterpreter.ProcessLine(f'#include "../include/ThermalModels.h"')
from ROOT import ThermalModels

def sample_blast_wave(mass:float, outfile:TFile, n_samples:int = 1_000_000, hist_name:str='h_blast_wave',
                      normalisation:float=1.):

    # values taken form the Blast-Wave fit in https://arxiv.org/abs/2405.19839
    T_kin = 0.103 # GeV
    beta_s = 0.694
    n = 0.498
    
    thermal_model = ThermalModels()
    blast_wave = thermal_model.GetBGBW(mass, beta_s, T_kin, n, normalisation)
    
    blast_wave.SetTitle(f'Blast-Wave, m = {mass} GeV/#it{{c}}^{{2}}; #it{{p}}_{{T}} (GeV/#it{{c}}); d^{{2}}#it{{N}}/d#it{{p}}_{{T}}d#it{{y}}')

    h_blast_wave = TH1F(hist_name, 
                        f'Blast-Wave, m = {mass} GeV/#it{{c}}^{{2}}; #it{{p}}_{{T}} (GeV/#it{{c}}); d^{{2}}#it{{N}}/d#it{{p}}_{{T}}d#it{{y}}',
                        2000, 0, 20)
    for isample in range(n_samples):
        h_blast_wave.Fill(blast_wave.GetRandom())
    h_blast_wave.Scale(normalisation)

    outfile.cd()
    blast_wave.Write()
    h_blast_wave.Write()

    return h_blast_wave

if __name__ == '__main__':

    outfile = TFile.Open('output/hypertriton_spectra.root', 'recreate')
    
    mass_hypertriton = 2.991134 # GeV
    mass_he3 = ParticleMasses['He']

    h_hypertriton = load_hist('input/HEPData-ins2791616-v1-root.root', 
                              '(Anti)hypertriton spectrum in 0-10% V0M centrality class/Hist1D_y1')
    h_helium3 = load_hist('input/He3_spectrum.root', 'hHe30_10')

    h_blast_wave_hypertriton = sample_blast_wave(mass_hypertriton, outfile, 10_000_000, 'hBlastWaveHypertriton',
                                                 normalisation=h_hypertriton.Integral())
    h_blast_wave_he3 = sample_blast_wave(mass_he3, outfile, 10_000_000, 'hBlastWaveHe3',
                                         normalisation=h_helium3.Integral())
    
    outfile.Close()
