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

def sample_boltzmann(mass:float, outfile:TFile, n_samples:int = 1_000_000):

    boltzmann = TF1('boltzmann', 
                    'std::sqrt([1]*[1] + x*x) * x * std::exp(- std::sqrt([1]*[1] + x*x) / [0] )',
                    0, 20, 2)
    boltzmann.SetTitle(f'Boltzmann, m = {mass} GeV/#it{{c}}^{{2}}; #it{{p}}_{{T}} (GeV/#it{{c}}); d^{{2}}#it{{N}}/d#it{{p}}_{{T}}d#it{{y}}')
    boltzmann.SetParameter(1, mass)
    
    # values taken form the Blast-Wave fit in https://arXiv.org/abs/2311.11758
    T_kin = 0.132 # GeV
    beta_s = 0.670
    T_slope = T_kin * np.sqrt( (1 + beta_s) / (1 - beta_s))
    boltzmann.SetParameter(0, T_slope)

    canvas = TCanvas('canavs', boltzmann.GetTitle())
    canvas.SetLogy()
    
    canvas.DrawFrame(0., 1e-20, 20., 1e-5, boltzmann.GetTitle())
    boltzmann.Draw('same')

    h_boltzmann = TH1F('h_boltzmann', 
                                      f'Boltzmann, m = {mass} GeV/#it{{c}}^{{2}}; #it{{p}}_{{T}} (GeV/#it{{c}}); d^{{2}}#it{{N}}/d#it{{p}}_{{T}}d#it{{y}}',
                                      200, 0, 20)
    for isample in range(n_samples):
        h_boltzmann.Fill(boltzmann.GetRandom())

    outfile.cd()
    boltzmann.Write()
    canvas.Write()
    h_boltzmann.Write()
    canvas.SaveAs('output/boltzmann.pdf')

    return h_boltzmann

def sample_blast_wave(mass:float, outfile:TFile, n_samples:int = 1_000_000):

    # values taken form the Blast-Wave fit in https://arXiv.org/abs/2311.11758
    T_kin = 0.132 # GeV
    beta_s = 0.670
    n = 0.55
    normalisation = 1.
    
    thermal_model = ThermalModels()
    blast_wave = thermal_model.GetBGBW(mass, beta_s, T_kin, n, normalisation)
    
    blast_wave.SetTitle(f'Blast-Wave, m = {mass} GeV/#it{{c}}^{{2}}; #it{{p}}_{{T}} (GeV/#it{{c}}); d^{{2}}#it{{N}}/d#it{{p}}_{{T}}d#it{{y}}')
    canvas = TCanvas('canavs', blast_wave.GetTitle())
    canvas.SetLogy()
    
    canvas.DrawFrame(0., 1e-20, 20., 1e-5, blast_wave.GetTitle())
    blast_wave.Draw('same')

    h_blast_wave = TH1F('h_blast_wave', 
                                      f'Blast-Wave, m = {mass} GeV/#it{{c}}^{{2}}; #it{{p}}_{{T}} (GeV/#it{{c}}); d^{{2}}#it{{N}}/d#it{{p}}_{{T}}d#it{{y}}',
                                      200, 0, 20)
    for isample in range(n_samples):
        h_blast_wave.Fill(blast_wave.GetRandom())

    outfile.cd()
    blast_wave.Write()
    canvas.Write()
    h_blast_wave.Write()
    canvas.SaveAs('output/blast_wave.pdf')

    return h_blast_wave

def correct_by_efficiency(h_dNdpt:TH1F, outfile:TFile):

    h_efficiency_matter = load_hist('/home/galucia/Lithium4/preparation/output/efficiency.root', 'hEfficiencyMatter')
    h_efficiency_antimatter = load_hist('/home/galucia/Lithium4/preparation/output/efficiency.root', 'hEfficiencyAntimatter')

    h_dNdpt_matter = h_dNdpt.Clone(h_dNdpt.GetName()+'_matter')
    h_dNdpt_antimatter = h_dNdpt.Clone(h_dNdpt.GetName()+'_antimatter')

    for ibin in range(1, h_dNdpt.GetNbinsX()+1):
        h_dNdpt_matter.SetBinContent(ibin, h_dNdpt.GetBinContent(ibin)*h_efficiency_matter.GetBinContent(ibin))
        h_dNdpt_antimatter.SetBinContent(ibin, h_dNdpt.GetBinContent(ibin)*h_efficiency_antimatter.GetBinContent(ibin))

    outfile.cd()
    h_efficiency_matter.Write()
    h_efficiency_antimatter.Write()
    h_dNdpt_matter.Write()
    h_dNdpt_antimatter.Write()

    return h_dNdpt_matter, h_dNdpt_antimatter




if __name__ == '__main__':

    outfile = TFile.Open('output/spectra.root', 'recreate')
    
    mass_li4 = 3.751 # GeV/c^2
    mass_he3 = ParticleMasses['He']
    mass_p = ParticleMasses['Pr']

    h_boltzmann = sample_boltzmann(mass_li4, outfile, 10_000_000)
    h_boltzmann_matter, h_boltzmann_antimatter = \
                            correct_by_efficiency(h_boltzmann, outfile)
    
    h_blast_wave = sample_blast_wave(mass_li4, outfile, 10_000_000)
    h_blast_wave_matter, h_blast_wave_antimatter = \
                            correct_by_efficiency(h_blast_wave, outfile)

    outfile.Close()

    RunTwoBodyDecaySimulation('output/spectra.root', 'h_boltzmann_matter',
                              'output/lifetime_li4.root', 'hDecayTimeLi4',
                              'output/phase_space_studies_matter.root',
                              mass_li4, mass_he3, mass_p)
    
    RunTwoBodyDecaySimulation('output/spectra.root', 'h_boltzmann_antimatter',
                              'output/lifetime_li4.root', 'hDecayTimeLi4',
                              'output/phase_space_studies_antimatter.root',
                              mass_li4, mass_he3, mass_p)
    
    RunTwoBodyDecaySimulation('output/spectra.root', 'h_blast_wave_matter',
                              'output/lifetime_li4.root', 'hDecayTimeLi4',
                              'output/phase_space_studies_matter.root',
                              mass_li4, mass_he3, mass_p)
    
    RunTwoBodyDecaySimulation('output/spectra.root', 'h_blast_wave_antimatter',
                              'output/lifetime_li4.root', 'hDecayTimeLi4',
                              'output/phase_space_studies_antimatter.root',
                              mass_li4, mass_he3, mass_p)


