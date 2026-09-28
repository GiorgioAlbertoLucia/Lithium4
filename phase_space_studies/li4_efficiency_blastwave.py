import numpy as np
from ROOT import TF1, TH1F, TFile, TCanvas, gInterpreter
from torchic.physics.simulations import RunTwoBodyDecaySimulation
from torchic.core.histogram import load_hist

import sys
sys.path.append('..')
from utils.particles import ParticleMasses

gInterpreter.ProcessLine(f'#include "../include/ThermalModels.h"')
from ROOT import ThermalModels

INPUT_FILE_PERIODS = {
    '2023': '/home/galucia/Lithium4/preparation/output/PbPb/LHC25g11_efficiency.root',
    '2024': '/home/galucia/Lithium4/preparation/output/PbPb/LHC26e5_efficiency.root',
    '2025': '/home/galucia/Lithium4/preparation/output/PbPb/LHC26e6_efficiency.root'
}

def sample_blast_wave(mass:float, outfile:TFile, n_samples:int = 1_000_000):

    # values taken form the Blast-Wave fit in https://arXiv.org/abs/2311.11758
    T_kin = 0.132 # GeV
    beta_t = 0.670
    n = 0.55
    normalisation = 1.
    beta_s = (2 + n) / 2 * beta_t
    
    thermal_model = ThermalModels()
    blast_wave = thermal_model.GetBGBW(mass, beta_s, T_kin, n, normalisation)
    
    blast_wave.SetTitle(f'Blast-Wave, m = {mass} GeV/#it{{c}}^{{2}}; #it{{p}}_{{T}} (GeV/#it{{c}}); d^{{2}}#it{{N}}/d#it{{p}}_{{T}}d#it{{y}}')
    canvas = TCanvas('canavs', blast_wave.GetTitle())
    canvas.SetLogy()
    
    canvas.DrawFrame(0., 1e-16, 10., 1e-13, blast_wave.GetTitle())
    blast_wave.Draw('same')

    h_blast_wave = TH1F('h_blast_wave', 
                                      f'Blast-Wave, m = {mass} GeV/#it{{c}}^{{2}}; #it{{p}}_{{T}} (GeV/#it{{c}}); d^{{2}}#it{{N}}/d#it{{p}}_{{T}}d#it{{y}}',
                                      100, 0, 10)
    for isample in range(n_samples):
        h_blast_wave.Fill(blast_wave.GetRandom())

    outfile.cd()
    blast_wave.Write()
    canvas.Write()
    h_blast_wave.Write()
    canvas.SaveAs('output/li4_blast_wave.pdf')

    return h_blast_wave

def compute_weighted_efficiency(h_dNdpt:TH1F, outfile:TFile):

    for period in ['2023', '2024', '2025']:
    
        h_efficiency_matter = load_hist(INPUT_FILE_PERIODS[period], 'Efficiency/hEfficiencyMatter')
        h_efficiency_antimatter = load_hist(INPUT_FILE_PERIODS[period], 'Efficiency/hEfficiencyAntimatter')

        efficiency_matter, efficiency_antimatter = 0., 0.
        total_weight_matter, total_weight_antimatter = 0., 0.

        for ibin in range(1, h_efficiency_matter.GetNbinsX() + 1):
            pt = h_efficiency_matter.GetBinCenter(ibin)
            spectrum_bin_content = h_dNdpt.GetBinContent(h_dNdpt.FindBin(pt))

            efficiency_matter += h_efficiency_matter.GetBinContent(ibin) * spectrum_bin_content
            total_weight_matter += spectrum_bin_content

            efficiency_antimatter += h_efficiency_antimatter.GetBinContent(ibin) * spectrum_bin_content
            total_weight_antimatter += spectrum_bin_content

        efficiency_matter /= total_weight_matter
        efficiency_antimatter /= total_weight_antimatter

        print(f"Period: {period}")
        print(f"Weighted efficiency matter: {efficiency_matter}")
        print(f"Weighted efficiency antimatter: {efficiency_antimatter}\n")
        
        outfile.cd()
        h_dNdpt.Write()
        h_efficiency_matter.Write()
        h_efficiency_antimatter.Write()

def compute_weighted_efficiency_centrality(h_dNdpt:TH1F, outfile:TFile):

    for period in ['2023', '2024', '2025']:
        for centrality in ['010', '1050']:

            h_efficiency = load_hist(INPUT_FILE_PERIODS[period], f'Efficiency/hEfficiency{centrality}')

            efficiency_matter, efficiency_antimatter = 0., 0.
            efficiency_matter_error, efficiency_antimatter_error = 0., 0.
            total_weight_matter, total_weight_antimatter = 0., 0.

            for ibin in range(1, h_efficiency.GetNbinsX() + 1):
                pt = h_efficiency.GetBinCenter(ibin)
                spectrum_bin_content = h_dNdpt.GetBinContent(h_dNdpt.FindBin(abs(pt)))

                if pt >= 0:  # matter
                    efficiency_matter += h_efficiency.GetBinContent(ibin) * spectrum_bin_content
                    efficiency_matter_error += (h_efficiency.GetBinError(ibin) * spectrum_bin_content)**2
                    total_weight_matter += spectrum_bin_content
                else:        # antimatter
                    efficiency_antimatter += h_efficiency.GetBinContent(ibin) * spectrum_bin_content
                    efficiency_antimatter_error += (h_efficiency.GetBinError(ibin) * spectrum_bin_content)**2
                    total_weight_antimatter += spectrum_bin_content

            efficiency_matter /= total_weight_matter
            efficiency_antimatter /= total_weight_antimatter
            efficiency_matter_error = np.sqrt(efficiency_matter_error) / total_weight_matter
            efficiency_antimatter_error = np.sqrt(efficiency_antimatter_error) / total_weight_antimatter

            print(f"Period: {period}, Centrality: {centrality}")
            print(f"Weighted efficiency matter: {efficiency_matter} ± {efficiency_matter_error}")
            print(f"Weighted efficiency antimatter: {efficiency_antimatter} ± {efficiency_antimatter_error}\n")

            outfile.cd()
            h_efficiency.Write()



if __name__ == '__main__':

    outfile = TFile.Open('output/li4_efficiency_blastwave.root', 'recreate')
    
    mass_li4 = 3.751 # GeV/c^2
    mass_he3 = ParticleMasses['He']
    mass_p = ParticleMasses['Pr']
    
    h_blast_wave = sample_blast_wave(mass_li4, outfile, 10_000_000)
    compute_weighted_efficiency(h_blast_wave, outfile)
    compute_weighted_efficiency_centrality(h_blast_wave, outfile)
    outfile.Close()
