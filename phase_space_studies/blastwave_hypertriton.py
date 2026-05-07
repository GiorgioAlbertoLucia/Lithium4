import numpy as np
from ROOT import TH1F, TFile, gInterpreter

from uncertainties import ufloat

import sys
sys.path.append('..')
from utils.particles import ParticleMasses

gInterpreter.ProcessLine(f'#include "../include/ThermalModels.h"')
from ROOT import ThermalModels

# values taken form the Blast-Wave fit in https://arxiv.org/abs/2405.19839
T_KIN = {   
    '0_10': ufloat(0.103, 0.005), # GeV
    '10_30': ufloat(0.132, 0.008),
    '30_50': ufloat(0.152, 0.010),
}
N = {   
    '0_10': ufloat(0.498, 0.009),
    '10_30': ufloat(0.507, 0.012),
    '30_50': ufloat(0.660, 0.022),
}
BETA_T = {   
    '0_10': ufloat(0.694, 0.003),
    '10_30': ufloat(0.666, 0.003),
    '30_50': ufloat(0.598, 0.005),
}

def fit_blast_wave(mass:float, outfile:TFile, hist_to_fit:TH1F, name:str, centrality:str = '0_10'):

    T_kin = T_KIN[centrality].n
    T_kin_error = T_KIN[centrality].s
    beta_t = BETA_T[centrality].n
    beta_t_error = BETA_T[centrality].s
    n = N[centrality].n
    n_error = N[centrality].s

    beta_s = (2 + n) / 2 * beta_t
    beta_s_error = (2 + n) / 2 * beta_t_error

    max_bin = hist_to_fit.GetMaximumBin()
    normalisation = hist_to_fit.GetBinContent(max_bin) * 50
    
    thermal_model = ThermalModels()
    thermal_model.SetVarType(0) # (dN / dpt dy)

    blast_wave = thermal_model.GetBGBW(mass, beta_s, T_kin, n, normalisation)

    if 'Hypertriton' not in name:
        blast_wave.FixParameter(1, beta_s)
        blast_wave.FixParameter(2, T_kin)
        blast_wave.FixParameter(3, n)
    else:
        blast_wave.SetParLimits(1, beta_s - beta_s_error, beta_s + beta_s_error)
        blast_wave.SetParLimits(2, T_kin - T_kin_error, T_kin + T_kin_error)
        blast_wave.SetParLimits(3, n - n_error, n + n_error)
    
    #blast_wave.SetParameter(4, normalisation)

    blast_wave.SetTitle(f'Blast-Wave, m = {mass} GeV/#it{{c}}^{{2}}; #it{{p}}_{{T}} (GeV/#it{{c}}); d^{{2}}#it{{N}}/d#it{{p}}_{{T}}d#it{{y}}')
    fit_result = hist_to_fit.Fit(blast_wave, 'RMS+')

    outfile.cd()
    hist_to_fit.Write(f'h{name}')
    blast_wave.Write(f'f{name}')
    fit_result.Write(f'result{name}')

HYPERTRITON_CENTRALITY_DICT = {
    '0_10': '0-10% V0M centrality class',
    '10_30': '10-30% V0M centrality class',
    '30_50': '30-50% V0M centrality class',
}

def load_hypertriton_spectrum_hepdata(centrality:str = '0_10'):

    centrality_hypertriton = HYPERTRITON_CENTRALITY_DICT[centrality]
    infile_hypertriton = TFile.Open('input/HEPData-ins2791616-v1-root.root')

    h_hypertriton = infile_hypertriton.Get(f'(Anti)hypertriton spectrum in {centrality_hypertriton}/Hist1D_y1')
    h_hypertriton.SetDirectory(0)

    h_hypertriton_stat = infile_hypertriton.Get(f'(Anti)hypertriton spectrum in {centrality_hypertriton}/Hist1D_y1_e1')
    h_hypertriton_syst = infile_hypertriton.Get(f'(Anti)hypertriton spectrum in {centrality_hypertriton}/Hist1D_y1_e2')

    for ibin in range(1, h_hypertriton.GetNbinsX()+1):
        stat_error = h_hypertriton_stat.GetBinContent(ibin)
        syst_error = h_hypertriton_syst.GetBinContent(ibin)

        error = np.sqrt(stat_error*stat_error + syst_error*syst_error)
        h_hypertriton.SetBinError(ibin, error)

    return h_hypertriton    


def load_helium3_spectrum_hepdata(centrality:str = '0_10'):

    infile_helium3 = TFile.Open('input/He3_spectrum.root')

    h_helium3 = infile_helium3.Get(f'hHe3{centrality}')
    h_helium3.SetDirectory(0)

    h_helium3_stat = infile_helium3.Get(f'hHe3Stat{centrality}')
    h_helium3_syst = infile_helium3.Get(f'hHe3Syst{centrality}')

    for ibin in range(1, h_helium3.GetNbinsX()+1):
        stat_error = h_helium3_stat.GetBinContent(ibin)
        syst_error = h_helium3_syst.GetBinContent(ibin)

        error = np.sqrt(stat_error*stat_error + syst_error*syst_error)
        h_helium3.SetBinError(ibin, error)

    return h_helium3    


if __name__ == '__main__':

    outfile = TFile.Open('output/hypertriton_spectra_new_fit.root', 'recreate')
    
    mass_hypertriton = 2.991134 # GeV
    mass_he3 = ParticleMasses['He']

    for centrality in ['0_10', '10_30', '30_50']:
        outdir = outfile.mkdir(f'Centrality_{centrality}')
        h_hypertriton = load_hypertriton_spectrum_hepdata(centrality)
        h_helium3 = load_helium3_spectrum_hepdata(centrality)

        fit_blast_wave(mass_hypertriton, outdir, h_hypertriton, 'BlastWaveHypertriton')
        fit_blast_wave(mass_he3, outdir, h_helium3, 'BlastWaveHe3')
    
    outfile.Close()
