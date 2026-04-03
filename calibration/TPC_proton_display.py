'''
    Refit the TPC caibration
'''

from ROOT import TFile, TF1, TCanvas
from torchic.physics import BetheBloch

from torchic.core.histogram import load_hist

import sys
sys.path.append('..')
from utils.particles import ParticleMasses

def momentum_Bethe_Bloch(mass):

    def func(x, p):
        original = x[0]
        x[0] = original / mass
        val = BetheBloch(x, p)
        x[0] = original   # restore value
        return val

    return func


if __name__ == '__main__':

    th2_file_name = '/home/galucia/Lithium4/calibration/input/ProtonsFromLambda_23zzh.root'

    th2_name_protons = 'ebye-maker/QA/tpcSignalPr'
    h2_tpc_protons = load_hist(th2_file_name, th2_name_protons)

    momentum_Bethe_Bloch_protons = momentum_Bethe_Bloch(ParticleMasses['Pr'])
    f_bb_protons = TF1('f_bb_protons', momentum_Bethe_Bloch_protons, 0.05, 10, 5)
    f_bb_protons.SetParameters(-13.9261, -2.4156, 0.9469, 1.2628, 3.4574)
    f_bb_protons.SetLineColor(2)

    th2_name_all = 'ebye-maker/QA/tpcSignal'
    h2_tpc_all = load_hist(th2_file_name, th2_name_all)

    momentum_Bethe_Bloch_pions = momentum_Bethe_Bloch(ParticleMasses['Pi'])
    f_bb_pions = TF1('f_bb_protons', momentum_Bethe_Bloch_pions, 0.05, 10, 5)
    f_bb_pions.SetParameters(-13.9261, -2.4156, 0.9469, 1.2628, 3.4574)
    f_bb_pions.SetLineColor(3)

    outfile = TFile.Open('output/TPC_pr_display.root', 'RECREATE')
    
    canvas_protons = TCanvas('canvas_protons', 'TPC dE/dx vs p for protons', 800, 600)
    h2_tpc_protons.Draw('colz')
    f_bb_protons.Draw('same')

    canvas_all = TCanvas('canvas_all', 'TPC dE/dx vs p for all particles', 800, 600)
    h2_tpc_all.Draw('colz')
    f_bb_protons.Draw('same')
    f_bb_pions.Draw('same')

    outfile.cd()
    canvas_protons.Write()
    canvas_all.Write()

    outfile.Close()

