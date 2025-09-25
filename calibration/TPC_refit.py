'''
    Refit the TPC caibration
'''

from ROOT import TFile, TF1
from torchic.physics import BetheBloch


if __name__ == '__main__':

    infile = TFile.Open('output/TPC_he_mc.root', 'UPDATE')
    graph = infile.Get('TPC_antimatter/g_mean_antimatter')
    graph.GetListOfFunctions().Clear()
    infile.Close()

    outfile = TFile.Open('output/TPC_he_mc_refit.root', 'RECREATE')
    f_mean = TF1('f_mean', '[0] / (x^[1]) + [2]', 0.8, 4.)
    f_mean.SetParameters(500, 1.5, 300)

    graph.Fit(f_mean, 'RMS+')
    graph.Write('g_mean_refit_3pFit')
    graph.GetListOfFunctions().Clear()

    f_mean = TF1('f_mean', BetheBloch, 0.8, 3.8, 5)
    f_mean.SetParameters(-241.4902, 0.374245, 1.397847, 1.0782504, 2.048336)
    graph.Fit(f_mean, 'RMS+')
    graph.Write('g_mean_refit_bethe_bloch')

    outfile.Close()

