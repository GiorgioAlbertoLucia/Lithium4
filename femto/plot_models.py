'''
    Plot individual lambda contributions and their sum for antimatter, 010 centrality.
'''

import numpy as np
from ROOT import TFile, TCanvas, TLegend, TH1F, TPaveText, \
                 kRed, kAzure, kGreen, kBlack

from torchic.core.histogram import load_hist

input_Ck_path = '/home/galucia/Lithium4/femto/models/CATS_square_well_CF_LS_plot.root'
input_Sigma_Ck_path = '/home/galucia/Lithium4/femto/models/he3_Sigma_plus_Coulomb_plot.root'

def convert_MeV_to_GeV(hist):
    """Clone histogram and rescale x-axis from MeV to GeV."""
    nbins = hist.GetNbinsX()
    xlow  = hist.GetXaxis().GetXmin() / 1000.
    xhigh = hist.GetXaxis().GetXmax() / 1000.
    h_out = TH1F(hist.GetName() + '_GeV', hist.GetTitle(), nbins, xlow, xhigh)
    for ibin in range(1, nbins + 1):
        h_out.SetBinContent(ibin, hist.GetBinContent(ibin))
        h_out.SetBinError(ibin, hist.GetBinError(ibin))
    return h_out

if __name__ == '__main__':

    sign       = 'Antimatter'
    centrality = '010'

    h_lambda_param       = load_hist('/home/galucia/Lithium4/calibration/output/lambda_parameters_smaller_tolerance_pr.root',
                                     f'{sign}/hLambdaParameters')
    h_lambda_Sigma_param = load_hist('/home/galucia/Lithium4/calibration/output/lambda_parameters_smaller_tolerance_pr.root',
                                     f'{sign}/hLambdaSigmaParameters')

    h_genuine_raw = load_hist(input_Ck_path,       'r=6.39_fm/hHe3_p_Coul_CF')
    h_sigma_raw   = load_hist(input_Sigma_Ck_path, 'r=6.39_fm/hhe3_Sigma_plus_CF')

    h_genuine = convert_MeV_to_GeV(h_genuine_raw)
    h_sigma   = convert_MeV_to_GeV(h_sigma_raw)

    nbins = h_genuine.GetNbinsX()

    h_genuine_weighted  = h_genuine.Clone('hGenuineWeighted')
    h_sigma_weighted    = h_sigma.Clone('hSigmaWeighted')
    h_baseline          = h_genuine.Clone('hBaseline')
    h_sum               = h_genuine.Clone('hSum')

    for ibin in range(1, nbins + 1):
        kstar = h_genuine.GetBinCenter(ibin)

        lam       = h_lambda_param.GetBinContent(h_lambda_param.FindBin(kstar))
        lam_sigma = h_lambda_Sigma_param.GetBinContent(h_lambda_Sigma_param.FindBin(kstar))

        ck_genuine = h_genuine.GetBinContent(ibin)
        ck_sigma   = h_sigma.GetBinContent(h_sigma.FindBin(kstar))

        genuine_contrib  = lam       * ck_genuine
        sigma_contrib    = lam_sigma * ck_sigma
        baseline_contrib = 1.0 - lam - lam_sigma

        h_genuine_weighted.SetBinContent(ibin, genuine_contrib)
        h_sigma_weighted.SetBinContent(ibin,   sigma_contrib)
        h_baseline.SetBinContent(ibin,          baseline_contrib)
        h_sum.SetBinContent(ibin, genuine_contrib + sigma_contrib + baseline_contrib)

    # --- styling ---
    h_genuine_weighted.SetLineColor(kAzure  + 1);  h_genuine_weighted.SetLineWidth(3)
    h_sigma_weighted.SetLineColor(kGreen    + 2);  h_sigma_weighted.SetLineWidth(3)
    h_baseline.SetLineColor(797                );  h_baseline.SetLineWidth(3)
    h_sum.SetLineColor(kBlack);                    h_sum.SetLineWidth(3)

    # --- draw ---
    canvas = TCanvas('cContributions', 'Lambda Contributions', 800, 600)
    canvas.SetLeftMargin(0.15)
    canvas.SetBottomMargin(0.15)

    hframe = canvas.DrawFrame(0.011, -0.05, 0.4, 1.3,
                              f'{sign} | centrality {centrality}; #it{{k}}* (GeV/c); C(#it{{k}}*)')
    hframe.GetXaxis().SetTitleSize(0.05)
    hframe.GetYaxis().SetTitleSize(0.05)
    hframe.GetXaxis().SetLabelSize(0.04)
    hframe.GetYaxis().SetLabelSize(0.04)

    h_genuine_weighted.Draw('HIST SAME')
    h_sigma_weighted.Draw('HIST SAME')
    h_baseline.Draw('HIST SAME')
    h_sum.Draw('HIST SAME')

    legend = TLegend(0.18, 0.75, 0.88, 0.88)
    legend.SetTextSize(0.038)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetNColumns(2)
    legend.AddEntry(h_genuine_weighted, '#lambda_{p-^{3}He} #it{C}_{genuine}',      'l')
    legend.AddEntry(h_sigma_weighted,   '#lambda_{p#leftarrow#Sigma-^{3}He} #it{C}_{#Sigma^{+}-^{3}He}', 'l')
    legend.AddEntry(h_baseline,         '1 #minus #lambda_{p-^{3}He} #minus #lambda_{p#leftarrow#Sigma^{+}-^{3}He}', 'l')
    legend.AddEntry(h_sum,              '#it{C}_{full model}',                                    'l')
    legend.Draw()

    text = TPaveText(0.15, 0.40, 0.42, 0.58, 'NDC')
    text.SetTextSize(0.04)
    text.SetBorderSize(0)
    text.SetFillStyle(0)
    text.SetTextAlign(12)
    text.SetTextFont(42)
    text.AddText('p - ^{3}He')
    text.AddText('#it{R}_{source} = 6.39 fm')
    #text.Draw()

    outfile = TFile.Open('models/lambda_contributions_antimatter_010.root', 'recreate')
    for h in [h_genuine_weighted, h_sigma_weighted, h_baseline, h_sum]:
        h.Write()
    canvas.Write()
    canvas.SaveAs('models/lambda_contributions_antimatter_010.pdf')
    outfile.Close()

    print("Done.")