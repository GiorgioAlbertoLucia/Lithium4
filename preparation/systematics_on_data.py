'''
    Script to load the data points, load the systematics, then create a plot showing the data with statistical uncertainties and systematic uncertainties separately
    and create a histogram to save in which the two are summed in quadrature.
'''
import numpy as np

from torchic.core.histogram import load_hist
from torchic.utils.root import set_root_object

from ROOT import TGraphMultiErrors, TCanvas, TLegend, gStyle, TFile, TDirectory, TH1F, TLatex, TGraph, TLine


def add_systematics_on_data(h_systematics:TH1F, outfile:TDirectory, sign:str, centrality:str):

    h_correlation_function = load_hist('/home/galucia/Lithium4/preparation/checks/correlation_hadronpid_pass1_pass4_refined_dca.root',
                                       f'Correlation{sign}/Default/hCorrelation{centrality if centrality != "1050" else "DirectComputation1050"}')
    
    h_correlation_stat_and_syst = h_correlation_function.Clone(h_correlation_function.GetName() + 'StatAndSyst')
    h_correlation_syst = h_correlation_function.Clone(h_correlation_function.GetName() + 'Syst')

    xs, ys, ex_low, ex_high, ey_stat_low, ey_stat_high, ey_syst_low, ey_syst_high = np.array([]), np.array([]), np.array([]), np.array([]), np.array([]), np.array([]), np.array([]), np.array([])

    for ibin in range(1, h_correlation_function.GetNbinsX() + 1):
        
        x = h_correlation_function.GetBinCenter(ibin)
        if x > 0.4:
            continue
        
        stat_error = h_correlation_function.GetBinError(ibin)
        syst_error = h_systematics.GetBinContent(ibin)

        xs = np.append(xs, x)
        ys = np.append(ys, h_correlation_function.GetBinContent(ibin))
        ex_low = np.append(ex_low, h_correlation_function.GetBinCenter(ibin) - h_correlation_function.GetBinLowEdge(ibin))
        ex_high = np.append(ex_high, h_correlation_function.GetBinLowEdge(ibin + 1) - h_correlation_function.GetBinCenter(ibin))
        ey_stat_low = np.append(ey_stat_low, stat_error)
        ey_stat_high = np.append(ey_stat_high, stat_error)
        ey_syst_low = np.append(ey_syst_low, syst_error)
        ey_syst_high = np.append(ey_syst_high, syst_error)

        total_error = (stat_error**2 + syst_error**2)**0.5
        h_correlation_stat_and_syst.SetBinError(ibin, total_error)
        h_correlation_syst.SetBinError(ibin, syst_error)
    
    g_correlation_stat_and_syst = TGraphMultiErrors(h_correlation_function.GetName().replace('h', 'g') + 'StatAndSyst',
                                                    f';{h_correlation_function.GetXaxis().GetTitle()};{h_correlation_function.GetYaxis().GetTitle()}',
                                                    len(xs), xs, ys, ex_low, ex_high, ey_stat_low, ey_stat_high)
    g_correlation_stat_and_syst.AddYError(len(xs), ey_syst_low, ey_syst_high)
    print(f'{ey_syst_low=}, {ey_syst_high=}')


    set_root_object(g_correlation_stat_and_syst, marker_style=20, marker_color=1, marker_size=1.2,
                     line_color=1, line_width=2)
    g_correlation_stat_and_syst.GetAttLine(0).SetLineColor(1)
    g_correlation_stat_and_syst.GetAttLine(0).SetLineWidth(2)
    g_correlation_stat_and_syst.GetAttLine(1).SetLineColor(1)
    g_correlation_stat_and_syst.GetAttLine(1).SetLineWidth(2)
    #g_correlation_stat_and_syst.GetAttFill(1).SetFillStyle(0)
    g_correlation_stat_and_syst.GetAttFill(1).SetFillColorAlpha(1, 0.3)

    canvas = TCanvas(f'cCorrelationFunctionWithSystematics{centrality}', 
                     f';{h_correlation_function.GetXaxis().GetTitle()};#it{{C}}(#it{{k}}*)', 800, 600)
    canvas.SetBottomMargin(0.13)
    canvas.SetLeftMargin(0.15)
    
    
    #canvas.SetLeftMargin(0.15)
    xmin = 0.0
    if centrality == '010' and sign == 'Antimatter':
        xmin = 0.02
    g_correlation_stat_and_syst.GetXaxis().SetRangeUser(xmin, 0.4)

    line = TLine(xmin, 1, 0.4, 1)
    set_root_object(line, line_color=15, line_style=2, line_width=2)
    
    hframe = canvas.DrawFrame(xmin, g_correlation_stat_and_syst.GetYaxis().GetXmin(), 
                              0.4, g_correlation_stat_and_syst.GetYaxis().GetXmax())
    hframe.SetTitle(';#it{k}* (GeV/#it{c});#it{C}(#it{k}*)')
    hframe.GetYaxis().SetTitleSize(0.05)
    hframe.GetXaxis().SetTitleSize(0.05)
    hframe.GetYaxis().SetLabelSize(0.045)
    hframe.GetXaxis().SetLabelSize(0.045)
    line.Draw('same')
    g_correlation_stat_and_syst.Draw('PS ; Z ; 2')

    latex = TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.045)
    latex.SetTextFont(42)
    latex.DrawLatex(0.22, 0.52, f'ALICE Preliminary')
    centrality_label = f'{centrality[:1]}-{centrality[1:]}' if centrality != '1050' else '10-50'
    sign_label = 'p#minus^{3}He' if sign == 'Matter' else '#bar{p}#minus^{3}#bar{He}'
    latex.DrawLatex(0.22, 0.46, f'Pb#minusPb #sqrt{{#it{{s}}_{{NN}}}} = 5.36 TeV')
    latex.DrawLatex(0.22, 0.40, f'FT0C Centrality: {centrality_label}%')
    latex.DrawLatex(0.22, 0.34, f'{sign_label}')

    dummy_stat = TGraph()
    dummy_stat.SetLineColor(1)
    dummy_stat.SetLineWidth(2)

    dummy_syst = TGraph()
    dummy_syst.SetFillColorAlpha(1, 0.3)
    dummy_syst.SetLineColor(0)

    legend = TLegend(0.6, 0.34, 0.75, 0.54)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.045)
    entry_data = legend.AddEntry(g_correlation_stat_and_syst, 'Data', 'p')
    entry_stat = legend.AddEntry(dummy_stat, 'Stat. uncert.', 'l')
    entry_syst = legend.AddEntry(dummy_syst, 'Syst. uncert.', 'f')
    legend.AddEntry(line, '#it{C}(#it{k}*) = 1', 'l')

    legend.Draw()
    canvas.SaveAs(f'/home/galucia/Lithium4/figures/systematics/correlation_with_systematics_{sign}_{centrality}.pdf')  

    canvas_comparison = TCanvas(f'cComparison{centrality}', 
                                     f';{h_correlation_function.GetXaxis().GetTitle()};{h_correlation_function.GetYaxis().GetTitle()}', 800, 600)
    set_root_object(h_correlation_stat_and_syst, marker_style=20, marker_color=797, line_color=797, line_width=3)
    set_root_object(h_correlation_function, marker_style=20, marker_color=797, line_color=797, line_width=3)
    h_correlation_stat_and_syst.Draw('p e0 same')
    h_correlation_function.Draw('p e0 same')

    outfile.cd()
    h_correlation_syst.Write()
    h_correlation_stat_and_syst.Write()
    g_correlation_stat_and_syst.Write()
    canvas.Write()
    canvas_comparison.Write()


if __name__ == "__main__":

    gStyle.SetOptStat(0)
    gStyle.SetPadTickX(1)
    gStyle.SetPadTickY(1)

    outfile  = TFile.Open('/home/galucia/Lithium4/preparation/output/correlation_with_systematics.root', 'RECREATE')

    for sign in ['Matter', 'Antimatter']:

        outdir = outfile.mkdir(sign)
        
        h_systematics = load_hist('/home/galucia/Lithium4/preparation/output/systematics_barlow_all_variables.root',
                                f'{sign}/hStd')
        
        for centrality in ['010', '1030', '3050', '050', '1050']:

            add_systematics_on_data(h_systematics, outdir, sign, centrality)

    outfile.Close()
