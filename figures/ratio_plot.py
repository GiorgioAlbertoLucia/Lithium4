from ROOT import TRatioPlot, TCanvas, gStyle, TLegend
from torchic.core.histogram import load_hist
from torchic.utils.root import set_root_object

if __name__ == '__main__':

    gStyle.SetOptStat(0)


    h23 = load_hist('/home/galucia/Lithium4/preparation/checks/correlation_all_pass1_pass4_23.root',
                     'CorrelationAntimatter/hCorrelation050')
    set_root_object(h23, name='2023', marker_color=601, marker_style=20, line_color=601)

    h24 = load_hist('/home/galucia/Lithium4/preparation/checks/correlation_all_pass1_pass4_24.root',
                         'CorrelationAntimatter/hCorrelation050')
    set_root_object(h24, name='2024', marker_color=797, marker_style=24, line_color=797)

    canvas = TCanvas('canvas', 'canvas', 800, 800)
    rp = TRatioPlot(h23, h24, 'divsym')

    rp.SetH1DrawOpt('p e0 same')
    rp.SetH2DrawOpt('p e0 same')

    rp.GetLowerPad().SetTitle('; #it{k}* (GeV/#it{c});Ratio (23/24)')
    rp.Draw('nogrid same')
    rp.GetUpperRefYaxis().SetRangeUser(0., 1.4)

    rp.GetUpperPad().cd()
    legend = TLegend(0.5, 0.3, 0.8, 0.5)
    legend.SetBorderSize(0)
    legend.AddEntry(h23, '2023', 'p')
    legend.AddEntry(h24, '2024', 'p')
    legend.Draw()

    canvas.SaveAs('/home/galucia/Lithium4/preparation/checks/ratio_23_24_all.pdf')

    
