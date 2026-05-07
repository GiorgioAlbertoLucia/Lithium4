from ROOT import TRatioPlot, TCanvas, gStyle, TLegend, TLine, TLatex
from torchic.core.histogram import load_hist
from torchic.utils.root import set_root_object

def perform_ratio_plot(h_numerator, h_denominator, output_path,
                       **kwargs):
    gStyle.SetOptStat(0)

    canvas = TCanvas('canvas', 'canvas', 800, 800)
    rp = TRatioPlot(h_numerator, h_denominator, 'divsym')

    rp.SetH1DrawOpt(kwargs.get('hist_draw_opt', 'p e0 same'))
    rp.SetH2DrawOpt(kwargs.get('hist_draw_opt', 'p e0 same'))

    rp.Draw('nogrid same')
    ratio_title = kwargs.get('ratio_title', '')
    if ratio_title:
        rp.GetLowerRefYaxis().SetTitle(ratio_title.split(';')[-1])
        rp.GetLowerRefXaxis().SetTitle(ratio_title.split(';')[1])
    
    x_max = kwargs.get('x_max', h_numerator.GetXaxis().GetXmax())
    rp.GetUpperRefXaxis().SetRangeUser(kwargs.get('x_min', 0.), x_max)
    rp.GetLowerRefXaxis().SetRangeUser(kwargs.get('x_min', 0.), x_max)

    rp.GetUpperRefYaxis().SetRangeUser(kwargs.get('y_min', 0.), kwargs.get('y_max', 1.4))
    reference_line = TLine(rp.GetLowerPad().GetUxmin(), kwargs.get('line_y', 1), x_max, kwargs.get('line_y',1))
    set_root_object(reference_line, line_color=15, line_style=2)
    rp.GetLowerPad().cd()
    reference_line.Draw('same')

    rp.GetUpperPad().SetRightMargin(kwargs.get('right_margin', 0.05))
    rp.GetLowerPad().SetRightMargin(kwargs.get('right_margin', 0.05))
    rp.GetUpperPad().SetLeftMargin(kwargs.get('left_margin', 0.15))
    rp.GetLowerPad().SetLeftMargin(kwargs.get('left_margin', 0.15))

    rp.GetUpperPad().cd()
    legend = TLegend(*kwargs.get('legend_position', (0.5, 0.3, 0.8, 0.5)))
    legend.SetBorderSize(0)
    legend.AddEntry(h_numerator, h_numerator.GetName(), 'p')
    legend.AddEntry(h_denominator, h_denominator.GetName(), 'p')
    legend.Draw()

    canvas.SaveAs(output_path)

if __name__ == '__main__':

    gStyle.SetOptStat(0)

    h_same_010_matter = load_hist('/home/galucia/Lithium4/preparation/checks/correlation_hadronpid_pass1_pass4_nohe3pcut_offlinetpc.root', 
                            'CorrelationMatter/Default/hSameEvent010')
    h_same_010_antimatter = load_hist('/home/galucia/Lithium4/preparation/checks/correlation_hadronpid_pass1_pass4_nohe3pcut_offlinetpc.root', 
                            'CorrelationAntimatter/Default/hSameEvent010')
    
    set_root_object(h_same_010_matter, marker_color=797, line_color=797, name='Matter Same Event')
    set_root_object(h_same_010_antimatter, marker_color=418, line_color=418, name='Matter Same Event')
    
    perform_ratio_plot(h_same_010_matter, h_same_010_antimatter, 'check_mixed/ratio_same_010.pdf',
                       hist_draw_opt='hist e1 same', y_max=2e4, line_y=1.0, ratio_title='; k* (GeV/c); Matter / Antimatter',
                       legend_position=(0.5, 0.25, 0.8, 0.45),
                       right_margin=0.03, left_margin=0.25)
    
    h_mixed_010_matter = load_hist('/home/galucia/Lithium4/preparation/checks/correlation_hadronpid_pass1_pass4_nohe3pcut_offlinetpc.root', 
                            'CorrelationMatter/Default/hMixedEvent010')
    h_mixed_010_antimatter = load_hist('/home/galucia/Lithium4/preparation/checks/correlation_hadronpid_pass1_pass4_nohe3pcut_offlinetpc.root', 
                            'CorrelationAntimatter/Default/hMixedEvent010')
    
    set_root_object(h_mixed_010_matter, marker_color=797, line_color=797, name='Matter Mixed Event')
    set_root_object(h_mixed_010_antimatter, marker_color=418, line_color=418, name='Antimatter Mixed Event')

    perform_ratio_plot(h_mixed_010_matter, h_mixed_010_antimatter, 'check_mixed/ratio_mixed_010.pdf',
                       hist_draw_opt='hist e1 same', y_max=3.5e5, line_y=1.0, ratio_title='; k* (GeV/c); Matter / Antimatter',
                       right_margin=0.03, left_margin=0.25)
    

    # Correlation with TOF at pt > 0.75 GeV/c
    h_same_010_default = load_hist('/home/galucia/Lithium4/preparation/checks/correlation_hadronpid_pass1_pass4_refined_dca.root', 
                            'CorrelationAntimatter/Default/hCorrelation010')
    h_same_010_low_tof = load_hist('/home/galucia/Lithium4/preparation/checks/correlation_hadronpid_pass1_pass4_refined_dca_lower_pt_tof.root', 
                            'CorrelationAntimatter/Default/hCorrelation010')
    
    set_root_object(h_same_010_low_tof, marker_color=797, marker_style=20, line_color=797, name='TOF for #it{p}_{T} > 0.75 GeV/#it{c}', title=';#it{k}* (GeV/#it{c}); C(#it{k}*)')
    set_root_object(h_same_010_default, marker_color=418, marker_style=20, line_color=418, name='TOF for #it{p}_{T} > 0.80 GeV/#it{c}', title=';#it{k}* (GeV/#it{c}); C(#it{k}*)')
    
    perform_ratio_plot(h_same_010_low_tof, h_same_010_default, 'check_mixed/ratio_Ck_tof_010.pdf',
                       hist_draw_opt='p e1 same', y_max=1.3, x_max=0.4, line_y=1.0, ratio_title='; k* (GeV/c); 0.75 / 0.80',
                       legend_position=(0.5, 0.25, 0.8, 0.45),
                       right_margin=0.03, left_margin=0.25)
    
    h_same_050_default = load_hist('/home/galucia/Lithium4/preparation/checks/correlation_hadronpid_pass1_pass4_refined_dca.root', 
                            'CorrelationAntimatter/Default/hCorrelation050')
    h_same_050_low_tof = load_hist('/home/galucia/Lithium4/preparation/checks/correlation_hadronpid_pass1_pass4_refined_dca_lower_pt_tof.root', 
                            'CorrelationAntimatter/Default/hCorrelation050')
    
    set_root_object(h_same_050_low_tof, marker_color=797, marker_style=20, line_color=797, name='TOF for #it{p}_{T} > 0.75 GeV/#it{c}', title=';#it{k}* (GeV/#it{c}); C(#it{k}*)')
    set_root_object(h_same_050_default, marker_color=418, marker_style=20, line_color=418, name='TOF for #it{p}_{T} > 0.80 GeV/#it{c}', title=';#it{k}* (GeV/#it{c}); C(#it{k}*)')
    
    perform_ratio_plot(h_same_050_low_tof, h_same_050_default, 'check_mixed/ratio_Ck_tof_050.pdf',
                       hist_draw_opt='p e1 same', y_max=1.3, x_max=0.4, line_y=1.0, ratio_title='; k* (GeV/c); 0.75 / 0.80',
                       legend_position=(0.5, 0.25, 0.8, 0.45),
                       right_margin=0.03, left_margin=0.25)
    

    
