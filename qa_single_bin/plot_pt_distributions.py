import numpy as np

from ROOT import TFile, TCanvas, TPad, TLine

from torchic.core.histogram import load_hist
from torchic.utils.root import set_root_object, init_legend, set_alice_global_style
from torchic.utils.colors import get_color

def produce_ratio_plot(hist1, hist2, canvas, h_reference):
    
    h_ratio = hist1.Clone(f'{hist1.GetName()}_ratio')
    h_ratio.Divide(hist2)
    set_root_object(h_ratio, line_color=get_color(2), line_width=2)
    
    canvas.cd()
    pad_ratio = TPad(f'cRatio_{hist1.GetName()}', f'cRatio_{hist1.GetName()}', 0, 0, 1, 0.3)
    pad_ratio.SetBottomMargin(0.3)
    pad_ratio.Draw()
    
    pad_ratio.cd()
    hframe_ratio = pad_ratio.DrawFrame(h_reference.GetXaxis().GetXmin(), 0, h_reference.GetXaxis().GetXmax(), 2,
                                      f';#it{{p}}_{{T}} (GeV/#it{{c}});SE / ME')
    hframe_ratio.GetXaxis().SetLabelSize(0.08)
    hframe_ratio.GetYaxis().SetLabelSize(0.08)
    hframe_ratio.GetXaxis().SetTitleSize(0.1)
    hframe_ratio.GetYaxis().SetTitleSize(0.1)
    hframe_ratio.GetYaxis().SetNdivisions(5)
    hframe_ratio.GetYaxis().SetTitleOffset(0.2)
    
    line = TLine(h_reference.GetXaxis().GetXmin(), 1, h_reference.GetXaxis().GetXmax(), 1)
    set_root_object(line, line_color=1, line_style=2, line_width=2)
    
    line.Draw('same')
    h_ratio.Draw('hist same')
    
    return h_ratio, pad_ratio, (hframe_ratio, line)


if __name__ == '__main__':
    
    set_alice_global_style()

    infile_path_se = '/home/galucia/Lithium4/preparation/output/PbPb/LHC25_PbPb_pass1_hadronpid_same_qa_kstar_bin.root'
    infile_path_me = '/home/galucia/Lithium4/preparation/checks/PbPb/LHC25_PbPb_pass1_hadronpid_event_mixing.root'
    
    h_pt_he3_vs_kstar_se = load_hist(infile_path_se, 'QA_dedicated_bin/hKstarVsPtHe1030')
    h_pt_he3_vs_kstar_me = load_hist(infile_path_me, 'QA_dedicated_bin/hKstarVsPtHe1030')
    
    h_pt_pr_vs_kstar_se = load_hist(infile_path_se, 'QA_dedicated_bin/hKstarVsPtPr1030')
    h_pt_pr_vs_kstar_me = load_hist(infile_path_me, 'QA_dedicated_bin/hKstarVsPtPr1030')
    
    outfile = TFile('output/pt_distributions.root', 'RECREATE')
    
    kstar_low = np.linspace(0.14, 0.24, 6)
    kstar_high = np.linspace(0.16, 0.26, 6)
    
    outdir_he = outfile.mkdir('he3')
    for i, (kstar_l, kstar_h) in enumerate(zip(kstar_low, kstar_high)):
        
        h_pt_he3_se = h_pt_he3_vs_kstar_se.ProjectionX(f'hPtHe3_SE_{i}', h_pt_he3_vs_kstar_se.GetYaxis().FindBin(kstar_l), h_pt_he3_vs_kstar_se.GetYaxis().FindBin(kstar_h))
        h_pt_he3_me = h_pt_he3_vs_kstar_me.ProjectionX(f'hPtHe3_ME_{i}', h_pt_he3_vs_kstar_me.GetYaxis().FindBin(kstar_l), h_pt_he3_vs_kstar_me.GetYaxis().FindBin(kstar_h))
        
        h_pt_he3_se.Scale(1.0 / h_pt_he3_se.Integral())
        h_pt_he3_me.Scale(1.0 / h_pt_he3_me.Integral())
        
        set_root_object(h_pt_he3_se, line_color=get_color(0), line_width=2)
        set_root_object(h_pt_he3_me, line_color=get_color(1), line_width=2, title=f'')
        
        canvas = TCanvas(f'cPtHe3_{i}', f'cPtHe3_{i}', 800, 600)
        upper_pad = TPad(f'cPtHe3_{i}_upper', f'cPtHe3_{i}_upper', 0, 0.3, 1, 1)
        upper_pad.SetBottomMargin(0.)
        upper_pad.Draw()
        
        upper_pad.cd()
        hframe = upper_pad.DrawFrame(-10, 0, 0, max(h_pt_he3_se.GetMaximum(), h_pt_he3_me.GetMaximum())*1.1, 
                                  f'10% < FT0C Centrality < 30%: {kstar_l:.2f} < #it{{k}}* < {kstar_h:.2f} GeV/#it{{c}};'
                                  f'#it{{p}}_{{T}} (GeV/#it{{c}});Normalized counts')
        h_pt_he3_me.Draw('hist same')
        h_pt_he3_se.Draw('hist same')
        
        legend = init_legend(0.2, 0.7, 0.36, 0.86)
        legend.AddEntry(h_pt_he3_se, f'SE', 'l')
        legend.AddEntry(h_pt_he3_me, f'ME', 'l')
        legend.Draw()
        
        h_ratio, pad_ratio, __ = produce_ratio_plot(h_pt_he3_se, h_pt_he3_me, canvas, hframe)
        
        outdir_he.cd()
        canvas.Write()
        
    outdir_pr = outfile.mkdir('pr')
    for i, (kstar_l, kstar_h) in enumerate(zip(kstar_low, kstar_high)):
        
        h_pt_pr_se = h_pt_pr_vs_kstar_se.ProjectionX(f'hPtPr_SE_{i}', h_pt_pr_vs_kstar_se.GetYaxis().FindBin(kstar_l), h_pt_pr_vs_kstar_se.GetYaxis().FindBin(kstar_h))
        h_pt_pr_me = h_pt_pr_vs_kstar_me.ProjectionX(f'hPtPr_ME_{i}', h_pt_pr_vs_kstar_me.GetYaxis().FindBin(kstar_l), h_pt_pr_vs_kstar_me.GetYaxis().FindBin(kstar_h))
        
        h_pt_pr_se.Scale(1.0 / h_pt_pr_se.Integral())
        h_pt_pr_me.Scale(1.0 / h_pt_pr_me.Integral())
        
        set_root_object(h_pt_pr_se, line_color=get_color(0), line_width=2)
        set_root_object(h_pt_pr_me, line_color=get_color(1), line_width=2)
        
        canvas = TCanvas(f'cPtpr_{i}', f'cPtpr_{i}', 800, 600)
        upper_pad = TPad(f'cPtpr_{i}_upper', f'cPtpr_{i}_upper', 0, 0.3, 1, 1)
        upper_pad.SetBottomMargin(0.02)
        upper_pad.Draw()
        
        upper_pad.cd()
        hframe = upper_pad.DrawFrame(-4, 0, 0, max(h_pt_pr_se.GetMaximum(), h_pt_pr_me.GetMaximum())*1.1, 
                                  f'10% < FT0C Centrality < 30%: {kstar_l:.2f} < #it{{k}}* < {kstar_h:.2f} GeV/#it{{c}};'
                                  f'#it{{p}}_{{T}} (GeV/#it{{c}});Normalized counts')
        h_pt_pr_me.Draw('hist same')
        h_pt_pr_se.Draw('hist same')
        
        legend = init_legend(0.2, 0.7, 0.36, 0.86)
        legend.AddEntry(h_pt_pr_se, f'SE', 'l')
        legend.AddEntry(h_pt_pr_me, f'ME', 'l')
        legend.Draw()
        
        h_ratio, pad_ratio, __ = produce_ratio_plot(h_pt_pr_se, h_pt_pr_me, canvas, hframe)
        
        outdir_pr.cd()
        canvas.Write()

        
    