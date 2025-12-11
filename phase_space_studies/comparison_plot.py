from ROOT import TFile, TCanvas, TH1F, TLegend
from torchic.core.histogram import load_hist
from torchic.utils.root import set_root_object

def compare_toy_data(h_toy:TH1F, h_data:TH1F, outpdf:str):

    canvas = TCanvas()

    h_toy.Scale(1./h_toy.Integral())
    set_root_object(h_toy, line_color=797, line_width=2)
    h_data.Scale(1./h_data.Integral())
    set_root_object(h_data, line_color=418, line_width=2)

    legend = TLegend(0.6, 0.4, 0.8, 0.6)
    legend.SetFillColor(0)
    legend.SetBorderSize(0)
    legend.AddEntry(h_data, 'data', 'l')
    legend.AddEntry(h_toy, 'mc', 'l')

    ymax = max(h_toy.GetMaximum(), h_data.GetMaximum())
    canvas.DrawFrame(0., 0.0001, 10., ymax * 1.1, ';#it{p}_{T} (GeV/#it{c}); Normalised counts')
    canvas.SetLogy()
    h_toy.Draw('hist same')
    h_data.Draw('hist same')
    legend.Draw()

    canvas.SaveAs(outpdf)


if __name__ == '__main__':

    h_datas = [load_hist('/home/galucia/Lithium4/preparation/checks/same_event_hadronpid_pass1_pass4_nohe3pcut.root', 'QA/hPtHe'),
               load_hist('/home/galucia/Lithium4/preparation/checks/same_event_hadronpid_pass1_pass4_nohe3pcut.root', 'QA/hPtPr'),
               load_hist('/home/galucia/Lithium4/preparation/checks/same_event_hadronpid_pass1_pass4_nohe3pcut.root', 'QA/hPtLi')
              ]
    h_toys = [load_hist('output/phase_space_studies_matter.root', 'hPtFirstDaughter'),
              load_hist('output/phase_space_studies_matter.root', 'hPtSecondDaughter'),
              load_hist('output/phase_space_studies_matter.root', 'hPtMother')
              ]
    outpdfs = ['output/pt_he_comparison.pdf',
               'output/pt_pr_comparison.pdf',
               'output/pt_li_comparison.pdf']

    for h_data, h_toy, outpdf in zip(h_datas, h_toys, outpdfs):
        compare_toy_data(h_toy, h_data, outpdf)