from dataclasses import dataclass

from ROOT import TCanvas, TPaveText, gStyle

from torchic.core.histogram import load_hist
from torchic.utils.root import set_root_object, set_alice_global_style, init_legend, set_alice_frame_style
from torchic.utils.colors import get_color

@dataclass
class nucleus:
    mass: float
    width: float
    l: int
    normalization: float

if __name__ == '__main__':

    set_alice_global_style()

    infile_path = '/home/galucia/Lithium4/femto/models/li4_contribution_proper_sill.root'
    hist_gs = load_hist(infile_path, 'hKstar_g.s.')
    hist_ex1 = load_hist(infile_path, 'hKstar_ex. 1')
    hist_ex2 = load_hist(infile_path, 'hKstar_ex. 2')
    hist_ex3 = load_hist(infile_path, 'hKstar_ex. 3')
    
    li4_states = {'g.s.': nucleus(mass=3.75073, width=0.006, l=2, normalization=7.8425),
                  'ex. 1': nucleus(mass=3.75105, width=0.00735, l=1, normalization=4.6964),
                  'ex. 2': nucleus(mass=3.75281, width=0.00935, l=0, normalization=1.5490),
                  'ex. 3': nucleus(mass=3.75358, width=0.01351, l=1, normalization=4.6255)}
    integral_gs = hist_gs.Integral() * li4_states['g.s.'].normalization
    
    legend = init_legend(0.52, 0.35, 0.85, 0.55, text_size=0.045)
    text = TPaveText(0.5, 0.65, 0.8, 0.8, 'NDC')
    text.SetFillStyle(0)
    text.SetBorderSize(0)
    text.SetTextSize(0.045)
    text.SetTextFont(42)
    text.AddText('ALICE Simulation')
    text.AddText('Pb-Pb #sqrt{#it{s}_{NN}} = 5.36 TeV')
    text.AddText('^{4}Li #rightarrow p + ^{3}He')
    
    for ihist, (hist, hist_name, nucleus_state) in enumerate(zip([hist_gs, hist_ex1, hist_ex2, hist_ex3], 
                                                  ['^{4}Li g.s.', '^{4}Li ex. s. 1', '^{4}Li ex. s. 2', '^{4}Li ex. s. 3'],
                                                  li4_states.values())):
        set_root_object(hist, line_color=get_color(ihist), line_width=2,
                        name=f'hist_{ihist}', title=f'; #it{{k}}* (GeV/#it{{c}}); d#it{{N}}/d#it{{k}}* (GeV/#it{{c}})^{{-1}}')
        normalization = nucleus_state.normalization
        hist.Scale(normalization / integral_gs)
        legend.AddEntry(hist, hist_name, 'l')

    canvas = TCanvas('canvas', 'canvas', 700, 800)
    canvas.SetLeftMargin(0.2)
    canvas.SetRightMargin(0.05)
    canvas.SetTopMargin(0.05)
    canvas.SetBottomMargin(0.15)
    
    set_alice_frame_style(hist_gs)

    hist_gs.Draw('hist')
    hist_ex1.Draw('hist same')
    hist_ex2.Draw('hist same')
    hist_ex3.Draw('hist same')
    legend.Draw()
    text.Draw()
    canvas.SaveAs('output/li4_signal_contributions.pdf')