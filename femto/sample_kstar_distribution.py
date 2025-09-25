from ROOT import TFile, TCanvas, \
                kOrange, \
                RooDataSet, RooRealVar
from torchic.roopdf import RooSillPdf

def sampling(n_samples: int = 1_000_000):

    M_LI4 = 3.7513 # GeV/c^2
    W_LI4 = 0.005 # GeV/c^2 (li4 width)
    E_TH = 3.7466 # GeV/c^2

    invmass = RooRealVar('invmass', 'invmass', 3.473, 4.473, 'GeV/#it{c}^{2}')
    sill_params = {
        'mass': RooRealVar('mass', 'mass', M_LI4),
        'gamma': RooRealVar('mass', 'mass', W_LI4),
        'eth': RooRealVar('eth', 'eth', E_TH)
    }
    sill_pdf = RooSillPdf('sill', 'sill', invmass, *sill_params.values())

    frame = invmass.frame()
    sill_pdf.plotOn(frame, LineColor=kOrange-3)

    canvas = TCanvas('canvas', 'canvas', 800, 600)
    frame.Draw()

    outfile = TFile.Open('output/sampling.root', 'recreate')
    canvas.Write()
    outfile.Close()

if __name__ == '__main__':

    sampling()
