from ROOT import TFile, RooRealVar, RooDataHist, TCanvas

from torchic.core.histogram import load_hist
from torchic.roopdf import RooGausExp, RooGausDExp

if __name__ == "__main__":
    
    kstar = RooRealVar('kstar', '#it{k}^{*} (GeV/#it{c})', 0.05, 0.14)

    h_kstar_distribution = load_hist('models/li4_contribution_proper_sill.root', 'hKstar')
    #h_kstar_distribution.Fit('gaus', 'RMS+', '', 0.055, 0.11)

    kstar_data = RooDataHist('kstar_data', 'kstar_data', [kstar], Import=h_kstar_distribution)

    mean = RooRealVar('mean', 'mean', 0.081, 0.073, 0.1)
    sigma = RooRealVar('sigma', 'sigma', 0.01, 0.001, 0.1)
    tau = RooRealVar('tau', 'tau', 2, 0.1, 10)
    tau_2 = RooRealVar('tau_2', 'tau_2', -5, -10, -0.1)
    
    #pdf = RooGausExp('pdf', 'pdf', kstar, mean, sigma, tau)
    pdf = RooGausDExp('pdf', 'pdf', kstar, mean, sigma, tau_2, tau)
    
    pdf.fitTo(kstar_data)

    frame = kstar.frame()
    kstar_data.plotOn(frame)
    pdf.plotOn(frame)
    pdf.paramOn(frame)
    canvas = TCanvas('canvas', 'canvas', 800, 600)
    frame.Draw()

    outfile = TFile('output/li4_peak_fit.root', 'recreate')
    #h_kstar_distribution.Write()
    canvas.Write()
    outfile.Close()
