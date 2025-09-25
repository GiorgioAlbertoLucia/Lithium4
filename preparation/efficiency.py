import argparse
import yaml
import numpy as np
import ROOT
from ROOT import TFile, TChain, TList, TTree

from utils.particles import ParticleMasses

ROOT.gROOT.LoadMacro('Common.h++')
from ROOT import NsigmaTpcHe, NsigmaITSHe, NsigmaITSPr, NsigmaTOFPr, averageClusterSize, CorrectPidTrkHe, \
  Kstar, ExpectedClusterSizeCosLambdaHe, ExpectedClusterSizeCosLambdaPr

ROOT.EnableImplicitMT(30)
ROOT.gROOT.SetBatch(True)

base_selection = '(fSignedPtHe3 > 0 && fSignedPtHad > 0) || (fSignedPtHe3 < 0 && fSignedPtHad < 0)'
selections = [
  'true',
  '((fDeltaPhi/0.006)*(fDeltaPhi/0.006) + (fDeltaEta/0.012)*(fDeltaEta/0.012) > 1)',

  '(std::abs(fEtaHe3) < 0.9)',
  '((fPtHe3 < 2.5) || (fPIDtrkHe3 == 7))',
  '((0.5 < fChi2TPCHe3) && (fChi2TPCHe3 < 4))',
  '((-1.5 < fNSigmaTPCHe3) && (fNSigmaTPCHe3 < 2.5))',
  #'((-0.9 < fNSigmaTPCHe3) && (fNSigmaTPCHe3 < 1.9))',
  '(std::abs(fDCAxyHe3) < 0.1)',
  '(std::abs(fDCAzHe3) < 1.0)',
  '(fClusterSizeCosLamHe3 > 4)',
  #'(fNSigmaITSHe3 > -1.5)',

  '(std::abs(fEtaHad) < 0.9)',
  '((0.5 < fChi2TPCHad) && (fChi2TPCHad < 4))',
  '(std::abs(fNSigmaTPCHad) < 2)',
  '((std::abs(fPtHad) < 0.8) || (std::abs(fNSigmaTOFHad) < 2))',
  #'(fNSigmaITSHad > -1.5)',
]
#selections = conf['selections']
selection = selections[0]
for sel in selections[1:]:
    selection += (' && ' + sel)
print(f'Selection: {selection}')

############################################################################################################

input_data = ['input/mc.root',
              ]
file_data_list = input_data if isinstance(input_data, list) else [input_data]

chainData = TChain("O2he3hadtable")
tree_name = "O2he3hadtable"

mode = 'DF'  # DataFrame mode

for fileName in file_data_list:
  fileData = TFile(fileName)

  if mode == 'DF':
    for key in fileData.GetListOfKeys():
      keyName = key.GetName()
      if 'DF_' in keyName :
          print(f'Adding {fileName}/{keyName}/{tree_name} to the chain')
          chainData.Add(f'{fileName}/{keyName}/{tree_name}')
  elif mode == 'tree':
    print(f'Adding {fileName}/{tree_name} to the chain')
    chainData.Add(f'{fileName}/{tree_name}')

############################################################################################################################################################################
rdf_gen = ROOT.ROOT.RDataFrame(chainData)
rdf_rec = rdf_gen.Filter('fPtHe3 > -900') \
        .Define('fSignedPtHad', 'fPtHad') \
        .Define('fSignHe3', 'fPtHe3/std::abs(fPtHe3)') \
        .Redefine('fPtHe3', 'std::abs(fPtHe3)') \
        .Redefine('fPtHad', 'std::abs(fPtHad)') \
        .Redefine('fPtHe3', '(fPIDtrkHe3 == 7) || (fPtHe3 > 2.5) ? fPtHe3 : CorrectPidTrkHe(fPtHe3)') \
        .Define('fSignedPtHe3', 'fPtHe3 * fSignHe3') \
        .Define(f'fEHe3', f'std::sqrt((fPtHe3 * std::cosh(fEtaHe3))*(fPtHe3 * std::cosh(fEtaHe3)) + {ParticleMasses["He"]}*{ParticleMasses["He"]})') \
        .Define(f'fEHad', f'std::sqrt((fPtHad * std::cosh(fEtaHad))*(fPtHad * std::cosh(fEtaHad)) + {ParticleMasses["Pr"]}*{ParticleMasses["Pr"]})') \
        .Define('fDeltaEta', 'fEtaHe3 - fEtaHad') \
        .Define('fDeltaPhi', 'fPhiHe3 - fPhiHad') \
        .Redefine('fInnerParamTPCHe3', 'fInnerParamTPCHe3 * 2') \
        .Redefine('fNSigmaTPCHe3', 'NsigmaTpcHe(std::abs(fInnerParamTPCHe3), fSignalTPCHe3)') \
        .Define('fNSigmaTOFHad', 'NsigmaTOFPr(fPtHad * std::cosh(fEtaHad), fMassTOFHad)') \
        .Define('fClusterSizeCosLamHe3', 'averageClusterSize(fItsClusterSizeHe3) / cosh(fEtaHe3)') \
        .Define('fClusterSizeCosLamHad', 'averageClusterSize(fItsClusterSizeHad) / cosh(fEtaHad)') \
        .Define('fExpectedClusterSizeHe3', 'ExpectedClusterSizeCosLambdaHe(fPtHe3)') \
        .Define('fExpectedClusterSizeHad', 'ExpectedClusterSizeCosLambdaPr(fPtHad)') \
        .Define('fNSigmaITSHe3', 'NsigmaITSHe(fPtHe3 * std::cosh(fEtaHe3), fClusterSizeCosLamHe3)') \
        .Define('fNSigmaITSHad', 'NsigmaITSPr(fPtHad * std::cosh(fEtaHad), fClusterSizeCosLamHad)') \
        .Filter(base_selection).Filter(selection) \
        .Define('fPxLi', 'fPtHe3 * std::cos(fPhiHe3) + fPtHad * std::cos(fPhiHad)') \
        .Define('fPyLi', 'fPtHe3 * std::sin(fPhiHe3) + fPtHad * std::sin(fPhiHad)') \
        .Define('fPtLi', 'sqrt(fPxLi*fPxLi + fPyLi*fPyLi)') \

rdf_gen = rdf_gen.Define('fSignedPtHad', 'fPtMCHad') \
        .Define('fSignedPtHe3', 'fPtMCHe3') \
        .Filter(base_selection) \
        .Define('fSignHe3', 'fPtMCHe3/std::abs(fPtMCHe3)') \
        .Define('fPtLiMC', 'std::abs(fSignedPtMC)') \
        

pt_bins, pt_min, pt_max = 100, 0, 10

hPthe3  = rdf_rec.Histo1D(("hPthe3", ";#it{p}_{T} (^{3}He) (GeV/#it{c});", 2*pt_bins, -10, pt_max), "fSignedPtHe3").GetValue()
hPtHe3MC = rdf_gen.Histo1D(("hPtHe3MC", ";#it{p}_{T} (^{3}He) (GeV/#it{c});", 2*pt_bins, -10, pt_max), "fPtMCHe3").GetValue()

hPtLiMatter = rdf_rec.Filter('fSignHe3 > 0').Histo1D(("hPtLiMatter", ";#it{p}_{T} (^{4}Li) (GeV/#it{c});", pt_bins, pt_min, pt_max), "fPtLi").GetValue()
hPtLiAntimatter = rdf_rec.Filter('fSignHe3 < 0').Histo1D(("hPtLiAntimatter", ";#it{p}_{T} (^{4}Li) (GeV/#it{c});", pt_bins, pt_min, pt_max), "fPtLi").GetValue()
hPtLiMCMatter = rdf_gen.Filter('fSignHe3 > 0').Histo1D(("hPtLiMCMatter", ";#it{p}_{T} (^{4}#bar{Li}) (GeV/#it{c});", pt_bins, pt_min, pt_max), "fPtLiMC").GetValue()
hPtLiMCAntimatter = rdf_gen.Filter('fSignHe3 < 0').Histo1D(("hPtLiMCAntimatter", ";#it{p}_{T} (^{4}#bar{Li}) (GeV/#it{c});", pt_bins, pt_min, pt_max), "fPtLiMC").GetValue()


def efficiency_histogram(hist_rec, hist_gen, name):
    """
    Calculate the efficiency histogram by dividing the reconstructed histogram by the generated histogram.
    """
    efficiency_hist = hist_rec.Clone(name)

    for ibin in range(1, hist_rec.GetNbinsX() + 1):
        if hist_gen.GetBinContent(ibin) > 0:
            efficiency = hist_rec.GetBinContent(ibin) / hist_gen.GetBinContent(ibin)
            efficiency_error = np.sqrt( efficiency * (1 - efficiency) / hist_gen.GetBinContent(ibin))
            efficiency_hist.SetBinContent(ibin, efficiency)
            efficiency_hist.SetBinError(ibin, efficiency_error)
        else:
            efficiency_hist.SetBinContent(ibin, 0)
            efficiency_hist.SetBinError(ibin, 0)

    return efficiency_hist

hEfficiencyMatter = efficiency_histogram(hPtLiMatter, hPtLiMCMatter, "hEfficiencyMatter")
hEfficiencyAntimatter = efficiency_histogram(hPtLiAntimatter, hPtLiMCAntimatter, "hEfficiencyAntimatter")

def get_mean_efficiency(hist_efficiency):
   
    mean = 0
    for ibin in range(1, hist_efficiency.GetNbinsX() + 1):
        mean += hist_efficiency.GetBinContent(ibin)

    mean /= hist_efficiency.GetNbinsX()
    return mean

average_efficiency_matter = get_mean_efficiency(hEfficiencyMatter)
average_efficiency_antimatter = get_mean_efficiency(hEfficiencyAntimatter)
print(f"Average efficiency for matter: {average_efficiency_matter:.4f}")
print(f"Average efficiency for antimatter: {average_efficiency_antimatter:.4f}")

output_file = 'output/efficiency.root'
outFile = ROOT.TFile(output_file, "RECREATE")
hPthe3.Write()
hPtHe3MC.Write()
hPtLiMatter.Write()
hPtLiAntimatter.Write()
hPtLiMCMatter.Write()
hPtLiMCAntimatter.Write()
hEfficiencyMatter.Write()
hEfficiencyAntimatter.Write()
print(f"Efficiency histograms saved to {output_file}")

outFile.Close()


# draw efficiencies
ROOT.gStyle.SetOptStat(0)
canvas_efficiency = ROOT.TCanvas("canvas_efficiency", "Efficiency Histograms", 800, 600)
hframe = canvas_efficiency.DrawFrame(0, 0, 10.5, 0.15)
hEfficiencyMatter.SetMarkerColor(ROOT.kBlue+2)
hEfficiencyMatter.SetMarkerStyle(22)
hEfficiencyMatter.SetLineColor(ROOT.kBlue+2)
hEfficiencyAntimatter.SetMarkerColor(ROOT.kOrange-3)
hEfficiencyAntimatter.SetMarkerStyle(23)
hEfficiencyAntimatter.SetLineColor(ROOT.kOrange-3)


legend = ROOT.TLegend(0.6, 0.7, 0.75, 0.85)
legend.AddEntry(hEfficiencyMatter, "^{4}Li", "p")
legend.AddEntry(hEfficiencyAntimatter, "^{4}#bar{Li}", "p")
legend.SetBorderSize(0)
legend.SetFillStyle(0)

text = ROOT.TPaveText(0.15, 0.55, 0.4, 0.85, 'ndc')
text.AddText("This work")
text.AddText("#bf{ALICE Run 3}")
text.AddText("#bf{Pb-Pb, #sqrt{s_{NN}} = 5.36 TeV}")
text.AddText(f"#bf{{#LT #varepsilon(^{{4}}Li) #GT = {average_efficiency_matter:.3f}}}")
text.AddText(f"#bf{{#LT #varepsilon(^{{4}}#bar{{Li}}) #GT = {average_efficiency_antimatter:.3f}}}")
text.SetFillStyle(0)
text.SetBorderSize(0)

hframe.SetTitle(";#it{p}_{T} (GeV/#it{c}); Efficiency")
hEfficiencyAntimatter.Draw("E1 same")
hEfficiencyMatter.Draw("E1 same")
legend.Draw()
text.Draw('same')

canvas_efficiency.SaveAs("output/efficiency_histograms.pdf")

