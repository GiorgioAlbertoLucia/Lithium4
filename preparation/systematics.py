import argparse
import yaml
import numpy as np
import ROOT
from ROOT import TFile, TChain, TList, TTree
from alive_progress import alive_bar

from utils.particles import ParticleMasses

ROOT.gROOT.LoadMacro('Common.h++')
from ROOT import NsigmaTpcHe, NsigmaITSHe, NsigmaITSPr, NsigmaTOFPr, averageClusterSize, CorrectPidTrkHe, \
  Kstar, ExpectedClusterSizeCosLambdaHe, ExpectedClusterSizeCosLambdaPr

ROOT.EnableImplicitMT(30)
ROOT.gROOT.SetBatch(True)

selections = [
    '((fSignedPtHe3 > 0 && fSignedPtHad > 0) || (fSignedPtHe3 < 0 && fSignedPtHad < 0))',
    '((fDeltaPhi/0.006)*(fDeltaPhi/0.006) + (fDeltaEta/0.012)*(fDeltaEta/0.012) > 1)',

    '(std::abs(fEtaHe3) < 0.9)',
    '((fPtHe3 < 2.5) || (fPIDtrkHe3 == 7))',
    '((0.5 < fChi2TPCHe3) && (fChi2TPCHe3 < 4))',

    '(std::abs(fEtaHad) < 0.9)',
    '((0.5 < fChi2TPCHad) && (fChi2TPCHad < 4))',
]
#selections = conf['selections']
base_selection = selections[0]
for sel in selections[1:]:
    base_selection += (' && ' + sel)
print(f'Selection: {base_selection}')


def load_same():
    
    input_data = ['input/LHC23_PbPb_pass4_long_same_lsus_merged.root',
                  'input/LHC24ar_pass1_same_merged.root',
                  'input/LHC24as_pass1_same_merged.root']
    file_data_list = input_data if isinstance(input_data, list) else [input_data]

    chainData = TChain("O2he3hadtable")
    tree_name = "O2he3hadtable"

    for fileName in file_data_list:
        fileData = TFile(fileName)

        for key in fileData.GetListOfKeys():
          keyName = key.GetName()
          if 'DF_' in keyName :
              print(f'Adding {fileName}/{keyName}/{tree_name} to the chain')
              chainData.Add(f'{fileName}/{keyName}/{tree_name}')
      
    rdf = ROOT.ROOT.RDataFrame(chainData) \
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
          .Redefine('fNSigmaTPCHe3', 'NsigmaTpcHe(fInnerParamTPCHe3, fSignalTPCHe3)') \
          .Define('fNSigmaTOFHad', 'NsigmaTOFPr(fPtHad * std::cosh(fEtaHad), fMassTOFHad)') \
          .Define('fClusterSizeCosLamHe3', 'averageClusterSize(fItsClusterSizeHe3) / cosh(fEtaHe3)') \
          .Define('fClusterSizeCosLamHad', 'averageClusterSize(fItsClusterSizeHad) / cosh(fEtaHad)') \
          .Define('fExpectedClusterSizeHe3', 'ExpectedClusterSizeCosLambdaHe(fPtHe3)') \
          .Define('fExpectedClusterSizeHad', 'ExpectedClusterSizeCosLambdaPr(fPtHad)') \
          .Define('fNSigmaITSHe3', 'NsigmaITSHe(fPtHe3 * std::cosh(fEtaHe3), fClusterSizeCosLamHe3)') \
          .Define('fNSigmaITSHad', 'NsigmaITSPr(fPtHad * std::cosh(fEtaHad), fClusterSizeCosLamHad)') \
          .Filter(base_selection) \
          .Define('fKstar', f'Kstar(fPtHe3, fEtaHe3, fPhiHe3, {ParticleMasses["He"]}, fPtHad, fEtaHad, fPhiHad, {ParticleMasses["Pr"]})') \
          
    return rdf, chainData

def load_mixed():

    input_data = ['/data/galucia/lithium_local/mixing/LHC23_PbPb_pass4_long_mixing_lsus.root',
                  '/data/galucia/lithium_local/mixing/LHC24ar_pass1_mixing_lsus.root',
                  '/data/galucia/lithium_local/mixing/LHC24as_pass1_mixing_lsus.root',]
    file_data_list = input_data if isinstance(input_data, list) else [input_data]

    chainData = TChain("O2he3hadtable")
    tree_name = "MixedTree"

    for fileName in file_data_list:
        print(f'Adding {fileName}/{tree_name} to the chain')
        chainData.Add(f'{fileName}/{tree_name}')

    rdf = ROOT.ROOT.RDataFrame(chainData) \
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
          .Redefine('fNSigmaTPCHe3', 'NsigmaTpcHe(fInnerParamTPCHe3, fSignalTPCHe3)') \
          .Define('fNSigmaTOFHad', 'NsigmaTOFPr(fPtHad * std::cosh(fEtaHad), fMassTOFHad)') \
          .Define('fClusterSizeCosLamHe3', 'averageClusterSize(fItsClusterSizeHe3) / cosh(fEtaHe3)') \
          .Define('fClusterSizeCosLamHad', 'averageClusterSize(fItsClusterSizeHad) / cosh(fEtaHad)') \
          .Define('fExpectedClusterSizeHe3', 'ExpectedClusterSizeCosLambdaHe(fPtHe3)') \
          .Define('fExpectedClusterSizeHad', 'ExpectedClusterSizeCosLambdaPr(fPtHad)') \
          .Define('fNSigmaITSHe3', 'NsigmaITSHe(fPtHe3 * std::cosh(fEtaHe3), fClusterSizeCosLamHe3)') \
          .Define('fNSigmaITSHad', 'NsigmaITSPr(fPtHad * std::cosh(fEtaHad), fClusterSizeCosLamHad)') \
          .Filter(base_selection) \
          .Filter('(fCentralityFT0C < 50)') \
          .Define('fKstar', f'Kstar(fPtHe3, fEtaHe3, fPhiHe3, {ParticleMasses["He"]}, fPtHad, fEtaHad, fPhiHad, {ParticleMasses["Pr"]})') \
          
    return rdf, chainData
          
def run_systematics():

    rdf_same, _chain_same = load_same()
    rdf_mixed, _chain_mixed = load_mixed()

    outFile = TFile("output/systematics.root", "RECREATE")
    h_corr_iters = []

    N_ITERATIONS = 200
    with alive_bar(N_ITERATIONS, title='Running systematics...') as bar:
        for iter in range(N_ITERATIONS):
            
            dcaxy_he3_max = np.random.uniform(0.05, 0.15)
            dcaz_he3_max = np.random.uniform(0.75, 1.0)
            n_sigma_tpc_he3_high = np.random.uniform(2., 3.)
            n_sigma_tpc_he3_low = np.random.uniform(-2.0, -1.0)
            n_sigma_tpc_had_max = np.random.uniform(1.5, 2.5)
            n_sigma_tof_had_max = np.random.uniform(1.5, 2.5)

            condition = f'((std::abs(fDCAxyHe3) < {dcaxy_he3_max}) && (std::abs(fDCAzHe3) < {dcaz_he3_max}) && ' \
                        f'((fNSigmaTPCHe3 < {n_sigma_tpc_he3_high}) && (fNSigmaTPCHe3 > {n_sigma_tpc_he3_low})) && ' \
                        f'((fNSigmaTPCHe3 < {n_sigma_tpc_had_max}) && (fNSigmaTPCHe3 > -{n_sigma_tpc_had_max})) && ' \
                        f'((fNSigmaTOFHad < {n_sigma_tof_had_max}) && (fNSigmaTOFHad > -{n_sigma_tof_had_max})) '

            h_same_iter = rdf_same.Filter(condition).Histo1D((f"hKstarSameIter{iter}", ";#it{k}^{*} (GeV/#it{c});", 100, 0, 1.), "fKstar").GetValue()
            h_mixed_iter = rdf_mixed.Filter(condition).Histo1D((f"hKstarMixedIter{iter}", ";#it{k}^{*} (GeV/#it{c});", 100, 0, 1.), "fKstar").GetValue()
            
            low_bin = h_same_iter.FindBin(0.25)
            high_bin = h_same_iter.FindBin(1)
            normalization_factor = h_same_iter.Integral(low_bin, high_bin) / h_mixed_iter.Integral(low_bin, high_bin)
            h_mixed_iter.Scale(normalization_factor)

            h_corr_iter = h_same_iter.Clone(f'hCorrelationIter{iter}')
            h_corr_iter.Divide(h_mixed_iter)

            h_corr_iters.append(h_corr_iter)

            bar()
    
    point_positions = {}
    outFile.cd()
    for ibin in range(1, h_corr_iters[0].GetNbinsX()+1):
        for hist in h_corr_iters:
            if ibin == 1:
                hist.Write()
            point_positions[ibin] = point_positions.get(ibin, []) + [hist.GetBinContent(ibin)]

    hist_systematics = h_corr_iters[0].Clone('hSystematics')
    hist_systematics.Reset()
    for ibin in range(1, hist_systematics.GetNbinsX()+1):
        point_systematics = np.std(point_positions[ibin])
        hist_systematics.SetBinContent(ibin, point_systematics)
    
    outFile.cd()
    hist_systematics.Write()
    
    outFile.Close()

def run_indivisual_systematics():

    rdf_same, _chain_same = load_same()
    rdf_mixed, _chain_mixed = load_mixed()

    outFile = TFile("output/systematics_individual.root", "RECREATE")
    h_corr_iters = {}

    # Nominal values for selections
    dcaxy_he3_max = 0.1
    dcaz_he3_max = 1.0
    n_sigma_tpc_he3_high = 2.5
    n_sigma_tpc_he3_low = -1.5
    n_sigma_tpc_had_max = 2.
    n_sigma_tof_had_max = 2.

    systematic_variables = ['dcaxy_he3_max', 'dcaz_he3_max', 'n_sigma_tpc_he3',
                          'n_sigma_tpc_had_max', 'n_sigma_tof_had_max']

    N_ITERATIONS = 50
    with alive_bar(N_ITERATIONS*len(systematic_variables), title='Running systematics...') as bar:
        for systematic_variable in systematic_variables:
            
            h_corr_iters_variable = []
            
            # Nominal values for selections
            dcaxy_he3_max = 0.1
            dcaz_he3_max = 1.0
            n_sigma_tpc_he3_high = 2.5
            n_sigma_tpc_he3_low = -1.5
            n_sigma_tpc_had_max = 2.
            n_sigma_tof_had_max = 2.

            for iter in range(N_ITERATIONS):

                if systematic_variable == 'dcaxy_he3_max':
                    dcaxy_he3_max = np.random.uniform(0.05, 0.15)
                elif systematic_variable == 'dcaz_he3_max':
                    dcaz_he3_max = np.random.uniform(0.75, 1.0)
                elif systematic_variable == 'n_sigma_tpc_he3':
                    n_sigma_tpc_he3_high = np.random.uniform(2., 3.)
                    n_sigma_tpc_he3_low = np.random.uniform(-2.0, -1.0)
                elif systematic_variable == 'n_sigma_tpc_had_max':
                    n_sigma_tpc_had_max = np.random.uniform(1.5, 2.5)
                elif systematic_variable == 'n_sigma_tof_had_max':
                    n_sigma_tof_had_max = np.random.uniform(1.5, 2.5)

                condition = f'((std::abs(fDCAxyHe3) < {dcaxy_he3_max}) && (std::abs(fDCAzHe3) < {dcaz_he3_max}) && ' \
                            f'((fNSigmaTPCHe3 < {n_sigma_tpc_he3_high}) && (fNSigmaTPCHe3 > {n_sigma_tpc_he3_low})) && ' \
                            f'((fNSigmaTPCHe3 < {n_sigma_tpc_had_max}) && (fNSigmaTPCHe3 > -{n_sigma_tpc_had_max})) && ' \
                            f'((fNSigmaTOFHad < {n_sigma_tof_had_max}) && (fNSigmaTOFHad > -{n_sigma_tof_had_max})) && ' \
                            f'(fCentralityFT0C < 50))'

                h_same_iter = rdf_same.Filter(condition).Histo1D((f"hKstarSameIter{iter}", ";#it{k}* (GeV/#it{c});", 100, 0, 1.), "fKstar").GetValue()
                h_mixed_iter = rdf_mixed.Filter(condition).Histo1D((f"hKstarMixedIter{iter}", ";#it{k}* (GeV/#it{c});", 100, 0, 1.), "fKstar").GetValue()

                low_bin = h_same_iter.FindBin(0.25)
                high_bin = h_same_iter.FindBin(1)
                normalization_factor = h_same_iter.Integral(low_bin, high_bin) / h_mixed_iter.Integral(low_bin, high_bin)
                h_mixed_iter.Scale(normalization_factor)

                h_corr_iter = h_same_iter.Clone(f'hCorrelationIter{iter}')
                h_corr_iter.Divide(h_mixed_iter)

                h_corr_iters_variable.append(h_corr_iter)

                bar()
            
            h_corr_iters[systematic_variable] = h_corr_iters_variable
    
    outFile.cd()
    hist_systematics_variables = []
    for variable in systematic_variables:
        outFile.mkdir(f'{variable}')
        point_positions = {}

        for ibin in range(1, h_corr_iters[variable][0].GetNbinsX()+1):
            for hist in h_corr_iters[variable]:
                if ibin == 1:
                    outFile.cd(f'{variable}')
                    hist.Write()
                point_positions[ibin] = point_positions.get(ibin, []) + [hist.GetBinContent(ibin)]

        hist_systematics = h_corr_iters[variable][0].Clone(f'hSystematics{variable}')
        hist_systematics.Reset()
        for ibin in range(1, hist_systematics.GetNbinsX()+1):
            point_systematics = np.std(point_positions[ibin])
            hist_systematics.SetBinContent(ibin, point_systematics)
        hist_systematics_variables.append(hist_systematics)
    
    outFile.cd()
    for hist_systematics in hist_systematics_variables:
        hist_systematics.Write()
    
    outFile.Close()

def display_systematics():

    ROOT.gStyle.SetOptStat(0)

    infile = TFile("output/systematics.root", "READ")
    h_systematics = infile.Get("hSystematics")
    h_systematics.SetDirectory(0)

    #infile_correlation = TFile("output/correlation.root", "READ")
    #h_correlation = infile_correlation.Get("Correlation/hCorrelation050")
    #h_correlation.SetDirectory(0)

    infile_correlation = TFile("/home/galucia/antiLithium4/analysis/output/PbPb/studies.root", "READ")
    h_correlation = infile_correlation.Get("CorrelationAnti/hCorrelation_kstar")
    h_correlation.SetDirectory(0)

    h_relative_systematics = h_systematics.Clone("hRelativeSystematics")
    h_relative_statistical = h_correlation.Clone("hRelativeStatistical")
    ratio = h_relative_systematics.Clone("hRatio")

    for ibin in range(1, h_relative_systematics.GetNbinsX() + 1):
        correlation = h_correlation.GetBinContent(ibin)
        systematic = h_systematics.GetBinContent(ibin)
        statistical = h_correlation.GetBinError(ibin)
        _ratio = systematic / statistical if statistical != 0 else 0

        h_relative_systematics.SetBinContent(ibin, systematic / correlation if correlation != 0 else 0)
        h_relative_statistical.SetBinContent(ibin, statistical / correlation if correlation != 0 else 0)
        ratio.SetBinContent(ibin, _ratio)
        #ratio.SetBinError(ibin, 1)

        print(f"Bin {ibin}: Corr = {correlation}, Stat = {statistical:.4f}, Syst = {systematic:.4f}, Ratio = {_ratio:.4f}")

    canvas = ROOT.TCanvas("canvas", "Systematics Analysis", 800, 600)

    pad = ROOT.TPad("pad", "pad", 0, 0.03, 1, 1)
    pad.SetBottomMargin(0.15)
    pad.Draw()

    pad_ratio = ROOT.TPad("pad_ratio", "pad_ratio", 0, 0, 1, 0.29)
    pad_ratio.SetTopMargin(0.0)
    pad_ratio.SetBottomMargin(0.3)
    pad_ratio.Draw()

    h_relative_statistical.SetLineColor(ROOT.kBlue+2)
    h_relative_statistical.SetLineWidth(2)
    h_relative_statistical.SetTitle(";#it{k}* (GeV/#it{c});#sigma/ C(#it{k}*)")

    h_relative_systematics.SetLineColor(ROOT.kOrange-3)
    h_relative_systematics.SetLineWidth(2)
    h_relative_systematics.SetTitle(";#it{k}* (GeV/#it{c});#sigma/ C(#it{k}*)")

    legend = ROOT.TLegend(0.25, 0.7, 0.55, 0.8)
    legend.SetFillColor(0)
    legend.SetBorderSize(0)
    legend.AddEntry(h_relative_statistical, "#sigma_{stat}/ C(#it{k}*)", "l")
    legend.AddEntry(h_relative_systematics, "#sigma_{syst}/ C(#it{k}*)", "l")

    text = ROOT.TPaveText(0.55, 0.7, 0.85, 0.85, "NDC")
    text.SetFillColor(0)
    text.SetBorderSize(0)
    text.AddText("This work")
    text.AddText("#bf{ALICE, Run 3}")
    text.AddText("#bf{Pb-Pb #sqrt{s_{NN}} = 5.36 TeV}")
    
    pad.cd()
    hframe = pad.DrawFrame(0, 1e-4, 0.39, 1, ";#it{k}* (GeV/#it{c});#sigma/ C(#it{k}*)")
    h_relative_statistical.Draw("HIST same")
    h_relative_systematics.Draw("HIST SAME")
    legend.Draw('same')
    text.Draw('same')
    pad.SetLogy()

    ratio.SetTitle(";#it{k}* (GeV/#it{c});#sigma_{syst}/#sigma_{stat}")
    ratio.SetMarkerStyle(20)
    ratio.SetMarkerSize(1.5)
    ratio.SetMarkerColor(ROOT.kOrange-3)
    ratio.SetLineColor(ROOT.kGreen+2)
    ratio.SetLineWidth(2)

    #line = ROOT.TLine(0, 1, 0.39, 1)
    #line.SetLineColor(ROOT.kGray+2)
    #line.SetLineStyle(2)
    #line.SetLineWidth(2)

    fit = ROOT.TF1("fit", "pol0", 0, 0.39)
    fit.SetLineColor(ROOT.kGray+2)
    fit.SetLineStyle(2)
    fit.SetLineWidth(2)
    ratio.Fit(fit, "rms+", "")

    pad_ratio.cd()
    hframe_ratio = pad_ratio.DrawFrame(0, 1e-2, 0.39, 2, ";#it{k}* (GeV/#it{c});#sigma_{syst}/#sigma_{stat}")
    hframe_ratio.GetYaxis().SetTitleOffset(0.5)
    hframe_ratio.GetYaxis().SetTitleSize(0.1)
    hframe_ratio.GetXaxis().SetTitleSize(0.1)
    hframe_ratio.GetXaxis().SetTitleOffset(0.9)
    hframe_ratio.GetXaxis().SetLabelSize(0.08)
    hframe_ratio.GetYaxis().SetLabelSize(0.08)

    #line.Draw("same")
    fit.Draw("same")
    ratio.Draw("hist e1 same")
    pad_ratio.SetLogy()

    canvas.SaveAs("output/systematics.pdf")

    #########################################################

    xs, ys, ex, eystat, eysyst = [], [], [], [], []
    for ibin in range(1, h_correlation.GetNbinsX() + 1):
        xs.append(h_correlation.GetBinCenter(ibin))
        ys.append(h_correlation.GetBinContent(ibin))
        ex.append(h_correlation.GetBinWidth(ibin) / 2)
        eystat.append(h_correlation.GetBinError(ibin))
        eysyst.append(h_systematics.GetBinContent(ibin))

    correlation_systematics = ROOT.TGraphMultiErrors('gme', ';#it{k}* (GeV/#it{c});C(#it{k}*)',
                                                    len(xs), np.array(xs), np.array(ys),
                                                    np.array(ex), np.array(ex), np.array(eystat), np.array(eystat))
    correlation_systematics.AddYError(len(xs), np.array(eysyst), np.array(eysyst))
    
    canvas_correlation = ROOT.TCanvas("canvas_correlation", "Correlation", 800, 600)
    correlation_systematics.SetTitle(";#it{k}* (GeV/#it{c});C(#it{k}*)")
    correlation_systematics.SetMarkerStyle(20)
    correlation_systematics.SetMarkerColor(ROOT.kBlue+2)
    correlation_systematics.SetLineColor(ROOT.kOrange-3)
    correlation_systematics.GetAttLine(0).SetLineColor(ROOT.kOrange-3)
    correlation_systematics.GetAttLine(1).SetLineColor(ROOT.kGreen+2)
    correlation_systematics.GetAttFill(1).SetFillStyle(0)
    correlation_systematics.Draw("APS;;5")

    text = ROOT.TPaveText(0.55, 0.35, 0.85, 0.55, "NDC")
    text.SetFillColor(0)
    text.SetBorderSize(0)
    text.AddText("This work")
    text.AddText("#bf{ALICE, Run 3}")
    text.AddText("#bf{Pb-Pb #sqrt{s_{NN}} = 5.36 TeV}")
    text.AddText("#bf{p-^{3}He #oplus #bar{p}-^{3}#bar{He}}")
    text.Draw('same')

    canvas_correlation.SaveAs("output/correlation_systematics.pdf")

def display_individual_systematics():

    ROOT.gStyle.SetOptStat(0)

    systematic_variables = {'DCAxy (^{3}He)':'dcaxy_he3_max', 
                            'DCAz (^{3}He)': 'dcaz_he3_max', 
                            'n#sigma_{TPC} (^{3}He)': 'n_sigma_tpc_he3',
                            'n#sigma_{TPC} (p)': 'n_sigma_tpc_had_max', 
                            'n#sigma_{TOF} (p)': 'n_sigma_tof_had_max'}
    colors = [ROOT.kRed+1, ROOT.kBlue+2, ROOT.kGreen+2, ROOT.kMagenta+2, ROOT.kCyan+2]

    infile = TFile("output/systematics_individual.root", "READ")

    #infile_correlation = TFile("output/correlation.root", "READ")
    #h_correlation = infile_correlation.Get("Correlation/hCorrelation050")
    #h_correlation.SetDirectory(0)

    infile_correlation = TFile("/home/galucia/antiLithium4/analysis/output/PbPb/studies.root", "READ")
    h_correlation = infile_correlation.Get("CorrelationAnti/hCorrelation_kstar")
    h_correlation.SetDirectory(0)

    hs_relative_systematics = []

    for variable in systematic_variables.values():

        h_systematics = infile.Get(f"hSystematics{variable}")
        h_systematics.SetDirectory(0)
        h_relative_systematics = h_systematics.Clone(f"hRelativeSystematics{variable}")

        for ibin in range(1, h_relative_systematics.GetNbinsX() + 1):
            correlation = h_correlation.GetBinContent(ibin)
            systematic = h_systematics.GetBinContent(ibin)

            h_relative_systematics.SetBinContent(ibin, systematic / correlation if correlation != 0 else 0)

        hs_relative_systematics.append(h_relative_systematics)

    canvas = ROOT.TCanvas("canvas", "Systematics Analysis", 800, 600)
    hframe = canvas.DrawFrame(0, 1e-5, 0.39, 1, ";#it{k}* (GeV/#it{c});#sigma_{syst}/ C(#it{k}*)")

    legend = ROOT.TLegend(0.55, 0.5, 0.85, 0.69)
    legend.SetFillColor(0)
    legend.SetBorderSize(0)
    legend.SetNColumns(2)

    canvas.cd()
    for h_relative_systematics, variable_title, color in zip(hs_relative_systematics, systematic_variables.keys(), colors):
        h_relative_systematics.SetLineColor(color)
        h_relative_systematics.SetLineWidth(2)
        h_relative_systematics.Draw("HIST SAME")
        legend.AddEntry(h_relative_systematics, variable_title, "l")
    
    text = ROOT.TPaveText(0.55, 0.7, 0.85, 0.85, "NDC")
    text.SetFillColor(0)
    text.SetBorderSize(0)
    text.AddText("This work")
    text.AddText("#bf{ALICE, Run 3}")
    text.AddText("#bf{Pb-Pb #sqrt{s_{NN}} = 5.36 TeV}")
    
    canvas.cd()
    legend.Draw('same')
    text.Draw('same')
    canvas.SetLogy()

    canvas.SaveAs("output/systematics_individual.pdf")

    #########################################################

    xs, ys, ex, eystat, eysyst = [], [], [], [], []
    for ibin in range(1, h_correlation.GetNbinsX() + 1):
        xs.append(h_correlation.GetBinCenter(ibin))
        ys.append(h_correlation.GetBinContent(ibin))
        ex.append(h_correlation.GetBinWidth(ibin) / 2)
        eystat.append(h_correlation.GetBinError(ibin))
        eysyst.append(h_systematics.GetBinContent(ibin))

    correlation_systematics = ROOT.TGraphMultiErrors('gme', ';#it{k}* (GeV/#it{c});C(#it{k}*)',
                                                    len(xs), np.array(xs), np.array(ys),
                                                    np.array(ex), np.array(ex), np.array(eystat), np.array(eystat))
    correlation_systematics.AddYError(len(xs), np.array(eysyst), np.array(eysyst))
    
    canvas_correlation = ROOT.TCanvas("canvas_correlation", "Correlation", 800, 600)
    correlation_systematics.SetTitle(";#it{k}* (GeV/#it{c});C(#it{k}*)")
    correlation_systematics.SetMarkerStyle(20)
    correlation_systematics.SetMarkerColor(ROOT.kBlue+2)
    correlation_systematics.SetLineColor(ROOT.kOrange-3)
    correlation_systematics.GetAttLine(0).SetLineColor(ROOT.kOrange-3)
    correlation_systematics.GetAttLine(1).SetLineColor(ROOT.kGreen+2)
    correlation_systematics.GetAttFill(1).SetFillStyle(0)
    correlation_systematics.Draw("APS;;5")

    text = ROOT.TPaveText(0.55, 0.35, 0.85, 0.55, "NDC")
    text.SetFillColor(0)
    text.SetBorderSize(0)
    text.AddText("This work")
    text.AddText("#bf{ALICE, Run 3}")
    text.AddText("#bf{Pb-Pb #sqrt{s_{NN}} = 5.36 TeV}")
    text.AddText("#bf{p-^{3}He #oplus #bar{p}-^{3}#bar{He}}")
    text.Draw('same')

    canvas_correlation.SaveAs("output/correlation_systematics.pdf")


if __name__ == "__main__":

    #run_systematics()
    #run_indivisual_systematics()
    display_systematics()
    display_individual_systematics()
    print("Systematics analysis completed successfully.")