/**
 * systematics_with_upper_limit.cpp
 *
 * C++ translation of systematics_with_upper_limit.py
 *
 * Dependencies: ROOT (TFile, TChain, RDataFrame, TH1D, TF1, RooWorkspace, etc.)
 *
 * Compile example:
 *   g++ -O2 -std=c++17 systematics_with_upper_limit.cpp \
 *       $(root-config --cflags --libs) -lRooFit -lRooFitCore \
 *       -I../include -o systematics_with_upper_limit
 */

#include <iostream>
#include <cmath>
#include <map>
#include <memory>
#include <random>
#include <string>
#include <vector>

#include "TFile.h"
#include "TChain.h"
#include "TKey.h"
#include "TH1D.h"
#include "TF1.h"
#include "TDirectory.h"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RDF/RVariationsDescription.hxx"
#include "RooWorkspace.h"
#include "RooMsgService.h"

#include "../include/Common.h"
#include "../include/Variation.h"
#include "../include/Systematics.h"

#include "SignalFitter.h"
#include "BkgFitter.h"
#include "ModelFitter.h"

static constexpr double kMassHe  = 2.80923; // GeV/c^2  (He-3 nucleus)
static constexpr double kMassPr  = 0.93827; // GeV/c^2  (proton)

static const std::string kBaseSelection =
    "((fSignedPtHe3 > 0 && fSignedPtHad > 0) || (fSignedPtHe3 < 0 && fSignedPtHad < 0))"
    " && ((fChi2TPCHe3 > 0.5) || (fIs23 == false))"
    " && (fChi2TPCHe3 < 4)"
    " && (fChi2TPCHad < 4)"
    " && (std::abs(fEtaHe3) < 0.9)"
    " && (std::abs(fEtaHad) < 0.9)"
    " && (fNClsTPCHe3 > 110)"
    " && ((fPtHe3 < 2.5) || (fPIDtrkHe3 == 7))"
    " && ((fPIDtrkHe3 == 6) || (fPIDtrkHe3 == 7) || (fPIDtrkHe3 == 8))"
    " && (std::abs(fNSigmaDCAxyHe3) < 3)"
    " && (std::abs(fNSigmaDCAzHe3) < 3)"
    " && (std::abs(fNSigmaDCAxyHad) < 3)"
    " && (std::abs(fNSigmaDCAzHad) < 3)";

struct ChainConfig {
    std::vector<std::string> input_data;
    std::string tree_name;
    std::string mode; // "DF" or "tree"
};

TChain* PrepareInputTChain(const ChainConfig& config)
{
    auto* chain = new TChain("tchain");

    for (const auto& file_name : config.input_data) {
        TFile* f = TFile::Open(file_name.c_str());
        if (!f || f->IsZombie()) {
            std::cerr << "Cannot open file: " << file_name << std::endl;
            continue;
        }
        if (config.mode == "DF") {
            TIter next(f->GetListOfKeys());
            TKey* key = nullptr;
            while ((key = dynamic_cast<TKey*>(next()))) {
                std::string key_name = key->GetName();
                if (key_name.find("DF_") != std::string::npos) {
                    std::string full_path = file_name + "/" + key_name + "/" + config.tree_name;
                    std::cout << "Adding " << full_path << " to chain\n";
                    chain->Add(full_path.c_str());
                }
            }
        } else { // "tree"
            std::string full_path = file_name + "/" + config.tree_name;
            std::cout << "Adding " << full_path << " to chain\n";
            chain->Add(full_path.c_str());
        }
        f->Close();
        delete f;
    }
    return chain;
}

ROOT::RDF::RNode PrepareRDataFrame(TChain* chain,
                                   const std::string& base_selection,
                                   const std::string& extra_selection)
{
    ROOT::RDataFrame rdf(*chain);

    // Conditionally define fNSigmaTPCHad (checked at runtime via column names)
    // In C++ we cannot easily branch at runtime on column existence in a pure
    // RDataFrame chain, so we use a lambda approach with HasColumn.
    auto col_names = rdf.GetColumnNames();
    auto has_col = [&col_names](const std::string& name) {
        for (auto& c : col_names) if (c == name) return true;
        return false;
    };

    ROOT::RDF::RNode rdf1 = has_col("fNSigmaTPCHadPr")
        ? rdf.Define("fNSigmaTPCHad", "fNSigmaTPCHadPr")
        : rdf; // already defined upstream

    ROOT::RDF::RNode rdf2 = has_col("fNSigmaTOFHadPr")
        ? rdf1.Define("fNSigmaTOFHad", "fNSigmaTOFHadPr")
        : rdf1.Define("fNSigmaTOFHad",
              "ComputeNsigmaTOFPr(std::abs(fPtHad), fMassTOFHad)");

    ROOT::RDF::RNode rdf3 = rdf2
        .Define  ("fSignedPtHad",   "fPtHad")
        .Define  ("fSignHe3",       "fPtHe3/std::abs(fPtHe3)")
        .Redefine("fPtHe3",         "std::abs(fPtHe3)")
        .Redefine("fPtHad",         "std::abs(fPtHad)")
        .Redefine("fPtHe3",
            "(fPIDtrkHe3 == 7) || (fPIDtrkHe3 == 8) || (fPtHe3 > 2.5)"
            " ? fPtHe3 : CorrectPidTrkHe(fPtHe3)")
        .Define  ("fSignedPtHe3",   "fPtHe3 * fSignHe3")
        .Define  ("fEHe3",
            "std::sqrt((fPtHe3*std::cosh(fEtaHe3))*(fPtHe3*std::cosh(fEtaHe3))"
            " + " + std::to_string(kMassHe) + "*" + std::to_string(kMassHe) + ")")
        .Define  ("fEHad",
            "std::sqrt((fPtHad*std::cosh(fEtaHad))*(fPtHad*std::cosh(fEtaHad))"
            " + " + std::to_string(kMassPr) + "*" + std::to_string(kMassPr) + ")")
        .Define  ("fDeltaEta",              "fEtaHe3 - fEtaHad")
        .Define  ("fDeltaPhi",              "fPhiHe3 - fPhiHad")
        .Redefine("fInnerParamTPCHe3",      "fInnerParamTPCHe3 * 2")
        .Define  ("fClusterSizeCosLamHe3",
            "ComputeAverageClusterSize(fItsClusterSizeHe3) / cosh(fEtaHe3)")
        .Define  ("fClusterSizeCosLamHad",
            "ComputeAverageClusterSize(fItsClusterSizeHad) / cosh(fEtaHad)")
        .Define  ("fNSigmaITSHe3",
            "ComputeNsigmaITSHe(fPtHe3 * std::cosh(fEtaHe3), fClusterSizeCosLamHe3)")
        .Define  ("fNSigmaITSHad",
            "ComputeNsigmaITSPr(fPtHad * std::cosh(fEtaHad), fClusterSizeCosLamHad)")
        .Redefine("fNSigmaTPCHe3",
            "ComputeNsigmaTPCHe(std::abs(fInnerParamTPCHe3), fSignalTPCHe3)")
        .Define  ("fNSigmaDCAxyHe3", "ComputeNsigmaDCAxyHe(fPtHe3, fDCAxyHe3)")
        .Define  ("fNSigmaDCAzHe3",  "ComputeNsigmaDCAzHe(fPtHe3, fDCAzHe3)")
        .Define  ("fNSigmaDCAxyHad", "ComputeNsigmaDCAxyPr(fPtHad, fDCAxyHad)")
        .Define  ("fNSigmaDCAzHad",  "ComputeNsigmaDCAzPr(fPtHad, fDCAzHad)")
        .Filter  (base_selection)
        .Filter  (extra_selection)
        .Filter  ("fCentralityFT0C < 50")
        .Define  ("fKstar",
            "ComputeKstar(fPtHe3, fEtaHe3, fPhiHe3, " + std::to_string(kMassHe) +
            ", fPtHad, fEtaHad, fPhiHad, " + std::to_string(kMassPr) + ")");

    return rdf3;
}

std::pair<ROOT::RDF::RNode, TChain*> LoadSame()
{
    ChainConfig cfg;
    cfg.input_data = {
        "/data/galucia/lithium_local/same_merged/LHC23_PbPb_pass4_hadronpid_same.root",
        "/data/galucia/lithium_local/same_merged/LHC24ar_pass1_hadronpid_same.root",
        "/data/galucia/lithium_local/same_merged/LHC24as_pass1_hadronpid_same.root",
    };
    cfg.tree_name = "O2he3hadtable";
    cfg.mode = "DF";

    TChain* chain = PrepareInputTChain(cfg);
    auto rdf = PrepareRDataFrame(chain, kBaseSelection, "true");
    return {rdf, chain};
}

std::pair<ROOT::RDF::RNode, TChain*> LoadMixed()
{
    ChainConfig cfg;
    cfg.input_data = {
        "/data/galucia/lithium_local/mixing/LHC23_PbPb_pass4_hadronpid_event_mixing_batch1995_refined_dca.root",
        "/data/galucia/lithium_local/mixing/LHC24ar_pass1_hadronpid_event_mixing_batch1995_refined_dca.root",
        "/data/galucia/lithium_local/mixing/LHC24as_pass1_hadronpid_event_mixing_batch1995_refined_dca.root",
        "/data/galucia/lithium_local/mixing/LHC23_PbPb_pass4_hadronpid_event_mixing_batch42_refined_dca.root",
        "/data/galucia/lithium_local/mixing/LHC24ar_pass1_hadronpid_event_mixing_batch42_refined_dca.root",
        "/data/galucia/lithium_local/mixing/LHC24as_pass1_hadronpid_event_mixing_batch42_refined_dca.root",
        "/data/galucia/lithium_local/mixing/LHC23_PbPb_pass4_hadronpid_event_mixing_batch256_refined_dca.root",
        "/data/galucia/lithium_local/mixing/LHC24ar_pass1_hadronpid_event_mixing_batch256_refined_dca.root",
        "/data/galucia/lithium_local/mixing/LHC24as_pass1_hadronpid_event_mixing_batch256_refined_dca.root",
        "/data/galucia/lithium_local/mixing/LHC23_PbPb_pass4_hadronpid_event_mixing_batch3112_refined_dca.root",
        "/data/galucia/lithium_local/mixing/LHC24ar_pass1_hadronpid_event_mixing_batch3112_refined_dca.root",
        "/data/galucia/lithium_local/mixing/LHC24as_pass1_hadronpid_event_mixing_batch3112_refined_dca.root",
    };
    cfg.tree_name = "MixedTree";
    cfg.mode = "tree";

    TChain* chain = PrepareInputTChain(cfg);
    auto rdf = PrepareRDataFrame(chain, kBaseSelection, "true");
    return {rdf, chain};
}

void PoissonSampling(TH1D* hist)
{
    static std::mt19937 rng(std::random_device{}());
    for (int ibin = 1; ibin <= hist->GetNbinsX(); ++ibin) {
        double content = hist->GetBinContent(ibin);
        std::poisson_distribution<long> dist(content > 0 ? content : 0);
        long new_content = dist(rng);
        hist->SetBinContent(ibin, static_cast<double>(new_content));
        hist->SetBinError(ibin, std::sqrt(static_cast<double>(new_content)));
    }
}

//using VarHistMap = std::map<std::string, ROOT::RDF::RResultPtr<TH1D>>;
using VarHistMap = std::map<std::string, TH1D>;

std::tuple<VarHistMap, VarHistMap, VarHistMap>
PrepareSystematicsHistograms(ROOT::RDF::RNode& rdf,
                             int n_variations,
                             const std::string& hist_name_suffix,
                             const std::string& variable)
{
    // Build the vary expression depending on 'variable'
    const std::string vary_cols =
        "fNSigmaITSHe3Systematic, fAbsNSigmaTPCHe3Systematic,"
        " fAbsNSigmaDCAxyHe3Systematic, fAbsNSigmaDCAzHe3Systematic,"
        " fNClsTPCHe3Systematic, fChi2TPCHe3Systematic,"
        " fNSigmaITSHadSystematic, fAbsNSigmaTPCHadSystematic, fAbsNSigmaTOFHadSystematic,"
        " fAbsNSigmaDCAxyHadSystematic, fAbsNSigmaDCAzHadSystematic, fAbsZVertexSystematic";

    std::string vary_expr;
    if (variable == "all") {
        vary_expr = "systematicCuts(" + std::to_string(n_variations) +
            ", fNSigmaITSHe3Systematic, fAbsNSigmaTPCHe3Systematic,"
            " fAbsNSigmaDCAxyHe3Systematic, fAbsNSigmaDCAzHe3Systematic,"
            " fNClsTPCHe3Systematic, fChi2TPCHe3Systematic,"
            " fNSigmaITSHadSystematic, fAbsNSigmaTPCHadSystematic, fAbsNSigmaTOFHadSystematic,"
            " fAbsNSigmaDCAxyHadSystematic, fAbsNSigmaDCAzHadSystematic, fChi2TPCHadSystematic,"
            " fAbsZVertexSystematic)";
    } else {
        vary_expr = "systematicCutsSingleVariable(" + std::to_string(n_variations) +
            ", \"" + variable + "\","
            " fNSigmaITSHe3Systematic, fAbsNSigmaTPCHe3Systematic,"
            " fAbsNSigmaDCAxyHe3Systematic, fAbsNSigmaDCAzHe3Systematic,"
            " fNClsTPCHe3Systematic, fChi2TPCHe3Systematic,"
            " fNSigmaITSHadSystematic, fAbsNSigmaTPCHadSystematic, fAbsNSigmaTOFHadSystematic,"
            " fAbsNSigmaDCAxyHadSystematic, fAbsNSigmaDCAzHadSystematic, fChi2TPCHadSystematic,"
            " fAbsZVertexSystematic)";
    }

    const std::string syst_filter =
        "(fNSigmaITSHe3Systematic > 0)"
        " && (fAbsNSigmaTPCHe3Systematic < 0)"
        " && (fNSigmaITSHadSystematic > 0)"
        " && (fAbsNSigmaTPCHadSystematic < 0)"
        " && ((std::abs(fPtHad) < 0.8) || (fAbsNSigmaTOFHadSystematic < 0))";

    ROOT::RDF::RNode vrdf = rdf
        .Define  ("fNSigmaITSHe3Systematic",        "fNSigmaITSHe3")
        .Define  ("fAbsNSigmaTPCHe3Systematic",      "std::abs(fNSigmaTPCHe3)")
        .Define  ("fAbsNSigmaDCAxyHe3Systematic",    "std::abs(fNSigmaDCAxyHe3)")
        .Define  ("fAbsNSigmaDCAzHe3Systematic",     "std::abs(fNSigmaDCAzHe3)")
        .Define  ("fNClsTPCHe3Systematic",           "static_cast<float>(fNClsTPCHe3)")
        .Define  ("fChi2TPCHe3Systematic",           "fChi2TPCHe3")
        .Define  ("fNSigmaITSHadSystematic",         "fNSigmaITSHad")
        .Define  ("fAbsNSigmaTPCHadSystematic",      "std::abs(fNSigmaTPCHad)")
        .Define  ("fAbsNSigmaTOFHadSystematic",      "std::abs(fNSigmaTOFHad)")
        .Define  ("fAbsNSigmaDCAxyHadSystematic",    "std::abs(fNSigmaDCAxyHad)")
        .Define  ("fAbsNSigmaDCAzHadSystematic",     "std::abs(fNSigmaDCAzHad)")
        .Define  ("fChi2TPCHadSystematic",           "fChi2TPCHad")
        .Define  ("fAbsZVertexSystematic",           "std::abs(fZVertex)")
        .Vary    ({  "fNSigmaITSHe3Systematic", "fAbsNSigmaTPCHe3Systematic",
                     "fAbsNSigmaDCAxyHe3Systematic", "fAbsNSigmaDCAzHe3Systematic",
                     "fNClsTPCHe3Systematic", "fChi2TPCHe3Systematic",
                     "fNSigmaITSHadSystematic", "fAbsNSigmaTPCHadSystematic",
                     "fAbsNSigmaTOFHadSystematic",
                     "fAbsNSigmaDCAxyHadSystematic", "fAbsNSigmaDCAzHadSystematic",
                     "fAbsZVertexSystematic" },
                  vary_expr, n_variations, "systematic")
        .Filter  (syst_filter);

    auto nominal_010 = vrdf.Filter("fCentralityFT0C < 10")
        .Histo1D({("hKstar010" + hist_name_suffix).c_str(),
                  ";#it{k}* (GeV/#it{c});", 20, 0., 0.4}, "fKstar");
    auto varHists010 = ROOT::RDF::Experimental::VariationsFor(nominal_010);

    auto nominal_1030 = vrdf.Filter("fCentralityFT0C >= 10 && fCentralityFT0C < 30")
        .Histo1D({("hKstar1030" + hist_name_suffix).c_str(),
                  ";#it{k}* (GeV/#it{c});", 20, 0., 0.4}, "fKstar");
    auto varHists1030 = ROOT::RDF::Experimental::VariationsFor(nominal_1030);

    auto nominal_3050 = vrdf.Filter("fCentralityFT0C >= 30 && fCentralityFT0C < 50")
        .Histo1D({("hKstar3050" + hist_name_suffix).c_str(),
                  ";#it{k}* (GeV/#it{c});", 20, 0., 0.4}, "fKstar");
    auto varHists3050 = ROOT::RDF::Experimental::VariationsFor(nominal_3050);

    // Convert RVariationsFor results to plain maps of TH1D result-ptrs.
    // (RVariationsFor returns an RResultMap<TH1D>; access via operator[])
    VarHistMap map010, map1030, map3050;
    for (int i = 0; i < n_variations; ++i) {
        std::string key = "systematic:" + std::to_string(i);
        map010[key]   = varHists010[key];
        map1030[key]  = varHists1030[key];
        map3050[key]  = varHists3050[key];
    }
    return {map010, map1030, map3050};
}

TH1D* NormaliseHistogram(TH1D* h_same, TH1D* h_mixed,
                          double norm_low, double norm_high)
{
    int low_bin  = h_same->FindBin(norm_low);
    int high_bin = h_same->FindBin(norm_high);
    double norm  = h_same->Integral(low_bin, high_bin)
                 / h_mixed->Integral(low_bin, high_bin);

    TH1D* h_norm = static_cast<TH1D*>(
        h_mixed->Clone((std::string(h_mixed->GetName()) + "_Normalised").c_str()));
    h_norm->Scale(norm);
    return h_norm;
}

TH1D* CorrelationFunctionCentralityIntegrated(
    const std::vector<TH1D*>& h_sames,
    const std::vector<TH1D*>& h_mixeds,
    const std::string& suffix)
{
    TH1D* h_same  = static_cast<TH1D*>(h_sames[0]->Clone());
    TH1D* h_mixed = static_cast<TH1D*>(h_mixeds[0]->Clone());

    for (size_t i = 1; i < h_sames.size(); ++i) {
        h_same->Add(h_sames[i]);
        h_mixed->Add(h_mixeds[i]);
    }

    TH1D* h_corr = static_cast<TH1D*>(
        h_same->Clone(("hCorrelation" + suffix).c_str()));
    h_corr->Divide(h_mixed);

    delete h_same;
    delete h_mixed;
    return h_corr;
}

double FittingRoutine(TDirectory* outdir,
                      const std::string& bkg_input_path,
                      const std::string& h_bkg_name,
                      TH1D* h_data,
                      TH1D* h_mixed_event,
                      const std::string& sign,
                      const std::string& centrality,
                      bool use_smoothening = false)
{
    const std::string signal_path =
        "/home/galucia/Lithium4/femto/models/li4_contribution_finer_binning.root";
    const std::string signal_hist_name = "hCkHist";

    bool is_antimatter_010 = (centrality.find("010") != std::string::npos && sign == "Antimatter");
    double kstar_min = is_antimatter_010 ? 0.02 : 0.01;
    double kstar_max = 0.4;

    // Load background and signal histograms
    TFile* f_bkg = TFile::Open(bkg_input_path.c_str(), "READ");
    TH1D*  h_bkg = static_cast<TH1D*>(f_bkg->Get(h_bkg_name.c_str()));

    TFile* f_sig = TFile::Open(signal_path.c_str(), "READ");
    TH1D*  h_sig = static_cast<TH1D*>(f_sig->Get(signal_hist_name.c_str()));

    RooWorkspace workspace("roows");

    // NOTE: SignalFitter, BkgFitter, ModelFitter are C++ classes from the
    // femto library – call their constructors/methods as defined there.
    // The logic below mirrors the Python code exactly.

    SignalFitter signal_fitter(std::string("signal"), std::string("kstar"), std::string(";#it{k}* (GeV/#it{c});"),
                               kstar_min, kstar_max, outdir, &workspace);
    std::string signal_init_mode = use_smoothening ? "from_kde" : "from_mc";
    signal_fitter.init_signal(signal_init_mode, h_sig);
    signal_fitter.set_title("^{4}Li");
    signal_fitter.save_to_workspace();

    double bkg_rho = (centrality.find("010") != std::string::npos) ? 0.1 : 0.05;
    BkgFitter bkg_fitter(std::string("bkg"), std::string("kstar"), std::string(";#it{k}* (GeV/#it{c});"),
                         kstar_min, kstar_max, outdir, &workspace);
    std::string bkg_init_mode = use_smoothening ? "from_kde" : "from_mc";
    bkg_fitter.init_bkg(bkg_init_mode, h_bkg, "bkg_pdf", /*extended=*/true, kstar_min, kstar_max, bkg_rho);
    bkg_fitter.set_title("Full model");
    bkg_fitter.save_to_workspace();

    ModelFitter model_fitter(std::string("model"), std::string("kstar"), std::string(";#it{k}* (GeV/#it{c});"),
                                kstar_min, kstar_max, outdir,
                                {"signal_pdf"}, {"bkg_pdf"}, &workspace, 
                                /*extended=*/true);
    model_fitter.load_data(h_data, h_data->GetName());
    model_fitter.prefit_background(h_data, 0.2, 0.4, "bkg_fit_range",
                                  /*save_norm=*/true);
    model_fitter.fractions["signal_pdf"]->setVal(0.3);
    model_fitter.fractions["signal_pdf"]->setRange(-1e4, 1e4);
    model_fitter.fit_model(h_data, "signal_pdf", "bkg_fit_range");
    model_fitter.save_to_workspace();
    model_fitter.compute_chi2(h_data);

    auto* xvar = dynamic_cast<RooRealVar*>(workspace.obj("kstar"));
    xvar->setRange(kstar_min, 0.39);
    double raw_yield = model_fitter.compute_raw_yield(h_mixed_event, "signal_pdf", "bkg_pdf");

    f_bkg->Close(); delete f_bkg;
    f_sig->Close(); delete f_sig;

    return raw_yield;
}

void PrepareHistograms()
{
    auto [rdf_same,  chain_same]  = LoadSame();
    auto [rdf_mixed, chain_mixed] = LoadMixed();

    const int N_ITERATIONS = 100;
    TFile* outFile = TFile::Open("output/hist_systematics_with_upper_limit.root", "RECREATE");

    for (const auto& [sign, condition] : std::map<std::string, std::string>{
            {"Matter",    "fSignedPtHe3 > 0"},
            {"Antimatter","fSignedPtHe3 < 0"}})
    {
        TDirectory* outDir = outFile->mkdir(sign.c_str());

        ROOT::RDF::RNode tmp_rdf_same  = rdf_same.Filter(condition);
        ROOT::RDF::RNode tmp_rdf_mixed = rdf_mixed.Filter(condition);

        constexpr double NORM_LOW  = 0.2;
        constexpr double NORM_HIGH = 0.4;

        // --- Reference (n_variations=1) ---
        auto [v010_s_ref, v1030_s_ref, v3050_s_ref] =
            PrepareSystematicsHistograms(tmp_rdf_same,  1, "Same",  "");
        auto [v010_m_ref, v1030_m_ref, v3050_m_ref] =
            PrepareSystematicsHistograms(tmp_rdf_mixed, 1, "Mixed", "");

        //TH1D* same_010_ref  = static_cast<TH1D*>(v010_s_ref.at("systematic:0").GetPtr());
        //TH1D* same_1030_ref = static_cast<TH1D*>(v1030_s_ref.at("systematic:0").GetPtr());
        //TH1D* same_3050_ref = static_cast<TH1D*>(v3050_s_ref.at("systematic:0").GetPtr());
        //TH1D* mixed_010_ref  = static_cast<TH1D*>(v010_m_ref.at("systematic:0").GetPtr());
        //TH1D* mixed_1030_ref = static_cast<TH1D*>(v1030_m_ref.at("systematic:0").GetPtr());
        //TH1D* mixed_3050_ref = static_cast<TH1D*>(v3050_m_ref.at("systematic:0").GetPtr());

        TH1D* same_010_ref  = &v010_s_ref.at("systematic:0");
        TH1D* same_1030_ref = &v1030_s_ref.at("systematic:0");
        TH1D* same_3050_ref = &v3050_s_ref.at("systematic:0");
        TH1D* mixed_010_ref  = &v010_m_ref.at("systematic:0");
        TH1D* mixed_1030_ref = &v1030_m_ref.at("systematic:0");
        TH1D* mixed_3050_ref = &v3050_m_ref.at("systematic:0");

        TH1D* mn_010  = NormaliseHistogram(same_010_ref,  mixed_010_ref,  NORM_LOW, NORM_HIGH);
        TH1D* mn_1030 = NormaliseHistogram(same_1030_ref, mixed_1030_ref, NORM_LOW, NORM_HIGH);
        TH1D* mn_3050 = NormaliseHistogram(same_3050_ref, mixed_3050_ref, NORM_LOW, NORM_HIGH);

        TH1D* h_corr_ref = CorrelationFunctionCentralityIntegrated(
            {same_010_ref, same_1030_ref, same_3050_ref},
            {mn_010, mn_1030, mn_3050}, "0");

        outDir->cd();
        h_corr_ref->Write("hCorrelationReference");
        delete mn_010; delete mn_1030; delete mn_3050; delete h_corr_ref;

        // --- Full systematics (N_ITERATIONS) ---
        auto [v010_s, v1030_s, v3050_s] =
            PrepareSystematicsHistograms(tmp_rdf_same,  N_ITERATIONS, "Same",  "all");
        auto [v010_m, v1030_m, v3050_m] =
            PrepareSystematicsHistograms(tmp_rdf_mixed, N_ITERATIONS, "Mixed", "all");

        struct IterHistSet {
            std::map<std::string, TH1D*> sames, mixeds, mixeds_norm;
        };
        std::vector<IterHistSet> iter_data;

        for (int iter = 0; iter < N_ITERATIONS; ++iter) {
            std::cout << "Processing iteration " << iter+1 << "/" << N_ITERATIONS
                      << " for " << sign << "...\n";

            std::string key = "systematic:" + std::to_string(iter);
            //TH1D* s010  = static_cast<TH1D*>(v010_s.at(key).GetPtr());
            //TH1D* s1030 = static_cast<TH1D*>(v1030_s.at(key).GetPtr());
            //TH1D* s3050 = static_cast<TH1D*>(v3050_s.at(key).GetPtr());
            //TH1D* m010  = static_cast<TH1D*>(v010_m.at(key).GetPtr());
            //TH1D* m1030 = static_cast<TH1D*>(v1030_m.at(key).GetPtr());
            //TH1D* m3050 = static_cast<TH1D*>(v3050_m.at(key).GetPtr());

            TH1D* s010  = &v010_s.at(key);
            TH1D* s1030 = &v1030_s.at(key);
            TH1D* s3050 = &v3050_s.at(key);
            TH1D* m010  = &v010_m.at(key);
            TH1D* m1030 = &v1030_m.at(key);
            TH1D* m3050 = &v3050_m.at(key);

            TH1D* mn010  = NormaliseHistogram(s010,  m010,  NORM_LOW, NORM_HIGH);
            TH1D* mn1030 = NormaliseHistogram(s1030, m1030, NORM_LOW, NORM_HIGH);
            TH1D* mn3050 = NormaliseHistogram(s3050, m3050, NORM_LOW, NORM_HIGH);

            iter_data.push_back({
                {{"010", s010}, {"1030", s1030}, {"3050", s3050}},
                {{"010", m010}, {"1030", m1030}, {"3050", m3050}},
                {{"010", mn010},{"1030", mn1030},{"3050", mn3050}}
            });
        }

        for (int iter = 0; iter < N_ITERATIONS; ++iter) {
            TDirectory* outDirIter = outDir->mkdir(("iter_" + std::to_string(iter)).c_str());
            outDirIter->cd();
            for (const auto& cent : {"010", "1030", "3050"}) {
                iter_data[iter].sames.at(cent)->Write(("hSame_"   + std::string(cent)).c_str());
                iter_data[iter].mixeds.at(cent)->Write(("hMixed_" + std::string(cent)).c_str());
                iter_data[iter].mixeds_norm.at(cent)->Write(
                    ("hMixedNormalised_" + std::string(cent)).c_str());
            }
            for (auto& [k,v] : iter_data[iter].mixeds_norm) delete v;
        }
    }

    outFile->Close();
    delete outFile;
    delete chain_same;
    delete chain_mixed;
}

void UpperLimitSystematicRoutine()
{
    const int N_ITERATIONS             = 1;
    const int N_UPPER_LIMIT_ITERATIONS = 1000;

    TFile* inFile  = TFile::Open("output/hist_systematics_with_upper_limit.root");
    TFile* outFile = TFile::Open("output/systematics_with_upper_limit.root", "RECREATE");

    for (const std::string& sign : {"Matter", "Antimatter"}) {

        TDirectory* inDir  = dynamic_cast<TDirectory*>(inFile->Get(sign.c_str()));
        TDirectory* outDir = outFile->mkdir(sign.c_str());

        TH1D h_upper_limits_010("hUpperLimits_010",
            ";[#it{N}_{^{4}Li}^{raw}]_{upper limit};", 600, 0, 600);

        for (int iter = 0; iter < N_ITERATIONS; ++iter) {

            std::cout << "Processing " << sign << " iter " << iter+1
                      << "/" << N_ITERATIONS << "\n";

            TDirectory* outDirIter = outDir->mkdir(("iter_" + std::to_string(iter)).c_str());

            TH1D* h_same_010 = dynamic_cast<TH1D*>(
                inDir->Get(("iter_" + std::to_string(iter) + "/hSame_010").c_str()));
            TH1D* h_mixed_010 = dynamic_cast<TH1D*>(
                inDir->Get(("iter_" + std::to_string(iter) + "/hMixed_010").c_str()));
            TH1D* h_mixed_norm_010 = dynamic_cast<TH1D*>(
                inDir->Get(("iter_" + std::to_string(iter) + "/hMixedNormalised_010").c_str()));

            TH1D* h_corr_010 = static_cast<TH1D*>(
                h_same_010->Clone(("hCorrelation_010_" + std::to_string(iter)).c_str()));
            h_corr_010->Divide(h_mixed_norm_010);

            TDirectory* outDirNominal = outDirIter->mkdir("nominal");
            double nominal_raw_yield = FittingRoutine(
                outDirNominal,
                "/home/galucia/Lithium4/femto/models/lambda_models.root",
                sign + "/010/hLambdaSigmaCorrectedCk",
                h_corr_010, h_mixed_norm_010,
                sign, "010", /*use_smoothening=*/true);

            TH1D h_raw_yields_010(
                ("hRawYields_010_Iter_" + std::to_string(iter)).c_str(),
                ";Raw yield;", 300, -200, 400);

            for (int iter_upper = 0; iter_upper < N_UPPER_LIMIT_ITERATIONS; ++iter_upper) {

                TH1D* h_same_up = static_cast<TH1D*>(
                    h_same_010->Clone(("hSame_010_Iter_Upper_"
                        + std::to_string(iter_upper)).c_str()));
                TH1D* h_mixed_norm_up = static_cast<TH1D*>(
                    h_mixed_norm_010->Clone(("hMixed_010_Iter_Upper_"
                        + std::to_string(iter_upper)).c_str()));

                PoissonSampling(h_same_up);
                PoissonSampling(h_mixed_norm_up);

                TH1D* h_corr_up = static_cast<TH1D*>(
                    h_same_up->Clone(("hCorrelation_010_Iter_Upper_"
                        + std::to_string(iter_upper)).c_str()));
                h_corr_up->Divide(h_mixed_norm_up);

                TDirectory* outDir_to_pass = (iter_upper == 0)
                    ? outDirIter->mkdir("inner_iter_0")
                    : nullptr;

                double raw_yield = FittingRoutine(
                    outDir_to_pass,
                    "/home/galucia/Lithium4/femto/models/lambda_models.root",
                    sign + "/010/hLambdaSigmaCorrectedCk",
                    h_corr_up, h_mixed_norm_up,
                    sign, "010", /*use_smoothening=*/true);

                h_raw_yields_010.Fill(raw_yield);

                delete h_same_up;
                delete h_mixed_norm_up;
                delete h_corr_up;
            }

            TF1 fit_func(("fit_func_010_" + std::to_string(iter)).c_str(),
                         "gaus", -200, 400);
            fit_func.SetParameters(h_raw_yields_010.GetMaximum(),
                                   h_raw_yields_010.GetMean(),
                                   h_raw_yields_010.GetRMS());
            h_raw_yields_010.Fit(&fit_func, "RMS+", "", -200, 400);

            outDirIter->cd();
            h_raw_yields_010.Write();

            // 95% confidence level upper limit
            double upper_limit_010 = nominal_raw_yield + 1.96 * fit_func.GetParameter(2);
            h_upper_limits_010.Fill(upper_limit_010);

            delete h_same_010;
            delete h_mixed_010;
            delete h_mixed_norm_010;
            delete h_corr_010;
        }

        outDir->cd();
        h_upper_limits_010.Write();
    }

    inFile->Close();
    outFile->Close();
    delete inFile;
    delete outFile;
}

int cutOptimization()
{
    RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL); // level 5

    PrepareHistograms();  // uncomment to (re)generate input histograms
    //UpperLimitSystematicRoutine();

    std::cout << "Systematics analysis completed successfully.\n";
    return 0;
}