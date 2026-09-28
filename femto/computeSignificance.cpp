/**
 * This code computes the chi2 disitrbution for the null hypothesis 
*/

#include <iostream>
#include <vector>
#include <deque>
#include <numeric>

#include <TGraphAsymmErrors.h>
#include <TFile.h>
#include <TH1F.h>
#include <TGraph.h>
#include <TRandom3.h>
#include <TMath.h>

namespace InputData {
    const char * inputMixedFile = "/home/galucia/Lithium4/preparation/output/PbPb/correlation_PbPb_hadronpid.root";
    const char * inputCorrectionFile = "models/lambda_models_LL_10.root";
    //const char * inputCorrectionFile = "models/lambda_models_coulomb_10.root";
    const char * inputChi2File = "/home/galucia/Lithium4/femto/output/PbPb_fit_correlation_function_hadronpid__smoothened_finer_binning_LL_10.root";
    const char * outputFile = "/home/galucia/Lithium4/femto/output/significance_LL_010_1050.root";

    const std::vector<std::string> bkgVariationHistNames = {
    //};
        "nominal/hLambdaSigmaCorrectedCk_Smeared_nominal",
        "nominal/hLambdaSigmaCorrectedCk_Smeared_nominal_higher",
        "nominal/hLambdaSigmaCorrectedCk_Smeared_nominal_lower",
        "upper/hLambdaSigmaCorrectedCk_Smeared_upper",
        "upper/hLambdaSigmaCorrectedCk_Smeared_upper_higher",
        "upper/hLambdaSigmaCorrectedCk_Smeared_upper_lower",
        "lower/hLambdaSigmaCorrectedCk_Smeared_lower",
        "lower/hLambdaSigmaCorrectedCk_Smeared_lower_higher",
        "lower/hLambdaSigmaCorrectedCk_Smeared_lower_lower   ",
    };
};

namespace {
    const float kstarMin = 0.0;
    const float kstarMax = 0.4;
    
    const int N_ITERATIONS = 1000000; // 1 million
    const int N_BINS = 40; // 40 bins (this has to match the binning of the input mixed event histogram)
    const int N_BINS_WINDOW = 4; // 6 bins for the window

    const int NBINS_CHI2 = 1600;
    const float CHI2_MAX_VALUE = 160;
};

enum class MatterMode { Matter, Antimatter, Both };

struct CentralityConfig {
    const char* name;
    bool directComputation;
};

const std::vector<CentralityConfig> kCentralities = {
    {"010",  false},
    //{"1030", false},
    //{"3050", false},
    //{"5080", false},    
    //{"050",  true },
    //{"080",  true },
    {"1050", true },
    //{"1080", true },
};

enum class BkgVariation { Nominal, Low, High };

const char* bkgVariationSuffix(BkgVariation var) {
    switch (var) {
        case BkgVariation::Low:  return "_bkg_low";
        case BkgVariation::High: return "_bkg_high";
        default: return "";
    }
}

TH1F* buildBkgEnvelope(TFile* fileCorrection, const std::string& baseDir, TH1F* hRef,
                       const std::vector<std::string>& variationNames, bool useMax) {

    TH1F* hEnvelope = (TH1F*)hRef->Clone(useMax ? "hBkgEnvelopeHigh" : "hBkgEnvelopeLow");
    hEnvelope->Reset();
    hEnvelope->SetDirectory(0);

    std::vector<TH1F*> variations;
    for (const auto& name : variationNames) {
        std::string fullName = baseDir + "/" + name;
        auto h = (TH1F*)fileCorrection->Get(fullName.c_str());
        if (!h) {
            std::cerr << "Warning: could not find variation histogram " << fullName << std::endl;
            continue;
        }
        variations.push_back(h);
    }

    for (int ibin = 1; ibin <= hEnvelope->GetNbinsX(); ++ibin) {
        double kstar = hEnvelope->GetBinCenter(ibin);
        double extremum = useMax ? -1e300 : 1e300;
        for (auto h : variations) {
            double value = h->Interpolate(kstar); // handles differing binning, as in build_bkg_envelope
            extremum = useMax ? std::max(extremum, value) : std::min(extremum, value);
        }
        hEnvelope->SetBinContent(ibin, extremum);
    }

    return hEnvelope;
}

void loadSameMixedSingle(TH1F *& hSame, TH1F *& hMixed, const char* centrality, const bool directComputation, 
                         const BkgVariation bkgVariation, const bool isMatter = false) {

    TFile *fileMixed = TFile::Open(InputData::inputMixedFile);
    std::string mixedHistName = std::string(isMatter ? "CorrelationMatter/Default/hMixedEvent" 
                                                     : "CorrelationAntimatter/Default/hMixedEvent");
    if (directComputation) mixedHistName += "DirectComputation";
    mixedHistName += centrality;
    hMixed = (TH1F*)fileMixed->Get(mixedHistName.c_str());
    
    TFile *fileCorrection = TFile::Open(InputData::inputCorrectionFile);
    std::string baseDir = std::string(isMatter ? "Matter/" : "Antimatter/") + centrality;
    std::string nominalName = baseDir + "/hLambdaSigmaCorrectedCk";
    auto hCorrectionNominal = (TH1F*)fileCorrection->Get(nominalName.c_str());

    TH1F* hCorrection = nullptr;
    if (bkgVariation == BkgVariation::Nominal) {
        hCorrection = (TH1F*)hCorrectionNominal->Clone("hCorrectionTmp");
    } else {
        hCorrection = buildBkgEnvelope(fileCorrection, baseDir, hCorrectionNominal,
                                        InputData::bkgVariationHistNames, bkgVariation == BkgVariation::High);
    }
    hCorrection->SetDirectory(0);
    fileCorrection->Close();

    hSame = (TH1F*)hMixed->Clone("hSame");
    hSame->SetDirectory(0);

    for (int ibin = 1; ibin <= hSame->GetNbinsX()+1; ++ibin) {
        double mixedValue = hMixed->GetBinContent(ibin);
        double binCenter = hMixed->GetBinCenter(ibin);
        double correctionValue = hCorrection->GetBinContent(hCorrection->FindBin(binCenter));
        hSame->SetBinContent(ibin, mixedValue * correctionValue);
    }

    delete hCorrection;
}

void loadSameMixed(TH1F *& hSame, TH1F *& hMixed, const char* centrality, const bool directComputation, 
                   const BkgVariation bkgVariation, const MatterMode matterMode) {

    TH1F *hSameMatter = nullptr, *hMixedMatter = nullptr;
    TH1F *hSameAnti  = nullptr, *hMixedAnti  = nullptr;

    if (matterMode == MatterMode::Matter || matterMode == MatterMode::Both) {
        loadSameMixedSingle(hSameMatter, hMixedMatter, centrality, directComputation, bkgVariation, true);
    }
    if (matterMode == MatterMode::Antimatter || matterMode == MatterMode::Both) {
        loadSameMixedSingle(hSameAnti, hMixedAnti, centrality, directComputation, bkgVariation, false);
    }

    if (matterMode == MatterMode::Both) {
        hSame  = (TH1F*)hSameMatter->Clone("hSame");
        hMixed = (TH1F*)hMixedMatter->Clone("hMixed");
        hSame->Add(hSameAnti);
        hMixed->Add(hMixedAnti);
        hSame->SetDirectory(0);
        hMixed->SetDirectory(0);
        delete hSameMatter; delete hMixedMatter;
        delete hSameAnti;   delete hMixedAnti;
    } else {
        hSame  = (matterMode == MatterMode::Matter) ? hSameMatter  : hSameAnti;
        hMixed = (matterMode == MatterMode::Matter) ? hMixedMatter : hMixedAnti;
    }
}

void computeCorrelationFunction(TH1F *hSame, TH1F *hMixed, TH1F* hCorrelation) {

    for  (int ibin = 1; ibin <= hSame->GetNbinsX()+1; ++ibin) {
        double sameValue = hSame->GetBinContent(ibin);
        double mixedValue = hMixed->GetBinContent(ibin);
        double sameError = std::sqrt(sameValue);
        double mixedError = std::sqrt(mixedValue);

        if (mixedValue > 0) {
            double correlationValue = sameValue / mixedValue;
            hCorrelation->SetBinContent(ibin, correlationValue);
            hCorrelation->SetBinError(ibin, correlationValue*std::sqrt((sameError/sameValue)*(sameError/sameValue) + (mixedError/mixedValue)*(mixedError/mixedValue)) );
        } else {
            hCorrelation->SetBinContent(ibin, 0.0);
        }
    }
}

void runMcChi2(TH1F *& hChi2, TH1F *& hChi2FarFromSignal,
               std::vector<TH1F *>& runningChi2Histograms,
               std::vector<TH1F *>& windowChi2Histograms,
               float kstarBinCenters[],
               TDirectory * outfile,
               const int N_BINS = 40 /* 40 bins */,
               const int N_ITERATIONS = 1000000 /* 1 mln */,
               const char* centrality = "010",
               const bool directComputation = false,
               const MatterMode matterMode = MatterMode::Antimatter,
               const BkgVariation bkgVariation = BkgVariation::Nominal) {

    TH1F* hSame, * hMixed;
    loadSameMixed(hSame, hMixed, centrality, directComputation, bkgVariation, matterMode);
    auto hCorrelation = (TH1F*)hSame->Clone("hCorrelation");
    std::cout << "Cloned histogram for correlation function." << std::endl;
    computeCorrelationFunction(hSame, hMixed, hCorrelation);

    std::cout << "Loaded histograms: " << hSame->GetName() << " and " << hMixed->GetName() << std::endl;

    auto hSameIter = (TH1F*)hSame->Clone("hSameIter");
    auto hMixedIter = (TH1F*)hMixed->Clone("hMixedIter");
    auto hCorrelationIter = (TH1F*)hSameIter->Clone("hCorrelationIter");

    const int N_WINDOW_BINS = runningChi2Histograms.size() - windowChi2Histograms.size();
    std::deque<float> chi2Deque;
    chi2Deque.resize(N_WINDOW_BINS);
    std::vector<std::vector<uint64_t>> runningChi2Counts(N_BINS, std::vector<uint64_t>(runningChi2Histograms[0]->GetNbinsX(), 0));

    for (int iter = 0; iter < N_ITERATIONS; ++iter) {
        if (iter % static_cast<int>(N_ITERATIONS / 100) == 0) {
            std::cout << "Processing iteration: " << iter << "/" << N_ITERATIONS << std::endl;
        }
        hSameIter->Reset();
        hMixedIter->Reset();

        for (int ibin = 1; ibin <= N_BINS; ++ibin) {

            hSameIter->SetBinContent(ibin, gRandom->Poisson(hSame->GetBinContent(ibin)));
            hMixedIter->SetBinContent(ibin, gRandom->Poisson(hMixed->GetBinContent(ibin)));
        }
        computeCorrelationFunction(hSameIter, hMixedIter, hCorrelationIter);

        double chi2 = 0.0, kstar = 0.0, expected = 0.0, observed = 0.0, error = 0.0;
        double chi2Cumulated = 0.0;
        double chi2Limited = 0;
        double chi2FarFromSignal = 0;
        //const int FIRST_BIN = hCorrelation->FindBin(0.01);
        const int FIRST_BIN = hCorrelation->FindBin(kstarMin);

        for (int ibin = FIRST_BIN; ibin <= N_BINS; ++ibin) {
            kstar = hSameIter->GetBinCenter(ibin);
            expected = hCorrelationIter->GetBinContent(ibin);  
            observed = hCorrelation->GetBinContent(ibin);
            error = hCorrelationIter->GetBinError(ibin);

            if (error > 0) {
                chi2 = (observed - expected) * (observed - expected) / (error * error);
            }
            chi2Cumulated += chi2;

            if (kstar < 0.15) {
                chi2Limited = chi2Cumulated;
                //chi2Limited = chi2Cumulated / (ibin+1); // reduced chi2
            } else {
                chi2FarFromSignal += chi2;
            }
            
            //runningChi2Histograms[ibin-1]->Fill(chi2Cumulated);
            int idx = std::min(NBINS_CHI2 - 1, std::max(0, (int)(chi2Cumulated / CHI2_MAX_VALUE * NBINS_CHI2)));
            runningChi2Counts[ibin-1][idx]++;
            
            chi2Deque.pop_front();
            chi2Deque.push_back(chi2);
            if (ibin > N_WINDOW_BINS) {
                double chi2Window = std::accumulate(chi2Deque.begin(), chi2Deque.end(), 0., std::plus<double>());
                windowChi2Histograms[ibin - N_WINDOW_BINS - 1]->Fill(chi2Window);
            }
            //runningChi2Histograms[ibin-1]->Fill(chi2 / (ibin+1)); // reduced chi2

            if (iter == 0) {
                kstarBinCenters[ibin-1] = hCorrelation->GetBinCenter(ibin);
            }
        }

        hChi2->Fill(chi2Limited);
        hChi2FarFromSignal->Fill(chi2FarFromSignal);
    }

    for (int ibin = 0; ibin < N_BINS; ++ibin) {
        for (int idx = 0; idx < NBINS_CHI2; ++idx) {
            runningChi2Histograms[ibin]->SetBinContent(idx+1, runningChi2Counts[ibin][idx]);
        }
    }

    outfile->cd();
    hSame->Write();
    hMixed->Write();
    hCorrelation->Write();
    hSameIter->Write();
    hMixedIter->Write();
    hCorrelationIter->Write();

    delete hSame;
    delete hMixed;
    delete hCorrelation;
    delete hSameIter;
    delete hMixedIter;
    delete hCorrelationIter;

}

void displayRunningResult(TH1F *& hChi2, std::vector<TH1F *> & runningChi2Histograms,
                          TDirectory * outfile, float kstarBinCenters[],
                          const char * centrality,
                          const MatterMode matterMode,
                          std::vector<double>& outSignificance,
                          std::vector<double>& outPvalue,
                          const char * suffix = "",
                          const BkgVariation bkgVariation = BkgVariation::Nominal,
                          const int N_BINS = 40 /* 40 bins */,
                          const int N_ITERATIONS = 1000000 /* 1 mln */) {


    auto infile = TFile::Open(InputData::inputChi2File);
    std::string matterLabel = (matterMode == MatterMode::Matter) ? "Matter" : (matterMode == MatterMode::Antimatter) ? "Antimatter" : "";
    std::string chi2Name = std::string(matterLabel) + "" + centrality + "/model/chi2" + bkgVariationSuffix(bkgVariation) + "_stat_only";
    auto hChi2Data = (TH1F*)infile->Get(chi2Name.c_str());
    
    std::vector<double> runningChi2(hChi2Data->GetNbinsX());
    for (size_t ichi2 = 0; ichi2 < runningChi2.size(); ichi2++) {
        runningChi2[ichi2] = hChi2Data->GetBinContent(ichi2+1);
    }

    outfile->cd();
    hChi2->Write();
    hChi2Data->Write(Form("hChi2Data%s", suffix));

    outfile->mkdir(Form("runningChi2%s", suffix));
    outfile->cd(Form("runningChi2%s", suffix));

    TGraph *gRunningPvalue = new TGraph(N_BINS);
    gRunningPvalue->SetTitle("Running P-value;#it{k}* (GeV/#it{c});P-value");
    TGraph *gRunningSignificance = new TGraph(N_BINS);
    gRunningSignificance->SetTitle("Running Significance;#it{k}* (GeV/#it{c});Significance");

    outSignificance.assign(N_BINS, 0.0);
    outPvalue.assign(N_BINS, 0.0);

    for (int ibin = 0; ibin < N_BINS; ++ibin) {
        
        const float chi2Value = runningChi2[ibin];
        const float kstar = kstarBinCenters[ibin];
        
        TF1 *fGamma = new TF1(Form("fGamma_%d_%s", ibin, suffix), "[0]*ROOT::Math::gamma_pdf(x,[1],[2],0)", 2, 200);
        fGamma->SetParameters(runningChi2Histograms[ibin]->Integral("width"), ibin/2., 0.5); // rough starting guess: shape~ndof, scale~2
        runningChi2Histograms[ibin]->Fit(fGamma, "RLQ"); // "L" = log-likelihood fit, much better behaved for tails than chi2 fit

        double k     = fGamma->GetParameter(1);
        double theta = fGamma->GetParameter(2);
        //double pvalue = 1.0 - ROOT::Math::gamma_cdf(chi2Value, k, theta);
        //double significance = TMath::NormQuantile(1. - pvalue/2.);
        double pvalue = ROOT::Math::gamma_cdf_c(chi2Value, k, theta);
        double significance = std::sqrt(2.0) * TMath::ErfcInverse(pvalue);
        std::cout << "Bin " << ibin << ", k* = " << kstar << ", chi2 = " << chi2Value << ", p-value = " << pvalue << ", significance = " << significance << std::endl;

        runningChi2Histograms[ibin]->Write();
        // //const float chi2Value = runningChi2[ibin] / (ibin+1); // reduced chi2
        // const float pvalue = runningChi2Histograms[ibin]->Integral(runningChi2Histograms[ibin]->FindBin(chi2Value), runningChi2Histograms[ibin]->GetNbinsX()+1) / N_ITERATIONS;
        // const float significance = TMath::NormQuantile(1. - pvalue/2.);
        
        gRunningPvalue->SetPoint(ibin, kstar, pvalue);
        gRunningSignificance->SetPoint(ibin, kstar, significance);

        outSignificance[ibin] = significance;
        outPvalue[ibin] = pvalue;

        //delete fGamma;
    }
    gRunningPvalue->SetMarkerStyle(20);
    gRunningPvalue->Write(Form("gRunningPvalue%s", suffix));
    gRunningSignificance->SetMarkerStyle(20);
    gRunningSignificance->SetMaximum(10.0);
    gRunningSignificance->Write(Form("gRunningSignificance%s", suffix));

    delete hChi2Data;
    infile->Close();

}

void computeSignificance() {

    auto outfile = TFile::Open(InputData::outputFile, "RECREATE");
    auto hChi2 = new TH1F("hChi2", "Chi2 Distribution;#chi^{2};Counts", NBINS_CHI2, 0, CHI2_MAX_VALUE);
    auto hChi2FarFromSignal = new TH1F("hChi2FarFromSignal", "Chi2 Distribution (far from signal);#chi^{2};Counts", NBINS_CHI2, 0, CHI2_MAX_VALUE);
    
    std::vector<TH1F*> runningChi2Histograms, windowChi2Histograms;
    runningChi2Histograms.reserve(N_BINS);
    windowChi2Histograms.reserve(N_BINS - N_BINS_WINDOW);
    for (int ibin = 0; ibin < N_BINS; ++ibin) {
        std::string name_running = "hRunningChi2_" + std::to_string(ibin);
        auto hRunningChi2 = new TH1F(name_running.c_str(), Form("Running Chi2 %d ;#chi^{2};Counts", ibin), NBINS_CHI2, 0, CHI2_MAX_VALUE);
        runningChi2Histograms.emplace_back(hRunningChi2);

        if (ibin > N_BINS_WINDOW / 2 && ibin <= N_BINS - (N_BINS_WINDOW / 2)) {
            std::string nameWindow = "hWindowChi2_" + std::to_string(ibin - N_BINS_WINDOW);
            auto hWindowChi2 = new TH1F(nameWindow.c_str(), Form("Window Chi2 %d ;#chi^{2};Counts", ibin - N_BINS_WINDOW), NBINS_CHI2, 0, CHI2_MAX_VALUE);
            windowChi2Histograms.emplace_back(hWindowChi2);
        }
    }
    float kstarBinCenters[N_BINS];

    for (const auto& cent : kCentralities) {
        auto runForMode = [&](MatterMode mode) {
            const char* label = (mode == MatterMode::Matter) ? "Matter" 
                            : (mode == MatterMode::Antimatter) ? "Antimatter" 
                            : "Both";
            
            std::vector<double> sigNominal, sigLow, sigHigh, pvalNominal, pvalLow, pvalHigh;
            for (auto bkgVariation : {BkgVariation::Nominal, BkgVariation::Low, BkgVariation::High}) {
                hChi2->Reset();
                hChi2FarFromSignal->Reset();
                for (auto& h : runningChi2Histograms) h->Reset();
                for (auto& h : windowChi2Histograms)  h->Reset();

                auto outdir = outfile->mkdir(Form("%s%s%s", label, cent.name, bkgVariationSuffix(bkgVariation)));
                runMcChi2(hChi2, hChi2FarFromSignal, runningChi2Histograms, windowChi2Histograms,
                        kstarBinCenters, outdir, N_BINS, N_ITERATIONS, cent.name, cent.directComputation, mode, bkgVariation);
                
                std::vector<double> sig, pval;
                displayRunningResult(hChi2, runningChi2Histograms, outdir, kstarBinCenters,
                                    cent.name, mode, sig, pval,
                                    Form("%s%s%s", label, cent.name, bkgVariationSuffix(bkgVariation)),
                                    bkgVariation, N_BINS, N_ITERATIONS);

                if (bkgVariation == BkgVariation::Nominal) { sigNominal = sig; pvalNominal = pval; }
                else if (bkgVariation == BkgVariation::Low)  { sigLow  = sig; pvalLow  = pval; }
                else if (bkgVariation == BkgVariation::High) { sigHigh = sig; pvalHigh = pval; }
            }

            auto outdirNominal = outfile->GetDirectory(Form("%s%s", label, cent.name));
            outdirNominal->cd();

            TGraphAsymmErrors *gRunningSignificanceWithUnc = new TGraphAsymmErrors(N_BINS);
            gRunningSignificanceWithUnc->SetTitle("Running Significance (with bkg variation uncertainty);#it{k}* (GeV/#it{c});Significance");
            gRunningSignificanceWithUnc->SetMarkerStyle(20);

            for (int ibin = 0; ibin < N_BINS; ++ibin) {
                double central = sigNominal[ibin];
                double a = sigLow[ibin] - central;
                double b = sigHigh[ibin] - central;
                double errLow  = std::max(0.0, -std::min(a, b)); // lower error: how far below nominal the smaller of low/high falls
                double errHigh = std::max(0.0,  std::max(a, b)); // upper error: how far above nominal the larger of low/high falls

                gRunningSignificanceWithUnc->SetPoint(ibin, kstarBinCenters[ibin], central);
                gRunningSignificanceWithUnc->SetPointError(ibin, 0., 0., errLow, errHigh);
            }
            gRunningSignificanceWithUnc->Write(Form("gRunningSignificanceWithUnc%s%s", label, cent.name));
        };

        runForMode(MatterMode::Antimatter);
        runForMode(MatterMode::Matter);
        runForMode(MatterMode::Both);
    }

    outfile->Close();

}
