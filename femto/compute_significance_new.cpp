/**
 * computeSignificance.cpp
 * 
 * This code computes the chi2 distribution for the null hypothesis 
 * and estimates the statistical significance
 */

#include <iostream>
#include <vector>
#include <TFile.h>
#include <TH1F.h>
#include <TGraph.h>
#include <TRandom3.h>

// Include unified framework headers
#include "FemtoConfig.h"
#include "HistogramUtils.h"
#include "Chi2Utils.h"

using namespace FemtoAnalysis;

// Configuration specific to significance calculation
namespace SignificanceConfig {
    const int N_ITERATIONS = 1000000; // 1 million
    const int N_BINS = 40; // 40 bins (must match input histogram binning)
    const int N_BINS_WINDOW = 4; // bins for the sliding window
    
    // Input file paths
    InputDataPaths paths = {
        .inputMixedFile = "/home/galucia/Lithium4/preparation/checks/correlation_all_pass1_pass4_nclstpc.root",
        .inputMixedNameMatter = "CorrelationAntimatter/hMixedEvent050",
        .inputMixedNameAntimatter = "CorrelationAntimatter/hMixedEvent050",
        .inputCorrectionFile = "models/lambda_models.root",
        .inputCorrectionNameAntimatter = "Antimatter/050/hLambdaCorrectedCk",
        .inputCorrectionNameMatter = "Matter/050/hLambdaCorrectedCk",
        .inputChi2File = "/home/galucia/Lithium4/femto/output/fit_correlation_function_all_pass1_pass4_nclstpc.root",
        .inputChi2NameAntimatter = "Antimatter050/model/chi2",
        .inputChi2NameMatter = "Matter050/model/chi2"
    };
}

/**
 * Run Monte Carlo chi2 calculation
 */
void runMcChi2(TH1F*& hChi2, TH1F*& hChi2FarFromSignal,
               std::vector<TH1F*>& runningChi2Histograms,
               std::vector<TH1F*>& windowChi2Histograms,
               float kstarBinCenters[],
               TFile* outfile,
               const bool isMatter = false) {

    using namespace SignificanceConfig;
    
    // Load histograms
    TH1F *hSame, *hMixed;
    const char* mixedPath = isMatter ? paths.inputMixedNameMatter : paths.inputMixedNameAntimatter;
    const char* correctionPath = isMatter ? paths.inputCorrectionNameMatter : paths.inputCorrectionNameAntimatter;
    
    if (!loadSameMixedWithCorrection(hSame, hMixed, paths.inputMixedFile, mixedPath,
                                     paths.inputCorrectionFile, correctionPath)) {
        std::cerr << "Failed to load histograms for " << (isMatter ? "Matter" : "Antimatter") << std::endl;
        return;
    }
    
    std::cout << "Loaded histograms: " << hSame->GetName() << " and " << hMixed->GetName() << std::endl;
    
    // Compute observed correlation function
    auto hCorrelation = (TH1F*)hSame->Clone("hCorrelation");
    computeCorrelationFunction(hSame, hMixed, hCorrelation);
    
    // Create iteration histograms
    auto hSameIter = (TH1F*)hSame->Clone("hSameIter");
    auto hMixedIter = (TH1F*)hMixed->Clone("hMixedIter");
    auto hCorrelationIter = (TH1F*)hSameIter->Clone("hCorrelationIter");
    
    // Initialize window chi2 tracking
    const int N_WINDOW_BINS = N_BINS_WINDOW;
    std::deque<float> chi2Deque(N_WINDOW_BINS, 0.0);
    
    // Main MC loop
    const int FIRST_BIN = hCorrelation->FindBin(0.0);
    
    for (int iter = 0; iter < N_ITERATIONS; ++iter) {
        if (iter % 10000 == 0) {
            std::cout << "Processing iteration: " << iter << "/" << N_ITERATIONS << std::endl;
        }
        
        // Sample histograms using Poisson distribution
        sampleHistogramPoisson(hSame, hSameIter);
        sampleHistogramPoisson(hMixed, hMixedIter);
        computeCorrelationFunction(hSameIter, hMixedIter, hCorrelationIter);
        
        // Compute chi2 values
        double chi2Cumulated = 0.0;
        double chi2Limited = 0.0;
        double chi2FarFromSignal = 0.0;
        
        for (int ibin = FIRST_BIN; ibin <= hSameIter->GetNbinsX(); ++ibin) {
            double kstar = hSameIter->GetBinCenter(ibin);
            double observed = hCorrelation->GetBinContent(ibin);
            double expected = hCorrelationIter->GetBinContent(ibin);
            double error = hCorrelationIter->GetBinError(ibin);
            
            double chi2 = 0.0;
            if (error > 0) {
                chi2 = std::pow((observed - expected) / error, 2);
            }
            chi2Cumulated += chi2;
            
            // Track chi2 in signal region (k* < 0.15)
            if (kstar < 0.15) {
                chi2Limited = chi2Cumulated;
            } else {
                chi2FarFromSignal += chi2;
            }
            
            // Fill running chi2 histogram
            if (ibin - FIRST_BIN < runningChi2Histograms.size()) {
                runningChi2Histograms[ibin - FIRST_BIN]->Fill(chi2Cumulated);
            }
            
            // Fill window chi2 histograms
            chi2Deque.pop_front();
            chi2Deque.push_back(chi2);
            if (ibin > FIRST_BIN + N_WINDOW_BINS - 1) {
                double chi2Window = std::accumulate(chi2Deque.begin(), chi2Deque.end(), 0.0);
                int windowIdx = ibin - FIRST_BIN - N_WINDOW_BINS;
                if (windowIdx >= 0 && windowIdx < windowChi2Histograms.size()) {
                    windowChi2Histograms[windowIdx]->Fill(chi2Window);
                }
            }
            
            // Store k* bin centers on first iteration
            if (iter == 0 && ibin - FIRST_BIN < N_BINS) {
                kstarBinCenters[ibin - FIRST_BIN] = hCorrelation->GetBinCenter(ibin);
            }
        }
        
        hChi2->Fill(chi2Limited);
        hChi2FarFromSignal->Fill(chi2FarFromSignal);
    }
    
    // Save histograms
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

/**
 * Display running chi2 results and compute significance
 */
void displayRunningResult(TH1F*& hChi2, 
                         std::vector<TH1F*>& runningChi2Histograms,
                         TFile* outfile, 
                         float kstarBinCenters[],
                         const char* inputChi2Name,
                         const char* suffix = "") {

    using namespace SignificanceConfig;
    
    // Load observed chi2 data
    auto infile = TFile::Open(paths.inputChi2File);
    if (!infile || infile->IsZombie()) {
        std::cerr << "Error: Cannot open chi2 file" << std::endl;
        return;
    }
    
    auto hChi2Data = (TH1F*)infile->Get(inputChi2Name);
    if (!hChi2Data) {
        std::cerr << "Error: Cannot load chi2 data histogram" << std::endl;
        infile->Close();
        return;
    }
    
    std::vector<double> runningChi2(hChi2Data->GetNbinsX());
    for (size_t ichi2 = 0; ichi2 < runningChi2.size(); ichi2++) {
        runningChi2[ichi2] = hChi2Data->GetBinContent(ichi2 + 1);
    }
    
    // Save overall chi2 distributions
    outfile->cd();
    hChi2->Write();
    hChi2Data->Write(Form("hChi2Data%s", suffix));
    
    // Create directory for running chi2 results
    outfile->mkdir(Form("runningChi2%s", suffix));
    outfile->cd(Form("runningChi2%s", suffix));
    
    // Create significance graphs
    TGraph *gRunningPvalue, *gRunningSignificance;
    createSignificanceGraphs(runningChi2Histograms, runningChi2, kstarBinCenters,
                            N_ITERATIONS, gRunningPvalue, gRunningSignificance);
    
    // Save individual chi2 histograms and graphs
    for (size_t ibin = 0; ibin < runningChi2Histograms.size(); ++ibin) {
        runningChi2Histograms[ibin]->Write();
    }
    
    gRunningPvalue->Write(Form("gRunningPvalue%s", suffix));
    gRunningSignificance->Write(Form("gRunningSignificance%s", suffix));
    
    delete hChi2Data;
    delete gRunningPvalue;
    delete gRunningSignificance;
    infile->Close();
}

/**
 * Main function
 */
void computeSignificance() {

    using namespace SignificanceConfig;
    
    auto outfile = TFile::Open("output/significance.root", "RECREATE");
    
    // Create chi2 distribution histograms
    auto hChi2 = new TH1F("hChi2", "Chi2 Distribution;#chi^{2};Counts", 
                         Config::NBINS_CHI2, 0, Config::CHI2_MAX_VALUE);
    auto hChi2FarFromSignal = new TH1F("hChi2FarFromSignal", 
                                       "Chi2 Distribution (far from signal);#chi^{2};Counts", 
                                       Config::NBINS_CHI2, 0, Config::CHI2_MAX_VALUE);
    
    // Create running and window chi2 histograms
    std::vector<TH1F*> runningChi2Histograms, windowChi2Histograms;
    runningChi2Histograms.reserve(N_BINS);
    windowChi2Histograms.reserve(N_BINS - N_BINS_WINDOW);
    
    for (int ibin = 0; ibin < N_BINS; ++ibin) {
        std::string name_running = "hRunningChi2_" + std::to_string(ibin);
        auto hRunningChi2 = new TH1F(name_running.c_str(), 
                                     Form("Running Chi2 %d ;#chi^{2};Counts", ibin), 
                                     Config::NBINS_CHI2, 0, Config::CHI2_MAX_VALUE);
        runningChi2Histograms.emplace_back(hRunningChi2);
        
        if (ibin > N_BINS_WINDOW / 2 && ibin <= N_BINS - (N_BINS_WINDOW / 2)) {
            std::string nameWindow = "hWindowChi2_" + std::to_string(ibin - N_BINS_WINDOW);
            auto hWindowChi2 = new TH1F(nameWindow.c_str(), 
                                       Form("Window Chi2 %d ;#chi^{2};Counts", ibin - N_BINS_WINDOW), 
                                       Config::NBINS_CHI2, 0, Config::CHI2_MAX_VALUE);
            windowChi2Histograms.emplace_back(hWindowChi2);
        }
    }
    
    float kstarBinCenters[N_BINS];
    
    // Process Antimatter
    std::cout << "\n=== Processing Antimatter ===" << std::endl;
    runMcChi2(hChi2, hChi2FarFromSignal, runningChi2Histograms, windowChi2Histograms,
              kstarBinCenters, outfile, false /*isMatter*/);
    displayRunningResult(hChi2, runningChi2Histograms, outfile, kstarBinCenters,
                        paths.inputChi2NameAntimatter, "Antimatter");
    
    // Reset histograms for Matter
    hChi2->Reset();
    hChi2FarFromSignal->Reset();
    for (auto& hist : runningChi2Histograms) {
        hist->Reset();
    }
    for (auto& hist : windowChi2Histograms) {
        hist->Reset();
    }
    
    // Process Matter
    std::cout << "\n=== Processing Matter ===" << std::endl;
    runMcChi2(hChi2, hChi2FarFromSignal, runningChi2Histograms, windowChi2Histograms,
              kstarBinCenters, outfile, true /*isMatter*/);
    displayRunningResult(hChi2, runningChi2Histograms, outfile, kstarBinCenters,
                        paths.inputChi2NameMatter, "Matter");
    
    outfile->Close();
    
    std::cout << "\n=== Analysis complete. Results saved to output/significance.root ===" << std::endl;
}
