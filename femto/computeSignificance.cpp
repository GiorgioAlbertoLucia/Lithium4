/**
 * This code computes the chi2 disitrbution for the null hypothesis 
*/

#include <iostream>
#include <vector>
#include <deque>
#include <numeric>
#include <TFile.h>
#include <TH1F.h>
#include <TGraph.h>
#include <TRandom3.h>
#include <TMath.h>

struct InputData {
    const char * inputMixedFile = "/home/galucia/Lithium4/preparation/checks/correlation_all_pass1_pass4_nclstpc.root";
    const char * inputMixedName = "CorrelationAntimatter/hMixedEvent050";
    const char * inputCorrectionFile = "models/CATS_converted.root";
    const char * inputCorrectionName = "hHe3_p_Coul_CF";
    const char * inputChi2File = "/home/galucia/Lithium4/femto/output/fit_correlation_function_all_pass1_pass4_nclstpc.root";
    const char * inputChi2NameAntimatter = "Antimatter050/model/chi2";
    const char * inputChi2NameMatter = "Matter050/model/chi2";
};

void loadSameMixed(TH1F *& hSame, TH1F *& hMixed) {

    InputData inputData;
    TFile *fileMixed = TFile::Open(inputData.inputMixedFile);
    hMixed = (TH1F*)fileMixed->Get(inputData.inputMixedName);
    hMixed->SetDirectory(0);
    fileMixed->Close();

    TFile *fileCorrection = TFile::Open(inputData.inputCorrectionFile);
    auto hCorrection = (TH1F*)fileCorrection->Get(inputData.inputCorrectionName);
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
               TFile * outfile,
               const int N_BINS = 40 /* 40 bins */,
               const int N_ITERATIONS = 1000000 /* 1 mln */) {

    TH1F* hSame, * hMixed;
    loadSameMixed(hSame, hMixed);
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

    for (int iter = 0; iter < N_ITERATIONS; ++iter) {
        if (iter % 10000 == 0) {
            std::cout << "Processing iteration: " << iter << "/" << N_ITERATIONS << std::endl;
        }
        hSameIter->Reset();
        hMixedIter->Reset();

        for (int ibin = 1; ibin <= hSame->GetNbinsX(); ++ibin) {

            hSameIter->SetBinContent(ibin, gRandom->Poisson(hSame->GetBinContent(ibin)));
            hMixedIter->SetBinContent(ibin, gRandom->Poisson(hMixed->GetBinContent(ibin)));
        }
        computeCorrelationFunction(hSameIter, hMixedIter, hCorrelationIter);

        double chi2 = 0.0, kstar = 0.0, expected = 0.0, observed = 0.0, error = 0.0;
        double chi2Cumulated = 0.0;
        double chi2Limited = 0;
        double chi2FarFromSignal = 0;
        //const int FIRST_BIN = hCorrelation->FindBin(0.01);
        const int FIRST_BIN = hCorrelation->FindBin(0.0);

        for (int ibin = FIRST_BIN; ibin <= hSameIter->GetNbinsX(); ++ibin) {
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
            runningChi2Histograms[ibin-1]->Fill(chi2Cumulated);
            
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
                          TFile * outfile, float kstarBinCenters[],
                          const char * inputChi2Name,
                          const char * suffix = "",
                          const int N_BINS = 40 /* 40 bins */,
                          const int N_ITERATIONS = 1000000 /* 1 mln */) {


    InputData inputData;
    auto infile = TFile::Open(inputData.inputChi2File);
    auto hChi2Data = (TH1F*)infile->Get(inputChi2Name);
    
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

    for (int ibin = 0; ibin < N_BINS; ++ibin) {
        runningChi2Histograms[ibin]->Write();
        const float kstar = kstarBinCenters[ibin];
        const float chi2Value = runningChi2[ibin];
        //const float chi2Value = runningChi2[ibin] / (ibin+1); // reduced chi2
        const float pvalue = runningChi2Histograms[ibin]->Integral(runningChi2Histograms[ibin]->FindBin(chi2Value), runningChi2Histograms[ibin]->GetNbinsX()+1) / N_ITERATIONS;
        const float significance = TMath::NormQuantile(1. - pvalue/2.);

        gRunningPvalue->SetPoint(ibin, kstar, pvalue);
        gRunningSignificance->SetPoint(ibin, kstar, significance);
    }
    gRunningPvalue->SetMarkerStyle(20);
    gRunningPvalue->Write(Form("gRunningPvalue%s", suffix));
    gRunningSignificance->SetMarkerStyle(20);
    gRunningSignificance->Write(Form("gRunningSignificance%s", suffix));

    delete hChi2Data;
    infile->Close();

}

void computeSignificance() {

    const int N_ITERATIONS = 1000000; // 1 million
    const int N_BINS = 40; // 40 bins (this has to match the binning of the input mixed event histogram)
    const int N_BINS_WINDOW = 4; // 6 bins for the window

    const int NBINS_CHI2 = 1600;
    const float CHI2_MAX_VALUE = 160;

    auto outfile = TFile::Open("output/significance.root", "RECREATE");
    
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
    runMcChi2(hChi2, hChi2FarFromSignal, runningChi2Histograms, windowChi2Histograms,
                kstarBinCenters, outfile, N_BINS, N_ITERATIONS);

    InputData inputData;
    displayRunningResult(hChi2, runningChi2Histograms, outfile, kstarBinCenters,
                        inputData.inputChi2NameAntimatter, "Antimatter",
                        N_BINS, N_ITERATIONS);
    displayRunningResult(hChi2, runningChi2Histograms, outfile, kstarBinCenters,
                        inputData.inputChi2NameMatter, "Matter",
                        N_BINS, N_ITERATIONS);

    outfile->Close();

}