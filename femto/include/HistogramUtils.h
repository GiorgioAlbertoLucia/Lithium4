/**
 * HistogramUtils.h
 * 
 * Common histogram manipulation utilities for femtoscopy analysis
 */

#ifndef HISTOGRAM_UTILS_H
#define HISTOGRAM_UTILS_H

#include <TH1F.h>
#include <TFile.h>
#include <TMath.h>
#include <iostream>

namespace FemtoAnalysis {

/**
 * Load same-event and mixed-event histograms from file
 * @param hSame Output histogram for same events
 * @param hMixed Output histogram for mixed events
 * @param filePath Path to the input file
 * @param sameEventPath Path to same-event histogram within file
 * @param mixedEventPath Path to mixed-event histogram within file
 * @return true if successful, false otherwise
 */
inline bool loadHistograms(TH1F*& hSame, TH1F*& hMixed, 
                          const char* filePath,
                          const char* sameEventPath,
                          const char* mixedEventPath) {
    TFile* file = TFile::Open(filePath);
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Cannot open file " << filePath << std::endl;
        return false;
    }
    
    hSame = (TH1F*)file->Get(sameEventPath);
    hMixed = (TH1F*)file->Get(mixedEventPath);
    
    if (!hSame || !hMixed) {
        std::cerr << "Error: Cannot load histograms from file" << std::endl;
        file->Close();
        return false;
    }
    
    hSame->SetDirectory(0);
    hMixed->SetDirectory(0);
    file->Close();
    
    return true;
}

/**
 * Load same-event histogram and compute it from mixed-event with correction
 * Used in significance calculation where same = mixed * correction
 * @param hSame Output histogram for same events (computed)
 * @param hMixed Output histogram for mixed events
 * @param mixedFilePath Path to the mixed-event file
 * @param mixedHistPath Path to mixed histogram within file
 * @param correctionFilePath Path to the correction file
 * @param correctionHistPath Path to correction histogram within file
 * @return true if successful, false otherwise
 */
inline bool loadSameMixedWithCorrection(TH1F*& hSame, TH1F*& hMixed,
                                       const char* mixedFilePath,
                                       const char* mixedHistPath,
                                       const char* correctionFilePath,
                                       const char* correctionHistPath) {
    TFile* fileMixed = TFile::Open(mixedFilePath);
    if (!fileMixed || fileMixed->IsZombie()) {
        std::cerr << "Error: Cannot open mixed event file" << std::endl;
        return false;
    }
    
    hMixed = (TH1F*)fileMixed->Get(mixedHistPath);
    if (!hMixed) {
        std::cerr << "Error: Cannot load mixed histogram" << std::endl;
        fileMixed->Close();
        return false;
    }
    hMixed->SetDirectory(0);
    fileMixed->Close();
    
    TFile* fileCorrection = TFile::Open(correctionFilePath);
    if (!fileCorrection || fileCorrection->IsZombie()) {
        std::cerr << "Error: Cannot open correction file" << std::endl;
        return false;
    }
    
    TH1F* hCorrection = (TH1F*)fileCorrection->Get(correctionHistPath);
    if (!hCorrection) {
        std::cerr << "Error: Cannot load correction histogram" << std::endl;
        fileCorrection->Close();
        return false;
    }
    hCorrection->SetDirectory(0);
    fileCorrection->Close();
    
    // Compute same = mixed * correction
    hSame = (TH1F*)hMixed->Clone("hSame");
    hSame->SetDirectory(0);
    
    for (int ibin = 1; ibin <= hSame->GetNbinsX() + 1; ++ibin) {
        double mixedValue = hMixed->GetBinContent(ibin);
        double binCenter = hMixed->GetBinCenter(ibin);
        double correctionValue = hCorrection->GetBinContent(hCorrection->FindBin(binCenter));
        hSame->SetBinContent(ibin, mixedValue * correctionValue);
    }
    
    delete hCorrection;
    return true;
}

/**
 * Compute correlation function C(k*) = Same/Mixed
 * Includes proper error propagation
 * @param hSame Input same-event histogram
 * @param hMixed Input mixed-event histogram
 * @param hCorrelation Output correlation function histogram
 */
inline void computeCorrelationFunction(const TH1F* hSame, const TH1F* hMixed, TH1F* hCorrelation) {
    const int nbins = hSame->GetNbinsX();
    
    for (int ibin = 1; ibin <= nbins + 1; ++ibin) {
        const double sameValue = hSame->GetBinContent(ibin);
        const double mixedValue = hMixed->GetBinContent(ibin);
        
        if (mixedValue <= 0) {
            hCorrelation->SetBinContent(ibin, 1e-12);
            hCorrelation->SetBinError(ibin, 1e-12);
            continue;
        }
        
        const double ratio = sameValue / mixedValue;
        const double sameError = std::sqrt(sameValue);
        const double mixedError = std::sqrt(mixedValue);
        const double relError = std::sqrt(
            std::pow(sameError / sameValue, 2) + 
            std::pow(mixedError / mixedValue, 2)
        );
        
        hCorrelation->SetBinContent(ibin, ratio);
        hCorrelation->SetBinError(ibin, ratio * relError);
    }
}

/**
 * Sample histogram using Poisson distribution for each bin
 * @param hInput Input histogram to sample
 * @param hOutput Output histogram with sampled values
 */
inline void sampleHistogramPoisson(const TH1F* hInput, TH1F* hOutput) {
    const int nbins = hInput->GetNbinsX();
    for (int ibin = 1; ibin <= nbins; ++ibin) {
        double value = hInput->GetBinContent(ibin);
        hOutput->SetBinContent(ibin, gRandom->Poisson(value));
    }
}

/**
 * Sample histogram using Gaussian distribution for each bin
 * @param hInput Input histogram to sample
 * @param hOutput Output histogram with sampled values
 */
inline void sampleHistogramGaussian(const TH1F* hInput, TH1F* hOutput) {
    const int nbins = hInput->GetNbinsX();
    for (int ibin = 1; ibin <= nbins; ++ibin) {
        double centralValue = hInput->GetBinContent(ibin);
        double error = hInput->GetBinError(ibin);
        double sampledValue = gRandom->Gaus(centralValue, error);
        hOutput->SetBinContent(ibin, sampledValue);
        hOutput->SetBinError(ibin, error);
    }
}

} // namespace FemtoAnalysis

#endif // HISTOGRAM_UTILS_H
