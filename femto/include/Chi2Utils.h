/**
 * Chi2Utils.h
 * 
 * Utilities for chi2 calculation and significance estimation
 */

#ifndef CHI2_UTILS_H
#define CHI2_UTILS_H

#include <TH1F.h>
#include <TGraph.h>
#include <TFile.h>
#include <TMath.h>
#include <vector>
#include <deque>
#include <numeric>
#include <iostream>

namespace FemtoAnalysis {

/**
 * Compute chi2 between observed and expected correlation functions
 * @param hObserved Observed correlation function (data)
 * @param hExpected Expected correlation function (from toy MC)
 * @param firstBin First bin to include in chi2
 * @param lastBin Last bin to include in chi2
 * @return chi2 value
 */
inline double computeChi2(const TH1F* hObserved, const TH1F* hExpected, 
                         int firstBin = 1, int lastBin = -1) {
    if (lastBin < 0) lastBin = hObserved->GetNbinsX();
    
    double chi2 = 0.0;
    for (int ibin = firstBin; ibin <= lastBin; ++ibin) {
        double observed = hObserved->GetBinContent(ibin);
        double expected = hExpected->GetBinContent(ibin);
        double error = hExpected->GetBinError(ibin);
        
        if (error > 0) {
            chi2 += std::pow((observed - expected) / error, 2);
        }
    }
    return chi2;
}

/**
 * Compute running chi2 for each bin (cumulative from first bin)
 * Fills histogram array with chi2 distributions
 * @param hObserved Observed correlation function
 * @param hExpected Expected correlation function (from toy MC)
 * @param runningChi2Histograms Vector of histograms to fill
 * @param firstBin First bin to start chi2 calculation
 */
inline void fillRunningChi2(const TH1F* hObserved, const TH1F* hExpected,
                           std::vector<TH1F*>& runningChi2Histograms,
                           int firstBin = 1) {
    double chi2Cumulated = 0.0;
    
    for (int ibin = firstBin; ibin <= hExpected->GetNbinsX(); ++ibin) {
        double observed = hObserved->GetBinContent(ibin);
        double expected = hExpected->GetBinContent(ibin);
        double error = hExpected->GetBinError(ibin);
        
        double chi2 = 0.0;
        if (error > 0) {
            chi2 = std::pow((observed - expected) / error, 2);
        }
        chi2Cumulated += chi2;
        
        if (ibin - firstBin < runningChi2Histograms.size()) {
            runningChi2Histograms[ibin - firstBin]->Fill(chi2Cumulated);
        }
    }
}

/**
 * Compute windowed chi2 (chi2 over a sliding window of bins)
 * @param hObserved Observed correlation function
 * @param hExpected Expected correlation function (from toy MC)
 * @param windowChi2Histograms Vector of histograms to fill
 * @param windowSize Number of bins in the window
 * @param firstBin First bin to start calculation
 */
inline void fillWindowChi2(const TH1F* hObserved, const TH1F* hExpected,
                          std::vector<TH1F*>& windowChi2Histograms,
                          int windowSize, int firstBin = 1) {
    std::deque<float> chi2Deque;
    chi2Deque.resize(windowSize, 0.0);
    
    for (int ibin = firstBin; ibin <= hExpected->GetNbinsX(); ++ibin) {
        double observed = hObserved->GetBinContent(ibin);
        double expected = hExpected->GetBinContent(ibin);
        double error = hExpected->GetBinError(ibin);
        
        double chi2 = 0.0;
        if (error > 0) {
            chi2 = std::pow((observed - expected) / error, 2);
        }
        
        chi2Deque.pop_front();
        chi2Deque.push_back(chi2);
        
        if (ibin > firstBin + windowSize - 1) {
            double chi2Window = std::accumulate(chi2Deque.begin(), chi2Deque.end(), 0.0);
            int windowIdx = ibin - firstBin - windowSize;
            if (windowIdx >= 0 && windowIdx < windowChi2Histograms.size()) {
                windowChi2Histograms[windowIdx]->Fill(chi2Window);
            }
        }
    }
}

/**
 * Create graphs of p-value and significance vs k*
 * @param runningChi2Histograms Vector of chi2 distribution histograms
 * @param observedChi2 Vector of observed chi2 values (one per bin)
 * @param kstarBinCenters Array of k* bin centers
 * @param nIterations Number of MC iterations used
 * @param gPvalue Output TGraph for p-values
 * @param gSignificance Output TGraph for significance
 */
inline void createSignificanceGraphs(const std::vector<TH1F*>& runningChi2Histograms,
                                    const std::vector<double>& observedChi2,
                                    const float* kstarBinCenters,
                                    int nIterations,
                                    TGraph*& gPvalue,
                                    TGraph*& gSignificance) {
    const int nBins = runningChi2Histograms.size();
    gPvalue = new TGraph(nBins);
    gSignificance = new TGraph(nBins);
    
    gPvalue->SetTitle("Running P-value;#it{k}* (GeV/#it{c});P-value");
    gSignificance->SetTitle("Running Significance;#it{k}* (GeV/#it{c});Significance");
    
    for (int ibin = 0; ibin < nBins; ++ibin) {
        const float kstar = kstarBinCenters[ibin];
        const float chi2Value = observedChi2[ibin];
        
        // Calculate p-value from chi2 distribution
        int chi2Bin = runningChi2Histograms[ibin]->FindBin(chi2Value);
        int lastBin = runningChi2Histograms[ibin]->GetNbinsX() + 1;
        float pvalue = runningChi2Histograms[ibin]->Integral(chi2Bin, lastBin) / nIterations;
        
        // Convert p-value to significance
        float significance = TMath::NormQuantile(1.0 - pvalue / 2.0);
        
        gPvalue->SetPoint(ibin, kstar, pvalue);
        gSignificance->SetPoint(ibin, kstar, significance);
    }
    
    gPvalue->SetMarkerStyle(20);
    gSignificance->SetMarkerStyle(20);
}

} // namespace FemtoAnalysis

#endif // CHI2_UTILS_H
