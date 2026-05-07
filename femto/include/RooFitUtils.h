/**
 * RooFitUtils.h
 * 
 * Utilities for RooFit model preparation and fitting
 */

#ifndef ROOFIT_UTILS_H
#define ROOFIT_UTILS_H

#include <TH1F.h>
#include <TFile.h>
#include <TCanvas.h>
#include <iostream>

#include <RooRealVar.h>
#include <RooDataHist.h>
#include <RooHistPdf.h>
#include <RooAddPdf.h>
#include <RooPlot.h>
#include <RooFitResult.h>

namespace FemtoAnalysis {

/**
 * Load and prepare signal model from histogram
 * @param kstar RooRealVar for k*
 * @param signalData Output RooDataHist
 * @param signalPdf Output RooHistPdf
 * @param signalFilePath Path to file containing signal histogram
 * @param signalHistName Name of signal histogram in file
 * @param outfile Output file for diagnostic plots
 * @return true if successful
 */
inline bool prepareSignalModel(RooRealVar& kstar, 
                              RooDataHist*& signalData, 
                              RooHistPdf*& signalPdf,
                              const char* signalFilePath,
                              const char* signalHistName,
                              TFile* outfile = nullptr) {
    using namespace RooFit;
    
    TFile* file = TFile::Open(signalFilePath);
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Cannot open signal file " << signalFilePath << std::endl;
        return false;
    }
    
    TH1F* hSignal = (TH1F*)file->Get(signalHistName);
    if (!hSignal) {
        std::cerr << "Error: Cannot load signal histogram " << signalHistName << std::endl;
        file->Close();
        return false;
    }
    hSignal->SetDirectory(0);
    file->Close();
    
    signalData = new RooDataHist("signal_data", "signal_data", kstar, Import(*hSignal));
    signalPdf = new RooHistPdf("signal_pdf", "signal_pdf", kstar, *signalData);
    
    // Save diagnostic plot if outfile provided
    if (outfile) {
        RooPlot* frame = kstar.frame();
        signalData->plotOn(frame, MarkerStyle(20), MarkerSize(0.5));
        signalPdf->plotOn(frame, LineColor(kBlue), LineWidth(2));
        frame->SetTitle(";#it{k}* (GeV/#it{c});C(#it{k}*)");
        
        TCanvas canvas("canvas_signal_fit", "Signal Fit", 800, 600);
        canvas.cd();
        frame->Draw();
        
        outfile->mkdir("signal_fit");
        outfile->cd("signal_fit");
        hSignal->Write("h_correlation_signal");
        canvas.Write();
        
        delete frame;
    }
    
    delete hSignal;
    return true;
}

/**
 * Load and prepare background model from histogram
 * @param kstar RooRealVar for k*
 * @param backgroundData Output RooDataHist
 * @param backgroundPdf Output RooHistPdf
 * @param backgroundFilePath Path to file containing background histogram
 * @param backgroundHistName Name of background histogram in file
 * @param outfile Output file for diagnostic plots
 * @return true if successful
 */
inline bool prepareBackgroundModel(RooRealVar& kstar,
                                  RooDataHist*& backgroundData,
                                  RooHistPdf*& backgroundPdf,
                                  const char* backgroundFilePath,
                                  const char* backgroundHistName,
                                  TFile* outfile = nullptr) {
    using namespace RooFit;
    
    TFile* file = TFile::Open(backgroundFilePath);
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Cannot open background file " << backgroundFilePath << std::endl;
        return false;
    }
    
    TH1F* hBackground = (TH1F*)file->Get(backgroundHistName);
    if (!hBackground) {
        std::cerr << "Error: Cannot load background histogram " << backgroundHistName << std::endl;
        file->Close();
        return false;
    }
    hBackground->SetDirectory(0);
    file->Close();
    
    backgroundData = new RooDataHist("background_data", "background_data", kstar, Import(*hBackground));
    backgroundPdf = new RooHistPdf("background_pdf", "background_pdf", kstar, *backgroundData);
    
    // Save diagnostic plot if outfile provided
    if (outfile) {
        RooPlot* frame = kstar.frame();
        backgroundData->plotOn(frame, MarkerStyle(20), MarkerSize(0.5));
        backgroundPdf->plotOn(frame, LineColor(kRed), LineWidth(2));
        frame->SetTitle(";#it{k}* (GeV/#it{c});C(#it{k}*)");
        
        TCanvas canvas("canvas_background_fit", "Background Fit", 800, 600);
        canvas.cd();
        frame->Draw();
        
        outfile->mkdir("background_fit");
        outfile->cd("background_fit");
        hBackground->Write("h_correlation_background");
        canvas.Write();
        
        delete frame;
    }
    
    delete hBackground;
    return true;
}

/**
 * Perform fit and return chi2
 * @param kstar RooRealVar for k*
 * @param data RooDataHist to fit
 * @param model RooAddPdf model
 * @param fitRange Fit range string (e.g., "fit_range")
 * @param outfile Output file for saving results
 * @param fitName Name for this fit
 * @param doDrawing Whether to create diagnostic plots
 * @return chi2/ndf value
 */
inline double performFit(RooRealVar& kstar, 
                        RooDataHist* data, 
                        RooAddPdf* model,
                        const char* fitRange,
                        TFile* outfile = nullptr,
                        const char* fitName = "fit_result",
                        bool doDrawing = false) {
    using namespace RooFit;
    
    RooFitResult* fitResult = model->fitTo(*data, Range(fitRange), Save(), PrintLevel(-1));
    
    double chi2 = 0.0;
    if (fitResult) {
        RooPlot* frame = kstar.frame(Range(fitRange));
        data->plotOn(frame, MarkerStyle(20), MarkerSize(0.5));
        model->plotOn(frame, LineColor(kBlue), LineWidth(2));
        
        chi2 = frame->chiSquare();
        
        if (doDrawing && outfile) {
            TCanvas canvas(Form("canvas_%s", fitName), "Fit Result", 800, 600);
            canvas.cd();
            frame->Draw();
            
            outfile->cd("fit_results");
            canvas.Write(Form("canvas_%s", fitName));
        }
        
        delete frame;
        delete fitResult;
    }
    
    return chi2;
}

} // namespace FemtoAnalysis

#endif // ROOFIT_UTILS_H
