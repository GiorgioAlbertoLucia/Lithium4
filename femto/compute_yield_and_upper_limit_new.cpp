/**
 * compute_yield_and_upper_limit.cpp
 * 
 * This code computes the chi2 distribution for the null hypothesis 
 * and estimates the yield uncertainty for Li4 from correlation function fits
 */

#include <iostream>
#include <unordered_map>

#include <TFile.h>
#include <TH1F.h>
#include <TF1.h>
#include <TRandom3.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TString.h>

#include <RooRealVar.h>
#include <RooDataHist.h>
#include <RooHistPdf.h>
#include <RooAddPdf.h>
#include <RooPlot.h>
#include <RooFitResult.h>

// Include unified framework headers
#include "FemtoConfig.h"
#include "HistogramUtils.h"
#include "Chi2Utils.h"
#include "RooFitUtils.h"

using namespace FemtoAnalysis;

// Configuration specific to yield calculation
namespace YieldConfig {
    const char* CORRELATION_FILE = "/home/galucia/Lithium4/preparation/checks/correlation_hadronpid_pass1_pass4_nohe3pcut_offlinetpc.root";
    const char* SIGNAL_FILE = "/home/galucia/Lithium4/femto/output/sampling_check.root";
    const char* BACKGROUND_FILE = "/home/galucia/Lithium4/femto/models/lambda_models.root";
    
    std::unordered_map<std::string, const char*> OUTPUT_FILE = {
        {"Antimatter", "/home/galucia/Lithium4/femto/output/yield_upper_limit_Antimatter.root"},
        {"Matter", "/home/galucia/Lithium4/femto/output/yield_upper_limit_Matter.root"}
    };
    
    const int N_ITERATIONS = 1000;
    const int PRINT_INTERVAL = 100;
    
    SamplingMethod SAMPLING_METHOD = SamplingMethod::POISSONIAN;
}

// Forward declarations
bool load_same_mixed(TH1F*& h_same, TH1F*& h_mixed, const char* sign, const char* centrality);
void prefit_background(RooRealVar& kstar, TH1F* h_data, RooDataHist*& data, RooAddPdf*& model, 
                      RooRealVar& nbkg, TFile* outfile, double& bkg_normalisation_at_reference, 
                      double& bkg_full_correction, const char* background_name, bool do_drawing);
double fit(RooRealVar& kstar, RooDataHist*& data, RooAddPdf*& model, TFile* outfile, 
          double& sig_normalisation_at_reference, const char* fit_name, bool do_drawing);
double compute_yield(RooAbsPdf* model_pdf, RooRealVar& nbkg, RooRealVar& nsig, RooRealVar& xvar, 
                    TH1F* h_mixed_event, TFile* outfile, double bkg_normalisation_at_reference, 
                    double sig_normalisation_at_reference, const char* yield_name, bool do_drawing);
void save_high_yield_fit(RooRealVar& kstar, RooDataHist*& data, RooAddPdf*& model, int iter, TFile* outfile);
void setup_output_directories(TFile* outfile);
IterationResult process_single_iteration(RooRealVar& kstar, RooAddPdf* model, RooRealVar& nsig, 
                                        RooRealVar& nbkg, TH1F*& h_same_iter, TH1F*& h_mixed_iter, 
                                        TH1F*& h_correlation_iter, const TH1F* h_same, const TH1F* h_mixed, 
                                        const TH1F* h_correlation, TFile* outfile, int iter, 
                                        bool do_drawing, SamplingMethod sampling_method);

/**
 * Load same and mixed event histograms
 */
bool load_same_mixed(TH1F*& h_same, TH1F*& h_mixed, const char* sign, const char* centrality) {
    using namespace YieldConfig;
    
    TString sameEventPath = Form("Correlation%s/Default/hSameEvent%s", sign, centrality);
    TString mixedEventPath = Form("Correlation%s/Default/hNormalisedMixedEvent%s", sign, centrality);
    
    bool success = loadHistograms(h_same, h_mixed, CORRELATION_FILE, 
                                  sameEventPath.Data(), mixedEventPath.Data());
    
    if (success) {
        h_same->Print();
        h_mixed->Print();
    }
    
    return success;
}

/**
 * Pre-fit background in the flat region
 */
void prefit_background(RooRealVar& kstar, TH1F* h_data, RooDataHist*& data, RooAddPdf*& model, 
                      RooRealVar& nbkg, TFile* outfile, double& bkg_normalisation_at_reference, 
                      double& bkg_full_correction, const char* background_name = "background_pdf", 
                      bool do_drawing = false) {
    using namespace RooFit;
    
    // Set range for pre-fit (flat region)
    kstar.setRange("prefit_range", Config::PREFIT_MIN, Config::PREFIT_MAX);
    
    // Create RooDataHist from histogram
    data = new RooDataHist(background_name, background_name, kstar, Import(*h_data));
    
    // Fit in the flat region
    RooFitResult* fitResult = model->fitTo(*data, Range("prefit_range"), Save(), PrintLevel(-1));
    
    // Calculate normalization factors
    RooAbsPdf* bkg_pdf = model->pdfList().at(1);
    int ref_bin = h_data->FindBin(Config::REFERENCE_KSTAR_FOR_BKG_NORMALIZATION);
    bkg_normalisation_at_reference = bkg_pdf->getVal(RooArgSet(kstar)) * nbkg.getVal();
    bkg_full_correction = h_data->GetBinContent(ref_bin) / bkg_normalisation_at_reference;
    
    // Save diagnostic plot if requested
    if (do_drawing) {
        RooPlot* frame = kstar.frame(Range("prefit_range"));
        data->plotOn(frame, MarkerStyle(20), MarkerSize(0.5));
        model->plotOn(frame, LineColor(kRed), LineWidth(2));
        frame->SetTitle(Form("Background Pre-fit;#it{k}* (GeV/#it{c});C(#it{k}*)"));
        
        TCanvas canvas(Form("canvas_%s", background_name), "Background Pre-fit", 800, 600);
        canvas.cd();
        frame->Draw();
        
        outfile->cd("prefit_background");
        canvas.Write(Form("canvas_%s", background_name));
        
        delete frame;
    }
    
    delete fitResult;
}

/**
 * Perform full fit with signal and background
 */
double fit(RooRealVar& kstar, RooDataHist*& data, RooAddPdf*& model, TFile* outfile, 
          double& sig_normalisation_at_reference, const char* fit_name = "fit_result", 
          bool do_drawing = false) {
    using namespace RooFit;
    
    // Set full fit range
    kstar.setRange("fit_range", Config::KSTAR_MIN, Config::KSTAR_MAX);
    
    // Perform fit
    RooFitResult* fitResult = model->fitTo(*data, Range("fit_range"), Save(), PrintLevel(-1));
    
    // Calculate chi2
    RooPlot* frame = kstar.frame(Range("fit_range"));
    data->plotOn(frame, MarkerStyle(20), MarkerSize(0.5));
    model->plotOn(frame, LineColor(kBlue), LineWidth(2));
    
    double chi2 = frame->chiSquare();
    
    // Calculate signal normalization
    RooAbsPdf* sig_pdf = model->pdfList().at(0);
    RooRealVar* nsig = (RooRealVar*)&model->coefList().at(0);
    kstar.setVal(Config::REFERENCE_KSTAR_FOR_SIG_NORMALIZATION);
    sig_normalisation_at_reference = sig_pdf->getVal(RooArgSet(kstar)) * nsig->getVal();
    
    // Save diagnostic plot if requested
    if (do_drawing) {
        TCanvas canvas(Form("canvas_%s", fit_name), "Fit Result", 800, 600);
        canvas.cd();
        frame->Draw();
        
        outfile->cd("fit_results");
        canvas.Write(Form("canvas_%s", fit_name));
    }
    
    delete frame;
    delete fitResult;
    
    return chi2;
}

/**
 * Compute yield from fitted model
 */
double compute_yield(RooAbsPdf* model_pdf, RooRealVar& nsig, RooRealVar& nbkg, RooRealVar& xvar, 
                    TH1F* h_mixed_event, TFile* outfile, double bkg_normalisation_at_reference, 
                    double sig_normalisation_at_reference, const char* yield_name = "yield_result", 
                    bool do_drawing = false) {
    
    // Store original values
    double nsig_stored = nsig.getVal();
    double nbkg_stored = nbkg.getVal();
    
    // Get PDFs
    RooAbsPdf* sig_pdf = ((RooAddPdf*)model_pdf)->pdfList().at(0);
    RooAbsPdf* bkg_pdf = ((RooAddPdf*)model_pdf)->pdfList().at(1);
    
    // Create histograms for signal and background
    auto h_signal_correlation = (TH1F*)h_mixed_event->Clone(Form("%s_sig", yield_name));
    auto h_background_correlation = (TH1F*)h_mixed_event->Clone(Form("%s_bkg", yield_name));
    auto h_same_event_signal = (TH1F*)h_mixed_event->Clone(yield_name);
    
    double yield = 0.0;
    
    // Fill histograms bin by bin
    for (int ibin = 1; ibin <= h_mixed_event->GetNbinsX(); ++ibin) {
        double kstar_val = h_mixed_event->GetBinCenter(ibin);
        xvar.setVal(kstar_val);
        
        // Signal correlation
        double sig_ck = sig_pdf->getVal(RooArgSet(xvar)) * nsig.getVal() / sig_normalisation_at_reference;
        h_signal_correlation->SetBinContent(ibin, sig_ck);
        
        // Background correlation  
        double bkg_ck = bkg_pdf->getVal(RooArgSet(xvar)) * nbkg.getVal() / bkg_normalisation_at_reference;
        h_background_correlation->SetBinContent(ibin, bkg_ck);
        
        // Same event yield
        double mixed_val = h_mixed_event->GetBinContent(ibin);
        double same_val = mixed_val * (sig_ck + bkg_ck);
        h_same_event_signal->SetBinContent(ibin, same_val);
        
        yield += same_val;
    }
    
    // Save results
    outfile->cd("yield_results");
    h_same_event_signal->Write(yield_name);
    if (do_drawing) {
        h_background_correlation->Write(Form("%s_bkg_correlation", yield_name));
    }
    
    // Restore original values
    nsig.setVal(nsig_stored);
    nbkg.setVal(nbkg_stored);
    
    delete h_signal_correlation;
    delete h_same_event_signal;
    delete h_background_correlation;
    
    return yield;
}

/**
 * Save high yield fit for inspection
 */
void save_high_yield_fit(RooRealVar& kstar, RooDataHist*& data, RooAddPdf*& model, 
                        int iter, TFile* outfile) {
    using namespace RooFit;
    
    RooPlot* frame = kstar.frame();
    data->plotOn(frame, MarkerStyle(20), MarkerSize(0.5));
    model->plotOn(frame, LineColor(kBlue), LineWidth(2));
    model->paramOn(frame, Layout(0.45, 0.85, 0.35), Format("NEU", AutoPrecision(2)));
    frame->SetTitle(Form("High Yield Fit Result - Iteration %d;#it{k}* (GeV/#it{c});C(#it{k}*)", iter));
    
    TCanvas canvas(Form("canvas_high_yield_%d", iter), "High Yield Fit Result", 800, 600);
    canvas.cd();
    frame->Draw();
    outfile->cd("high_yield_results");
    canvas.Write();
    
    delete frame;
}

/**
 * Setup output directory structure
 */
void setup_output_directories(TFile* outfile) {
    outfile->mkdir("prefit_background");
    outfile->mkdir("fit_results");
    outfile->mkdir("yield_results");
    outfile->mkdir("high_yield_results");
}

/**
 * Process a single MC iteration
 */
IterationResult process_single_iteration(RooRealVar& kstar, RooAddPdf* model, RooRealVar& nsig, 
                                        RooRealVar& nbkg, TH1F*& h_same_iter, TH1F*& h_mixed_iter, 
                                        TH1F*& h_correlation_iter, const TH1F* h_same, 
                                        const TH1F* h_mixed, const TH1F* h_correlation, 
                                        TFile* outfile, int iter, bool do_drawing,
                                        SamplingMethod sampling_method) {
    
    const int nbins = h_same->GetNbinsX();
    
    // Sample histograms based on method
    if (sampling_method == SamplingMethod::POISSONIAN) {
        sampleHistogramPoisson(h_same, h_same_iter);
        sampleHistogramPoisson(h_mixed, h_mixed_iter);
        computeCorrelationFunction(h_same_iter, h_mixed_iter, h_correlation_iter);
    } 
    else if (sampling_method == SamplingMethod::GAUSSIAN) {
        sampleHistogramGaussian(h_correlation, h_correlation_iter);
        for (int ibin = 1; ibin <= nbins; ++ibin) {
            h_mixed_iter->SetBinContent(ibin, h_mixed->GetBinContent(ibin));
        }
    }
    
    // Pre-fit background
    RooDataHist* data_prefit = nullptr;
    nsig.setVal(0);
    nsig.setConstant(true);
    double bkg_normalisation_at_reference = 1.0;
    double bkg_full_correction = 1.0;
    prefit_background(kstar, h_correlation_iter, data_prefit, model, nbkg, outfile, 
                     bkg_normalisation_at_reference, bkg_full_correction,
                     Form("background_pdf_%d", iter), do_drawing);
    delete data_prefit;
    
    // Full fit with signal
    nsig.setConstant(false);
    nbkg.setConstant(true);
    double sig_normalisation_at_reference = 1.0;
    
    RooDataHist* data_iter = new RooDataHist("data_iter", "data_iter", kstar, 
                                            RooFit::Import(*h_correlation_iter));
    
    double chi2 = fit(kstar, data_iter, model, outfile, sig_normalisation_at_reference,
                     Form("fit_result_%d", iter), do_drawing);
    double yield = compute_yield(model, nsig, nbkg, kstar, h_mixed_iter, outfile, 
                                bkg_normalisation_at_reference, sig_normalisation_at_reference,
                                Form("yield_result_%d", iter), do_drawing);
    
    // Save high yield fits for inspection
    if (yield > Config::HIGH_YIELD_THRESHOLD) {
        save_high_yield_fit(kstar, data_iter, model, iter, outfile);
    }
    
    delete data_iter;
    
    return {yield, chi2};
}

/**
 * Main function
 */
void compute_yield_and_upper_limit() {
    
    using namespace YieldConfig;
    
    const char* sign = "Antimatter";
    const char* centrality = "010";
    
    // Load data
    TH1F *h_same, *h_mixed;
    if (!load_same_mixed(h_same, h_mixed, sign, centrality)) {
        std::cerr << "Failed to load histograms. Exiting." << std::endl;
        return;
    }
    
    // Compute correlation function
    auto h_correlation = (TH1F*)h_same->Clone("h_correlation");
    computeCorrelationFunction(h_same, h_mixed, h_correlation);
    std::cout << "Loaded histograms: " << h_same->GetName() << " and " << h_mixed->GetName() << std::endl;
    
    // Create iteration histograms
    auto h_same_iter = (TH1F*)h_same->Clone("h_same_iter");
    auto h_mixed_iter = (TH1F*)h_mixed->Clone("h_mixed_iter");
    auto h_correlation_iter = (TH1F*)h_same_iter->Clone("h_correlation_iter");
    std::cout << "Cloned histograms for iterations." << std::endl;
    
    // Setup output
    auto outfile = TFile::Open(OUTPUT_FILE[sign], "RECREATE");
    setup_output_directories(outfile);
    
    // Setup RooFit variables
    RooRealVar kstar("kstar", "kstar", Config::KSTAR_MIN, Config::KSTAR_MAX);
    
    // Prepare models
    RooDataHist* signal_data = nullptr;
    RooHistPdf* signal_pdf = nullptr;
    prepareSignalModel(kstar, signal_data, signal_pdf, SIGNAL_FILE, "hCkHist", outfile);
    
    RooDataHist* background_data = nullptr;
    RooHistPdf* background_pdf = nullptr;
    TString backgroundHistName = Form("%s/%s/hLambdaCorrectedCk", sign, centrality);
    prepareBackgroundModel(kstar, background_data, background_pdf, BACKGROUND_FILE, 
                          backgroundHistName.Data(), outfile);
    
    // Create result histograms
    auto h_raw_yield = new TH1F("h_raw_yield", "Raw yield distribution;#it{N}_{raw}(^{4}Li);Counts", 
                                400, -200, 200);
    auto h_chi2_fit = new TH1F("h_chi2_fit", "Chi2 distribution;#chi^{2};Counts", 500, 0, 100);
    
    // Main iteration loop
    for (int iter = 0; iter < N_ITERATIONS; ++iter) {
        bool do_drawing = (iter % PRINT_INTERVAL == 0);
        
        RooRealVar nsig("nsig", "Signal Yield", 0.2, -1, 1);
        RooRealVar nbkg("nbkg", "Background Yield", 20, 0, 1e6);
        RooAddPdf model("model", "Signal + Background", 
                       RooArgList(*signal_pdf, *background_pdf), 
                       RooArgList(nsig, nbkg));
        
        if (do_drawing) {
            std::cout << "Processing iteration: " << iter << "/" << N_ITERATIONS << std::endl;
        }
        
        auto result = process_single_iteration(kstar, &model, nsig, nbkg, h_same_iter, h_mixed_iter, 
                                              h_correlation_iter, h_same, h_mixed, h_correlation, 
                                              outfile, iter, do_drawing, SAMPLING_METHOD);
        
        h_raw_yield->Fill(result.yield);
        h_chi2_fit->Fill(result.chi2);
    }
    
    // Fit yield distribution
    TF1 f_gaus("f_gaus", "gaus", -200, 200);
    h_raw_yield->Fit(&f_gaus, "RMS+");
    gStyle->SetOptFit(1);
    
    // Save results
    outfile->cd();
    h_raw_yield->Write();
    h_chi2_fit->Write();
    outfile->Close();
    
    // Cleanup
    delete h_same;
    delete h_mixed;
    delete h_correlation;
    delete h_same_iter;
    delete h_mixed_iter;
    delete h_correlation_iter;
    delete background_data;
    delete background_pdf;
    
    std::cout << "Analysis complete. Results saved to " << OUTPUT_FILE[sign] << std::endl;
}
