/**
 * This code computes the chi2 distribution for the null hypothesis 
 * and estimates the yield uncertainty for Li4 from correlation function fits
*/

#include <iostream>
#include <unordered_map>

#include <TFile.h>
#include <TH1F.h>
#include <TF1.h>
#include <TRandom3.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TString.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TStyle.h>

#include <RooRealVar.h>
#include <RooCrystalBall.h>
#include <RooDataHist.h>
#include <RooPlot.h>
#include <RooFitResult.h>
#include <RooHistPdf.h>
#include <RooAddPdf.h>
#include <RooCurve.h>
#include <RooKeysPdf.h>
#include <RooDataSet.h>

enum class SamplingMethod {
    POISSONIAN,
    GAUSSIAN
};

namespace Config {

    const char* CENTRALITY = "010";
    
    const char* DATA_FILE = "/home/galucia/Lithium4/preparation/checks/correlation_hadronpid_pass1_pass4_refined_dca.root";
    //const char* CORRELATION_FILE = "/home/galucia/Lithium4/preparation/checks/correlation_hadronpid_pass1_pass4_refined_dca.root";
    //std::unordered_map<std::string, const char*> CORRELATION_NAME = {
    //    {"Antimatter", "CorrelationAntimatter/Default/hCorrelation010"},
    //    {"Matter", "CorrelationMatter/Default/hCorrelation010"},
    //    //{"Antimatter", "CorrelationAntimatter/Default/hCorrelation050"},
    //    //{"Matter", "CorrelationMatter/Default/hCorrelation050"},
    //};
    
    const char* CORRELATION_FILE = "/home/galucia/Lithium4/preparation/output/correlation_with_systematics.root";
    std::unordered_map<std::string, const char*> CORRELATION_NAME = {
        {"Antimatter", "Antimatter/hCorrelation010Syst"},
        {"Matter", "Matter/hCorrelation010Syst"},
        //{"Antimatter", "Antimatter/hCorrelation050Syst"},
        //{"Matter", "Matter/hCorrelation050Syst"},
    };
    const char* CORRELATION_OPTION = "load"; // "load" or "compute"
    
    
    //const char* SIGNAL_FILE = "/home/galucia/Lithium4/femto/models/li4_contribution_finer_binning.root";
    //const char* SIGNAL_FILE = "/home/galucia/Lithium4/femto/models/li4_contribution.root";
    const char* SIGNAL_FILE = "/home/galucia/Lithium4/femto/models/li4_contribution_proper_sill.root";
    const char* BACKGROUND_FILE = "/home/galucia/Lithium4/femto/models/lambda_models.root";
    
    std::unordered_map<std::string, const char*> OUTPUT_FILE = {
        // ------------------------ 010 ------------------------
        // Poisson - stat
        //{"Antimatter", "/home/galucia/Lithium4/femto/output/yield_upper_limit_Antimatter_010.root"},
        //{"Matter", "/home/galucia/Lithium4/femto/output/yield_upper_limit_Matter_010.root"},
        //{"Both", "/home/galucia/Lithium4/femto/output/yield_upper_limit_Both_010.root"}

        // Gaus - stat
        //{"Antimatter", "/home/galucia/Lithium4/femto/output/yield_upper_limit_Antimatter_gaus_010.root"},
        //{"Matter", "/home/galucia/Lithium4/femto/output/yield_upper_limit_Matter_gaus_010.root"},
        
        // Gaus - syst
        {"Antimatter", "/home/galucia/Lithium4/femto/output/yield_upper_limit_Antimatter_Syst_gaus_010.root"},
        {"Matter", "/home/galucia/Lithium4/femto/output/yield_upper_limit_Matter_Syst_gaus_010.root"},
        //{"Both", "/home/galucia/Lithium4/femto/output/yield_upper_limit_Both_gaus_010.root"}

        // ------------------------ 050 ------------------------
        // Poisson - stat
        //{"Antimatter", "/home/galucia/Lithium4/femto/output/yield_upper_limit_Antimatter_050.root"},
        //{"Matter", "/home/galucia/Lithium4/femto/output/yield_upper_limit_Matter_050.root"},

        // Gaus - stat
        //{"Antimatter", "/home/galucia/Lithium4/femto/output/yield_upper_limit_Antimatter_gaus_050.root"},
        //{"Matter", "/home/galucia/Lithium4/femto/output/yield_upper_limit_Matter_gaus_050.root"},
        
        // Gaus - syst
        //{"Antimatter", "/home/galucia/Lithium4/femto/output/yield_upper_limit_Antimatter_Syst_gaus_050.root"},
        //{"Matter", "/home/galucia/Lithium4/femto/output/yield_upper_limit_Matter_Syst_gaus_050.root"}


        // ----------------------------- Hist fit -----------------------------
        // ------------------------ 010 ------------------------
        // Poisson - stat
        //{"Antimatter", "/home/galucia/Lithium4/femto/output/yield_upper_limit_Antimatter_010_hist_fit.root"},
        //{"Matter", "/home/galucia/Lithium4/femto/output/yield_upper_limit_Matter_010_hist_fit.root"},
        //{"Both", "/home/galucia/Lithium4/femto/output/yield_upper_limit_Both_010_hist_fit.root"}

        // Gaus - stat
        //{"Antimatter", "/home/galucia/Lithium4/femto/output/yield_upper_limit_Antimatter_gaus_010_hist_fit.root"},
        //{"Matter", "/home/galucia/Lithium4/femto/output/yield_upper_limit_Matter_gaus_010_hist_fit.root"},
        
        // Gaus - syst
        //{"Antimatter", "/home/galucia/Lithium4/femto/output/yield_upper_limit_Antimatter_Syst_gaus_010_hist_fit.root"},
        //{"Matter", "/home/galucia/Lithium4/femto/output/yield_upper_limit_Matter_Syst_gaus_010_hist_fit.root"},
        //{"Both", "/home/galucia/Lithium4/femto/output/yield_upper_limit_Both_gaus_010_hist_fit.root"}

        // ------------------------ 050 ------------------------
        // Poisson - stat
        //{"Antimatter", "/home/galucia/Lithium4/femto/output/yield_upper_limit_Antimatter_050_hist_fit.root"},
        //{"Matter", "/home/galucia/Lithium4/femto/output/yield_upper_limit_Matter_050_hist_fit.root"},

        // Gaus - stat
        //{"Antimatter", "/home/galucia/Lithium4/femto/output/yield_upper_limit_Antimatter_gaus_050_hist_fit.root"},
        //{"Matter", "/home/galucia/Lithium4/femto/output/yield_upper_limit_Matter_gaus_050_hist_fit.root"},
        
        // Gaus - syst
        //{"Antimatter", "/home/galucia/Lithium4/femto/output/yield_upper_limit_Antimatter_Syst_gaus_050_hist_fit.root"},
        //{"Matter", "/home/galucia/Lithium4/femto/output/yield_upper_limit_Matter_Syst_gaus_050_hist_fit.root"}
    };
    const int N_ITERATIONS = 10000;
    const int PRINT_INTERVAL = 100;
    const double HIGH_YIELD_THRESHOLD = 1000.0;
    
    const double KSTAR_MIN = 0.02;
    const double KSTAR_MAX = 0.4;
    const double PREFIT_MIN = 0.2;
    const double PREFIT_MAX = 0.4;
    const double REFERENCE_KSTAR_FOR_BKG_NORMALIZATION = 0.3; // in the flat region
    const double REFERENCE_KSTAR_FOR_SIG_NORMALIZATION = 0.07; // in the flat region

    //SamplingMethod SAMPLING_METHOD = SamplingMethod::POISSONIAN;
    SamplingMethod SAMPLING_METHOD = SamplingMethod::GAUSSIAN;
}

// Iteration result structure
struct IterationResult {
    double yield;
    double chi2;
};

// Function declarations
bool load_same_mixed(TH1F* &h_same, TH1F* &h_mixed, const char* sign = "Antimatter", const char* centrality = "010");
void compute_correlation_function(const TH1F *h_same, const TH1F *h_mixed, TH1F* &h_correlation);
void prepare_signal_model(RooRealVar &kstar, RooDataHist* &signal_data, RooHistPdf* &signal_pdf, TFile *outfile);
void prepare_background_model(RooRealVar &kstar, RooDataHist* &background_data, RooHistPdf* &background_pdf, TFile *outfile, 
                                const char* sign = "Antimatter", const char* centrality = "010");
void prefit_background(RooRealVar &kstar, TH1F *h_data, RooDataHist* &data, RooAddPdf* &model, RooRealVar& nbkg,
                        TFile *outfile, 
                       double& bkg_normalisation_at_reference, double& bkg_full_correction,
                       const char* background_name = "background_pdf", bool do_drawing = false);
double fit(RooRealVar &kstar, RooDataHist* &data, RooAddPdf* &model, TFile *outfile, 
            double& sig_normalisation_at_reference,
            const char* fit_name = "fit_result", bool do_drawing = false);
double compute_yield(RooAbsPdf* model_pdf, RooRealVar& nbkg, RooRealVar& nsig, RooRealVar& xvar, 
                    TH1F* h_mixed_event, TFile *outfile, 
                    double bkg_normalisation_at_reference, double sig_normalisation_at_reference,
                    const char * yield_name = "yield_result", bool do_drawing = false);
void save_high_yield_fit(RooRealVar& kstar, RooDataHist* &data, RooAddPdf* &model, int iter, TFile* outfile);
void setup_output_directories(TFile* outfile);
IterationResult process_single_iteration(RooRealVar& kstar, RooAddPdf* model,
                                        RooRealVar& nsig, RooRealVar& nbkg,
                                        TH1F* &h_same_iter, TH1F* &h_mixed_iter, TH1F* &h_correlation_iter,
                                        const TH1F* h_same, const TH1F* h_mixed, TFile* outfile, int iter, bool do_drawing,
                                        SamplingMethod sampling_method);
void compute_yield_error();

// Function implementations

bool load_same_mixed(TH1F* &h_same, TH1F* &h_mixed, const char* sign, const char * centrality) {
    TFile *file = TFile::Open(Config::DATA_FILE);
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Cannot open file " << Config::DATA_FILE << std::endl;
        return false;
    }
    
    const char * sign_suffix = (strcmp(sign, "Both") == 0) ? "" : sign;
    h_same = (TH1F*)file->Get(Form("Correlation%s/Default/hSameEvent%s", sign_suffix, centrality));
    h_mixed = (TH1F*)file->Get(Form("Correlation%s/Default/hNormalisedMixedEvent%s", sign_suffix, centrality));
    
    if (!h_same || !h_mixed) {
        std::cerr << "Error: Cannot load histograms from file" << std::endl;
        file->Close();
        return false;
    }
    
    h_same->SetDirectory(0);
    h_mixed->SetDirectory(0);
    file->Close();

    h_same->Print();
    h_mixed->Print();
    
    return true;
}

void prepare_correlation_function(const TH1F *h_same, const TH1F *h_mixed, TH1F* &h_correlation, const char* sign) {
    if (strcmp(Config::CORRELATION_OPTION, "load") == 0)
    {   
        std::cout << "Loading correlation function from file: " << Config::CORRELATION_FILE << std::endl;
        std::cout << "Histogram name: " << Config::CORRELATION_NAME[sign] << std::endl;
        TFile *file = TFile::Open(Config::CORRELATION_FILE);
        if (!file || file->IsZombie()) {
            std::cerr << "Error: Cannot open file " << Config::CORRELATION_FILE << std::endl;
            return;
        }
        h_correlation = (TH1F*)file->Get(Config::CORRELATION_NAME[sign]);
        h_correlation->SetDirectory(0);
        file->Close();
    } 
    else if (strcmp(Config::CORRELATION_OPTION, "compute") == 0) 
    {
        compute_correlation_function(h_same, h_mixed, h_correlation);
    } 
    else 
    {
        std::cerr << "Error: Invalid CORRELATION_OPTION in Config" << std::endl;
    }
}

void compute_correlation_function(const TH1F *h_same, const TH1F *h_mixed, TH1F* &h_correlation) {
    const int nbins = h_same->GetNbinsX();
    
    for (int ibin = 1; ibin <= nbins + 1; ++ibin) {
        const double same_val = h_same->GetBinContent(ibin);
        const double mixed_val = h_mixed->GetBinContent(ibin);
        
        if (mixed_val <= 0) {
            h_correlation->SetBinContent(ibin, 1e-12);
            h_correlation->SetBinError(ibin, 1e-12);
            continue;
        }
        
        const double ratio = same_val / mixed_val;
        const double same_err = std::sqrt(same_val);
        const double mixed_err = std::sqrt(mixed_val);
        const double rel_err = std::sqrt(
            std::pow(same_err / same_val, 2) + 
            std::pow(mixed_err / mixed_val, 2)
        );
        
        h_correlation->SetBinContent(ibin, ratio);
        h_correlation->SetBinError(ibin, ratio * rel_err);
    }
}

//void prepare_signal_model(RooRealVar &kstar, RooDataHist* &signal_data, RooHistPdf* &signal_pdf, TFile *outfile) {
//    using namespace RooFit;
//
//    TFile *file = TFile::Open(Config::SIGNAL_FILE);
//    if (!file || file->IsZombie()) {
//        std::cerr << "Error: Cannot open CATS file" << std::endl;
//        return;
//    }
//    auto h_correlation_signal = (TH1F*)file->Get("hCkHist");
//    if (!h_correlation_signal) {
//        std::cerr << "Error: Cannot load signal histogram" << std::endl;
//        file->Close();
//        return;
//    }
//    h_correlation_signal->SetDirectory(0);
//    file->Close();
//
//    signal_data = new RooDataHist("signal_data", "signal_data", kstar, Import(*h_correlation_signal));
//    signal_pdf = new RooHistPdf("signal_pdf", "signal_pdf", kstar, *signal_data);
//
//    RooPlot* frame = kstar.frame();
//    signal_data->plotOn(frame, MarkerStyle(20), MarkerSize(0.5));
//    signal_pdf->plotOn(frame, LineColor(kBlue), LineWidth(2));
//    frame->SetTitle(";#it{k}* (GeV/#it{c});C(#it{k}*)");
//    TCanvas canvas("canvas_signal_fit", "Signal Fit", 800, 600);
//    canvas.cd();
//    frame->Draw();
//
//    outfile->mkdir("signal_fit");
//    outfile->cd("signal_fit");
//    h_correlation_signal->Write("h_correlation_signal");
//    canvas.Write();
//
//    delete h_correlation_signal;
//    delete frame;
//}

void prepare_signal_model(RooRealVar &kstar, RooDataSet* &signal_data, RooKeysPdf* &signal_pdf, TFile *outfile) {
    using namespace RooFit;

    TFile *file = TFile::Open(Config::SIGNAL_FILE);
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Cannot open CATS file" << std::endl;
        return;
    }
    auto h_correlation_signal = (TH1F*)file->Get("hCkHist");
    if (!h_correlation_signal) {
        std::cerr << "Error: Cannot load signal histogram" << std::endl;
        file->Close();
        return;
    }
    h_correlation_signal->SetDirectory(0);
    file->Close();

    kstar.setRange(0.01, 0.42);
    RooRealVar weight_var("weight", "weight", 0, 1e6);
    RooArgSet vars(kstar, weight_var);
    signal_data = new RooDataSet("signal_data", "signal_data", vars, RooFit::WeightVar(weight_var));

    for (int ibin = 1; ibin <= h_correlation_signal->GetNbinsX(); ++ibin) {
        double x_val = h_correlation_signal->GetBinCenter(ibin);
        double y_val = h_correlation_signal->GetBinContent(ibin);
        if (x_val < kstar.getMin() || x_val > kstar.getMax()) continue;
        if (y_val <= 0) continue;

        kstar.setVal(x_val);
        weight_var.setVal(y_val);
        signal_data->add(vars, y_val);
    }

    double rho = 2.;
    signal_pdf = new RooKeysPdf("signal_pdf", "signal_pdf",
                                    kstar, *signal_data,
                                    RooKeysPdf::NoMirror, rho);

    RooPlot* frame = kstar.frame();
    signal_data->plotOn(frame, MarkerStyle(20), MarkerSize(0.5));
    signal_pdf->plotOn(frame, LineColor(kBlue), LineWidth(2));
    frame->SetTitle(";#it{k}* (GeV/#it{c});C(#it{k}*)");
    TCanvas canvas("canvas_signal_fit", "Signal Fit", 800, 600);
    canvas.cd();
    frame->Draw();

    outfile->mkdir("signal_fit");
    outfile->cd("signal_fit");
    h_correlation_signal->Write("h_correlation_signal");
    canvas.Write();

    kstar.setRange(Config::KSTAR_MIN, Config::KSTAR_MAX); // Ensure range is set correctly for later use

    delete h_correlation_signal;
    delete frame;
}

//void prepare_background_model(RooRealVar &kstar, RooDataHist* &background_data, RooHistPdf* &background_pdf, TFile *outfile, 
//                                const char* sign, const char* centrality) {
//    using namespace RooFit;
//
//    TFile *file = TFile::Open(Config::BACKGROUND_FILE);
//    if (!file || file->IsZombie()) {
//        std::cerr << "Error: Cannot open BACKGROUND file" << std::endl;
//        return;
//    }
//    auto h_correlation_background = (TH1F*)file->Get(Form("%s/%s/hLambdaCorrectedCk", sign, centrality));
//    if (!h_correlation_background) {
//        std::cerr << "Error: Cannot load background histogram" << std::endl;
//        file->Close();
//        return;
//    }
//    h_correlation_background->SetDirectory(0);
//    file->Close();
//
//    background_data = new RooDataHist("background_data", "background_data", kstar, Import(*h_correlation_background));
//    background_pdf = new RooHistPdf("background_pdf", "background_pdf", kstar, *background_data);
//
//    RooPlot* frame = kstar.frame();
//    background_data->plotOn(frame, MarkerStyle(20), MarkerSize(0.5));
//    background_pdf->plotOn(frame, LineColor(kBlue), LineWidth(2));
//    frame->SetTitle(";#it{k}* (GeV/#it{c});C(#it{k}*)");
//    TCanvas canvas("canvas_background_fit", "Background Fit", 800, 600);
//    canvas.cd();
//    frame->Draw();
//
//    outfile->mkdir("background_fit");
//    outfile->cd("background_fit");
//    canvas.Write("background_fit");
//
//    delete h_correlation_background;
//    delete frame;
//}

void prepare_background_model(RooRealVar &kstar, RooDataSet* &background_data, RooKeysPdf* &background_pdf, TFile *outfile, 
                                const char* sign, const char* centrality) {
    using namespace RooFit;

    TFile *file = TFile::Open(Config::BACKGROUND_FILE);
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Cannot open BACKGROUND file" << std::endl;
        return;
    }
    //auto h_correlation_background = (TH1F*)file->Get(Form("%s/%s/hLambdaSigmaCorrectedCk", sign, centrality));
    auto h_correlation_background = (TH1F*)file->Get(Form("%s/%s/hLambdaSigmaCorrectedCk_Smeared", sign, centrality));
    if (!h_correlation_background) {
        std::cerr << "Error: Cannot load background histogram" << std::endl;
        file->Close();
        return;
    }
    h_correlation_background->SetDirectory(0);
    file->Close();

    kstar.setRange(0.01, 0.42);
    RooRealVar weight_var("weight", "weight", 0, 1e6);
    RooArgSet vars(kstar, weight_var);
    background_data = new RooDataSet("background_data", "background_data", vars, RooFit::WeightVar(weight_var));

    for (int ibin = 1; ibin <= h_correlation_background->GetNbinsX(); ++ibin) {
        double x_val = h_correlation_background->GetBinCenter(ibin);
        double y_val = h_correlation_background->GetBinContent(ibin);
        if (x_val < kstar.getMin() || x_val > kstar.getMax()) continue;
        if (y_val <= 0) continue;

        kstar.setVal(x_val);
        weight_var.setVal(y_val);
        background_data->add(vars, y_val);
    }

    double rho = 0.1; // KDE bandwidth — tune this (smaller = less smooth, larger = more smooth)
    background_pdf = new RooKeysPdf("background_pdf", "background_pdf",
                                    kstar, *background_data,
                                    RooKeysPdf::NoMirror, rho);

    RooPlot* frame = kstar.frame();
    background_data->plotOn(frame, MarkerStyle(20), MarkerSize(0.5));
    background_pdf->plotOn(frame, LineColor(kBlue), LineWidth(2));
    frame->SetTitle(";#it{k}* (GeV/#it{c});C(#it{k}*)");
    TCanvas canvas("canvas_background_fit", "Background Fit", 800, 600);
    canvas.cd();
    frame->Draw();

    outfile->mkdir("background_fit");
    outfile->cd("background_fit");
    canvas.Write("background_fit");

    kstar.setRange(Config::KSTAR_MIN, Config::KSTAR_MAX); // Ensure range is set correctly for later use

    delete h_correlation_background;
    delete frame;
}

void prefit_background(RooRealVar &kstar, TH1F* h_data, RooDataHist* &data, RooAddPdf* &model, RooRealVar& nbkg, TFile *outfile, 
                       double& bkg_normalisation_at_reference, double& bkg_full_correction,
                       const char* background_name, bool do_drawing) {
    using namespace RooFit;

    // Store the old range
    const double old_min = kstar.getMin();
    const double old_max = kstar.getMax();
    
    // Set ranges for integral calculation
    kstar.setRange("full", old_min, old_max);
    kstar.setRange("prefit", Config::PREFIT_MIN, Config::PREFIT_MAX);
    kstar.setRange(Config::PREFIT_MIN, Config::PREFIT_MAX);

    data = new RooDataHist("data_prefit", "data_prefit", kstar, Import(*h_data));
    
    model->chi2FitTo(*data, PrintLevel(-2), SumW2Error(true), 
                     Range("prefit"), 
                     Extended(true), Save());

    // Get the background PDF component - cast RooAbsArg to RooAbsPdf
    RooAbsArg* bkg_arg = model->pdfList().at(1);
    RooAbsPdf* bkg_pdf = dynamic_cast<RooAbsPdf*>(bkg_arg);
    
    if (!bkg_pdf) {
        std::cerr << "Error: Could not cast background component to RooAbsPdf" << std::endl;
        return;
    }
    
    RooPlot* frame = kstar.frame();
    data->plotOn(frame, MarkerStyle(20), MarkerSize(0.5));
    model->plotOn(frame, LineColor(kBlue), LineWidth(2), Name("background_model"),
                    Normalization(1.0, RooAbsReal::RelativeExpected));
    frame->SetTitle(";#it{k}* (GeV/#it{c});C(#it{k}*)");
    
    TCanvas canvas("canvas_prefit_background", "Prefit Background", 800, 600);
    frame->Draw();
    

    auto integral_full = bkg_pdf->createIntegral(kstar, NormSet(kstar), Range("full"));
    auto integral_prefit = bkg_pdf->createIntegral(kstar, NormSet(kstar), Range("prefit"));
    
    const double I_full = integral_full->getVal();
    const double I_prefit = integral_prefit->getVal();
    bkg_full_correction = I_full / I_prefit;

    const double norm = nbkg.getVal();
    nbkg.setVal(norm * bkg_full_correction);
    nbkg.setConstant(true);
    
    delete integral_full;
    delete integral_prefit;

    const auto bkg_curve = (RooCurve*)frame->findObject("background_model");
    // Get the normalization AFTER applying the correction
    bkg_normalisation_at_reference = bkg_curve->interpolate(Config::REFERENCE_KSTAR_FOR_BKG_NORMALIZATION);

    if (do_drawing) {

        TLatex latex;
        latex.SetNDC(); 
        latex.SetTextSize(0.04);
        
        canvas.cd();
        latex.DrawLatex(0.15, 0.85, Form("N_{bkg} (reference) = %.3f", bkg_normalisation_at_reference));
        
        outfile->cd("prefit_background");
        canvas.Write(background_name);
    }

    kstar.setRange(old_min, old_max); // Restore original range

    delete frame;
}

double fit(RooRealVar &kstar, RooDataHist* &data, RooAddPdf* &model,
            TFile *outfile, 
            double& sig_normalisation_at_reference,
            const char* fit_name, bool do_drawing) {
    using namespace RooFit;

    RooFitResult* result = model->chi2FitTo(*data, PrintLevel(-2), SumW2Error(true),
                                            //NormRange("prefit"), 
                                           Range("full"), 
                                           Extended(true), Save());

    RooPlot* frame = kstar.frame(kstar.getMin(), kstar.getMax());
    data->plotOn(frame, MarkerStyle(20), MarkerSize(0.5));
    model->plotOn(frame, LineColor(kBlue), LineWidth(2), Name("total_model"),
                    Normalization(1.0, RooAbsReal::RelativeExpected));
    model->plotOn(frame, Components("background_pdf"), LineStyle(kDashed), LineColor(kRed), LineWidth(2),
                  Normalization(1.0, RooAbsReal::RelativeExpected));
    model->plotOn(frame, Components("signal_pdf"), LineStyle(kDashed), LineColor(kGreen+2), LineWidth(2),
                  Name("signal_model"),
                  Normalization(1.0, RooAbsReal::RelativeExpected));
    model->paramOn(frame, Layout(0.45, 0.85, 0.35), Format("NEU", AutoPrecision(2)));
    frame->SetTitle(";#it{k}* (GeV/#it{c});C(#it{k}*)");

    RooAbsArg* sig_arg = model->pdfList().at(0);
    RooAbsPdf* sig_pdf = dynamic_cast<RooAbsPdf*>(sig_arg);
    const auto sig_curve = (RooCurve*)frame->findObject("signal_model");
    sig_normalisation_at_reference = sig_curve->interpolate(Config::REFERENCE_KSTAR_FOR_SIG_NORMALIZATION);

    double chi2 = frame->chiSquare();

    if (do_drawing) {
        TCanvas canvas("canvas_fit", "Fit Result", 800, 600);
        canvas.cd();
        frame->Draw();

        outfile->cd("fit_results");
        canvas.Write(fit_name);
    }
    
    delete frame;
    delete result;

    return chi2;
}

double compute_yield(RooAbsPdf* model_pdf, RooRealVar& nsig, RooRealVar& nbkg, RooRealVar& xvar, 
                    TH1F* h_mixed_event, TFile *outfile, 
                    double bkg_normalisation_at_reference, double sig_normalisation_at_reference,
                    const char * yield_name, bool do_drawing) {
    
    // Store original values
    const double nsig_stored = nsig.getVal();
    const double nbkg_stored = nbkg.getVal();

    TH1F* h_signal_correlation = (TH1F*)h_mixed_event->Clone("h_signal_correlation");
    TH1F* h_same_event_signal = (TH1F*)h_mixed_event->Clone("h_same_event_signal");
    TH1F* h_background_correlation = (TH1F*)h_mixed_event->Clone("h_background_correlation");
    const int nbins = h_same_event_signal->GetNbinsX();
    
    // Get the normalization correction for the background AT THE REFERENCE POINT
    xvar.setVal(Config::REFERENCE_KSTAR_FOR_BKG_NORMALIZATION);
    nsig.setVal(0.);
    nbkg.setVal(nbkg_stored);
    const double bkg_value_at_reference = model_pdf->getVal(xvar);
    const double bkg_correction = bkg_normalisation_at_reference / bkg_value_at_reference;

    // Get the noramalization correction for the signal AT THE REFERENCE POINT
    xvar.setVal(Config::REFERENCE_KSTAR_FOR_SIG_NORMALIZATION);
    nsig.setVal(nsig_stored);
    nbkg.setVal(0.);
    const double sig_value_at_reference = model_pdf->getVal(xvar);
    const double sig_correction = sig_normalisation_at_reference / sig_value_at_reference;
    
    for (int ibin = 1; ibin <= nbins; ++ibin) {
        const double x = h_same_event_signal->GetBinCenter(ibin);
        xvar.setVal(x);
        
        // Evaluate signal only (background = 0)
        nsig.setVal(nsig_stored);
        nbkg.setVal(0.);
        const double signal_value = model_pdf->getVal(xvar) * sig_correction;
        
        // Evaluate background only (signal = 0) 
        nsig.setVal(0.);
        nbkg.setVal(nbkg_stored);
        const double bkg_value = model_pdf->getVal(xvar) * bkg_correction;
        
        h_signal_correlation->SetBinContent(ibin, signal_value);
        h_background_correlation->SetBinContent(ibin, bkg_value);
        h_same_event_signal->SetBinContent(ibin, signal_value * h_mixed_event->GetBinContent(ibin));
    }

    const int last_bin = h_same_event_signal->FindBin(Config::KSTAR_MAX);
    double yield = h_same_event_signal->Integral(1, last_bin);

    if (do_drawing) {
        TCanvas canvas_corr("canvas_signal_correlation", "Signal Correlation", 800, 600);
        canvas_corr.cd();
        
        // Calculate the total model for comparison
        TH1F* h_total_correlation = (TH1F*)h_mixed_event->Clone("h_total_correlation");
        for (int ibin = 1; ibin <= nbins; ++ibin) {
            h_total_correlation->SetBinContent(ibin, 
                h_signal_correlation->GetBinContent(ibin) + 
                h_background_correlation->GetBinContent(ibin));
        }
        
        h_total_correlation->SetLineColor(kBlack);
        h_total_correlation->SetMarkerColor(kBlack);
        h_total_correlation->SetMarkerStyle(20);
        h_total_correlation->SetTitle(";#it{k}* (GeV/#it{c});C(#it{k}*)");
        
        h_signal_correlation->SetLineColor(kGreen+2);
        h_signal_correlation->SetMarkerColor(kGreen+2);
        h_signal_correlation->SetMarkerStyle(20);
        h_signal_correlation->SetLineStyle(kDashed);
        
        h_background_correlation->SetLineColor(kRed);
        h_background_correlation->SetMarkerColor(kRed);
        h_background_correlation->SetMarkerStyle(20);
        h_background_correlation->SetLineStyle(kDashed);
        
        h_total_correlation->Draw("hist");
        h_signal_correlation->Draw("hist same");
        h_background_correlation->Draw("hist same");
        
        TLegend legend(0.6, 0.6, 0.85, 0.85);
        legend.AddEntry(h_total_correlation, "Total", "l");
        legend.AddEntry(h_signal_correlation, "Signal", "l");
        legend.AddEntry(h_background_correlation, "Background", "l");
        legend.SetBorderSize(0);
        legend.Draw();
        
        outfile->cd("yield_results");
        canvas_corr.Write(Form("%s_components", yield_name));
        
        delete h_total_correlation;
    }

    if (do_drawing) {
        outfile->cd("yield_results");
        h_signal_correlation->Write(Form("%s_correlation", yield_name));
        h_same_event_signal->Write(yield_name);
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

void save_high_yield_fit(RooRealVar& kstar, RooDataHist* &data, RooAddPdf* &model, 
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

void setup_output_directories(TFile* outfile) {
    outfile->mkdir("prefit_background");
    outfile->mkdir("fit_results");
    outfile->mkdir("yield_results");
    outfile->mkdir("high_yield_results");
}

IterationResult process_single_iteration(RooRealVar& kstar, RooAddPdf* model, RooRealVar& nsig, RooRealVar& nbkg,
                                            TH1F* &h_same_iter, TH1F* &h_mixed_iter, TH1F* &h_correlation_iter,
                                            const TH1F* h_same, const TH1F* h_mixed, const TH1F* h_correlation,
                                            TFile* outfile, int iter, bool do_drawing,
                                            SamplingMethod sampling_method) {

    const int nbins = h_same->GetNbinsX();

    if (sampling_method == SamplingMethod::POISSONIAN) {

        for (int ibin = 1; ibin <= nbins; ++ibin) {
            h_same_iter->SetBinContent(ibin, gRandom->Poisson(h_same->GetBinContent(ibin)));
            h_mixed_iter->SetBinContent(ibin, gRandom->Poisson(h_mixed->GetBinContent(ibin)));
        }
        compute_correlation_function(h_same_iter, h_mixed_iter, h_correlation_iter);
    } 
    else if (sampling_method == SamplingMethod::GAUSSIAN) {

        for (int ibin = 1; ibin <= nbins; ++ibin) {
            const double central_value = h_correlation->GetBinContent(ibin);
            const double error = h_correlation->GetBinError(ibin);
            
            const double sampled_value = gRandom->Gaus(central_value, error);
            h_correlation_iter->SetBinContent(ibin, sampled_value);
            h_correlation_iter->SetBinError(ibin, error);
        }
        
        for (int ibin = 1; ibin <= nbins; ++ibin) {
            h_mixed_iter->SetBinContent(ibin, h_mixed->GetBinContent(ibin));
        }
    }

    RooDataHist* data_prefit = nullptr;
    
    nsig.setVal(0);
    nsig.setConstant(true);
    double bkg_normalisation_at_reference = 1.0;
    double bkg_full_correction = 1.0;
    prefit_background(kstar, h_correlation_iter, data_prefit, model, nbkg,
                        outfile, bkg_normalisation_at_reference, bkg_full_correction,
                        Form("background_pdf_%d", iter), do_drawing);

    delete data_prefit;
    
    nsig.setConstant(false);
    nbkg.setConstant(true);
    double sig_normalisation_at_reference = 1.0;

    RooDataHist* data_iter = new RooDataHist("data_iter", "data_iter", kstar, RooFit::Import(*h_correlation_iter));
    
    double chi2 = fit(kstar, data_iter, model, outfile, sig_normalisation_at_reference,
                        Form("fit_result_%d", iter), do_drawing);
    double yield = compute_yield(model, nsig, nbkg, kstar, h_mixed_iter, outfile, 
                            bkg_normalisation_at_reference, sig_normalisation_at_reference,
                            Form("yield_result_%d", iter), do_drawing);

    if (yield > Config::HIGH_YIELD_THRESHOLD) {
        save_high_yield_fit(kstar, data_iter, model, iter, outfile);
    }
    
    delete data_iter;

    return {yield, chi2};
}

void compute_yield_and_upper_limit() {

    const char* centrality = Config::CENTRALITY;

    for (const auto& sign : {"Antimatter", "Matter"}) {
        
        TH1F *h_same, *h_mixed;
        if (!load_same_mixed(h_same, h_mixed, sign, centrality)) {
            std::cerr << "Failed to load histograms. Exiting." << std::endl;
            return;
        }

        h_same->Print();
        h_mixed->Print();
        
        auto h_correlation = (TH1F*)h_same->Clone("h_correlation");
        std::cout << "Cloned histogram for correlation function." << std::endl;
        prepare_correlation_function(h_same, h_mixed, h_correlation, sign);
        std::cout << "Loaded histograms: " << h_same->GetName() << " and " << h_mixed->GetName() << std::endl;
        
        auto h_same_iter = (TH1F*)h_same->Clone("h_same_iter");
        auto h_mixed_iter = (TH1F*)h_mixed->Clone("h_mixed_iter");
        auto h_correlation_iter = (TH1F*)h_correlation->Clone("h_correlation_iter");
        std::cout << "Cloned histograms for iterations." << std::endl;
        
        auto outfile = TFile::Open(Config::OUTPUT_FILE[sign], "RECREATE");
        h_correlation_iter->Write("h_correlation_iter");
        setup_output_directories(outfile);

        
        RooRealVar kstar("kstar", "kstar", Config::KSTAR_MIN, Config::KSTAR_MAX);

        //RooDataHist* signal_data = nullptr;
        //RooHistPdf* signal_pdf = nullptr;
        RooDataSet* signal_data = nullptr;
        RooKeysPdf* signal_pdf = nullptr;
        prepare_signal_model(kstar, signal_data, signal_pdf, outfile);
        
        //RooDataHist* background_data = nullptr;
        //RooHistPdf* background_pdf = nullptr;
        RooDataSet* background_data = nullptr;
        RooKeysPdf* background_pdf = nullptr;
        prepare_background_model(kstar, background_data, background_pdf, outfile,
                                sign, centrality);
        
        auto h_raw_yield = new TH1F("h_raw_yield", "Raw yield distribution;#it{N}_{raw}(^{4}Li);Counts", 
                                    400, -800, 800);
        auto h_chi2_fit = new TH1F("h_chi2_fit", "Chi2 distribution;#chi^{2};Counts", 500, 0, 100);
        
        for (int iter = 0; iter < Config::N_ITERATIONS; ++iter) {
            bool do_drawing = (iter % Config::PRINT_INTERVAL == 0);
            
            RooRealVar nsig("nsig", "Signal Yield", 0.2, -1, 1);
            RooRealVar nbkg("nbkg", "Background Yield", 20, 0, 1e6);
            RooAddPdf model("model", "Signal + Background", 
                            RooArgList(*signal_pdf, *background_pdf), 
                            RooArgList(nsig, nbkg));
            
            if (do_drawing) {
                std::cout << "Processing iteration: " << iter << "/" << Config::N_ITERATIONS << std::endl;
            }
            
            auto result = process_single_iteration(
                kstar, &model, nsig, nbkg,
                h_same_iter, h_mixed_iter, h_correlation_iter,
                h_same, h_mixed, h_correlation,
                outfile, iter, do_drawing,
                Config::SAMPLING_METHOD
            );
            
            h_raw_yield->Fill(result.yield);
            h_chi2_fit->Fill(result.chi2);
        }
        
        TF1 f_gaus("f_gaus", "gaus", -800, 800);
        h_raw_yield->Fit(&f_gaus, "RMS+");
        gStyle->SetOptFit(1);

        outfile->cd();
        h_raw_yield->Write();
        h_chi2_fit->Write();

        std::cout << "Finished processing for sign: " << sign << std::endl;
        
        //delete h_same;
        //delete h_mixed;
        //delete h_correlation;
        //delete h_same_iter;
        //delete h_mixed_iter;
        //delete h_correlation_iter;
        //delete background_data;
        //delete background_pdf;
        
        outfile->Close();
        
        std::cout << "Analysis complete. Results saved to " << Config::OUTPUT_FILE[sign] << std::endl;

    }
}