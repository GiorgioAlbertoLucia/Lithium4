#pragma once

#include "Fitter.h"

#include "TDirectory.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "RooRealVar.h"
#include "RooFit.h"
#include "RooHistPdf.h"
#include "RooDataHist.h"
#include "RooWorkspace.h"
#include "RooAddPdf.h"
#include "RooArgList.h"
#include "RooPlot.h"
#include "RooAbsReal.h"

#include <string>
#include <vector>
#include <map>
#include <stdexcept>
#include <cmath>

/**
 * ModelFitter.h
 * Equivalent to femto/core/model_fitter.py
 */
class ModelFitter : public Fitter {

public:
    ModelFitter(const std::string& name,
                const std::string& xvar_name,
                const std::string& xvar_title,
                double xmin, double xmax,
                TDirectory* outfile,
                const std::vector<std::string>& signal_func_names,
                const std::vector<std::string>& bkg_func_names,
                RooWorkspace* workspace = nullptr,
                bool extended = false)
        : Fitter(name, xvar_name, xvar_title, xmin, xmax, outfile, workspace)
        , _model_pdf(nullptr)
        , _roo_data_hist(nullptr)
        , reference_kstar_value_for_normalisation(0.31)
        , bkg_normalisations_at_reference_kstar(0.)
        , reference_kstar_value_for_signal_normalisation(0.07)
        , signal_normalisations_at_reference_kstar(0.)
    {
        if (signal_func_names.empty() || bkg_func_names.empty())
            throw std::invalid_argument("At least one signal and one bkg name must be provided");

        _outdir = (outfile != nullptr) ? outfile->mkdir("model") : nullptr;
        _init_model(name, signal_func_names, bkg_func_names, extended);
    }

    // -----------------------------------------------------------------------
    // Properties
    // -----------------------------------------------------------------------
    RooAddPdf* model_pdf() const { return _model_pdf; }

    // Publicly accessible fractions (mirroring Python's self.fractions dict)
    std::map<std::string, RooRealVar*> fractions;

    // Reference values (public, like Python attributes)
    double reference_kstar_value_for_normalisation;
    double bkg_normalisations_at_reference_kstar;
    double reference_kstar_value_for_signal_normalisation;
    double signal_normalisations_at_reference_kstar;

    // -----------------------------------------------------------------------
    // Pre-fit: fix signal to 0 and fit the background region
    // -----------------------------------------------------------------------
    void prefit_background(TH1D* h_data,
                            double range_low = -1., double range_high = -1.,
                            const std::string& range_name = "",
                            bool use_chi2_method = true,
                            bool save_normalisation_value = true)
    {
        auto* xvar = dynamic_cast<RooRealVar*>(_roo_workspace->obj(_xvar_name.c_str()));
        const double old_min = xvar->getMin();
        const double old_max = xvar->getMax();
        if (range_low  < 0.) range_low  = old_min;
        if (range_high < 0.) range_high = old_max;

        xvar->setRange(range_low, range_high);
        if (!range_name.empty())
            xvar->setRange(range_name.c_str(), range_low, range_high);

        // Fix signal fractions to 0
        for (auto& [sig_name, _] : _signal_pdfs) {
            if (fractions.count(sig_name)) {
                fractions[sig_name]->setVal(0.);
                fractions[sig_name]->setConstant(true);
            }
        }

        RooDataHist datahist((std::string(h_data->GetName()) + "_datahist").c_str(),
                             (std::string(h_data->GetName()) + "_datahist").c_str(),
                             *xvar, RooFit::Import(*h_data));

        if (use_chi2_method)
            _model_pdf->chi2FitTo(datahist,
                RooFit::Save(), RooFit::Extended(true),
                !range_name.empty() ? RooFit::Range(range_name.c_str()) : RooFit::Save(), // dummy reuse
                RooFit::PrintLevel(-1), RooFit::Verbose(false));
        else
            _model_pdf->fitTo(datahist,
                RooFit::Save(), RooFit::Extended(true),
                !range_name.empty() ? RooFit::Range(range_name.c_str()) : RooFit::Save(),
                RooFit::PrintLevel(-1), RooFit::Verbose(false));

        RooPlot* frame = xvar->frame();
        _plot_model(frame, &datahist, "prefit_bkg");

        xvar->setRange("full",   old_min,    old_max);
        xvar->setRange("prefit", range_low,  range_high);

        // Correct bkg fractions from the sideband to the full range
        for (auto& [bkg_name, bkg_pdf] : _bkg_pdfs) {
            if (!fractions.count(bkg_name)) continue;
            RooAbsReal* I_full = bkg_pdf->createIntegral(*xvar,
                RooFit::NormSet(*xvar), RooFit::Range("full"));
            RooAbsReal* I_side = bkg_pdf->createIntegral(*xvar,
                RooFit::NormSet(*xvar), RooFit::Range("prefit"));
            double bkg_full_correction = I_full->getVal() / I_side->getVal();
            double norm = fractions[bkg_name]->getVal();
            fractions[bkg_name]->setVal(norm * bkg_full_correction);
            fractions[bkg_name]->setConstant(true);
            delete I_full; delete I_side;
        }

        // Release signal fractions
        for (auto& [sig_name, _] : _signal_pdfs)
            if (fractions.count(sig_name))
                fractions[sig_name]->setConstant(false);

        // Store reference background level for chi2 / raw yield calculation
        RooCurve* bkg_curve = dynamic_cast<RooCurve*>(frame->findObject(_model_pdf->GetName()));
        if (bkg_curve)
            bkg_normalisations_at_reference_kstar =
                bkg_curve->interpolate(reference_kstar_value_for_normalisation);

        xvar->setRange(old_min, old_max);
        delete frame;
    }

    // -----------------------------------------------------------------------
    // Full model fit
    // -----------------------------------------------------------------------
    void fit_model(TH1D* h_data,
                   const std::string& signal_name,
                   bool use_chi2_fit_method = true,
                   const std::string& norm_range = "")
    {
        auto* xvar = dynamic_cast<RooRealVar*>(_roo_workspace->obj(_xvar_name.c_str()));
        _roo_data_hist = new RooDataHist((std::string(h_data->GetName()) + "_datahist").c_str(),
                                          (std::string(h_data->GetName()) + "_datahist").c_str(),
                                          *xvar, RooFit::Import(*h_data));

        if (use_chi2_fit_method)
            _model_pdf->chi2FitTo(*_roo_data_hist,
                RooFit::Save(), RooFit::SumW2Error(true), RooFit::Extended(true),
                !norm_range.empty() ? RooFit::NormRange(norm_range.c_str()) : RooFit::Save(),
                RooFit::PrintLevel(-1), RooFit::Verbose(false));
        else
            _model_pdf->fitTo(*_roo_data_hist,
                RooFit::Save(), RooFit::SumW2Error(true), RooFit::Extended(true),
                !norm_range.empty() ? RooFit::NormRange(norm_range.c_str()) : RooFit::Save(),
                RooFit::PrintLevel(-1), RooFit::Verbose(false));

        RooPlot* frame = xvar->frame(xvar->getMin(), xvar->getMax());
        _plot_model(frame, _roo_data_hist, "fit_signal");

        // Store reference signal level
        if (_signal_pdfs.count(signal_name)) {
            RooCurve* sig_curve = dynamic_cast<RooCurve*>(
                frame->findObject(_signal_pdfs[signal_name]->GetName()));
            if (sig_curve)
                signal_normalisations_at_reference_kstar =
                    sig_curve->interpolate(reference_kstar_value_for_signal_normalisation);
        }
        delete frame;
    }

    // -----------------------------------------------------------------------
    void save_to_workspace() override {
        _roo_workspace->import(*_model_pdf);
    }

    // -----------------------------------------------------------------------
    // Chi2 between data and background-only model
    // -----------------------------------------------------------------------
    void compute_chi2(TH1D* h_data) {
        auto* xvar = dynamic_cast<RooRealVar*>(_roo_workspace->obj(_xvar_name.c_str()));

        int nbins = h_data->GetNbinsX();
        TH1D h_chi2("chi2",     ";#it{k}* (GeV/#it{c}); #chi^{2}",       nbins, h_data->GetBinLowEdge(1), h_data->GetBinLowEdge(nbins+1));
        TH1D h_nsigma("nsigma", ";#it{k}* (GeV/#it{c}); n#sigma",        nbins, h_data->GetBinLowEdge(1), h_data->GetBinLowEdge(nbins+1));
        TH1D h_bkg_check("bkg_check", ";#it{k}* (GeV/#it{c}); C(k*)",   nbins, h_data->GetBinLowEdge(1), h_data->GetBinLowEdge(nbins+1));
        TH1D h_chi2_ndf("chi2_ndf", ";#it{k}* (GeV/#it{c}); #chi^{2}/NDF", nbins, h_data->GetBinLowEdge(1), h_data->GetBinLowEdge(nbins+1));

        // Save and zero signal fractions
        std::map<std::string, double> stored_signal_fractions;
        for (auto& [sig_name, _] : _signal_pdfs) {
            stored_signal_fractions[sig_name] = fractions[sig_name]->getVal();
            fractions[sig_name]->setVal(0.);
        }

        double chi2 = 0., ndf = 0.;
        for (int ibin = 1; ibin <= nbins; ++ibin) {
            double kstar_value = h_data->GetBinCenter(ibin);
            double data_value  = h_data->GetBinContent(ibin);
            double data_error  = h_data->GetBinError(ibin);

            xvar->setVal(reference_kstar_value_for_normalisation);
            double bkg_at_ref = _model_pdf->getVal(RooArgSet(*xvar));
            double correction = bkg_normalisations_at_reference_kstar / bkg_at_ref;

            xvar->setVal(kstar_value);
            double bkg_value = _model_pdf->getVal(RooArgSet(*xvar)) * correction;
            double bkg_error = 0.;

            double difference  = data_value - bkg_value;
            double uncertainty = std::sqrt(data_error * data_error + bkg_error * bkg_error);
            double nsigma      = (uncertainty > 0.) ? difference / uncertainty : 0.;
            chi2 += nsigma * nsigma;
            ndf  += 1.;

            h_chi2.SetBinContent    (ibin, chi2);
            h_chi2_ndf.SetBinContent(ibin, chi2 / ndf);
            h_nsigma.SetBinContent  (ibin, nsigma);
            h_bkg_check.SetBinContent(ibin, bkg_value);
            h_bkg_check.SetBinError  (ibin, bkg_error);
        }

        // Restore signal fractions
        for (auto& [sig_name, val] : stored_signal_fractions)
            fractions[sig_name]->setVal(val);

        TCanvas canvas("data_bkg_comparison", "");
        h_data->Draw("e1");
        h_bkg_check.Draw("e1 same");
        TLegend* legend = canvas.BuildLegend(0.5, 0.3, 0.8, 0.5);
        legend->SetBorderSize(0);

        if (_outdir) {
            _outdir->cd();
            h_chi2.Write();
            h_nsigma.Write();
            h_chi2_ndf.Write();
            h_bkg_check.Write();
            canvas.Write();
        }
    }

    // -----------------------------------------------------------------------
    // Raw yield extraction
    // -----------------------------------------------------------------------
    double compute_raw_yield(TH1D* h_mixed_event,
                              const std::string& signal_pdf_name,
                              const std::string& bkg_pdf_name)
    {
        auto* xvar = dynamic_cast<RooRealVar*>(_roo_workspace->obj(_xvar_name.c_str()));

        double nsig_stored  = fractions[signal_pdf_name]->getVal();
        double nbkg_stored  = fractions[bkg_pdf_name]->getVal();
        int    nbins        = h_mixed_event->GetNbinsX();

        TH1D* h_signal_correlation  = static_cast<TH1D*>(h_mixed_event->Clone("h_signal_correlation"));
        TH1D* h_background_correlation = static_cast<TH1D*>(h_mixed_event->Clone("h_background_correlation"));
        TH1D* h_same_event_signal   = static_cast<TH1D*>(h_mixed_event->Clone("h_same_event_signal"));

        // Compute bkg correction factor at reference kstar
        xvar->setVal(reference_kstar_value_for_normalisation);
        fractions[signal_pdf_name]->setVal(0.);
        fractions[bkg_pdf_name]->setVal(nbkg_stored);
        double bkg_at_ref  = _model_pdf->getVal(RooArgSet(*xvar));
        double bkg_correction = bkg_normalisations_at_reference_kstar / bkg_at_ref;

        // Compute signal correction factor at reference kstar
        xvar->setVal(reference_kstar_value_for_signal_normalisation);
        fractions[signal_pdf_name]->setVal(nsig_stored);
        fractions[bkg_pdf_name]->setVal(0.);
        double sig_at_ref  = _model_pdf->getVal(RooArgSet(*xvar));
        double signal_correction = signal_normalisations_at_reference_kstar / sig_at_ref;

        for (int ibin = 1; ibin <= nbins; ++ibin) {
            double kstar_value = h_mixed_event->GetBinCenter(ibin);
            xvar->setVal(kstar_value);

            fractions[signal_pdf_name]->setVal(0.);
            fractions[bkg_pdf_name]->setVal(nbkg_stored);
            double bkg_value = _model_pdf->getVal(RooArgSet(*xvar)) * bkg_correction;

            fractions[signal_pdf_name]->setVal(nsig_stored);
            fractions[bkg_pdf_name]->setVal(0.);
            double signal_value = _model_pdf->getVal(RooArgSet(*xvar)) * signal_correction;

            h_signal_correlation->SetBinContent(ibin, signal_value);
            h_background_correlation->SetBinContent(ibin, bkg_value);
            h_same_event_signal->SetBinContent(ibin, signal_value * h_mixed_event->GetBinContent(ibin));
        }

        int last_bin  = h_same_event_signal->FindBin(xvar->getMax());
        double yield  = h_same_event_signal->Integral(1, last_bin);

        TPaveText text(0.5, 0.5, 0.8, 0.7, "ndc");
        text.SetFillStyle(0);
        text.SetBorderSize(0);
        text.SetTextSize(0.04);
        text.AddText(Form("Raw yield = %.2f", yield));

        TCanvas canvas("yield_extraction", "");
        h_same_event_signal->Draw("hist");
        text.Draw("same");

        if (_outdir) {
            _outdir->cd();
            h_same_event_signal->Write();
            canvas.Write();
        }

        // Restore fractions
        fractions[signal_pdf_name]->setVal(nsig_stored);
        fractions[bkg_pdf_name]->setVal(nbkg_stored);

        delete h_signal_correlation;
        delete h_background_correlation;
        delete h_same_event_signal;

        return yield;
    }

private:
    // -----------------------------------------------------------------------
    void _init_model(const std::string& name,
                     const std::vector<std::string>& signal_func_names,
                     const std::vector<std::string>& bkg_func_names,
                     bool extended)
    {
        RooArgList pdf_list, fraction_list;

        for (const auto& sig_name : signal_func_names) {
            auto* pdf = dynamic_cast<RooAbsPdf*>(_roo_workspace->obj(sig_name.c_str()));
            _signal_pdfs[sig_name] = pdf;

            std::string title = pdf->GetTitle();
            auto pos = title.find(';');
            if (pos != std::string::npos) title = title.substr(0, pos);

            auto* frac = new RooRealVar((sig_name + "_frac").c_str(),
                                        ("#it{f}_{" + title + "}").c_str(),
                                        0.5, 0., 1.);
            if (extended) {
                frac->setRange(0., 1e4);
                frac->setVal(0.2);
                frac->SetTitle(("#it{N}_{" + title + "}").c_str());
            }
            fractions[sig_name] = frac;
            pdf_list.add(*pdf);
            fraction_list.add(*frac);
        }

        for (size_t ibkg = 0; ibkg < bkg_func_names.size(); ++ibkg) {
            const auto& bkg_name = bkg_func_names[ibkg];
            auto* pdf = dynamic_cast<RooAbsPdf*>(_roo_workspace->obj(bkg_name.c_str()));
            _bkg_pdfs[bkg_name] = pdf;
            pdf_list.add(*pdf);

            std::string title = pdf->GetTitle();
            auto pos = title.find(';');
            if (pos != std::string::npos) title = title.substr(0, pos);

            // Last bkg has no explicit fraction in non-extended mode
            if (ibkg == bkg_func_names.size() - 1 && !extended) continue;

            auto* frac = new RooRealVar((bkg_name + "_frac").c_str(),
                                        ("#it{f}_{" + title + "}").c_str(),
                                        0.5, 0., 1.);
            if (extended) {
                frac->setRange(0., 1e4);
                frac->SetTitle(("#it{N}_{" + title + "}").c_str());
            }
            fractions[bkg_name] = frac;
            fraction_list.add(*frac);
        }

        _model_pdf = new RooAddPdf(name.c_str(), name.c_str(), pdf_list, fraction_list);
    }

    // -----------------------------------------------------------------------
    void _plot_model(RooPlot* frame, RooDataHist* roodatahist, const std::string& canvas_name) {
        int line_color = 2;
        roodatahist->plotOn(frame);
        _model_pdf->plotOn(frame,
            RooFit::Name(_model_pdf->GetName()),
            RooFit::Title(_model_pdf->GetTitle()),
            RooFit::Normalization(1.0, RooAbsReal::RelativeExpected),
            RooFit::LineColor(line_color++));

        for (auto& [sig_name, sig_pdf] : _signal_pdfs) {
            RooArgSet comp_set(*sig_pdf);
            _model_pdf->plotOn(frame,
                RooFit::Name(sig_pdf->GetName()),
                RooFit::Title(sig_pdf->GetTitle()),
                RooFit::Normalization(1.0, RooAbsReal::RelativeExpected),
                RooFit::Components(comp_set),
                RooFit::LineColor(line_color++));
        }
        for (auto& [bkg_name, bkg_pdf] : _bkg_pdfs) {
            RooArgSet comp_set(*bkg_pdf);
            _model_pdf->plotOn(frame,
                RooFit::Name(bkg_pdf->GetName()),
                RooFit::Title(bkg_pdf->GetTitle()),
                RooFit::Normalization(1.0, RooAbsReal::RelativeExpected),
                RooFit::Components(comp_set),
                RooFit::LineColor(line_color++));
        }

        // Text box with fit parameters
        TPaveText text(0.58, 0.43, 0.81, 0.58, "NDC");
        text.SetFillColor(0);
        text.SetBorderSize(0);
        text.SetTextSize(0.044);
        if (fractions.count("signal_pdf") && fractions.count("bkg_pdf")) {
            text.AddText(Form("#bf{%s = %.2f}",
                fractions["signal_pdf"]->GetTitle(), fractions["signal_pdf"]->getVal()));
            text.AddText(Form("#bf{%s = %.2f (fixed)}",
                fractions["bkg_pdf"]->GetTitle(), fractions["bkg_pdf"]->getVal()));
        }
        text.AddText(Form("#bf{#chi^{2} / ndf = %.2f}", frame->chiSquare()));

        TLegend legend(0.58, 0.23, 0.81, 0.38);
        legend.SetBorderSize(0);
        legend.SetTextSize(0.045);
        if (frame->findObject(_model_pdf->GetName()))
            legend.AddEntry(frame->findObject(_model_pdf->GetName()), "Total", "l");
        for (auto& [_, sig_pdf] : _signal_pdfs)
            if (frame->findObject(sig_pdf->GetName()))
                legend.AddEntry(frame->findObject(sig_pdf->GetName()), sig_pdf->GetTitle(), "l");
        for (auto& [_, bkg_pdf] : _bkg_pdfs)
            if (frame->findObject(bkg_pdf->GetName()))
                legend.AddEntry(frame->findObject(bkg_pdf->GetName()), bkg_pdf->GetTitle(), "l");

        TCanvas canvas(canvas_name.c_str());
        frame->Draw();
        text.Draw("same");
        legend.Draw();

        if (_outdir) {
            _outdir->cd();
            canvas.Write();
        }
    }

    // -----------------------------------------------------------------------
    // Members
    // -----------------------------------------------------------------------
    std::map<std::string, RooAbsPdf*> _signal_pdfs;
    std::map<std::string, RooAbsPdf*> _bkg_pdfs;
    RooAddPdf*   _model_pdf;
    RooDataHist* _roo_data_hist;
    TDirectory*  _outdir;
};
