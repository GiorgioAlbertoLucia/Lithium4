#pragma once

#include "Fitter.h"

#include "TDirectory.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TTree.h"
#include "RooRealVar.h"
#include "RooFit.h"
#include "RooHistPdf.h"
#include "RooDataHist.h"
#include "RooWorkspace.h"
#include "RooDataSet.h"
#include "RooKeysPdf.h"
#include "RooExtendPdf.h"
#include "RooAbsReal.h"
#include "RooPlot.h"

#include <string>
#include <map>
#include <vector>
#include <stdexcept>
#include <cmath>

/**
 * BkgFitter.h
 * Equivalent to femto/core/bkg_fitter.py
 */
class BkgFitter : public Fitter {

public:
    BkgFitter(const std::string& name,
              const std::string& xvar_name,
              const std::string& xvar_title,
              double xmin, double xmax,
              TDirectory* outfile = nullptr,
              RooWorkspace* workspace = nullptr)
        : Fitter(name, xvar_name, xvar_title, xmin, xmax, outfile, workspace)
        , _bkg_pdf(nullptr)
        , _bkg_datahist(nullptr)
        , _bkg_dataset(nullptr)
        , _bkg_normalisation(nullptr)
        , _title("")
    {
        _outdir = (_outfile != nullptr) ? _outfile->mkdir("bkg") : nullptr;
    }

    // -----------------------------------------------------------------------
    // Properties
    // -----------------------------------------------------------------------
    const std::map<std::string, RooRealVar*>& bkg_pars() const { return _bkg_pars; }
    RooRealVar* bkg_normalisation() const { return _bkg_normalisation; }
    RooAbsReal*  bkg_pdf()           const { return _bkg_pdf; }

    const std::string& title() const { return _title; }
    void set_title(const std::string& t) {
        _title = t;
        if (_bkg_pdf) _bkg_pdf->SetTitle(_title.c_str());
    }

    // -----------------------------------------------------------------------
    // Initialisation dispatcher
    // -----------------------------------------------------------------------
    void init_bkg(const std::string& mode, TH1D* h_bkg,
                  const std::string& name = "bkg_pdf",
                  bool extended = false,
                  double xmin = 0.01, double xmax = 0.42, double rho = 0.05)
    {
        if      (mode == "from_mc")  _init_bkg_from_mc (h_bkg, name, extended);
        else if (mode == "from_kde") _init_bkg_from_kde(h_bkg, name, xmin, xmax, rho);
        else throw std::invalid_argument("Supported modes: from_mc, from_kde");
    }

    // -----------------------------------------------------------------------
    // Fit the background in a given range
    // -----------------------------------------------------------------------
    void fit_bkg(TH1D* hist,
                 double range_low = -1., double range_high = -1.,
                 const std::string& range_name = "",
                 bool save_normalisation_value = true)
    {
        auto* xvar = dynamic_cast<RooRealVar*>(_roo_workspace->obj(_xvar_name.c_str()));
        const double old_min = xvar->getMin();
        const double old_max = xvar->getMax();

        if (range_low < 0.)  range_low  = old_min;
        if (range_high < 0.) range_high = old_max;

        xvar->setRange(range_low, range_high);
        if (!range_name.empty())
            xvar->setRange(range_name.c_str(), range_low, range_high);

        RooDataHist datahist((std::string(hist->GetName()) + "_datahist").c_str(),
                             (std::string(hist->GetName()) + "_datahist").c_str(),
                             *xvar, RooFit::Import(*hist));

        if (!range_name.empty())
            _bkg_pdf->chi2FitTo(datahist, RooFit::Save(), RooFit::Range(range_name.c_str()), RooFit::PrintLevel(-1), RooFit::Verbose(false));
        else
            _bkg_pdf->chi2FitTo(datahist, RooFit::Save(), RooFit::PrintLevel(-1), RooFit::Verbose(false));

        // Normalisation over the full range
        xvar->setRange("full", range_low, range_high);
        RooAbsReal* integral_full = _bkg_pdf->createIntegral(*xvar, RooFit::Range("full"));
        if (save_normalisation_value && _bkg_normalisation)
            _bkg_normalisation->setVal(integral_full->getVal());
        delete integral_full;

        RooPlot* frame = xvar->frame();
        datahist.plotOn(frame);
        _bkg_pdf->plotOn(frame);

        if (_outdir) {
            _outdir->cd();
            TCanvas canvas("fit_bkg");
            frame->Draw();
            canvas.Write();
        }
        delete frame;

        xvar->setRange(old_min, old_max);
    }

    // -----------------------------------------------------------------------
    void save_to_workspace() override {
        _roo_workspace->import(*_bkg_pdf);
        for (auto& [key, par] : _bkg_pars)
            _roo_workspace->import(*par);
        if (_bkg_datahist)
            _roo_workspace->import(*_bkg_datahist);
        if (_bkg_normalisation)
            _roo_workspace->import(*_bkg_normalisation);
    }

private:
    // -----------------------------------------------------------------------
    void _init_bkg_from_mc(TH1D* h_bkg, const std::string& name, bool extended) {
        auto* xvar = dynamic_cast<RooRealVar*>(_roo_workspace->obj(_xvar_name.c_str()));

        // Zero bins outside the xvar range
        for (int ibin = 1; ibin <= h_bkg->GetNbinsX(); ++ibin) {
            double center = h_bkg->GetBinCenter(ibin);
            if (center < xvar->getMin() || center > xvar->getMax())
                h_bkg->SetBinContent(ibin, 0.);
        }

        _bkg_datahist = new RooDataHist("bkg_dh", "bkg_dh", *xvar, RooFit::Import(*h_bkg));

        if (extended) {
            auto* tmp = new RooHistPdf("tmp", "tmp", *xvar, *_bkg_datahist);
            _bkg_normalisation = new RooRealVar("bkg_normalisation", "#it{N}_{bkg}", 1., 0., 1e4);
            _bkg_pdf = new RooExtendPdf(name.c_str(), name.c_str(), *tmp, *_bkg_normalisation);
        } else {
            _bkg_pdf = new RooHistPdf(name.c_str(), name.c_str(), *xvar, *_bkg_datahist);
        }

        RooPlot* frame = xvar->frame();
        _bkg_datahist->plotOn(frame);
        _bkg_pdf->plotOn(frame);

        if (_outdir) {
            _outdir->cd();
            TCanvas canvas("bkg_pdf");
            frame->Draw();
            canvas.Write();
        }
        delete frame;
    }

    // -----------------------------------------------------------------------
    void _init_bkg_from_kde(TH1D* h_bkg, const std::string& name,
                             double xmin, double xmax, double rho) {
        auto* xvar = dynamic_cast<RooRealVar*>(_roo_workspace->obj(_xvar_name.c_str()));
        const double old_min = xvar->getMin();
        const double old_max = xvar->getMax();
        xvar->setRange(xmin, xmax);

        std::vector<double> x_data, weights;
        for (int ibin = 1; ibin <= h_bkg->GetNbinsX(); ++ibin) {
            double x_val = h_bkg->GetBinCenter(ibin);
            if (x_val < xmin || x_val > xmax) continue;
            double y_val = h_bkg->GetBinContent(ibin);
            if (y_val <= 0.) continue;
            x_data.push_back(x_val);
            weights.push_back(y_val);
        }

        TTree tree("tree", "tree");
        double x_branch{0.}, w_branch{0.};
        tree.Branch("kstar",  &x_branch, "kstar/D");
        tree.Branch("weight", &w_branch, "weight/D");
        for (size_t i = 0; i < x_data.size(); ++i) {
            x_branch = x_data[i];
            w_branch = weights[i];
            tree.Fill();
        }

        RooRealVar weight_var("weight", "weight", 0., 1e6);
        std::string ds_name = std::string(h_bkg->GetName()) + "_roodata";
        _bkg_dataset = new RooDataSet(ds_name.c_str(), ds_name.c_str(),
                                       &tree, RooArgSet(*xvar, weight_var), "", "weight");
        _roo_workspace->import(*_bkg_dataset);

        _bkg_pdf = new RooKeysPdf(name.c_str(), name.c_str(),
                                   *xvar, *_bkg_dataset,
                                   RooKeysPdf::NoMirror, rho);

        RooPlot* frame = xvar->frame(static_cast<int>(x_data.size()));
        _bkg_dataset->plotOn(frame, RooFit::MarkerStyle(20), RooFit::MarkerSize(0.8), RooFit::LineColor(kBlack));
        _bkg_pdf->plotOn(frame,     RooFit::LineColor(kRed),  RooFit::LineWidth(2));

        TCanvas canvas(("cKeysPdf_" + std::string(h_bkg->GetName())).c_str(),
                        ("cKeysPdf_" + std::string(h_bkg->GetName())).c_str(), 800, 600);
        frame->Draw();

        if (_outdir) {
            _outdir->cd();
            h_bkg->Write((std::string(h_bkg->GetName()) + "_original").c_str());
            _bkg_pdf->Write((std::string(h_bkg->GetName()) + "_keyspdf").c_str());
            canvas.Write();
        }
        delete frame;

        xvar->setRange(old_min, old_max);
    }

    // -----------------------------------------------------------------------
    // Members
    // -----------------------------------------------------------------------
    std::map<std::string, RooRealVar*> _bkg_pars;
    RooAbsReal*   _bkg_pdf;
    RooDataHist* _bkg_datahist;
    RooDataSet*  _bkg_dataset;
    RooRealVar*  _bkg_normalisation;
    std::string  _title;
    TDirectory*  _outdir;
};
