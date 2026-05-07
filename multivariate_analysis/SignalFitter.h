#pragma once

#include "Fitter.h"

#include "TDirectory.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TTree.h"
#include "RooRealVar.h"
#include "RooCrystalBall.h"
#include "RooFit.h"
#include "RooHistPdf.h"
#include "RooDataHist.h"
#include "RooWorkspace.h"
#include "RooDataSet.h"
#include "RooKeysPdf.h"
#include "RooPlot.h"

#include <string>
#include <map>
#include <vector>
#include <stdexcept>
#include <cmath>

/**
 * SignalFitter.h
 * Equivalent to femto/core/signal_fitter.py
 */
class SignalFitter : public Fitter {

public:
    SignalFitter(const std::string& name,
                 const std::string& xvar_name,
                 const std::string& xvar_title,
                 double xmin, double xmax,
                 TDirectory* outfile        = nullptr,
                 RooWorkspace* workspace = nullptr)
        : Fitter(name, xvar_name, xvar_title, xmin, xmax, outfile, workspace)
        , _signal_pdf(nullptr)
        , _signal_datahist(nullptr)
        , _signal_dataset(nullptr)
        , _title("")
    {
        _outdir = (_outfile != nullptr) ? _outfile->mkdir("signal") : nullptr;
    }

    // -----------------------------------------------------------------------
    // Properties
    // -----------------------------------------------------------------------
    const std::map<std::string, RooRealVar*>& signal_pars() const { return _signal_pars; }
    RooAbsPdf* signal_pdf() const { return _signal_pdf; }

    const std::string& title() const { return _title; }
    void set_title(const std::string& t) {
        _title = t;
        if (_signal_pdf) _signal_pdf->SetTitle(_title.c_str());
    }

    // -----------------------------------------------------------------------
    // Initialisation
    // -----------------------------------------------------------------------
    void init_signal(const std::string& mode, TH1D* h_signal = nullptr,
                     const std::string& name = "signal_pdf",
                     double xmin = 0., double xmax = 0.42, double rho = 2.)
    {
        if      (mode == "from_mc")       _init_signal_from_mc(h_signal, name);
        else if (mode == "from_kde")      _init_signal_from_kde(h_signal, name, xmin, xmax, rho);
        else if (mode == "crystal_ball")  _init_signal_crystal_ball(name);
        else throw std::invalid_argument("Supported modes: from_mc, from_kde, crystal_ball");
    }

    // -----------------------------------------------------------------------
    // Fit
    // -----------------------------------------------------------------------
    void fit_signal(double range_low, double range_high, TH1D* h_signal) {
        auto* xvar = dynamic_cast<RooRealVar*>(_roo_workspace->obj(_xvar_name.c_str()));
        RooDataHist datahist("signal_dh", "signal_dh", *xvar, RooFit::Import(*h_signal));

        _signal_pdf->fitTo(datahist,
            RooFit::Save(),
            RooFit::Range(range_low, range_high),
            RooFit::SumW2Error(true),
            RooFit::Extended(true),
            RooFit::PrintLevel(-1),
            RooFit::Verbose(false));

        RooPlot* frame = xvar->frame();
        datahist.plotOn(frame);
        _signal_pdf->plotOn(frame);

        if (_outdir) {
            _outdir->cd();
            TCanvas canvas("fit_signal");
            frame->Draw();
            canvas.Write();
        }
        delete frame;
    }

    // -----------------------------------------------------------------------
    void save_to_workspace() override {
        _roo_workspace->import(*_signal_pdf);
        for (auto& [key, par] : _signal_pars)
            _roo_workspace->import(*par);
    }

private:
    // -----------------------------------------------------------------------
    void _init_signal_from_mc(TH1D* h_signal, const std::string& name) {
        // Zero the unphysical first bin (acts as underflow)
        h_signal->SetBinContent(h_signal->FindBin(0.), 0.);

        auto* xvar = dynamic_cast<RooRealVar*>(_roo_workspace->obj(_xvar_name.c_str()));
        _signal_datahist = new RooDataHist("signal_dh", "signal_dh", *xvar, RooFit::Import(*h_signal));
        _signal_pdf = new RooHistPdf(name.c_str(), name.c_str(), *xvar, *_signal_datahist);

        RooPlot* frame = xvar->frame();
        _signal_datahist->plotOn(frame);
        _signal_pdf->plotOn(frame);

        if (_outdir) {
            _outdir->cd();
            TCanvas canvas("signal_pdf");
            frame->Draw();
            canvas.Write();
        }
        delete frame;
    }

    // -----------------------------------------------------------------------
    void _init_signal_from_kde(TH1D* h_signal, const std::string& name,
                                double xmin, double xmax, double rho) {
        auto* xvar = dynamic_cast<RooRealVar*>(_roo_workspace->obj(_xvar_name.c_str()));
        const double old_min = xvar->getMin();
        const double old_max = xvar->getMax();
        xvar->setRange(xmin, xmax);

        // Collect (x, weight) pairs inside [xmin, xmax]
        std::vector<double> x_data, weights;
        for (int ibin = 1; ibin <= h_signal->GetNbinsX(); ++ibin) {
            double x_val = h_signal->GetBinCenter(ibin);
            if (x_val < xmin || x_val > xmax) continue;
            double y_val = h_signal->GetBinContent(ibin);
            if (y_val <= 0.) continue;
            x_data.push_back(x_val);
            weights.push_back(y_val);
        }

        // Fill a TTree so RooDataSet can be constructed with weights
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
        std::string ds_name = std::string(h_signal->GetName()) + "_roodata";
        _signal_dataset = new RooDataSet(ds_name.c_str(), ds_name.c_str(),
                                         &tree, RooArgSet(*xvar, weight_var), "", "weight");
        _roo_workspace->import(*_signal_dataset);

        _signal_pdf = new RooKeysPdf(name.c_str(), name.c_str(),
                                     *xvar, *_signal_dataset,
                                     RooKeysPdf::NoMirror, rho);

        RooPlot* frame = xvar->frame(static_cast<int>(x_data.size()));
        _signal_dataset->plotOn(frame, RooFit::MarkerStyle(20), RooFit::MarkerSize(0.8), RooFit::LineColor(kBlack));
        _signal_pdf->plotOn(frame,     RooFit::LineColor(kRed),  RooFit::LineWidth(2));

        TCanvas canvas(("cKeysPdf_" + std::string(h_signal->GetName())).c_str(),
                        ("cKeysPdf_" + std::string(h_signal->GetName())).c_str(), 800, 600);
        frame->Draw();

        if (_outdir) {
            _outdir->cd();
            h_signal->Write((std::string(h_signal->GetName()) + "_original").c_str());
            _signal_pdf->Write((std::string(h_signal->GetName()) + "_keyspdf").c_str());
            canvas.Write();
        }
        delete frame;

        xvar->setRange(old_min, old_max);
    }

    // -----------------------------------------------------------------------
    void _init_signal_crystal_ball(const std::string& name) {
        auto* xvar = dynamic_cast<RooRealVar*>(_roo_workspace->obj(_xvar_name.c_str()));
        _signal_pars["mean"]  = new RooRealVar("sig_mean",  "#mu",    0.081, 0.06,  0.1);
        _signal_pars["sigma"] = new RooRealVar("sig_sigma", "#sigma", 0.01,  0.001, 0.1);
        _signal_pars["aL"]    = new RooRealVar("sig_aL",    "a_{L}",  1.3,   0.1,   10.);
        _signal_pars["nL"]    = new RooRealVar("sig_nL",    "n_{L}",  2.7,   0.1,   10.);
        _signal_pars["aR"]    = new RooRealVar("sig_aR",    "a_{R}",  1.1,   0.1,   10.);
        _signal_pars["nR"]    = new RooRealVar("sig_nR",    "n_{R}",  5.7,   0.1,   10.);

        _signal_pdf = new RooCrystalBall(name.c_str(), name.c_str(),
                                          *xvar,
                                          *_signal_pars["mean"],
                                          *_signal_pars["sigma"],
                                          *_signal_pars["aL"],
                                          *_signal_pars["nL"],
                                          *_signal_pars["aR"],
                                          *_signal_pars["nR"]);
    }

    // -----------------------------------------------------------------------
    // Members
    // -----------------------------------------------------------------------
    std::map<std::string, RooRealVar*> _signal_pars;
    RooAbsPdf*   _signal_pdf;
    RooDataHist* _signal_datahist;
    RooDataSet*  _signal_dataset;
    std::string  _title;
    TDirectory*  _outdir;
};
