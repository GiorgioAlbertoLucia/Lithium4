#pragma once

#include "TH1D.h"
#include "TDirectory.h"
#include "RooDataHist.h"
#include "RooWorkspace.h"
#include "RooRealVar.h"

#include <string>
#include <stdexcept>

/**
 * Fitter.h
 * Abstract base class — equivalent to femto/core/fitter.py
 */
class Fitter {

public:
    Fitter(const std::string& name,
           const std::string& xvar_name,
           const std::string& xvar_title,
           double xmin, double xmax,
           TDirectory* outfile = nullptr,
           RooWorkspace* workspace = nullptr)
        : _hist_data(nullptr)
        , _outfile(outfile)
        , _xvar_name(xvar_name)
        , _roo_data_hist_name("")
    {
        if (workspace == nullptr) {
            _roo_workspace = new RooWorkspace(name.empty() ? "roows" : name.c_str());
            RooRealVar xvar(xvar_name.c_str(), xvar_title.c_str(), xmin, xmax);
            _roo_workspace->import(xvar);
        } else {
            _roo_workspace = workspace;
            if (!_roo_workspace->obj(xvar_name.c_str())) {
                RooRealVar xvar(xvar_name.c_str(), xvar_title.c_str(), xmin, xmax);
                _roo_workspace->import(xvar);
            }
        }
    }

    virtual ~Fitter() = default;

    // Accessors
    RooWorkspace* roo_workspace() const { return _roo_workspace; }
    const std::string& xvar_name()  const { return _xvar_name; }

    /**
     * Import histogram as RooDataHist into the workspace.
     */
    void load_data(TH1D* hist_data, const std::string& name = "datahist") {
        auto* xvar = dynamic_cast<RooRealVar*>(_roo_workspace->obj(_xvar_name.c_str()));
        _hist_data = hist_data;
        RooDataHist rdh(name.c_str(), name.c_str(), *xvar, RooFit::Import(*hist_data));
        _roo_workspace->import(rdh);
        _roo_data_hist_name = name;
    }

    virtual void save_to_workspace() = 0;

protected:
    TH1D*         _hist_data;
    TDirectory*   _outfile;
    std::string   _xvar_name;
    RooWorkspace* _roo_workspace; // non-owning if provided from outside
    std::string   _roo_data_hist_name;
};
