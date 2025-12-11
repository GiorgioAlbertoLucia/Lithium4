#include <string>
#include <vector>
#include <TH1D.h>
#include <TString.h>
#include <ROOT/RResultPtr.hxx>

void Normalise(TH1D* h_same, TH1D* h_mixed, double low, double high)
{
    int low_bin  = h_same->FindBin(low);
    int high_bin = h_same->FindBin(high);

    double same_int  = h_same->Integral(low_bin, high_bin);
    double mixed_int = h_mixed->Integral(low_bin, high_bin);

    if (mixed_int <= 0) return;

    double norm = same_int / mixed_int;
    h_mixed->Scale(norm);
}

TH1D* ComputeCorrelation(const std::vector<TH1D*>& same, const std::vector<TH1D*>& mixed, int suffix)
{
    TH1D* h_same  = (TH1D*) same[0]->Clone();
    TH1D* h_mixed = (TH1D*) mixed[0]->Clone();

    for (size_t i = 1; i < same.size(); ++i) {
        h_same->Add(same[i]);
        h_mixed->Add(mixed[i]);
    }

    TString name = TString::Format("hCorrelation%d", suffix);
    TH1D* h_corr = (TH1D*) h_same->Clone(name);
    h_corr->Divide(h_mixed);

    return h_corr;
}

using RMap = ROOT::RDF::Experimental::RResultMap<TH1D>;
std::vector<TH1D*> ComputeAllSystematics(std::vector<RMap>& same,
                                        std::vector<RMap>& mixed,
                                        int N_ITER, double low, double high)
{
    std::vector<TH1D*> h_correlations;
    h_correlations.reserve(N_ITER);

    for (int iter = 0; iter < N_ITER; ++iter) {

        TString key = TString::Format("systematic:%d", iter);
        auto h010_same = same[0][std::string(key.Data())];
        auto h1030_same = same[1][std::string(key.Data())];
        auto h3050_same = same[2][std::string(key.Data())];

        auto h010_mixed = (TH1D*) mixed[0][std::string(key.Data())].Clone();
        auto h1030_mixed = (TH1D*) mixed[1][std::string(key.Data())].Clone();
        auto h3050_mixed = (TH1D*) mixed[2][std::string(key.Data())].Clone();

        // normalisation
        Normalise(&h010_same, h010_mixed, low, high);
        Normalise(&h1030_same, h1030_mixed, low, high);
        Normalise(&h3050_same, h3050_mixed, low, high);

        // correlation
        TH1D* h_correlation = ComputeCorrelation({&h010_same, &h1030_same, &h3050_same},
                                                 {h010_mixed, h1030_mixed, h3050_mixed}, iter);

        h_correlations.push_back(h_correlation);
    }

    return h_correlations;
}
