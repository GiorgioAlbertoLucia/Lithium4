#include <string>

#include <TFile.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH1D.h>
#include <ROOT/RDataFrame.hxx>

constexpr int kHePdg = 1000020030;
constexpr int isPhysicalPrimary = 0x4;

void checkMcParticle() {


    const char * inputFileName = "/data/galucia/lithium_local/raw/LHC25g12_001_raw.root";
    auto infile = TFile::Open(inputFileName);
    
    TChain chain("mcParticleChain");
    for (const auto& key: *infile->GetListOfKeys()) {
        if(std::string(key->GetName()).rfind("DF_", 0) == 0) 
            chain.Add((std::string(inputFileName)+"/"+std::string(key->GetName())+"/O2mcparticle_001").c_str());
    }

    ROOT::RDataFrame rdf(chain);
    std::cout << "[\n";
    for (const auto& column: rdf.GetColumnNames())
        std::cout << column << ", ";
    std::cout << "\n]\n";

    auto rdfNew = rdf.Define("fY", "0.5 * std::log((fE + fPz)/(fE - fPz))").
                      Define("fP", "std::hypot(fPx, fPy, fPz)").
                      Define("fEta", "0.5 * std::log((fP + fPz)/(fP - fPz))").
                      Filter("std::abs(fPdgCode) == kHePdg && std::abs(fY) < 1").
                      Define("fPt", "std::hypot(fPx, fPy)");

    auto histPrimaries = rdfNew.Filter("(fFlags & isPhysicalPrimary) == isPhysicalPrimary").Histo1D("fPt");
    auto histSecondaries = rdfNew.Filter("(fFlags & isPhysicalPrimary) != isPhysicalPrimary").Histo1D("fPt");
    auto histEtaSecondaries = rdfNew.Filter("(fFlags & isPhysicalPrimary) != isPhysicalPrimary").Histo1D("fEta");

    TCanvas canvas;
    histPrimaries->Draw("hist");
    canvas.SaveAs("checkMcParticlePrimaries.pdf");
    canvas.Clear();
    histSecondaries->Draw("hist");
    canvas.SaveAs("checkMcParticleSecondaries.pdf");
    canvas.Clear();
    histEtaSecondaries->Draw("hist");
    canvas.SaveAs("checkMcParticleEtaSecondaries.pdf");

}