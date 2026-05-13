
#include <iostream>
#include <TGenPhaseSpace.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TStyle.h>
#include <TPad.h>
#include <TFile.h>


void he4l_decay() {
    
    TFile outfile("output/he4l_decay.root", "RECREATE");
    
    // Particle masses in GeV/c^2
    const Double_t M_hyperHe4 = 3.921;  // Hyperhelium-4 mass
    const Double_t M_proton = 0.938272; // Proton mass
    const Double_t M_He3 = 2.809;       // Helium-3 mass
    const Double_t M_pion = 0.13957;    // Pion mass
    
    // Number of events to generate
    const Int_t nEvents = 1000000;
    
    // Create TGenPhaseSpace object
    TGenPhaseSpace event;
    
    // Set up the decay: hyperhelium-4 at rest
    TLorentzVector W(0.0, 0.0, 0.0, M_hyperHe4);
    
    // Daughter masses
    Double_t masses[3] = {M_proton, M_He3, M_pion};
    
    // Initialize the decay
    event.SetDecay(W, 3, masses);
    
    // Create histograms for invariant masses
    TH1F *h_invmass_all = new TH1F("h_invmass_all", 
        "Invariant Mass of All Daughters;#it{M}(p + He3 + #pi^{-}) (GeV/#it{c}^{2});Counts", 
        200, 3.743, 3.943);
    
    TH1F *h_invmass_pHe3 = new TH1F("h_invmass_pHe3", 
        "Invariant Mass of Proton-Helium3 Pair;#it{M}(p + He3) (GeV/#it{c}^{2});Counts", 
        200, 3.743, 3.943);

    TH1F *h_kstar_pHe3 = new TH1F("h_kstar_pHe3", 
        "Relative Momentum k* of Proton-Helium3 Pair;#it{k}* (GeV/#it{c});Counts", 
        50, 0.0, 0.5);

    TH1F *h_correlation_pHe3 = new TH1F("h_correlation_pHe3", 
        "Contribution to the correlation function of p-^{3}He from the ^{4}_{#Lambda}He weak decay;#it{k}* (GeV/#it{c});C_{^{4}_{#Lambda}He}(#it{k}*)", 
        50, 0.0, 0.5);

    TFile infileMixedEvent("/home/galucia/Lithium4/preparation/checks/correlation_hadronpid_pass1_pass4_refined_dca.root");
    TH1F* h_mixed_event_kstar_pHe3 = (TH1F*)infileMixedEvent.Get("Correlation/Default/hNormalisedMixedEvent050");
    h_mixed_event_kstar_pHe3->Print();
    
    // Generate events
    for (Int_t i = 0; i < nEvents; i++) {

        Double_t weight = event.Generate();
        
        TLorentzVector *pProton = event.GetDecay(0);
        TLorentzVector *pHe3 = event.GetDecay(1);
        TLorentzVector *pPion = event.GetDecay(2);
        
        TLorentzVector sum_all = *pProton + *pHe3 + *pPion;
        h_invmass_all->Fill(sum_all.M());
        
        TLorentzVector sum_pHe3 = *pProton + *pHe3;
        h_invmass_pHe3->Fill(sum_pHe3.M());

        TLorentzVector pProtonStar(*pProton);
        TLorentzVector pHe3Star(*pHe3);

        auto boostVector = -sum_pHe3.BoostVector();
        pProtonStar.Boost(boostVector);
        pHe3Star.Boost(boostVector);

        const double kstar = 0.5 * (pProtonStar - pHe3Star).P();
        h_kstar_pHe3->Fill(kstar);
    }

    for (int ibin = 1; ibin < h_kstar_pHe3->GetNbinsX() + 1; ibin++) {
        double kstar = h_kstar_pHe3->GetBinCenter(ibin);
        double mixed_count = h_mixed_event_kstar_pHe3->GetBinContent(ibin);
        double signal_count = h_kstar_pHe3->GetBinContent(ibin);
        
        if (mixed_count > 0) {
            double ratio = signal_count / mixed_count;
            h_correlation_pHe3->SetBinContent(ibin, ratio);
        }
    }
    
    // Create canvas and draw histograms
    TCanvas *c1 = new TCanvas("c1", "Hyperhelium-4 Decay", 1200, 600);
    c1->Divide(2, 1);
    
    // Plot invariant mass of all daughters
    c1->cd(1);
    gPad->SetLeftMargin(0.12);
    h_invmass_all->SetLineColor(kBlue+1);
    h_invmass_all->SetLineWidth(2);
    h_invmass_all->Draw();
    gPad->SetGrid();
    
    // Add statistics
    gStyle->SetOptStat(1111);
    
    // Plot invariant mass of proton-helium3 pair
    c1->cd(2);
    gPad->SetLeftMargin(0.12);
    h_invmass_pHe3->SetLineColor(kRed+1);
    h_invmass_pHe3->SetLineWidth(2);
    h_invmass_pHe3->Draw();
    gPad->SetGrid();

    outfile.cd();
    h_invmass_all->Write();
    h_invmass_pHe3->Write();
    h_kstar_pHe3->Write();
    h_mixed_event_kstar_pHe3->Write();
    h_correlation_pHe3->Write();
    c1->Write();
    
    // Print summary information
    std::cout << "\n=== Hyperhelium-4 Decay Simulation ===" << std::endl;
    std::cout << "Mother particle mass: " << M_hyperHe4 << " GeV/c^2" << std::endl;
    std::cout << "Decay products: p + He3 + pi-" << std::endl;
    std::cout << "Number of events generated: " << nEvents << std::endl;
    std::cout << "\nInvariant mass of all daughters:" << std::endl;
    std::cout << "  Mean: " << h_invmass_all->GetMean() << " GeV/c^2" << std::endl;
    std::cout << "  RMS:  " << h_invmass_all->GetRMS() << " GeV/c^2" << std::endl;
    std::cout << "\nInvariant mass of p-He3 pair:" << std::endl;
    std::cout << "  Mean: " << h_invmass_pHe3->GetMean() << " GeV/c^2" << std::endl;
    std::cout << "  RMS:  " << h_invmass_pHe3->GetRMS() << " GeV/c^2" << std::endl;
}
