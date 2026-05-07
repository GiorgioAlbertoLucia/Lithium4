#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <iostream>
#include <chrono>

void filterTree(const char* inputFile,
                const char* outputFile,
                const char* treeName = "your_tree_name")  // <-- change this
{
    // --- Open input ---
    TFile* fIn = TFile::Open(inputFile, "READ");
    if (!fIn || fIn->IsZombie()) {
        std::cerr << "Cannot open input file\n";
        return;
    }

    TTree* tIn = dynamic_cast<TTree*>(fIn->Get(treeName));
    if (!tIn) {
        std::cerr << "Cannot find tree: " << treeName << "\n";
        return;
    }

    Long64_t nEntries = tIn->GetEntries();
    std::cout << "Total entries: " << nEntries << "\n";

    // Optimize read speed: cache 100 MB in memory
    tIn->SetCacheSize(100 * 1024 * 1024);
    tIn->AddBranchToCache("*", true);

    // --- Open output ---
    // Compression: 101 = ZLIB level 1 (fast), 404 = LZ4 level 4 (faster on modern CPUs)
    TFile* fOut = TFile::Open(outputFile, "RECREATE", "", ROOT::RCompressionSetting::EAlgorithm::kLZ4 + 4);
    if (!fOut || fOut->IsZombie()) {
        std::cerr << "Cannot open output file\n";
        return;
    }

    // Clone structure (no entries), set large basket buffer for throughput
    TTree* tOut = tIn->CloneTree(0);
    tOut->SetAutoSave(1000000000LL);  // AutoSave every ~1 GB written
    tOut->SetBasketSize("*", 64000);  // 64 kB baskets

    // --- Loop ---
    Long64_t nGood = 0, nBad = 0;
    auto t0 = std::chrono::steady_clock::now();

    for (Long64_t i = 0; i < nEntries; ++i) {

        // Progress every 5M entries
        if (i % 5000000 == 0 && i > 0) {
            auto dt = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - t0).count();
            double rate = i / dt / 1e6;
            std::printf("  [%lld / %lld]  good=%lld  bad=%lld  rate=%.1f Mev/s\n",
                        i, nEntries, nGood, nBad, rate);
        }

        try {
            Long64_t bytes = tIn->GetEntry(i);
            if (bytes <= 0) {          // 0 = no data, -1 = I/O error
                ++nBad;
                continue;
            }
            tOut->Fill();
            ++nGood;
        } catch (...) {
            ++nBad;
        }
    }

    // --- Write ---
    fOut->cd();
    tOut->Write("", TObject::kOverwrite);
    fOut->Close();
    fIn->Close();

    auto dt = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - t0).count();

    std::printf("\nDone in %.1f s\n", dt);
    std::printf("Good: %lld   Bad: %lld   (%.4f%% corrupted)\n",
                nGood, nBad, 100.0 * nBad / nEntries);
}

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Usage: filterTree <input.root> <output.root> [treeName]\n";
        return 1;
    }
    const char* tree = (argc >= 4) ? argv[3] : "your_tree_name";
    filterTree(argv[1], argv[2], tree);
    return 0;
}