/*
 * merge_trees.C
 *
 * Horizontal merge of trees with different names within the same DF_ folders.
 * Uses TTree::CopyEntries with "fast" option for basket-level copying
 * (no decompression/recompression for existing branches).
 * Adds a constant fIs23 branch.
 *
 * Usage:
 *   root -l -b -q merge_trees.C
 */

#include <TFile.h>
#include <TTree.h>
#include <TKey.h>
#include <TDirectory.h>
#include <TBranch.h>
#include <TFileMerger.h>
#include <iostream>
#include <vector>
#include <string>

// ---------------------------------------------------------------------------
// Process one input file: for each DF_ folder, merge tree_names horizontally.
// Strategy:
//   1. CloneTree(0) from tree[0] to get the output structure
//   2. CopyEntries("fast") from each input tree into the output tree
//      The "fast" option copies baskets without decompression.
//   3. Fill the constant fIs23 branch separately.
// ---------------------------------------------------------------------------
void process_file(const std::string& infile_path,
                  const std::string& outfile_path,
                  const std::vector<std::string>& tree_names,
                  bool is23)
{
    TFile* infile = TFile::Open(infile_path.c_str(), "READ");
    if (!infile || infile->IsZombie()) {
        std::cerr << "[ERROR] Cannot open " << infile_path << std::endl;
        return;
    }

    TFile* outfile = TFile::Open(outfile_path.c_str(), "RECREATE", "",
                                 ROOT::RCompressionSetting::EDefaults::kUseGeneralPurpose);
    if (!outfile || outfile->IsZombie()) {
        std::cerr << "[ERROR] Cannot create " << outfile_path << std::endl;
        infile->Close();
        return;
    }

    TIter next(infile->GetListOfKeys());
    TKey* key = nullptr;

    while ((key = (TKey*)next())) {
        std::string key_name = key->GetName();
        if (key_name.find("DF_") == std::string::npos)
            continue;

        TDirectory* indir = (TDirectory*)infile->Get(key_name.c_str());
        if (!indir) continue;

        // Load all input trees
        std::vector<TTree*> trees;
        for (const auto& tname : tree_names) {
            TTree* t = (TTree*)indir->Get(tname.c_str());
            if (!t || t->IsZombie()) {
                std::cerr << "[WARN] Tree " << tname << " not found in " << key_name << ", skipping folder." << std::endl;
                trees.clear();
                break;
            }
            trees.push_back(t);
        }
        if (trees.empty()) continue;

        Long64_t nentries = trees[0]->GetEntries();
        for (size_t i = 1; i < trees.size(); i++) {
            if (trees[i]->GetEntries() != nentries) {
                std::cerr << "[WARN] Entry count mismatch in " << key_name << ", skipping." << std::endl;
                trees.clear();
                break;
            }
        }
        if (trees.empty()) continue;

        std::cout << "[INFO] Processing " << key_name << "  (" << nentries << " entries)" << std::endl;

        TDirectory* outdir = outfile->mkdir(key_name.c_str());
        outdir->cd();

        // Step 1: clone structure of all trees into one output tree
        // CloneTree(0) copies only the structure, not the data
        TTree* outtree = trees[0]->CloneTree(0);
        outtree->SetName(tree_names[0].c_str());
        for (size_t i = 1; i < trees.size(); i++) {
            // Temporarily add friend to expose branches during CloneTree structure copy
            TList friend_list;
            TObjArray* src_branches = trees[i]->GetListOfBranches();
            for (int b = 0; b < src_branches->GetEntries(); b++) {
                TBranch* src_br = (TBranch*)src_branches->At(b);
                // Use branch title (leaf list) to recreate branch in outtree
                outtree->Branch(src_br->GetName(), nullptr, src_br->GetTitle());
            }
        }

        // Step 2: CopyEntries with "fast" for each source tree
        // "fast" = raw basket copy, no decompression
        for (size_t i = 0; i < trees.size(); i++) {
            // Activate only the branches belonging to trees[i]
            outtree->SetBranchStatus("*", 0);
            TObjArray* src_branches = trees[i]->GetListOfBranches();
            for (int b = 0; b < src_branches->GetEntries(); b++) {
                TBranch* src_br = (TBranch*)src_branches->At(b);
                outtree->SetBranchStatus(src_br->GetName(), 1);
            }
            trees[i]->CopyEntries(outtree, -1, "fast");
        }
        outtree->SetBranchStatus("*", 1);

        // Step 3: add constant fIs23 branch
        Bool_t fIs23_val = (Bool_t)is23;
        TBranch* is23_branch = outtree->Branch("fIs23", &fIs23_val, "fIs23/O");
        for (Long64_t i = 0; i < nentries; i++)
            is23_branch->Fill();

        outtree->Write("", TObject::kOverwrite);
    }

    outfile->Close();
    infile->Close();
}

// ---------------------------------------------------------------------------
// Vertically merge per-file tmp outputs into the final file via TFileMerger
// ---------------------------------------------------------------------------
void merge_files(const std::vector<std::string>& tmp_files,
                 const std::string& outfile_path)
{
    std::cout << "\n[INFO] Final merge into " << outfile_path << std::endl;

    TFileMerger merger(kFALSE);
    merger.OutputFile(outfile_path.c_str(), "RECREATE");
    merger.SetFastMethod(kTRUE);
    for (const auto& f : tmp_files)
        merger.AddFile(f.c_str());

    if (!merger.Merge())
        std::cerr << "[ERROR] Merge failed." << std::endl;
    else
        std::cout << "[INFO] Done: " << outfile_path << std::endl;
}

void merge_trees()
{
    const std::vector<std::string> tree_names = { "O2he3hadtable", "O2he3hadmult" };

    std::vector<std::string> inputs = {
        "/data/galucia/lithium_local/same/LHC23_PbPb_pass4_hadronpid_same.root",
        //"/data/galucia/lithium_local/same/LHC24ar_pass1_hadronpid_same.root",
        //"/data/galucia/lithium_local/same/LHC24as_pass1_hadronpid_same.root",
    };

    const std::string outfile_path = "/data/galucia/lithium_local/same_merged/LHC23_PbPb_pass4_hadronpid_same_merged.root";

    std::vector<std::string> tmp_files;
    for (const auto& inpath : inputs) {
        bool is23 = (inpath.find("LHC23") != std::string::npos);
        std::string tmp_path = inpath + ".tmp_merged.root";
        std::cout << "\n[INFO] Processing: " << inpath << "  (fIs23=" << is23 << ")" << std::endl;
        process_file(inpath, tmp_path, tree_names, is23);
        tmp_files.push_back(tmp_path);
    }

    merge_files(tmp_files, outfile_path);

    for (const auto& f : tmp_files) {
        std::remove(f.c_str());
        std::cout << "[INFO] Removed tmp: " << f << std::endl;
    }
}