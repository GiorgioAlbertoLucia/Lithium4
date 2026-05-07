"""
recover_tree_fast.py

Recovers a corrupted ROOT TTree by:
  1. Scanning for bad entries using a single branch (fast)
  2. Bulk-copying good entry ranges with CopyTree (no entry-by-entry loop)
  3. Merging all chunks with TTree::MergeTrees

Usage:
    python recover_tree_fast.py <input_file> <output_file> [tree_name] [scan_branch]
"""

import sys
import ROOT
from tqdm import tqdm
from ROOT import TFile, TTree


def find_bad_regions(tree, n_total, scan_branch):
    """
    Scans the tree entry-by-entry using only one branch for speed.
    Returns a list of (first_bad, last_bad) inclusive ranges.
    """

    tree.SetBranchStatus("*", 0)
    tree.SetBranchStatus(scan_branch, 1)

    bad_regions = []
    in_bad = False
    bad_start = -1

    print(f"  Scanning with branch '{scan_branch}' only...")
    for i in tqdm(range(n_total)):

        ok = tree.GetEntry(i) > 0
        if not ok and not in_bad:
            bad_start = i
            in_bad = True
        elif ok and in_bad:
            bad_regions.append((bad_start, i - 1))
            print(f"  Bad region found: entries [{bad_start}, {i-1}] ({i - bad_start} entries)")
            in_bad = False

    if in_bad:
        bad_regions.append((bad_start, n_total - 1))
        print(f"  Bad region found: entries [{bad_start}, {n_total-1}] ({n_total - bad_start} entries)")

    tree.SetBranchStatus("*", 1)
    return bad_regions


def recover_tree_fast(input_path: str, output_path: str, tree_name: str = "MixedTree", scan_branch: str = "fEtaHad"):

    ROOT.gROOT.SetBatch(True)
    ROOT.gErrorIgnoreLevel = ROOT.kFatal

    print(f"Opening input file: {input_path}")
    input_file = TFile.Open(input_path, "READ")
    if not input_file or input_file.IsZombie():
        raise RuntimeError(f"Cannot open input file: {input_path}")

    tree = input_file.Get(tree_name)
    if not tree:
        raise RuntimeError(f"Tree '{tree_name}' not found in {input_path}")

    n_total = tree.GetEntries()
    print(f"Tree '{tree_name}' has {n_total} entries total.\n")

    print("Step 1: Locating bad entry regions...")
    bad_regions = find_bad_regions(tree, n_total, scan_branch)

    if not bad_regions:
        print("No bad entries found — tree appears clean.")
        input_file.Close()
        return

    print(f"\nFound {len(bad_regions)} bad region(s).\n")

    # Build good ranges as the complement of bad regions
    good_ranges = []
    cursor = 0
    for bad_start, bad_end in bad_regions:
        if cursor < bad_start:
            good_ranges.append((cursor, bad_start - 1))
        cursor = bad_end + 1
    if cursor < n_total:
        good_ranges.append((cursor, n_total - 1))

    print("Step 2: Bulk-copying good ranges...")
    output_file = TFile.Open(output_path, "RECREATE")
    if not output_file or output_file.IsZombie():
        raise RuntimeError(f"Cannot open output file: {output_path}")
    output_file.cd()

    chunks = []
    n_good = 0
    for (start, end) in good_ranges:
        count = end - start + 1
        print(f"  Copying entries [{start}, {end}] ({count} entries)...")
        chunk = tree.CopyTree("", "", count, start)
        chunks.append(chunk)
        n_good += count

    print(f"\nStep 3: Merging {len(chunks)} chunk(s)...")
    merge_list = ROOT.TList()
    for chunk in tqdm(chunks):
        merge_list.Add(chunk)

    out_tree = ROOT.TTree.MergeTrees(merge_list)
    out_tree.SetName(tree_name)

    ROOT.gErrorIgnoreLevel = ROOT.kInfo

    n_bad = n_total - n_good
    print(f"\nDone.")
    print(f"  Good entries copied : {n_good}")
    print(f"  Bad entries skipped : {n_bad}")
    print(f"  Bad regions         : {bad_regions}")

    output_file.cd()
    out_tree.Write("", ROOT.TObject.kOverwrite)
    output_file.Close()
    input_file.Close()

    print(f"\nRecovered tree written to: {output_path}")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python recover_tree_fast.py <input_file> <output_file> <scan_branch> [tree_name]")
        sys.exit(1)

    input_path = sys.argv[1]
    output_path = sys.argv[2]
    scan_branch = sys.argv[3]
    tree_name = sys.argv[4] if len(sys.argv) > 4 else "MixedTree"

    recover_tree_fast(input_path, output_path, tree_name, scan_branch)