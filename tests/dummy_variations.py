import ROOT
from ROOT import RDataFrame, gRandom, TTree, TFile, gInterpreter
from array import array

gInterpreter.ProcessLine(f'#include "include/Variation.h"')

if __name__ == '__main__':

    file = TFile.Open("dummy.root", "recreate")
    tree = TTree("tree", "example")

    branch_a = array('d', [0.])
    branch_b = array('d', [0.])
    tree.Branch("a", branch_a, "a/D")
    tree.Branch("b", branch_b, "b/D")

    for _ in range(1_000_000):
        branch_a[0] = gRandom.Gaus(0, 1)
        branch_b[0] = gRandom.Gaus(5, 2)
        tree.Fill()

    tree.Write()

    print("Testing direct function call...")
    result = ROOT.systematicCutsDummyTest(5.0)
    print(f"Result: {result}")
    print(f"Result[0]: {result[0]}")
    
    rdf = RDataFrame(tree)

    # variations will be:    a > 0, a > -0.5;
    # they will be performed requiring a_dummy > 0
    
    # inline definition
    #nominal_hist = rdf.Define('a_dummy', 'a') \
    #                  .Define('b_dummy', 'b') \
    #                  .Vary(['a_dummy', 'b_dummy'], 'ROOT::RVec<ROOT::RVecD>{{a_dummy, a_dummy+0.5}, {b_dummy, b_dummy+0.5}}', 2, "ab") \
    #                  .Filter('a_dummy > 0') \
    #                  .Filter('b_dummy > 0') \
    #                  .Histo1D('a', )

    # importing function
    n_variations = 100
    #systematicCutsDummy = buildSystematicCutsDummy(n_variations)
    nominal_hist = rdf.Define('a_dummy', 'a') \
                      .Vary(['a_dummy'], f'systematicCutsDummy({n_variations}, a_dummy)', n_variations, 'a_variations') \
                      .Filter('a_dummy > 0') \
                      .Histo1D(('a_dummy', '', 100, -5, 5), 'a_dummy')

    #nominal_hist = rdf.Define('a_dummy', 'a') \
    #              .Define('b_dummy', 'b') \
    #              .Vary(['a_dummy'], 'ROOT::RVec<ROOT::RVecD>{{a_dummy + 0.5, a_dummy + 1.0}}', 2, 'a_variations') \
    #              .Filter('a_dummy > 0') \
    #              .Histo1D('a')

    hists = ROOT.RDF.Experimental.VariationsFor(nominal_hist)

    print(f'{hists=}')
    for hist_name in hists.GetKeys():
        hists[hist_name].Write()
    
    file.Close()    
