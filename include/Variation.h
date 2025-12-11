/*
    Functions for the Vary() method of the RDataFrames used in systematics evaluation
*/

#include <stdio.h>
#include <array>
#include <ROOT/RVec.hxx>
#include <TRandom3.h>

namespace VariationRanges 
{
    std::array<float, 2> nSigmaItsHe3Min = {-1.5, -1.};
    std::array<float, 2> absNsigmaTpcHe3Max = {1.5, 2.5};
    
    std::array<float, 2> nSigmaItsPrMin = {-2., -1.5};
    std::array<float, 2> absNsigmaTpcPrMax = {1.5, 2.5};
    std::array<float, 2> absNsigmaTofPrMax = {1.5, 2.5};

    std::array<float, 2> absAmax = {-2., 2};
};

ROOT::RVec<ROOT::RVecF> systematicCuts(const int nCuts, const float nSigmaItsHe3, const float absNsigmaTpcHe3,
                                       const float nSigmaItsPr, const float absNsigmaTpcPr, const float absNsigmaTofPr)
{
    thread_local TRandom3 rng;
    
    ROOT::RVec<ROOT::RVecF> variationsVectors(5);
    for (auto& variableVector : variationsVectors)
        variableVector.reserve(nCuts);

    // ordered variation
    //for (int icut = 0; icut < nCuts; ++icut) {
    //    variationsVectors[0].push_back(nSigmaItsHe3 - (VariationRanges::nSigmaItsHe3Min[1] -
    //                                    icut*(VariationRanges::nSigmaItsHe3Min[1] - VariationRanges::nSigmaItsHe3Min[0])/nCuts) );
    //    variationsVectors[1].push_back(absNsigmaTpcHe3 - (VariationRanges::absNsigmaTpcHe3Max[1] -
    //                                    icut*(VariationRanges::absNsigmaTpcHe3Max[1] - VariationRanges::absNsigmaTpcHe3Max[0])/nCuts) );
    //    variationsVectors[2].push_back(nSigmaItsPr - (VariationRanges::nSigmaItsPrMin[1] -
    //                                    icut*(VariationRanges::nSigmaItsPrMin[1] - VariationRanges::nSigmaItsPrMin[0])/nCuts) );
    //    variationsVectors[3].push_back(absNsigmaTpcPr - (VariationRanges::absNsigmaTpcPrMax[1] -
    //                                    icut*(VariationRanges::absNsigmaTpcPrMax[1] - VariationRanges::absNsigmaTpcPrMax[0])/nCuts) );
    //    variationsVectors[4].push_back(absNsigmaTofPr - (VariationRanges::absNsigmaTofPrMax[1] -
    //                                    icut*(VariationRanges::absNsigmaTofPrMax[1] - VariationRanges::absNsigmaTofPrMax[0])/nCuts) );
    //}

    // random variation
    for (int icut = 0; icut < nCuts; ++icut) {
        variationsVectors[0].push_back(nSigmaItsHe3 - gRandom->Uniform(VariationRanges::nSigmaItsHe3Min[0], VariationRanges::nSigmaItsHe3Min[1]));
        variationsVectors[1].push_back(absNsigmaTpcHe3 - gRandom->Uniform(VariationRanges::absNsigmaTpcHe3Max[0], VariationRanges::absNsigmaTpcHe3Max[1]));
        variationsVectors[2].push_back(nSigmaItsPr - gRandom->Uniform(VariationRanges::nSigmaItsPrMin[0], VariationRanges::nSigmaItsPrMin[1]));
        variationsVectors[3].push_back(absNsigmaTpcPr - gRandom->Uniform(VariationRanges::absNsigmaTpcPrMax[0], VariationRanges::absNsigmaTpcPrMax[1]));
        variationsVectors[4].push_back(absNsigmaTofPr - gRandom->Uniform(VariationRanges::absNsigmaTofPrMax[0], VariationRanges::absNsigmaTofPrMax[1]));
    }

    return variationsVectors;
}


// functions used to check the method
ROOT::RVec<ROOT::RVecD> systematicCutsDummyTest(const double a) {
    
    const size_t nVariables = 1;
    ROOT::RVec<ROOT::RVecD> variationsVector(nVariables);
    variationsVector[0] = {a + VariationRanges::absAmax[0], a + VariationRanges::absAmax[1]};
    return variationsVector;
}

ROOT::RVec<ROOT::RVecD> systematicCutsDummy(const int nCuts, const double a)
{
    thread_local TRandom3 rng;
    
    ROOT::RVec<ROOT::RVecD> variationsVector(1);
    variationsVector[0].reserve(nCuts);  // Reserve space first

    for (int icut = 0; icut < nCuts; ++icut) {
        //variationsVector[0].push_back(a + rng.Uniform(VariationRanges::absAmax[0], VariationRanges::absAmax[1]));
        variationsVector[0].push_back(a + icut*(VariationRanges::absAmax[1] - VariationRanges::absAmax[0])/nCuts );
    }

    return variationsVector;
}