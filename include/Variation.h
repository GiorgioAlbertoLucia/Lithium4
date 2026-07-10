/*
    Functions for the Vary() method of the RDataFrames used in systematics evaluation
*/

#include <stdio.h>
#include <array>
#include <ROOT/RVec.hxx>
#include <TRandom3.h>

namespace VariationRanges 
{
    std::array<float, 2> nSigmaItsHe3Min = {-2.5, -1.5};
    std::array<float, 2> nSigmaTpcHe3Max = {2.5, 3.5};
    std::array<float, 2> nSigmaTpcHe3Min = {-2.5, -1.5};
    std::array<float, 2> absNsigmaDcaXyHe3Max = {3., 7.};
    std::array<float, 2> absNsigmaDcaZHe3Max = {3., 7.};
    std::array<float, 2> nClsTpcHe3Min = {100., 120.};
    std::array<float, 2> chi2TpcHe3Max = {2., 4.};
    
    std::array<float, 2> nSigmaItsPrMin = {-2.5, -1.5};
    std::array<float, 2> absNsigmaTpcPrMax = {1.5, 2.5};
    std::array<float, 2> absNsigmaTofPrMax = {1.5, 2.5};
    std::array<float, 2> absNsigmaDcaXyPrMax = {3., 7.};
    std::array<float, 2> absNsigmaDcaZPrMax = {3., 7.};
    std::array<float, 2> chi2TpcPrMax = {2., 4.}; 
    
    std::array<float, 2> absZVertexMax = {8., 10.};

    const float nominalNsigmaItsHe3 = -1.5;
    const float nominalNsigmaTpcMinHe3 = -2.;
    const float nominalNsigmaTpcMaxHe3 = 3.;
    const float nominalAbsNsigmaDcaXyHe3 = 3;
    const float nominalAbsNsigmaDcaZHe3 = 3;
    const float nominalNclsTpcHe3 = 110;
    const float nominalChi2TpcHe3 = 4.;
    
    const float nominalNsigmaItsPr = -2.;
    const float nominalAbsNsigmaTpcPr = 2.;
    const float nominalAbsNsigmaTofPr = 2.;
    const float nominalAbsNsigmaDcaXyPr = 3;
    const float nominalAbsNsigmaDcaZPr = 3;
    const float nominalChi2TpcPr = 4.;
    
    const float nominalAbsZVertex = 10.;

    // testing
    std::array<float, 2> absAmax = {-2., 2};

    bool debug = true;
};

ROOT::RVec<ROOT::RVecF> systematicCuts(const int nCuts, const float nSigmaItsHe3, const float nSigmaTpcMinHe3, const float nSigmaTpcMaxHe3,
                                       const float absNsigmaDcaXyHe3, const float absNsigmaDcaZHe3,
                                       const float nClsTpcHe3, const float chi2TpcHe3,
                                       const float nSigmaItsPr, const float absNsigmaTpcPr, const float absNsigmaTofPr,
                                       const float absNsigmaDcaXyPr, const float absNsigmaDcaZPr, const float chi2TpcPr,
                                       const float absZVertex)
{
    thread_local TRandom3 rng;
    
    ROOT::RVec<ROOT::RVecF> variationsVectors(14);
    for (auto& variableVector : variationsVectors)
        variableVector.reserve(nCuts);

        // random variation
        for (int icut = 0; icut < nCuts; ++icut) {
        const float offsetNsigmaItsHe3 = rng.Uniform(VariationRanges::nSigmaItsHe3Min[0], VariationRanges::nSigmaItsHe3Min[1]);
        const float offsetNsigmaTpcMinHe3 = rng.Uniform(VariationRanges::nSigmaTpcHe3Min[0], VariationRanges::nSigmaTpcHe3Min[1]);
        const float offsetNsigmaTpcMaxHe3 = rng.Uniform(VariationRanges::nSigmaTpcHe3Max[0], VariationRanges::nSigmaTpcHe3Max[1]);
        const float offsetAbsNsigmaDcaXyHe3 = rng.Uniform(VariationRanges::absNsigmaDcaXyHe3Max[0], VariationRanges::absNsigmaDcaXyHe3Max[1]);
        const float offsetAbsNsigmaDcaZHe3 = rng.Uniform(VariationRanges::absNsigmaDcaZHe3Max[0], VariationRanges::absNsigmaDcaZHe3Max[1]);
        const float offsetNclsTpcHe3 = rng.Uniform(VariationRanges::nClsTpcHe3Min[0], VariationRanges::nClsTpcHe3Min[1]);
        const float offsetChi2TpcHe3 = rng.Uniform(VariationRanges::chi2TpcHe3Max[0], VariationRanges::chi2TpcHe3Max[1]);
        const float offsetNsigmaItsPr = rng.Uniform(VariationRanges::nSigmaItsPrMin[0], VariationRanges::nSigmaItsPrMin[1]);
        const float offsetAbsNsigmaTpcPr = rng.Uniform(VariationRanges::absNsigmaTpcPrMax[0], VariationRanges::absNsigmaTpcPrMax[1]);
        const float offsetAbsNsigmaTofPr = rng.Uniform(VariationRanges::absNsigmaTofPrMax[0], VariationRanges::absNsigmaTofPrMax[1]);
        const float offsetAbsNsigmaDcaXyPr = rng.Uniform(VariationRanges::absNsigmaDcaXyPrMax[0], VariationRanges::absNsigmaDcaXyPrMax[1]);
        const float offsetAbsNsigmaDcaZPr = rng.Uniform(VariationRanges::absNsigmaDcaZPrMax[0], VariationRanges::absNsigmaDcaZPrMax[1]);
        const float offsetChi2TpcPr = rng.Uniform(VariationRanges::chi2TpcPrMax[0], VariationRanges::chi2TpcPrMax[1]);
        const float offsetAbsZVertex = rng.Uniform(VariationRanges::absZVertexMax[0], VariationRanges::absZVertexMax[1]);

        variationsVectors[0].push_back(nSigmaItsHe3 - offsetNsigmaItsHe3);
        variationsVectors[1].push_back(nSigmaTpcMinHe3 - offsetNsigmaTpcMinHe3);
        variationsVectors[2].push_back(nSigmaTpcMaxHe3 - offsetNsigmaTpcMaxHe3);
        variationsVectors[3].push_back(absNsigmaDcaXyHe3 - offsetAbsNsigmaDcaXyHe3);
        variationsVectors[4].push_back(absNsigmaDcaZHe3 - offsetAbsNsigmaDcaZHe3);
        variationsVectors[5].push_back(nClsTpcHe3 - offsetNclsTpcHe3);
        variationsVectors[6].push_back(chi2TpcHe3 - offsetChi2TpcHe3);
        variationsVectors[7].push_back(nSigmaItsPr - offsetNsigmaItsPr);
        variationsVectors[8].push_back(absNsigmaTpcPr - offsetAbsNsigmaTpcPr);
        variationsVectors[9].push_back(absNsigmaTofPr - offsetAbsNsigmaTofPr);
        variationsVectors[10].push_back(absNsigmaDcaXyPr - offsetAbsNsigmaDcaXyPr);
        variationsVectors[11].push_back(absNsigmaDcaZPr - offsetAbsNsigmaDcaZPr);
        variationsVectors[12].push_back(chi2TpcPr - offsetChi2TpcPr);
        variationsVectors[13].push_back(absZVertex - offsetAbsZVertex);
        
        if (VariationRanges::debug) { // DEBUGGING
                printf("nsigmaItsHe3: %f, offset: %f\n", nSigmaItsHe3, offsetNsigmaItsHe3);
                printf("nsigmaTpcMinHe3: %f, offset: %f\n", nSigmaTpcMinHe3, offsetNsigmaTpcMinHe3);
                printf("nsigmaTpcMaxHe3: %f, offset: %f\n", nSigmaTpcMaxHe3, offsetNsigmaTpcMaxHe3);
                printf("absNsigmaDcaXyHe3: %f, offset: %f\n", absNsigmaDcaXyHe3, offsetAbsNsigmaDcaXyHe3);
                printf("absNsigmaDcaZHe3: %f, offset: %f\n", absNsigmaDcaZHe3, offsetAbsNsigmaDcaZHe3);
                printf("nClsTpcHe3: %f, offset: %f\n", nClsTpcHe3, offsetNclsTpcHe3);
                printf("chi2TpcHe3: %f, offset: %f\n",  chi2TpcHe3, offsetChi2TpcHe3);
                printf("nsigmaItsPr: %f, offset: %f\n", nSigmaItsPr, offsetNsigmaItsPr);
                printf("absNsigmaTpcPr: % f, offset: %f\n", absNsigmaTpcPr, offsetAbsNsigmaTpcPr);
                printf("absNsigmaTofPr: % f, offset: %f\n", absNsigmaTofPr, offsetAbsNsigmaTofPr);
                printf("absNsigmaDcaXyPr: % f, offset: %f\n", absNsigmaDcaXyPr, offsetAbsNsigmaDcaXyPr);
                printf("absNsigmaDcaZPr: % f, offset: %f\n", absNsigmaDcaZPr, offsetAbsNsigmaDcaZPr);
                printf("chi2TpcPr: % f, offset: %f\n", chi2TpcPr, offsetChi2TpcPr);
                printf("absZVertex: % f, offset: %f\n", absZVertex, offsetAbsZVertex);
                VariationRanges::debug = false;
            }
    }

    return variationsVectors;
}

ROOT::RVec<ROOT::RVecF> systematicCutsSingleVariable(const int nCuts, const char * variable,
                                            const float nSigmaItsHe3, const float nSigmaTpcMinHe3, const float nSigmaTpcMaxHe3,
                                            const float absNsigmaDcaXyHe3, const float absNsigmaDcaZHe3,
                                            const float nClsTpcHe3,  const float chi2TpcHe3,
                                            const float nSigmaItsPr, const float absNsigmaTpcPr, const float absNsigmaTofPr,
                                            const float absNsigmaDcaXyPr, const float absNsigmaDcaZPr, const float chi2TpcPr,
                                            const float absZVertex)
{
    thread_local TRandom3 rng;
    
    ROOT::RVec<ROOT::RVecF> variationsVectors(14);
    for (auto& variableVector : variationsVectors)
        variableVector.reserve(nCuts);

    // random variation
    for (int icut = 0; icut < nCuts; ++icut) {

        float offsetNsigmaItsHe3 = VariationRanges::nominalNsigmaItsHe3, 
            offsetNsigmaTpcMinHe3 = VariationRanges::nominalNsigmaTpcMinHe3,
            offsetNsigmaTpcMaxHe3 = VariationRanges::nominalNsigmaTpcMaxHe3,
            offsetAbsNsigmaDcaXyHe3 = VariationRanges::nominalAbsNsigmaDcaXyHe3,
            offsetAbsNsigmaDcaZHe3 = VariationRanges::nominalAbsNsigmaDcaZHe3,
            offsetNclsTpcHe3 = VariationRanges::nominalNclsTpcHe3,
            offsetChi2TpcHe3 = VariationRanges::nominalChi2TpcHe3,
            offsetNsigmaItsPr = VariationRanges::nominalNsigmaItsPr, 
            offsetAbsNsigmaTpcPr = VariationRanges::nominalAbsNsigmaTpcPr, 
            offsetAbsNsigmaTofPr = VariationRanges::nominalAbsNsigmaTofPr,
            offsetAbsNsigmaDcaXyPr = VariationRanges::nominalAbsNsigmaDcaXyPr,
            offsetAbsNsigmaDcaZPr = VariationRanges::nominalAbsNsigmaDcaZPr,
            offsetChi2TpcPr = VariationRanges::nominalChi2TpcPr,
            offsetAbsZVertex = VariationRanges::nominalAbsZVertex;

        if (strcmp(variable, "fNSigmaITSHe3") == 0) {
            offsetNsigmaItsHe3 = rng.Uniform(VariationRanges::nSigmaItsHe3Min[0], VariationRanges::nSigmaItsHe3Min[1]);
        } else if (strcmp(variable, "fNSigmaTPCHe3") == 0) {
            offsetNsigmaTpcMinHe3 = rng.Uniform(VariationRanges::nSigmaTpcHe3Min[0], VariationRanges::nSigmaTpcHe3Min[1]);
            offsetNsigmaTpcMaxHe3 = rng.Uniform(VariationRanges::nSigmaTpcHe3Max[0], VariationRanges::nSigmaTpcHe3Max[1]);
        } else if (strcmp(variable, "fNSigmaDCAxyHe3") == 0) {
            offsetAbsNsigmaDcaXyHe3 = rng.Uniform(VariationRanges::absNsigmaDcaXyHe3Max[0], VariationRanges::absNsigmaDcaXyHe3Max[1]);
        } else if (strcmp(variable, "fNSigmaDCAzHe3") == 0) {
            offsetAbsNsigmaDcaZHe3 = rng.Uniform(VariationRanges::absNsigmaDcaZHe3Max[0], VariationRanges::absNsigmaDcaZHe3Max[1]);
        } else if (strcmp(variable, "fNClsTPCHe3") == 0) {
            offsetNclsTpcHe3 = rng.Uniform(VariationRanges::nClsTpcHe3Min[0], VariationRanges::nClsTpcHe3Min[1]);
        } else if (strcmp(variable, "fChi2TPCHe3") == 0) {
            offsetChi2TpcHe3 = rng.Uniform(VariationRanges::chi2TpcHe3Max[0], VariationRanges::chi2TpcHe3Max[1]);
        }
        
        else if (strcmp(variable, "fNSigmaITSHad") == 0) {
            offsetNsigmaItsPr = rng.Uniform(VariationRanges::nSigmaItsPrMin[0], VariationRanges::nSigmaItsPrMin[1]); 
        } else if (strcmp(variable, "fNSigmaTPCHad") == 0) {
            offsetAbsNsigmaTpcPr = rng.Uniform(VariationRanges::absNsigmaTpcPrMax[0], VariationRanges::absNsigmaTpcPrMax[1]);
        } else if (strcmp(variable, "fNSigmaTOFHad") == 0) {
            offsetAbsNsigmaTofPr = rng.Uniform(VariationRanges::absNsigmaTofPrMax[0], VariationRanges::absNsigmaTofPrMax[1]);
        } else if (strcmp(variable, "fNSigmaDCAxyHad") == 0) {
            offsetAbsNsigmaDcaXyPr = rng.Uniform(VariationRanges::absNsigmaDcaXyPrMax[0], VariationRanges::absNsigmaDcaXyPrMax[1]);
        } else if (strcmp(variable, "fNSigmaDCAzHad") == 0) {
            offsetAbsNsigmaDcaZPr = rng.Uniform(VariationRanges::absNsigmaDcaZPrMax[0], VariationRanges::absNsigmaDcaZPrMax[1]);
        } else if (strcmp(variable, "fChi2TPCHad") == 0) {
            offsetChi2TpcPr = rng.Uniform(VariationRanges::chi2TpcPrMax[0], VariationRanges::chi2TpcPrMax[1]);
        } else if (strcmp(variable, "fZVertex") == 0) {
            offsetAbsZVertex = rng.Uniform(VariationRanges::absZVertexMax[0], VariationRanges::absZVertexMax[1]);
        }

        if (VariationRanges::debug) { // DEBUGGING
            printf("nsigmaItsHe3: %f, offset: %f\n", nSigmaItsHe3, offsetNsigmaItsHe3);
            printf("nsigmaTpcMinHe3: %f, offset: %f\n", nSigmaTpcMinHe3, offsetNsigmaTpcMinHe3);
            printf("nsigmaTpcMaxHe3: %f, offset: %f\n", nSigmaTpcMaxHe3, offsetNsigmaTpcMaxHe3);
            printf("absNsigmaDcaXyHe3: %f, offset: %f\n", absNsigmaDcaXyHe3, offsetAbsNsigmaDcaXyHe3);
            printf("absNsigmaDcaZHe3: %f, offset: %f\n", absNsigmaDcaZHe3, offsetAbsNsigmaDcaZHe3);
            printf("nClsTpcHe3: %f, offset: %f\n", nClsTpcHe3, offsetNclsTpcHe3);
            printf("chi2TpcHe3: %f, offset: %f\n",  chi2TpcHe3, offsetChi2TpcHe3);
            printf("nsigmaItsPr: %f, offset: %f\n", nSigmaItsPr, offsetNsigmaItsPr);
            printf("absNsigmaTpcPr: % f, offset: %f\n", absNsigmaTpcPr, offsetAbsNsigmaTpcPr);
            printf("absNsigmaTofPr: % f, offset: %f\n", absNsigmaTofPr, offsetAbsNsigmaTofPr);
            printf("absNsigmaDcaXyPr: % f, offset: %f\n", absNsigmaDcaXyPr, offsetAbsNsigmaDcaXyPr);
            printf("absNsigmaDcaZPr: % f, offset: %f\n", absNsigmaDcaZPr, offsetAbsNsigmaDcaZPr);
            printf("chi2TpcPr: % f, offset: %f\n", chi2TpcPr, offsetChi2TpcPr);
            printf("absZVertex: % f, offset: %f\n", absZVertex, offsetAbsZVertex);
            VariationRanges::debug = false;
        }

        variationsVectors[0].push_back(nSigmaItsHe3 - offsetNsigmaItsHe3);
        variationsVectors[1].push_back(nSigmaTpcMinHe3 - offsetNsigmaTpcMinHe3);
        variationsVectors[2].push_back(nSigmaTpcMaxHe3 - offsetNsigmaTpcMaxHe3);
        variationsVectors[3].push_back(absNsigmaDcaXyHe3 - offsetAbsNsigmaDcaXyHe3);
        variationsVectors[4].push_back(absNsigmaDcaZHe3 - offsetAbsNsigmaDcaZHe3);
        variationsVectors[5].push_back(nClsTpcHe3 - offsetNclsTpcHe3);
        variationsVectors[6].push_back(chi2TpcHe3 - offsetChi2TpcHe3);
        
        variationsVectors[7].push_back(nSigmaItsPr - offsetNsigmaItsPr);
        variationsVectors[8].push_back(absNsigmaTpcPr - offsetAbsNsigmaTpcPr);
        variationsVectors[9].push_back(absNsigmaTofPr - offsetAbsNsigmaTofPr);
        variationsVectors[10].push_back(absNsigmaDcaXyPr - offsetAbsNsigmaDcaXyPr);
        variationsVectors[11].push_back(absNsigmaDcaZPr - offsetAbsNsigmaDcaZPr);
        variationsVectors[12].push_back(chi2TpcPr - offsetChi2TpcPr);

        variationsVectors[13].push_back(absZVertex - offsetAbsZVertex);
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