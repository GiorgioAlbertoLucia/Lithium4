#pragma once

#include "ROOT/RDataFrame.hxx"
#include "Math/Boost.h"
#include "Math/Vector4D.h"

#include <iostream>
#include <map>
#include <string>
#include <vector>

using std::map;
using std::string;
using std::vector;

enum class species {
  kPr = 0,
  kHe = 1,
  kPi = 2,
  kNspecies
};

namespace parametrisation
{
    /// Default values
    std::array<float, 4> kDCAxyResolutionParams[static_cast<int>(species::kNspecies)] = {
      // mean, par0, par1, par2
      {8.19e-5, 0.004, 0.0026, 1.1741},    // Pr
      {1.09e-4, 0.0011, 0.0065, 1.0399}       // He
    };
    std::array<float, 4> kDCAzResolutionParams[static_cast<int>(species::kNspecies)] = {
      // mean, par0, par1, par2
      {1.18e-4, 0.0020, 0.0025, 1.3460},   // Pr
      {9.36e-5, 0.0019, 0.0080, 1.416},  // He
    };

    std::array<float, 4> kDCAxyResolutionParamsMC[static_cast<int>(species::kNspecies)] = {
      // mean, par0, par1, par2
      {-3.60e-5, 0.0013, 0.0021, 1.4722},    // Pr
      {-4.54e-5, 0.0015, 0.0081, 1.6000}       // He
    };
    std::array<float, 4> kDCAzResolutionParamsMC[static_cast<int>(species::kNspecies)] = {
      // mean, par0, par1, par2
      {-1.86e-5, 0.0016, 0.0022, 1.5610},   // Pr
      {-2.06e-5, 0.0018, 0.0080, 1.6000},  // He
    };
    
    // pp 24 and 25 pass1
    std::array<float, 5> kHeTPCParams = {-261.62, 0.1859, 1.1907, 0.5000, 3.0748};
    float kHeTPCResolution = 0.045;
    std::array<float, 6> kHeTPCParamsResiduals = {44.784, 0.6374, 0.0715, -13.329, 0.8686, 0.1638};

    std::array<float, 5> kHeTPCParamsMC = {-283.35, 0.2261, 1.2371, 0.6793, 2.6862};
    float kHeTPCResolutionMC = 0.050;

    std::array<float, 5> kPrTPCParams = {-13.9261, -2.4156, 0.9469, 1.2628, 3.4574};
    float kPrTPCResolution = 0.076;

    // not computed yet, using the same as data
    std::array<float, 5> kPrTPCParamsMC = {-13.9261, -2.4156, 0.9469, 1.2628, 3.4574};
    float kPrTPCResolutionMC = 0.076;

    // pp 24 and 25 pass1
    std::array<float, 3> kITSParams[static_cast<int>(species::kNspecies)] = {
      {1.0228, 1.9634, 2.2081},  // Pr for pp 23
      {2.6692, 1.5112, 5.0994}  // He 
    };
    std::array<float, 3> kITSParamsMC[static_cast<int>(species::kNspecies)] = {
      {1.5273, 1.8027, 2.8959},  // Pr 
      {2.4882, 1.7105, 7.3087}   // He
    };
    
    std::array<float, 3> kITSResolutionParams[static_cast<int>(species::kNspecies)] = {
      {0.1575, 0., 0.},          // Pr, constant
      {0.1132, -0.0135, 0.0001}  // He
    };

    std::array<float, 3> kITSResolutionParamsMC[static_cast<int>(species::kNspecies)] = {
      {0.1300, 0., 0.},          // Pr, constant
      {0.0885, 0.0015, -0.0010}  // He
    };

    std::array<float, 3> kPrTOFParams = {0.9452, -0.0098, 0.0039};
    std::array<float, 2> kPrTOFResolutionParams = {-0.0069, 0.0289};

    std::array<float, 2> kHePidTrkParamsPt = {0.1593, -0.0445};
    std::array<float, 2> kHePidTrkParamsP = {0.1113, -0.0208};

    /// Load parameters

    void SetDCAxyResolutionParams(const int iSpecies, const float mean, const float p0, const float p1, const float p2) {
        kDCAxyResolutionParams[iSpecies] = {mean, p0, p1, p2};
    }
    void SetDCAzResolutionParams(const int iSpecies, const float mean, const float p0, const float p1, const float p2) {
        kDCAzResolutionParams[iSpecies] = {mean, p0, p1, p2};
    }
    void SetDCAxyResolutionParamsMC(const int iSpecies, const float mean, const float p0, const float p1, const float p2) {
        kDCAxyResolutionParamsMC[iSpecies] = {mean, p0, p1, p2};
    }
    void SetDCAzResolutionParamsMC(const int iSpecies, const float mean, const float p0, const float p1, const float p2) {
        kDCAzResolutionParamsMC[iSpecies] = {mean, p0, p1, p2};
    }
    void SetHeTPCParams(const float p0, const float p1, const float p2, const float p3, const float p4) {
        kHeTPCParams = {p0, p1, p2, p3, p4};
    }
    void SetHeTPCResolution(const float res) { kHeTPCResolution = res; }
    void SetHeTPCParamsResiduals(const float p0, const float p1, const float p2, const float p3, const float p4, const float p5) {
        kHeTPCParamsResiduals = {p0, p1, p2, p3, p4, p5};
    }
    void SetHeTPCParamsMC(const float p0, const float p1, const float p2, const float p3, const float p4) {
        kHeTPCParamsMC = {p0, p1, p2, p3, p4};
    }
    void SetHeTPCResolutionMC(const float res) { kHeTPCResolutionMC = res; }
    void SetPrTPCParams(const float p0, const float p1, const float p2, const float p3, const float p4) {
        kPrTPCParams = {p0, p1, p2, p3, p4};
    }
    void SetPrTPCResolution(const float res) { kPrTPCResolution = res; }
    void SetITSParams(const int iSpecies, const float p0, const float p1, const float p2) {
        kITSParams[iSpecies] = {p0, p1, p2};
    }
    void SetITSParamsMC(const int iSpecies, const float p0, const float p1, const float p2) {
        kITSParamsMC[iSpecies] = {p0, p1, p2};
    }
    void SetITSResolutionParams(const int iSpecies, const float p0, const float p1, const float p2) {
        kITSResolutionParams[iSpecies] = {p0, p1, p2};
    }
    void SetITSResolutionParamsMC(const int iSpecies, const float p0, const float p1, const float p2) {
        kITSResolutionParamsMC[iSpecies] = {p0, p1, p2};
    }
    void SetPrTOFParams(const float p0, const float p1, const float p2) {
        kPrTOFParams = {p0, p1, p2};
    }
    void SetPrTOFResolutionParams(const float p0, const float p1) {
        kPrTOFResolutionParams = {p0, p1};
    }
    void SetHePidTrkParamsPt(const float p0, const float p1) {
        kHePidTrkParamsPt = {p0, p1};
    }
    void SetHePidTrkParamsP(const float p0, const float p1) {
        kHePidTrkParamsP = {p0, p1};
    }
}