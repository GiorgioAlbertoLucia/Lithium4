#pragma once

#include "ROOT/RDataFrame.hxx"
#include "Math/Boost.h"
#include "Math/Vector4D.h"

#include <iostream>
#include <map>
#include <string>
#include <vector>

#include <TList.h>
#include <TF1.h>
#include <TObject.h>
#include <TMath.h>
#include <TSystem.h>

#include "Parametrisation.h"

using std::map;
using std::string;
using std::vector;

namespace constant
{
  const float kMass[static_cast<int>(species::kNspecies)] = {0.938272, 2.8089, 0.13957}; // GeV
}


// -------------------------------------------- DCA ----------------------------------------------------

float ComputeNsigmaDCA(const float pt, const float dca, const int iSpecies, const char * dcaType = "xy", const bool isMC = false) {
  
  //std::array<float, 3> parameters;
  std::array<float, 4> parameters;
  if (std::strcmp(dcaType, "xy") == 0) {
    parameters = isMC ? parametrisation::kDCAxyResolutionParamsMC[iSpecies] : parametrisation::kDCAxyResolutionParams[iSpecies];
  } else if (std::strcmp(dcaType, "z") == 0) {
    parameters = isMC ? parametrisation::kDCAzResolutionParamsMC[iSpecies] : parametrisation::kDCAzResolutionParams[iSpecies];
  } else {
    std::cout << "Invalid dcaType. Accepted types are 'xy' 'z'" << std::endl;
    //parameters = {0., 0., 0.};
    parameters = {0., 0., 0., 0.};
  }

  //const float sigma = parameters[0] *
  //                    std::exp(- std::abs(pt) * parameters[1]) +
  //                    parameters[2];
  //return dca / sigma;

  //const float sigma = std::sqrt(parameters[1] * parameters[1] +
  //                             parameters[2] * parameters[2] / (pt * pt) +
  //                             parameters[3] * parameters[3] / (pt * pt * pt * pt));
  
  const float sigma = parameters[1] +
                    parameters[2]/std::pow(std::abs(pt), parameters[3]);
  const float mean = parameters[0];
  return (dca - mean) / sigma;
  
}

float ComputeNsigmaDCAxyHe(const float pt, const float dcaxy, const bool isMC = false) {
  return ComputeNsigmaDCA(pt, dcaxy, static_cast<int>(species::kHe), "xy", isMC);
}

float ComputeNsigmaDCAzHe(const float pt, const float dcaz, const bool isMC = false) {
  return ComputeNsigmaDCA(pt, dcaz, static_cast<int>(species::kHe), "z", isMC);
}

float ComputeNsigmaDCAxyPr(const float pt, const float dcaxy, const bool isMC = false) {
  return ComputeNsigmaDCA(pt, dcaxy, static_cast<int>(species::kPr), "xy", isMC);
}

float ComputeNsigmaDCAzPr(const float pt, const float dcaz, const bool isMC = false) {
  return ComputeNsigmaDCA(pt, dcaz, static_cast<int>(species::kPr), "z", isMC);
}

// -------------------------------------------- TPC ----------------------------------------------------

double BetheBlochParametrisation(double bg, double kp1, double kp2, double kp3, double kp4, double kp5) {
  double beta = bg / std::sqrt(1. + bg * bg);
  double aa = std::pow(beta, kp4);
  double bb = std::pow(1. / bg, kp5);
  bb = std::log(kp3 + bb);
  return (kp2 - aa - bb) * kp1 / aa;
}

float BetheBlochHe(const float momentum, const bool isMC=false, const bool applyResiduals=false) {
  float betagamma = std::abs(momentum) / constant::kMass[static_cast<int>(species::kHe)];
  std::array<float, 5> params = isMC ? parametrisation::kHeTPCParamsMC : parametrisation::kHeTPCParams;

  double dEdxExpected = BetheBlochParametrisation(betagamma, params[0], params[1], params[2], params[3], params[4]);

  if (applyResiduals) {
    double residualCorrection = (parametrisation::kHeTPCParamsResiduals[0] * std::exp(- (betagamma - parametrisation::kHeTPCParamsResiduals[1]) * (betagamma - parametrisation::kHeTPCParamsResiduals[1]) / (2 * parametrisation::kHeTPCParamsResiduals[2] * parametrisation::kHeTPCParamsResiduals[2]))
                           + parametrisation::kHeTPCParamsResiduals[3] * std::exp(- (betagamma - parametrisation::kHeTPCParamsResiduals[4]) * (betagamma - parametrisation::kHeTPCParamsResiduals[4]) / (2 * parametrisation::kHeTPCParamsResiduals[5] * parametrisation::kHeTPCParamsResiduals[5])));
    dEdxExpected += residualCorrection;
  }

  return dEdxExpected;
}

float ComputeNsigmaTPCHe(const float momentum, const float tpcSignal, const bool isMC=false, const bool applyResiduals=false) {
  const float resolution = isMC ? parametrisation::kHeTPCResolutionMC : parametrisation::kHeTPCResolution;
  return (tpcSignal / BetheBlochHe(std::abs(momentum), isMC, applyResiduals) - 1.) / resolution;
}

float BetheBlochPi(const float momentum, const bool isMC=false) {
  float betagamma = std::abs(momentum) / constant::kMass[static_cast<int>(species::kPi)];
  std::array<float, 5> params = isMC ? parametrisation::kPrTPCParamsMC : parametrisation::kPrTPCParams;

  return BetheBlochParametrisation(betagamma, params[0], params[1], params[2], params[3], params[4]);
}

float ComputeNsigmaTPCPi(const float momentum, const float tpcSignal, const bool isMC=false) {
  const float resolution = isMC ? parametrisation::kPrTPCResolutionMC : parametrisation::kPrTPCResolution;
  return (tpcSignal / BetheBlochPi(std::abs(momentum), isMC) - 1.) / resolution;
}

// -------------------------------------------- ITS ----------------------------------------------------

float ComputeAverageClusterSize(const uint32_t itsClusterSizes, const bool useTruncatedMean = false)  {
  float sum = 0;
  int nclusters = 0;
  int max = 0;
  for (int layer = 0; layer < 7; layer++) {
    int clsize = (itsClusterSizes >> (layer * 4)) & 0xf;
    if (clsize > 0) {
      nclusters++;
      sum += clsize;
      if (clsize > max) {
        max = clsize;
      }
    }
  }
  
  if (nclusters == 0) {
    return 0;
  }
  
  if (useTruncatedMean && nclusters > 1) {
    return (sum - max) / (nclusters - 1);
  }
  return sum / nclusters;
};

float ComputeExpectedClusterSizeCosLambda(const float momentum, const int iSpecies, const bool isMC=false) { 
  const float mass = constant::kMass[iSpecies];
  const std::array<float, 3>& parameters = isMC ? parametrisation::kITSParamsMC[iSpecies]:
                                               parametrisation::kITSParams[iSpecies];
  const float betagamma = std::abs(momentum) / mass;
  return parameters[0] / std::pow(betagamma, parameters[1]) + parameters[2];
}

float ComputeExpectedClusterSizeCosLambdaHe(const float momentum, const bool isMC=false) {
  return ComputeExpectedClusterSizeCosLambda(momentum, static_cast<int>(species::kHe), isMC);
}

float ComputeExpectedClusterSizeCosLambdaPr(const float momentum) {
  return ComputeExpectedClusterSizeCosLambda(momentum, static_cast<int>(species::kPr));
}


float ComputeClusterSizeResolutionHe(const float momentum, const bool isMC=false)  {
  const float mass = constant::kMass[static_cast<int>(species::kHe)];
  const float betagamma = std::abs(momentum) / mass;
  const std::array<float, 3>& parameters = isMC ? parametrisation::kITSResolutionParamsMC[static_cast<int>(species::kHe)]:
                                               parametrisation::kITSResolutionParams[static_cast<int>(species::kHe)];
  return parameters[0] + betagamma * parameters[1] + betagamma * betagamma * parameters[2];
}

float ComputeClusterSizeResolutionPr(const float momentum, const bool isMC=false)  {
  const std::array<float, 3>& parameters = isMC ? parametrisation::kITSResolutionParamsMC[static_cast<int>(species::kPr)]:
                                               parametrisation::kITSResolutionParams[static_cast<int>(species::kPr)];
  //const float mass = constant::kMass[static_cast<int>(species::kPr)];
  //const float betagamma = std::abs(momentum) / mass;
  //return parameters[0] * TMath::Erf((betagamma - parameters[1]) / parameters[2]);
  return parameters[0]; // constant
}

float ComputeClusterSizeResolution(const float momentum, const int iSpecies, const bool isMC=false)  { 
  if (iSpecies == static_cast<int>(species::kPr)) {
    return ComputeClusterSizeResolutionPr(momentum, isMC);
  } else if (iSpecies == static_cast<int>(species::kHe)) {
    return ComputeClusterSizeResolutionHe(momentum, isMC);
  } else {
    std::cout << "Invalid species" << std::endl;
  }
  return 1.;
}


float ComputeNsigmaITS(const float momentum, const float averageClusterSizeCosLambda, const int iSpecies, const bool isMC=false) { 
  const float mass = constant::kMass[iSpecies];
  
  const float betagamma = std::abs(momentum) / mass;
  const float expected = ComputeExpectedClusterSizeCosLambda(momentum, iSpecies, isMC);
  const float resolution = ComputeClusterSizeResolution(momentum, iSpecies, isMC);
  return (averageClusterSizeCosLambda - expected) / (resolution * expected);
}

float ComputeNsigmaITSHe(const float momentum, const float averageClusterSizeCosLambda, const bool isMC=false) {
  return ComputeNsigmaITS(momentum, averageClusterSizeCosLambda, static_cast<int>(species::kHe), isMC);
}

float ComputeNsigmaITSPr(const float momentum, const float averageClusterSizeCosLambda, const bool isMC=false) {
  return ComputeNsigmaITS(momentum, averageClusterSizeCosLambda, static_cast<int>(species::kPr), isMC);
}

// -------------------------------------------- TOF ----------------------------------------------------

float ComputeNsigmaTOFPr(const float pt, const float tofMass) {
  const float expected = parametrisation::kPrTOFParams[0] + parametrisation::kPrTOFParams[1]*std::abs(pt) + parametrisation::kPrTOFParams[2]*std::abs(pt)*std::abs(pt);
  const float resolution = parametrisation::kPrTOFResolutionParams[0] + std::abs(pt) * parametrisation::kPrTOFResolutionParams[1];
  return (tofMass - expected) / (resolution * expected);
}

// -------------------------------------- PID in Tracking ----------------------------------------------

float CorrectPidTrkHe(const float momentum, const bool isPt = true) {
    const auto& params = isPt ? parametrisation::kHePidTrkParamsPt : parametrisation::kHePidTrkParamsP;
    return momentum * (1. - params[0] - params[1] * momentum);
}

// ------------------------------------------- Femto ---------------------------------------------------

float ComputeKstar(const double pt1, const double eta1, const double phi1, const double m1, const double pt2, const double eta2, const double phi2, const double m2)  {
    using namespace ROOT::Math;

    PtEtaPhiMVector p1mu(pt1, eta1, phi1, m1);
    PtEtaPhiMVector p2mu(pt2, eta2, phi2, m2);
    auto P_beta_vector = (p1mu + p2mu).BoostToCM();
    double P_bx(0), P_by(0), P_bz(0);
    P_beta_vector.GetCoordinates(P_bx, P_by, P_bz);
    auto P_boost = Boost(P_bx, P_by, P_bz);

    PtEtaPhiMVector p1mu_to_boost(pt1, eta1, phi1, m1);
    PtEtaPhiMVector p2mu_to_boost(pt2, eta2, phi2, m2);

    auto p1mu_star = P_boost(p1mu_to_boost);
    auto p2mu_star = P_boost(p2mu_to_boost);

    PtEtaPhiMVector kmu_star = p1mu_star - p2mu_star;
    const double kstar = 0.5 * kmu_star.P();

    return kstar;
}

inline int ReadPidTrkFromFlags(const int flags) { return (flags >> 12) & 0x1F; }

inline int ReadBitFromFlags(const int flags, const int indexBit) { return (flags >> indexBit) & 0b1; }
