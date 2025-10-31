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

using std::map;
using std::string;
using std::vector;

enum class species {
  kPr = 0,
  kHe = 1,
  kNspecies
};

namespace constant
{
  const float kMass[static_cast<int>(species::kNspecies)] = {0.938272, 2.8089}; // GeV
}

namespace parametrisation
{
    std::array<float, 3> kDCAxyResolutionParams[static_cast<int>(species::kNspecies)] = {
      {0.0032, 0.5206, 0.0012}, // Pr
      {0.0118, 0.6889, 0.0017}     // He
    };
    std::array<float, 3> kDCAzResolutionParams[static_cast<int>(species::kNspecies)] = {
      {0.0021, 1.1122, 0.0021},    // Pr
      {0.1014, 1.7512, 0.0024}     // He
    };
    
    std::array<float, 5> kHeTPCParams = {-178.17, 0.2942, 2.0095, 1.6669, 3.4239};
    float kHeTPCResolution = 0.063;

    std::array<float, 3> kITSParams[static_cast<int>(species::kNspecies)] = {
      {1.0228, 1.9634, 2.2081},  // Pr from Ka fitting
      //{0.6389, 3.4378, 2.3707},  // Pr
      {2.6916, 1.2630, 4.8939}   // He
    };
    
    std::array<float, 3> kITSResolutionParams[static_cast<int>(species::kNspecies)] = {
      //{0.1564, 0.6588, 0.5633}, // Pr
      {0.1575, 0., 0.},          // Pr, constant
      {0.1132, -0.0135, 0.0001}  // He
    };

    std::array<float, 3> kPrTOFParams = {0.9443, -0.0101, 0.0037};
    std::array<float, 2> kPrTOFResolutionParams = {-0.0059, 0.0302};

    std::array<float, 2> kHePidTrkParams = {0.1593, -0.0445};
}

// -------------------------------------------- DCA ----------------------------------------------------

float ComputeNsigmaDCA(const float pt, const float dca, const int iSpecies, const char * dcaType = "xy") {
  
  std::array<float, 3> parameters;
  if (std::strcmp(dcaType, "xy") == 0) {
    parameters = parametrisation::kDCAxyResolutionParams[iSpecies];
  } else if (std::strcmp(dcaType, "z") == 0) {
    parameters = parametrisation::kDCAxyResolutionParams[iSpecies];
  } else {
    std::cout << "Invalid dcaType. Accepted types are 'xy' 'z'" << std::endl;
    parameters = {0., 0., 0.};
  }
  const float sigma = parameters[0] *
                      std::exp(- std::abs(pt) * parameters[1]) +
                      parameters[2];
  return dca / sigma;
}

float ComputeNsigmaDCAxyHe(const float pt, const float dcaxy) {
  return ComputeNsigmaDCA(pt, dcaxy, static_cast<int>(species::kHe), "xy");
}

float ComputeNsigmaDCAzHe(const float pt, const float dcaxy) {
  return ComputeNsigmaDCA(pt, dcaxy, static_cast<int>(species::kHe), "z");
}

float ComputeNsigmaDCAxyPr(const float pt, const float dcaxy) {
  return ComputeNsigmaDCA(pt, dcaxy, static_cast<int>(species::kPr), "xy");
}

float ComputeNsigmaDCAzPr(const float pt, const float dcaxy) {
  return ComputeNsigmaDCA(pt, dcaxy, static_cast<int>(species::kPr), "z");
}

// -------------------------------------------- TPC ----------------------------------------------------

double BetheBlochParametrisation(double bg, double kp1, double kp2, double kp3, double kp4, double kp5) {
  double beta = bg / std::sqrt(1. + bg * bg);
  double aa = std::pow(beta, kp4);
  double bb = std::pow(1. / bg, kp5);
  bb = std::log(kp3 + bb);
  return (kp2 - aa - bb) * kp1 / aa;
}

float BetheBlochHe(const float momentum)  {
  float betagamma = std::abs(momentum) / constant::kMass[static_cast<int>(species::kHe)];
  return BetheBlochParametrisation(betagamma, parametrisation::kHeTPCParams[0], parametrisation::kHeTPCParams[1],
                                   parametrisation::kHeTPCParams[2], parametrisation::kHeTPCParams[3],
                                   parametrisation::kHeTPCParams[4]);
}

float ComputeNsigmaTPCHe(const float momentum, const float tpcSignal) {
  return (tpcSignal / BetheBlochHe(std::abs(momentum)) - 1.) / parametrisation::kHeTPCResolution;
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


float ComputeExpectedClusterSizeCosLambda(const float momentum, const int iSpecies) { 
  const float mass = constant::kMass[iSpecies];
  const std::array<float, 3>& parameters = parametrisation::kITSParams[iSpecies];
  const float betagamma = std::abs(momentum) / mass;
  return parameters[0] / std::pow(betagamma, parameters[1]) + parameters[2];
}

float ComputeExpectedClusterSizeCosLambdaHe(const float momentum) {
  return ComputeExpectedClusterSizeCosLambda(momentum, static_cast<int>(species::kHe));
}

float ComputeExpectedClusterSizeCosLambdaPr(const float momentum) {
  return ComputeExpectedClusterSizeCosLambda(momentum, static_cast<int>(species::kPr));
}


float ComputeClusterSizeResolutionHe(const float momentum)  {
  const float mass = constant::kMass[static_cast<int>(species::kHe)];
  const float betagamma = std::abs(momentum) / mass;
  const std::array<float, 3>& parameters = parametrisation::kITSResolutionParams[static_cast<int>(species::kHe)];
  return parameters[0] + betagamma * parameters[1] + betagamma * betagamma * parameters[2];
}

float ComputeClusterSizeResolutionPr(const float momentum)  {
  const std::array<float, 3>& parameters = parametrisation::kITSResolutionParams[static_cast<int>(species::kPr)];
  //const float mass = constant::kMass[static_cast<int>(species::kPr)];
  //const float betagamma = std::abs(momentum) / mass;
  //return parameters[0] * TMath::Erf((betagamma - parameters[1]) / parameters[2]);
  return parameters[0]; // constant
}

float ComputeClusterSizeResolution(const float momentum, const int iSpecies)  { 
  if (iSpecies == static_cast<int>(species::kPr)) {
    return ComputeClusterSizeResolutionPr(momentum);
  } else if (iSpecies == static_cast<int>(species::kHe)) {
    return ComputeClusterSizeResolutionHe(momentum);
  } else {
    std::cout << "Invalid species" << std::endl;
  }
  return 1.;
}


float ComputeNsigmaITS(const float momentum, const float averageClusterSizeCosLambda, const int iSpecies) { 
  const float mass = constant::kMass[iSpecies];
  const std::array<float, 3>& parameters = parametrisation::kITSParams[iSpecies];

  const float betagamma = std::abs(momentum) / mass;
  const float expected = ComputeExpectedClusterSizeCosLambda(momentum, iSpecies);
  const float resolution = ComputeClusterSizeResolution(momentum, iSpecies);
  return (averageClusterSizeCosLambda - expected) / (resolution * expected);
}

float ComputeNsigmaITSHe(const float momentum, const float averageClusterSizeCosLambda) {
  return ComputeNsigmaITS(momentum, averageClusterSizeCosLambda, static_cast<int>(species::kHe));
}

float ComputeNsigmaITSPr(const float momentum, const float averageClusterSizeCosLambda) {
  return ComputeNsigmaITS(momentum, averageClusterSizeCosLambda, static_cast<int>(species::kPr));
}

// -------------------------------------------- TOF ----------------------------------------------------

float ComputeNsigmaTOFPr(const float pt, const float tofMass) {
  const float expected = parametrisation::kPrTOFParams[0] + parametrisation::kPrTOFParams[1]*std::abs(pt) + parametrisation::kPrTOFParams[2]*std::abs(pt)*std::abs(pt);
  const float resolution = parametrisation::kPrTOFResolutionParams[0] + std::abs(pt) * parametrisation::kPrTOFResolutionParams[1];
  return (tofMass - expected) / (resolution * expected);
}

// -------------------------------------- PID in Tracking ----------------------------------------------

float CorrectPidTrkHe(const float momentum) {
    return momentum * (1. - parametrisation::kHePidTrkParams[0] -
           parametrisation::kHePidTrkParams[1] * momentum);
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
