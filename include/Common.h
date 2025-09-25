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

namespace constant
{
    const float kMassHe = 2.8089; // GeV
    const float kMassPr = 0.938272; // GeV

}

namespace parametrisation
{
    std::array<float, 3> kHeDCAxyResolutionParams = {8.214e-3, 1.253, 1.041e-3};
    std::array<float, 3> kPrDCAxyResolutionParams = {3.755e-3, 7.035e-1, 1.264e-4};
    std::array<float, 3> kHeDCAzResolutionParams = {6.420e-3, 3.485, 5.118e-3};
    std::array<float, 3> kPrDCAzResolutionParams = {2.937e-4, 2.322, 5.121e-3};

    std::array<float, 5> kHeTpcParams = {-251.9, 0.3223, 1.355, 0.6456, 2.675};
    float kHeTpcResolution = 0.07;

    std::array<float, 3> kHeITSParams = {2.35117, 1.80347, 5.14355};
    std::array<float, 3> kHeITSResolutionParams = {8.74371e-02, -1.82804, 5.06449e-01};

    std::array<float, 2> kPrTofResolutionParams = {1.22204e-02, 7.48467e-01};

    std::array<float, 3> kPrITSParams = {1.18941, 1.53792, 1.69961,};
    std::array<float, 3> kPrITSResolutionParams = {1.94669e-01, -2.08616e-01, 1.30753,};

    std::array<float, 2> kHePidTrkParams = {0.1593, -0.0445};
}

// -------------------------------------------- DCA ----------------------------------------------------

float NsigmaDCAxyHe(const float pt, const float dcaxy) {
  const float sigma = parametrisation::kHeDCAxyResolutionParams[0] /
                      std::pow(pt, parametrisation::kHeDCAxyResolutionParams[1]) +
                      parametrisation::kHeDCAxyResolutionParams[2];
  return dcaxy / sigma;
}

float NsigmaDCAxyPr(const float pt, const float dcaxy) {
  const float sigma = parametrisation::kPrDCAxyResolutionParams[0] /
                      std::pow(pt, parametrisation::kPrDCAxyResolutionParams[1]) +
                      parametrisation::kPrDCAxyResolutionParams[2];
  return dcaxy / sigma;
}

float NsigmaDCAzHe(const float pt, const float dcaz) {
  const float sigma = parametrisation::kHeDCAzResolutionParams[0] /
                      std::pow(pt, parametrisation::kHeDCAzResolutionParams[1]) +
                      parametrisation::kHeDCAzResolutionParams[2];
  return dcaz / sigma;
}

float NsigmaDCAzPr(const float pt, const float dcaz) {
  const float sigma = parametrisation::kPrDCAzResolutionParams[0] /
                      std::pow(pt, parametrisation::kPrDCAzResolutionParams[1]) +
                      parametrisation::kPrDCAzResolutionParams[2];
  return dcaz / sigma;
}

// -------------------------------------------- TPC ----------------------------------------------------

double BetheBlochParametrisation(double bg, double kp1, double kp2, double kp3, double kp4, double kp5)
{
  double beta = bg / std::sqrt(1. + bg * bg);
  double aa = std::pow(beta, kp4);
  double bb = std::pow(1. / bg, kp5);
  bb = std::log(kp3 + bb);
  return (kp2 - aa - bb) * kp1 / aa;
}

float BetheBlochHe(const float momentum)
{
  float betagamma = std::abs(momentum) / constant::kMassHe;
  return BetheBlochParametrisation(betagamma, parametrisation::kHeTpcParams[0], parametrisation::kHeTpcParams[1],
                                   parametrisation::kHeTpcParams[2], parametrisation::kHeTpcParams[3],
                                   parametrisation::kHeTpcParams[4]);
}

float NsigmaTpcHe(const float momentum, const float tpcSignal)
{
  return (tpcSignal / BetheBlochHe(std::abs(momentum)) - 1.) / parametrisation::kHeTpcResolution;
}

// -------------------------------------------- ITS ----------------------------------------------------

float averageClusterSize(const uint32_t itsClusterSizes, const bool useTruncatedMean = true)
{
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

float ExpectedClusterSizeCosLambdaHe(const float momentum)
{
    const float betagamma = std::abs(momentum) / constant::kMassHe;
    return parametrisation::kHeITSParams[0] / std::pow(betagamma, parametrisation::kHeITSParams[1]) +
           parametrisation::kHeITSParams[2];
}

float NsigmaITSHe(const float momentum, const float averageClusterSizeCosLambda)
{
    const float betagamma = std::abs(momentum) / constant::kMassHe;
    const float expected = ExpectedClusterSizeCosLambdaHe(momentum);
    const float resolution = parametrisation::kHeITSResolutionParams[0] * 
                  TMath::Erf((betagamma - parametrisation::kHeITSResolutionParams[1]) / parametrisation::kHeITSResolutionParams[2]);
    return (averageClusterSizeCosLambda - expected) / (resolution * expected);
}

float ExpectedClusterSizeCosLambdaPr(const float momentum)
{
    const float betagamma = std::abs(momentum) / constant::kMassPr;
    return parametrisation::kPrITSParams[0] / std::pow(betagamma, parametrisation::kPrITSParams[1]) +
           parametrisation::kPrITSParams[2];
}

float NsigmaITSPr(const float momentum, const float averageClusterSizeCosLambda)
{
    const float betagamma = std::abs(momentum) / constant::kMassPr;
    const float expected = ExpectedClusterSizeCosLambdaPr(momentum);
    const float resolution = parametrisation::kPrITSResolutionParams[0] *
                  TMath::Erf((betagamma - parametrisation::kPrITSResolutionParams[1]) / parametrisation::kPrITSResolutionParams[2]);
    return (averageClusterSizeCosLambda - expected) / (resolution * expected);
}

// -------------------------------------------- TOF ----------------------------------------------------

float NsigmaTOFPr(const float momentum, const float tofMass)
{
  const float expected = constant::kMassPr;
  const float resolution = parametrisation::kPrTofResolutionParams[0] * std::exp(std::abs(momentum) + parametrisation::kPrTofResolutionParams[1]);
  return (tofMass - expected) / (resolution * expected);
}

// -------------------------------------- PID in Tracking ----------------------------------------------

float CorrectPidTrkHe(const float momentum)
{
    return momentum * (1. - parametrisation::kHePidTrkParams[0] -
           parametrisation::kHePidTrkParams[1] * momentum);
}

// ------------------------------------------- Femto ---------------------------------------------------

float Kstar(const double pt1, const double eta1, const double phi1, const double m1, const double pt2, const double eta2, const double phi2, const double m2)
{
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
