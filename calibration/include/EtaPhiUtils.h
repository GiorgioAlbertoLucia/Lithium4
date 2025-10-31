#pragma once

#include <cmath>
#include <TMath.h>

double ComputePhiStar(const double phi, const double p, const double radius) {
    const double B = 0.5 ; // magnetic field in Tesla
    double phiStar = phi + std::asin(0.3 * radius * B / (2 * p));
    phiStar = std::fmod(phiStar,  2. * TMath::Pi());
    return phiStar;
}
