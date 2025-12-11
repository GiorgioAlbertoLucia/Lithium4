#include <TMath.h>

double RooGausExpEvaluate(double * xs, double * pars) {

    double x = xs[0];
    double norm = pars[0], mu = pars[1], sigma = pars[2], tau = pars[3];
    double u = (x - mu) / sigma;
    double abstau = tau;
    if (tau < 0) {
        u = -u;
        abstau = -tau;
    }
    
    if (u <= abstau)
        return TMath::Exp(-u * u * 0.5);
    else
        return TMath::Exp(-abstau * (u - 0.5 * abstau));
}