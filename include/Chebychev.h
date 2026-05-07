#include "TF1.h"
#include "TH1.h"

#include <vector>


class Chebychev {
public: 
   Chebychev(int n, double xmin, double xmax) : 
      fA(xmin), fB(xmax),
      fT(std::vector<double>(n) )  {}

   double operator() (const double * xx, const double *p) { 
      double x = (xx[0] - fA -fB)/(fB-fA);
      int order = fT.size(); 
      if (order == 1) return p[0]; 
      if (order == 2) return p[0] + x*p[1]; 
      // build the polynomials
      fT[0] = 1;
      fT[1] = x; 
      for (int i = 1; i< order; ++i) { 
         fT[i+1] =  2 *x * fT[i] - fT[i-1]; 
      }
      double sum = p[0]*fT[0]; 
      for (int i = 1; i<= order; ++i) { 
         sum += p[i] * fT[i]; 
      }
      return sum; 
   }

private: 
   double fA; 
   double fB; 
   std::vector<double> fT; // polynomial
   std::vector<double> fC; // coefficients
};
