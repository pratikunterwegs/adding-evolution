#ifndef VONMISES_H
#define VONMISES_H

/// copied from the rotations R package
// this is an expensive function but no better implementations
// appear to exist

#include <vector>
#include <Rcpp.h>

/// helper function for vonmises
int sign(float x)
{
  if (x < 0)
    return -1;

  return 1;
}

/// Sample from the circular von Mises distribution
float vonMisesAngle(float kappa)
{
  Rcpp::RNGScope scope;
  Rcpp::NumericVector u(3);
  float theta = 10.f;

  u = Rcpp::runif(3, 0, 1);
  float a = 1.0f + std::sqrt(1.0f + 4.0f * std::pow(kappa, 2.0f));
  float b = (a - std::sqrt(2.0f * a)) / (2.0f * kappa);
  float r = (1.0f + std::pow(b, 2.0f)) / (2.0f * b);
  float z = 0.0f;
  float f = 0.0f;
  float c = 0.0f;

  while (theta > 4)
  {
    // Step 1
    u = Rcpp::runif(3, 0, 1);
    z = std::cos(M_PI * u[0]);
    f = (1.0 + r * z) / (r + z);
    c = kappa * (r - f);

    // Step 2
    u = Rcpp::runif(3, 0, 1);
    if ((c * (2.0 - c) - u[1]) > 0)
      theta = (sign(u[2] - 0.5)) * std::acos(f);
    else
    {
      if ((std::log(c / u[1]) + 1.0 - c) < 0)
        u = Rcpp::runif(3, 0, 1);
      else {
        u = Rcpp::runif(3, 0, 1);
        theta = (sign(u[2] - 0.5)) * std::acos(f);
      }
    }
  }
  return theta;
}

#endif // VONMISES_H
