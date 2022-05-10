#ifndef PARAMETERS_H
#define PARAMETERS_H

// Enable C++14 via this plugin to suppress 'long long' errors
// [[Rcpp::plugins("cpp14")]]
// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::depends(RcppGSL)]]
#include <random>
#include <RcppGSL.h>
#include <chrono>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

extern std::mt19937 rng;
extern gsl_rng *r;

// landscape
const double foodEnergy = 1.0;

#endif // PARAMETERS_H
