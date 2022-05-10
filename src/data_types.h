#ifndef DATA_TYPES_H
#define DATA_TYPES_H
/// data types and related functions
#include <vector>
#include <random>
#include <iostream>
#include <fstream>
#include <algorithm>
#include "agents.h"
#include "network.h"

struct moveData {
public:

    moveData(const int tmax, const int popsize) :
        tmax(tmax),
        popsize(popsize),
        timesteps(tmax, std::vector<int>(popsize, 0)),
        x(tmax, std::vector<float>(popsize, 0)),
        y(tmax, std::vector<float>(popsize, 0))
    {}
    ~moveData() {}

    const int tmax;
    const int popsize;
    std::vector<std::vector<int> > timesteps;
    std::vector<std::vector<float> > x;
    std::vector<std::vector<float> > y;

    void updateMoveData (Population &pop, const int t_);
    Rcpp::List getMoveData();
};

#endif  //
