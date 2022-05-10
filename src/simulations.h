#ifndef SIMULATIONS_H
#define SIMULATIONS_H

#include "data_types.h"
#include "landscape.h"
#include "agents.h"

class simulation {
public:
    simulation(const int scenario,
                const int popsize,
                const int nItems, 
                const float landsize,
                const int nClusters,
                const float clusterSpread,
                const int regen_time,
                const int tmax,
                const int genmax,
                const float paramBallisticGammaA,
                const float paramBallisticGammaB,
                const float paramBallisticNormalSD,
                const float paramSearchGammaA,
                const float paramSearchGammaB,
                const float paramSearchNormalSD,
                const float range_perception,
                const float costMove,
                const int tSearch,
                const float pSearchSlow,
                const float pSearchFast,
                const float pStrategy,
                const int nThreads,
                const float dispersal,
                const float mProb,
                const float mSize):
        // population, food, and data structures
        pop (popsize, paramBallisticGammaA, paramBallisticGammaB,
            paramBallisticNormalSD, paramSearchGammaA,
            paramSearchGammaB, paramSearchNormalSD,
            range_perception,
            costMove, tSearch,
            pSearchSlow, pSearchFast, pStrategy,
            scenario
        ),
        food(nItems, landsize, nClusters, clusterSpread, regen_time),
        
        // eco-evolutionary parameters
        scenario(scenario),
        tmax(tmax),
        genmax(genmax),

        // // agent perception and behaviour, food growth
        // range_perception(range_perception),
        // tSearch(tSearch),
        // regen_time(regen_time),

        // parallelisation
        nThreads (nThreads),

        // natal dispersal
        dispersal(dispersal),
        // mutation probability and step size
        mProb(mProb),
        mSize(mSize)

        // movement data
        // mdPre(tmax, popsize),
        // mdPost(tmax, popsize)
    {}
    ~simulation() {}

    Population pop;
    Resources food;
    const int scenario, tmax, genmax;
    // const float range_perception;
    // const int tSearch;
    // const int regen_time;
    int nThreads;
    const float dispersal;
    
    const float mProb, mSize;

    // moveData mdPre, mdPost;

    // funs
    Rcpp::List do_simulation_eco();

};

#endif // SIMULATIONS_H
