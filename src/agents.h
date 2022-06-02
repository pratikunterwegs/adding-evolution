#ifndef AGENTS_H
#define AGENTS_H

#define _USE_MATH_DEFINES
/// code to make agents
#include <vector>
#include <cassert>
#include <algorithm>
#include <iostream>
#include <boost/foreach.hpp>
#include "landscape.h"
#include "network.h"


// Agent class
struct Population {
public:
    Population(const int popsize, 
        const float paramMu,
        const float paramKappa,
        const float pMove,
        const float range_perception,
        const float costMove,
        const int scenario
    ) :
        // agents, positions, energy and traits
        nAgents (popsize),
        coordX (popsize, 0.0f),
        coordY (popsize, 0.0f),
        initX (popsize, 0.0f),
        initY (popsize, 0.0f),
        intake (popsize, 0.001f),
        energy (popsize, 0.001f),

        // parameters for the initial random walk
        paramMu(popsize, paramMu),
        paramKappa(popsize, paramKappa),
        pMove(popsize, pMove),

        range_perception(range_perception),
        costMove(costMove),
        
        // counters for handling and social metrics
        counter (popsize, 0),
        associations(popsize, 0),

        // vectors for agent order
        order(popsize, 1),
        forageItem(popsize, -1),
        
        // distance moved
        moved(popsize, 0.f),

        // a network object
        pbsn(popsize)
    {}
    ~Population() {}

    // agent count, coords, and energy
    const int nAgents;
    std::vector<float> coordX;
    std::vector<float> coordY;
    std::vector<float> initX;
    std::vector<float> initY;
    std::vector<float> intake;
    std::vector<float> energy;
    
    std::vector<float> paramMu;
    std::vector<float> paramKappa;
    std::vector<float> pMove;

    const float range_perception;
    const float costMove;

    // counter and metrics
    std::vector<int> counter;
    std::vector<int> associations; // number of total interactions

    // shuffle vector and transmission
    std::vector<int> order;
    std::vector<int> forageItem;

    // movement distances
    std::vector<float> moved;

    // network object
    Network pbsn;

    // position rtree
    bgi::rtree< value, bgi::quadratic<16> > agentRtree;

    /// functions for the population ///
    // population order, trait and position randomiser
    void shufflePop();
    void initPos(Resources food);

    int countFood (const Resources &food, const float xloc, const float yloc);
    
    std::vector<int> getFoodId (
        const Resources &food,
        const float xloc, const float yloc
    );
    
    std::pair<int, int > countAgents (
        const float xloc, const float yloc);
    
    // functions to move and forage on a landscape
    void move_random(const Resources &food);
    void pickForageItem(const Resources &food, const int nThreads);
    void doForage(Resources &food);
    
    // funs to handle fitness and reproduce
    std::vector<float> handleFitness();
    void Reproduce(const Resources food, 
        const float dispersal,
        const float mProb,
        const float mSize
    );
    
    // counting proximity based interactions
    // make rtree and get nearest agents and food
    void updateRtree();
    std::vector<int> getNeighbourId(const float xloc, const float yloc);
    void countAssoc(const int nThreads);

    // return population data
    Rcpp::DataFrame returnPopData();
};

// a dinky function for distance and passed to catch test
float get_distance(float x1, float x2, float y1, float y2);

#endif // AGENTS_H
