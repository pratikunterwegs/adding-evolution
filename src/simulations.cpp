#include <vector>
#include <random>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <Rcpp.h>
#include "simulations.h"

using namespace Rcpp;


/// simulation for random movement
Rcpp::List simulation::do_simulation_eco() {
    unsigned seed = static_cast<unsigned> (std::chrono::system_clock::now().time_since_epoch().count());
    rng.seed(seed);
    
    // prepare landscape and pop
    food.initResources();
    food.countAvailable();
    Rcpp::Rcout << "landscape with " << food.nClusters << " clusters\n";

    pop.setTrait(mSize);
    Rcpp::Rcout << "pop with " << pop.nAgents << " agents for " << genmax << " gens " << tmax << " timesteps\n";

    // prepare scenario
    // return scenario as string
    Rcpp::Rcout << "this is scenario " << scenario << ": ecological\n";

    // agent random position in first gen
    pop.initPos(food);

    Rcpp::Rcout << "initialised population positions\n";
    Rcpp::DataFrame edgeList;

    Rcpp::Rcout << "created single edge list object\n";

    // all ecological dynamics
    food.countAvailable();
    // reset counter and positions
    pop.counter = std::vector<int> (pop.nAgents, 0);

    // timesteps start here
    for (size_t t = 0; t < static_cast<size_t>(tmax); t++)
    {
        // resources regrow
        food.regenerate();
        pop.updateRtree();
        // movement section
        pop.move_random(food);

        // mdPre.updateMoveData(pop, t);

        // foraging -- split into parallelised picking
        // and non-parallel exploitation
        pop.pickForageItem(food, nThreads);
        pop.doForage(food);

        // count associations
        pop.countAssoc(nThreads);
        // timestep ends here
    }
    pop.energy = pop.intake;
    
    edgeList = pop.pbsn.getNtwkDf();

    Rcpp::Rcout << "gen: " << 1 << " --- logged edgelist\n";
    Rcpp::Rcout << "data prepared\n";

    return Rcpp::List::create(
        Named("gen_data") = pop.returnPopData(),
        Named("edge_list") = edgeList
        // Named("move_data") = mdPre.getMoveData()
    );
}

