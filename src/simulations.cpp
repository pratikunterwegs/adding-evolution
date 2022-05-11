#include <vector>
#include <random>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <Rcpp.h>
#include "simulations.h"

using namespace Rcpp;

/// simulation for evolved probability of switching to search
Rcpp::List simulation::do_simulation_evo() {
    unsigned seed = static_cast<unsigned> (std::chrono::system_clock::now().time_since_epoch().count());
    rng.seed(seed);
    
    // prepare landscape and pop
    food.initResources();
    food.countAvailable();
    Rcpp::Rcout << "landscape with " << food.nClusters << " clusters\n";

    pop.setTrait();
    Rcpp::Rcout << "pop with " << pop.nAgents << " agents for " << genmax << " gens " << tmax << " timesteps\n";

    // prepare scenario
    // return scenario as string
    Rcpp::Rcout << "this is scenario " << scenario << ": evolutionary\n";

    // agent random position in first gen
    pop.initPos(food);

    Rcpp::Rcout << "initialised population positions\n";
    Rcpp::DataFrame edgeList;
    Rcpp::DataFrame pop_trait_data;

    Rcpp::Rcout << "created single edge list object\n";

    // all ecological dynamics
    food.countAvailable();
    // reset counter and positions
    pop.counter = std::vector<int> (pop.nAgents, 0);

    // go over gens
    for(int gen = 0; gen < genmax; gen ++) {
        food.countAvailable();
        pop.counter = std::vector<int> (pop.nAgents, 0);

        // timesteps start here
        for (size_t t = 0; t < static_cast<size_t>(tmax); t++)
        {
            // resources regrow
            food.regenerate();
            pop.updateRtree();
            // movement section
            pop.move_random(food);

            md.updateMoveData(pop, t);

            // foraging -- split into parallelised picking
            // and non-parallel exploitation
            pop.pickForageItem(food, nThreads);
            pop.doForage(food);

            // count associations
            pop.countAssoc(nThreads);
            // timestep ends here
        }
        
        // log data in the last generation
        if (gen == (genmax - 1)) {
            pop_trait_data = pop.returnPopData();
            edgeList = pop.pbsn.getNtwkDf();
        }

        // reproduce
        pop.Reproduce(food, dispersal, mProb, mSize);
        // generation ends here

    }
    Rcpp::Rcout << "gen: " << (genmax - 1) << " --- logged edgelist\n";
    Rcpp::Rcout << "data prepared\n";

    return Rcpp::List::create(
        Named("gen_data") = pop_trait_data,
        Named("edge_list") = edgeList,
        Named("move_data") = md.getMoveData()
    );
}

/// simulation for random movement
Rcpp::List simulation::do_simulation_eco() {
    unsigned seed = static_cast<unsigned> (std::chrono::system_clock::now().time_since_epoch().count());
    rng.seed(seed);
    
    // prepare landscape and pop
    food.initResources();
    food.countAvailable();
    Rcpp::Rcout << "landscape with " << food.nClusters << " clusters\n";

    pop.setTrait();
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

        md.updateMoveData(pop, t);

        // foraging -- split into parallelised picking
        // and non-parallel exploitation
        pop.pickForageItem(food, nThreads);
        pop.doForage(food);

        // count associations
        pop.countAssoc(nThreads);
        // timestep ends here
    }
    
    edgeList = pop.pbsn.getNtwkDf();

    Rcpp::Rcout << "gen: " << 1 << " --- logged edgelist\n";
    Rcpp::Rcout << "data prepared\n";

    return Rcpp::List::create(
        Named("gen_data") = pop.returnPopData(),
        Named("edge_list") = edgeList,
        Named("move_data") = md.getMoveData()
    );
}

//' Run a simulation of fast and slow movement.
//' @description Run a simulation in which agents move on a landscape with
//' discrete food items and switch to area restricted search with some
//' probability, after encountering a food item.
//'
//' @param scenario The scenario: 0 for ecological scenario only, or 
//' 1 for evolved probability of switching to area restricted search.
//' @param popsize The population size.
//' @param nItems How many food items on the landscape.
//' @param landsize The size of the landscape as a numeric (double).
//' @param nClusters Number of clusters around which food is generated.
//' @param clusterSpread How dispersed food is around the cluster centre.
//' @param regen_time The item regeneration time.
//' @param tmax The number of timesteps per generation.
//' @param genmax The maximum number of generations per simulation.
//' @param paramBallisticGammaA Alpha parameter for a step length distribution
//' used to draw steps for agent ballistic movement.
//' @param paramBallisticGammaB Beta parameter for a step length distribution
//' used to draw steps for agent ballistic movement.
//' @param paramBallisticKappa Concentration parameter of a von Mises distribution from
//' which turning angles are drawn in radians, for ballistic movement.
//' Should be smaller than `searchAngle`.
//' @param paramSearchGammaA Alpha parameter for a step length distribution
//' used to draw steps for agent searching movement.
//' @param paramSearchGammaB Beta parameter for a step length distribution
//' used to draw steps for agent searching movement.
//' @param paramSearchKappa Concentration parameter of a von Mises distribution from
//' which turning angles are drawn in radians, for searching movement.
//' Should be greater than `ballisticAngle`.
//' @param range_perception The range at which agents detect items.
//' @param costMove The energetic cost per distance moved.
//' @param tSearch The duration of area restricted search; this is fixed, while
//' the probability of switching to this mode varies.
//' @param pSearchSlow The probability of switching to search mode for so-called slow
//' agents.
//' @param pSearchFast The probability of switching to search mode for so-called fast
//' agents.
//' @param pStrategy The initial proportion of fast individuals in the population.
//' @param nThreads How many threads to parallelise over. Set to 1 to run on
//' the HPC Peregrine cluster.
//' @param dispersal A float value; the standard deviation of a normal
//' distribution centred on zero, which determines how far away from its parent
//' each individual is initialised. The standard value is 5 percent of the
//' landscape size (\code{landsize}), and represents local dispersal.
//' Setting this to 10 percent is already almost equivalent to global dispersal.
//' @param mProb The probability of mutation. The suggested value is 0.01.
//' While high, this may be more appropriate for a small population; change this
//' value and \code{popsize} to test the simulation's sensitivity to these values.
//' @param mSize Controls the mutational step size, and represents the scale
//' parameter of a Cauchy distribution. 
//' @return An S4 class, `simulation_output`, with simulation outcomes.
// [[Rcpp::export]]
S4 run_model(const int scenario,
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
                const float paramBallisticKappa,
                const float paramSearchGammaA,
                const float paramSearchGammaB,
                const float paramSearchKappa,
                const float range_perception,
                const float costMove,
                const int tSearch,
                const float pSearchSlow,
                const float pSearchFast,
                const float pStrategy,
                const int nThreads,
                const float dispersal,
                const float mProb,
                const float mSize) {

    // make simulation class with input parameters                            
    simulation this_sim(scenario, popsize, nItems,
        landsize, nClusters, clusterSpread, regen_time, tmax,
        genmax, 
        paramBallisticGammaA, paramBallisticGammaB,
        paramBallisticKappa,
        paramSearchGammaA, paramSearchGammaB,
        paramSearchKappa,
        range_perception, costMove, tSearch,
        pSearchSlow, pSearchFast, pStrategy,
        nThreads, dispersal,
        mProb, mSize
    );

    // return scenario as string
    std::string scenario_str;
    Rcpp::List simOutput;
    if(scenario == 0) {
        scenario_str = std::string("ecological");
        // do the simulation using the simulation class function                        
        simOutput = this_sim.do_simulation_eco();
    } else if(scenario == 1) {
        scenario_str = std::string("evolutionary movement");
        // do the simulation using the simulation class function                        
        simOutput = this_sim.do_simulation_evo();
    }

    // parameter list
    Rcpp::List param_list = Rcpp::List::create(
            Named("scenario") = scenario_str,
            Named("generations") = genmax,
            Named("timesteps") = tmax,
            Named("pop_size") = popsize,
            Named("pop_density") = static_cast<float>(popsize) / landsize,
            Named("item_density") = static_cast<float>(nItems) / landsize,
            Named("pSearchSlow") = pSearchSlow,
            Named("pSearchFast") = pSearchFast,
            Named("pStrategy") = pStrategy,
            Named("dispersal") = dispersal
        );

    // create S4 class pathomove output and fill slots
    S4 x("simulation_output");
    x.slot("parameters") = Rcpp::wrap(param_list);
    x.slot("trait_data") = Rcpp::wrap(simOutput["gen_data"]);
    x.slot("edge_list") = Rcpp::wrap(simOutput["edge_list"]);
    x.slot("move_data") = Rcpp::wrap(simOutput["move_data"]);

    return(x);
}
