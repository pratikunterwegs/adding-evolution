#define _USE_MATH_DEFINES
/// code to make agents
#include <vector>
#include <cassert>
#include <algorithm>
#include <iostream>
#include <random>

#include <boost/foreach.hpp>
#include <RcppGSL.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <Rcpp.h>
#include <RcppParallel.h>

#include "network.h"
#include "landscape.h"
#include "agents.h"

/// random number generator
std::mt19937 rng;
gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);

// to shuffle pop id
void Population::shufflePop() {
    if (order[0] == order[nAgents - 1])
    {
        for (size_t i = 0; i < static_cast<size_t>(nAgents); i++)
        {
            order[i] = i;
        }
        std::random_shuffle ( order.begin(), order.end() );
    }
    else {
        std::random_shuffle ( order.begin(), order.end() );
    }
    
}

// to update agent Rtree
void Population::updateRtree () {
    // initialise rtree
    bgi::rtree< value, bgi::quadratic<16> > tmpRtree;
    for (int i = 0; i < nAgents; ++i)
    {
        point p = point(coordX[i], coordY[i]);
        tmpRtree.insert(std::make_pair(p, i));
    }
    std::swap(agentRtree, tmpRtree);
    tmpRtree.clear();
}

// uniform distribution for agent position
std::uniform_real_distribution<float> agent_ran_pos(0.0f, 1.f);

// function for initial positions
void Population::initPos(Resources food) {
    for (size_t i = 0; i < static_cast<size_t>(nAgents); i++) {
        coordX[i] = agent_ran_pos(rng) * food.dSize;
        initX[i] = coordX[i];
        coordY[i] = agent_ran_pos(rng) * food.dSize;
        initY[i] = coordY[i];
    }
    updateRtree();
}

// set agent probability of switching to area restricted search
void Population::setTrait() {

    // create a cauchy distribution, mSize is the scale
    std::bernoulli_distribution fastOrSlow(pStrategy);

    for(int i = 0; i < nAgents; i++) {
        if(fastOrSlow(rng)) {
            pArs[i] = pArsFast;
        } else {
            pArs[i] = pArsSlow;
        }
    }
}

float get_distance(float x1, float x2, float y1, float y2) {
    return std::sqrt(std::pow((x1 - x2), 2) + std::pow((y1 - y2), 2));
}

// general function for agents within distance
int Population::countAgents (
    const float xloc, const float yloc) {
    
    int ag = 0;
    std::vector<value> near_agents;
    // query for a simple box
    agentRtree.query(bgi::satisfies([&](value const& v) {
        return bg::distance(v.first, point(xloc, yloc)) < range_perception;}),
        std::back_inserter(near_agents));

    BOOST_FOREACH(value const& v, near_agents) {
        
        ag++;
    }
    near_agents.clear();
    // first element is number of near entities
    // second is the identity of entities
    return ag;
}

// function for near agent ids
std::vector<int> Population::getNeighbourId (
    const float xloc, const float yloc) {
    
    std::vector<int> agent_id;
    std::vector<value> near_agents;
    // query for a simple box
    // neighbours for associations are counted over the MOVEMENT RANGE
    agentRtree.query(bgi::satisfies([&](value const& v) {
        return bg::distance(v.first, point(xloc, yloc)) < range_perception;}),
        std::back_inserter(near_agents));

    BOOST_FOREACH(value const& v, near_agents) {
        agent_id.push_back(v.second);
    }
    near_agents.clear();
    // first element is number of near entities
    // second is the identity of entities
    return agent_id;
}

// general function for items within distance
int Population::countFood (
    const Resources &food,
    const float xloc, const float yloc) {

    int nFood = 0;
    std::vector<value> near_food;

    // check any available
    if (food.nAvailable > 0) {
        // query for a simple box
        food.rtree.query(bgi::satisfies([&](value const& v) {
            return bg::distance(v.first, point(xloc, yloc)) < range_perception;}),
            std::back_inserter(near_food));

        BOOST_FOREACH(value const& v, near_food) {
            // count only which are available!
            if (food.available[v.second]) {
                nFood++;
            }
        }
        near_food.clear();
    }

    return nFood;
}

// function for the nearest available food item
std::vector<int> Population::getFoodId (
    const Resources &food,
    const float xloc, const float yloc) {
        
    std::vector<int> food_id;
    std::vector<value> near_food;
    // check any available
    if (food.nAvailable > 0) {
        // query for a simple box
        // food is accessed over the MOVEMENT RANGE
        food.rtree.query(bgi::satisfies([&](value const& v) {
            return bg::distance(v.first, point(xloc, yloc)) < range_perception;}), 
            std::back_inserter(near_food));

        BOOST_FOREACH(value const& v, near_food) {
            // count only which are available!
            if (food.available[v.second]) {
                food_id.push_back(v.second);
            }
        }
        near_food.clear();
    }

    // first element is number of near entities
    // second is the identity of entities
    return food_id;
}

/// rng for suitability
std::normal_distribution<float> noise(0.f, 0.01f);
std::cauchy_distribution<float> noise_cauchy(0.f, 0.001f);

/// function to wrap location given a maximum size
float wrapLoc(float l, float maxl) {
    if(l > maxl) {
        return (l - maxl);
    } else if(l < 0.f) {
        return(maxl + l);
    } else {
        return l;
    }
}

// function to move after drawing step lengths from a distribution
// and also turning angles
void Population::move_random(const Resources &food) {
    float twopi = 2.f * M_PI;

    // set up distributions
    std::gamma_distribution<float> distanceBallistic (paramBallisticGammaA, paramBallisticGammaB);
    std::normal_distribution<float> angleBallistic (0.f, paramBallisticNormalSD);

    std::gamma_distribution<float> distanceArs (paramArsGammaA, paramArsGammaB);
    std::normal_distribution<float> angleArs (0.f, paramArsNormalSD);
    
    for (int i = 0; i < nAgents; ++i) {
        // check if locked into search mode
        if (counter[i] > 0) {

            // agent is searching and moves with brownian motion
            float distance = distanceArs(rng);
            float angle = angleArs(rng);

            float t1_ = static_cast<float>(cos(angle));
            float t2_ = static_cast<float>(sin(angle));

            coordX[i] = coordX[i] + (distance * t1_);
            coordY[i] = coordY[i] + (distance * t2_);

            coordX[i] = wrapLoc(coordX[i]);
            coordY[i] = wrapLoc(coordY[i]);

            // movement and cost of movement
            moved[i] += distance;
            energy[i] -= (distance * costMove);

            // agent counter for search reduces
            counter[i] --;
        }
        else {
            // agent is searching and moves with brownian motion
            float distance = distanceBallistic(rng);
            float angle = angleBallistic(rng);

            float t1_ = static_cast<float>(cos(angle));
            float t2_ = static_cast<float>(sin(angle));

            coordX[i] = coordX[i] + (distance * t1_);
            coordY[i] = coordY[i] + (distance * t2_);

            coordX[i] = wrapLoc(coordX[i]);
            coordY[i] = wrapLoc(coordY[i]);
            
            // movement and cost of movement
            moved[i] += distance;
            energy[i] -= (distance * costMove);   
        }
    }
}

// function to paralellise choice of forage item
void Population::pickForageItem(const Resources &food, const int nThreads){
    shufflePop();
    // nearest food
    std::vector<int> idTargetFood (nAgents, -1);

    if (nThreads > 1)
    {
        // loop over agents --- no shuffling required here
        tbb::task_scheduler_init _tbb(tbb::task_scheduler_init::automatic); // automatic for now
        // try parallel foraging --- agents pick a target item
        tbb::parallel_for(
            tbb::blocked_range<unsigned>(1, order.size()),
                [&](const tbb::blocked_range<unsigned>& r) {
                for (unsigned i = r.begin(); i < r.end(); ++i) {
                    if ((counter[i] > 0) | (food.nAvailable == 0)) { 
                        // nothing -- agent cannot forage or there is no food
                    }
                    else {
                        // find nearest item ids
                        std::vector<int> theseItems = getFoodId(food, coordX[i], coordY[i]);
                        int thisItem = -1;

                        // check near items count
                        if(theseItems.size() > 0) {
                            // take first item by default
                            thisItem = theseItems[0];
                            idTargetFood[i] = thisItem;
                        }
                    }
                }
            }
        );
    } else if (nThreads == 1)
    {
        for (int i = 0; i < nAgents; ++i) {
            if ((counter[i] > 0) | (food.nAvailable == 0)) { 
                // nothing -- agent cannot forage or there is no food
            }
            else {
                // find nearest item ids
                std::vector<int> theseItems = getFoodId(food, coordX[i], coordY[i]);
                int thisItem = -1;

                // check near items count
                if(theseItems.size() > 0) {
                    // take first item by default
                    thisItem = theseItems[0];
                    idTargetFood[i] = thisItem;
                }
            }
        }
    }

    forageItem = idTargetFood;
}

// function to exploitatively forage on picked forage items
void Population::doForage(Resources &food) {
    // all agents have picked a food item if they can forage
    // now forage in a serial loop --- this cannot be parallelised
    // this order is randomised
    for (size_t i = 0; i < static_cast<size_t>(nAgents); i++)
    {
        int id = order[i];
        if ((food.nAvailable == 0)) {
            // do nothing if there is no food available
        } else {
            int thisItem = forageItem[id]; //the item picked by this agent
            // check selected item is available
            if (thisItem != -1)
            {
                intake[id] += 1.0; // increased here --- not as described --- okay for now.
                energy[id] += 1.0;
                // individuals have some probability of shifting to ars
                if(gsl_ran_bernoulli(r, pArs[i]) == 1) {
                    counter[id] = tArs;
                }

                // reset food availability
                food.available[thisItem] = false;
                food.counter[thisItem] = food.regen_time;
                food.nAvailable --;
            }
        }
    }
}

void Population::countAssoc(const int nThreads) {
    for (int i = 0; i < nAgents; ++i) {
        // count nearby agents and update raw associations
        std::vector<int> nearby_agents = getNeighbourId(coordX[i], coordY[i]);
        associations[i] += nearby_agents.size();
        // subtract 1 to exclude self

        // loop over nearby agents and update association matrix
        for (size_t j = 0; j < nearby_agents.size(); j++)
        {
            int target_agent = nearby_agents[j];
            pbsn.adjMat (i, target_agent) += 1;
        }
    }
}

/// minor function to normalise vector
std::vector<float> Population::handleFitness() {
    // sort vec fitness
    std::vector<float> vecFitness = energy;
    std::sort(vecFitness.begin(), vecFitness.end()); // sort to to get min-max
    // scale to max fitness
    float maxFitness = vecFitness[vecFitness.size()-1];
    float minFitness = vecFitness[0];

    // reset to energy
    vecFitness = energy;
    // rescale copied energy vector by min anx max fitness
    for(size_t i = 0; i < static_cast<size_t>(nAgents); i++) {
        vecFitness[i] = ((vecFitness[i]  - minFitness) / (maxFitness - minFitness)) + 0.0001f;
    }
    
    return vecFitness;
}

// fun for replication
void Population::Reproduce(const Resources food,
    const float dispersal, const float mProb, const float mSize) 
{
    // mutation probability and size distribution --- inefficient but oh well
    std::bernoulli_distribution mutation_happens(mProb);
    std::cauchy_distribution<float> mutation_size(0.0, mSize);

    // choose the range over which individuals are dispersed
    std::normal_distribution<float> sprout(0.f, dispersal);
    std::vector<float> vecFitness;
    vecFitness = handleFitness();

    // set up weighted lottery
    std::discrete_distribution<> weightedLottery(vecFitness.begin(), vecFitness.end());

    // get parent trait based on weighted lottery
    std::vector<float> tmp_pArs (nAgents, 0.5f);
    
    // reset associations
    associations = std::vector<int> (nAgents, 0);

    // reset distance moved
    moved = std::vector<float> (nAgents, 0.f);

    // reset adjacency matrix
    pbsn.adjMat = Rcpp::NumericMatrix(nAgents, nAgents);

    // reset positions
    std::vector<float> coord_x_2 (nAgents, 0.f);
    std::vector<float> coord_y_2 (nAgents, 0.f);
    
    for (int a = 0; a < nAgents; a++) {
        size_t parent_id = static_cast<size_t>(weightedLottery(rng));

        tmp_pArs[a] = pArs[parent_id];

        // inherit positions from parent
        coord_x_2[a] = coordX[parent_id] + sprout(rng);
        coord_y_2[a] = coordY[parent_id] + sprout(rng);

        // robustly wrap positions
        coord_x_2[a] = wrapLoc(coord_x_2, food.dSize);
        coord_y_2[a] = wrapLoc(coord_y_2, food.dSize);
    }

    // swap coords --- this initialises individuals near their parent's position
    std::swap(coordX, coord_x_2);
    std::swap(coordY, coord_y_2);
    coord_x_2.clear(); coord_y_2.clear();

    // update initial positions!
    initX = coordX;
    initY = coordY;

    // reset counter
    counter = std::vector<int> (nAgents, 0);
    assert(static_cast<int>(counter.size()) == nAgents && "counter size wrong");

    // mutate trait: trait shifts up or down with an equal prob
    // trait mutation prob is mProb, in a two step process
    for (int a = 0; a < nAgents; a++) {
        if(mutation_happens(rng)) {
            tmp_pArs[a] = tmp_pArs[a] + mutation_size(rng);
            
            if(tmp_pArs[a] < 0.f) tmp_pArs = 0.00001f;
            if(tmp_pArs[a] > 1.f) tmp_pArs = 0.99f;
        }
    }
    
    // swap trait matrices
    std::swap(pArs, tmp_pArs);

    tmp_pArs.clear();
    
    // swap energy
    std::vector<float> tmpEnergy (nAgents, 0.001f);
    std::swap(energy, tmpEnergy);
    tmpEnergy.clear();

    // swap intake
    std::vector<float> tmpIntake (nAgents, 0.001f);
    std::swap(intake, tmpIntake);
    tmpIntake.clear();
}
