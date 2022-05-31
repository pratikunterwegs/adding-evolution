
#include "data_types.h"

// function to return gen data as an rcpp list
Rcpp::DataFrame Population::returnPopData() {
    Rcpp::DataFrame this_data = Rcpp::DataFrame::create(
        Named("intake") = intake,
        Named("energy") = energy,
        Named("mu") = paramMu,
        Named("kappa") = paramKappa,
        Named("x") = initX,
        Named("y") = initY,
        Named("xn") = coordX,
        Named("yn") = coordY,
        Named("assoc") = associations,
        Named("moved") = moved
    );
    return this_data;
}

void moveData::updateMoveData(Population &pop, const int t_) {
    assert(t_ <= tmax && "too many timesteps logged");

    timesteps[t_] = std::vector<int> (popsize, t_);
    x[t_] = pop.coordX;
    y[t_] = pop.coordY;

}

Rcpp::List moveData::getMoveData() {
    Rcpp::List mDataList (tmax);
    std::vector<int> id(popsize, 0);
    std::iota(std::begin(id), std::end(id), 0);

    for (int i = 0; i < tmax; i++)
    {
        mDataList[i] = DataFrame::create(
            Named("time") = timesteps[i],
            Named("x") = x[i],
            Named("y") = y[i],
            Named("id") = id
        );
    }

    return mDataList;
}