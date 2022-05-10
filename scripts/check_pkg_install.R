# check function
Rcpp::compileAttributes()
devtools::build()
devtools::load_all()
devtools::install(build = T, upgrade = "never")
devtools::document()

library(ggplot2)

# test case 0
a = run_model(
  scenario = 0,
  popsize = 250,
  nItems = 1800,
  landsize = 200,
  nClusters = 100,
  clusterSpread = 1,
  regen_time = 10,
  tmax = 100,
  genmax = 1,
  paramBallisticGammaA = 5.0,
  paramBallisticGammaB = 1.0,
  paramBallisticKappa = 10.0,
  paramSearchGammaA = 0.05,
  paramSearchGammaB = 1.0,
  paramSearchKappa = 0.1,
  range_perception = 1.0,
  costMove = 0.2,
  tSearch = 5,
  pSearchSlow = 0.8,
  pSearchFast = 0.2,
  pStrategy = 0.8,
  nThreads = 1,
  dispersal = 3.0,
  mProb = 0.01,
  mSize = 0.01
)

str(a)

# plot displacement
ggplot(a@trait_data)+
  geom_segment(
    aes(
      x = x, y = y,
      xend = xn, yend = yn
    )
  )

b = make_network(a, 1)

plot_network(b, pSearch) +
  scale_fill_viridis_c()

# look at trait in relation to movement and intake
ggplot(a@trait_data)+
  geom_jitter(
    aes(
      moved, intake, col = pSearch
    )
  )
