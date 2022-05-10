# check function
Rcpp::compileAttributes()
devtools::build()
devtools::load_all()
devtools::install(build = T, upgrade = "never")
devtools::document()

library(ggplot2)

# test landscapes
l = ecoevomove1::get_test_landscape(
  nItems = 1800,
  landsize = 60,
  nClusters = 60,
  clusterSpread = 1,
  regen_time = 100
)

l |>
  ggplot()+
    geom_point(
      aes(x, y)
    )

# test gamm distributions
ballistic = rgamma(1000, shape = 0.5, scale = 1)
ballistic |>
  hist()

# test case 0
a = ecoevomove1::run_model(
  scenario = 1,
  popsize = 500,
  nItems = 1800,
  landsize = 60,
  nClusters = 60,
  clusterSpread = 1.0,
  regen_time = 100,
  tmax = 100,
  genmax = 1000,
  paramBallisticGammaA = 0.25,
  paramBallisticGammaB = 1.0,
  paramBallisticKappa = 10.0,
  paramSearchGammaA = 0.25,
  paramSearchGammaB = 1.0,
  paramSearchKappa = 0.1,
  range_perception = 1.0,
  costMove = 0.075,
  tSearch = 5,
  pSearchSlow = 0.5,
  pSearchFast = 0.5,
  pStrategy = 0.5,
  nThreads = 2,
  dispersal = 5.0,
  mProb = 0.01,
  mSize = 0.01
)

str(a)

b = ecoevomove1::make_network(a, 5)

# ecoevomove1::
plot_network(b, pSearch) +
  scale_fill_distiller(
    palette = "RdBu",
    limits = c(0, 1)
  )

# look at trait in relation to movement and intake
ggplot(a@trait_data)+
  geom_jitter(
    aes(
      moved, assoc, col = pSearch
    )
  )+
  scale_colour_viridis_c(
    option = "H",
    limits = c(0, 1)
  )

# get movement data and plot
d = get_move_data(a)

ggplot(d)+
  geom_path(
    aes(
      x, y, col = as.factor(id)
    ),
    show.legend = FALSE
  )+
  geom_point(
    aes(
      x, y, col = as.factor(id)
    ),
    show.legend = FALSE
  )+
  coord_fixed(
    xlim = c(0, 200),
    ylim = c(0, 200)
  )
