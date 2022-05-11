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
  nClusters = 180,
  clusterSpread = 0.3,
  regen_time = 100
)

l |>
  ggplot()+
    geom_point(
      aes(x, y)
    )

# test gamm distributions
rgamma(1000, shape = 0.25, scale = 1) |>
  hist()

# test case 0
a = ecoevomove1::run_model(
  scenario = 1,
  popsize = 500,
  nItems = 1800,
  landsize = 180,
  nClusters = 180,
  clusterSpread = 0.3,
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
  costMove = 0.0,
  tSearch = 5,
  pSearchSlow = 0.5,
  pSearchFast = 0.5,
  pStrategy = 0.5,
  nThreads = 2,
  dispersal = 5.0,
  mProb = 0.01,
  mSize = 0.01
)

b = ecoevomove1::make_network(a, 1)

# ecoevomove1::
plot_network(b, pSearch) +
  scico::scale_fill_scico(
    palette = "romaO",
    direction = 1,
    limits = c(0, 1)
  )

# look at trait in relation to movement and intake
ggplot(a@trait_data)+
  geom_jitter(
    aes(
      moved, assoc, col = pSearch
    )
  )+
  scico::scale_colour_scico(
    palette = "vik",
    direction = 1,
    limits = c(0, 1)
  )+
  theme_test()

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
