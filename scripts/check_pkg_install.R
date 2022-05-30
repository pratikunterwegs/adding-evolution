# check function
Rcpp::compileAttributes()
devtools::build()
# to find RCPPGSL
Sys.setenv("LIB_GSL" = "C:/local323")
devtools::install(build = T, upgrade = "never",
  args = c("--no-multiarch")  
)
devtools::document()

library(ggplot2)

# test landscapes
l = ecoevomove1::get_test_landscape(
  nItems = 450,
  landsize = 30,
  nClusters = 30,
  clusterSpread = 0.5,
  regen_time = 100
)

l |>
  ggplot()+
    geom_point(
      aes(x, y)
    )

# test gamm distributions
rgamma(1000, shape = 0.5, scale = 0.2) |>
  hist()

# test case 0
a = ecoevomove1::run_model(
  scenario = 1,
  popsize = 500,
  nItems = 1800,
  landsize = 30,
  nClusters = 30,
  clusterSpread = 0.5,
  regen_time = 100,
  tmax = 400,
  genmax = 500,
  paramGammaA = 0.5,
  paramGammaB = 0.2,
  paramKappa = 0.1,
  range_perception = 1.0,
  costMove = 0.1,
  nThreads = 2,
  dispersal = 10.0,
  mProb = 0.01,
  mSize = 0.01
)

# check movement data
d = ecoevomove1::get_move_data(a)

ggplot(d[id == 5])+
  geom_path(
    aes(x, y)
  )+
  coord_equal()

d = a@trait_data

ggplot(d)+
  annotate(
    geom = "point",
    x = 0.5, y = 0.1,
    col = "red"
  )+
  geom_jitter(
    aes(gammaA, kappa)
  )+
  labs(
    x = "alpha", y = "kappa"
  )

ggplot(d)+
  geom_jitter(
    aes(gammaA, intake)
  )

b = ecoevomove1::make_network(a, 1)

# ecoevomove1::
plot_network(b, kappa)+
  scale_fill_viridis_c(
    direction = -1
    # limits = c(1e-3, 10)
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
