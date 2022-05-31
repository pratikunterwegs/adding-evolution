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
    )+
    coord_equal()

# test gamm distributions
mu = 1
rexp(n = 1000, rate = 1/(mu)) |>
  hist()

rotations::rvmises(n = 1000, kappa = 1) |>
  # as.vector() |>
  hist()

# test case 0
a = ecoevomove1::run_model(
  scenario = 1,
  popsize = 250,
  nItems = 450,
  landsize = 30,
  nClusters = 30,
  clusterSpread = 0.5,
  regen_time = 10,
  tmax = 400,
  genmax = 500,
  paramMu = 1,
  paramKappa = 1,
  range_perception = 1,
  costMove = 0.05,
  nThreads = 2,
  dispersal = 10.0,
  mProb = 0.05,
  mSize = 0.01
)

d = ecoevomove1::get_trait_data(a)

ggplot(d[gen %in% c(min(gen), max(gen))])+
  geom_histogram(
    aes(moved, col = as.factor(gen)),
    alpha = 0.2
  )

ggplot(d[gen %in% c(min(gen), max(gen))])+
  geom_histogram(
    aes(mu, col = as.factor(gen)),
    alpha = 0.2
  )

ggplot(d[gen %in% c(min(gen), max(gen))])+
  geom_histogram(
    aes(kappa, col = as.factor(gen)),
    alpha = 0.2
  )

ggplot(d)+
  geom_bin_2d(
    aes(
      gen, mu
    ),
    binwidth = c(1, 0.01)
  )+
  scale_y_log10()+
  scale_fill_viridis_c(
    option = "A",
    direction = -1
  )+
  coord_cartesian(
    # ylim = c(1e-2, NA)
  )

ggplot(d[gen %in% c(max(gen))])+
  annotate(
    geom = "point",
    x = 0.5, y = 0.1,
    col = "red"
  )+
  geom_jitter(
    aes(mu, kappa, col = gen)
  )+
  labs(
    x = "mu", y = "kappa"
  )

ggplot(d[gen %in% c(max(gen))])+
  geom_jitter(
    aes(mu, energy)
  )

ggplot(d[gen %in% c(max(gen))])+
  geom_jitter(
    aes(kappa, intake)
  )

b = ecoevomove1::make_network(, 12)

# ecoevomove1::
plot_network(b, gammaA)+
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
