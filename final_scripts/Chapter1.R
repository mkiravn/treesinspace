# Loading libraries
devtools::load_all("~/treesinspace")
# library(tidyverse)
# library(ape)
# library(sf)
# library(MetBrewer)
# library(ggrepel)
# library(ggpubr)
# library(ggtree)
# library(moments)
# library(phangorn)
# library(apTreeshape)
# library(VGAM)
# library(scales)
# library(beeswarm)
# library(ggbeeswarm)
devtools::load_all("~/slendr")
setup_env()
check_env()
set.seed(206)

######### Here we define parameters

# defining parameters
Ns <- 100 # number of individuals in population
mat_dists <- c(1, 5, 10, 100) # mating distance
comp_dists <- c(0, 0.2, 20, 40) # competition distance
disp_dists <- c(1:3) # dispersal distance (sigma)
disp_funs <-
  c("brownian", "normal", "cauchy", "exponential", "uniform") # dispersal kernel
reps <-
  c(1:3) # number of simulation runs for each parameter combination
ngens <- 50 # number of generations

# generating a parameter table
pars <- expand.grid(
  N = Ns,
  mat_dist = mat_dists,
  comp_dist = comp_dists,
  disp_fun = disp_funs,
  sigma = disp_dists,
  rep = reps,
  ngen = ngens
) %>%
  rownames_to_column("sim")
print(paste("Running" , dim(pars)[1], "simulations."))
write_delim(
  x = pars,
  file = paste0(as.character(Sys.Date()), "-Chapter1.tsv"),
  delim = "\t"
)



######### Here we run the simulations


# these will hold our data
trees <- c()
trees_us <- c()


# defining a world
map <- world(
  xrange = c(-25, 25),
  # min-max longitude
  yrange = c(-25, 25),
  # min-max latitude
  landscape = "blank"
)

i <- 1
for (row in c(1:dim(pars)[1])) {
  # define a population
  pop <- population(
    name="POP",
    time = 1,
    N = pars[row, "N"],
    center = c(0, 0),
    radius = 50,
    map = map,
    mating = pars[row, "mat_dist"],
    competition = pars[row, "comp_dist"],
    dispersal_fun = as.character(pars[row, "disp_fun"]),
    dispersal = pars[row, "sigma"]
  )

  # compile and run the model
  model <- compile_model(
    populations = pop,
    generation_time = 1,
    simulation_length = ngens,
    resolution = 1,
    # resolution in "distance units per pixel"
    path = "./logs/Chapter1",
    overwrite = TRUE,
    force=TRUE
  )

  # sampling strategy
  samples <-
    schedule_sampling(model, times = ngens, list(pop, pars[row, "N"] / 2))

  # simulate
  # simulate
  try_again(
    10,
    slim(
      model,
      samples = samples,
      # simulate only a single locus
      sequence_length = 1,
      recombination_rate = 0,
      method = "batch",
      # change to "gui" to execute the model in SLiMgui
      random_seed = sample(c(1:100)),
      verbose = FALSE,
      coalescent_only = FALSE,
      burnin = 0,
    )
  )


  # extract simplified trees
  ts <- ts_load(model) %>%
    ts_recapitate(
      Ne = pars[row, "N"],
      recombination_rate = 0,
      random_seed = 206
    ) %>%
    ts_simplify()

  # extract unsimplified tree
  tsu <- ts_load(model, simplify = FALSE)

  # collect trees
  trees <- c(trees, ts)
  trees_us <- c(trees_us, tsu)

  print(paste("Simulation number:", i))
  i <- i + 1
}

######### Here we convert our data to something more usable

# getting out ape trees
trees_phylo <-
  lapply(trees, ts_phylo, i = 1, quiet = TRUE) # convert
class(trees_phylo) <-  "multiPhylo" # convert list

# merge all the spatial tree data into a dataframe
# first for simplified trees
tree_data <- data.frame()
for (row in c(1:dim(pars)[1])) {
  tree <- trees_phylo[[row]]
  tree_data_i <-
    tree %>% ts_data() %>% rowwise() %>% mutate(
      N = pars[row, "N"],
      fun = pars[row, "disp_fun"],
      rep = as.factor(pars[row, "rep"]),
      mat_dist =
        pars[row, "mat_dist"],
      comp_dist =
        pars[row, "comp_dist"],
      sigma =
        pars[row, "sigma"],
      sim = as.factor(row),
      x = unlist(location)[1],
      y = unlist(location)[2],
      simplified =
        "simplified"
    ) %>% ungroup()

  tree_data <- rbind(tree_data, tree_data_i)
  print(paste("Converted tree no", row))
}
# now for unsimplified trees
for (row in c(1:dim(pars)[1])) {
  tree <- trees_us[[row]]
  tree_data_i <-
    tree %>% ts_data() %>% rowwise() %>% mutate(
      N = pars[row, "N"],
      fun = pars[row, "disp_fun"],
      rep = as.factor(pars[row, "rep"]),
      mat_dist =
        pars[row, "mat_dist"],
      comp_dist =
        pars[row, "comp_dist"],
      sigma =
        pars[row, "sigma"],
      sim = as.factor(row),
      x = unlist(location)[1],
      y = unlist(location)[2],
      simplified =
        "unsimplified",
      phylo_id =
        NA
    ) %>% ungroup()

  tree_data <- rbind(tree_data, tree_data_i)
  print(paste("Converted unsimplified tree no", row))
}

# save the tree data
tree_data <- tree_data %>%
  group_by(sim) %>%
  # adjust so that time starts at 0 (for recapitated trees)
  mutate(timeoriginal = time,
         time = time - min(time))

# the connections (spatial lines between nodes)
conns_simp <- lapply(trees, ts_connect)
conns <- lapply(trees_us, ts_connect)

all_connections_l <- data.frame() # this will hold the connections

for (row in c(1:dim(pars)[1])) {
  print(paste("Connecting tree number:", row))
  all_connections_i <- conns[[row]][[1]] %>%
    as.data.frame() %>%
    mutate(simplified = "unsimplified") %>%
    cbind(pars[row,])
  all_connections_i <- conns_simp[[row]][[1]] %>%
    as.data.frame() %>%
    mutate(simplified = "simplified") %>%
    cbind(pars[row,]) %>% rbind(all_connections_i)
  all_connections_l <- rbind(all_connections_l, all_connections_i)
}


######### Plotting
global_labeller <- labeller(
  .default = label_value,
  `mating distance` = label_both,
  `competition distance` = label_both,
  sigma = label_both
)

all_connections_l <- all_connections_l  %>%
  rename(
    "competition distance" = "comp_dist",
    "mating distance" = "mat_dist",
    "distance" = "dist",
    "sigma" = "sigma",
    "dispersal function" = "disp_fun",
    "x displacement" = "x_dist",
    "y displacement" = "y_dist"
  )


## Plotting curves of parent-offspring distance
# 1. Mating distance and competition
all_connections_l %>%
  filter(simplified == "unsimplified",
         sigma == 1) %>%
  ggplot(aes(x = distance,
             col = as.factor(`competition distance`),)) +
  geom_line(stat = "density", size = 0.8) +
  lims(x = c(0, 10)) +
  theme_minimal() +
  geom_vline(
    aes(xintercept = `mating distance`),
    lty = 2,
    size = 0.5,
    alpha = 0.8
  ) +
  facet_grid(
    cols = vars(`mating distance`),
    rows = vars(`dispersal function`),
    labeller = global_labeller,
    scales = "free_y"
  ) + theme(strip.text.x = element_text(size = 6)) +
  scale_color_manual(values = met.brewer("Tsimshian", length(comp_dists)))  +
  labs(col = "competition \ndistance",
       title = "Mating and competition distance",
       x = "distance",
       y = "density") -> curves_plot1



# 2. Dispersal function
all_connections_l %>%
  filter(
    `competition distance` == 0,
    `mating distance` == min(mat_dists),
    sigma == 1,
    simplified == "unsimplified"
  ) %>%
  ggplot(aes(x = distance, col = `dispersal function`), size = 1) +
  geom_line(stat = "density") +
  lims(x = c(0, 8)) +
  theme_minimal() +
  stat_function(
    fun = drayleigh,
    args = list(scale = 1),
    alpha = 0.8,
    lty = 2,
    col = "black"
  ) +
  theme(strip.text.x = element_text(size = 6)) +
  labs(
    col = "dispersal function",
    y = "density",
    x = "distance",
    title = "Dispersal function",
    subtitle = "Simulated"
  ) -> curves_plot2


# theoretical plots
ggplot() +
  theme_minimal() +
  stat_function(
    fun = dfoldnorm,
    args = list(mean = 0, sd = 1),
    alpha = 1,
    lty = 1,
    aes(col = "normal")
  ) + stat_function(
    fun = drayleigh,
    args = list(scale = 1),
    alpha = 1,
    lty = 1,
    aes(col = "brownian")
  ) +  stat_function(
    fun = dcauchy,
    args = list(scale = 1),
    alpha = 1,
    lty = 1,
    aes(col = "cauchy")
  ) +   stat_function(
    fun = dunif,
    args = list(min = 0, max = 1),
    alpha = 1,
    lty = 1,
    aes(col = "uniform")
  ) + stat_function(
    fun = dexp,
    args = list(rate = 1),
    alpha = 1,
    lty = 1,
    aes(col = "exponential")
  ) +
  stat_function(
    fun = drayleigh,
    args = list(scale = 1),
    alpha = 1,
    lty = 2,
    col = "black"
  )  +
  lims(x = c(0, 8)) +
  labs(
    col = "distribution",
    y = "density",
    x = "distance",
    title = "Dispersal function",
    subtitle = "Theoretical"
  ) +
  scale_color_manual(
    values = c(
      "brownian" = hue_pal()(5)[1],
      "normal" = hue_pal()(5)[2],
      "cauchy" = hue_pal()(5)[3],
      "exponential" = hue_pal()(5)[4],
      "uniform" = hue_pal()(5)[5]
    )
  ) -> theory_plot

ggarrange(curves_plot2,
          theory_plot,
          ncol = 2,
          common.legend = T)  -> dispfun_plot

curves_plot2 + coord_cartesian(xlim = c(4, 8), ylim = c(0, 0.01)) -> curves_plot2.tails
theory_plot + coord_cartesian(xlim = c(4, 8), ylim = c(0, 0.01)) -> theory_plot.tails

ggarrange(curves_plot2.tails,
          theory_plot.tails,
          ncol = 2,
          common.legend = T)  -> dispfun_plot.tails

# dispersal distance plot
all_connections_l %>% filter(`competition distance` == 0,
                             `mating distance` == 1,
                             simplified == "unsimplified") %>%
  ggplot() +
  geom_line(stat = "density",
            aes(x = distance, col = as.factor(sigma)),
            size = 0.9) +
  facet_grid(
    cols = vars(`dispersal function`),
    scales = "free",
    labeller = global_labeller
  ) +
  theme_minimal() +
  scale_color_manual(values = met.brewer("Tsimshian",
                                         length(disp_dists)))  +
  labs(colour = "sigma", title = "Dispersal distance (sigma)") +
  lims(x = c(0, 15)) -> dispdist_plot

ggarrange(
  curves_plot1,
  dispdist_plot,
  dispfun_plot,
  ncol = 1,
  heights = c(2.5, 1, 1.5),
  labels = "auto"
) -> curves_plots

## Plotting statistics of these distributions
# Calculate variance and kurtosis
all_connections_stats <-
  all_connections_l %>% filter(simplified == "simplified") %>%
  group_by(
    `mating distance`,
    `competition distance`,
    sigma,
    `dispersal function`,
    simplified,
    rep,
    sim
  ) %>%
  summarise(
    `mean distance` = mean(`distance`),
    `variance (x displacement)` = var(`x displacement`),
    `kurtosis (x displacement)` = kurtosis(`x displacement`)
  ) %>%
  pivot_longer(
    cols = c(
      "mean distance",
      "variance (x displacement)",
      "kurtosis (x displacement)"
    ),
    names_to = "statistic",
    values_to = "value"
  )
# order factor
all_connections_stats$statistic <-
  factor(
    all_connections_stats$statistic,
    levels = c(
      "mean distance",
      "variance (x displacement)",
      "kurtosis (x displacement)"
    )
  )

# Effect of mating and competition distance
all_connections_stats %>% filter(sigma == 1) %>%
  ggplot(aes(x = `mating distance`, y = value, col = `dispersal function`)) +
  geom_line(
    stat = "smooth",
    method = "loess",
    alpha = 0.5,
    size = 0.5,
    se = F
  ) +
  geom_point(size = 0.5) + scale_x_log10() +
  facet_grid(
    rows = vars(statistic),
    cols = vars(`competition distance`, sigma),
    labeller = global_labeller,
    scales = "free_y"
  ) +
  theme_minimal() +
  theme(strip.text.x = element_text(size = 6)) +
  labs(title = "Mating and competition distance") -> mating_competition

# Effect of dispersal distance
all_connections_stats %>% filter(`mating distance` == min(mat_dists),
                                 `competition distance` == min(comp_dists)) %>%
  ggplot(aes(x = sigma, y = value, col = `dispersal function`)) +
  geom_line(
    stat = "smooth",
    method = "lm",
    size = 0.5,
    se = F,
    alpha = 0.5
  ) +
  geom_point(size = 0.5) +
  facet_grid(
    rows = vars(statistic),
    cols = vars(`competition distance`, `mating distance`),
    labeller = global_labeller,
    scales = "free_y"
  ) +
  theme_minimal() +
  theme(strip.text.x = element_text(size = 6)) +
  labs(title = "Dispersal distance") -> dispersal_dist

ggarrange(
  mating_competition,
  dispersal_dist,
  nrow = 1,
  labels = "auto",
  widths = c(2, 1),
  common.legend = T
) -> statistics_plots

### Plotting locations of individuals
## 1. Across a mating/competition distance grid
tree_data %>%
  rename(
    "competition distance" = "comp_dist",
    "mating distance" = "mat_dist",
    "dispersal function" = "fun"
  ) %>%
  filter(
    simplified == "unsimplified",
    `dispersal function` == "brownian",
    sigma == 1,
    time == max(time)
  ) %>%
  ggplot(aes(col = rep)) +
  geom_sf(show.legend = F) +
  facet_grid(
    cols = vars(`competition distance`),
    rows = vars(`mating distance`),
    labeller = global_labeller
  ) +
  labs(colour = "simulation repeat",
       x = "Eastings",
       y = "Northings") +
  theme_minimal() +
  theme(strip.text.x = element_text(size = 6),
        strip.text.y = element_text(size = 6)) -> points_plot1

# showing the connections between clusters according to mating distance
all_connections_l %>%
  filter(
    rep == 1,
    `competition distance` == max(comp_dists),
    `dispersal function` == "brownian",
    sigma == 1,
    child_time == max(child_time)
  ) %>%
  ggplot(aes()) +
  geom_sf(aes(geometry = child_location), col = "grey") +
  geom_sf(aes(geometry = connection), col = "lightpink4") +
  facet_grid(
    rows = vars(`competition distance`),
    cols = vars(`mating distance`),
    labeller = global_labeller
  ) +
  theme(legend.position = "none") +
  labs(x = "Eastings",
       y = "Northings") +
  theme_minimal() +
  theme(strip.text.x = element_text(size = 6),
        strip.text.y = element_text(size = 6))  -> mating_distance_clusters

## 2. In time groups
# a) Mating distance
tree_data %>%
  ungroup() %>%
  rename(
    "competition distance" = "comp_dist",
    "mating distance" = "mat_dist",
    "dispersal function" = "fun"
  ) %>%
  dplyr::filter(
    simplified == "unsimplified",
    `dispersal function` == "brownian",
    sigma == 1,
    `competition distance` == 0.2,
    `mating distance` == min(`mating distance`) |
      `mating distance` == max(`mating distance`),
    rep == 1
  ) %>% group_by(sim) %>%
  mutate(timegroup = cut(time, breaks = 4)) %>%
  ggplot(aes(col = timeoriginal)) +
  geom_sf() +
  facet_grid(
    cols = vars(timegroup),
    rows = vars(`mating distance`),
    labeller = global_labeller
  ) +
  theme_minimal() +
  scale_colour_gradientn(colours = met.brewer(name = "Hokusai1")) +
  theme_minimal() +
  theme(strip.text.x = element_text(size = 6),
        strip.text.y = element_text(size = 6)) +
  labs(
    x = "Eastings",
    y = "Northings",
    colour = "time (from root)",
    caption = sprintf(
      'Competition distance: %s; Sigma: %s; Dispersal function: %s',
      0.2,
      1,
      '"brownian"'
    )
  ) -> timegroups_md
# b) competition distance
tree_data %>%
  ungroup() %>%
  rename(
    "competition distance" = "comp_dist",
    "mating distance" = "mat_dist",
    "dispersal function" = "fun"
  ) %>%
  dplyr::filter(
    simplified == "unsimplified",
    `dispersal function` == "brownian",
    sigma == 1,
    `mating distance` == min(mat_dists),
    `competition distance` == min(comp_dists) |
      `competition distance` == max(comp_dists),
    rep == 1
  ) %>%
  mutate(timegroup = cut(timeoriginal, breaks = 4)) %>%
  ggplot(aes(col = timeoriginal)) +
  geom_sf() +
  facet_grid(
    cols = vars(timegroup),
    rows = vars(`competition distance`),
    labeller = global_labeller
  ) +
  theme_minimal() +
  scale_colour_gradientn(colours = met.brewer(name = "Hokusai1")) +
  theme_minimal() +
  theme(strip.text.x = element_text(size = 6),
        strip.text.y = element_text(size = 6)) +
  labs(
    x = "Eastings",
    y = "Northings",
    colour = "time (from root)",
    caption = sprintf(
      'Mating distance: %s; Sigma: %s; Dispersal function: %s',
      min(mat_dists),
      1,
      '"brownian"'
    )
  ) -> timegroups_cd

ggarrange(
  timegroups_md,
  timegroups_cd,
  common.legend = T,
  ncol = 1,
  labels = "auto",
  heights = c(1, 1)
) -> points_plot2

######### Exporting all the plots

curves_plots %>%
  ggsave(
    filename = paste0("figs/", as.character(Sys.Date()) , "curves_plot.pdf"),
    device = "pdf",
    height = 9,
    width = 7
  )

dispfun_plot.tails %>%
  ggsave(
    filename = paste0("figs/", as.character(Sys.Date()) , "curves_plot_tails.pdf"),
    device = "pdf",
    height = 2,
    width = 7
  )

statistics_plots %>%
  ggsave(
    filename = paste0("figs/", as.character(Sys.Date()), "statistics_plot.pdf"),
    device = "pdf",
    height = 7,
    width = 7
  )

ggarrange(
  points_plot1,
  mating_distance_clusters,
  ncol = 1,
  heights = c(3, 1),
  labels = "auto"
) %>%
  ggsave(
    filename = paste0("figs/", as.character(Sys.Date()), "points.pdf"),
    device = "pdf",
    height = 7,
    width = 7
  )

points_plot2 %>%
  ggsave(
    filename = paste0("figs/", as.character(Sys.Date()), "points_2.pdf"),
    device = "pdf",
    height = 9,
    width = 7
  )
