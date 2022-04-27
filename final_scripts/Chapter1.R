# Loading libraries
devtools::load_all("~/treesinspace")
library(tidyverse)
library(ape)
library(sf)
library(MetBrewer)
library(ggrepel)
library(ggpubr)
library(ggtree)
library(moments)
library(phangorn)
library(apTreeshape)
library(VGAM)
devtools::load_all("~/slendr")
setup_env()
check_env()
set.seed(206)

######### Here we define parameters

# defining parameters
Ns <- 100 # number of individuals in population
mat_dists <- c(1, 5, 10, 1000) # mating distance
comp_dists <- c(0, 0.2, 20, 40) # competition distance
disp_dists <- c(1:3) # dispersal distance (sigma)
disp_funs <-
  c("brownian", "normal", "cauchy", "uniform") # dispersal kernel
reps <-
  c(1:3) # number of simulation runs for each parameter combination
ngens <- 500 # number of generations

# generating a parameter table
pars <- expand.grid(
  N = Ns,
  mat_dist = mat_dists,
  comp_dist = comp_dists,
  disp_fun = disp_funs,
  disp_dist = disp_dists,
  rep = reps,
  ngen = ngens
) %>%
  rownames_to_column("sim")


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
    "POP",
    time = 1,
    N = pars[row, "N"],
    center = c(0, 0),
    radius = 50,
    map = map,
    mate_dist = pars[row, "mat_dist"],
    competition_dist = pars[row, "comp_dist"],
    dispersal_fun = as.character(pars[row, "disp_fun"]),
    dispersal_dist = pars[row, "disp_dist"]
  )

  # compile and run the model
  model <- compile(
    populations = pop,
    generation_time = 1,
    sim_length = ngens,
    resolution = 1,
    # resolution in "distance units per pixel"
    path = "~/Desktop/test-model",
    overwrite = TRUE
  )

  # sampling strategy
  samples <- sampling(model, times = ngens, list(pop, 25))

  # simulate
  slim(
    model,
    sampling = samples,
    # simulate only a single locus
    sequence_length = 1,
    recombination_rate = 0,
    method = "batch",
    # change to "gui" to execute the model in SLiMgui
    random_seed = i,
    verbose = FALSE,
    retainall = TRUE,
    burnin = 0,
  )

  # extract simplified trees
  ts <- ts_load(model, simplify = TRUE) %>%
    ts_recapitate(Ne = pars[row, "N"],
                  recombination_rate = 0) %>%
    ts_simplify() %>%
    ts_mutate(mutation_rate = 1e-4)

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
      disp_dist =
        pars[row, "disp_dist"],
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
      disp_dist =
        pars[row, "disp_dist"],
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

all_connections <- data.frame() # this will hold the connections

for (row in c(1:dim(pars)[1])) {
  print(paste("Connecting tree number:", row))
  all_connections_i <- conns[[row]][[1]] %>%
    as.data.frame() %>%
    mutate(simplified = "unsimplified") %>%
    cbind(pars[row, ])
  all_connections_i <- conns_simp[[row]][[1]] %>%
    as.data.frame() %>%
    mutate(simplified = "simplified") %>%
    cbind(pars[row, ]) %>% rbind(all_connections_i)
  all_connections <- rbind(all_connections, all_connections_i)
}


######### Plotting

## Plotting curves of parent-offspring distance
# 1. Mating distance
parameter_curves_plot <- all_connections %>%
  filter(comp_dist == 0,
         disp_dist == 1,
         disp_fun == "brownian",
         simplified == "unsimplified") %>%
  ggplot(aes(x = dist, col = as.factor(mat_dist))) +
  geom_density(alpha = 0.5) +
  lims(x = c(0, 30)) +
  theme_minimal() +
  labs(col = "mating distance") +
  # 2. Competition distance
  all_connections %>%
  filter(
    mat_dist == min(mat_dists),
    disp_dist == 1,
    disp_fun == "brownian",
    simplified == "unsimplified"
  ) %>%
  ggplot(aes(x = dist, col = as.factor(comp_dist))) +
  geom_density(alpha = 0.5) +
  lims(x = c(0, 30)) +
  theme_minimal() +
  labs(col = "competition distance") +
  # 3. Dispersal distance
  all_connections %>%
  filter(
    comp_dist == 0,
    mat_dist == min(mat_dists),
    disp_fun == "brownian",
    simplified == "unsimplified"
  ) %>%
  ggplot(aes(x = dist, col = as.factor(disp_dist))) +
  geom_density(alpha = 0.5) +
  lims(x = c(0, 30)) +
  theme_minimal() +
  labs(col = "dispersal distance") -> curves_plot1
  # 4. Dispersal function
all_connections %>%
  filter(comp_dist == 0,
         mat_dist == min(mat_dists),
         disp_dist == 1,
         simplified == "unsimplified") %>%
  ggplot(aes(x = dist, col = disp_fun)) +
  geom_density(alpha = 0.5) +
  lims(x = c(0, 20)) +
  theme_minimal() +
  labs(col = "dispersal function") -> curves_plot2

## Plotting statistics of these distributions
# Calculate variance and kurtosis
all_connections_stats <- all_connections %>%
  group_by(mat_dist, comp_dist, disp_dist, disp_fun, simplified, rep, sim) %>%
  summarise(
    mean = mean(dist),
    variance = var(dist),
    kurtosis = kurtosis(dist)
  ) %>%
  pivot_longer(
    cols = c("mean", "variance", "kurtosis"),
    names_to = "statistic",
    values_to = "value"
  )

# 1. Mating distance
all_connections_stats %>%
  filter(comp_dist == 0,
         disp_dist == 1,
         disp_fun == "brownian",
         simplified == "unsimplified") %>%
  ggplot(aes(x = mat_dist, y = value)) +
  geom_point() +
  geom_line() +
  facet_wrap( ~ statistic, labeller = label_both,scales="free") +
  theme_minimal() +
  # 2. Competition distance
  all_connections_stats %>%
  filter(
    mat_dist == min(mat_dists),
    disp_fun == "brownian",
    simplified == "unsimplified"
  ) %>%
  ggplot(aes(x = comp_dist, y = value)) +
  geom_point() +
  geom_line() +
  facet_wrap( ~ statistic, labeller = label_both,scales="free") +
  theme_minimal() +
  # 3. Dispersal distance
  all_connections_stats %>%
  filter(
    comp_dist == 0,
    mat_dist == min(mat_dists),
    disp_fun == "brownian",
    simplified == "unsimplified"
  ) %>%
  ggplot(aes(x = disp_dist, y = value)) +
  geom_point() +
  geom_line() +
  facet_wrap( ~ statistic, labeller = label_both,scales="free") +
  theme_minimal() +
  # 4. Dispersal function
  all_connections_stats %>%
  filter(comp_dist == 0,
         mat_dist == min(mat_dists),
         disp_dist == 1,
         simplified == "unsimplified") %>%
  ggplot(aes(x = disp_fun, y = value)) +
  geom_point() +
  facet_wrap( ~ statistic, labeller = label_both,scales="free") +
  theme_minimal() -> curves_plot3

### Plotting locations of individuals
## 1. Across a mating/competition distance grid
tree_data %>%
  filter(simplified == "unsimplified",
         fun == "brownian",
         disp_dist == 1,
         time == max(time)) %>%
  ggplot(aes(col = rep)) +
  geom_sf() +
  facet_grid(
    cols = vars(comp_dist),
    rows = vars(mat_dist),
    labeller = label_both
  ) +
  theme_minimal() -> points_plot1
## 2. In time groups
# a) Mating distance
tree_data %>% ungroup() %>%
  filter(
    simplified == "unsimplified",
    fun == "brownian",
    disp_dist == 1,
    comp_dist == 0.2,
    mat_dist == min(mat_dist) |
      mat_dist == max(mat_dist),
    rep == 1
  ) %>%
  mutate(timegroup = cut(time, breaks = 4)) %>%
  ggplot(aes(col = time)) +
  geom_sf() +
  facet_grid(
    cols = vars(timegroup),
    rows = vars(mat_dist),
    labeller = labeller(mat_dist = label_both, timegroup = label_value
  )) +
  theme_minimal() +
  scale_colour_gradientn(colours = met.brewer(name = "Hokusai2")) +
  theme_minimal() +
  labs(x = "", y = "") -> timegroups_md
# b) competition distance
tree_data %>% ungroup() %>%
  filter(
    simplified == "unsimplified",
    fun == "brownian",
    disp_dist == 1,
    mat_dist == min(mat_dists),
    comp_dist == min(comp_dists) |
      comp_dist == max(comp_dists),
    rep == 1
  ) %>%
  mutate(timegroup = cut(timeoriginal, breaks = 4)) %>%
  ggplot(aes(col = time)) +
  geom_sf() +
  facet_grid(
    cols = vars(timegroup),
    rows = vars(comp_dist),
    labeller = labeller(comp_dist = label_both, timegroup = label_value)
  ) +
  theme_minimal() +
  scale_colour_gradientn(colours = met.brewer(name = "Hokusai2")) +
  theme_minimal() +
  labs(x = "", y = "") -> timegroups_cd

ggarrange(timegroups_md, timegroups_cd, common.legend = T,ncol=1,labels="auto") -> points_plot2

######### Exporting all the plots

ggarrange(curves_plot1, curves_plot2, curves_plot3, ncol = 1,labels="auto") %>%
  ggsave(
    filename = paste0("figs/", as.character(Sys.Date()), "-curves_plot.pdf"),
    device = "pdf",
    height = 8,
    width = 10
  )

ggarrange(points_plot1, points_plot2, ncol = 2,labels="auto") %>%
  ggsave(
    filename = paste0("figs/", as.character(Sys.Date()), "-points_plot1.pdf"),
    device = "pdf",
    height = 6,
    width = 13
  )





