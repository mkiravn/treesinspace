#### Illustrating sampling differences


# libraries
library(ape)
library(tidyverse)
library(sf)
devtools::load_all("../slendr")
devtools::load_all("../treesinspace")


n <- 20

# Tree based simulation
tr <- rtree(n)
tr <- move_nodes(edge_calculator(tr,tree_mode="coalescent",Ne=n),disp_dist = 1)
plot(tr)


# slendr simulation:
# setting some parameters
disp_f <- "brownian"
comp_d <- 3
mate_d <- 100
disp_d <- 1

# how many gens?
gens <- 1000
# how large a pop?
N <- n

# making a landscape...
map <-
  world(
    xrange = c(1, 100),
    yrange = c(1, 100),
    landscape = "blank"
  )
# defining the population
pop <- population(
  "POP",
  time = 1,
  N = N,
  center = c(50, 50),
  radius = 50,
  map = map,
  mate_dist = mate_d,
  competition_dist = comp_d,
  dispersal_fun = disp_f,
  dispersal_dist = disp_d
)

# make the model
model <- compile(
  populations = pop,
  generation_time = 1,
  sim_length = gens,
  resolution = 0.5
)

# who should we sample?
samples <- sampling(model, times = gens, list(pop, 10))
# run model
slim(
  model,
  sequence_length = 1,
  recombination_rate = 0, # one tree locus
  method = "batch",
  random_seed = 50,
  retainall=TRUE, # a new option
  sampling = samples
)

# the full tree
tsu <- ts_load(model)
# its data
data <- tsu %>% ts_data() %>% as.data.frame()
# the connections and branches
res <- ts_connect(tsu)
# the simplified tree
ts <- ts_load(model) %>%
  ts_recapitate(Ne = N, recombination_rate = 0) %>%
  ts_simplify()
# the phylogenetic tree
tree <- ts %>% ts_phylo(i=1)
class(tree) <- "phylo"

#1) simplify and show


#2) don't simplify, and sample descendants of a node in the past


# compute tree-based statistics of genetic diversity
