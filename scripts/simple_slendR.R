
# setting some parameters
disp_f <- "brownian"
comp_d <- 1
mate_d <- 1
disp_d <- 1

# how many gens?
gens <- 1000
# how large a pop?
N <- 100

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
samples <- sampling(model, times = gens, list(pop, 25))
# run model
slim(
  model,
  sequence_length = 1,
  recombination_rate = 0, # one tree locus
  method = "batch",
  random_seed = 50,
  retainall=FALSE, # a new option
  sampling = samples
)

# the full tree
tsu <- ts_load(model)
# its data
data <- tsu %>% ts_data() %>% as.data.frame()
# plot node locations
data %>% ggplot(aes(col=time))+geom_sf(aes(geometry=location))
# extract the phylogenetic tree
tree <- tsu %>% ts_simplify() %>% ts_phylo(i=1)
# plot the tree
plot(tree)

