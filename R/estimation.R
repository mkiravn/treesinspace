set.seed(206)

# defining parameters
n <- 2000
Ns <- n
mat_dists <- c(0.2, 0.5, 1, 2, 5) # mate distance
comp_dists <- c(0.2) # competition distance
disp_dists <- c(1) # mean dispersal distance
disp_funs <- c("brownian", "cauchy")
reps <- c(1:5)
ngens <- 8000
fileout <- "part2a"
pars <- set_slendr_pars(Ns,mat_dists,comp_dists,disp_dists,reps,ngens,fileout)

# defining a world
map <- world(
  xrange = c(-25, 25),
  # min-max longitude
  yrange = c(-25, 25),
  # min-max latitude
  landscape = "blank"
)


res <- run_slendr_simulation(pars,map=map,pathout="part2")

trees <- res$simplified
trees_us <- res$unsimplified

tree_data <- convert_trees(trees,pars) %>% rbind(convert_trees(trees_us,pars))

# save the tree data
tree_data <- tree_data %>%
  group_by(sim) %>%
  # adjust so that time starts at 0 (for recapitated trees)
  mutate(timeoriginal = time,
         time = time - min(time))

all_connections <- retrieve_connections(trees,trees_us,pars)
