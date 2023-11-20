######## script for treesinspace
library(latex2exp)
library(spatstat)
set.seed(206)

##### Part 1: exploration of parameters. Illustrate positions of individuals.
# defining parameters
Ns <- 200 # number of individuals in population
mat_dists <- c(.2, 2, 20, 200) # mating distance
comp_dists <- c(0, .2, 2, 20, 200) # competition distance
disp_dists <- c(1) # dispersal distance (sigma)
disp_funs <-
  c("brownian") # dispersal kernel
reps <-
  c(1) # number of simulation runs for each parameter combination
ngens <- 10 # number of generations
fileout <- "illustration"
pars <- set_slendr_pars(Ns,mat_dists,comp_dists,disp_dists,reps,ngens,fileout)

# defining a world
map <- world(
  xrange = c(-25, 25),
  # min-max longitude
  yrange = c(-25, 25),
  # min-max latitude
  landscape = "blank"
)

res <- run_slendr_simulation(pars,map=map,pathout=fileout)

trees <- res$simplified
trees_us <- res$unsimplified

tree_data <- convert_trees(trees,pars) %>% rbind(convert_trees(trees_us,pars))

# save the tree data
tree_data <- tree_data %>%
  group_by(sim) %>%
  # adjust so that time starts at 0 (for recapitated trees)
  mutate(timeoriginal = time,
         time = time - min(time))

all_connections_mc <- retrieve_connections(trees,trees_us,pars,nice_names=T)

all_connections_mc <- all_connections_mc %>% mutate(parent_time = - (max(parent_time) - parent_time))

all_connections_mc %>% write_delim("illustration_conns.tsv",delim="\t")

######### Plotting
global_labeller <- labeller(
  .default = label_value,
  `mating distance` = label_both,
  `competition distance` = label_both,
  sigma = label_both
)

all_connections_mc <-  all_connections_mc %>% mutate(parent_time = - parent_time)

all_connections_mc %>%
  dplyr::filter(simplified == "unsimplified",
         `dispersal function` == "brownian",
         sigma == 1) %>%
  ggplot(aes(col = parent_time,geometry=parent_location)) +
  geom_sf() +
  facet_grid(
    cols = vars(`competition distance`),
    rows = vars(`mating distance`),
    labeller = label_both
  ) +
  #geom_sf(aes(geometry=connection),alpha=0.4,col="grey",size=0.5)+
  scale_colour_gradientn(colours = met.brewer(name = "Hokusai2")) +
  theme_minimal()+
  labs(color="time in past",x="Eastings",y="Northings")+
  theme(strip.text.x = element_text(size = 6),
        strip.text.y = element_text(size = 6))-> points_plot1

points_plot1 %>%
  ggsave(
    filename = paste0("figs/", as.character(Sys.Date()), "mdcd_plot.pdf"),
    device = "pdf",
    height = 5,
    width = 7
  )


