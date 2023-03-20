######## script for treesinspace
### load libraries and set up
# devtools::load_all("~/treesinspace")
# library(slendr)
# library(testthat)
# library(cowplot)
# setup_env()
# check_env()
library(latex2exp)
set.seed(206)

##### Part 1: exploration of parameters
# defining parameters
Ns <- 50 # number of individuals in population
mat_dists <- c(1, 5, 10, 100) # mating distance
comp_dists <- c(0, 0.2, 20, 40) # competition distance
disp_dists <- c(1,5,10) # dispersal distance (sigma)
disp_funs <-
  c("brownian") # dispersal kernel
reps <-
  c(1) # number of simulation runs for each parameter combination
ngens <- 10 # number of generations
fileout <- "part1"
pars <- set_slendr_pars(Ns,mat_dists,comp_dists,disp_dists,reps,ngens,fileout)

# defining a world
map <- world(
  xrange = c(-25, 25),
  # min-max longitude
  yrange = c(-25, 25),
  # min-max latitude
  landscape = "blank"
)

res <- run_slendr_simulation(pars,map=map,pathout="part1")

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

######### Plotting
global_labeller <- labeller(
  .default = label_value,
  `mating distance` = label_both,
  `competition distance` = label_both,
  sigma = label_both
)


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
  geom_sf(aes(geometry=connection),alpha=0.4,col="grey",size=0.5)+
  scale_colour_gradientn(colours = met.brewer(name = "Hokusai2")) +
  theme_minimal()+
  labs(color="time")+
  theme(strip.text.x = element_text(size = 6),
        strip.text.y = element_text(size = 6))-> points_plot1

points_plot1 %>%
  ggsave(
    filename = paste0("figs/", as.character(Sys.Date()), "mdcd_plot.pdf"),
    device = "pdf",
    height = 5,
    width = 7
  )
