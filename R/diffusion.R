# defining parameters
Ns <- 2000 # number of individuals in population
mat_dists <- 0.1 # mating distance
comp_dists <- 0.2 # competition distance
disp_dists <- 1 # dispersal distance (sigma)
disp_funs <-
  "brownian" # dispersal kernel
reps <-
  1 # number of simulation runs for each parameter combination
ngens <- 100 # number of generations
fileout <- "part2"
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

# main plots
all_connections %>%
  dplyr::filter(simplified=="unsimplified") %>%
  ggplot(aes(x=distance,lty="simulated",col="simulated")) +
  geom_density(show.legend = F) +
  stat_function(aes(lty="theoretical"),fun = drayleigh, args = list(scale=1),alpha=0.8,col="black") +
  theme_minimal() +
  labs(y="density",x="distance",lty="") -> absdist

all_connections %>%
  dplyr::filter(simplified == "unsimplified") %>%
  select(x = `x displacement`, y =`y displacement`) %>%
  pivot_longer(cols = c("x", "y"),
               names_to = "dimension",
               values_to = "displacement") %>%
  ggplot(aes(x = displacement, lty = "simulated", col = dimension)) +
  geom_density() +
  stat_function(aes(lty = "theoretical"),
                fun = dnorm,
                args = list(mean = 0, sd = 1),
                alpha = 1,
                col = "black"
  ) +
  theme_minimal() +
  labs(y = "density", x = "displacement", lty = "") -> xandy

# adding qqplots
all_connections %>%
  dplyr::filter(simplified=="unsimplified") %>%
  ggplot() +
  stat_qq(aes(sample=distance,col="distance"),distribution = qrayleigh,alpha=0.5)  +
  stat_qq_line(aes(sample=distance),
               distribution = qrayleigh,
               alpha=0.8,
               lty=2,
               col="black",
               dparams = list(scale=1))  +
  labs(y="simulated",x="theoretical") +
  theme_minimal() -> absdist_qq
all_connections %>%
  dplyr::filter(simplified=="unsimplified") %>%
  ggplot() +
  stat_qq(aes(sample=`x displacement`,col="x displacement"),distribution = qnorm,alpha=0.5)  +
  stat_qq(aes(sample=`y displacement`,col="y displacement"),distribution = qnorm,alpha=0.5)  +
  stat_qq_line(aes(sample=`y displacement`),
               distribution = qnorm,
               alpha=0.8,
               lty=2,
               col="black",
               dparams = list(mean=0,sd=1))  +
  labs(y="simulated",x="theoretical") +
  theme_minimal() -> xandy_qq
ggarrange(absdist,absdist_qq,xandy,xandy_qq,common.legend = T) -> brownian_plots

