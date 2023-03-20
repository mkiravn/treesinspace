######## script for treesinspace
### load libraries and set up
devtools::load_all("~/treesinspace")
library(slendr)
library(testthat)
library(cowplot)
setup_env()
check_env()
set.seed(206)

##### Part 1: exploration of parameters
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

res <- run_slendr_simulation(pars,map=map,pathout="trial")

trees <- res$simplified
trees_us <- res$unsimplified

tree_data <- convert_trees(trees,pars) %>% rbind(convert_trees(trees_us,pars))

# save the tree data
tree_data <- tree_data %>%
  group_by(sim) %>%
  # adjust so that time starts at 0 (for recapitated trees)
  mutate(timeoriginal = time,
         time = time - min(time))

all_connections_l <- retrieve_connections(trees,trees_us,pars)

######### Plotting
global_labeller <- labeller(
  .default = label_value,
  `mating distance` = label_both,
  `competition distance` = label_both,
  sigma = label_both
)

## Plotting curves of parent-offspring distance
# 1. Mating distance and competition
all_connections_l %>%
  dplyr::filter(simplified == "unsimplified",
                sigma == 1,
                `dispersal function`=="brownian") %>%
  ggplot(aes(x = distance,
             col = as.factor(`competition distance`))) +
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
  dplyr::filter(
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

curves_plot2 + coord_cartesian(xlim = c(4, 8), ylim = c(0, 0.01)) -> curves_plot2.tails
theory_plot + coord_cartesian(xlim = c(4, 8), ylim = c(0, 0.01)) -> theory_plot.tails

ggarrange(curves_plot2,
          theory_plot,
          curves_plot2.tails,
          theory_plot.tails,
          common.legend = T)  -> dispfun_plot

## Plotting statistics of these distributions
# Calculate variance and kurtosis
all_connections_stats <-
  all_connections_l %>%
  dplyr::filter(simplified == "unsimplified") %>%
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
    `variance (distance)` = var(`distance`),
    `kurtosis (distance)` = kurtosis(`distance`)
  ) %>%
  pivot_longer(
    cols = c(
      "variance (distance)",
      "kurtosis (distance)"
    ),
    names_to = "statistic",
    values_to = "value"
  )
# order factor
all_connections_stats$statistic <-
  factor(
    all_connections_stats$statistic,
    levels = c("variance (distance)",
      "kurtosis (distance)"
    )
  )

# Effect of mating and competition distance
all_connections_stats %>%
  dplyr::filter(sigma == 1,`dispersal function`=="brownian") %>%
  ggplot(aes(x = `mating distance`, y = value, col = as.factor(`competition distance`))) +
  geom_boxplot(aes(group=`mating distance`+`competition distance`, y = value),alpha=0.8,width=0.1) +
  geom_line(
    stat = "smooth",
    method = "loess",
    alpha = 0.5,
    size = 0.5,
    se = F
  ) +
  geom_point(size = 0.5) +
  scale_x_log10() +
  facet_grid(
    rows = vars(statistic),
    labeller = global_labeller,
    scales = "free_y"
  ) +
  theme_minimal() +
  theme(strip.text.x = element_text(size = 6)) +
  scale_color_manual(values = met.brewer("Tsimshian", length(comp_dists)))  +
  labs(title = "Mating and competition distance",
       col = "competition \ndistance",y="") -> mating_competition

ggarrange(curves_plot1,
          mating_competition,
          nrow=2,
          common.legend = T) -> md_plot

dispfun_plot %>%
  ggsave(
    filename = paste0("figs/", as.character(Sys.Date()), "dispfun_plot.pdf"),
    device = "pdf",
    height = 7,
    width = 7
  )

md_plot %>%
  ggsave(
    filename = paste0("figs/", as.character(Sys.Date()), "md_plot.pdf"),
    device = "pdf",
    height = 7,
    width = 7
  )
