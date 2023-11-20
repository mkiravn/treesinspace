######## script for treesinspace
library(latex2exp)
library(phyloTop)
set.seed(206)

##### Part 1: exploration of parameters
# defining parameters
Ns <- 100 # number of individuals in population
mat_dists <- c(.2, 2, 20, 200)# mating distance
comp_dists <- c(0, .2, 2, 20, 200) # competition distance
disp_dists <- c(1:3) # dispersal distance (sigma)
disp_funs <-
  c("brownian", "normal", "cauchy", "exponential", "uniform") # dispersal kernel
reps <-
  c(1:20) # number of simulation runs for each parameter combination
ngens <- 50 # number of generations
fileout <- "exploration"
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

all_connections_l <- retrieve_connections(trees,trees_us,pars,nice_names=T)

######### Plotting
global_labeller <- labeller(
  .default = label_value,
  `mating distance` = label_both,
  `competition distance` = label_both,
  sigma = label_both
)
# rename normal to half-normal
all_connections_l <- all_connections_l %>%
  mutate(`dispersal function` = str_replace(`dispersal function`, "normal", "half-normal"))

all_connections_l %>% write_delim("exploration_conns.tsv",delim="\t")
# in case you need to read in the table all_connections_l <- read.delim("exploration_conns.tsv",check.names=FALSE)

# Plot of dispersal function
all_connections_l %>%
  dplyr::filter(
    `competition distance` == min(comp_dists),
    `mating distance` == min(mat_dists),
    sigma == 1,
    simplified == "unsimplified"
  ) %>%
  ggplot(aes(x = distance, col = `dispersal function`), size = 1) +
  geom_line(stat = "density",bw=0.1) +
  lims(x = c(0, 8)) +
  theme_minimal() +
  theme(strip.text.x = element_text(size = 6)) +
  labs(
    col = "dispersal function",
    y = "density",
    x = "distance",
    subtitle = "Sampled from spatial simulations"
  ) -> curves_plot2


# plot of theoretical curves
ggplot() +
  theme_minimal() +
  stat_function(
    fun = dfoldnorm,
    args = list(mean = 0, sd = 1),
    alpha = 1,
    lty = "longdash",
    aes(col = "half-normal")
  ) +  stat_function(
    fun = dcauchy,
    args = list(scale = 1),
    alpha = 1,
    lty = "longdash",
    aes(col = "cauchy")
  ) + stat_function(
    fun = drayleigh,
    args = list(scale = 1),
    alpha = 1,
    lty = "longdash",
    aes(col = "brownian")
  ) +   stat_function(
    fun = dunif,
    args = list(min = 0, max = 1),
    alpha = 1,
    lty = "longdash",
    aes(col = "uniform")
  ) + stat_function(
    fun = dexp,
    args = list(rate = 1),
    alpha = 1,
    lty = "longdash",
    aes(col = "exponential")
  ) +
  lims(x = c(0, 8)) +
  labs(
    col = "distribution",
    y = "density",
    x = "distance",
    subtitle = "Theoretical"
  ) +
  scale_color_manual(
    values = c(
      "brownian" = hue_pal()(5)[1],
      "half-normal" = hue_pal()(5)[4],
      "cauchy" = hue_pal()(5)[2],
      "exponential" = hue_pal()(5)[3],
      "uniform" = hue_pal()(5)[5]
    )
  ) -> theory_plot
# plot of tails
curves_plot2 +
  coord_cartesian(xlim = c(4, 8), ylim = c(0, 0.01)) -> curves_plot2.tails
theory_plot + coord_cartesian(xlim = c(4, 8), ylim = c(0, 0.01)) -> theory_plot.tails
# arranging into one plot
ggarrange(curves_plot2,
          theory_plot,
          curves_plot2.tails,
          theory_plot.tails,
          heights = c(4,3),
          common.legend = T)  -> dispfun_plot
# writing out
dispfun_plot %>%
  ggsave(
    filename = paste0("figs/", as.character(Sys.Date()), "dispfun_plot.pdf"),
    device = "pdf",
    height = 4,
    width = 7
  )



## Plotting statistics of these distributions
# Calculate variance
all_connections_stats <-
  all_connections_l %>%
  dplyr::filter(simplified == "unsimplified") %>%
  group_by(
    `mating distance`,
    `competition distance`,
    sigma,
    `dispersal function`,
    simplified,
  ) %>%
  summarise(
    `mean (distance)` = mean(`distance`),
    `variance (distance)` = var(`distance`),
  ) %>%
  pivot_longer(
    cols = c(
      "mean (distance)",
      "variance (distance)",
    ),
    names_to = "statistic",
    values_to = "value"
  )
# output a table of statistics for these distributions
stats <- all_connections_stats %>%
  ungroup() %>%
  select(`dispersal function`) %>%
  distinct()
# check these
stats$expvariance <- NA
stats[stats$`dispersal function`==
        "brownian","expvariance"] <- (4-pi)/2
stats[stats$`dispersal function`==
        "half-normal","expvariance"] <-  1 - 2/pi
# undefined for cauchy- set to zero
stats[stats$`dispersal function`==
"cauchy","expvariance"] <- 0
stats[stats$`dispersal function`==
        "exponential","expvariance"] <- 1
stats[stats$`dispersal function`==
        "uniform","expvariance"] <- 1/12

stats$`theoretical variance` <- NA
stats[stats$`dispersal function`==
        "brownian","theoretical variance"] <- "$\\sigma^2(4-\\pi)/2$"
stats[stats$`dispersal function`==
        "half-normal","theoretical variance"] <-  "$\\sigma^2-\\frac{2}{\\pi}\\sigma^2$"
# undefined for cauchy- set to zero
stats[stats$`dispersal function`==
        "cauchy","theoretical variance"] <- "undefined"
stats[stats$`dispersal function`==
        "exponential","theoretical variance"] <- "$\\sigma^2$"
stats[stats$`dispersal function`==
        "uniform","theoretical variance"] <- "$\\sigma^2/12$"

stats$`theoretical mean` <- NA
stats[stats$`dispersal function`==
        "brownian","theoretical mean"] <- "$\\sigma\\sqrt{\\pi/2}$"
stats[stats$`dispersal function`==
        "half-normal","theoretical mean"] <-  "$\\sigma\\sqrt{2/\\pi}$"
# undefined for cauchy- set to zero
stats[stats$`dispersal function`==
        "cauchy","theoretical mean"] <- "undefined"
stats[stats$`dispersal function`==
        "exponential","theoretical mean"] <- "$\\sigma$"
stats[stats$`dispersal function`==
        "uniform","theoretical mean"] <- "$\\sigma/2$"

stats$`parametrisation` <- NA
stats[stats$`dispersal function`==
        "brownian","parametrisation"] <- "Distance in x and y dimensions drawn independently from N(0,$\\sigma^2$). Distance follows Rayleigh($\\sigma$)"
stats[stats$`dispersal function`==
        "half-normal","parametrisation"] <- "Angle drawn uniformly, distance drawn from N(0,$\\sigma^2$). Distance follows folded normal distribution"
# undefined for cauchy- set to zero
stats[stats$`dispersal function`==
        "cauchy","parametrisation"] <- "Angle drawn uniformly, distance drawn from Cauchy(scale=$\\sigma$,location=0)"
stats[stats$`dispersal function`==
        "exponential","parametrisation"] <- "Angle drawn uniformly, distance drawn from Exp($1/\\sigma$)"
stats[stats$`dispersal function`==
        "uniform","parametrisation"] <- "Angle drawn uniformly, distance drawn from U(0,$\\sigma$)"
library(kableExtra)

stats %>% select(`dispersal function`,parametrisation,`theoretical mean`,`theoretical variance`) %>%
  knitr::kable( escape = FALSE, "latex") %>%
  column_spec(2, width = "1in") %>%
  kableExtra::save_kable(paste0("figs/", as.character(Sys.Date()), "table.tex"))

all_connections_stats_evar <- all_connections_stats %>%
  dplyr::filter(sigma==1) %>%
  pivot_wider(names_from="statistic",values_from="value") %>%
  full_join(stats,by="dispersal function") %>%
  mutate(`excess variance`=`variance (distance)`-expvariance)

# plot of excess variance
all_connections_stats_evar %>%
  ggplot(aes(x=as.factor(`mating distance`),
             y=as.factor(`competition distance`),
             fill=log(`excess variance`+1)))+
  geom_tile() +
  facet_wrap(~`dispersal function`,nrow = 1) +
  scale_fill_gradient2(
    low = "rosybrown1",
    mid = "mintcream",
    high = "lightslateblue",
    space = "Lab",
    na.value = "mediumpurple1",
    guide = "none",
    aesthetics = "fill"
  ) +
  geom_text(aes(label = round(`excess variance`, 2)),angle = 45,size=2.5) +
  labs(x="mating distance",y="competition distance") +
  theme_minimal() -> evar_plot

evar_plot %>%
  ggsave(
    filename = paste0("figs/", as.character(Sys.Date()), "evar_plot.pdf"),
    device = "pdf",
    height = 2,
    width = 7
  )

# plot of dispersal distance curves
all_connections_l %>%
  dplyr::filter(`mating distance`==min(mat_dists),`competition distance`==min(comp_dists),simplified=="unsimplified") %>%
  ggplot(aes(x = distance, col=as.factor(sigma)), size = 1) +
  facet_wrap(~`dispersal function`,nrow=1) +
  geom_line(stat = "density",bw=0.1) +
  lims(x = c(0, 8)) +
  theme_minimal() +
  theme(strip.text.x = element_text(size = 10)) +
  scale_color_met_d("Hokusai2")+
  labs(
    y = "density",
    x = "distance",
    col=TeX("$sigma$")
  ) -> disp_dist_curves_plot

disp_dist_curves_plot %>%
  ggsave(
    filename = paste0("figs/", as.character(Sys.Date()), "sigma_plot.pdf"),
    device = "pdf",
    height = 3,
    width = 7
  )
