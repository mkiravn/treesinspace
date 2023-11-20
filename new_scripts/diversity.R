######## script for treesinspace
library(latex2exp)
library(phyloTop)
library(pegas)
library(apTreeshape)
set.seed(210)

##### Assessing the effect of ecological parameters on tree-based statistics
# defining parameters
Ns <- 100 # number of individuals in population
mat_dists <- c(.2, 2, 20, 200)# mating distance
comp_dists <- c(0, .2, 2, 20, 200) # competition distance
disp_dists <- c(1,5,10) # dispersal distance (sigma)
disp_funs <-
  c("brownian", "normal", "cauchy", "exponential", "uniform") # dispersal kernel
reps <-
  c(1:20) # number of simulation runs for each parameter combination
ngens <- 500 # number of generations
fileout <- "diversity"
pars <- set_slendr_pars(Ns,mat_dists,comp_dists,disp_dists,reps,ngens,fileout)
# defining a world
map <- world(
  xrange = c(-25, 25),
  # min-max longitude
  yrange = c(-25, 25),
  # min-max latitude
  landscape = "blank"
)

res <- run_slendr_simulation(pars,map=map,pathout="diversity")

trees <- res$simplified

######### Plotting
global_labeller <- labeller(
  .default = label_value,
  `mating distance` = label_both,
  `competition distance` = label_both,
  sigma = label_both
)


# Define a function to apply to each tree to get the diversity
get_diversity <- function(tree_element) {
  result <- ts_mutate(tree_element, mutation_rate = 1e-1) %>%
    ts_diversity(sample_sets = list("all" = 1:50), mode = "branch")
  return(result)
}

# Use lapply to apply the function to each element of the list
div_list <- lapply(trees, get_diversity)
div_df <- do.call(rbind, div_list) %>% cbind(pars)

# Define a function to apply to each tree to get the sackin index
get_sackin <- function(tree){
  tr <- tree %>% ts_recapitate(recombination_rate = 0,Ne = 50) %>% ts_phylo(1)
  class(tr) <- "phylo"
  return(sackin.phylo(tr,normalise = T))
}

# Apply to trees
sackin_list <- lapply(trees, get_sackin)
sackin_df <- do.call(rbind, sackin_list) %>% cbind(pars)

# Define a function to apply to each tree to get Tajima's D
get_watt <- function(tree_element) {
  result <- ts_segregating(tree_element,list(p=c(1:50)),mode="branch")
  return(result)
}

# Apply to list
watt_list <- lapply(trees, get_watt)
watt_df <- do.call(rbind, watt_list) %>% cbind(pars)



# Combine data
combined_df <-
  sackin_df %>%
  full_join(div_df, by = colnames(pars)) %>%
  full_join(watt_df, by = colnames(pars)) %>% # modify from here
  rename(
    "mating distance" = "mat_dist",
    "competition distance" = "comp_dist",
    "dispersal function" = "disp_fun",
    "Sackin index" = ".",
    "num. segregating sites" = segsites
  ) %>%
  mutate(`dispersal function` = str_replace(`dispersal function`, "normal", "half-normal"))

combined_df %>% write_delim("tree_statistics.tsv",delim="\t")

combined_df_longer <- combined_df %>%
  pivot_longer(cols=c(`Sackin index`, diversity, `num. segregating sites`),
               names_to = "statistic",
               values_to = "value")

# visualising results
# boxplot
combined_df_longer %>%
  dplyr::filter(`competition distance`==1,`mating distance`==1) %>%
  ggboxplot(x = "sigma",
            y = "value") +
    stat_compare_means(method = "t.test",
                       comparisons = list(c("1","5"),c("5","10"),c("1","10")),
                       label = "p.signif",
                       hide.ns = T,
                       vjust = 1.5) +
    facet_grid(cols=vars(`dispersal function`),rows=vars(statistic),scales = "free_y") +
    theme_bw() +
  labs(x=TeX("$sigma$"),y="") -> sig

combined_df_longer %>%
  dplyr::filter(`competition distance`==1,`sigma`==1) %>%
  ggboxplot(x = "mating distance",
            y = "value") +
  stat_compare_means(method = "t.test",
                     comparisons = list(c("1","5"),c("5","10"),c("1","10"),c("1","100"),c("5","100"),c("10","100")),
                     label = "p.signif",
                     hide.ns = T,
                     vjust = 1.5
                     ) +
  facet_grid(cols=vars(`dispersal function`),rows=vars(statistic),scales = "free_y") +
  theme_bw() +
  labs(y="") -> md

combined_df_longer %>%
  dplyr::filter(`mating distance`==1,`sigma`==1) %>%
  ggboxplot(x = "competition distance",
            y = "value") +
  stat_compare_means(method = "t.test",
                     comparisons = list(c("1","5"),c("5","10"),c("1","10"),c("1","100"),c("5","100"),c("10","100")),
                     label = "p.signif",
                     hide.ns = T,
                     vjust = 1.5
  ) +
  facet_grid(cols=vars(`dispersal function`),rows=vars(statistic),scales = "free_y") +
  theme_bw() +
  labs(y="") -> cd

ggarrange(sig,md,cd,ncol=1) %>%
  ggsave(
    filename = paste0("figs/",
                      as.character(Sys.Date()),
                      "statistics_plot.pdf"),
    device = "pdf",
    height = 9,
    width = 7
  )

combined_df_longer <- combined_df_longer %>% dplyr::filter(statistic != "Tajima's D")

combined_df_longer %>%
  dplyr::filter(`competition distance`==1,`mating distance`==1,
                `dispersal function` == "brownian") %>%
  ggboxplot(x = "sigma",
            y = "value") +
  stat_compare_means(method = "t.test",
                     comparisons = list(c("1","5"),c("5","10"),c("1","10")),
                     label = "p.signif",
                     hide.ns = T,
                     vjust = 1.5) +
  facet_grid(rows=vars(statistic),scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(x=TeX("$sigma$"),y="") -> sigs


combined_df_longer %>%
  dplyr::filter(`competition distance`==1,`sigma`==1,
                `dispersal function` == "brownian") %>%
  mutate(`mating distance` = case_when(`mating distance`==100 ~ "random",
                                            .default = as.character(`mating distance`))) %>%
  ggboxplot(x = "mating distance",
            y = "value") +
  stat_compare_means(method = "t.test",
                     comparisons = list(c("1","5"),c("5","10"),c("1","10"),c("1","random"),c("5","random"),c("10","random")),
                     label = "p.signif",
                     hide.ns = T,
                     vjust = 1.5
  ) +
  facet_grid(rows=vars(statistic),scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(size=10)) +
  labs(y="") -> mds

combined_df_longer %>%
  dplyr::filter(`mating distance`==1,`sigma`==1,
                `dispersal function` == "brownian") %>%
  mutate(`competition distance` = case_when(`competition distance`==100 ~ "no comp.",
                                            .default = as.character(`competition distance`))) %>%
  ggboxplot(x = "competition distance",
            y = "value") +
  stat_compare_means(method = "t.test",
                     comparisons = list(c("1","5"),c("5","10"),c("1","10"),c("1","no comp."),c("5","no comp."),c("10","no comp.")),
                     label = "p.signif",
                     hide.ns = T,
                     vjust = 1.5
  ) +
  facet_grid(rows=vars(statistic),scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(size=10)) +
  labs(y="") -> cds

ggarrange(sigs,mds,cds,nrow = 1) %>%
  ggsave(
    filename = paste0("figs/",
                      as.character(Sys.Date()),
                      "statistics_plot_small.pdf"),
    device = "pdf",
    height = 7,
    width = 8
  )

