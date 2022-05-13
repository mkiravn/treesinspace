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
library(gt)
library(webshot)
library(scattermore)
library(ggbeeswarm)
library(ggrastr)


devtools::load_all("~/slendr")
setup_env()
check_env()
set.seed(206)

######### Here we define parameters

# defining parameters
n <- 200
Ns <- n
mat_dists <- c(0.2, 0.5, 1, 2, 5) # mate distance
comp_dists <- c(0.2) # competition distance
disp_dists <- c(1) # mean dispersal distance
disp_funs <- c("brownian", "cauchy")
reps <- c(1:5)
ngens <- 1000

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
print(paste("Running" , dim(pars)[1], "simulations."))
write_delim(
  x = pars,
  file = paste0(as.character(Sys.Date()), "-Chapter2.tsv"),
  delim = "\t"
)


######### Here we run the simulations


# these will hold our data
trees <- c()
trees_us <- c()
pwds <- data.frame()


# defining a world
map <- world(
  xrange = c(-50, 50),
  # min-max longitude
  yrange = c(-50, 50),
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
    radius = 25,
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
    path = "~/Desktop/Chapter2c",
    overwrite = TRUE
  )

  # sampling strategy
  samples <-
    sampling(model, times = ngens, list(pop, pars[row, "N"] / 4))

  # simulate
  try_again(
    10,
    slim(
      model,
      sampling = samples,
      # simulate only a single locus
      sequence_length = 1,
      recombination_rate = 0,
      method = "batch",
      # change to "gui" to execute the model in SLiMgui
      random_seed = sample(c(1:100)),
      verbose = FALSE,
      retainall = TRUE,
      burnin = 0,
    )
  )


  # extract simplified trees
  ts <- ts_load(model) %>%
    ts_recapitate(Ne = pars[row, "N"],
                  recombination_rate = 0,
                  random_seed = 206) %>%
    ts_simplify() %>%
    ts_mutate(mutation_rate = 1e-4)

  # extract unsimplified tree
  tsu <- ts_load(model, simplify = FALSE)

  # collect trees
  trees <- c(trees, ts)
  trees_us <- c(trees_us, tsu)

  # get pairwise distances
  pwd <- pairwise_distance(ts_phylo(ts, 1))
  pwds <- pwd %>% cbind(pars[row, ]) %>% rbind(pwds)

  print(paste("Simulation number", i, "done."))
  i <- i + 1
}

######### Here we convert our data to something more usable

# merge all the spatial tree data into a dataframe
# the connections (spatial lines between nodes)
conns_simp <- lapply(trees, ts_connect)
conns <- lapply(trees_us, ts_connect)

all_connections <- data.frame() # this will hold the connections

for (row in c(1:dim(pars)[1])) {
  print(paste("Connecting tree number:", row))
  all_connections_i <- conns[[row]][[1]] %>%
    as.data.frame() %>%
    mutate(simplified = "unsimplified") %>%
    cbind(pars[row,])
  all_connections_i <- conns_simp[[row]][[1]] %>%
    as.data.frame() %>%
    mutate(simplified = "simplified") %>%
    cbind(pars[row,]) %>% rbind(all_connections_i)
  all_connections <- rbind(all_connections, all_connections_i)
}


pwds %>%
  mutate(simplified = "tips only") %>%
  distinct(.keep_all = T, edge_gens, dist) %>%
  select(mat_dist,
         comp_dist,
         disp_fun,
         disp_dist,
         simplified,
         sim,
         rep,
         dist,
         edge_gens) %>%
  rbind(
    all_connections %>%
      select(
        mat_dist,
        comp_dist,
        disp_fun,
        disp_dist,
        simplified,
        sim,
        rep,
        dist,
        edge_gens
      )
  ) -> all_dists

#  Estimates from unsimplified trees, simplified trees and tips only
estimates <- all_dists %>%
  group_by(disp_dist, mat_dist, disp_fun, simplified, rep, comp_dist) %>%
  summarise(n = n(),
            sigma_d = sqrt((1 / (2 * n(
            ))) *
              sum((dist / sqrt(
                sqrt(edge_gens)
              )) ^ 2))) %>% ungroup()

# can we estimate mating distance?
all_dists %>% ggplot(aes(x = dist / sqrt(edge_gens), col = simplified)) +
  geom_density() + scale_x_log10() +
  #geom_vline(data=estimates,aes(xintercept=sigma_d,col=simplified))+
  geom_vline(
    data = estimates,
    aes(xintercept = mat_dist),
    col = "purple",
    lty = 3,
    alpha = 1
  ) +
  #geom_vline(aes(xintercept=1),col="purple",lty=3,alpha=0.5)+
  facet_grid(
    rows = vars(mat_dist),
    cols = vars(disp_fun),
    labeller = label_both,
    scales = "free_y"
  ) +
  #stat_function(fun = drayleigh, args = list(scale=1),alpha=0.8,lty=2,col="black") +
  theme_minimal() +
  labs(y = "density") -> estimates.md.plot_


# pairwise distances
all_dists  %>% filter(mat_dist == 0.2, disp_fun == "brownian") %>%
  ggplot(aes(y = dist, x = edge_gens)) +
  #ggrastr::rasterize(geom_point(alpha = 0.5, size = 0.5,col="grey")) +
  geom_bin2d(alpha = 0.9, size = 0.5) +
  stat_function(
    fun = sqrt,
    alpha = 0.8,
    lty = 2,
    col = "black"
  ) +
  facet_grid(
    rows = vars(mat_dist),
    cols = vars(disp_fun, comp_dist),
    labeller = label_both,
    scales = "free_y"
  ) +
  geom_smooth(se = F, method = "loess") +
  theme_minimal() +
  facet_wrap( ~ simplified, nrow = 1) +
  lims(y=c(0,13)) +
  scale_fill_gradientn(colours = met.brewer(name = "Hokusai2"))+
  scale_x_sqrt() -> pairwise_plot

# how does relatedness decay with distance?
all_dists  %>% filter(mat_dist == 0.2, disp_fun == "brownian") %>%
  ggplot(aes(x = dist, y = edge_gens)) +
  #ggrastr::rasterize(geom_point(alpha = 0.5, size = 0.5,col="grey")) +
  geom_bin2d(alpha = 0.9, size = 0.5) +
  stat_function(
    fun = sqrt,
    alpha = 0.8,
    lty = 2,
    col = "black"
  ) +
  facet_grid(
    rows = vars(mat_dist),
    cols = vars(disp_fun, comp_dist),
    labeller = label_both,
    scales = "free_y"
  ) +
  geom_smooth(se = F, method = "loess") +
  theme_minimal() +
  facet_wrap( ~ simplified, nrow = 1) +
  scale_fill_gradientn(colours = met.brewer(name = "Hokusai2")) -> pairwise_plot2

pairwise_plot2 %>%  ggsave(
  filename = paste0("figs/", as.character(Sys.Date()), "-pairs.pdf"),
  device = "pdf",
  height = 3,
  width = 10
)

#  Estimates from unsimplified trees, simplified trees and tips only
estimates_cut <- all_dists %>%
  group_by(mat_dist, disp_fun, simplified, rep, comp_dist, shortedge = edge_gens <
             100) %>%
  summarise(n = n(),
            sigma_d = sqrt((1 / (2 * n(
            ))) *
              sum((dist / sqrt(
                sqrt(edge_gens)
              )) ^ 2))) %>% ungroup()

all_dists %>%  filter(mat_dist == 0.2,disp_fun=="brownian") %>%
  ggplot(aes(x = dist / (sqrt(edge_gens*pi/2)), col = simplified)) +
  geom_density(aes(lty="uncut")) +
  geom_density(data=all_dists %>%  filter(mat_dist == 0.2,disp_fun=="brownian",edge_gens<100),aes(lty="cut")) +
  geom_vline(data = estimates %>%  filter(mat_dist == 0.2,disp_fun=="brownian"), aes(
    xintercept = sigma_d,
    col = simplified, lty="uncut"
  ),alpha=0.5) +
  geom_vline(data = estimates_cut %>% filter(mat_dist == 0.2,disp_fun=="brownian",shortedge==T), aes(
    xintercept = sigma_d,
    col = simplified, lty="cut"
  ),alpha=0.5) +
  geom_vline(aes(xintercept = 1),
             col = "black",
             lty = 3,
             alpha = 0.5) +
  facet_grid(
    rows = vars(mat_dist),
    cols = vars(disp_fun, comp_dist),
    labeller = label_both,
    scales = "free_y"
  ) +
  stat_function(
    fun = drayleigh,
    args = list(scale = 1),
    alpha = 0.8,
    lty = 3,
    col = "black"
  ) +
  theme_minimal() +
  labs(y = "density") -> estimates.plot.cut

# histograms of edge generations
all_dists %>%  filter(mat_dist == 0.2,
                      disp_fun == "brownian",
                      simplified != "unsimplified") %>%
  ggplot(aes(x = edge_gens)) +
  geom_histogram(alpha = 0.6,fill="#00BFC4") +
  facet_wrap( ~ simplified, scales = "free") +
  geom_vline(xintercept = 100,
             lty = 2,
             col = "black") +
  theme_minimal() + labs(title = "uncut branches") -> h1

all_dists %>%  filter(mat_dist == 0.2,
                      disp_fun == "brownian",
                      simplified != "unsimplified",
                      edge_gens < 100) %>%
  ggplot(aes(x = edge_gens)) +
  geom_histogram(alpha = 0.6,fill="#F8766D") +
  facet_wrap( ~ simplified, scales = "free") +
  geom_vline(xintercept = 100,
             lty = 2,
             col = "black") +
  theme_minimal() + labs(title = "cut branches") -> h2
ggarrange(h1,h2,common.legend = T,ncol=1) -> hist.plot


all_dists %>%  filter(mat_dist == 0.2,disp_fun=="brownian") %>%
  ggplot(aes(y = dist / (sqrt(edge_gens*pi/2)),x=simplified),alpha=0.1,size=0.5) +
  geom_violin(aes(col="uncut"),alpha=0.1,size=0.5) +
  geom_violin(data=all_dists %>%  filter(mat_dist == 0.2,disp_fun=="brownian",edge_gens<100),
                   aes(x=simplified,col="cut"),alpha=0.1,size=0.5) +
  geom_point(data = estimates %>%  filter(mat_dist == 0.2,disp_fun=="brownian"), aes(
    y = sigma_d,
    shape=as.factor(rep),
    x=simplified,col="uncut"
  ),size=2) +
  geom_point(data = estimates_cut %>% filter(mat_dist == 0.2,disp_fun=="brownian",shortedge==T), aes(
    y = sigma_d,
    shape=as.factor(rep),
    x=simplified,col="cut"
  ),size=2) +
  facet_grid(
    rows = vars(mat_dist),
    cols = vars(disp_fun, comp_dist),
    labeller = label_both,
    scales = "free_y"
  ) +
  geom_line(data = estimates %>%  filter(mat_dist == 0.2,disp_fun=="brownian"), aes(
    y = sigma_d,
    x=simplified,group=rep,col="uncut"),lty=2)+
  geom_line(data = estimates_cut %>% filter(mat_dist == 0.2,disp_fun=="brownian",shortedge==T), aes(
    y = sigma_d,
    x=simplified,group=rep,col="cut"),lty=2)+
  theme_minimal() + labs(shape="rep")+
  labs(y = "estimated sigma") -> estimates.plot.cut.2

ggplot(data = filter(estimates),
       aes(col = simplified, x = mat_dist, y = sigma_d)) +
  geom_quasirandom(size = 2,width = 0.1) +
  geom_smooth(se = F, size = 0.5,method="lm") +
  theme_minimal() +
  facet_grid(
    cols = vars(disp_fun, comp_dist),
    labeller = label_both,
    scales = "free_y"
  ) +
  labs(y = "estimated sigma") -> estimates.plot.matdist

estimates.plot.matdist  %>%   ggsave(
  filename = paste0("figs/", as.character(Sys.Date()), "-estimate_plots_md.pdf"),
  device = "pdf",
  height = 3,
  width = 10
)

ggarrange(
  pairwise_plot,
  hist.plot,
  estimates.plot.cut.2,
  labels = "auto",
  ncol = 1)  %>%   ggsave(
  filename = paste0("figs/", as.character(Sys.Date()), "-estimate_plots_cut.pdf"),
  device = "pdf",
  height = 12,
  width = 10
)

### output table of results
all_dists %>%
  filter(disp_dist==1) %>%
  cbind(sigma = 1) %>%
  group_by(mat_dist, disp_fun, simplified, sigma) %>%
  summarise(sigma_d = sqrt((1 / (2 * n(
  ))) *
    sum((dist / sqrt(
      sqrt(edge_gens)
    )) ^ 2))) %>%
  arrange(disp_fun, mat_dist, simplified) %>%
  pivot_wider(names_from = c("simplified"),
              values_from = sigma_d) %>%
  ungroup() %>%
  select(simplified,`tips only`) -> t1

all_dists %>%
  filter(edge_gens<100,disp_dist==1) %>%
  cbind(sigma = 1) %>%
  group_by(mat_dist, disp_fun, simplified, sigma) %>%
  summarise(sigma_d = sqrt((1 / (2 * n(
  ))) *
    sum((dist / sqrt(
      sqrt(edge_gens)
    )) ^ 2))) %>%
  arrange(disp_fun, mat_dist, simplified) %>%
  pivot_wider(names_from = c("simplified"), values_from = sigma_d) %>%
  relocate(unsimplified, .before = simplified) %>%
  relocate(sigma, .before = mat_dist) %>%
  rename("mating distance" = mat_dist,
         "dispersal function" = disp_fun,
         "simplified (cut)"= simplified,
         "tips only (cut)"= `tips only`
         ) %>%
  cbind(t1) %>%
  relocate(simplified, .after = unsimplified) %>%
  relocate(`tips only (cut)`, .after = `simplified (cut)`) %>%
  as_tibble() %>%
  gt() %>%
  tab_spanner(label = "estimated sigma (cut)",
              columns = c(`simplified (cut)`,`tips only (cut)`)) %>%
  tab_spanner(label = "estimated sigma",
              columns = c(unsimplified,`simplified`,`tips only`)) %>%
  tab_spanner(label = "dispersal parameters",
              columns = c(sigma, `mating distance`, `dispersal function`)) %>%
  tab_header(title = "Estimating sigma from trees",
              subtitle = "'Cut' means branches longer than 100 generations excluded") %>%
  data_color(
    columns = c(unsimplified,simplified,`tips only`,
                  `simplified (cut)`, `tips only (cut)`),
    colors = scales::col_numeric(
      palette = c(
        "lightpink",
        "white",
        "lightpink",
        "lightpink1",
        "lightpink2",
        "lightpink3",
        "lightpink4"
      ),
      domain = c(0, 6)
    )
  ) %>% gtsave(paste0("figs/",
                      as.character(Sys.Date()),
                      "results_estimatingsigma_cut.png"))
