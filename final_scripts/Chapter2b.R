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
devtools::load_all("~/slendr")
setup_env()
check_env()
set.seed(206)

######### Here we define parameters

# defining parameters
Ns <- 2000 # number of individuals in population
mat_dists <- c(0.5,1,2,5,10) # mating distance
comp_dists <- 0.2 # competition distance
disp_dists <- 1 # dispersal distance (sigma)
disp_funs <-
  "brownian" # dispersal kernel
reps <-
  1 # number of simulation runs for each parameter combination
ngens <- 50 # number of generations

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
write_delim(x = pars,file =paste0(as.character(Sys.Date()), "-Chapter2b.tsv"),delim="\t" )

######### Here we run the simulations


# these will hold our data
trees <- c()
trees_us <- c()


# defining a world
map <- world(
  xrange = c(-25, 25),
  # min-max longitude
  yrange = c(-25, 25),
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
    path = "~/Desktop/Chapter2b",
    overwrite = TRUE
  )

  # sampling strategy
  samples <- sampling(model, times = ngens, list(pop, pars[row, "N"]/2))

  # simulate
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
    ts_simplify()

  # extract unsimplified tree
  tsu <- ts_load(model, simplify = FALSE)

  # collect trees
  trees <- c(trees, ts)
  trees_us <- c(trees_us, tsu)

  print(paste("Simulation number:", i))
  i <- i + 1
}

######### Here we convert our data to something more usable

# getting out ape trees
trees_phylo <-
  lapply(trees, ts_phylo, i = 1, quiet = TRUE) # convert
class(trees_phylo) <-  "multiPhylo" # convert list

# merge all the spatial tree data into a dataframe
# first for simplified trees
tree_data <- data.frame()
for (row in c(1:dim(pars)[1])) {
  tree <- trees_phylo[[row]]
  tree_data_i <-
    tree %>% ts_data() %>% rowwise() %>% mutate(
      N = pars[row, "N"],
      fun = pars[row, "disp_fun"],
      rep = as.factor(pars[row, "rep"]),
      mat_dist =
        pars[row, "mat_dist"],
      comp_dist =
        pars[row, "comp_dist"],
      disp_dist =
        pars[row, "disp_dist"],
      sim = as.factor(row),
      x = unlist(location)[1],
      y = unlist(location)[2],
      simplified =
        "simplified"
    ) %>% ungroup()

  tree_data <- rbind(tree_data, tree_data_i)
  print(paste("Converted tree no", row))
}
# now for unsimplified trees
for (row in c(1:dim(pars)[1])) {
  tree <- trees_us[[row]]
  tree_data_i <-
    tree %>% ts_data() %>% rowwise() %>% mutate(
      N = pars[row, "N"],
      fun = pars[row, "disp_fun"],
      rep = as.factor(pars[row, "rep"]),
      mat_dist =
        pars[row, "mat_dist"],
      comp_dist =
        pars[row, "comp_dist"],
      disp_dist =
        pars[row, "disp_dist"],
      sim = as.factor(row),
      x = unlist(location)[1],
      y = unlist(location)[2],
      simplified =
        "unsimplified",
      phylo_id =
        NA
    ) %>% ungroup()

  tree_data <- rbind(tree_data, tree_data_i)
  print(paste("Converted unsimplified tree no", row))
}

# save the tree data
tree_data <- tree_data %>%
  group_by(sim) %>%
  # adjust so that time starts at 0 (for recapitated trees)
  mutate(timeoriginal = time,
         time = timeoriginal - min(timeoriginal))

# the connections (spatial lines between nodes)
conns_simp <- lapply(trees, ts_connect)
conns <- lapply(trees_us, ts_connect)

all_connections <- data.frame() # this will hold the connections

for (row in c(1:dim(pars)[1])) {
  print(paste("Connecting tree number:", row))
  all_connections_i <- conns[[row]][[1]] %>%
    as.data.frame() %>%
    mutate(simplified = "unsimplified") %>%
    cbind(pars[row, ])
  all_connections_i <- conns_simp[[row]][[1]] %>%
    as.data.frame() %>%
    mutate(simplified = "simplified") %>%
    cbind(pars[row, ]) %>% rbind(all_connections_i)
  all_connections <- rbind(all_connections, all_connections_i)
}

all_connections <- all_connections %>% mutate(angle=atan(x_dist/y_dist))
all_connections$mat_dist <- factor(all_connections$mat_dist, levels = mat_dists)

# reading in theoretical results
probsm <- read_delim("ThirdSidePDF.tsv.dat",col_names = seq(1,10,0.2))
colnames(probsm) <- seq(1,10,0.2)
probsm <- probsm %>% as.data.frame() %>% cbind(data.frame(dist=seq(0,20,0.5)))
probsm <- probsm  %>% pivot_longer(cols = -dist,names_to = "mat_dist",values_to = "p")
ggplot(probsm) +
  geom_line(aes(
    x = dist,
    y = p,
    group = mat_dist,
    col = as.numeric(mat_dist)
    ),alpha=0.8) + scale_color_viridis_c() +
  theme_minimal() +
  stat_function(col="black",fun = drayleigh, args = list(scale=1),alpha=0.8,lty=2) +
  labs(col="Mating distance r_b",y="g_y(y)",x="y") -> monster.plot


# probability if both sides are Rayleigh
ThirdSideRR <- function(y,sig,rb){
  0.5*drayleigh(y,scale=sig) + 0.5*drayleigh(y,scale=sqrt(sig^2+rb^2))
}
probsm <- probsm %>% rowwise() %>% mutate(prr=ThirdSideRR(dist,1,as.numeric(mat_dist)))

ggplot(probsm) + geom_line(aes(
  x = dist,
  y = prr,
  group = mat_dist,
  col = as.numeric(mat_dist)
  ),alpha=0.8) + scale_color_viridis_c() +
  theme_minimal() +
  stat_function(col="black",fun = drayleigh, args = list(scale=1),alpha=0.8,lty=2) +
  labs(col="Standard deviation s_b",y="h_y(y)",x="y") -> rr.plot

ggarrange(monster.plot,rr.plot,ncol=1,labels="auto") %>%
  ggsave(
    filename = paste0("figs/", as.character(Sys.Date()), "-theoretical_plots.pdf"),
    device = "pdf",
    height = 7,
    width = 10
  )

probsm <- probsm %>% filter(mat_dist %in% mat_dists)
probsm$mat_dist <- factor(probsm$mat_dist, levels = mat_dists)
probsm <- probsm %>% group_by(mat_dist) %>% mutate(cp=cumsum(p))

######### Plotting


# Checking that dispersal is normal in x and y and the norm is rayleigh


all_connections <- all_connections %>% mutate(mat_dist=as.factor(mat_dist))
# main plots
all_connections %>% filter(simplified=="unsimplified") %>%
  ggplot(aes(x=dist,col=mat_dist)) +
  geom_density() +
  stat_function(col="black",fun = drayleigh, args = list(scale=1),alpha=0.8,lty=2) +
  theme_minimal() +
  labs(y="density") -> brownian_dispersal_simulated
probsm %>%
  ggplot(aes(x=dist,col=mat_dist,y=p)) +
  geom_line() +
  stat_function(col="black",fun = drayleigh, args = list(scale=1),alpha=0.8,lty=2) +
  theme_minimal() +
  labs(y="density") -> brownian_dispersal_theoretical
ggarrange(brownian_dispersal_simulated,
          brownian_dispersal_theoretical,
          ncol=2,common.legend = T) -> brownian_dispersal


# adding qqplots
all_connections %>% filter(simplified=="unsimplified") %>%
  ggplot() +
  stat_qq(aes(sample=dist,col=mat_dist),distribution = qrayleigh,alpha=1,geom="path")  +
  geom_abline(slope=1,intercept=0,lty=2,col="black")  +
  labs(y="simulated",x="theoretical") +
  theme_minimal() -> brownian_qqplots

# comparison to theoretical distribution
all_connections %>% filter(simplified=="unsimplified") %>%
  ggplot(aes(x=dist,col=mat_dist)) +
  geom_line(data=probsm,aes(x=dist,y=p,col=mat_dist,lty="theoretical"))+
  #geom_line(data=probsm,aes(x=dist,y=prr,col=mat_dist,lty="h_y(y)"))+
  geom_density(aes(lty="simulated")) +
  stat_function(col="black",fun = drayleigh, args = list(scale=1),alpha=0.8,lty=2) +
  theme_minimal() +
  facet_wrap(~mat_dist,labeller=label_both,nrow = 1) +
  labs(y="density")  -> sim_theo

ggarrange(sim_theo,brownian_dispersal,brownian_qqplots,ncol=1,labels="auto",common.legend = T) %>%
  ggsave(
    filename = paste0("figs/", as.character(Sys.Date()), "-mating_distance.pdf"),
    device = "pdf",
    height = 7,
    width = 10
  )
