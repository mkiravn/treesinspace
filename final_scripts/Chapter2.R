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
library(gridExtra)
library(patchwork)
devtools::load_all("~/slendr")
setup_env()
check_env()
set.seed(206)

######### Here we define parameters

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
write_delim(x = pars,file =paste0(as.character(Sys.Date()), "-Chapter2.tsv"),delim="\t" )



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
    radius = 50, # (bigger than world)
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
    path = "~/Desktop/Chapter2",
    overwrite = TRUE
  )

  # sampling strategy
  samples <- sampling(model, times = ngens, list(pop, pars[row, "N"]/2))

  # simulate
  slim(
    model,
    sampling = samples,
    # simulate only a single locus
    sequence_length = 1,
    recombination_rate = 0,
    method = "batch",
    # change to "gui" to execute the model in SLiMgui
    random_seed = i,
    verbose = FALSE,
    retainall = TRUE,
    burnin = 0,
  )

  # extract simplified trees
  ts <- ts_load(model, simplify = TRUE) %>%
    ts_recapitate(Ne = pars[row, "N"],
                  recombination_rate = 0) %>%
    ts_simplify() %>%
    ts_mutate(mutation_rate = 1e-4)

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


######### Plotting

# main plots
all_connections %>% filter(simplified=="unsimplified") %>%
  ggplot(aes(x=dist,lty="simulated",col="simulated")) +
  geom_density(show.legend = F) +
  stat_function(aes(lty="theoretical"),fun = drayleigh, args = list(scale=1),alpha=0.8,col="black") +
  theme_minimal() +
  labs(y="density",x="distance",lty="") -> absdist
all_connections %>%
  filter(simplified == "unsimplified") %>%
  select(x = x_dist, y =
             y_dist) %>%
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
all_connections %>% filter(simplified=="unsimplified") %>%
  ggplot() +
  stat_qq(aes(sample=dist,col="distance"),distribution = qrayleigh,alpha=0.5)  +
  stat_qq_line(aes(sample=dist),distribution = qrayleigh,alpha=0.8,lty=2,col="black",dparams = list(scale=1))  +
  labs(y="simulated",x="theoretical") +
  theme_minimal() ->absdist_qq
all_connections %>% filter(simplified=="unsimplified") %>%
  ggplot() +
  stat_qq(aes(sample=x_dist,col="x displacement"),distribution = qnorm,alpha=0.5)  +
  stat_qq(aes(sample=y_dist,col="y displacement"),distribution = qnorm,alpha=0.5)  +
  stat_qq_line(aes(sample=y_dist),distribution = qnorm,alpha=0.8,lty=2,col="black",dparams = list(mean=0,sd=1))  +
  labs(y="simulated",x="theoretical") +
  theme_minimal() -> xandy_qq
ggarrange(absdist,absdist_qq,xandy,xandy_qq,common.legend = T) -> brownian_plots

# and checking that x and y are independent of each other
pcor <- cor.test(all_connections$x_dist, all_connections$y_dist, method="pearson")

all_connections %>% filter(simplified=="unsimplified") %>%
  ggplot(aes(x=x_dist,y=y_dist,col="simulated")) +
  geom_point(alpha=1,size=0.3,show.legend = F) +
  theme_minimal() +
  labs(x="x displacement",y="y displacement")+
  geom_density_2d(alpha=1,show.legend = F,col="grey") +
  stat_cor(col="black",geom="label") +
  coord_fixed()-> brownian_dispersal_angle
ggarrange(brownian_plots,
          brownian_dispersal_angle,
          ncol = 2,
          widths = c(2,1)) -> brownian_plots_panel

# Checking that increments are stationary
all_connections %>%
  mutate(timegroup=cut_width(parent_time,width=25,boundary=0)) %>%
  filter(simplified=="unsimplified") %>%
  ggplot(aes(x=dist,col=timegroup)) +
  geom_density() +
  stat_function(aes(col="Rayleigh"),
                fun = drayleigh,
                args = list(scale=1),
                alpha=0.8,
                lty=2, col="black") +
  theme_minimal() +
  labs(y="density",
       col="time grouping",
       x="distance") -> stationary_increments

all_connections %>%
  mutate(timegroup=cut_width(parent_time,width=25,boundary=0)) %>%
  filter(simplified=="unsimplified") %>%
  ggplot(aes(sample=dist,col=timegroup)) +
  geom_qq(distribution = qrayleigh,dparams=list(scale=1),alpha=0.8,geom="path") +
  stat_qq_line(aes(sample=dist),distribution = qrayleigh,alpha=0.8,lty=2,col="black",dparams = list(scale=1))  +
  theme_minimal() +
  labs(y="simulated",
       x="theoretical",
       col="time grouping") -> stationary_qqplot

all_connections %>% filter(simplified=="unsimplified") %>%
  mutate(timegroup=cut_width(parent_time,width=25,boundary=0)) %>%
  ggplot(aes(x=x_dist,y=y_dist,col=timegroup)) +
  geom_point(alpha=0.5,size=0.1,show.legend = F) +
  theme_minimal() +
  geom_density_2d(show.legend = F) +
  #stat_cor(geom="label",aes(group=timegroup),fill="white",r.accuracy = 0.01,p.accuracy = 0.01,size=4)+
  labs(x="x displacement",y="y displacement",col="time grouping")+
  coord_equal()-> stationary_dispersal_angle
# table of results
all_connections %>% filter(simplified=="unsimplified") %>%
  mutate(timegroup=cut_width(parent_time,width=25,boundary=0)) %>%
  group_by(`time grouping` = timegroup) %>%
  summarise(n=n(),
            `R`= round(cor.test(x_dist,y_dist,method="pearson")$statistic,3),
            `P-value`=round(cor.test(x_dist,y_dist,method="pearson")$p.value,3)) %>%
  ggtexttable(rows=NULL,theme=ttheme("minimal",base_size=8) )-> stationary_table

ggarrange(stationary_table,stationary_dispersal_angle,ncol = 1) -> stationary_rhs

ggarrange(
  stationary_increments,
  stationary_qqplot,
  stationary_table,
  common.legend = T,
  nrow = 1
) -> stationary

# checking for independence
cp_jumps <- data.frame()

while (dim(cp_jumps)[1]<500){
  row <- sample(x=c(1:dim(filter(all_connections,simplified=="unsimplified"))[1]),size=1)
  focal <- filter(all_connections,simplified=="unsimplified")[row,"parent_node_id"]
  if (focal %in% cp_jumps$focal_node){
    print("seen already")
  } else{
    downstream_connection <- filter(all_connections,parent_node_id==focal,simplified=="unsimplified")[1,]
    upstream_connection <- filter(all_connections,child_node_id==focal,simplified=="unsimplified")
    toadd <- data.frame(focal_node=focal,
                        downstream_distance=downstream_connection$dist,
                        downstream_angle=downstream_connection$angle,
                        upstream_distance=upstream_connection$dist,
                        upstream_angle=upstream_connection$angle)
    cp_jumps <- rbind(cp_jumps,toadd) %>% filter(is.na(downstream_distance)==F)
  }
}

cp_jumps %>%
  ggplot(aes(y=downstream_distance,x=upstream_distance)) +
  geom_point(show.legend=F, aes(col=as.factor(focal_node)),size=0.5,alpha=0.5) +
  labs(x="upstream",y="downstream",title="Distance") +
  theme_minimal() +
  stat_cor(geom="label") +

  cp_jumps %>%
  ggplot(aes(y=downstream_angle,x=upstream_angle)) +
  geom_point(show.legend=F, aes(col=as.factor(focal_node)),size=0.5,alpha=0.5) +
  labs(x="upstream",y="downstream",title="Angle") +
  theme_minimal() +
  stat_cor(geom="label") -> independence_plot

ggarrange(
  brownian_plots_panel,
  stationary,
  independence_plot,
  labels = "auto",
  ncol = 1,
  heights = c(2, 1, 1)
) -> panel_plot

panel_plot %>%   ggsave(
  filename = paste0("figs/", as.character(Sys.Date()), "Brownian_plots.pdf"),
  device = "pdf",
  height = 8,
  width = 7
)

write_delim(pars,file=paste0("figs/", as.character(Sys.Date()), "-Ch2-pars.tsv"),
            delim="\t",quote = "none")




