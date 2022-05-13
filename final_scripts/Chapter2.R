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


# Checking that dispersal is normal in x and y and the norm is rayleigh


# ks statistics
p <- all_connections %>% filter(simplified=="unsimplified") %>% select(dist) %>% ks.test("prayleigh", 1)
xp <- all_connections %>% filter(simplified=="unsimplified") %>% select(x_dist) %>% ks.test("pnorm", 1,0)
yp <- all_connections %>% filter(simplified=="unsimplified") %>% select(y_dist) %>% ks.test("pnorm", 1,0)


# main plots
all_connections %>% filter(simplified=="unsimplified") %>%
  ggplot(aes(x=dist,col="simulated")) +
  geom_density() +
  stat_function(aes(col="Rayleigh"),fun = drayleigh, args = list(scale=1),alpha=0.8,lty=2) +
  theme_minimal() +
  labs(y="density",caption = paste("P =",p$p.value))+
all_connections %>% filter(simplified=="unsimplified") %>%
  ggplot(aes(x=x_dist,col="simulated")) +
  geom_density() +
  stat_function(aes(col="Normal"),fun = dnorm, args = list(mean=0,sd=1),alpha=0.8,lty=2) +
  theme_minimal() +
  labs(y="density",caption = paste("P =",xp$p.value)) +
all_connections %>% filter(simplified=="unsimplified") %>%
  ggplot(aes(x=y_dist,col="simulated")) +
  geom_density() +
  stat_function(aes(col="Normal"),fun = dnorm, args = list(mean=0,sd=1),alpha=0.8,lty=2) +
  theme_minimal() +
  labs(y="density",caption = paste("P =",yp$p.value))-> brownian_dispersal

# adding qqplots
all_connections %>% filter(simplified=="unsimplified") %>%
  ggplot() +
  stat_qq(aes(sample=dist,col="dist"),distribution = qrayleigh,alpha=0.5)  +
  stat_qq_line(aes(sample=dist),distribution = qrayleigh,alpha=0.8,lty=2,col="black",dparams = list(scale=1))  +
  labs(y="simulated",x="theoretical") +
  theme_minimal() +
all_connections %>% filter(simplified=="unsimplified") %>%
  ggplot() +
  stat_qq(aes(sample=x_dist,col="x_dist"),distribution = qnorm,alpha=0.5)  +
  stat_qq(aes(sample=y_dist,col="y_dist"),distribution = qnorm,alpha=0.5)  +
  stat_qq_line(aes(sample=y_dist),distribution = qnorm,alpha=0.8,lty=2,col="black",dparams = list(mean=0,sd=1))  +
  labs(y="simulated",x="theoretical") +
  theme_minimal() -> brownian_qqplots


# and checking that x and y are independent of each other

pcor <- cor.test(all_connections$x_dist, all_connections$y_dist, method="pearson")

all_connections %>% filter(simplified=="unsimplified") %>%
  ggplot(aes(x=angle,col="simulated"),fill=NA) +
  geom_histogram(alpha=0.2) +
  geom_vline(xintercept = -pi/2,alpha=0.8,lty=2) +
  geom_vline(xintercept = pi/2,alpha=0.8,lty=2) +
  theme_minimal() +
  labs(x="angle",y="density",title="Angle distribution") +
  labs(y="density")+
all_connections %>% filter(simplified=="unsimplified") %>%
  ggplot(aes(x=x_dist,y=y_dist,col="simulated")) +
  geom_point(alpha=0.5,size=0.5) +
  theme_minimal() +
  labs(x="",y="",title="Individual dispersals",caption=paste("P =",round(pcor$p.value,3)))+
  coord_fixed()-> brownian_dispersal_angle

# Checking that increments are stationary
all_connections %>%
  mutate(timegroup=cut(child_time,4)) %>%
  filter(simplified=="unsimplified") %>%
  ggplot(aes(x=dist,col=timegroup)) +
  geom_density() +
  stat_function(aes(col="Rayleigh"),
                fun = drayleigh,
                args = list(scale=1),
                alpha=0.8,
                lty=2) +
  theme_minimal() +
  labs(y="density") -> stationary_increments
# qqplot
all_connections %>%
  mutate(timegroup=cut(child_time,4)) %>%
  filter(simplified=="unsimplified") %>%
  ggplot(aes(sample=dist,col=timegroup)) +
  geom_qq(distribution = qrayleigh,dparams=list(scale=1),alpha=0.8,geom="path") +
  stat_qq_line(aes(sample=dist),distribution = qrayleigh,alpha=0.8,lty=2,col="black",dparams = list(scale=1))  +
  theme_minimal() +
  labs(y="simulated",x="theoretical") -> stationary_qqplot
all_connections %>% filter(simplified=="unsimplified") %>%
  mutate(timegroup=cut(child_time,4)) %>%
  ggplot(aes(x=angle,fill=timegroup)) +
  geom_histogram(alpha=0.2) +
  geom_vline(xintercept = -pi/2,alpha=0.8,lty=2) +
  geom_vline(xintercept = pi/2,alpha=0.8,lty=2) +
  facet_wrap(~timegroup) +
  theme_minimal() +
  labs(x="angle",y="density",title="Angle distribution") +
  labs(y="density")+
  all_connections %>% filter(simplified=="unsimplified") %>%
  mutate(timegroup=cut(child_time,4)) %>%
  ggplot(aes(x=x_dist,y=y_dist,col=timegroup)) +
  geom_point(alpha=0.5,size=0.5) +
  facet_wrap(~timegroup) +
  theme_minimal() +
  labs(x="",y="",title="Individual dispersals")+
  coord_fixed() -> stationary_dispersal_angle


# Checking that increments are independent
samp <- c()
while (length(unique(samp$node_id))<50){
  print(paste(50-length(unique(samp$node_id)),"left to sample."))
  i <- sample(c(1:Ns),1)
  sampi <-
    ts_ancestors(tsu, i) %>% as.data.frame() %>%
    filter(!(parent_id %in% samp$parent_id))
  if (dim(sampi)[1]>40){samp <- rbind(samp,sampi)}
}


all_connections %>%
  select(dist,parent_id=parent_node_id) %>%
  right_join(samp,by="parent_id")%>%
  ggplot(aes(x = dist, col = as.factor(node_id)),size=0.1) +
  geom_line(aes(color=as.factor(node_id)), stat="density", size=0.5, alpha=0.3, show.legend = FALSE) +
  stat_function(
    aes(col = "Rayleigh"),
    fun = drayleigh,
    args = list(scale = 1),
    alpha = 1,
    lty = 2,
    col="black",
    show.legend = FALSE
  ) +
  theme_minimal()  +
  labs(col = "sample",lty="simulated") +
  labs(y="density") -> independent_increments
# qqplot
all_connections %>%
  select(dist,parent_id=parent_node_id) %>%
  right_join(samp,by="parent_id")%>%
  ggplot(aes(sample = dist, col = as.factor(node_id)),alpha=0.1,size=1) +
  geom_qq(distribution = qrayleigh,dparams=list(scale=1),alpha=0.5, show.legend = FALSE,geom="path") +
  stat_qq_line(aes(sample=dist),distribution = qrayleigh,alpha=0.8,lty=2,col="black",dparams = list(scale=1))  +
  theme_minimal()  +
  labs(col = "sampled tip",lty="simulated") +
  labs(y="simulated",x="theoretical") -> independent_qqplot
all_connections %>% filter(simplified=="unsimplified") %>%
  select(dist,parent_id=parent_node_id,angle) %>%
  right_join(samp,by="parent_id") %>%
  ggplot(aes(sample=angle,col=as.factor(node_id))) +
  geom_qq(distribution = qunif,dparams=list(min=-pi/2,max=pi/2),alpha=0.5, show.legend = FALSE,geom="path") +
  stat_qq_line(aes(sample=angle),distribution = qunif,dparams=list(min=-pi/2,max=pi/2),alpha=0.8,lty=2,col="black") +
  theme_minimal() +
  labs(x="theoretical",y="simulated",title="Angle distribution") +
  all_connections %>% filter(simplified=="unsimplified") %>%
  select(dist,parent_id=parent_node_id,x_dist,y_dist) %>%
  right_join(samp,by="parent_id") %>%
  ggplot(aes(x=x_dist,y=y_dist,col=as.factor(node_id))) +
  geom_point(alpha=0.5,size=0.5,show.legend=FALSE) +
  theme_minimal() +
  labs(x="",y="",title="Individual dispersals")+
  coord_fixed() -> independent_dispersal_angle


# Scaling branch lengths
all_connections %>%
  ggplot(aes(x=dist/sqrt(edge_gens),col=simplified)) +
  geom_density() +
  stat_function(fun = drayleigh,
                args = list(scale=1),
                alpha=0.8,
                lty=2,col="black") +
  theme_minimal() +
  labs(y="density") -> scaling
all_connections %>%
  ggplot(aes(sample=dist/sqrt(edge_gens),col=simplified)) +
  geom_qq(distribution = qrayleigh,dparams=list(scale=1),alpha=0.5) +
  stat_qq_line(aes(sample=dist),distribution = qrayleigh,alpha=0.8,lty=2,col="black",dparams = list(scale=1))  +
  theme_minimal() +
  labs(y="simulated",x="theoretical") -> scaling_qqplot


ggarrange(brownian_dispersal,
          stationary_increments+independent_increments,
          ncol = 1, labels="auto") -> brownian_plots
brownian_plots %>%   ggsave(
  filename = paste0("figs/", as.character(Sys.Date()), "-Brownian_plots.pdf"),
  device = "pdf",
  height = 5,
  width = 10
)
ggarrange(brownian_qqplots,
          stationary_qqplot+independent_qqplot,
          ncol = 1, labels="auto") -> brownian_qplots
brownian_qplots %>%   ggsave(
  filename = paste0("figs/", as.character(Sys.Date()), "-Brownian_qqplots.pdf"),
  device = "pdf",
  height = 5,
  width = 10
)

ggarrange(brownian_dispersal_angle,
          stationary_dispersal_angle,independent_dispersal_angle,
          ncol = 1, labels="auto") %>%  ggsave(
  filename = paste0("figs/", as.character(Sys.Date()), "-Brownian_plots_angle.pdf"),
  device = "pdf",
  height = 8,
  width = 10
)

ggarrange(scaling,scaling_qqplot)  %>%  ggsave(
  filename = paste0("figs/", as.character(Sys.Date()), "-Brownian_scaling.pdf"),
  device = "pdf",
  height = 2,
  width = 10
)

write_delim(pars,file=paste0("figs/", as.character(Sys.Date()), "-Ch2-pars.tsv"),delim="\t",quote = "none")
