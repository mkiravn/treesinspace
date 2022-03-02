### Simulating some trees to look at the spread of nodes over time

# setting up slendr:
devtools::load_all("~/slendr")
setup_env()
check_env()

# defining a region
map <- world(
  xrange = c(1, 100), # min-max longitude
  yrange = c(1, 100), # min-max latitude
  landscape = "blank"
)

# defining parameters
Ns <- 300
mat_dists <- 3 # the mate distance- let's keep this fixed here
comp_dists <- 3
disp_dists <- c(1:3) # mean dispersal distance
disp_funs <- c("normal","brownian") # trying a few distributions
reps <- c(1:3)
ngens <- 1000

# let's also make a parameter table
pars <- expand.grid(N=Ns,
                    mat_dist=mat_dists,
                    comp_dist=comp_dists,
                    disp_fun=disp_funs,
                    disp_dist=disp_dists,
                    rep=reps,
                    ngen=ngens)

# running the model
trees <- c()
i <- 1

for (row in c(1:dim(pars)[1])){

  # define a population
  pop <- population("POP",
                    time = 1,
                    N = pars[row,"N"],
                    center = c(50, 50),
                    radius = 50,
                    map = map,
                    mate_dist = pars[row,"mat_dist"],
                    competition_dist = pars[row,"comp_dist"],
                    dispersal_fun = as.character(pars[row,"disp_fun"]),
                    dispersal_dist = pars[row,"disp_dist"])

  # compile and run the model
  model <- compile(
    populations = pop,
    generation_time = 1,
    sim_length = ngens, # a forward in time model in slendr needs length of the sim
    resolution = 1, # resolution in "distance units per pixel"
    # how can we control this distribution?
    path = "~/Desktop/test-model",
    overwrite = TRUE
  )

  # let's sample 10 individuals in the "present"
  samples <- sampling(model, times = ngens, list(pop, 100))

  # simulate
  slim(
    model, sampling = samples,
    sequence_length = 1,
    recombination_rate = 0, # simulate only a single locus
    method = "batch", # change to "gui" to execute the model in SLiMgui
    #random_seed = 314159,
    verbose = FALSE,
    burnin = 0
  )

  # extract trees
  ts <- ts_load(model, simplify = TRUE) %>%
    ts_recapitate(Ne = pars[row,"N"],
                  recombination_rate = 0) %>%
    ts_simplify() %>%
    ts_mutate(mutation_rate = 1e-4)

  # collect trees
  trees <- c(trees,ts)

  print(paste("Simulation number:",i))
  i <- i+1
}

### Exploring the trees ------------------------
library(treespace)
library(adegraphics)
library(adagenet)
library(pheatmap)

# getting out the tree objects
trees_phylo <- lapply(trees, ts_phylo, i=1) # convert
class(trees_phylo) <-  "multiPhylo" # convert list

# Spatial plots
library(sf)
library(tidyverse)
library(MetBrewer)
library(ggpubr)

tree_data <- data.frame()

# merge all the spatial tree data into a dataframe to plot
for (row in c(1:dim(pars)[1])){

  tree <- trees_phylo[[row]]

  tree_data_i <- tree %>% ts_data() %>% mutate(N=pars[row,"N"],
                                               fun=pars[row,"disp_fun"],
                                               rep=as.factor(pars[row,"rep"]),
                                               mat_dist=pars[row,"mat_dist"],
                                               comp_dist=pars[row,"comp_dist"],
                                               disp_dist=pars[row,"disp_dist"])

  tree_data <- rbind(tree_data,tree_data_i)
  print(paste("Converting tree no",row))
  #print(paste(colnames(pars),pars[row,]))
}

bks <- c(-Inf,seq(ngens-100,ngens,10))
tree_data <- tree_data %>%
  mutate(timegroup=cut(time, breaks=bks))

tp <- tree_data %>%
  group_by(timegroup,
           N,
           fun,
           rep=as.factor(rep),
           mat_dist,
           comp_dist,
           disp_dist) %>%
  summarize(geometry = st_union(location)) %>%
  st_convex_hull() %>%
  mutate(area = st_area(geometry))

# let's make lots of summary plots!
d_dist <- disp_dists[1]
m_dist <- mat_dists[1]
c_dist <- comp_dists[1]
n <- Ns[1]
f <- disp_funs[1]

pt <- ggplot(data=filter(tree_data,N==n,disp_dist==d_dist,comp_dist==c_dist),aes(col=rep,fill=rep,lty=fun)) +
  lims(x=c(1,100),y=c(1,100)) +
  facet_wrap(~timegroup) +
  ggtitle("Node Distribution") +
  theme(axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank())+
  theme_minimal()

pt <- pt + geom_sf(data=filter(tp,N==n,fun==f,disp_dist==d_dist,comp_dist==c_dist),
                   alpha=0.1) +
  scale_fill_manual(values=met.brewer("Hokusai1",n = length(reps)))+
  scale_colour_manual(values=met.brewer("Hokusai1",n=length(reps)))+
  geom_sf(size=0.1)

pt2 <- ggplot(data=filter(tree_data,N==n,disp_dist==d_dist,comp_dist==c_dist),aes(col=timegroup,fill=timegroup)) +
  lims(x=c(1,100),y=c(1,100)) +
  facet_wrap(~fun+rep) +
  ggtitle("Node Distribution") +
  theme(axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank())+
  theme_minimal()

pt2 <- pt2 + geom_sf(data=filter(tp,N==n,disp_dist==d_dist,comp_dist==c_dist,rep==max(reps)),
                   alpha=0.1)

pt2 <- pt2+geom_sf(size=0.1)
pt2 <- pt2 +
  scale_fill_manual(values=met.brewer("Hokusai2",n = length(bks)))+
  scale_colour_manual(values=met.brewer("Hokusai2",n=length(bks)))


tp_f <- tp %>% filter(N==n,disp_dist==d_dist,comp_dist==c_dist,is.na(timegroup)==FALSE,!grepl("Inf",timegroup))
mdd <- ggplot(tp_f,aes(x=timegroup,col=fun,y=sqrt(area),group=interaction(timegroup,fun))) +
  geom_smooth(aes(x=timegroup,y=sqrt(area),group=fun),method="loess",se=FALSE,alpha=0.1)+
  geom_smooth(aes(x=timegroup,y=sqrt(area),group=interaction(fun,rep),col=fun),method="loess",se=FALSE,alpha=0.1,size=0.1)+
  geom_point(size=0.5) +
  theme_classic() +
  scale_colour_manual(values=met.brewer("Hokusai1",n=length(reps)))+
  ggtitle("Bounding Area- time")

ggarrange(pt2,mdd,ncol=1)
ggsave(filename="figs/clouds-spread.pdf",device="pdf",width=10,height=20)

tp_f <- tp %>% filter(N==n,disp_dist==disp_dists[2],comp_dist==comp_dists[1],is.na(timegroup)==FALSE,!grepl("Inf",timegroup))

tp_f <- tp_f %>% mutate(timen=as.numeric(timegroup)*10)
l <- lm(sqrt(area) ~timen,data=tp_f)
attributes(l)

