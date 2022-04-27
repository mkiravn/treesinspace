devtools::load_all()
library(tidyverse)
library(ape)
library(sf)
library(MetBrewer)
library(ggrepel)
library(ggpubr)
library(ggtree)
library(moments)
library(phangorn)
library(ape)
library(gganimate)
library(transformr)
devtools::load_all("~/slendr")
setup_env()
check_env()
set.seed(206)





### Demonstration


### set up simulation
mating_distance <- c(1)
competition_distance <- c(0.2)
sigma <- 1
dispersal_function <- "brownian"


map <- world(
  xrange = c(-20, 20), # min-max longitude
  yrange = c(-20, 20), # min-max latitude
  landscape = "blank"
)

pop <- population("POP",
                  time = 1,
                  N = 100,
                  center = c(0,0),
                  radius = 20,
                  map = map,
                  mate_dist = mating_distance,
                  competition_dist = competition_distance,
                  dispersal_fun = dispersal_function,
                  dispersal_dist = sigma)


# compile the model
model <- compile(
  populations = pop, # population to include
  generation_time = 1, # generation time
  sim_length = 1000, # length of the simulation
  resolution = 1, # resolution in "distance units per pixel"
  path = "~/Desktop/test-model",
  overwrite = TRUE
)

# who to sample
samples <- sampling(model, times = 1000, list(pop, 5))

# simulate
slim(
  model, sampling = samples,
  sequence_length = 1, # simulate only a single site
  recombination_rate = 0, # which is non-recombining
  method = "batch", # change to "gui" to execute the model in SLiMgui
  random_seed = 7,
  verbose = FALSE,
  retainall = TRUE, # remember all individuals- not just nodes
  burnin = 0,
)

# extract trees
ts <- ts_load(model) %>%
  ts_recapitate(Ne = 100,
                recombination_rate = 0) %>%
  ts_simplify()
# extract unsimplified trees
tsu <- ts_load(model, simplify = FALSE)
# extract connections between individuals
connections <-ts_connect(ts)[[1]]
connections_unsimplified <-ts_connect(tsu)[[1]]
# extract tree
tree <- ts_phylo(ts,i=1)
data <- ts_data(tree)
class(tree) <- "phylo" # convert class to phylo
# extract data of unsimplified tree
data_unsimplified <- ts_data(tsu)

data <- data %>% rowwise() %>% mutate(x=unlist(location)[[1]],y=unlist(location)[[2]]) %>% ungroup()

### some plots!

# only individuals in the present
ggplot(filter(data,remembered==T)) +
  geom_sf() +
  geom_label_repel(aes(x=x,y=y,label=phylo_id)) +
  lims(x=c(-20,20),y=c(-20,20)) + theme_minimal()


# adding a ggtree tree
ggtree(tree,col="grey") +
  geom_point(aes(col=x)) +
  geom_label(aes(label=node,col=x)) +
  theme_minimal() +
  theme(legend.position="none")  +
  scale_color_gradientn(colors=met.brewer("Hokusai1")) +
# plot of simplified and unsimplified data
ggplot() +
  geom_sf(data=data_unsimplified,aes(geometry=location,col=time),alpha=0.2,size=0.5) +
  geom_sf(data=connections_unsimplified,aes(geometry=connection,col=child_time),alpha=0.1,size=0.5) +
  geom_sf(data=connections,aes(geometry=connection,col=child_time))+
  geom_sf(data=data,aes(geometry=location,col=time),size=2) +
  geom_sf_label(data=data,aes(col=time,label=phylo_id)) +
  scale_colour_gradientn(colors = met.brewer(name = "Hokusai1",n = 50)) +
  theme_minimal() +
# another plot of points
ggplot() +
  geom_sf(data=data_unsimplified,aes(geometry=location,col=time)) +
  scale_colour_gradientn(colors = met.brewer(name = "Hokusai1",n = 50)) +
  theme_minimal()


all_connections <- connections %>% mutate(simplified="yes") %>%
  rbind(mutate(connections_unsimplified,simplified="no"))

connections_unsimplified %>% ggplot(aes(x=disp)) +geom_density() +theme_minimal()

ancestors <- find_ancestor_rec_ts_ut(tree,data=data)
ancestors$nodes <- ancestors$nodes %>% mutate(phylo_id=node_id) %>% select(-node_id)

ancestors_centroid <- find_ancestors(tree,data=data)
ancestors_centroid <- ancestors_centroid %>% mutate(phylo_id=MRCA) %>% select(-MRCA)

ggplot() +
  geom_sf_label(data=filter(data,time<1000),aes(col=time,label=phylo_id)) +
  scale_colour_gradientn(colors = met.brewer(name = "Hokusai1",n = 50)) +
  geom_sf(data=filter(ancestors$nodes,time<1000),aes(geometry=inf_loc,col=time),size=2) +
  geom_sf_label(data=filter(ancestors$nodes,time<1000),aes(col=time,label=phylo_id,geometry=inf_loc),fill="lightgrey") +
  theme_minimal() +
ggplot() +
  scale_colour_gradientn(colors = met.brewer(name = "Hokusai1",n = 50)) +
  geom_point(data=filter(ancestors$nodes,time<1000),aes(x=time,col=time,y=var),size=2) +
  geom_label_repel(data=filter(ancestors$nodes,time<1000),aes(x=time,col=time,y=var,label=phylo_id)) +
  theme_minimal()


