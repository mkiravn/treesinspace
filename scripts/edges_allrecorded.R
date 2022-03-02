devtools::load_all("~/slendr")
setup_env()
check_env()
library(ape)
library(ggtree)
library(cowplot)
library(tidyverse)
library(sf)

disp_fun <- "brownian" # trying a few distributions
map <- world(xrange = c(1, 100), yrange = c(1, 100), landscape = "blank")
pop <- population("POP",
                  time = 1,
                  N = 500,
                  dispersal_fun = disp_fun,
                  competition_dist = 10, mate_dist = 1, dispersal_dist = 1,
                  center = c(50, 50),
                  radius = 50,
                  map = map)

# the simulation runs quite long to force the coalescence of the tree
model <- compile(
  populations = pop,
  generation_time = 1,
  sim_length = 100,
  resolution = 0.5
)

# let's get a table of all sampled individuals (because we didn't specify an
# explicit sampling schedule, this table has 100 individuals -- everyone living
# at the end of the simulation)
samples <- sampling(model, times = 100, list(pop, 20))

# note we don't have to specify the sampling schedule explicitly -- in that
# case, slendr will automatically sample everyone living at the very end of
# the simulation (later we will simplify to only a couple of individuals
# ourselves)
slim(
  model,
  sequence_length = 1, recombination_rate = 0, # one tree locus
  method = "batch",
  random_seed = 42,
  retainall=TRUE,
  sampling = samples
)

# load the tree sequence and simplify it (this shouldn't give a recapitation
# warning because the locus should be fully coalesced)
ts <- ts_load(model) %>%
  ts_recapitate(Ne = 500,recombination_rate = 0) %>%
  ts_simplify()

# now let's play with all the remembered nodes
tsu <- ts_load(model)
nds <- ts_load(model) %>% ts_nodes()
edg <- ts_load(model) %>% ts_edges()
dat <- ts_load(model) %>% ts_data()

ggplot(dat)+geom_sf(aes(col=time))

dat <- dat %>%
  rowwise() %>%
  mutate(x=unlist(location)[1],
         y=unlist(location)[2]) %>%
  ungroup()

ggplot(dat)+geom_density(aes(x=y,col=as.factor(time)))

parent_nodes <- dat %>%
  dplyr::as_tibble() %>%
  dplyr::filter(ind_id %in% edg$parent) %>%
  dplyr::select(parent_pop = pop,
                parent_phylo_id = ind_id, parent_node_id = node_id,
                parent_time = time, parent_location = location) %>%
  dplyr::left_join(edges, by = c("parent_phylo_id" = "parent")) %>%
  dplyr::arrange(parent_phylo_id)

# take the `parent_nodes` able above and do another join operation, this time
# with the table of child nodes' times/locations
branch_nodes <- dat %>%
  dplyr::as_tibble() %>%
  dplyr::filter(ind_id %in% edges$child) %>%
  dplyr::select(child_pop  = pop,
                child_phylo_id = ind_id, child_node_id = node_id,
                child_time = time, child_location = location) %>%
  dplyr::inner_join(parent_nodes, by = c("child_phylo_id" = "child")) %>%
  dplyr::arrange(child_phylo_id)

connections <- branch_nodes %>%
  rowwise() %>%
  mutate(pts=st_union(child_location,parent_location)) %>%
  mutate(connection=st_cast(pts,"LINESTRING")) %>% ungroup()

branches <- connections %>% sf::st_set_geometry("connection")
branches %>% ggplot()+geom_sf()
connections <- connections %>% rowwise() %>%
  mutate(edge_gens=child_time-parent_time,
         # calculate the x and y displacements
         disp=st_length(connection),
         x_disp=unlist(child_location)[1]-unlist(parent_location)[1],
         y_disp=unlist(child_location)[2]-unlist(parent_location)[2],
         # and scale everything by generations
         disp_pergen=disp/edge_gens,
         x_disp_pergen=x_disp/edge_gens,
         y_disp_pergen=y_disp/edge_gens) %>%
  ungroup() %>%
  sf::st_set_geometry("child_location")


### diagnostic plots
connections %>% ggplot()+geom_sf(aes(col=child_time)) # this looks good
connections %>% ggplot(aes(x=child_time,y=edge_gens))+geom_point() # this looks really bad
connections %>% ggplot(aes(x=child_time,y=parent_time))+geom_point() # this looks really bad
connections %>% ggplot(aes(x=y_disp_pergen,col=as.factor(child_time))) + geom_density() # this looks really bad
ggplot(connections)+geom_sf(aes(col=edge_gens)) # this looks really bad

### plotting the paths of a lineage
anc <- ts_ancestors(tsu,c(1)) %>%
  select(child_phylo_id=child_id,parent_time,child_time)

anc <- anc %>%
  rowwise() %>%
  # some displacement calculations again
  mutate(xd=unlist(connection)[2]-unlist(connection)[1],
         yd=unlist(connection)[4]-unlist(connection)[3]) %>%
  ungroup()
# let's see it walk in space...
ggplot(anc)+geom_sf(aes(col=child_time)) # this looks good

# and let's see how the dispersals are going!
brown <- data.frame(norm=rnorm(1000,0,1))
anc_longer <- anc %>% pivot_longer(cols=c("xd","yd"),
                     names_to = "coord",
                     values_to = "disp")
ggplot(anc_longer)+geom_density(aes(x=disp,col=coord))+
  geom_density(data=brown,aes(x=norm,col="Expected normal"),lty=2) # this looks great!


