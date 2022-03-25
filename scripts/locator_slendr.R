### Locating ancestors by recursive mid-point estiamation

# The idea of this algorith is to recursively walk up a tree and estimate the loactions of internal nodes.

# - Compute age of every internal node
# - While all nodes not solved:
#   - Pick youngest unsolved node
#   - Get children
#   - Infer position based on weighted* child position

# Weighting should happen by branch length. Ie. if we have
#       l1 ___|____ l2
#         |   n3   |
#         n1       n2
# then the location of n3 will be at l1/(l1+l2) p1 + l2/(l1+l2) p2
# where l1 and l2 are the branch lengths and p1 and p2 are the vectors of coordinates of n1 and n2 respectively.
devtools::load_all("~/treesinspace")
devtools::load_all("~/slendr")
library(ape)
library(tidyverse)
library(sf)
library(MetBrewer)
library(ggpubr)

n <- 50


# try on slendr data

# setting some parameters
disp_f <- "brownian"
comp_d <- 1
mate_d <- 100
disp_d <- 1

# how many gens?
gens <- 1000
# how large a pop?
N <- n

# making a landscape...
map <-
  world(
    xrange = c(1, 100),
    yrange = c(1, 100),
    landscape = "blank"
  )

# defining the population
pop <- population(
  "POP",
  time = 1,
  N = N,
  center = c(50, 50),
  radius = 50,
  map = map,
  mate_dist = mate_d,
  competition_dist = comp_d,
  dispersal_fun = disp_f,
  dispersal_dist = disp_d
)

# make the model
model <- compile(
  populations = pop,
  generation_time = 1,
  sim_length = gens,
  resolution = 0.5
)
# who should we sample?
samples <- sampling(model, times = gens, list(pop, 10))
# run model
slim(
  model,
  sequence_length = 1,
  recombination_rate = 0, # one tree locus
  method = "batch",
  random_seed = 7,
  retainall=TRUE, # a new option
  sampling = samples
)

# the full tree
tsu <- ts_load(model)
# the connections and branches
res <- ts_connect(tsu)
# the simplified tree
ts <- ts_load(model) %>%
  #ts_recapitate(Ne = N, recombination_rate = 0) %>%
  ts_simplify()
# the phylogenetic tree
tree <- ts %>% ts_phylo(i=1)
#tree$tip.label <- c(1:length(tree$tip.label))
class(tree) <- "phylo"
# the data
data <- ts %>% ts_data() %>% as.data.frame()


# this isn't working and I think it's because of node labelling...
find_ancestor_rec_ts <- function(tree,data){
  # initialise some stuff
  nodes <- data %>%
    as.data.frame() %>%
    mutate(node_id=node_id+1) %>%
    arrange(node_id,decreasing=FALSE) # node information

  dists <- dist.nodes(tree) # distance matrix
  # initialise column
  nodes$inf_loc <- NA
  nodes$inf_loc[1:length(tree$tip.label)] <- nodes$location[1:length(tree$tip.label)]

  while (sum(is.na(nodes$inf_loc))!=0){ # termination condition
    unsolved <- is.na(nodes$inf_loc)==TRUE & nodes$node_id>length(tree$tip.label)
    # pick which node to solve
    if (sum(is.na(nodes$inf_loc))>1){ # doesnt work- look, root is not the right thing
      tosolve <- nodes[nodes$time==max(nodes[unsolved,"time"]),"node_id"] # need to add a tiebreak if they are the same age
    }
    else tosolve <- length(tree$tip.label)+1 # in case it is the root
    children <- phangorn::Children(tree,tosolve) # find children
    edges <- c(dists[children[1],tosolve],dists[children[2],tosolve]) # get edge lengths
    weights <- c(edges[1]/sum(edges),edges[2]/sum(edges)) # weightings
    locations <- c(nodes$inf_loc[children[1]],nodes$inf_loc[children[2]]) # get children locations
    inf_locn <- weights[1]*unlist(locations[1])+weights[2]*unlist(locations[2]) # place the parent
    nodes$inf_loc[tosolve] <- list(inf_locn) # add to dataframe
    print(paste(sum(is.na(nodes$inf_loc)),"nodes left."))
    print(tosolve)
  }
  print("done inferring!")
  # ugly dataframe stuff
  nodes <- nodes %>%
    rowwise() %>%
    mutate(inf_x=inf_loc[[1]],inf_y=inf_loc[[2]]) %>%
    st_as_sf(coords=c("inf_x","inf_y")) %>%
    select(-inf_loc) %>%
    rename(inf_loc=geometry) %>%
    st_set_geometry("location") %>% ungroup()
  # return in nodes
  tree$nodes <- nodes
  return(tree)
}

tree <- find_ancestor_rec_ts(tree,data = data)
