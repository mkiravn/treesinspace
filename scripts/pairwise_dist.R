devtools::load_all()
devtools::load_all("~/slendr")
library(ape)
library(tidyverse)
library(sf)


tr <- rtree(200)
tr <- edge_calculator(tr, tree_mode="coalescent", Ne=10)
tr <- move_nodes(tr,disp_dist = 1)
nodes <- tr$nodes
plot(tr)

pairs <- pairwise_distance(tr,nodes)

ggplot(pairs,aes(x=timedist,y=geodist)) +
  geom_point(aes(col=as.factor(n1))) +
  geom_smooth() +
  theme(legend.position="none")
