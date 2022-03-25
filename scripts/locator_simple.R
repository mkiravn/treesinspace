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


tr <- rtree(50)
tr <- edge_calculator(tr, tree_mode="coalescent", Ne=100)
plot(tr)

bias <- 1

tr <- move_nodes(tr,disp_dist = 1, bias_x=bias,bias_y=0)

tr <- find_ancestor_rec(tr)
s <- find_ancestors(tr,tr$nodes)

locs <- s %>%
  select(node_id=MRCA,
         simple_loc=multiple) %>%
  distinct(node_id,simple_loc) %>%
  full_join(tr$nodes,by="node_id") %>%
  rowwise() %>%
  mutate(mid_line=st_cast(st_union(inf_loc,location),"LINESTRING"),
         mid_error=st_length(mid_line)) %>%
  mutate(simple_line=st_cast(st_union(simple_loc,location),"LINESTRING"),
         simple_error=st_length(simple_line)) %>% ungroup()

# evaluation
p1 <- ggplot(filter(locs,node_id>tr$Nnode+1),aes(col=time)) +
  geom_sf(aes(geometry=location,shape="true")) +
  geom_sf(aes(geometry=simple_loc,shape="simple")) +
  geom_sf(aes(geometry=inf_loc,shape="midpoint")) +
  geom_sf(aes(geometry=simple_line,lty="simple")) +
  geom_sf(aes(geometry=mid_line,lty="midpoint")) +
  scale_colour_gradientn(colours = met.brewer(name = "Hokusai1"))+
  theme_minimal()

p2 <- ggplot(locs) +
  geom_point(aes(x=mid_error,y=simple_error,col=time))+
  scale_colour_gradientn(colours = met.brewer(name = "Hokusai1")) +
  theme_minimal()


ggarrange(p1,p2,common.legend = TRUE)
