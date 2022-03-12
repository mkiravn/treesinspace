devtools::load_all()
library(ape)
library(tidyverse)
library(ggrepel)
library(ggtree)
library(castor)
library(sf)

# let's generate a random tree with 50 tips
n <- 50
tree <- rtree(n)
# and let's have a look at it
plot(tree)
# let's make an Ne which scales the rate of coalescence
Ne <- 500
# edge_calculator takes the tree above with its topology, and sets branch lengths
tree <- edge_calculator(tree,
                        tree_mode = "coalescent", #stretches the edges so that T_n~exp(nC2)/2Ne
                        # if we had tree_mode = "brown", T_n~exp(1)
                        Ne=Ne) # (this is not needed if tree_mode="brown")
plot(tree) # now we see that the branch lengths are different and all the tips are in the present
hist(tree$edge.length) # we can also take a look at a histogram of branch lengths to see what's going on

# Now we take the nodes and move them in space. The distance in x and y from the parent node ~N(0,t)
# where t is the length of the edge separating them.
tree <- move_nodes(tree)
# the now has an attribute "nodes" which is a dataframe containing node locations
# get the nodes
nodes <- tree$nodes
# plot the nodes, coloured by their time
nodes %>% ggplot(aes(x=x,y=y,col=time,label=node_id))+
  geom_point()+
  geom_label_repel()
# some other modes of exploration!
# calculate the convex hull of the nodes at certain time points and its area
nodes.g <- nodes %>% ungroup() %>%
  mutate(timegroup=cut(time,breaks =50))%>%
  group_by(timegroup,
           ) %>%
  summarize(geometry = st_union(location)) %>%
  st_convex_hull() %>%
  mutate(area = st_area(geometry))
# let's plot the hulls!
ggplot(nodes.g,aes(col=timegroup))+geom_sf(alpha=0.5)
# this is a function something which gets the line connecting parent and child as a spatial object
links <- treesim_connect(tree)
# links is a list of two elements. The second is the interesting one for the edges
links[[2]]
# we can use this to visualise the tree edges in space:
ggplot(links[[2]]) + geom_sf(aes(col=parent_time))

### let's put these plots together:
ggplot() + geom_sf(data=links[[2]],aes(col=child_time)) +
  geom_sf(data=nodes.g,alpha=0.5,fill=NA,size=0.1) +
  geom_point(data=nodes,aes(x=x,y=y,col=time,label=node_id))+
  geom_label_repel(data=nodes,aes(x=x,y=y,col=time,label=node_id))

