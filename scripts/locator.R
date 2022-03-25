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

n <- 100

finding_exploration <- function(n,bias=0){
  # compares mid-point estimation
  tree <- rtree(n)
  tree <- edge_calculator(tree, tree_mode="coalescent", Ne=100)
  tree <- move_nodes(tree ,disp_dist = 1, bias_x = bias,bias_y = 0)
  tree$tip.label <- c(1:n)
  tree <- find_ancestor_rec(tree)
  s <- find_ancestors(tree,tree$nodes,stripped = TRUE)

  locs <- s %>%
    select(node_id=MRCA,
           simple_loc=multiple) %>%
    distinct(node_id,simple_loc) %>%
    full_join(tree$nodes,by="node_id") %>%
    rowwise() %>%
    mutate(mid_line=st_cast(st_union(inf_loc,location),"LINESTRING"),
           mid_error=st_length(mid_line)) %>%
    mutate(simple_line=st_cast(st_union(simple_loc,location),"LINESTRING"),
           simple_error=st_length(simple_line)) %>% ungroup()
  return(locs)
}

ress <- data.frame()

for (bias in c(0,1,10,100)){
  for (rep in c(1:10)) {
    res <- finding_exploration(n=n,bias=bias) %>% mutate(rep=rep,bias=bias)
    ress <- rbind(ress,res)
  }
}

# evaluation
p1 <- ggplot(filter(ress,node_id>n),aes(col=time)) +
  geom_sf(aes(geometry=location,shape="true")) +
  geom_sf(aes(geometry=simple_loc,shape="simple")) +
  geom_sf(aes(geometry=inf_loc,shape="midpoint")) +
  geom_sf(aes(geometry=simple_line,lty="simple")) +
  geom_sf(aes(geometry=mid_line,lty="midpoint")) +
  facet_wrap(~bias,nrow = 1) +
  scale_colour_gradientn(colours = met.brewer(name = "Hokusai1"))+
  theme_minimal()

p2 <- ggplot(filter(ress,node_id>n)) +
  geom_point(aes(x=mid_error,y=simple_error,col=time))+
  facet_wrap(~bias)+
  geom_abline(slope=1,intercept=0,lty=2,col="grey") +
  scale_x_log10() + scale_y_log10()  +
  scale_colour_gradientn(colours = met.brewer(name = "Hokusai1")) +
  theme_minimal()

ggarrange(p1,p2,common.legend = TRUE,ncol = 1) %>%
  ggsave(filename = paste0("figs/",as.character(Sys.Date()),"-locator.pdf"),device="pdf")

