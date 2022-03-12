library(ape)
library(tidyverse)

n <- 6
tree <- rtree(n)

plot(tree)
nodelabels()

# Compute the times of all internal nodes
# (I could imagine this being called `node_ranks`, with the couple of lines
# of code at the end of the script actually going at the end of this function,
# which would then return the ranks of those nodes, instead of the times simulated
# by rtree)
node_times <- function(tree) {
  root <- length(tree$tip.label) + 1

  # create an R environment object -- it's basically a list() with some special
  # properties but otherwise the same data structure as list()
  env <- new.env()
  # create a new variable (a vector of NAs) in this environment called times
  # -- this is basically what happens when you do stuff like x <- 42 in your R
  # session (that creates a variable x in a "global environment")
  #
  # for fun, you can note that there's a hidden object .GlobalEnv in your R
  # session (just try typing it into R) which contains all variables you've
  # created in that R session
  #
  # you can inspect the contents of an environment by calling ls(envir =
  # .GlobalEnv), which in this case is the same as running just ls()
  env$times <- rep(NA, tree$Nnode + length(tree$tip.label))

  # set the root time in the vector to 0
  env$times[root] <- 0

  # dive down along the tree, starting at the root -- we're passing three
  # parameters through all layers of the recursion: the tree (always the same
  # object), the currently visited node (changes with each recursive dive), and
  # the environment object, in which we're collecting the times of each ndoe
  recursion(tree, root, env)

  # return the final vector of internal node times
  env$times
}

# At a given node in a given tree, inspect left and right children, calculate
# their times based on the time of the current node
recursion <- function(tree, node, env) {
  # indices of parent-child branches in the tree
  edge <- tree$edge

  # take the children of the current node
  children <- phangorn::Children(tree, node)
  left <- children[1]
  right <- children[2]

  # a "tip" is a node in an R phylo tree which has an index 1...nleaves
  # (so any node which has an ID > nleaves is an internal node)
  nleaves <- length(tree$tip.label)

  # recursively dive down to the left subtree (if left isn't a tip)
  if (left > nleaves) {
    # compute the time of the `left` node as the time of the parent plus
    # the branch length connecting those two
    left_branch <- tree$edge.length[edge[, 1] == node & edge[, 2] == left]
    env$times[left] <- env$times[node] + left_branch
    # call this function for node = left
    recursion(tree, left, env)
  }

  # recursively dive down to the right subtree (if right isn't a tip)
  if (right > nleaves) {
    # compute the time of the `right` node as the time of the parent plus
    # the branch length connecting those two
    right_branch <- tree$edge.length[edge[, 1] == node & edge[, 2] == right]
    env$times[right] <- env$times[node] + right_branch
    # call this function for node = right
    recursion(tree, right, env)
  }
}

# get the internal node times
times <- node_times(tree)
# create a vector of indices, one index for each node (internal nodes and tips)
indices <- seq_along(times)

# order the "ranks" internal nodes based on computed times, starting from the
# root all the way to the last coalescent event (or the first coalescent event
# if we're looking backwards in time from the present)
ranks <- indices[order(times, na.last = TRUE)][1:tree$Nnode]
ranks

# Now we want to assign times to all the nodes.
# To do this, we need to know the ranks. For rank 1 (the root), the time is zero
# For each rank k, we draw a time from an exponential (kC2)-
# this is a waiting time, or the time Tk that the tree spends with k nodes.
# We will sum these cumulatively to get the actual time a node existed
# And later use this to calculate branch lengths.
waiting_times <- rep(0,length(times)) # initialise waiting times
waiting_times[1] <- 0 # root has time 0
for (i in c(1:length(ranks)+1)){ # set the waiting times as kC2
  waiting_times[i] <- rexp(n=1,rate=choose(node_index,2))
}
node_times <- cumsum(waiting_times) # cumulatively sum to get node times
ranks_c <- c(ranks,c(1:n)) # assign the missing ranks of tips (they all have the same time)
node_times_ordered <- node_times[order(ranks_c)] # now in order of node index
# calculate edge lengths
tree_edge <- tree$edge %>%
  as.tibble() %>%
  rename("parent"="V1","child"="V2") %>%
  mutate(edge_length=node_times_ordered[child]-node_times_ordered[parent])
# set edge lengths
tree$edge.length <- tree_edge$edge_length
# let's
plot(tree)
nodelabels()
