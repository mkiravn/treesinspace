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
  # the environment object, in which we're collecting the times of each node
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


# Now we want to assign times to all the nodes.
# To do this, we need to know the ranks. For rank 1 (the root), the time is zero
# For each rank k, we draw a time from an exponential (kC2)-
# this is a waiting time, or the time Tk that the tree spends with k nodes.
# We will sum these cumulatively to get the actual time a node existed
# And later use this to calculate branch lengths.
edge_calculator <- function(tree,tree_mode="coalescent",Ne=length(tree$tip.label)){

  # get the internal node times
  times <- node_times(tree)
  # create a vector of indices, one index for each node (internal nodes and tips)
  indices <- seq_along(times)

  # order the "ranks" internal nodes based on computed times, starting from the
  # root all the way to the last coalescent event (or the first coalescent event
  # if we're looking backwards in time from the present)
  ranks <- indices[order(times, na.last = TRUE)][1:tree$Nnode]
  waiting_times <- rep(0,length(times)) # initialise waiting times
  waiting_times[1] <- 0 # root has time 0
  if (tree_mode=="coalescent"){
    for (i in c(1:length(ranks)+1)){ # set the waiting times as kC2/2Ne
      node_index <- i
      waiting_times[i] <- rexp(n=1,rate=choose(node_index,2)/(2*Ne))
    }
  }
  else if (tree_mode=="brown"){
    for (i in c(1:length(ranks)+1)){ # set the waiting times as kC2
      node_index <- i
      waiting_times[i] <- rexp(n=1,rate=node_index)
    }
  }

  node_times <- cumsum(waiting_times) # cumulatively sum to get node times
  n <- length(tree$tip.label)
  ranks_c <- c(ranks,c(1:n)) # assign the missing ranks of tips (they all have the same time)
  node_times_ordered <- node_times[order(ranks_c)] # now in order of node index, rather than rank

  # calculate edge lengths
  tree_tbl <- tree$edge %>%
    as_tibble() %>%
    rename("parent"="V1","child"="V2") %>%
    mutate(edge_length=node_times_ordered[child]-node_times_ordered[parent])

  # set edge lengths
  tree$edge.length <- tree_tbl$edge_length
  tree$table <- tree_tbl # add a new attribute to the tree with the edge lengths in a table
  tree$node.times <- node_times_ordered
  return(tree)
}

# Now let's simulate some displacement.
# for Brownian displacement, an edge of length t displaces a node by N(0,t)
# from its parent. I will happily cannibalise some of the functions above.

node_space <- function(tree,dim) {
  root <- length(tree$tip.label) + 1
  env <- new.env()
  env$`dim` <- rep(NA, tree$Nnode + length(tree$tip.label))

  # set the root position to 0
  env$`dim`[root] <- 0

  # dive down along the tree, starting at the root -- we're passing three
  # parameters through all layers of the recursion: the tree (always the same
  # object), the currently visited node (changes with each recursive dive), and
  # the environment object, in which we're collecting the times of each node
  recursion_space(tree, root, env,dim)

  # return the final vector of internal node positions
  env$dim
}

# At a given node in a given tree, inspect left and right children, calculate
# their x or y positions based on the position of the current node
recursion_space <- function(tree, node, env , dim) {

  # indices of parent-child branches in the tree
  edge <- tree$table

  # take the children of the current node
  children <- phangorn::Children(tree, node)
  left <- children[1]
  right <- children[2]

  nleaves <- length(tree$tip.label)
  # recursively dive down to the left subtree (if left isn't a tip)
  if (left > nleaves) {
    # compute the time of the `left` node as the position of the parent plus
    # the branch length connecting those two
    left_branch <- edge[edge[, 1] == node & edge[, 2] == left,dim]
    env$`dim`[left] <- env$`dim`[node] + left_branch
    # call this function for node = left
    recursion_space(tree, left, env,dim)
  }
  # recursively dive down to the right subtree (if right isn't a tip)
  if (right > nleaves) {
    # compute the time of the `right` node as the position of the parent plus
    # the branch length connecting those two
    right_branch <- edge[edge[, 1] == node & edge[, 2] == right,dim]
    env$`dim`[right] <- env$`dim`[node] + right_branch
    # call this function for node = right
    recursion_space(tree, right, env,dim)
  }

}

move_nodes <- function(tree,disp_dist=1,bias_x=0,bias_y=0){
  # this function will displace the nodes. It uses the edge lengths and draws from
  # a normal distribution in x and y
  # the tree, after its nodes are moved, will have a new attribute called "nodes"
  # which will hold this spatial information
  tree_tbl <- tree$table %>%
    rowwise() %>%
    mutate(x=rnorm(n=1,mean=bias_x,sd=edge_length*disp_dist),
           y=rnorm(n=1,mean=bias_y,sd=edge_length*disp_dist))
  tree$table <- tree_tbl
  n <- length(tree$tip.label)
  # get the internal node coordinates
  x <- unlist(node_space(tree,"x"))
  y <- unlist(node_space(tree,"y"))
  # get the coordinates of tips
  for (tip in c(1:n)){
    parent <- phangorn::Ancestors(tree,tip,"parent")
    x[tip] <- unlist(x[parent]+tree_tbl[tree_tbl$child==tip & tree_tbl$parent==parent,"x"])
    y[tip] <- unlist(y[parent]+tree_tbl[tree_tbl$child==tip & tree_tbl$parent==parent,"y"])
  }
  node_tbl <- data.frame(node=c(1:(tree$Nnode+length(tree$tip.label))),
                         time=tree$node.times,
                         x=x,
                         y=y) %>%
    rowwise() %>%
    st_as_sf(coords = c("x","y"), remove = FALSE) %>%
    rename("location" = "geometry","node_id"="node")
  tree$nodes <- node_tbl # add a new attribute with node positions
  return(tree)
}

treesim_connect <- function(ts){
  dat <- ts$nodes
  print(dat)
  edg <- ts$edge %>% as.data.frame() %>% rownames_to_column("edge_id") %>% rename("parent"="V1","child"="V2")
  #dat <- ts %>% ts_data()
  parent_nodes <- dat %>%
    dplyr::as_tibble() %>%
    dplyr::filter(node_id %in% edg$parent) %>%
    dplyr::select(parent_node_id = node_id,
                  parent_time = time, parent_location = location) %>%
    dplyr::left_join(edg, by = c("parent_node_id" = "parent")) %>%
    dplyr::arrange(parent_node_id)
  # next we get all the nodes which are children
  branch_nodes <- dat %>%
    dplyr::as_tibble() %>%
    dplyr::filter(node_id %in% edg$child) %>%
    dplyr::select(child_node_id = node_id,
                  child_time = time, child_location = location) %>%
    dplyr::inner_join(parent_nodes, by = c("child_node_id" = "child")) %>%
    # and smash the two dataframes together
    dplyr::arrange(child_node_id)
  # get all the lines connecting the parent and child locations
  connections <- branch_nodes %>%
    rowwise() %>%
    mutate(pts=st_union(child_location,parent_location)) %>%
    mutate(connection=st_cast(pts,"LINESTRING")) %>% ungroup()
  # let the branches hold the line information
  branches <- connections %>%
    sf::st_set_geometry("connection")
  # let the connections hold the point position of the individuals
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
  return(list(connections, branches))
}

find_ancestors <- function(tree, data, stripped=TRUE){
  # a function which estimates the location of internal nodes by getting the midpoint of tips
  # this is done pairwise- so you get many estimates for a simgle node if it has more than
  # two descendants.
  ts_MRCA <- data %>%
    rename("MRCA_location"="location",
           "MRCA_time"="time",
           "MRCA"="node_id") %>%
    select(MRCA,MRCA_time,MRCA_location)

  # renaming node ids to agree between tree and data
  tree$tip.label <- c(1:length(tree$tip.label))

  ts_n1 <- data %>%
    select(location,time,node_id) %>%
    rename("n1_location"="location","n1_time"="time","n1"="node_id") %>%
    mutate(n1=as.factor(n1))

  ts_n2 <- data %>%
    select(location,time,node_id) %>%
    rename("n2_location"="location","n2_time"="time","n2"="node_id") %>%
    mutate(n2=as.factor(n2))

  # now for all pairs, add data
  if (stripped==FALSE){
  pairs <- expand.grid(tree$tip.label,tree$tip.label) %>% # get all pairs of tips
    filter(Var1!=Var2) %>%
    mutate(Var1=as.factor(Var1),Var2=as.factor(Var2)) %>%
    rename("n1"="Var1","n2"="Var2") %>%
    rowwise() %>%
    mutate(MRCA=getMRCA(tree,c(n1,n2))) %>% # find their mrca
    # add data for nodes and mrca
    left_join(ts_MRCA,by="MRCA") %>%
    left_join(ts_n1,by="n1") %>%
    left_join(ts_n2,by="n2") %>%
    #filter(MRCA_time>0) %>%
    # draw lines, find midpoints and errors
    mutate(line=st_cast(st_union(n1_location,n2_location),"LINESTRING"),
           midpoint=st_centroid(st_union(n1_location,n2_location)),
           midtomrca=st_cast(st_union(midpoint,MRCA_location),"LINESTRING"),
           error=st_length(midtomrca),
           MRCA_x=unlist(MRCA_location)[1],
           MRCA_y=unlist(MRCA_location)[2]) %>%
    ungroup() %>%
    group_by(MRCA) %>%
    mutate(multiple=st_centroid(st_union(midpoint))) %>%
    ungroup() %>%
    # check geometry validity
    filter(st_is_valid(line)==TRUE,st_is_valid(midtomrca)==TRUE,st_is_valid(midpoint)==TRUE,st_is_valid(multiple)==TRUE)
  }
  else {
    pairs <- expand.grid(tree$tip.label,tree$tip.label) %>% # get all pairs of tips
      filter(Var1!=Var2) %>%
      mutate(Var1=as.factor(Var1),Var2=as.factor(Var2)) %>%
      rename("n1"="Var1","n2"="Var2") %>%
      rowwise() %>%
      mutate(MRCA=getMRCA(tree,c(n1,n2))) %>% # find their mrca
      # add data for nodes and mrca
      left_join(ts_MRCA,by="MRCA") %>%
      left_join(ts_n1,by="n1") %>%
      left_join(ts_n2,by="n2") %>%
      #filter(MRCA_time>0) %>%
      # draw lines, find midpoints and errors
      mutate(midpoint=st_centroid(st_union(n1_location,n2_location))) %>%
      ungroup() %>%
      group_by(MRCA) %>%
      mutate(multiple=st_centroid(st_union(midpoint))) %>%
      ungroup()
  }
  return(pairs)
}



find_ancestor_rec <- function(tree){
  # initialise some stuff
  nodes <- tree$nodes %>% as.data.frame() # node information
  dists <- dist.nodes(tree) # distance matrix
  # initialise column
  nodes$inf_loc <- NA
  nodes$inf_loc[1:length(tree$tip.label)] <- nodes$location[1:length(tree$tip.label)]

  while (sum(is.na(nodes$inf_loc))!=0){ # termination condition
    unsolved <- is.na(nodes$inf_loc)==TRUE & nodes$node_id>length(tree$tip.label)+1
    # pick which node to solve
    if (sum(is.na(nodes$inf_loc))>1){
      tosolve <- nodes[nodes$time==max(nodes[unsolved,"time"]),"node_id"]
    }
    else tosolve <- length(tree$tip.label)+1 # in case it is the root
    children <- phangorn::Children(tree,tosolve) # find children
    edges <- c(dists[children[1],tosolve],dists[children[2],tosolve]) # get edge lengths
    weights <- c(edges[1]/sum(edges),edges[2]/sum(edges)) # weightings
    locations <- c(nodes$inf_loc[children[1]],nodes$inf_loc[children[2]]) # get children locations
    inf_locn <- weights[1]*unlist(locations[1])+weights[2]*unlist(locations[2]) # place the parent
    nodes$inf_loc[tosolve] <- list(inf_locn) # add to dataframe
  }
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



find_ancestor_rec_ts <- function(tree,data){
  # initialise some stuff
  nodes <- data %>% select(node_id=phylo_id,location,time) %>%
    as.data.frame() %>%
    #mutate(node_id=node_id+1) %>%
    arrange(node_id,decreasing=FALSE) # node information

  dists <- dist.nodes(tree) # distance matrix

  # initialise column
  nodes$inf_loc <- NA
  nodes$inf_loc[1:length(tree$tip.label)] <- nodes$location[1:length(tree$tip.label)]
  nodes$var <- 0

  while (sum(is.na(nodes$inf_loc))!=0){ # termination condition
    unsolved <- is.na(nodes$inf_loc)==TRUE #& nodes$node_id>length(tree$tip.label)
    # pick which node to solve
    if (sum(is.na(nodes$inf_loc))>1){
      tosolve <- nodes[nodes$time==max(nodes[unsolved,"time"]) ,"node_id"] # need to add a tiebreak if they are the same age
      if (length(tosolve)>1) tosolve <- sample(tosolve,1)
    } else tosolve <- length(tree$tip.label)+1 # in case it is the root
    print(tosolve)
    children <- phangorn::Children(tree,tosolve) # find children
    edges <- c(dists[children[1],tosolve],dists[children[2],tosolve]) # get edge lengths
    weights <- c(edges[1]/sum(edges),edges[2]/sum(edges)) # weightings
    locations <- c(nodes$inf_loc[children[1]],nodes$inf_loc[children[2]]) # get children locations
    inf_locn <- weights[1]*unlist(locations[1])+weights[2]*unlist(locations[2]) # place the parent
    tms <- list("p"=nodes[nodes$node_id==tosolve,"time"], # relative times pretending that time goes in one direction
                "c1"=nodes[nodes$node_id==children[2],"time"]-sum(edges), # walking backwards in time from one child
                "c2"=nodes[nodes$node_id==children[2],"time"])
    var <- (tms$c2-tms$p)*(tms$p-tms$c1)/(tms$c2-tms$c1) + # standard deviation of inferred location assuming brownian motion
      nodes[nodes$node_id==children[1],"var"] + # sum uncertainty of children
      nodes[nodes$node_id==children[2],"var"]
    nodes$inf_loc[tosolve] <- list(inf_locn) # add to dataframe
    nodes$var[tosolve] <- var # add to dataframe
    print(paste(sum(is.na(nodes$inf_loc)),"nodes left."))

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


pairwise_distance <- function(tr,data){
  # gets the pairwise distance between nodes, and the time separating them
  nodes <- data
  # make two dataframes with node information
  n1 <- nodes %>%
    select(location,time,node_id) %>%
    rename("n1_location"="location","n1_time"="time","n1"="node_id") %>%
    mutate(n1=as.factor(n1))
  n2 <- nodes %>%
    select(location,time,node_id) %>%
    rename("n2_location"="location","n2_time"="time","n2"="node_id") %>%
    mutate(n2=as.factor(n2))
  # matrix of distances
  dists <- dist.nodes(tr)
  # do some dataframe joining to get distances
  pairs <- expand.grid(c(1:dim(nodes)[1]),c(1:dim(nodes)[1])) %>% # get all pairs of tips
    filter(Var1!=Var2) %>%
    mutate(Var1=as.factor(Var1),Var2=as.factor(Var2)) %>%
    rename("n1"="Var1","n2"="Var2") %>%
    rowwise() %>%
    mutate(timedist=dists[n1,n2]) %>% # distance in time along tree
    left_join(n1,by="n1") %>%
    left_join(n2,by="n2") %>%
    mutate(geodist=st_length(st_cast(st_union(n1_location,n2_location),"LINESTRING")), # distance in space
           timediff=abs(n2_time-n1_time)) %>% # time disparity (not along tree)
    ungroup()

  return(pairs)
}


relcomp <- function(a, b) {

  comp <- vector()

  for (i in a) {
    if (i %in% a && !(i %in% b)) {
      comp <- append(comp, i)
    }
  }

  return(comp)
}















