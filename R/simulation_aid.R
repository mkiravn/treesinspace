retrieve_connections <- function(trees=NULL,trees_us=NULL,pars,nice_names=F){
  # extract connections from the results of a simulation
  if (is.null(trees_us)){
    print("Only simplified trees provided.")
  }
  if (is.null(trees)){
    print("Only unsimplified trees provided.")
  }
  if (!is.null(trees)&!is.null(trees_us)){
    print("Both simplified and unsimplified trees provided")
  }
  print("Connecting...")
  # apply `ts_edges` to extract connections
  conns_simp <- lapply(trees, ts_edges)
  conns <- lapply(trees_us, ts_edges)
  all_connections_l <- data.frame() # this will hold the connections
  for (row in c(1:dim(pars)[1])) {
    print(paste("Connecting tree number:", row))
    if (!is.null(trees_us)){
      all_connections_i <- conns[[row]] %>%
        as.data.frame() %>%
        mutate(simplified = "unsimplified",
               distance=st_length(connection)) %>%
        cbind(pars[row,])
    }
    if (!is.null(trees)){
      all_connections_i <- conns_simp[[row]] %>%
        as.data.frame() %>%
        mutate(simplified = "simplified",
               distance=st_length(connection)) %>%
        cbind(pars[row,]) %>%
        rbind(all_connections_i)
    }
    all_connections_l <- rbind(all_connections_l, all_connections_i)
  }
  # get x and y displacements (missing in `ts_edges`)
  all_connections_l <- all_connections_l %>%
    rowwise() %>%
    mutate(x_dist=unlist(child_location)[1]-unlist(parent_location)[1],
           y_dist=unlist(child_location)[2]-unlist(parent_location)[2],
           angle=atan(x_dist/y_dist)) %>%
    ungroup()
  if (nice_names==T){
    all_connections_l <- all_connections_l  %>%
      rename(
        "competition distance" = "comp_dist",
        "mating distance" = "mat_dist",
        "sigma" = "sigma",
        "dispersal function" = "disp_fun",
        "x displacement" = "x_dist",
        "y displacement" = "y_dist"
      )}
  return(all_connections_l)
}

set_slendr_pars <- function(Ns,mat_dists,comp_dists,disp_dists,reps,ngens,fileout=NULL){
  # function to generate and save a parameter table.
  # generating a parameter table
  pars <- expand.grid(
    N = Ns,
    mat_dist = mat_dists,
    comp_dist = comp_dists,
    disp_fun = disp_funs,
    sigma = disp_dists,
    rep = reps,
    ngen = ngens
  ) %>%
    rownames_to_column("sim")
  print(paste("Running" , dim(pars)[1], "simulations."))
  dir.create("./logs")
  write_delim(
    x = pars,
    file = paste0("./logs/",as.character(Sys.Date()),
                  ifelse(is.null(fileout),"pars.tsv",paste0(fileout,".tsv"))),
    delim = "\t"
  )
  return(pars)
}

convert_trees <- function(trees,pars){
  tree_data <- data.frame()
  for (row in c(1:dim(pars)[1])) {
    tree <- trees[[row]]
    tree_data_i <-
      tree %>% ts_nodes() %>%
      rowwise() %>%
      mutate(
        N = pars[row, "N"],
        fun = pars[row, "disp_fun"],
        rep = as.factor(pars[row, "rep"]),
        mat_dist =
          pars[row, "mat_dist"],
        comp_dist =
          pars[row, "comp_dist"],
        sigma =
          pars[row, "sigma"],
        sim = as.factor(row),
        x = unlist(location)[1],
        y = unlist(location)[2],
        simplified =
          "simplified"
      ) %>% ungroup()

    tree_data <- rbind(tree_data, tree_data_i)
    print(paste("Converted tree no", row))
  }
  return(tree_data)
}


run_slendr_simulation <- function(pars,map,pathout,samplenum=NULL,pairs=F){
  # these will hold our data
  trees <- c()
  trees_us <- c()
  i <- 1
  for (row in c(1:dim(pars)[1])) {
    # define a population
    pop <- population(
      name="POP",
      time = 1,
      N = pars[row, "N"],
      center = c(0, 0),
      radius = 50,
      map = map,
      mating = pars[row, "mat_dist"],
      competition = pars[row, "comp_dist"],
      dispersal_fun = as.character(pars[row, "disp_fun"]),
      dispersal = pars[row, "sigma"]
    )

    # compile and run the model
    model <- compile_model(
      populations = pop,
      generation_time = 1,
      simulation_length = ngens,
      resolution = 1,
      # resolution in "distance units per pixel"
      path = paste0("./logs/",pathout),
      overwrite = TRUE,
      force=TRUE
    )

    # sampling strategy
    samples <-
      schedule_sampling(model, times = ngens,
                        list(pop, ifelse(is.null(samplenum),
                                         pars[row, "N"] / 2,
                                         samplenum)))
    # simulate
    try_again(
      10,
      slim(
        model,
        samples = samples,
        # simulate only a single locus
        sequence_length = 1,
        recombination_rate = 0,
        method = "batch",
        # change to "gui" to execute the model in SLiMgui
        random_seed = sample(c(1:100)),
        verbose = FALSE,
        coalescent_only = FALSE,
        burnin = 0,
      ) -> tsu
    )


    # extract simplified trees
    ts <- tsu %>%
      ts_simplify()

    # collect trees
    trees <- c(trees, ts)
    trees_us <- c(trees_us, tsu)

    if (pairs==T){
      pwds <- data.frame()
      # get pairwise distances
      pwd <- ts %>% ts_recapitate(Ne = pars[row, "N"],
                                  recombination_rate = 0,
                                  random_seed = 206) %>%
        ts_phylo(i=1) %>% pairwise_distance()
      pwds <- pwd %>% cbind(pars[row, ]) %>% rbind(pwds)
    }

    print(paste("Simulation number:", i))
    i <- i + 1
  }
  if (pairs==T){return(list(simplified=trees,unsimplified=trees_us,pairs=pwds))}
      else{return(list(simplified=trees,unsimplified=trees_us))}
}

