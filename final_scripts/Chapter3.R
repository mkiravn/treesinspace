# Loading libraries
devtools::load_all("~/treesinspace")
library(tidyverse)
library(ape)
library(sf)
library(MetBrewer)
library(ggrepel)
library(ggpubr)
library(ggtree)
library(moments)
library(phangorn)
library(apTreeshape)
library(VGAM)
library(gt)
library(ggstatsplot)
devtools::load_all("~/slendr")
setup_env()
check_env()
set.seed(206)

######### Here we define parameters

# defining parameters
Ns <- 1000 # number of individuals in population
mat_dists <- c(1) # mating distance
comp_dists <- c(0.5) # competition distance
disp_dists <- c(1) # dispersal distance (sigma)
disp_funs <-
  c("brownian") # dispersal kernel
reps <-
  c(1:5) # number of simulation runs for each parameter combination
ngens <- 4000 # number of generations
sampling <- c("present", "ancient", "staggered")
migration <- c("stationary", "migrated")

# generating a parameter table
pars <- expand.grid(
  N = Ns,
  mat_dist = mat_dists,
  comp_dist = comp_dists,
  disp_fun = disp_funs,
  disp_dist = disp_dists,
  rep = reps,
  ngen = ngens,
  sampling = sampling,
  migration = migration
) %>%
  rownames_to_column("sim")
print(paste("Running" , dim(pars)[1], "simulations."))
write_delim(
  x = pars,
  file = paste0(as.character(Sys.Date()), "-Chapter3.tsv"),
  delim = "\t"
)

######### Here we run the simulations


# these will hold our data
trees <- c()
trees_us <- c()
all_locations <- c()


# defining a world
map <- world(
  xrange = c(-50, 150),
  # min-max longitude
  yrange = c(-50, 50),
  # min-max latitude
  landscape = "blank"
)

i <- 1

for (row in c(1:dim(pars)[1])) {
  print(paste("Simulation number:", i))
  i <- i + 1
  # define a population
  pop <- population(
    "POP",
    time = 1,
    N = pars[row, "N"],
    center = c(0, 0),
    radius = 50,
    map = map,
    mate_dist = pars[row, "mat_dist"],
    competition_dist = pars[row, "comp_dist"],
    dispersal_fun = as.character(pars[row, "disp_fun"]),
    dispersal_dist = pars[row, "disp_dist"]
  )

  # adding a migration
  if (pars[row, "migration"] == "migrated") {
    pop <- pop %>% move(
      trajectory = list(c(0, 0), c(100, 0)),
      start = 1,
      end = ngens,
      snapshots = 5
    )
  }


  # compile and run the model
  model <- compile(
    populations = pop,
    generation_time = 1,
    sim_length = ngens,
    resolution = 1,
    # resolution in "distance units per pixel"
    path = "~/Desktop/Chapter3/",
    overwrite = TRUE
  )



  # sampling strategy
  if (pars[row, "sampling"] == "present") {
    modern_samples <-
      sampling(model, times = ngens, list(pop, 50))
  } else{
    modern_samples <- sampling(model, times = ngens, list(pop, 40))
  }

  if (pars[row, "sampling"] == "ancient") {
    ancient_samples <-
      sampling(model, times = c(ngens - 1000), list(pop, 10))
  } else{
    ancient_samples <- data.frame()
  }

  if (pars[row, "sampling"] == "staggered") {
    ancient_samples_staggered <- sampling(model,
                                          times = runif(
                                            n = 10,
                                            min = ngens - 2000,
                                            max = ngens
                                          ), list(pop, 1))
  } else{
    ancient_samples_staggered <- data.frame()
  }

  # simulate
  try_again(10,
            {
              slim(
                model,
                sampling = rbind(
                  modern_samples,
                  ancient_samples,
                  ancient_samples_staggered
                ),
                # simulate only a single locus
                sequence_length = 1,
                recombination_rate = 0,
                method = "batch",
                # change to "gui" to execute the model in SLiMgui
                random_seed = sample(c(1:100)),
                verbose = FALSE,
                retainall = TRUE,
                burnin = 0,
              )


              # extract simplified trees
              ts <- ts_load(model) %>%
                ts_recapitate(
                  Ne = pars[row, "N"],
                  recombination_rate = 0,
                  random_seed = 207
                ) %>%
                ts_simplify()

              # extract unsimplified tree
              tsu <- ts_load(model, simplify = FALSE)

              # collect trees
              trees <- c(trees, ts)
              trees_us <- c(trees_us, tsu)

              print("locating")
              # need to get data via ts_phylo
              tree <- ts %>% ts_phylo(1)
            })

  data <- tree %>% ts_data() %>% as.data.frame()

  rec <- find_ancestor_recursive(tree, data = data)
  r <- rec$nodes %>%
    rename("recursive_inferred" = "inf_loc")
  s <- find_ancestor_centroid(tree, data = data)
  # comparing the two methods
  locs <- s %>%
    select(node_id = MRCA,
           centroid_inferred = multiple) %>%
    distinct(node_id, centroid_inferred) %>%
    full_join(r, by = "node_id") %>%
    rowwise() %>%
    filter(is.na(var) == F) %>%
    mutate(
      location = st_cast(location, "POINT"),
      centroid_line = st_cast(st_union(centroid_inferred, location), "LINESTRING"),
      centroid_error = st_length(centroid_line)
    ) %>%
    mutate(
      recursive_line = st_cast(st_union(recursive_inferred, location), "LINESTRING"),
      recursive_error = st_length(recursive_line)
    ) %>%
    mutate(true_x = unlist(location)[[1]],
           true_y = unlist(location)[[2]]) %>%
    ungroup() %>% cbind(pars[row,])
  # getting nodes in the past
  locs <- locs %>% filter(time != max(time))

  all_locations <- rbind(all_locations, locs)
}

######### Plotting
all_locations <-
  all_locations %>% filter(recursive_error > 0, centroid_error > 0)

# how does recursive do?
all_locations %>%
  #select(-centroid_error) %>%
  filter(time > 0) %>%
  ggplot(aes(x = 1000 - time, y = recursive_error)) +
  geom_point(alpha = 0.5, size = 0.5) +
  facet_grid(cols = vars(sampling, migration), labeller = label_both) +
  labs(x = "years in past", y = "error (recursive method)") +
  geom_smooth() +
  theme_minimal() -> rec.time

rec.time %>% ggsave(
  filename = paste0(
    "figs/",
    as.character(Sys.Date()),
    "inferring_location_rec.pdf"
  ),
  device = "pdf",
  height = 5,
  width = 10
)


all_locations %>% filter(time > 0) %>%
  ggplot(aes(x = recursive_error, y = centroid_error, col = time)) +
  geom_abline(slope = 1,
              intercept = 0,
              lty = 2) +
  geom_point(size = 0.7,alpha=0.7) +
  geom_line(
    stat="smooth",
    method = "loess",
    col = met.brewer(name = "Austria",n = 1),
    se = F,
    size = 0.8
  ) +
  facet_grid(cols = vars(sampling),rows=vars( migration), labeller = label_value,scales="free_y") +
  scale_colour_gradientn(colours = met.brewer(name = "Hokusai1")) +
  labs(x = "error (recursive method)", y = "error (centroid method)") +
  theme_minimal() -> methods.comp

all_locations %>% filter(time > 0) %>%
  rename("centroid" = "centroid_error",
         "recursive" = "recursive_error") %>%
  pivot_longer(
    cols = c("recursive",
             "centroid"),
    names_to = "method",
    values_to = "error"
  ) %>%
  ggplot(aes(x = time, y = error, col = method)) +
  geom_line(aes(group=node_id+time),col="grey",alpha=0.8) +
  facet_grid(cols = vars(sampling),rows=vars(migration), labeller = label_value,scales="free_y") +
  geom_point(alpha = 1, size = 0.8) +
  geom_line(stat="smooth",method="gam",size=0.5,se=F,col="black",aes(lty=method)) +
  scale_x_reverse() +
  theme_minimal() -> methods.time

# add grey dots in background to show where the population was
all_locations  %>% filter(time > 0, disp_fun == "brownian") %>%
  filter(st_is_valid(recursive_line) ==T) %>%
  ggplot() +
  #geom_sf(data=filter(all_connections,rep==1,migration=="migrated",simplified=="simplified"),
  #        aes(geometry=parent_location),col="grey") +
  geom_sf(
    aes(geometry = recursive_line),
    col = "red",
    lty = 1,
    alpha = 0.2
  ) +
  geom_sf(
    aes(geometry = centroid_line),
    col = "blue",
    lty = 1,
    alpha = 0.2
  ) +
  geom_sf(
    aes(geometry = recursive_inferred),
    col = "red",
    shape = 21,
    fill = "white"
  ) +
  geom_sf(
    aes(geometry = centroid_inferred),
    col = "blue",
    shape = 21,
    fill = "white"
  ) +
  geom_sf(aes(geometry = location,shape=as.factor(rep), col = time), size = 0.8) +
  facet_grid(cols = vars(sampling),rows=vars( migration), labeller = label_value) +
  scale_colour_gradientn(colours = met.brewer(name = "Hokusai1")) +
  labs(x = "Eastings", y = "Northings") + guides(shape="none") +
  theme_minimal() -> points.plot

points.plot %>%
  ggsave(
    filename = paste0(
      "figs/",
      as.character(Sys.Date()),
      "-inferring_location_points.pdf"
    ),
    device = "pdf",
    height = 5,
    width = 7
  )

ggarrange(
  methods.time,
  methods.comp,
  labels = "auto",
  ncol = 1
) %>%
  ggsave(
    filename = paste0(
      "figs/",
      as.character(Sys.Date()),
      "-inferring_location_methods.pdf"
    ),
    device = "pdf",
    height = 7,
    width = 7
  )

# which conditions of sampling and migration are best?
all_locations %>% filter(time > 0, disp_fun == "brownian") %>%
  rename("centroid" = "centroid_error",
         "recursive" = "recursive_error") %>%
  pivot_longer(
    cols = c("recursive",
             "centroid"),
    names_to = "method",
    values_to = "error"
  ) %>%
  ggplot(aes(
    y = error,
    x = migration,
    col = method,
    fill = sampling
  )) +
  geom_boxplot(alpha = 0.6) +
  theme(axis.text.x = element_text(angle = 45)) +
  stat_summary(
    aes(col = method, fill = sampling),
    fun.y = mean,
    geom = "point",
    shape = 21,
    size = 2,
    position = position_dodge(width = 0.75)
  ) +
  theme_minimal() +
  scale_colour_grey(end = 0.6) -> boxplot

boxplot %>% ggsave(
  filename = paste0(
    "figs/",
    as.character(Sys.Date()),
    "inferring_location_boxplots.pdf"
  ),
  device = "pdf",
  height = 5,
  width = 7
)

# ggpaired(all_locations %>% filter(time>0,disp_fun=="brownian"),
#         cond1="recursive_error",cond2="centroid_error",col="sampling",
#         line.size=0.2,line.color="grey",xlab="method",ylab="error")+
#   stat_compare_means(method = "t.test", paired = TRUE,vjust =3,hide.ns = T)+
#   facet_grid(cols=vars(migration),rows=vars(sampling))
#
# grouped_ggbetweenstats(
#   data=all_locations,
#   x=sampling,
#   y=recursive_error,
#   grouping.var=migration,
#   type="p"
# ) -> pp

# all_locations %>% filter(time>0,sampling=="present") %>%
#   pivot_longer(cols=c("recursive_error",
#                       "centroid_error"),
#                names_to="method",
#                values_to="error") %>%
#   grouped_ggwithinstats(
#     x               = method,
#     y               = error,
#     type            = "p",
#     grouping.var    = migration) ->p1
# all_locations %>% filter(time>0,sampling=="ancient") %>%
#   pivot_longer(cols=c("recursive_error",
#                       "centroid_error"),
#                names_to="method",
#                values_to="error") %>%
#   grouped_ggwithinstats(
#     x               = method,
#     y               = error,
#     type            = "p",
#     grouping.var    = migration) ->p2

# boring tables to put in supplement
all_locations %>% filter(time > 0) %>%
  as.data.frame() -> all_locs_

compare_means(
  c(recursive_error,
    centroid_error) ~ sampling,
  all_locs_,
  method = "t.test",
  group.by = "migration"
) %>% select(-method, -p.format, -p.adj) %>% as_tibble() %>%
  gt() %>%
  tab_header(title = "Contrasts") %>% gtsave(paste0(
    "figs/",
    as.character(Sys.Date()),
    "results_findingancestors2.png"
  ))
compare_means(
  c(recursive_error,
    centroid_error) ~ migration,
  all_locs_,
  method = "t.test",
  group.by = "sampling"
) %>%
  select(-method, -p.format, -p.adj) %>%
  as_tibble() %>%
  gt() %>%
  tab_header(title = "Contrasts") %>%
  gtsave(paste0(
    "figs/",
    as.character(Sys.Date()),
    "results_findingancestors3.png"
  ))

# which method does best?
# all_locations %>% filter(time>0) %>%
#   pivot_longer(cols=c("recursive_error","centroid_error"),names_to="method",values_to="error") %>%
#   ggplot(aes(y=error,x=method,col=sampling)) +
#   geom_boxplot() +
#   stat_compare_means(method="t.test",hide.ns = F) +
#   facet_grid(rows=vars(migration),cols=vars(sampling))+
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45)) +
#   stat_summary(fun.y=mean, geom="point", shape=21, size=3,fill=NA) -> boxplot.2

all_locations %>% filter(time > 0) %>%
  pivot_longer(
    cols = c("recursive_error", "centroid_error"),
    names_to = "method",
    values_to = "error"
  ) %>%
  as.data.frame() -> all_locs_pivot
compare_means(
  error ~ method,
  all_locs_pivot,
  method = "t.test",
  group.by = c("sampling", "migration")
) %>%
  select(-method, -p.format, -p.adj) %>% as_tibble() %>%
  gt() %>%
  tab_header(title = "Contrasts") %>%
  gtsave(paste0(
    "figs/",
    as.character(Sys.Date()),
    "results_findingancestors.png"
  ))
