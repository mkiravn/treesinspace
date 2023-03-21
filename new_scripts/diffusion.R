# defining parameters
Ns <- 2000 # number of individuals in population
mat_dists <- c(1,2,5,10) # mating distance
comp_dists <- 0.2 # competition distance
disp_dists <- 1 # dispersal distance (sigma)
disp_funs <-
  "brownian" # dispersal kernel
reps <-
  1 # number of simulation runs for each parameter combination
ngens <- 100 # number of generations
fileout <- "diffusion"
pars <- set_slendr_pars(Ns,mat_dists,comp_dists,disp_dists,reps,ngens,fileout)

# defining a world
map <- world(
  xrange = c(-25, 25),
  # min-max longitude
  yrange = c(-25, 25),
  # min-max latitude
  landscape = "blank"
)

pars <- set_slendr_pars(Ns,mat_dists,comp_dists,disp_dists,reps,ngens,fileout)

res <- run_slendr_simulation(pars,map=map,pathout="diffusion2",pairs = F)

trees <- res$simplified
trees_us <- res$unsimplified

tree_data <- convert_trees(trees,pars) %>% rbind(convert_trees(trees_us,pars))

# save the tree data
tree_data <- tree_data %>%
  group_by(sim) %>%
  # adjust so that time starts at 0 (for recapitated trees)
  mutate(timeoriginal = time,
         time = time - min(time))

all_connections <- all_connections %>% rbind(retrieve_connections(trees,trees_us,pars))

# reading in theoretical results
probsm <- read_delim("data/ThirdSidePDF.tsv",col_names = seq(1,10,0.2))
colnames(probsm) <- seq(1,10,0.2)
probsm <- probsm %>% as.data.frame() %>% cbind(data.frame(dist=seq(0,20,0.5)))
probsm <- probsm  %>% pivot_longer(cols = -dist,names_to = "mat_dist",values_to = "p")

# probability if both sides are Rayleigh
ThirdSideRR <- function(y,sig,rb){
  0.5*drayleigh(y,scale=sig) + 0.5*drayleigh(y,scale=sqrt(sig^2+rb^2))
}
probsm <- probsm %>% rowwise() %>% mutate(prr=ThirdSideRR(dist,1,as.numeric(mat_dist)))

ggarrange(monster.plot,rr.plot,ncol=1,labels="auto") %>%
  ggsave(
    filename = paste0("figs/", as.character(Sys.Date()), "theoretical_plots.pdf"),
    device = "pdf",
    height = 6,
    width = 7
  )


probsm <- probsm %>%
  dplyr::filter(mat_dist %in% mat_dists)
probsm$mat_dist <- factor(probsm$mat_dist, levels = mat_dists)
probsm <- probsm %>% group_by(mat_dist) %>% mutate(cp=cumsum(p))

all_connections <- all_connections %>%
  mutate(mat_dist=as.factor(mat_dist))

# the plot to output properly
probsm %>% rename("mating distance"="mat_dist") %>%
  ggplot() +
  stat_function(col="orchid",fun = drayleigh, args = list(scale=1),alpha=0.8,lty=1) +
  geom_line(aes(x=dist,y=p,lty="theoretical")) +
  geom_density(data=dplyr::filter(all_connections,
                                  simplified=="unsimplified",
                                  mat_dist %in% factor(mat_dists)),
               aes(x=distance,lty="simulated")) +
  theme_minimal() + lims(x=c(0,15)) +
  facet_grid(cols=vars(`mating distance`),labeller=label_both) +
  labs(y="density",
       lty="",
       x="distance") -> brownian_dispersal_output

ggsave(brownian_dispersal_output,
       filename = paste0("figs/", as.character(Sys.Date()) ,"sim_theor_distance.pdf"),
       device = "pdf",
       height = 3,
       width = 7
)



