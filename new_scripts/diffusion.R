# defining parameters
Ns <- 2000 # number of individuals in population
mat_dists <- 0.1 # mating distance
comp_dists <- 0.2 # competition distance
disp_dists <- 1 # dispersal distance (sigma)
disp_funs <-
  "brownian" # dispersal kernel
reps <-
  1 # number of simulation runs for each parameter combination
ngens <- 100 # number of generations
fileout <- "part2a"
pars <- set_slendr_pars(Ns,mat_dists,comp_dists,disp_dists,reps,ngens,fileout)

# defining a world
map <- world(
  xrange = c(-25, 25),
  # min-max longitude
  yrange = c(-25, 25),
  # min-max latitude
  landscape = "blank"
)


res <- run_slendr_simulation(pars,map=map,pathout="part2")

trees <- res$simplified
trees_us <- res$unsimplified

tree_data <- convert_trees(trees,pars) %>% rbind(convert_trees(trees_us,pars))

# save the tree data
tree_data <- tree_data %>%
  group_by(sim) %>%
  # adjust so that time starts at 0 (for recapitated trees)
  mutate(timeoriginal = time,
         time = time - min(time))

all_connections <- retrieve_connections(trees,trees_us,pars)

# main plots
all_connections %>%
  dplyr::filter(simplified=="unsimplified") %>%
  ggplot(aes(x=distance,lty="simulated",col="simulated")) +
  geom_density(show.legend = F) +
  stat_function(aes(lty="theoretical"),fun = drayleigh, args = list(scale=1),alpha=0.8,col="black") +
  theme_minimal() +
  labs(y="density",x="distance",lty="") -> absdist

all_connections %>%
  dplyr::filter(simplified == "unsimplified") %>%
  select(x = `x displacement`, y =`y displacement`) %>%
  pivot_longer(cols = c("x", "y"),
               names_to = "dimension",
               values_to = "displacement") %>%
  ggplot(aes(x = displacement, lty = "simulated", col = dimension)) +
  geom_density() +
  stat_function(aes(lty = "theoretical"),
                fun = dnorm,
                args = list(mean = 0, sd = 1),
                alpha = 1,
                col = "black"
  ) +
  theme_minimal() +
  labs(y = "density", x = "displacement", lty = "") -> xandy

# adding qqplots
all_connections %>%
  dplyr::filter(simplified=="unsimplified") %>%
  ggplot() +
  stat_qq(aes(sample=distance,col="distance"),distribution = qrayleigh,alpha=0.5)  +
  stat_qq_line(aes(sample=distance),
               distribution = qrayleigh,
               alpha=0.8,
               lty=2,
               col="black",
               dparams = list(scale=1))  +
  labs(y="simulated",x="theoretical") +
  theme_minimal() -> absdist_qq

all_connections %>%
  dplyr::filter(simplified=="unsimplified") %>%
  ggplot() +
  stat_qq(aes(sample=`x displacement`,col="x displacement"),distribution = qnorm,alpha=0.5)  +
  stat_qq(aes(sample=`y displacement`,col="y displacement"),distribution = qnorm,alpha=0.5)  +
  stat_qq_line(aes(sample=`y displacement`),
               distribution = qnorm,
               alpha=0.8,
               lty=2,
               col="black",
               dparams = list(mean=0,sd=1))  +
  labs(y="simulated",x="theoretical") +
  theme_minimal() -> xandy_qq

ggarrange(absdist,absdist_qq,xandy,xandy_qq,common.legend = T) -> brownian_plots

all_connections %>%
  dplyr::filter(simplified=="unsimplified") %>%
  ggplot(aes(x=`x displacement`,y=`y displacement`,col="simulated")) +
  geom_point(alpha=1,size=0.3,show.legend = F) +
  theme_minimal() +
  labs(x="x displacement",y="y displacement")+
  geom_density_2d(alpha=1,show.legend = F,col="grey") +
  stat_cor(col="black",geom="label") +
  coord_fixed() -> brownian_dispersal_angle

ggarrange(brownian_plots,
          brownian_dispersal_angle,
          ncol = 2,
          widths = c(2,1)) -> brownian_plots_panel

# Checking that increments are stationary
all_connections %>%
  mutate(timegroup=cut_width(parent_time,width=25,boundary=0)) %>%
  dplyr::filter(simplified=="unsimplified") %>%
  ggplot(aes(x=distance,col=timegroup)) +
  geom_density() +
  stat_function(aes(col="Rayleigh"),
                fun = drayleigh,
                args = list(scale=1),
                alpha=0.8,
                lty=2, col="black") +
  theme_minimal() +
  labs(y="density",
       col="time grouping",
       x="distance") -> stationary_increments

all_connections %>%
  mutate(timegroup=cut_width(parent_time,width=25,boundary=0)) %>%
  dplyr::filter(simplified=="unsimplified") %>%
  ggplot(aes(sample=distance,col=timegroup)) +
  geom_qq(distribution = qrayleigh,dparams=list(scale=1),alpha=0.8,geom="path") +
  stat_qq_line(aes(sample=distance),distribution = qrayleigh,alpha=0.8,lty=2,col="black",dparams = list(scale=1))  +
  theme_minimal() +
  labs(y="simulated",
       x="theoretical",
       col="time grouping") -> stationary_qqplot


all_connections %>%
  dplyr::filter(simplified=="unsimplified") %>%
  mutate(timegroup=cut_width(parent_time,width=25,boundary=0)) %>%
  ggplot(aes(x=`x displacement`,y=`y displacement`,col=timegroup)) +
  geom_point(alpha=0.5,size=0.1,show.legend = F) +
  theme_minimal() +
  geom_density_2d(show.legend = F) +
  #stat_cor(geom="label",aes(group=timegroup),fill="white",r.accuracy = 0.01,p.accuracy = 0.01,size=4)+
  labs(x="x displacement",y="y displacement",col="time grouping")+
  coord_equal()-> stationary_dispersal_angle

# checking for independence
cp_jumps <- data.frame()

while (dim(cp_jumps)[1]<500){
  row <- sample(x=c(1:dim(dplyr::filter(all_connections,simplified=="unsimplified"))[1]),size=1)
  focal <- dplyr::filter(all_connections,simplified=="unsimplified")[row,"parent_node_id"]
  if (focal %in% cp_jumps$focal_node){
    print("seen already")
  } else{
    downstream_connection <- dplyr::filter(all_connections,parent_node_id==focal$parent_node_id,simplified=="unsimplified")[1,]
    upstream_connection <- dplyr::filter(all_connections,child_node_id==focal$parent_node_id,simplified=="unsimplified")
    toadd <- data.frame(focal_node=focal$parent_node_id,
                        downstream_distance=downstream_connection$distance,
                        downstream_angle=downstream_connection$angle,
                        upstream_distance=upstream_connection$distance,
                        upstream_angle=upstream_connection$angle)
    cp_jumps <- rbind(cp_jumps,toadd) %>% dplyr::filter(!is.na(downstream_distance))
  }
}

cp_jumps %>%
  ggplot(aes(y=downstream_distance,x=upstream_distance)) +
  geom_point(show.legend=F, aes(col=as.factor(focal_node)),size=0.5,alpha=0.5) +
  labs(x="upstream",y="downstream",title="Distance") +
  theme_minimal() +
  stat_cor(geom="label") +

  cp_jumps %>%
  ggplot(aes(y=downstream_angle,x=upstream_angle)) +
  geom_point(show.legend=F, aes(col=as.factor(focal_node)),size=0.5,alpha=0.5) +
  labs(x="upstream",y="downstream",title="Angle") +
  theme_minimal() +
  stat_cor(geom="label") -> independence_plot



mat_dists <- c(0.5,1,2,5,10) # mating distance
fileout <- "part2b"
pars <- set_slendr_pars(Ns,mat_dists,comp_dists,disp_dists,reps,ngens,fileout)

res <- run_slendr_simulation(pars,map=map,pathout="part2b")

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
ggplot(probsm) +
  geom_line(aes(
    x = dist,
    y = p,
    group = mat_dist,
    col = as.numeric(mat_dist)
  ),alpha=0.8) + scale_color_viridis_c() +
  theme_minimal() +
  stat_function(col="black",fun = drayleigh, args = list(scale=1),alpha=0.8,lty=2) +
  labs(col="Mating distance",y="g_y(y)",x="y") -> monster.plot

# probability if both sides are Rayleigh
ThirdSideRR <- function(y,sig,rb){
  0.5*drayleigh(y,scale=sig) + 0.5*drayleigh(y,scale=sqrt(sig^2+rb^2))
}
probsm <- probsm %>% rowwise() %>% mutate(prr=ThirdSideRR(dist,1,as.numeric(mat_dist)))

# plotting theoretical distributions
ggplot(probsm) + geom_line(aes(
  x = dist,
  y = prr,
  group = mat_dist,
  col = as.numeric(mat_dist)
),alpha=0.8) + scale_color_viridis_c() +
  theme_minimal() +
  stat_function(col="black",fun = drayleigh, args = list(scale=1),alpha=0.8,lty=2) +
  labs(col="Standard deviation s_b",y="h_y(y)",x="y") -> rr.plot

garrange(monster.plot,rr.plot,ncol=1,labels="auto") %>%
  ggsave(
    filename = paste0("figs/", as.character(Sys.Date()), "theoretical_plots.pdf"),
    device = "pdf",
    height = 6,
    width = 7
  )


probsm <- probsm %>%
  dplyr::filter(mat_dist %in% mat_dists) %>%
  rename("distance"="dist")
probsm$mat_dist <- factor(probsm$mat_dist, levels = mat_dists)
probsm <- probsm %>% group_by(mat_dist) %>% mutate(cp=cumsum(p))

# Checking that dispersal is normal in x and y and the norm is rayleigh

all_connections <- all_connections %>%
  mutate(mat_dist=as.factor(mat_dist))
# main plots
all_connections %>%
  dplyr::filter(simplified=="unsimplified") %>%
  ggplot(aes(x=distance,col=mat_dist)) +
  geom_density() +
  stat_function(col="black",fun = drayleigh, args = list(scale=1),alpha=0.8,lty=2) +
  theme_minimal() + lims(x=c(0,15)) +
  labs(y="density",
       col="mating distance",
       x="distance") -> brownian_dispersal_simulated
probsm %>%
  ggplot(aes(x=distance,col=mat_dist,y=p)) +
  geom_line() +
  stat_function(col="black",fun = drayleigh, args = list(scale=1),alpha=0.8,lty=2) +
  theme_minimal() + lims(x=c(0,15)) +
  labs(y="density",
       col="mating distance",
       x="distance") -> brownian_dispersal_theoretical
ggarrange(brownian_dispersal_simulated,
          brownian_dispersal_theoretical,
          ncol=2,common.legend = T) -> brownian_dispersal


# adding qqplots
all_connections %>%
  dplyr::filter(simplified=="unsimplified") %>%
  ggplot() +
  stat_qq(aes(sample=distance,col=mat_dist),distribution = qrayleigh,alpha=1,geom="path")  +
  geom_abline(slope=1,intercept=0,lty=2,col="black")  +
  labs(y="simulated",
       x="theoretical",
       col="mating distance") +
  theme_minimal() -> brownian_qqplots

# comparison to theoretical distribution
all_connections %>%
  dplyr::filter(simplified=="unsimplified") %>%
  ggplot(aes(x=distance,col=mat_dist)) +
  geom_line(data=probsm,aes(x=distance,y=p,col=mat_dist,lty="theoretical"))+
  #geom_line(data=probsm,aes(x=dist,y=prr,col=mat_dist,lty="h_y(y)"))+
  geom_density(aes(lty="simulated")) +
  stat_function(col="black",fun = drayleigh, args = list(scale=1),alpha=0.8,lty=2) +
  theme_minimal() +
  facet_wrap(~mat_dist,nrow = 1) +
  labs(y="density",
       x="distance",
       col="mating distance",
       lty="",
       strip.text = element_blank())  -> sim_theo

ggarrange(sim_theo,brownian_qqplots,ncol=1,labels="auto",common.legend = T) %>%
  ggsave(
    filename = paste0("figs/", as.character(Sys.Date()) ,"mating_distance.pdf"),
    device = "pdf",
    height = 5,
    width = 7
  )
