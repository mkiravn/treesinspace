set.seed(206)
library(latex2exp)

# defining parameters
n <- 2000
Ns <- n
mat_dists <- c(0.2) # mate distance
comp_dists <- c(0.2) # competition distance
disp_dists <- c(0.2, 0.5, 1, 2, 5) # mean dispersal distance
disp_funs <- c("brownian")
reps <- c(1:5)
ngens <- 8000
fileout <- "estimation_othersig"
pars <- set_slendr_pars(Ns,mat_dists,comp_dists,disp_dists,reps,ngens,fileout)

# defining a world
map <- world(
  xrange = c(-25, 25),
  # min-max longitude
  yrange = c(-25, 25),
  # min-max latitude
  landscape = "blank"
)

res <- run_slendr_simulation(pars,map=map,pathout="estimation_othersig",pairs=T,samplenum = 50)

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
all_connections <- all_connections %>% mutate(edge_gens=child_time-parent_time)

all_connections %>% write_delim("all_connections_1507.tsv",delim="\t")

pwds <- res$pairs
pwds %>%
  mutate(simplified = "tips only",parent_time=ngens-edge_gens/2) %>%
  distinct(.keep_all = T, edge_gens, dist) %>%
  select(mat_dist,
         comp_dist,
         disp_fun,
         sigma,
         simplified,
         sim,
         rep,
         distance=dist,
         edge_gens,
         parent_time) %>%
  rbind(
    all_connections %>%
      select(
        mat_dist,
        comp_dist,
        disp_fun,
        sigma,
        simplified,
        sim,
        rep,
        distance,
        edge_gens,
        parent_time
      )
  ) -> all_dists

all_dists %>% write_delim(file="all_dists_othersig.tsv",delim="\t")


#  Estimates from unsimplified trees, simplified trees and tips only
estimates <- all_dists %>%
  group_by(sigma, mat_dist, disp_fun, simplified, rep, comp_dist) %>%
  summarise(n = n(),
            sigma_d = sqrt((1 / (2 * n(
            ))) *
              sum((distance /
                     sqrt(edge_gens)
              ) ^ 2))) %>% ungroup() %>%
  mutate(recent=FALSE)

estimates_recent <- all_dists %>%
  group_by(sigma, mat_dist, disp_fun, simplified, rep, comp_dist, recent = parent_time > ngens-100) %>%
  summarise(n = n(),
            sigma_d = sqrt((1 / (2 * n(
            ))) *
              sum((distance /
                     sqrt(edge_gens)
              ) ^ 2))) %>% ungroup() %>%
  dplyr::filter(recent)

estimates <- rbind(estimates_recent,estimates)


estimates %>% write_delim(file="estimates_othersig.tsv",delim="\t")

estimates <- estimates %>%
  mutate(simplified_fct = factor(simplified, levels = c('unsimplified', 'simplified', 'tips only')))
all_dists <- all_dists %>%
  mutate(simplified_fct = factor(simplified, levels = c('unsimplified', 'simplified', 'tips only')))


all_dists %>%
  ggplot(aes(y = distance / (sqrt(edge_gens*pi/2)),x=sigma,shape=simplified_fct),
         alpha=0.1,size=0.5) +
  geom_abline(slope=1,intercept=0,col="grey",lty=2)+
  #add points for estimates
  geom_point(data = estimates %>%
               dplyr::filter(recent),
             aes( y = sigma_d,
                  col="tree cut at \n100 generations \nin past"),
             size=2) +
  geom_point(data = estimates %>%
               dplyr::filter(recent==FALSE),
             aes( y = sigma_d,
                  col="all branches"),
             size=2) +
  theme_minimal() +
  lims(y=c(0,5)) +
  scale_y_log10() +
  scale_x_log10() +
  #facet_wrap(~sigma,nrow=1) +
  labs(col="",x=TeX("$\\sigma$"),y = TeX("$\\hat{\\sigma}_{ML}$"),shape="") -> estimates.plot.cut

estimates.plot.cut %>% ggsave(filename=paste0("figs/",
                                                as.character(Sys.Date()),
                                                "estimate_plots_cut_othersig.pdf"),
                                device = "pdf",
                                height = 3,
                                width = 7
)
