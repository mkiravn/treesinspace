set.seed(206)
library(latex2exp)

### A script to simulate trees and get the ML estimate of sigma
# defining parameters
n <- 2000
Ns <- n
mat_dists <- c(0.2, 0.5, 1, 2, 5) # mate distance
comp_dists <- c(0.2) # competition distance
disp_dists <- c(1) # mean dispersal distance
disp_funs <- c("brownian", "cauchy")
reps <- c(1:5)
ngens <- 8000
fileout <- "estimation"
pars <- set_slendr_pars(Ns,mat_dists,comp_dists,disp_dists,reps,ngens,fileout)

# defining a world
map <- world(
  xrange = c(-25, 25),
  # min-max longitude
  yrange = c(-25, 25),
  # min-max latitude
  landscape = "blank"
)

res <- run_slendr_simulation(pars,map=map,pathout="estimation",pairs=T,samplenum = 50)

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

all_connections %>% write_delim("estimation_conns.tsv",delim="\t")
all_dists %>% write_delim(file="estimation_dists.tsv",delim="\t")
# all_dists <- read_delim(file="estimation_dists.tsv",delim="\t")



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

estimates %>% write_delim(file="estimates.tsv",delim="\t")
# pairwise distances
all_dists  %>%
  dplyr::filter(disp_fun == "brownian",simplified!="unsimplified",mat_dist ==0.2) %>%
  ggplot(aes(y = distance, x = edge_gens)) +
  #ggrastr::rasterize(geom_point(alpha = 0.5, size = 0.5,col="grey")) +
  geom_bin2d(alpha = 0.9, size = 0.5) +
  stat_function(
    fun = sqrt,
    alpha = 0.8,
    lty = 2,
    col = "black"
  ) +
  geom_smooth(se = F, method = "loess",col="indianred") +
  theme_minimal() +
  facet_grid( cols=vars(simplified)) +
  lims(y=c(0,15)) +
  labs(x="branch length (generations)",y="geographic distance") +
  scale_fill_gradientn(colours = met.brewer(name = "Hokusai2"))+
  theme(axis.text.x = element_text(angle = 45))+
  scale_x_sqrt() -> pairwise_plot

# how does distance increase with relatedness?
all_dists  %>%
  dplyr::filter(mat_dist == 0.2, disp_fun == "brownian", simplified!="unsimplified") %>%
  ggplot(aes(y = edge_gens, x = distance)) +
  #ggrastr::rasterize(geom_point(alpha = 0.5, size = 0.5,col="grey")) +
  geom_bin2d(alpha = 0.9, size = 0.5) +
  stat_function(
    fun = sqrt,
    alpha = 0.8,
    lty = 2,
    col = "black"
  ) +
  facet_grid(
    rows = vars(mat_dist),
    cols = vars(disp_fun, comp_dist),
    labeller = label_both,
    scales = "free_y"
  ) +
  geom_smooth(se = F, method = "loess") +
  scale_y_sqrt() +
  theme_minimal() +
  facet_wrap( ~ simplified, nrow = 1) +
  labs(x="genetic distance (generations)",y="geographic distance") +
  scale_fill_gradientn(colours = met.brewer(name = "Hokusai2")) -> pairwise_plot2

pairwise_plot %>%  ggsave(
  filename = paste0("figs/", as.character(Sys.Date()), "pairs.pdf"),
  device = "pdf",
  height = 3,
  width = 7
)

 # estimates_cut <- all_dists %>%
 #  group_by(mat_dist, disp_fun, simplified, rep, comp_dist, shortedge = edge_gens <
 #             100, recent) %>%
 #  summarise(n = n(),
 #            sigma_d = sqrt((1 / (2 * n(
 #            ))) *
 #              sum((distance / sqrt(
 #                sqrt(edge_gens)
 #              )) ^ 2))) %>% ungroup()
 #
 # estimates_cut_recent <- all_dists %>% filter(recent) %>%
 #   group_by(mat_dist, disp_fun, simplified, rep, comp_dist) %>%
 #   summarise(n = n(),
 #             sigma_d = sqrt((1 / (2 * n(
 #             ))) *
 #               sum((distance / sqrt(
 #                 sqrt(edge_gens)
 #               )) ^ 2))) %>% ungroup()
 #


ggplot(data = dplyr::filter(estimates,recent==FALSE),
       aes(col = simplified, x = mat_dist, y = sigma_d)) +
  geom_hline(yintercept = 1,col="grey") +
  geom_quasirandom(size = 2,width = 0.1) +
  geom_smooth(se = F, size = 0.5) +
  theme_minimal() +
  facet_grid(
    cols = vars(disp_fun),
    scales = "free_y"
  ) +
  labs(y = TeX("$\\hat{\\sigma}_{ML}$"),
       x="mating distance") -> estimates.plot.matdist


# histograms of edge generations
all_dists %>%
  dplyr::filter(mat_dist == 0.2,
         disp_fun == "brownian",
         simplified != "unsimplified") %>%
  ggplot(aes(x = edge_gens)) +
  geom_histogram(alpha = 0.6,fill="#00BFC4") +
  facet_wrap( ~ simplified, scales = "free") +
  geom_vline(xintercept = 100,
             lty = 2,
             col = "black") +
  theme_minimal() +
  labs(title = "uncut branches",
       x="branch length (generations)") -> h1

all_dists %>%
  dplyr::filter(mat_dist == 0.2,
                disp_fun == "brownian",
                simplified != "unsimplified",
                edge_gens < 100) %>%
  ggplot(aes(x = edge_gens)) +
  geom_histogram(alpha = 0.6,fill="#F8766D") +
  facet_wrap( ~ simplified, scales = "free") +
  geom_vline(xintercept = 100,
             lty = 2,
             col = "black") +
  theme_minimal() +
  labs(title = "cut branches",
       x="branch length (generations)") -> h2
ggarrange(h1,h2,common.legend = T,ncol=1) -> hist.plot


estimates <- estimates %>%
  mutate(simplified_fct = factor(simplified, levels = c('unsimplified', 'simplified', 'tips only')))
all_dists <- all_dists %>%
  mutate(simplified_fct = factor(simplified, levels = c('unsimplified', 'simplified', 'tips only')))


dplyr::filter(all_dists,mat_dist == 0.2,disp_fun=="brownian") %>%
  ggplot(aes(y = distance / (sqrt(edge_gens*pi/2)),x=simplified_fct),
         alpha=0.1,size=0.5) +
  geom_hline(yintercept=1,col="grey")+
  # add violins of all the ensemble of branch-wise contributions
  geom_violin(data=all_dists %>%
                dplyr::filter(mat_dist == 0.2,disp_fun=="brownian",parent_time > ngens-100),
              aes(col="tree cut at \n100 generations \nin past"),
              alpha=0.2,size=0.5,lty=1,fill="aliceblue") +
  geom_violin(aes(col="all branches"),alpha=0.2,size=0.5,fill="mistyrose") +
  # add points for estimates
  geom_point(data = estimates %>%
               dplyr::filter(mat_dist == 0.2,disp_fun=="brownian",recent),
             aes( y = sigma_d,
                  col="tree cut at \n100 generations \nin past"),
             size=2) +
  geom_point(data = estimates %>%
               dplyr::filter(mat_dist == 0.2,disp_fun=="brownian",recent==FALSE),
             aes( y = sigma_d,
                  col="all branches"),
             size=2) +
  theme_minimal() +
  scale_y_log10() +
  labs(shape="simulation replicate",
       col="",x="",y = TeX("$\\hat{\\sigma}_{ML}$")) -> estimates.plot.cut.2

#save.image(file = "estimation.RData")


estimates.plot.matdist %>% ggsave(filename=paste0("figs/",
                                                  as.character(Sys.Date()),
                                                  "estimate_plots_md.pdf"),
                                  device = "pdf",
                                  height = 3,
                                  width = 7
                                  )

estimates.plot.cut.2 %>% ggsave(filename=paste0("figs/",
                                                  as.character(Sys.Date()),
                                                  "estimate_plots_cut.pdf"),
                                  device = "pdf",
                                  height = 5,
                                  width = 7
)
