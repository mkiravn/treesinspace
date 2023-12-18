# treesinspace

## Exploring the Effects of Ecological Parameters on the Spatial Structure of Genealogies

Welcome to the `treesinspace` package. This will supply everything needed to reproduce the analyses and figures included in the manuscript. The code was written by M. Ianni-Ravn (mianniravn@uchicago.edu), with abundant help from Martin Petr and the rest of the Racimo lab at UCPH. The project extensively makes use of the R package `slendr` (current version at https://cran.r-project.org/web/packages/slendr/index.html). This will be installed, if not already present, when opening this R project. 

### Installation

To get going, please clone this repository:

```
git clone https://github.com/mkiravn/treesinspace.git
```

Once this has installed, please open the `treesinspace.Rproj` project. This will load all necessary modules.

### Navigating the package

The `new_scripts` directory contains files needed for the analysis. Namely:

* `exploration.R`: exploring the effects of dispersal distribution, competition distance and mating distance.
* `diffusion.R`: comparing a theoretical model of mate choice to simulated results. This also uses data from `data/ThirdSidePDF.tsv`, produced by `thirdsidepdf.nb`.
* `estimation.R`: investigating the accuracy of ML estimation of $\sigma$ under a brownian motion model
* `illustration.R`: a visualisation of the effects of competition and mate choice scale.

These scripts will output figures into a `figs` directory, and label them with the current date. 

