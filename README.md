# treesinspace

## Exploring the effects of ecological parameters on the spatial structure of genetic tree sequences

Welcome to the `treesinspace` package. This will supply everything needed to reproduce the analyses and figures included in the manuscript. The code was written by M. Ianni-Ravn (mianniravn@uchicago.edu), with abundant help from Martin Petr and the rest of the Racimo lab at UCPH. The scripts run on a previous version of `slendr` (current version at https://cran.r-project.org/web/packages/slendr/index.html). Please also make sure to install this version.

### Installation

To get going, please clone this repository:

```
git clone https://github.com/mkiravn/treesinspace.git
```

Please also install the August 2022 version of slendr:

``` 
git clone https://github.com/bodkan/slendr.git slendr_treesinspace 
```

Once this has installed, please open the `treesinspace.Rproj` project. This will load all necessary modules.

If there are troubles loading `slendr`, please check the path to the package

### Navigating the package

The `new_scripts` directory contains files needed for the analysis. Namely:

* `exploration.R`: exploring the effects of dispersal distribution, competition distance and mating distance.
* `diffusion.R`: comparing a theoretical model of mate choice to simulated results. This also uses data from `data/ThirdSidePDF.tsv`, produced by `thirdsidepdf.nb`.
* `estimation.R`: investigating the accuracy of ML estimation of $\sigma$ under a brownian motion model
* `illustration.R`: a visualisation of the effects of competition and mate choice scale.

These scripts will output figures into a `figs` directory, and label them with the current date. 


>>>>>>> Stashed changes
