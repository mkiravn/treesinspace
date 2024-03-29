% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tree_simulation.R
\name{find_ancestor_centroid}
\alias{find_ancestor_centroid}
\title{Find ancestor by centroid method}
\usage{
find_ancestor_centroid(tree, data, stripped = TRUE)
}
\arguments{
\item{tree}{a tree of ape class phylo}

\item{data}{the data associated with the tree}

\item{stripped}{give all pairwise estimates or only the consensus}
}
\value{
The cenotrid estimates of all internal nodes of the tree
}
\description{
Find ancestor by centroid method
}
