% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model.helpers.R
\name{tomatrix}
\alias{tomatrix}
\title{Converts an input graph object to an adjacency matrix and identifies its characteristics}
\usage{
tomatrix(graph)
}
\arguments{
\item{graph}{A graph object of class "matrix", "sparseMatrix", \link{igraph}, matrix or dataframe edgelist, or \link[network]{network}}
}
\value{
a list(summary, G)
   `summary` is a dataframe containing characteristics of the supplied object
   `G` is an adjacency matrix
}
\description{
Converts an input graph object to an adjacency matrix and identifies its characteristics
}
\examples{
M <- matrix(rbinom(5*5,1,.5),5,5)
test <- backbone:::tomatrix(M)
}
\keyword{internal}
