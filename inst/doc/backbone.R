## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(backbone)

## -----------------------------------------------------------------------------
data(davis) #load the dataset
op <- options(width = 100)
davis #view the dataset
options(op)

## -----------------------------------------------------------------------------
davis%*%t(davis) #The projected davis dataset

## -----------------------------------------------------------------------------
G <- davis%*%t(davis) #projected davis dataset, a weighted graph
universal_bb <- universal(G, upper = 0)
universal_bb$backbone
universal_bb$summary

## -----------------------------------------------------------------------------
graph <- igraph::graph_from_adjacency_matrix(universal_bb$backbone, mode = "undirected")
op <- par(mar=c(0,0,0,0))
lo <- igraph::layout_(graph, igraph::with_fr())
plot(graph, vertex.label = 1:18, layout = lo)
par(op)

## -----------------------------------------------------------------------------
universal_bb <- universal(davis, upper = 0, bipartite = TRUE)
universal_bb$summary
graph <- igraph::graph_from_adjacency_matrix(universal_bb$backbone, mode = "undirected")
op <- par(mar=c(0,0,0,0))
plot(graph, vertex.label = 1:18, layout = lo)
par(op)

## -----------------------------------------------------------------------------
universal_bb <- universal(davis, upper = 4, lower = 2, bipartite = TRUE)
universal_bb$backbone

## -----------------------------------------------------------------------------
universal_bb <- universal(davis, 
                          upper = function(x)mean(x)+sd(x), 
                          lower=function(x)mean(x)-sd(x), 
                          bipartite = TRUE)

## -----------------------------------------------------------------------------
hyperg <- hyperg(davis)
hyperg_bb <- backbone.extract(hyperg, signed = TRUE)

## ----echo=T, results='hide'---------------------------------------------------
sdsm <- sdsm(davis)

## -----------------------------------------------------------------------------
sdsm_bb <- backbone.extract(sdsm, signed = FALSE, alpha = 0.1, fwer = "bonferroni") 
sdsm_bb

## ----echo=T, results='hide'---------------------------------------------------
fdsm <- fdsm(davis, trials = 100, dyad=c(1,5))

## -----------------------------------------------------------------------------
fdsm$dyad_values
fdsm_bb <- backbone.extract(fdsm, signed = TRUE, alpha = 0.1)
fdsm_bb

