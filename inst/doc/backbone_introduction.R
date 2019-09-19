## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup---------------------------------------------------------------
library(backbone)

## ------------------------------------------------------------------------
data(davis) #load the dataset
op <- options(width = 100)
davis #view the dataset
options(op)

## ------------------------------------------------------------------------
davis%*%t(davis) #The projected davis dataset

## ------------------------------------------------------------------------
G <- davis%*%t(davis) #projected davis dataset, a weighted graph
universal_bb <- universal(G, upper = 0)
graph <- igraph::graph_from_adjacency_matrix(universal_bb, mode = "undirected")
op <- par(mar=c(0,0,0,0))
lo <- igraph::layout_(graph, igraph::with_fr())
plot(graph, vertex.label = 1:18, layout = lo)
par(op)

## ------------------------------------------------------------------------
universal_bb <- universal(davis, upper = 0, bipartite = TRUE)
graph <- igraph::graph_from_adjacency_matrix(universal_bb, mode = "undirected")
op <- par(mar=c(0,0,0,0))
plot(graph, vertex.label = 1:18, layout = lo)
par(op)

## ------------------------------------------------------------------------
universal_bb <- universal(davis, upper = 4, lower = 2, bipartite = TRUE)
universal_bb

## ------------------------------------------------------------------------
universal_bb <- universal(davis, 
                          upper = function(x)mean(x)+sd(x), 
                          lower=function(x)mean(x)-sd(x), 
                          bipartite = TRUE)

## ------------------------------------------------------------------------
hyperg_probs <- hyperg(davis)
hyperg_bb <- backbone.extract(hyperg_probs$positive, hyperg_probs$negative)

## ----echo=T, results='hide'----------------------------------------------
fdsm_props <- fdsm(davis, trials = 100, sparse = TRUE, dyad=c(1,5))

## ------------------------------------------------------------------------
fdsm_props$dyad_values
fdsm_bb <- backbone.extract(fdsm_props$positive, fdsm_props$negative, alpha = 0.05)
fdsm_bb

## ----echo=T, results='hide'----------------------------------------------
sdsm_props <- sdsm(davis, trials = 100, dyad = c("EVELYN", "CHARLOTTE")) 

## ------------------------------------------------------------------------
sdsm_props$dyad_values
sdsm_bb <- backbone.extract(sdsm_props$positive, alpha = 0.05) 
sdsm_bb

