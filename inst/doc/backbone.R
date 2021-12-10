## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
knitr::opts_knit$set(global.par = TRUE)

## ---- echo = FALSE------------------------------------------------------------
set.seed(5)
par(mar = c(0, 0, 0, 0) + 0.1)

## ----setup--------------------------------------------------------------------
library(backbone)

## -----------------------------------------------------------------------------
dat <- matrix(runif(100),10,10)  #Some data

backbone.suggest(dat)  #What should I do?

backbone <- backbone.suggest(dat, s = 0.05)  #Or, just do it

## -----------------------------------------------------------------------------
B <- rbind(cbind(matrix(rbinom(250,1,.8),10),
                 matrix(rbinom(250,1,.2),10),
                 matrix(rbinom(250,1,.2),10)),
           cbind(matrix(rbinom(250,1,.2),10),
                 matrix(rbinom(250,1,.8),10),
                 matrix(rbinom(250,1,.2),10)),
           cbind(matrix(rbinom(250,1,.2),10),
                 matrix(rbinom(250,1,.2),10),
                 matrix(rbinom(250,1,.8),10)))

## -----------------------------------------------------------------------------
B[1:5,1:5]

## -----------------------------------------------------------------------------
rowSums(B)
colSums(B)

## -----------------------------------------------------------------------------
P <- B%*%t(B)
plot(igraph::graph_from_adjacency_matrix(P, mode = "undirected", diag = FALSE, weighted = TRUE), vertex.label = NA)

## -----------------------------------------------------------------------------
bb <- sdsm(B, alpha = 0.075, narrative = TRUE, class = "igraph")

## -----------------------------------------------------------------------------
plot(bb, vertex.label = NA)

## -----------------------------------------------------------------------------
W <- matrix(c(0,10,10,10,10,75,0,0,0,0,
              10,0,1,1,1,0,0,0,0,0,
              10,1,0,1,1,0,0,0,0,0,
              10,1,1,0,1,0,0,0,0,0,
              10,1,1,1,0,0,0,0,0,0,
              75,0,0,0,0,0,100,100,100,100,
              0,0,0,0,0,100,0,10,10,10,
              0,0,0,0,0,100,10,0,10,10,
              0,0,0,0,0,100,10,10,0,10,
              0,0,0,0,0,100,10,10,10,0),10)

## -----------------------------------------------------------------------------
weighted <- igraph::graph_from_adjacency_matrix(W, mode = "undirected", weighted = TRUE, diag = FALSE)
plot(weighted, edge.width = sqrt(igraph::E(weighted)$weight), vertex.label = NA)

## -----------------------------------------------------------------------------
bb <- global(W, upper = 0, class = "igraph")
plot(bb, vertex.label = NA)

## -----------------------------------------------------------------------------
bb <- global(W, upper = function(x)mean(x), class = "igraph")
plot(bb, vertex.label = NA)

## -----------------------------------------------------------------------------
bb <- disparity(W, alpha = 0.05, narrative = TRUE, class = "igraph")
plot(bb, vertex.label = NA)

## -----------------------------------------------------------------------------
U.with.communities <- igraph::sbm.game(60, matrix(c(.75,.25,.25,.25,.75,.25,.25,.25,.75),3,3), c(20,20,20))
plot(U.with.communities, vertex.label = NA)

## -----------------------------------------------------------------------------
bb <- sparsify.with.lspar(U.with.communities, s = 0.6, narrative = TRUE)
plot(bb, vertex.label = NA)

## -----------------------------------------------------------------------------
U.with.hubs <- igraph::as.undirected(igraph::sample_pa(60, m = 3), mode = "collapse")
plot(U.with.hubs, vertex.size = igraph::degree(bb), vertex.label = NA) #A hairball

## -----------------------------------------------------------------------------
bb <- sparsify.with.localdegree(U.with.hubs, s = 0.3, narrative = TRUE)
plot(bb, vertex.size = igraph::degree(bb), vertex.label = NA)

## -----------------------------------------------------------------------------
B <- matrix(rbinom(100*1000,1,0.5),100,1000)
fdsm.trials(B, riskyp = .75, fwer = FALSE)

## -----------------------------------------------------------------------------
fdsm.trials(B, riskyp = .75, fwer = TRUE)

## -----------------------------------------------------------------------------
B <- rbind(cbind(matrix(rbinom(250,1,.8),10),
                 matrix(rbinom(250,1,.2),10),
                 matrix(rbinom(250,1,.2),10)),
           cbind(matrix(rbinom(250,1,.2),10),
                 matrix(rbinom(250,1,.8),10),
                 matrix(rbinom(250,1,.2),10)),
           cbind(matrix(rbinom(250,1,.2),10),
                 matrix(rbinom(250,1,.2),10),
                 matrix(rbinom(250,1,.8),10)))
bb.object <- fdsm(B)  #Backbone object containing edgewise p-values

## -----------------------------------------------------------------------------
bb1 <- backbone.extract(bb.object, alpha = 0.5, class = "igraph")  #Backbone extracted at alpha = 0.05
plot(bb1)

## -----------------------------------------------------------------------------
bb2 <- backbone.extract(bb.object, alpha = 0.05, class = "igraph")  #Backbone extracted at alpha = 0.05
plot(bb2)

## -----------------------------------------------------------------------------
B <- bipartite.from.probability(R = 5, C = 10, P = .5)
B

## -----------------------------------------------------------------------------
B <- bipartite.from.sequence(R = c(1,1,2), C = c(1,1,2))
B

## -----------------------------------------------------------------------------
par(mfrow=c(1,2), pty="s", mar=c(1,1,1,1))
B <- bipartite.from.distribution(R = 100, C = 100, P = 0.1, rowdist = c(1,10), coldist = c(1,10))
hist(rowSums(B), main = "Agent Deg.", xlab = "", yaxt='n')
hist(colSums(B), main = "Artifact Deg.", xlab = "", yaxt='n')

## -----------------------------------------------------------------------------
B <- bipartite.from.probability(R = 10, C = 10, P = .4)
B <- bipartite.add.blocks(B, blocks = 2, density = .75)
B  #Contains blocks

## -----------------------------------------------------------------------------
mat <- rbind(c(1,0,0), c(0,0,1), c(0,1,1))
mat
fastball(mat)

## -----------------------------------------------------------------------------
mat <- rbind(c(1,0,0), c(0,0,1), c(0,1,1))
bicm(mat)

