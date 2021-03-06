test_that("tomatrix classes match", {
  library(Matrix)
  library(igraph)
  library(network)

  ### unipartite matrix --> matrix
  M <- matrix(rbinom(5*5,1,.5),5,5)
  rownames(M) <- LETTERS[1:5]
  colnames(M) <- LETTERS[1:5]
  test <- tomatrix(M)
  expect_equal(class(test$G), c("matrix","array"))

  ### bipartite matrix --> matrix
  M <- matrix(rbinom(5*5,1,.5),5,5)
  rownames(M) <- LETTERS[1:5]
  colnames(M) <- letters[1:5]
  test <- tomatrix(M)
  expect_equal(class(test$G), c("matrix","array"))

  ### Matrix --> matrix
  M <- matrix(rbinom(5*5,1,.5),5,5)
  rownames(M) <- LETTERS[1:5]
  colnames(M) <- LETTERS[1:5]
  M <- Matrix::Matrix(M)
  test <- tomatrix(M)
  expect_equal(class(test$G), c("matrix","array"))

  ### sparse matrix --> matrix
  M <- matrix(rbinom(5*5,1,.5),5,5)
  rownames(M) <- LETTERS[1:5]
  colnames(M) <- LETTERS[1:5]
  M <- Matrix::Matrix(M, sparse = TRUE)
  test <- tomatrix(M)
  expect_equal(class(test$G), c("matrix","array"))

  ### edgelist (unweighted bipartite) --> matrix
  M <- data.frame(v1 = c("Senator A", "Senator A", "Senator B"),
                  v2 = c("Bill X", "Bill Y", "Bill X"))
  test <- tomatrix(M)
  expect_equal(class(test$G), c("matrix","array"))

  ### edgelist (weighted bipartite) --> matrix
  M <- data.frame(v1 = c("Senator A", "Senator A", "Senator B"),
                  v2 = c("Bill X", "Bill Y", "Bill X"),
                  v3 = c(1,2,3))
  test <- tomatrix(M)
  expect_equal(class(test$G), c("matrix","array"))

  ### edgelist (unweighted unipartite) --> matrix
  M <- data.frame(v1 = c("Senator A", "Senator A", "Senator C"),
                  v2 = c("Senator B", "Senator C", "Senator A"))
  test <- tomatrix(M)
  expect_equal(class(test$G), c("matrix","array"))

  ### edgelist (weighted unipartite) --> matrix
  M <- data.frame(v1 = c("Senator A", "Senator A", "Senator C"),
                  v2 = c("Senator B", "Senator C", "Senator A"),
                  v3 = c(2,7,1))
  test <- tomatrix(M)
  expect_equal(class(test$G), c("matrix","array"))

  ### igraph (unweighted bipartite) --> matrix
  M <- igraph::sample_bipartite(3,5,p=.5)
  test <- tomatrix(M)
  expect_equal(class(test$G), c("matrix","array"))

  ### igraph (weighted bipartite) --> matrix
  M <- igraph::sample_bipartite(3,5,p=.5)
  igraph::E(M)$weight <- 4
  test <- tomatrix(M)
  expect_equal(class(test$G), c("matrix","array"))

  ### igraph (unweighted unipartite) --> matrix
  M <- igraph::erdos.renyi.game(n = 5, p = .5, directed = TRUE)
  test <- tomatrix(M)
  expect_equal(class(test$G), c("matrix","array"))

  ### igraph (weighted unipartite) --> matrix
  M <- igraph::erdos.renyi.game(n = 5, p = .5, directed = TRUE)
  igraph::E(M)$weight <- 4
  test <- tomatrix(M)
  expect_equal(class(test$G), c("matrix","array"))

  ### statnet (unweighted bipartite) --> matrix
  M <- network::network(matrix(rbinom(5*5,1,.5),5,5), bipartite = TRUE)
  test <- tomatrix(M)
  expect_equal(class(test$G), c("matrix","array"))

  ### statnet (weighted bipartite) --> matrix
  M <- network::network(matrix(rbinom(5*5,1,.5),5,5), bipartite = TRUE)
  network::set.edge.attribute(M,"weight",4)
  test <- tomatrix(M)
  expect_equal(class(test$G), c("matrix","array"))

  ### statnet (unweighted unipartite) --> matrix
  M <- network::network(matrix(rbinom(5*5,1,.5),5,5), bipartite = FALSE)
  test <- tomatrix(M)
  expect_equal(class(test$G), c("matrix","array"))

  ### statnet (weighted bipartite) --> matrix
  M <- network::network(matrix(rbinom(5*5,1,.5),5,5), bipartite = FALSE)
  network::set.edge.attribute(M,"weight",4)
  test <- tomatrix(M)
  expect_equal(class(test$G), c("matrix","array"))

})
