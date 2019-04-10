## Correlated ER
rg.sample.correlated.gnp <- function(P,Q,sigma){
  require(igraph)
  n <-  nrow(P)
  U <- matrix(0, nrow = n, ncol = n)
  U[col(U) > row(U)] <- runif(n*(n-1)/2)
  U <- (U + t(U))
  diag(U) <- runif(n)
  A <- (U < P) + 0 ;
  diag(A) <- 0

  avec <- A[col(A) > row(A)]
  pvec <- P[col(P) > row(P)]
  qvec <- Q[col(Q) > row(Q)]
  bvec <- numeric(n*(n-1)/2)

  uvec <- runif(n*(n-1)/2)

  idx1 <- which(avec == 1)
  idx0 <- which(avec == 0)

  var.vec <- sqrt(pvec*(1-pvec)*qvec*(1 - qvec))

  bvec[idx1] <- (uvec[idx1] < (sigma*var.vec/pvec + qvec)[idx1]) + 0
  bvec[idx0] <- (uvec[idx0] < (qvec - sigma*var.vec/(1-pvec))[idx0]) + 0

  B <- matrix(0, nrow = n, ncol = n)
  B[col(B) > row(B)] <- bvec
  B <- B + t(B)
  diag(B) <- 0

  return(list(A = graph.adjacency(A,"undirected"), B = graph.adjacency(B,"undirected")))
}

# non-igraph version of correlated SBM
#gg <- rg.sample.SBM.correlated(n = 100, B = matrix(c(0.5,0.5,0.2,0.5), nrow = 2), rho = c(0.4,0.6), sigma = 0.2)
#cor(as.vector(gg$adjacency$A[]), as.vector(gg$adjacency$B[]))
rg.sample.SBM.correlated <- function(n, B, rho, sigma, conditional = FALSE){
  if(!conditional){
    tau <- sample(c(1:length(rho)), n, replace = TRUE, prob = rho)
  }
  else{
    tau <- unlist(lapply(1:2,function(k) rep(k, rho[k]*n)))
  }
  P <- B[tau,tau]
  return(list(adjacency=rg.sample.correlated.gnp(P, sigma),tau=tau))
}

