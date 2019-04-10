library(ggplot2)

## First setup the model for the simulated data
## Here we assume a simple linear regression model of the form Y = 0.5 + 2*X + epsilon where
## where epsilon is a standard normal r.v.
X <- seq(from = 0, to = 1, by = 0.02)
n <- length(X)
EY <- 0.5 + 2*X

nmc <- 1000
naive.error <- WH.error <- numeric(nmc)
for(i in 1:nmc){

Y <- EY + rnorm(n)

## We next obtain the estimates Yhat and the MSE
mod1 <- lm(Y ~ X)
Yhat <- mod1$fitted
mse <- summary(mod1)$sigma^2
txx <- sum((X - mean(X))^2)



## We first compute the individual 90% confidence intervals
alpha <- 0.1
s2.Yhat <- mse*(1/n + (X - mean(X))^2/txx)
upper.Yhat <- Yhat + sqrt(s2.Yhat)* qt( 1 - alpha/2, df = n -2)
lower.Yhat <- Yhat - sqrt(s2.Yhat)* qt( 1 - alpha/2, df = n -2)

## We next compute the 90% Working-Hotelling confidence bands
W <- sqrt(2*qf(1 - alpha, df1 = 2, df2 = n - 2))
upper.Yhat.WH <- Yhat + sqrt(s2.Yhat)*W
lower.Yhat.WH <- Yhat - sqrt(s2.Yhat)*W

## Compute the number of times the stated confidence intervals for X in [0,1] (in increments of 0.02)
## doesn't contain the mean
naive.error[i] <- sum((EY > upper.Yhat) | (EY < lower.Yhat))
WH.error[i] <- sum((EY > upper.Yhat.WH) | (EY < lower.Yhat.WH))

# Now plot the results
if(nmc == 1){
df.new <- data.frame(Yhat, upper.Yhat, lower.Yhat, upper.Yhat.WH,
                     lower.Yhat.WH, Y, EY)
ggplot(df.new, aes(x = X, y = Y)) + geom_point(alpha = 0.1) + theme_bw() +
  geom_line(aes(x = X, y = upper.Yhat), linetype = "dashed", color = "blue") +
  geom_line(aes(x = X, y = lower.Yhat), linetype = "dashed", color = "blue") +
  geom_line(aes(x = X, y = Yhat), linetype="dashed", colour = "black") +
  geom_line(aes(x = X, y = upper.Yhat.WH), linetype= "dashed", color = "green") +
  geom_line(aes(x = X, y = lower.Yhat.WH), linetype= "dashed", color = "green") +
  geom_line(aes(x = X, y = EY), color = "red") +
  ggtitle(paste("naive = ",naive.error[i], "; WH =", WH.error[i]))
}
}

table(naive.error)
table(WH.error)

mean(naive.error)

nmc <- 10000
coverage1 <- coverage2 <- numeric(nmc)
for(i in 1:nmc){
theta <- rnorm(20, sd = 0.04)
z <- rnorm(20, mean = theta, sd = 1)
CIz.lower <- z - qnorm(0.95, lower.tail = TRUE)
CIz.upper <- z + qnorm(0.95, lower.tail = TRUE)
theta.in.interval <- (theta <= CIz.upper) & (theta >= CIz.lower)
zero.not.in.interval <- (CIz.upper < 0) | (CIz.lower > 0)

coverage1[i] <- sum(theta.in.interval)/20
coverage2[i] <- sum(theta.in.interval & zero.not.in.interval)/sum(zero.not.in.interval)
}

mean(coverage1)
mean(coverage2,na.rm =TRUE)


procrustes <- function(X,Y, type = "I"){
  if(type == "C"){
    X <- X/norm(X, type = "F")*sqrt(nrow(X))
    Y <- Y/norm(Y, type = "F")*sqrt(nrow(Y))
  }
  if(type == "D"){
    tX <- rowSums(X^2)
    tX[tX <= 1e-15] <- 1
    tY <- rowSums(Y^2)
    tY[tY <= 1e-15] <- 1
    X <- X/sqrt(tX)
    Y <- Y/sqrt(tY)
  }

  tmp <- t(X) %*% Y
  tmp.svd <- svd(tmp)
  W <- tmp.svd$u %*% t(tmp.svd$v)
  return(list(error = norm(X%*%W - Y, type = "F"), W = W))
}

ase <- function(A,d) {
  A.svd <- irlba(A,d)
  Xhat <- A.svd$u %*% diag(sqrt(A.svd$d))
  return(Matrix(Xhat))
}

full.ase <- function(A, d, vec.in) {
  A.svd <- irlba(A,d)
  Xhat <- A.svd$u %*% diag(sqrt(A.svd$d))
  return(Matrix(Xhat))
}

## first n.in is insample! FIX THIS!
## => vec.in: index vector for "in" samples!
oos.ase <- function(A11, d, A21, cpp=FALSE)
{
  ## Xhat.insample <- ase(A11,d)

  A11.svd <- irlba(A11,d)
  Xhat.insample <- A11.svd$u %*% diag(sqrt(A11.svd$d))
  Xhat.oos <- A21 %*% A11.svd$u%*%diag(1/sqrt(A11.svd$d))

  # t.in <- t(Xhat.insample)
  #
  # if (cpp) {
  #   tmp1 <- solve(eigenMapMatMult(as.matrix(t.in), as.matrix(Xhat.insample)))
  #   tmp2 <- eigenMapMatMult(as.matrix(tmp1, t.in))
  #   Xhat.oos <- t(eigenMapMatMult(tmp2, as.matrix(A[vec.in,-vec.in])))
  # } else {
  #   tmp1 <- solve(crossprod(Xhat.insample))
  #   tmp2 <- tmp1 %*% t.in
  #   Xhat.oos <- t(tmp2 %*% A12)
  # }

  return(list(Xhat.insample=Xhat.insample, Xhat.oos=Xhat.oos))
}

bench <- function(n=1000, m=200, d=2, nmc=10, cpp=FALSE)
{
  set.seed(123)

  vec.in <- 1:m
  Tmat <- foreach (mc = 1:nmc, .combine='rbind') %dopar% {
    X <- matrix(runif(d*n, max = 1/sqrt(d)), nrow = d)
    g <- sample_dot_product(X, directed=FALSE)

    A <- g[]
    Tn <- system.time(Xhat <- full.ase(A, d))[3]
    Sm <- system.time(Xhat.oos <- oos.ase(A[vec.in,vec.in], d, A[-vec.in,vec.in], cpp))[3]
    save(Xhat, Xhat.oos, file=paste0("oosout-rdpg-n",n,"-m",m,"-d",d,"-mc",mc,".Rbin"))
    Err.is1 <- procrustes(t(X)[vec.in,], Xhat[vec.in,])$error
    Err.is2 <- procrustes(t(X)[vec.in,], Xhat.oos$Xhat.insample)$error
    Err.oos1 <- procrustes(t(X)[-vec.in,], Xhat[-vec.in,])$error
    Err.oos2 <- procrustes(t(X)[-vec.in,], Xhat.oos$Xhat.oos)$error
    #       Err.full1 <- procrustes(t(X), Xhat)$error
    #       Err.full2 <- procrustes(t(X), rbind(Xhat.oos$Xhat.in,Xhat.oos$Xhat.oos))$error
    cat("mc = ", mc, ", Tn = ", Tn, ", Sm = ", Sm, ", err.is1 = ", Err.is1, ", err.is2 = ", Err.is2, ", Err.oos1 = ", Err.oos1, ", Err.oos2 = ", Err.oos2, "\n")
    c(Tn, Sm, Err.is1, Err.is2, Err.oos1, Err.oos2)#, Err.full1, Err.full2)
  }
  colnames(Tmat) <- c("Tn","Tm+S(n-m)","Err.is1","Err.is2","Err.oos1","Err.oos2")
  return(Tmat)
}

plotBench <- function(out)
{
  out.t <- out[,c("Tn","Tm+S(n-m)")]
  df.t <- gather(as.data.frame(out.t), key="method", value="time")

  pt <- ggplot(df.t, aes(x=method, y=time, fill=method)) +
    geom_violin(trim=FALSE, alpha=0.3) +
    #           geom_boxplot(notch=TRUE, width=0.1, alpha=0.5) +
    scale_y_log10() +
    labs(y="time (seconds)") + theme(legend.position = "none")
  print(pt)

  # require(ggridges)
  #   ggplot(df.t, aes(x=time, y=factor(method), fill=..x..)) +
  # #     geom_density_ridges(alpha=0.5) + theme(legend.position = 'none')
  #       scale_y_discrete(limits = rev(levels(df.t$method))) +
  # #     labs(y="scale", x="ARI(ASE,LSE)", fill="ARI") +
  #       geom_density_ridges_gradient(alpha=0.5) +
  #       scale_fill_gradient(low="orange", high="navy") +
  #       theme(axis.title.y = element_blank())

  out.e <- as.tibble(out[,c("Err.is1","Err.is2","Err.oos1","Err.oos2")]);
  #   names(out.e) <- c("alt","null")
  #   cv <- nth(sort(out.e$null, decreasing = TRUE), nrow(out.e)*0.05)
  #   pval <- sum(out.e$alt > cv)/nrow(out.e)
  #   df.e <- gather(out.e, key="type", value="error")

  df.e <- melt(out.e, variable.name="etype", value.name="error")

  pe1 <- ggplot(subset(df.e, etype %in% c("Err.is1","Err.is2")),
                #           aes(x=error, color=etype, fill=etype)) + geom_density(alpha=0.5) +
                aes(x=error, y=..scaled.., color=etype, fill=etype)) + geom_density(alpha=0.5) +
    labs(x="|| Xhat - X ||_F") +
    scale_fill_hue(labels = c("Xhat.ase - X", "Xhat.oos - X")) +
    guides(color=FALSE, fill=guide_legend("Error_in")); pe1
  #       geom_vline(xintercept=cv, col="black", linetype="dashed") +
  #       geom_text(aes(label=paste0("pval = ", pval), x=cv-15, y=0.5), hjust=0)
  #       annotate("text", x=cv-15, y=0.1, label=paste0("pval = ", pval), angle=90)
  #       geom_text(aes(x=cv-15, y=0.5, label=paste0("pval = ", sum(out.e$alt>cv))), color="black")
  pe2 <- ggplot(subset(df.e, etype %in% c("Err.oos1","Err.oos2")),
                aes(x=error, y=..scaled.., color=etype, fill=etype)) + geom_density(alpha=0.5) +
    #            aes(x=error,  color=etype, fill=etype)) + geom_density(alpha=0.5) +
    labs(x="|| Xhat - X ||_F") +
    scale_fill_hue(labels = c("Xhat.ase - X", "Xhat.oos - X")) +
    guides(color=FALSE, fill=guide_legend("Error_out")); pe2

  ggarrange(pe1, pe2, heights=c(1,1))
}

library("igraph")
library("irlba")
nmc <- 10
err <- matrix(0,nmc,3)
for(i in 1:nmc){
B <- matrix(c(0.5,0.3,0.3,0.5),nrow = 2)
n <- 16000

A <- sbm.game(n, n^{-3/4}*B, n*c(0.5,0.5))
A.svd <- irlba(A[],2)
Xhat <- A.svd$u %*% diag(sqrt(A.svd$d))
idx1 <- 1:(n*0.5)
Xhat1 <- colMeans(Xhat[idx1,])
Xhat2 <- colMeans(Xhat[-idx1,])
err[i,1] <- sum(Xhat1*Xhat1) - n^{-3/4}*B[1,1]
err[i,2] <- sum(Xhat1*Xhat2) - n^{-3/4}*B[1,2]
err[i,3] <- sum(Xhat2*Xhat2) - n^{-3/4}*B[2,2]
}

