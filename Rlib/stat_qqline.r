require(proto)

stat_qqline <- function (mapping = NULL, data = NULL, geom = "abline", position = "identity", 
distribution = qnorm, dparams = list(), na.rm = FALSE, ...) { 
  StatQqline$new(mapping = mapping, data = data, geom = geom, position = position, 
  distribution = distribution, dparams = dparams, na.rm = na.rm, ...)
}
  
StatQqline <- proto(ggplot2:::Stat, {
  objname <- "qqline"

  default_geom <- function(.) GeomAbline
  default_aes <- function(.) aes(slope = ..slope.., intercept = ..int..)
  required_aes <- c("sample")

  calculate <- function(., data, scales, quantiles = NULL, distribution = qnorm, dparams = list(), na.rm = FALSE) {
    data <- remove_missing(data, na.rm, "sample", name = "stat_qqline")    
    
    y <- quantile(data$sample, c(0.25, 0.75))
    x <- distribution(c(0.25, 0.75))
    slope <- diff(y)/diff(x)
    int <- y[1L] - slope * x[1L]
    data.frame(slope = slope, int = int)
  }
  
})