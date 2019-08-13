#' bmmmcmc.
#' 
#' Package to sample from Bayesian Bernoulli Mixture Models
#'
#' @name bmmmcmc
#' @useDynLib bmmmcmc
#' @docType package
NULL

#' Simulated mixture dataset.
#'
#' Has 100 observations of 5 Bernoulli variables from 2 clusters.
#' The clusters are split in the ratio 0.7:0.3.
#' The true theta Bernoulli parameters are:
#' 
#'  [[0.7, 0.8, 0.2, 0.1, 0.1],
#'   [0.2, 0.2, 0.9, 0.8, 0.6]]
#'
#' @docType data
#' @name K2_N100_P5
#' @format A matrix with 100 rows and 5 columns
NULL

#' Simulated mixture dataset.
#'
#' Has 1000 observations of 5 Bernoulli variables from 2 clusters.
#' 
#' It is exactly the same as K2_N100_P5 but with 1000 
#' observations rather than 100.
#' 
#' The clusters are split in the ratio 0.7:0.3.
#' The true theta Bernoulli parameters are:
#' 
#'  [[0.7, 0.8, 0.2, 0.1, 0.1],
#'   [0.2, 0.2, 0.9, 0.8, 0.6]]
#'
#' @docType data
#' @name K2_N1000_P5
#' @format A matrix with 1000 rows and 5 columns
NULL

#' Simulated mixture dataset.
#'
#' Has 1000 observations of 5 Bernoulli variables from 3 clusters.
#' 
#' The clusters are split in the ratio 0.6:0.2:0.2.
#' The true theta Bernoulli parameters are:
#' [[0.7, 0.8, 0.2, 0.1, 0.1],
#'  [0.3, 0.5, 0.9, 0.8, 0.6],
#'  [0.1, 0.2, 0.5, 0.4, 0.9]]
#' 
#' @docType data
#' @name K3_N1000_P5
#' @format A matrix with 1000 rows and 5 columns
NULL
