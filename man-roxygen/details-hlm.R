#' @details The heteroscedastic linear model (HLM) is defined as
#' \deqn{
#' y_i \mid \boldsymbol{x}_i, \boldsymbol{z}_i \stackrel{\mathrm{ind}}{\sim} \mathcal N\big(\boldsymbol{x}_i'\boldsymbol{\beta}, \exp(\boldsymbol{z}_i'\boldsymbol{\gamma})\big),
#' }{
#' y_i | x_i, z_i ~ind N(x_i'\beta, exp(z_i'\gamma)),
#' }
#' where for each subject \eqn{i}, \eqn{y_i} is the response, and \eqn{\boldsymbol{x}_i \in \mathbb{R}^p}{x_i \in R^p} and \eqn{\boldsymbol{z}_i \in \mathbb{R}^q}{z_i \in R^q} are mean and variance covariate vectors, respectively.
