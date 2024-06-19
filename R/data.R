#' Simulated data under a historical functional Cox data generating mechanism
#'
#' A dataset containing the survival outcomes (survival or censoring time, event indicator)
#' and longitudinal data for a time-varying covariate for 1000 individuals. Data are in long format and thus survival times and event
#' indicators are repeated within individuals.
#'
#' @format A data frame with 42619 rows and 5 variables:
#' \describe{
#'   \item{id}{subject identifier}
#'   \item{T}{event or censoring time}
#'   \item{E}{event indicator (1 if event observed, 0 if event censored)}
#'   \item{Z}{longitidunal/functional predictor}
#'   \item{tind}{observation time for the functional predictor, Z}
#' }
#' @source simulated data
"hfcm_sim_data"
