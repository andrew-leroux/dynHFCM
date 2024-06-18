
#' Simulate survival times using a survival curve evaluated on a discrete grid
#' @param survfunction vector containing estimates of the survival function \eqn{S(t) \in [0,1]}
#' @param tind vector indicating the time points associated with survfunction input
#'
#' @export
getSurvT <- function(survfunction, tind){
  U <- stats::runif(1)
  inx <- which(survfunction <= U)
  ifelse(length(inx) > 0,tind[min(inx)],max(tind)+1)
}


#' Create a historical functional coefficient using an exponential decay model for a fixed landmark time
#' @param t vector containing the historical time for evaluating the historical functional coefficient
#' @param s scalar containing the "landmark" time for evaluating the historical functional coefficient
#' @param d numeric indicating the mean of the exponential random variable whose density function determins the rate of decay for the coefficient function
#' @param scale scaling factor for the coefficient
#'
#' @export
gamma_exp_decay = function (t, s, d = 1/5, scale = 1) {
  u <- (s - t)
  scale*stats::dexp(u, 1/d)
}

#' Create a historical functional coefficient where risk is of the form \eqn{c \times cos(2 \times \pi \times u/p)} with \eqn{u=s-t} for a fixed landmark time
#' @param t vector containing the historical time for evaluating the historical functional coefficient
#' @param s scalar containing the "landmark" time for evaluating the historical functional coefficient
#' @param period period of the cosine function (\eqn{p})
#' @param scale scaling factor for the coefficient (\eqn{c})
#'
#' @export
gamma_cos = function (t, s, period = 1/2, scale=10) {
  u <- (s - t)
  cos(2 * pi * u / period)  * scale
}

#' Function for evaluating the historical functional coefficient. Called by the genSurv function
#' @param gammafun function which evaluates the landmark historical functional coefficient. Should evaluate for a fixed landmark time and a vector of historical times using arguments passed down tot he function via ...
#' @param S vector containing the "landmark" times for evaluating the historical functional coefficient
#' @param period period of the cosine function (\eqn{p})
#' @param scale scaling factor for the coefficient (\eqn{c})
#' @param Tind vector containing the historical time for evaluating the historical functional coefficient
#' @param ... arguments to be passed to gammafun
#'
#' @export
gen_gamma_mat <- function(gammafun, S, Tind, ...){
  ret <- matrix(NA, nrow=length(S),ncol=length(Tind))
  for(i in seq_along(S)){
    tinx <- which(Tind <= S[i])
    ret[i,tinx] <- gammafun(Tind[tinx],S[i],...)
  }
  ret
}


#' Function for simulating survival curves for the historical functional Cox model using a single time-fixed and time-varying covariate
#' @param lambda0 function which takes time as an input and returns the baseline hazard
#' @param beta function which takes time as an input and returns a potentially time varying effect for a single scalar covariate
#' @param gamma function which takes time and historical time as input and returns a historical functional coefficient
#' @param tind grid on which to evaluate the survival curve
#' @param Z matrix of the longitudinal/time varying predictor to be used in creating the survival curve for individuals, evaluated on the grid supplied to the "tind" argument (should have N rows)
#' @param X vector of time-fixed covariate (should be of length N)
#' @param N integer, number of functions observed. Should be equal to the length of X and number of rows of Z
#' @param cens_fun function which generates censoring times
#' @param returnS logical indicating whether to return the simulated survival functions along with survival times
#' @param ... arguments to be passed to gen_gamma_mat
#'
#' @export
genSurv <- function(lambda0, beta, gamma=gamma_exp_decay,
                    tind, Z, X, N, cens_fun=NULL,returnS=FALSE,
                    ...){
  tlen = length(tind)
  stopifnot(N == nrow(Z))
  stopifnot(tlen == ncol(Z))
  stopifnot(all(is.function(lambda0), is.function(beta), is.function(gamma)))
  if(is.null(cens_fun)){
    cens_fun <- function(N, tind){runif(N, min(tind), max(tind))}
  }
  surv <- matrix(NA, nrow=N, ncol=tlen)
  lambda0_t <- lambda0(tind)
  beta_t    <- beta(tind)
  gamma_t     <- gen_gamma_mat(gamma, Tind = tind, S=tind, ...)
  delta_t     <- diff(tind)
  for(i in 1:N){
    ## get log hazard for all individuals at all time points
    ll <- vapply(seq_along(tind[-tlen]), function(x){
      sum(X[i] * beta_t[x]) + sum(Z[i,1:x] * gamma_t[x,1:x] * delta_t[1:x])
    }, numeric(1))
    ## approximate survival function using Riemann integration
    surv[i,] <- c(1,exp(-cumsum(exp(ll + log(lambda0_t)))*delta_t))
  }

  ## get survival and censoring time
  E  <- apply(surv,1,getSurvT,tind=tind)
  C  <- cens_fun(N, tind)
  ET <-  pmin(E,C)

  if(returnS){
    ret <- data.frame("E" = as.numeric(E < C),
                      "T" = ET,
                      "surv" = I(surv))
  } else {
    ret <- data.frame("E" = as.numeric(E < C),
                      "T" = ET)
  }

}
