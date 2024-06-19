#' @importFrom dplyr group_by ungroup mutate filter select slice %>% n
#' @importFrom reshape2 melt acast
#' @importFrom stats as.formula stepfun dexp runif
#' @importFrom rlang .data
NULL

#' Function for creating the landmark dataset
#' @param data frame of data (in long format)
#' @param vars_tv character vector of time-varying predictors
#' @param vars_tf character vector of time-fixed predictors (must not vary within grouping variable)
#' @param id character vector of length 1 indicating grouping ("id") variable
#' @param event_time character vector of length 1 indicating event time variable ("event_time")
#' @param time character vector of length 1 indicating time of observation of the time-varying predictors ("time"). Currently only supports a common grid of observation times if there are multiple time-varying predictors.
#' @param S numeric vector of landmark times (S)
#' @param wide logical indicating whether data supplied to the function are in wide format (FALSE) or long format (TRUE)
#' @param censor_var character vector of length 1 indicating event indicator variable. This variable must be coded as 1=event, 0=censor
#' @param w numeric vector of prediction windows (w). Must be either length 1 or same length as the vector of landmark times (s)
#' @param center logical indicating whether to also return landmark mean-centered functional predictors (only works for regularly observed data currently)
#'
#' @export
make_lm_data <- function(data, vars_tv, vars_tf, id, event_time, time, S,
                         wide=FALSE, censor_var, w=NULL){
        ## does not accept character time varying variables
        nms_sp <- c("J", "s","S")
        if(any(c("J", "s") %in% names(data))){
                stop(paste0(paste0(nms_sp, collapse=", "), " are reserved column names. Please re-name the columns of your input data matrix accordingly."))
        }
        if(is.null(S)){
                stop("One or more landmark times must be specified")
        }
        if(is.null(w)){
                warning("No prediction windows specified, assume to be infity")
                w <- Inf
        }
        if(is.null(censor_var)){
                stop("censoring variable must be specified")
        }
        if(!all(data[[censor_var]] %in% c(0,1))) stop("Censoring indicator must be either 0/1")

        ## add in time variable to vector of time-varying covariates
        vars_tv <- unique(c(time,vars_tv))
        ## add in event time to time-fixed variable
        vars_tf <- unique(c(event_time, censor_var, vars_tf))

        ##
        nS <- length(S)

        if(length(w)==1) w = rep(w, nS)

        ## acast() will re-order based on factor level, re-code id variable as a factor
        data[[id]] <- factor(data[[id]], levels=unique(data[[id]]))

        if(wide){
                nt   <- ncol(data[[vars_tv[1]]])
                N    <- nrow(data)
                data <- data.frame(data[c(id,vars_tf)][rep(1:N, each=nt),],
                                   lapply(data[vars_tv], function(x) as.vector(t(x)))
                                   )
        }
        ## create visit number variable
        data <- data %>%
            group_by(.data[[id]]) %>%
            dplyr::mutate(J = 1:n()) %>%
            ungroup() %>%
            filter(.data[[time]] < .data[[event_time]] & .data[[time]] <= max(S))

        ## make long version of the data using reshape2 functions
        ## then make an array of the "wide" data
        data_long    <- melt(data[,c(vars_tv,id,"J")], id=c(id,"J"))
        data_tv_wide <- acast(data_long, as.formula(paste0(id, " ~ J ~ variable")))


        ## function credit: https://stackoverflow.com/questions/14500707/select-along-one-of-n-dimensions-in-array
        index_array <- function(x, dim, value, drop = FALSE) {
                # Create list representing arguments supplied to [
                # bquote() creates an object corresponding to a missing argument
                indices <- rep(list(bquote()), length(dim(x)))
                indices[[dim]] <- value

                # Generate the call to [
                call <- as.call(c(
                        list(as.name("["), quote(x)),
                        indices,
                        list(drop = drop)))

                # Finally, evaluate it
                eval(call)
        }

        # uid <- unique(data[[id]])
        uid <- attributes(data_tv_wide)$dimnames[[1]]
        nid <- length(uid)
        # S_wide       <- rep(S, nid)
        data_tv_wide <- index_array(data_tv_wide, 1, value=rep(1:nrow(data_tv_wide), each=nS))

        # dim_smat <- dim(data_tv_wide)[3]+1
        # data_tv_wide[,,dim_smat] <- matrix(rep(S, each=dim(data_tv_wide)[2]), byrow=TRUE,ncol=dim(data_tv_wide)[2])
        smat <- matrix(rep(S, each=dim(data_tv_wide)[2]), byrow=TRUE,ncol=dim(data_tv_wide)[2])
        smat <- do.call(rbind, rep(list(smat), nid))

        na_mat <- smat < data_tv_wide[,,1]

        ## replace NAs where time of observation greater than landmark time
        for(p in seq_along(vars_tv)){
                data_tv_wide[,,p][na_mat] <- NA
        }

        ## combine time-invariant and time-varying variables into a single landmarked data frame
        data_lm     <- dplyr::select(data, !!id, !!vars_tf)  %>% group_by(.data[[id]]) %>% slice(rep(1, nS)) %>% mutate(s = S) %>% ungroup()
        data_lm$w_s <- rep(w, nid)
        data_lm <- data_lm %>%
                mutate("event_time_lm" = pmin(.data[[event_time]], s + w_s),
                       "event_lm" = .data[[censor_var]]*(event_time_lm >= .data[[event_time]]) )
        ## need to remove tbl and tbl_df classes from the landmarked data frame
        ## otherwise subsetting the the tv-matrices using `[` will fail
        data_lm <- data.frame(data_lm)

        var_nms <- attributes(data_tv_wide)$dimnames[[3]]
        for(p in seq_along(var_nms)){
                data_lm[[var_nms[p]]] <- I(data_tv_wide[,,p])
        }

        data_lm$smat <- I(smat)


        inx_keep <- which(rowSums(!na_mat , na.rm=TRUE) >= 1)
        data_lm <- data_lm[inx_keep,,drop=FALSE]
        data_lm <- filter(data_lm, s < event_time_lm)

          # center the  functional predictor at each landmark time in the training data and
          nT <- ncol(data_lm$Z)
          uS <- sort(unique(data_lm$s[data_lm$id]))
          nS <- length(uS)
          tind_obs <- sort(unique(data[[time]]))

          # divide by a constant corresponding to Riemann integration
          # gap between observations of the functional predictor
          dt <- diff(tind_obs[tind_obs <= max(uS)])
          # small mass we add to be able to evaluate the integral of the historical term at the baseline
          c0 <- mean(dt)

          data_lm$Z_int_cn <-
            I(matrix(NA_real_, ncol = nT, nrow = nrow(data_lm)))
          for(j in seq_along(uS)){
            inx_j  <- which(data_lm$s == uS[j])
            n_j    <- length(inx_j)
            # calculate mean at each landmark time
            inx_j  <- which(data_lm$s == uS[j])
            mu_Z_j       <- matrix(colMeans(data_lm$Z[inx_j,]), nrow=1) %x% matrix(1, nrow=n_j, ncol=1)
            # center and set up the integration matrix
            c_j    <- sum(tind_obs <= uS[j])
            Z_j_cn <- (data_lm$Z[inx_j,] - mu_Z_j) * (matrix(1, nrow=n_j, ncol=1) %x% matrix(c(c0, dt), nrow=1))
            data_lm$Z_int_cn[inx_j,] <- Z_j_cn
          }
        # replace missing values in the time and functional predictor matrices with 0s
        # the choice here is not entirely arbitrary as the model will attempt to estimate the coefficient at these points
        inx_na <- is.na(data_lm$tind)
        data_lm$tind[inx_na] <- data_lm$Z_int_cn[inx_na] <- 0


        data_lm
}


#' Function for creating the landmark dataset
#' @param fit fitter LM-HFLCM
#' @param s length 1 numeric vector indicating landmark time
#' @param tpred numeric vector of predictor times to evaluate the cumulative hazard function
#'
#' @export
predict_lm_cumulative_hazard <- function(fit, s, tpred){
  inx_s   <- which(fit$family$data$tr.strat == s)
  tind    <- fit$family$data$tr[inx_s]
  Lambda0 <- fit$family$data$h[inx_s]

  est <- data.frame(tind=rev(tind), Lambda0=rev(Lambda0))
  est_fn <- stats::stepfun(est$tind, c(0,est$Lambda0))

  est_fn(tpred)
}

