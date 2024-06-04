ID_AUC <- function(marker, Stime, status, predict.time, entry = NULL,...){
        if (length(entry) == 0) {
                entry = rep(0, NROW(Stime))
        }
        at.risk <- ((Stime >= predict.time) & (entry <= predict.time))
        eta     <- marker[at.risk]
        status2 <- status
        status2[Stime > predict.time] <- 0
        status2 <- status2[at.risk]
        
        C_num <- 0
        n_case_t    <- sum(status2)
        n_control_t <- sum(1-status2)
        
        if(n_case_t == 0 | n_control_t == 0) return(0.5)
        
        inx_ti <- which(status2 == 1)
        inx_tj <- which(status2 == 0)
        for(id in inx_ti){
                C_num <- C_num + sum(eta[id] < eta[inx_tj]) + 0.5*sum(eta[id] == eta[inx_tj])
        }
        out <- C_num/(n_case_t*n_control_t)
        
        out
        
}

PE <- function(marker, Stime, status, predict.time, SCfun, start.time){
        PE <- 0
        N <- length(marker)
        for(i in 1:N){
               PE <- PE +
                       (status[i]*(Stime[i] <= predict.time)*marker[i]^2)/(SCfun(Stime[i]-0.0000001)/SCfun(start.time)) +
                       ((1-status[i])*(Stime[i] >= predict.time)*(1-marker[i])^2)/(SCfun(predict.time)/SCfun(start.time))
                       
        }
        
        PE/N
}

predict_lm_cumulative_hazard <- function(fit, s, tpred){
        inx_s   <- which(fit$family$data$tr.strat == s)
        tind    <- fit$family$data$tr[inx_s]
        Lambda0 <- fit$family$data$h[inx_s]
        
        est <- data.frame(tind=rev(tind), Lambda0=rev(Lambda0))
        est_fn <- stepfun(est$tind, c(0,est$Lambda0))
        
        est_fn(tpred)
}


getSurvT <- function(survfunction, tind){
        U <- runif(1)
        inx <- which(survfunction <= U)
        ifelse(length(inx) > 0,tind[min(inx)],max(tind)+1)
}


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


gen_gamma_mat <- function(gammafun, S, Tind, ...){
        ret <- matrix(NA, nrow=length(S),ncol=length(Tind))
        for(i in seq_along(S)){
                tinx <- which(Tind <= S[i])
                ret[i,tinx] <- gammafun(Tind[tinx],S[i],...)
        }
        ret
}
