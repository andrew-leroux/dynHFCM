## Main Function: make_lm_data()
## Inputs: 
##  * Reqired
##    - data frame of data (in long format)
##    - vars_tv: character vector of time-varying predictors 
##    - vars_tf: character vector of time-fixed predictors (must not vary within grouping variable)
##    - id: character vector of length 1 indicating grouping ("id") variable
##    - event_time: character vector of length 1 indicating event time variable ("event_time")
##    - time: character vector of length 1 indicating time of observation of the time-varying predictors ("time").
##      for the initial implementation assume a common grid of observation times if there are multiple time-varying predictors!
##    - S: numeric vector of landmark times (S)
##    - w: numeric vector of prediction windows (w). Must be either length 1 or same length as the vector of landmark times (s)
##
## Function steps:
##    1) Input checking 
##       Make sure names in data frame supplied to the data argument do not conflict with any names generated in the function
##    2) Insert pseudo observations for each group at each landmark time
##        a. For each unit of the grouping variable create an index varibale j = 1,...,J denoting first, second, ..., Jth observation
##           of the time-varying variable(s) 
##        b. Insert pseudo observations at each landmark time, using only longitudinal data available AT OR PRIOR to the landmark time.
##    3) For each group, for each landmark time separately, convert time-varying predictors AND time of observation to wide format matrices.
##       At the end of this step, each row of the data should correspond to 1 row per subject-landmark time. 
##        a. separate out the time-varying data
##        b. reshape (tidyr::spread?) the time varying predictors separately
##        c. merge the time-varying data back with time-invariant data
##    4) Fill NA's in the time-varying predictor/time of observation matrices according to (user supplied?) preferences
make_lm_data <- function(data, vars_tv, vars_tf, id, event_time, time, S, 
                         wide=FALSE, censor_var, w=NULL){
        ## does not accept character time varying variables!
        nms_sp <- c("J", "s","S")
        if(any(c("J", "s") %in% names(data))){
                stop(paste0(paste0(nms_sp, collapse=", "), " are reserved column names. Please re-name the columns of your input data matrix accordingly."))
        }
        if(is.null(s)){
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
        data <- data %>% group_by(.data[[id]]) %>% dplyr::mutate(J = 1:n()) %>% ungroup() %>% filter(.data[[time]] < .data[[event_time]] & .data[[time]] <= max(S))
                
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
        
        
        data_lm
}







make_dyn_pred_data <- function(data_lm, vars_tv, id, time, argvals, sm_fns = NULL){
        vars_tv <- unique(c(time,vars_tv))
        N <- nrow(data_lm)
        t_new <- c()
        for(i in 1:N){
                t_new[[i]] <- sort(unique(na.omit(c(data_lm[[time]][i,], argvals))))
        }
        tmax <- max(vapply(t_new,length,numeric(1)))
        
        for(p in vars_tv){
                data_lm[[paste0(p,"_dyn")]] <- I(matrix(NA, ncol=tmax, nrow=N))
        }
        for(i in 1:N){
                t_i       <- na.omit(data_lm[[time]][i,])
                inx_col_i <- na.omit(match(t_i, t_new[[i]]))
                for(p in vars_tv){
                        if(p == time){
                                data_lm[[paste0(p,"_dyn")]][i,1:length(t_new[[i]])] <- t_new[[i]]        
                        } else {
                                data_lm[[paste0(p,"_dyn")]][i,inx_col_i] <-   data_lm[[p]][i,1:length(t_i)]     
                        }
                }
        }
        Nt_dyn <- ncol(data_lm[[paste0(vars_tv[1],"_dyn")]])
        data_lm$smat_dyn <- I(data_lm$smat[,1] %o% rep(1,Nt_dyn))
        
        if(!is.null(sm_fns)){
                vars_tv_d   <- paste0(vars_tv, "_dyn")
                vars_tv_p <- paste0(vars_tv, "_dyn_pred")
                for(p in seq_along(vars_tv_d)){
                        if(vars_tv[p] == time) next
                        data_lm[[vars_tv_p[p]]] <- I(matrix(NA, ncol=ncol(data_lm[[vars_tv_d[p] ]]), nrow=N))       
                }
                for(i in 1:N){
                        df_i <- data.frame("argvals"=as.vector(data_lm_pred[[paste0(time, "_dyn")]][i,]))
                        for(p in vars_tv_d){
                                if(p == paste0(time,"_dyn")) next
                                df_i[[p]] <-    as.vector(data_lm_pred[[p]][i,])     
                        }
                        df_i <- df_i[!rowSums(!is.na(df_i)) == 0,]
                        df_i <- df_i[df_i$argvals <= data_lm_pred$s[i],]
                        
                        
                        for(p in seq_along(vars_tv_p)){
                                if(vars_tv[p] == time) next
                                fn_p <-  sm_fns[[vars_tv[p]]]
                                df_i[[vars_tv_p[p]]] <- fn_p(df_i, y_nm=vars_tv_d[p], argvals_nm="argvals")
                                data_lm[[vars_tv_p[p]]][i,1:nrow(df_i)] <- df_i[[vars_tv_p[p]]]
                        }
                      
                if(i %% 100 == 0) print(i)
                }
                
        }
        
        data_lm
}

