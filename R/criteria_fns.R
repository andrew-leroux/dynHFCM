#' Function for calculating "local" AUC(s,t)
#' @param marker vector of predicted survival probability at the time indicated by the predict.time argument
#' @param Stime vector of observed survival times for participants
#' @param status vector of event indicators for participants
#' @param predict.time scalar indicating the time (t) AUC(s,t) should be calcuated at
#' @param entry exclude participants with events/censoring before this time
#'
#' @export
ID_AUC <- function(marker, Stime, status, predict.time, entry = NULL){
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

#'  Function for calculating "local" PE(s,t)
#' @param marker vector of predicted survival probability at the time indicated by the predict.time argument
#' @param Stime vector of observed survival times for participants
#' @param status vector of event indicators for participants
#' @param predict.time scalar indicating the time (t) AUC(s,t) should be calcuated at
#' @param SCfun for evaluating the survival function of censoring distribution
#' @param start.time landmark time associated with PE(s,t)
#'
#' @export
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
