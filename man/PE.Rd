% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/criteria_fns.R
\name{PE}
\alias{PE}
\title{Function for calculating "local" PE(s,t)}
\usage{
PE(marker, Stime, status, predict.time, SCfun, start.time)
}
\arguments{
\item{marker}{vector of predicted survival probability at the time indicated by the predict.time argument}

\item{Stime}{vector of observed survival times for participants}

\item{status}{vector of event indicators for participants}

\item{predict.time}{scalar indicating the time (t) AUC(s,t) should be calcuated at}

\item{SCfun}{for evaluating the survival function of censoring distribution}

\item{start.time}{landmark time associated with PE(s,t)}
}
\description{
Function for calculating "local" PE(s,t)
}
