\name{C_wcellwise}
\alias{C_wcellwise}

\title{Native C Interface Wrapper for the Function \code{wcellwise}}

\description{
  These compiled symbols represent the direct entry points to the
  underlying C source code, masked behind the \code{C_} namespace fix.

  The function \code{wcellwise} performs a weighted cellwise anomaly detection for stratum‑level wide datasets.
}

%\usage{
%.C(C_wcellwise, s, z, h, r, t, Xc, Xp, dimX)
%}
%
%\arguments{
%  \item{s}{An empty vector where to store the final scores.}
%  \item{z}{An empty vector for data‑format anomaly scores.}
%  \item{h}{An empty vector for historical anomaly scores.}
%  \item{r}{An empty vector for relational anomaly scores.}
%  \item{t}{An empty vector for tail anomaly scores.}
%  \item{Xc}{A matrix of current data.}
%  \item{Xp}{A matrix of previously reported data.}
%  \item{dimX}{A vector with the dimensions of `Xc` and `Xp`.}
%  \item{wt}{A vector of record-level anomaly scores (used as weights).}
%}
%
\keyword{internal}

