\name{C_wbayeswise}
\alias{C_wbayeswise}

\title{Native C Interface Wrapper for the Function \code{wbayeswise}}

\description{
  These compiled symbols represent the direct entry points to the
  underlying C source code, masked behind the \code{C_} namespace fix.

  The function \code{wbayeswise} performs a weighted cellwise anomaly detection for stratum‑level wide datasets using Bayesian testing.
}

%\usage{
%.C(C_wbayeswise, s, z, h, r, t, Xc, Xp, dimX, epochs)
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
%  \item{dimX}{A vector with the dimensions of \code{Xc} and \code{Xp}.}
%  \item{epochs}{Pointer to the number of gradient‑descent epochs.}
%  \item{wt}{A vector of record-level anomaly scores (used as weights).}
%}
%
\keyword{internal}

