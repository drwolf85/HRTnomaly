\name{C_pif}
\alias{C_pif}

\title{Native C Interface Wrapper for the Function \code{pif}}

\description{
  These compiled symbols represent the direct entry points to the
  underlying R code, masked behind the \code{C_} namespace fix.

  The function \code{pif} is an R‑wrapper for performing Proximity Isolation Forest (PIF) computations. It returns a generic S3 object containing the record-level anomaly scores from a PIF.
}

%\usage{
%.Call(C_pif, dta, _prx, _nt, _nss, max_depth, dst_fun, Rnv)
%}
%
%\arguments{
%  \item{dta}{A list of data.}
%  \item{\_prx}{The type of proximity isolation forest.}
%  \item{\_nt}{The total number of trees in the trained forest.}
%  \item{\_nss}{The number of subsamples.}
%  \item{max_depth}{The maximum depth of a tree in the forest.}
%  \item{dst_fun}{A metric function of choice.}
%  \item{Rnv}{An R environment used.}
%}
%
\keyword{internal}

