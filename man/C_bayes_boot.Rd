\name{C_bayes_boot}
\alias{C_bayes_boot}

\title{Native C Interface Wrapper for the Function \code{bayes_boot}}

\description{
  These compiled symbols represent the direct entry points to the
  underlying C source code, masked behind the \code{C_} namespace fix.

  The function \code{bayes_boot} performs the Bayesian bootstrap algorithm.
}

% \usage{
%   .C(C_bayes_boot, th, B, s, nn, theta)
% }
%
% \arguments{
%   \item{th}{A double vector that will receive the results.}
%   \item{B}{An integer number of bootstrap replicates.}
%   \item{s}{A numerical vector with the observed data values.}
%   \item{nn}{An integer corresponding to the length of \code{s}.}
%   \item{theta}{A numeric threshold in \eqn{(0, 1)}{(0, 1)}.}
% }
\keyword{internal}

