\name{C_post_results}
\alias{C_post_results}

\title{Native C Interface Wrapper for the Function \code{post_results}}

\description{
  These compiled symbols represent the direct entry points to the
  underlying C source code, masked behind the \code{C_} namespace fix.

  The function \code{post_results} does posterior computations over the whole data set.
}

%\usage{
%.C(C_post_results, H_err, R_err, T_err, prior_mat, G, dimX)
%}
%
%\arguments{
%  \item{H_err}{Matrix of errors/residuals from historical analyses.}
%  \item{R_err}{Matrix of errors/residuals from relational analyses.}
%  \item{T_err}{Matrix of errors/residuals from distribution‑tail analyses.}
%  \item{prior_mat}{Matrix of prior probabilities for the regular cases.}
%  \item{G}{Empty matrix of integers where to store the results.}
%  \item{dimX}{Vector of dimensions of the matrices listed above.}
%}
%
\keyword{internal}

