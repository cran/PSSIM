% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{trt_position}
\alias{trt_position}
\title{Starting and ending position in a vector}
\usage{
trt_position(i1, n)
}
\arguments{
\item{i1}{an integer between 1 and length(n).}

\item{n}{the vector of sample sizes.}
}
\description{
Function trt_position() gives the starting and ending
index of covariate values in the i1th group if all the
covariate values from all treatment groups are together in a vector.
E.g., covariate values in group 1 start from 1st value to the
n1 th value; those in group 2 start from n1+1 and end at (n1+n2)th
value.
This function is for retrieving the position of an observations
when the covariate values from all treatments are stored together
in one vector.
}
\examples{
i = 2; n=c(7, 8); trt_position(i, n)

}
