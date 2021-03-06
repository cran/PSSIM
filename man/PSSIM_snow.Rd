% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{PSSIM_snow}
\alias{PSSIM_snow}
\title{Image structural similarity measure PSSIM based on hypothesis test}
\usage{
PSSIM_snow(
  A,
  A1,
  nprocess = min(8, parallel::detectCores()),
  b = 64,
  a = 2,
  vs = 32,
  wavecoeff = FALSE,
  cs = 2,
  dyn = FALSE
)
}
\arguments{
\item{A}{a grayscale image stored as a matrix.}

\item{A1}{grayscale image stored as a matix. Same dimension as A.}

\item{nprocess}{number of cores (workers) to use for parallel computation.
Note:
 In personal computer, nprocess =detectCores() is good to use.
 On cluster machine, nprocess need to be specified to a number that is
 no more than its number of cores (for courtesy)}

\item{b}{Number of columns in each block. Suggest to use default value 64.}

\item{a}{Number of rows in each block. Suggest to use default value 2.}

\item{vs}{Block shift size. Suggest to use default value 32.}

\item{wavecoeff}{logical of whether the input matrices are wavelet coefficients.
Currently, wavelet version is not implemented.
This parameter is a placeholder for future implementation.}

\item{cs}{dividing factor to split index.}

\item{dyn}{logical, whether dynamic scheduling should be used.}
}
\value{
: Image structural similarity based on PSSIM. The value is in [0,1]
         with values close to 0 meaning the two images are different
         and values close to 1 meaning the two iamges are similar.
}
\description{
PSSIM_snow computes image structural similarity PSSIM of Wang, Maldonado and Silwal (2011) using
            parallel programming.
}
\examples{
  A=miniimagematrix$A
  B=miniimagematrix$B
  # see it with image(A, axes=FALSE, col  = gray((0:255)/256) )
  PSSIM_snow(A, B, nprocess=2)

}
\references{
Haiyan Wang, Diego Maldonado, and Sharad Silwal  (2011).   A Nonparametric-Test-Based
Structural Similarity Measure for Digital Images.  Computational Statistics and Data Analysis.  55: 2925-2936. Doi:10.1016/j.csda.2011.04.021
}
