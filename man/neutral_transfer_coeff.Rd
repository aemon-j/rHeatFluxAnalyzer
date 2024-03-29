% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/neutral_transfer_coeff.R
\name{neutral_transfer_coeff}
\alias{neutral_transfer_coeff}
\title{Calculate the neutral transfer coefficient}
\usage{
neutral_transfer_coeff(Uz, hu)
}
\arguments{
\item{Uz}{vector of wind speed in (m/s)}

\item{hu}{height of the wind speed measurement in m}
}
\value{
A data.frame containing ustar, u10, C_D10N, C_E10N, C_H10N, and C_DN
}
\description{
Calculate the neutral transfer coefficient
}
\examples{
neutral_transfer_coeff(6,10)
}
