% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clearSkySW.R
\name{clearSkySW}
\alias{clearSkySW}
\title{Calculate the clear sky short wave radiation}
\usage{
clearSkySW(time, lat, pressure = NULL, temperature = NULL, RH = NULL)
}
\arguments{
\item{time}{Vector of times as POSIXct}

\item{lat}{latitude for which to calculate clear sky radiation}

\item{pressure}{Air pressure in (hPa) if available}

\item{temperature}{Air temperature in (degrees C) if available}

\item{RH}{relative humidity in percentage if available}
}
\value{
A numerical vecot of clear sky short wave radiation
}
\description{
Calculate the clear sky short wave radiation
}
\references{
Crawford and Duchon, 1991
}
