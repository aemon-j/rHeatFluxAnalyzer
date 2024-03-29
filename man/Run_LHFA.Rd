% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Run_rLHFA.R
\name{Run_LHFA}
\alias{Run_LHFA}
\title{run rLakeHeatFluxAnalyzer}
\usage{
Run_LHFA(LakeName, Folder, skipLoad = FALSE)
}
\arguments{
\item{LakeName}{Name of the Lake. All forcing and config files must have this name}

\item{Folder}{Folder the forcing and config files are stored in}

\item{skipLoad}{Use own config file (must then be placed in Folder)}
}
\description{
run rLakeHeatFluxAnalyzer
}
\examples{
Run_LHFA(LakeName="Esthwaite",Folder="data",skipLoad=TRUE)
}
