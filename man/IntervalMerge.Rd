% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helper_functions.R
\name{IntervalMerge}
\alias{IntervalMerge}
\title{map MSMS data onto camera object}
\usage{
IntervalMerge(
  cameraObject,
  MSMSdata,
  PlusTime,
  MinusTime,
  ppm,
  listOfMS2Mapped = list(),
  listOfUnMapped = list()
)
}
\arguments{
\item{cameraObject}{A CAMERA object including addcuts and isotopes}

\item{MSMSdata}{a list of MSnbase objects}

\item{PlusTime}{The higher boundary of the retention time window}

\item{MinusTime}{The lower boundary of the retention time window}

\item{ppm}{ppm mass devation}

\item{listOfMS2Mapped}{An empty or full list of previously unmapped MS2s (the new mapped data will be added to this)}
}
\value{
A list of mapped and/or unmapped MS2s
}
\description{
This function accepts a CAMERA object, MSMS data and map MSMS data onto the features
}
\keyword{internal}
