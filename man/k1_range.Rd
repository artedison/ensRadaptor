% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/k1_range.R
\name{k1_range}
\alias{k1_range}
\title{k1 range production}
\usage{
k1_range(list.reac.addon = NULL, list.kine = NULL, varn = 1)
}
\arguments{
\item{list.reac.addon}{list. the list of reactions. must be provided}

\item{list.kine}{kinetic list. parameter list. must be provided}

\item{varn}{numeric. how variance the range can be, the range is defined as mean +- varn*sigma
default 1}
}
\value{
list. formualted k1 range
}
\description{
this is for producing a constraint enough k1 range to start the training
based on the assumption that each enzyme will have similar region size for kinetic parameter
depricated soon
}
