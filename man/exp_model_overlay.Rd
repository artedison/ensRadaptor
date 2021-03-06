% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/exp_model_overlay.R
\name{exp_model_overlay}
\alias{exp_model_overlay}
\title{overlay experiment measurement with model predicted result in o01}
\usage{
exp_model_overlay(o01 = NULL, species = NULL, addon = "",
  exp = NULL, dir.res = NULL, logtrans = FALSE)
}
\arguments{
\item{o01}{string. file path. must be provided}

\item{species}{string. the species that need overlay. must be provided}

\item{addon}{string. addon on output file name. default ""}

\item{exp}{numeric. the experiment number that need to be plot. must be provided}

\item{dir.res}{string. the output folder. must be provided}

\item{logtrans}{bool. whether log transformaed is performed on dat. default FALSE}
}
\value{
figure. return the figure itself
}
\description{
measurement is actually "scaled" measurement
the range [low high] will also be plot
}
