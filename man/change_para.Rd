% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/change.R
\name{change_para}
\alias{change_para}
\title{Editing or showing ensemble input files (ens.i01) for single term}
\usage{
change_para(term = NULL, value = NULL, infile = NULL,
  outfile = NULL, type = "show", indrep = 1)
}
\arguments{
\item{term}{string. the term that need to be changed. must be provided}

\item{value}{array, dim(array)=1. the set value. default NULL}

\item{infile}{string. input file path. must be provided}

\item{outfile}{string. output file path. default NULL}

\item{type}{string. "show"(show the parameter, default), "edit"(edit the value)
if there is mulitple locations for the corresponding term, only show value is possible
and both location/index and value are returned.}

\item{indrep}{array. the index of the one that need to be changed if there are multiple (edit). default 1}
}
\value{
If type="edit", no value is returned and just change the file correspondingly
        If type="show" OR length(uncomm)!=1&&indrep==1, a list of value and line number will be returned
}
\description{
The function will edit or show the value of specific parameters
}
\seealso{
[change_block()] for modifying one block
}
