% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilis.R
\name{as_data_frame_rowname}
\alias{as_data_frame_rowname}
\title{modified version of as.data.frame for rownames}
\usage{
as_data_frame_rowname(list)
}
\arguments{
\item{list}{list. list to be converted}
}
\value{
dataframe. will return t(as.data.frame(list))
}
\description{
modified version of as.data.frame to remove the changed caused by make.name
}