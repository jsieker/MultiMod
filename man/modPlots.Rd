
\name{modPlots}
\alias{modPlots}
\title{Plot Module Eigengenes of each module for each sample}
\usage{
modPlots(expr, modOrOut)
}
\arguments{
\item{expr}{Sample expression data in a dataframe with rows as
samples and columns as genes. }

\item{modOrOut}{This is the output from the modOr (Module-Oriented) function in this package.}
}
\value{
A list the same length as the number of modules. If a module has no assigned genes (and thus no expression data), it will have a zero instead of a plot in the output list.
}
\description{
Plot Module Eigengenes across samples
}
\references{
  Langfelder P and Horvath S, WGCNA: an R package for weighted
  correlation network analysis. BMC Bioinformatics 2008, 9:559
  doi:10.1186/1471-2105-9-559

  Peter Langfelder, Steve Horvath (2012). Fast R Functions for Robust
  Correlations and Hierarchical Clustering. Journal of Statistical
  Software, 46(11), 1-17. URL http://www.jstatsoft.org/v46/i11/.
}
