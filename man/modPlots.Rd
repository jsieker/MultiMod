
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

\item{ord}{This determines the order of the samples in the plots. If no argument is supplied or if ord=="intensity" the order will be based on module eigengene intensity. If ord=="alpha", the order will be alphabetical by sample name. For a custom order, supply a concatenated vector of names that match your samples. All samples must be given a location and all spelling and spacing must be copied perfectly or the function is likely to fail. If providing custom order, write it very carefully. }

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
