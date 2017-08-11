
\name{modOr}
\alias{modOr}
\title{Module-oriented multiple module membership clustering}
\usage{
modOr(threshold, kME, cap)
}
\arguments{
\item{threshold}{ A number between zero and one that serves as a cutoff of how well a Module must correlate with a certain module before it is considered part of that module in the analysis. A typical value is around 0.75 (the default), with more stringent analyses occuring around 0.90.}

\item{kME}{The module eigenvalues outputted by WGCNA's 'blockwiseModules' or 'moduleEigenModules' functions. This is a required input.}

\item{cap}{ Sets a maximum number of modules to which a gene may be assigned. The default is five. }
}
\references{
  Langfelder P and Horvath S, WGCNA: an R package for weighted
  correlation network analysis. BMC Bioinformatics 2008, 9:559
  doi:10.1186/1471-2105-9-559

  Peter Langfelder, Steve Horvath (2012). Fast R Functions for Robust
  Correlations and Hierarchical Clustering. Journal of Statistical
  Software, 46(11), 1-17. URL http://www.jstatsoft.org/v46/i11/.
}
\description{
Allows for a module-oriented, as opposed to a gene-oriented, output format of multiple module clustering. The analysis of these two methods are extremely similar, but the outputs are useful for very different types of analyses.
}
