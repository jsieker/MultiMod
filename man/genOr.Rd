
\name{genOr}
\alias{genOr}
\title{Gene-oriented multiple module membership clustering}
\usage{
genOr(threshold, kME, cutVal)
}
\arguments{
\item{threshold}{ A number between zero and one that serves as a cutoff of how well a gene must correlate with a certain module before it is considered part of that module in the analysis. A typical value is around 0.75, with more stringent analyses occuring around 0.90.}

\item{kME}{The module eigenvalues outputted by WGCNA's 'blockwiseModules' or 'moduleEigengenes' functions.}

\item{cutVal}{ Allows the reduction of meaningless parts of the output. Defaults to true but can be turned to false for the sake of standardizing the output dataframe dimensions no matter the number of modules in the sample.}
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
Allows for a gene-oriented, as opposed to a module-oriented, output format of multiple module clustering. The analysis of these two methods are extremely similar, but the outputs are useful for very different types of analyses.
}
