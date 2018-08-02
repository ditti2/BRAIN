\name{useBRAIN2}
\alias{useBRAIN2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function computing probabilities of aggregated isotopic variants using heuristics.
}
\description{
Function computing probabilities of isotopic variants  using heuristics,  for chemical components composed of carbon, hydrogen, oxygen, nitrogen 
and sulfur (e.g. peptides). Additionally the function returns also the monoisotopic mass and the average mass of given chemical component.
}
\usage{
useBRAIN2(aC, stopOption = "nrPeaks", nrPeaks, approxStart = 1, approxParam = NULL))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{aC}{
List with fields C, H, N, O, S of integer non-negative values (if any field is ommited, then its value is set to 0).
} 
\item{stopOption}{
only option "nrPeaks" allowed 
}
    \item{nrPeaks}{
Integer indicating the number of consecutive isotopic variants to be calculated,
starting from the monoisotopic one. This value can always be provided, even if <stop.option>
is not a default setting. In the latter case it is a hard stopping criterion.
}
    \item{approxStart}{
Integer indicating the number of first isotopic peak to be calculated
}
    \item{approxParam}{
Integer indicating the length of recurrence (see RCL in [Ditt2])
}
}
\details{
Function uses RCL and LSP heuristics from [Ditt2].
}
\value{
Named list with the following fields:
\itemize{
\item{isoDistr}{Probabilities of aggregated isotopic variants (numeric vector)}
}
}
\references{
[Ditt2] Dittwald P., Valkenborg D. BRAIN 2.0: time and memory complexity improvements in the algorithm for calculating the isotope distribution.
JASMS, 2014, doi: 10.1007/s13361-013-0796-5.
}
\author{
Piotr Dittwald <piotr.dittwald@mimuw.edu.pl>
}
\note{
Remember that the isotopic variants starts from monoisotopic one. For large chemical molecules, first masses may have very low abundances. So sufficient number of peaks should be calculated to reach most abundant isotopic variant.

If only isotopic probabilities are needed, then the function calculateIsotopicProbabilities should be used.
}

\seealso{
 \code{\link{calculateIsotopicProbabilities}}
}
\examples{
  nrPeaks = 1000  
  aC <-  list(C=23832, H=37816, N=6528, O=7031, S=170)  # Human dynein heavy chain  
  res <-  useBRAIN(aC = aC, stopOption="nrPeaks", nrPeaks = nrPeaks)
  res2 <-  useBRAIN2(aC = aC, stopOption="nrPeaks", nrPeaks = nrPeaks, approxStart = 10)
  old = res$iso[10:109]/res$iso[11:110]
  new  = res2$iso[1:100]/res2$iso[2:101]
  max(old - new)
  max((old - new)/old)
  res3 <-  useBRAIN2(aC = aC, stopOption="nrPeaks", nrPeaks = nrPeaks, approx=TRUE, approxParam = 10)
  max((res3$iso - res$iso)/res$iso)
}

