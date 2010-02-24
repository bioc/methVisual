\name{plotAbsMethyl}
\alias{plotAbsMethyl}

\title{Plot number mathylation on CpGs positions}
\description{
Plot of absolute/relative number of mathylation on aligned CpGs positions   
}
\usage{
plotAbsMethyl(methData,real)
}
\arguments{
  \item{methData}{List; contains information on the pairwise alignments, and methylated CpG motifs.}
  \item{real}{real position (real=TRUE) or relative position (real=FLASE) on reference sequence}
}

\details{
This function generates a plot of the absolute number of methylation in all CpG sites over all examined sequences. It returns a vector with the absolute number of methylation sites of all examined clone sequences. The user supplies a list
object containing information about the pairwise alignments, and methylated CpG sites as calculated by MethAlignNW() or makeDataMethGFF(). 
The user can display the absolute number of methylation over the CpG sites according
to their relative position on the reference sequence.
}

\value{
A vector with the absolute/relative number of methylation on genomic sequence CpG positions over all examined sequences

}

\author{Arie Zackay <zackay@molgen.mpg.de>, Christine Steinhoff <steinhof@molgen.mpg.de>}


\examples{
## using methData
data(methData)
plotAbsMethyl(methData,real=FALSE)
}
\keyword{graphs}

