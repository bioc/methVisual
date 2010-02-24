\name{cgMethFinder}
\alias{cgMethFinder}

\title{Methylation status}
\description{
CpGs methylation status on clone sequence}
\usage{
cgMethFinder(ref,str)
}
\arguments{
  \item{ref}{String, genomic sequence, see selectRefSeq()}
  \item{str}{String, Single sequence under study after alignment to ref}
}

\details{
The function determines the methylation status of each CpG site by comparing TpG and CpG sites within the clone sequence to corresponding CpG sites in the reference sequence. The input values are the reference sequence and one of the clone sequences which is explored. It returns a (0,1) vector. 1 stands for methylated and 0 for non methylated state.
This function is used in the methVisual package as internal function for the calculation of the methylation profiles.
}

\value{Returns a (0,1) vector. 1 stands for methylation and 0 for non methylation status.}

\author{Arie Zackay <zackay@molgen.mpg.de>, ,Christine Steinhoff <steinhof@molgen.mpg.de>}

\examples{
ref <- "TTCGGGATCGTTTTTTTAGTAGGTCGGAAGTTTCGTTATGGATTCGTTTTTC"
str <- "TTCGGGATCGTTTTTTTAGTAGGTTGGAAGTTTTGTTATGGATTCGTTTTTC"
cgMethFinder(ref,str)
}
\keyword{graphs}

