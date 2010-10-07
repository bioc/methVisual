\name{cgInAlign}
\alias{cgInAlign}

\title{Amount of CpGs}
\description{Calculating amount of CpGs between alignments border}
\usage{
cgInAlign(methData)
}
\arguments{
\item{methData}{An object of type list that contains information on the pairwise alignments and methylation status of all CpG motifs under study. Created by applying MethAlignNW()}
}

\details{
This function computes the amount of CpGs on the positions corresponds to genomic sequence over all sequences under study between the alignment borders.
}

\value{
Integer vector of CpG amount
}

\author{Arie Zackay <arie.zackay@mail.huji.ac.il>, Christine Steinhoff <steinhof@molgen.mpg.de>}

\examples{
data(methData)
cgInAlign(methData)
}
\keyword{graphs}
