\name{findNonAligned}
\alias{findNonAligned}

\title{Aligned CpG positions}

\description{Determination of non-aligned CpG positions}

\usage{
findNonAligned(methData)
}
\arguments{
 \item{methData}{List; contains information on the pairwise alignments, and methylated CpG motifs}
}

\details{
Determination of aligned and not-aligned positions of CpGs of examined sequences in relation to genomic sequence.
}

\value{
Integer matrix, where columns=CpG positions, row=clone sequences. 0 = not methylated, 1 = methylated, 2 = not aligned
}

\author{Arie Zackay <arie.zackay@mail.huji.ac.il>, Christine Steinhoff <steinhof@molgen.mpg.de>}

\examples{
data(methData)
findNonAligned(methData)
}
\keyword{graphs}

