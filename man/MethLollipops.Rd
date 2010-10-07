\name{MethLollipops}
\alias{MethLollipops}

\title{"Lollipops" Methylation plot}
\description{Visualization of methylation patterns in terms of "Lollipops" plot}
\usage{
MethLollipops(methData)
}
\arguments{
  \item{methData}{List; contains information on the pairwise alignments, and methylated CpG motifs.}
}


\details{
The lollipops plot allows the user to study the genomic localization and states of CpG sites. Each circle marks a CpG site under study. Full circles display methylated CpGs and the non filled ones stand for non methylated CpG states. The examined sequences are aligned with respect to the genomic sequence in order to allow for an intuitive visualization of methylation states according to their genomic position.
}

\value{Returns a "lollipop" plot}

\author{Arie Zackay <arie.zackay@mail.huji.ac.il>, Christine Steinhoff <steinhof@molgen.mpg.de>}

\examples{
data("methData")
MethLollipops(methData)
}
\keyword{graphs}

