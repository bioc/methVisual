\name{coversionGenom}
\alias{coversionGenom}

\title{Sequence conversion}
\description{Bisulfite conversion of genomic Sequence}

\usage{
coversionGenom(genomicSeq)
}
\arguments{
  \item{genomicSeq}{String; genomic sequence}
}

\details{
This method simulates the bisulfite reaction by converting all Cs outside from CpG sites into Ts. In doing so, the percent of identity between genomic sequence and sequence under study can be determined. The input is the reference sequence that was selected by the function selectRefSeq().
The method returns the reference sequence as a string object with Cs converted into Ts.
}

\value{
Returns String with C converted into T
}

\author{Arie Zackay <arie.zackay@mail.huji.ac.il>, Christine Steinhoff <steinhof@molgen.mpg.de>}


\examples{
genomicSeq <- "ACCGTTTGGCC"
coversionGenom(genomicSeq)
}	
\keyword{graphs}
