\name{methFisherTest}
\alias{methFisherTest}

\title{Fisher exact Test on methylation Data}
\description{
Fisher exact Test on two subsets of experiments over matched CpG sites
}
\usage{
methFisherTest(methData,set1,set2)
}
\arguments{
  \item{methData}{List; contains information on the pairwise alignments, and methylated CpG motifs.}
 \item{set1}{First subset - Integer vector of experiments due to there order at methData}
 \item{set2}{Second subset - Integer vector of indexes of experiments due to there order at methData}
}

\details{
Given two clone sequences groups A and B for each CpG site the user can investigate whether there is a difference of methylation status between the two groups at each of the CpG sites. In order to calculate this difference at each CpG site, the two-tailed p-value of Fisher's exact test is calculated from the 2*2 tables at each CpG site. This p-value indicates the level of difference at every single CpG in those two groups of clone sequences. A p-value smaller than 0.05 is an indication for the independence of the methylation state in a certain CpG site when comparing two groups of clone sequences.
}

\value{
P- Values vector and a Plot
}

\author{Arie Zackay <arie.zackay@mail.huji.ac.il>, Christine Steinhoff <steinhof@molgen.mpg.de>}


\examples{
data(methData)
methFisherTest(methData,c(1,2,3),c(4,5,6))

}
\keyword{graphs}

