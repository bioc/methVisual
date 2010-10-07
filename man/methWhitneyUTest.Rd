\name{methWhitneyUTest}
\alias{methWhitneyUTest}

\title{Mann Whitney U-Test on methylation data}
\description{
Mann Whitney U-Test on entire sets of CpG sites
}
\usage{
methWhitneyUTest(methData,set1,set2)
}
\arguments{
  \item{methData}{List; contains information on the pairwise alignments, and methylated CpG motifs.}
 \item{set1}{First subset - Integer vector of experiments due to there order at methData}
 \item{set2}{Second subset - Integer vector of indexes of experiments due to there order at methData}
}

\details{
Mann Whitney U-Test (known also as wilcoxon-rank-sum test), is a non parametric test for assessing whether two independent samples of observations come from the same distribution. The null hypothesis in the Mann-Whitney U test is that the two samples are from a single population, which means that their probability distributions are equal. In the methVisual package the Mann-Whitney U test is applied in order to test if the distribution of methylated and non methylated sites in the profile under study between the given experimental sub groups of clone sequences is different. In order to calculate it, the two-tailed p-value of the Mann-Whitney U test is computed by ranking the ratios of methylated CpG sites to all CpG sites in a given clone sequence. The p-value in this case indicates if the distribution of ratio over all methylation states differ between the two groups.
}

\value{
Returns a p-value
}

\author{Arie Zackay <arie.zackay@mail.huji.ac.il>, Christine Steinhoff <steinhof@molgen.mpg.de>}


\examples{
## using methData
data(methData)
methWhitneyUTest(methData,c(1,2,3),c(4,5,6))
}
\keyword{graphs}

