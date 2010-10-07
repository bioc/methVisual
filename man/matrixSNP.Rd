\name{matrixSNP}
\alias{matrixSNP}

\title{Correlation between methylation states}

\description{Correlation between methylation states for each CpG position and each sequence under study}

\usage{
matrixSNP(methData,correlation)
}
\arguments{
  \item{methData}{List; contains information on the pairwise alignments, and methylated CpG motifs.}
 \item{correlation}{calculation of corraltion on Matrix, if TRUE: calculation of spearman correlation matrix. if FALSE: computes a matrix with absolut numberes of methylation over all CpG position given methylation in a certain CpG sites.}
}

\details{
The function enable the user to explore cooccurrence between non neighbored CpG sites, which can be made by calculating all pairwise cooccurrences, due to correlation over all methylation CpG sites. Furthermore it is possible to display the computed distance cooccurrences in a plot containing all CpG sites and their correlation with other CpG sites.
The correlation values are color coded (gray levels) and the color coding bar is given beside the graph. The red numbers in the diagonal give the genomic position of each displayed CpG site.
}

\value{
Returns a correlation matrix displaying dependencies values for all CpG positions
}

\author{Arie Zackay <arie.zackay@mail.huji.ac.il>,  Christine Steinhoff <steinhof@molgen.mpg.de>}


\examples{
## using methData
data(methData)
matrixSNP(methData)
}
\keyword{graphs}

