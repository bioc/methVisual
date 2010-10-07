\name{plotMatrixSNP}
\alias{plotMatrixSNP}

\title{Plot of methylation states dependencies}
\description{
Visulisation of dependencies between methylation states over explored bisulafite sequences
}
\usage{
plotMatrixSNP(summeryMatrix,methData,file)
}

\arguments{
  \item{summeryMatrix}{see matrixSNP()}
   \item{methData}{List; contains information on the pairwise alignments, and methylated CpG motifs.}
 \item{file}{String; quoted character string for specification of path and file name for saving the result. The result file is in .pdf format}
}

\details{
The SNP Plot produce a visualisation of the cooccurrence between methylation states on CpGs over all explored sequences. The index of the CpGs is based on there position on genomic sequence.
}

\value{
SNP Plot as pdf file saved in given path and name. 
}

\author{Arie Zackay <arie.zackay@mail.huji.ac.il>, Christine Steinhoff <steinhof@molgen.mpg.de>}


\examples{
## using methData

data(methData)
summeryMatrix <- matrixSNP(methData)

## using methData, file is the path to R home directory.
## In order to save plotMatrixSNP.pdf, make sure that you have writing 
## permission under R.home() directory. If you do not have permission
## choose your own path. 
dir.create(file.path(R.home(component="home"),"/BiqAnalyzer")) 
data(methData)
summary <- matrixSNP(methData)
plotMatrixSNP(summary,methData,file=file.path(R.home(component="home"),
			"/BiqAnalyzer/","plotMatrixSNP.pdf")) 
}
\keyword{graphs}

