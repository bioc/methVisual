\name{selectRefSeq}
\alias{selectRefSeq}

\title{Uploading genomic sequence}
\description{Uploading genomic sequence}
\usage{
selectRefSeq(sFileName)
}
\arguments{
  \item{sFileName}{path and name of genomic sequence to be uploaded}
}

\details{
Uploading genomic sequence
}

\value{
This function is used in order to read the reference sequence into the R environment. The reference sequences must be in fasta format. The user must give the path and name of the reference sequence in a character string.
}

\author{Arie Zackay <arie.zackay@mail.huji.ac.il>, Christine Steinhoff <steinhof@molgen.mpg.de>}


\examples{
## make sure that you have reading 
## permission under R.home() directory. If you do not have permission
## choose your own path.

dir.create(file.path(R.home(component="home"),"/BiqAnalyzer")) 
makeLocalExpDir(dataPath="/examples/BiqAnalyzer",
	localDir=file.path(R.home(component="home"),"/BiqAnalyzer"))
refseq <- selectRefSeq(file.path(R.home(component="home"),
			"/BiqAnalyzer/Master_Sequence.txt"))
}
\keyword{graphs}

