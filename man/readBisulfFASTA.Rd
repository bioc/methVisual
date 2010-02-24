\name{readBisulfFASTA}
\alias{readBisulfFASTA}

\title{Read multiple FASTA file}
\description{
Read multiple FASTA file and write bisulfite sequences in separate FASTA files for every bisulfite sequence}
\usage{
readBisulfFASTA(sFileName,sDirName)
}
\arguments{
  \item{sFileName}{String; path to multiple FASTA file}
  \item{sDirName}{String;  path to directory of created separated FASTA files}
}


\author{Arie Zackay <zackay@molgen.mpg.de>, Christine Steinhoff <steinhof@molgen.mpg.de>}


\examples{
## In order to use the following example 
## make sure that you have writing permission under R.home()
## directory. If you do not have permission choose your own path. 
dir.create(file.path(R.home(component="home"),"/BiqAnalyzer")) 
makeLocalExpDir(dataPath="/examples/BiqAnalyzer",
	localDir=file.path(R.home(component="home"),"/BiqAnalyzer"))

	readBisulfFASTA(sFileName=file.path(R.home(component="home"),"/BiqAnalyzer/multiFASTA.fasta"), sDirName= file.path(R.home(component="home"),"/BiqAnalyzer/"))
}

\keyword{graphs}




