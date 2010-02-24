\name{MethDataInput}
\alias{MethDataInput}

\title{Sequences match control}
\description{
Control function in order to check existence of sequences .faste files}
\usage{
MethDataInput(sFileName)
}

\arguments{
\item{sFileName}{String; path and name of the Tab delimited text file which includes the names and paths of sequences under study}
}
\details{
This procedure controls whether the *.fasta files written in a tabdelimited file are exact match es to the *.fasta in given paths. It is useful if the tab delimited text file was created by hand. This can prevent non existing clone sequence files stopping the workflow in the later analysis steps. The input is a character string which is the path and name of the tab delimited text file which includes the names and paths of sequences under study. It returns a data frame object with names and paths of existing clone sequences.
}

\value{Returns a data frame object with names and paths of existing analysed sequences}

\author{Arie Zackay <zackay@molgen.mpg.de>, Christine Steinhoff <steinhof@molgen.mpg.de>}


\examples{
## This pipeline will save the package example data  
## in tests directory under your R.home() directory 
## make sure that you have writing 
## permission under R.home() directory. If you do not have permission
## choose your own path.
 
dir.create(file.path(R.home(component="home"),"/BiqAnalyzer")) 
makeLocalExpDir(dataPath="/examples/BiqAnalyzer",
	localDir=file.path(R.home(component="home"),"/BiqAnalyzer"))
datameth <-MethDataInput(file.path(R.home(component="home"),
				"/BiqAnalyzer/PathFileTab.txt"))
}
\keyword{graphs}

