\name{makeLocalExpDir}
\alias{makeLocalExpDir}

\title{Saving example data}

\description{Saving example data as provided by the package}

\usage{
makeLocalExpDir(dataPath,localDir)
}

\arguments{
  \item{dataPath}{String; path to the location of sequences under study and genomic reference sequence}

   \item{localDir}{string path to local directory for transferring the sequence files and Tab delimited text file}
}	

\details{
Help function for saving data provided along with the package followed by creation of a tab delimited text file with information on PATHs and FILEs
}

\author{Arie Zackay <zackay@molgen.mpg.de>, Christine Steinhoff <steinhof@molgen.mpg.de>}


\examples{
## saving example data under R.home() directory. 
## make sure that you have writing permission under 
## R.home() directory. If you do not have permission
## choose your own path (localDir=YOUR_OWN_PATH/).
dir.create(file.path(R.home(component="home"),"/BiqAnalyzer")) 
makeLocalExpDir(dataPath="/examples/BiqAnalyzer",
	localDir=file.path(R.home(component="home"),"/BiqAnalyzer"))
}
\keyword{graphs}

