\name{makeTabFilePath}
\alias{makeTabFilePath}

\title{Tab delimited text file}

\description{create tab delimited text file with PATH and FILE columns}

\usage{
makeTabFilePath(localDir)
}
\arguments{
  \item{localDir}{String; path to local directory where sequence files and tab delimited text file are saved}
 
}

\details{
The function creates a tab delimited text file with information on path and file names for files which are stored in a given directory. This tab-delimited text file is the input file that contains the clone sequences. The makeTabFilePath() function was implemented in order to create a control step before starting the analysis. The best way of using the function is by collecting all clone sequence files which need to be analyzed into one empty directory and apply this function to the directory.
}

\author{Arie Zackay <arie.zackay@mail.huji.ac.il>, Christine Steinhoff <steinhof@molgen.mpg.de>}


\examples{
## on tests directory under R home directory
## make sure that you have writing permission under 
## R.home() directory. If you do not have permission
## choose your own path (localDir=YOUR_OWN_PATH/).  

makeTabFilePath(file.path(R.home(component="home"),"/BiqAnalyzer"))

}
\keyword{graphs}
