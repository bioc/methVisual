\name{methCA}
\alias{methCA}

\title{CA on methylation}
\description{
Correspondence Analysis over methylation Data
}
\usage{
methCA(methData,file)
}
\arguments{
    \item{methData}{List; contains information on the pairwise alignments, and methylated CpG motifs.}
\item{file}{String; optionally, quoted character string for specification of path and file  name for saving the result. By default, the result file's format is .pdf. If argument is omitted, screen output is provided only.}
}

\details{
Correspondence analysis (CA) is a multivariate statistical technique which is applicable to tables of categorical data. CA can be useful in understanding the data in a inter relationship manner, meaning the dependencies between categories. Unlike conventional statistical methods that try to prove a hypothesis, the CA is an exploratory technique that can reveal data content. With the help of a graphical application, which displays each category as a point in two dimensional plot, it is easier to locate special characteristic in the numerical data. It used in order to study the association between CpG methylation states and clone sequences based on the aligned clone sequences under study.
}

\value{
CA diagram as postscript file saved in given path and name. 
}

\author{Arie Zackay <arie.zackay@mail.huji.ac.il>, Christine Steinhoff <steinhof@molgen.mpg.de>}


\examples{
## using methData, file is the path to R home directory.
## In order to save methCA.pdf, make sure that you have writing 
## permission under R.home() directory. If you do not have permission
## choose your own path.
dir.create(file.path(R.home(component="home"),"/BiqAnalyzer"))  
data(methData)
methCA(methData,file=file.path(R.home(component="home"),
			"/BiqAnalyzer/methCA.pdf")) 
}
\keyword{graphs}

