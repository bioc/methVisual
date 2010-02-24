\name{Cooccurrence}
\alias{Cooccurrence}

\title{Visualization of binary methylation data}
\description{Visualization of binary methylation data including neighboured cooccurrence}

\usage{
Cooccurrence(methData,file,real,lolli)
}
\arguments{
  \item{methData}{List; contains information on the pairwise alignments, and methylated CpG motifs.}
\item{file}{String; path and file name for saving the result. Default format is .pdf}
\item{real}{logical; real position (real=TRUE) or relative position (real=FLASE) according to the reference sequence}
\item{lolli}{Integer; size of lollipops}
}


\details{
Visualization of methylation states using lollipop graphs, percentage of methylation across experiments and value of neighboured cooccurrence due to calculation of spearman correlation.
Every single CpG site is marked with a circle with the following characteristic, a filled circle represents a methylated CpG and an empty circle a non methylated one.
The modification in methVisual lollipop plot is the calculation and visualization of dependencies (also known as cooccurrence) between neighbored methylated CpG site and non methylated CpG site. That means, given a set of bisulfite sequenced clones one would like to detect subgroups where specific CpG sites always occur coordinately either methylated or non methylated.
}

\value{
summary plot that will be saved by default as pdf file under the given path and name}

\author{Arie Zackay <zackay@molgen.mpg.de>, Christine Steinhoff <steinhof@molgen.mpg.de>}


\examples{
## using methData, file is the path to R home directory.
## In order to save Cooccurrence.pdf, make sure that you have writing 
## permission under R.home() directory. If you do not have permission
## choose your own path. 
#dir.create(file.path(R.home(component="home"),"/BiqAnalyzer")) 
data(methData)
Cooccurrence(methData,file=paste(tempdir(),"/Cooccurrence.pdf",sep="")) 
}
\keyword{graphs}


