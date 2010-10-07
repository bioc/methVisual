\name{MethAlignNW}
\alias{MethAlignNW}

\title{Summary of methylation states}
\description{
Summarize the methylation states after calculating pairwise alignments of each examined sequences and the genomic sequence}
\usage{
MethAlignNW(refSeq, QCdata, alignment)
}
\arguments{
  \item{refSeq}{String; Genomic sequence as String format}
  \item{QCdata}{Data frame;  Names and paths of analysed sequences after quality control}
  \item{alignment}{ If TRUE, alignments are included in summery, else not included}
}

\details{
Given aligned sequences after quality control, the function returns a list object with the following data: sequences name, methylation state on CpG position, start and end position of alignments and length of genomic sequence. The data includes the core information for the exploratory analysis and visualizations.
}

\value{
Returns a summery on sequence alignments and methylation states
}

\author{Arie Zackay <arie.zackay@mail.huji.ac.il>, Christine Steinhoff <steinhof@molgen.mpg.de>}


\examples{
## In order to use the following example 
## make sure that you have writing permission under R.home()
## directory. If you do not have permission choose your own path. 
dir.create(file.path(R.home(component="home"),"/BiqAnalyzer")) 
makeLocalExpDir(dataPath="/examples/BiqAnalyzer",
	localDir=file.path(R.home(component="home"),"/BiqAnalyzer"))
datameth <-MethDataInput(file.path(R.home(component="home"),
			"/BiqAnalyzer/PathFileTab.txt"))
refseq <- selectRefSeq(file.path(R.home(component="home"),
			"/BiqAnalyzer/Master_Sequence.txt"))
QCdata <- MethylQC(refseq, datameth)
methData <- MethAlignNW( refseq , QCdata)
}
\keyword{graphs}




