\name{MethylQC}
\alias{MethylQC}

\title{Quality controle (QC) on Methylation Data}
\description{
Processing a quality control (QC) procedure on bisulafite sequences
}

\usage{
MethylQC(refSeq, methFileDataFrame, makeChange, identity,conversion)
}
\arguments{
  \item{refSeq}{String; genomic sequence, see selectRefSeq()}
  \item{methFileDataFrame}{Data frame; sequences names and their paths, see MethDataInput()}
  \item{makeChange}{Logical; if TRUE changes take place automatically, by default TRUE}
  \item{identity}{min. identity value , by default 80 percent}
  \item{conversion}{min conversion rate, by default 90 percent}
}

\details{
In order to avoid bad qualitative data entering this methylation analysis three measurements are made: 1) alignment check: if reverse, complement or reverse-complement 2) sequence identity between genomic sequence and every examined sequence 3) bisulfite conversion
}

\value{
Returns data frame of sequences names after QC and their paths saves QCINFO.Rdata under sequence files directory
}

\author{Arie Zackay <zackay@molgen.mpg.de>, Christine Steinhoff <steinhof@molgen.mpg.de>}


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
}
\keyword{graphs}

