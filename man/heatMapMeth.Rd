\name{heatMapMeth}
\alias{heatMapMeth}

\title{HeatMap diagram over methylation data}
\description{
HeatMap diagram of methylation classes over CpG position and sequence level 
}
\usage{
heatMapMeth(methData,file)
}
\arguments{
   \item{methData}{List; contains information on the pairwise alignments, and methylated CpG motifs}

\item{file}{optionally, quoted character string for specification of path and file  name for saving the result. By default, the result file's format is .pdf. If argument is omitted, screen output is provided only.}
}

\details{
Clustering is a prominent method for visualizing and studying groups of similar features, and is also widely used for the analysis of microarrays. In the case of analyzing DNA methylation datasets the matrix to be explored is a I*J binary matrix, were I are the clone sequences and J are the CpG positions. Every index in this matrix has the value 1 for methylated or 0 for non methylated CpG sites.
Providing a hierarchical clustering option based on the quality checked methylation data can be useful in finding clusters in the explored data in two dimensions, of the methylation state of CpG sites (J) and for distribution of methylation state over explored clone sequences (I). Very importantly, one has to keep in mind, that this method does not take
into account the genomic ordering of CpG sites.
Clustering methods can be also useful in quality control check. Observing clones from the same PCR product in different clusters can point to bad quality clone sequences. The clustering method which was used in methVisual R package is the heatmap() function from stat package which is a color image with two dendrograms added to the sides of the columns and rows. Since the data is binary, the distances that are calculated among the columns and among the rows are computed with a binary distance function which is the asymmetric binary function. This method assumes that non zero elements are ON and
zero elements are OFF. The distance is the proportion of bits in which only one is ON amongst those in which at least one is ON.
Based on aligned sequences under study a heatmap is created that displays two way clustering of methylation status of all sequences and all aligned CpG positions.
}

\value{
 Heat-Map image is displayed and optionally saved as postscript in given path and name. 
}

\author{Arie Zackay <arie.zackay@mail.huji.ac.il>,Christine Steinhoff <steinhof@molgen.mpg.de>}


\examples{
## using methData, file is the path to R home directory.
## In order to save heatMapMeth.pdf, make sure that you have writing 
## permission under R.home() directory. If you do not have permission
## choose your own path.  
dir.create(file.path(R.home(component="home"),"/BiqAnalyzer")) 
data(methData)
heatMapMeth(methData,file=file.path(R.home(component="home"),
			"/BiqAnalyzer/heatMapMeth.pdf"))
}
\keyword{graphs}

