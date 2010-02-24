#######################################################################
#######################################################################
###
### This Programms checks if a tab delimited of *.fasta files are exact match to the
### *.fasta in given Directories										
###																 
### Author     : Arie Zackay											
###																	
### Last change: 27. Oktober. 2008										
###																	
### Change log														
###																	
#######################################################################
#######################################################################

MethDataInput <- function(sFileName) {

	if (sFileName != "") {

	dPathAndFileNames <- read.delim(paste(sFileName), header = TRUE, sep = "\t", quote="\"")
   		
		if (is.data.frame(dPathAndFileNames)==TRUE) {

			iPathPos <- grep ("PATH",ignore.case = TRUE, names(dPathAndFileNames))
			iFilePos <- grep ("FILE", ignore.case = TRUE ,names(dPathAndFileNames))

			if (length(iPathPos) > 0 && length(iFilePos) > 0) { 

				vbFileExists <- rep(FALSE, dim(dPathAndFileNames)[1])

				for (i in 1 : length(vbFileExists)) {

					sFileName <- paste( paste(dPathAndFileNames [i, iPathPos]), '/', paste(dPathAndFileNames [i, iFilePos]), sep="")
					vbFileExists [i] <- file.exists(sFileName)
				} 

        
		vFileChk<-c()      
		vPathChk<-c()

			if(sum(vbFileExists)>0) {
                            
				for(i in 1:sum(vbFileExists)) {

					vFileChk[i]<-paste(dPathAndFileNames[grep(TRUE,vbFileExists)[i], iFilePos])
					vPathChk[i]<-paste(dPathAndFileNames[grep(TRUE,vbFileExists)[i], iPathPos])
				}

			dPathAndFileChk<-data.frame(FILE=vFileChk, PATH=vPathChk)
         		}         
         
			if (sum(vbFileExists) != length(vbFileExists)) {
	
				cat("The following files do not exist.",'\n')

				for (i in 1 : length(vbFileExists)) {

					if (vbFileExists[i] == FALSE) { 
						cat("Row:", i + 1, " Path: \"", paste(dPathAndFileNames[i,iPathPos]), "\" Filename: \"", paste(dPathAndFileNames[i,iFilePos],"\"",sep=""),".",'\n',sep="")
					}
				}
			}   
	}   


	else {

		if (length(iPathPos) == 0) {

		cat ("Did not find column \"PATH\" in header of file \"", sFileName, "\"") ## was stop
		}
		
		if (length(iFilePos) == 0) {

		cat ("Did not find column \"FILE\" in header of file \"", sFileName, "\"")
		}

		} 
  	 }
	}

	else {
		cat("No file selected. Program stopped.")  ## was stop
	}

   	
	return(dPathAndFileChk)

gc(FALSE)

}


########################################################################
########################################################################
## The function reads reference sequence for the Methylation analysis
########################################################################
########################################################################

selectRefSeq <- function(sFileName) {


if (sFileName != "") {

	refSeq	<- toupper(scan(sFileName,what='character',skip=1))
}

else {
	cat("No file selected. Program stopped.")
}

return(refSeq)

}


##############################################################################
##############################################################################
## creating local data directory from examples data which saved under Package
##############################################################################
##############################################################################

makeLocalExpDir <- function(dataPath,localDir) {

ref <- readLines(paste(system.file(paste(dataPath, sep = ""), package = "methVisual"),"/Master_Sequence.txt", sep = ""), warn = FALSE)

writeLines(ref,paste(localDir,"/Master_Sequence.txt",sep=""))


fasta <- readLines(paste(system.file(paste(dataPath, sep = ""), package = "methVisual"),"/multiFASTA.fasta", sep = ""), warn = FALSE)

writeLines(ref,paste(localDir,"/multiFASTA.fasta",sep=""))


seq_dirpackage <- list.files(system.file(paste(dataPath,"/Demo_data",sep=""), package = "methVisual"))

	for (i in 1:length(seq_dirpackage)) {

		biseq <- readLines(paste(system.file(paste(dataPath,"/Demo_data",sep=""),package="methVisual"),"/",seq_dirpackage[i],sep=""),warn = FALSE)
		print(seq_dirpackage[i])
		writeLines(biseq,paste(localDir,"/",seq_dirpackage[i],sep=""))
	}

FILE <- seq_dirpackage
PATH <- rep(localDir,length(seq_dirpackage))
filepath <- data.frame(FILE=FILE, PATH=PATH)

write.table(filepath,paste(localDir,"/","PathFileTab.txt",sep=""),row.names = FALSE,quote=FALSE,sep='\t')
 
return(filepath)

}


####################################################################################
####################################################################################
## creating tab delimeted text file with FILE and PATH columns from a given directory  
####################################################################################
####################################################################################

makeTabFilePath <- function(localDir) {

seq_dir <- list.files(localDir)
FILE <- seq_dir
PATH <- rep(localDir,length(seq_dir))
filepath <- data.frame(FILE=FILE, PATH=PATH)

write.table(filepath,paste(localDir,"/","PathFileTab.txt",sep=""),row.names = FALSE,quote=FALSE,sep='\t')
 
return(filepath)

}


#################################################
#################################################
##
## Read clone sequences from multiple FASTA file
##
#################################################
#################################################


readBisulfFASTA <- function (sFileName, sDirName){

seqDNA <- c() 
seqName <- c()

if (sFileName != "") {

  bis <- readFASTA(sFileName, strip.descs=TRUE)
  seqDNA <- sapply(bis, function(x) x$seq)
  seqName <-  sapply(bis, function(x) x$desc)

for (i in 1:length(seqDNA)) {

write.table(paste(">",seqName[i],'\n', seqDNA[i],sep=""), file=paste(sDirName,'/',seqName[i],".fasta",sep=""), col.names=FALSE, row.names=FALSE,quote = FALSE)

}

}

else {
	cat("No file selected. Program stopped.")
}

cat(paste("bisulafite sequences in FASTA file format in",sDirName),'\n')

}



