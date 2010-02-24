#############################################################
#############################################################
### Alignment of bisulfite sequences to a referance sequence
#############################################################
#############################################################

MethAlignNW <- function (refSeq,  QCdata, alignment ) {

seqName <- c()
alignment <- c()
files <- paste( QCdata$FILE) 
paths <- QCdata$PATH 
mat <- matrix(0,length(DNA_ALPHABET),length(DNA_ALPHABET))
mat[1:4,1:4] <- c(1,0,0,0,0,1,0,1,0,0,1,0,0,0,0,1)
rownames(mat) <- colnames(mat) <- DNA_ALPHABET[1:length(DNA_ALPHABET)]
lengthRef <-nchar(refSeq)
cg_ref <- gregexpr("CG",refSeq)
positionCGIRef <- sort(cg_ref[[1]] [1: length(cg_ref[[1]])] )
methPos <- matrix(nrow=length( QCdata$PATH),ncol=length(cg_ref[[1]]))
startEnd <- matrix(nrow=length( QCdata$PATH),ncol=2)

for (i in 1: length (paths)) {

	seqFileTemp <- paste( paths[i],'/',files[i], sep="")
	seqTemp	<- toupper(scan(seqFileTemp,what='character', sep="",quiet=TRUE))
	align <- pairwiseAlignment(refSeq, seqTemp,substitutionMatrix=mat,type="local-global")
	startRef <- start(pattern(align))
	endRef <- end(pattern(align))
#	lengthAlignRef <- nchar(pattern(align))
	startEnd_temp <- c(startRef,endRef)
	startEnd[i,] <- startEnd_temp
	newRef <- paste(substr(refSeq,1,startRef-1),pattern(align), substr(refSeq,endRef+1,lengthRef),sep="")
	newSeq <- paste( paste(rep('-',startRef-1),collapse=""), subject(align), paste(rep('-',lengthRef-endRef)),sep="")
	alignment[i] <- paste(subject(align))

	methPos_Temp <- cgMethFinder (newRef,newSeq)
	methPos[i,] <- methPos_Temp
	seqName [i] <- files [i]
	cat('Alignment with ',seqName[i],' done','\n')  
}


if (missing(alignment) || alignment==FALSE) {
mathylAlign <- list(seqName=seqName, methPos=methPos, positionCGIRef=positionCGIRef, startEnd=startEnd, lengthRef=nchar(refSeq))

}

else {
mathylAlign <- list(seqName=seqName,alignment=alignment, methPos=methPos, positionCGIRef=positionCGIRef, startEnd=startEnd, lengthRef=nchar(refSeq))

}


return (mathylAlign)

}


################################################################
################################################################
## cgMethFinder() finds the positions of TGs against CGs motivs	
## between  two aligned sequences. 							
################################################################
################################################################

cgMethFinder <- function (ref,str) {

ref <- paste(ref)
str <- paste(str)
cg_ref <- gregexpr("C-{0,}G",ref,perl=TRUE)	
cg_str <- gregexpr("C-{0,}G",str,perl=TRUE) 	
indexcg_ref <- cg_ref[[1]] [1: length(cg_ref[[1]])] 
indexcg_str <- cg_str[[1]] [1: length(cg_str[[1]])]	
intersect_cg_cg	<- as.numeric(indexcg_ref%in% indexcg_str)

return(intersect_cg_cg)

}


##################################################################
##################################################################
## Processing and extractiong CpG methylation data from .gff files
##################################################################
##################################################################

makeDataMethGFF <- function(dir,chr,start,end,meth_value) {

meth_value<- 1-meth_value
lengthRef <- c()
positionCGIRefTemp <- list()
methPosList <- list()
methPosListTemp <- list()
seqName <- c()
methData <- list()
positionCGIRef <-c()
files <- dir(path=dir,pattern=".gff")
startEnd <- matrix( ncol=2 , nrow=length(files))

for (i in 1:length(files)) {

	cat("processing ", paste(files[i]),'\n')
	gff <- read.table(paste(dir,'/',files[i],sep=""))
	gff <- data.frame(gff)
	colV1 <- gsub("chr","",gff$V1)
	colV1 <- gsub("X","30",colV1)
	colV1 <- gsub("Y","40",colV1)
	if (chr=='X') {chr<-"30"}
	if (chr=='Y') {chr<-"40"}
	cat("processing ", paste(files[i])," DONE",'\n')
	gff$V1 <- colV1  
	dataframe <- fn$sqldf("select *  from gff where V1=$chr and V4>$start and V4 < $end order by V4")
	cat("extracting data from ", paste(files[i])," DONE",'\n')
	seqName[i] <- files[i]
	startEnd[i,] <- c(dataframe$V4[1] , dataframe$V4[length(dataframe$V4)])
	positionCGIRefTemp[[i]] <- dataframe$V4
	whatIn <- positionCGIRefTemp[[i]] %in% positionCGIRef
	positionCGIRef <- sort(c( positionCGIRef ,positionCGIRefTemp[[i]][which(whatIn==FALSE)]))
	methPosList[[i]] <- dataframe$V6
	methPosListTemp[[i]] <- dataframe$V6
	methPosList[[i]] <- replace(methPosList[[i]],which(methPosListTemp[[i]]>=meth_value),0)
	methPosList[[i]] <- replace(methPosList[[i]],which(methPosListTemp[[i]]<meth_value),1)
}

n <- length(positionCGIRef)
methPos <- matrix(2,ncol=n , nrow=length(files))

for (i in 1:length(files)) {

	whatIn <- which( positionCGIRef %in% positionCGIRefTemp[[i]] )
	methPos[i,whatIn] <- methPosList[[i]]
}

lengthRef <- startEnd[,2][1]-startEnd[,1][1]+10
methData <- list(seqName=seqName, methPos=methPos, positionCGIRef=positionCGIRef, startEnd=startEnd, lengthRef=lengthRef)

return(methData)

}



