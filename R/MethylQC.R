#############################################
#############################################
###																			
###  Quality Control (QC) for mehtylation				
###									
#############################################
#############################################

MethylQC <- function(refSeq, methFileDataFrame,makeChange,identity,conversion) {

seqName	<- c()
alignment <- c()
score_normal <- c()
score_complement <- c() 
score_reverse <- c()
score_reverseComplement <- c()
conversion_ratio <- c()
seq_identity <- c()
message <- c()
FILE <- c()
PATH <- c()
files <- paste(methFileDataFrame$FILE) 
paths <- paste(methFileDataFrame$PATH)
mat <- matrix(0,length(DNA_ALPHABET),length(DNA_ALPHABET))
mat[1:4,1:4] <- c(1,0,0,0,0,1,0,1,0,0,1,0,0,0,0,1)
rownames(mat) <- colnames(mat) <- DNA_ALPHABET[1:length(DNA_ALPHABET)]

if (missing(identity)) {

	identity <- 80
}

if (missing(conversion)) {

	conversion <- 90
}

pf <-0	

	for (i in 1: length (paths)) {
		pf <- pf+1	
		seqFileTemp <- paste( paths[i],'/',files[i], sep="")
		seqTemp	<- toupper(scan(seqFileTemp,what='character', sep="",quiet=TRUE,skip=1))
		seqName[i] <- files [i]
		cat("checking", seqName [i],'\n' )

		align_normal <- pairwiseAlignment(refSeq, seqTemp,substitutionMatrix=mat,type="local-global")
		score_normal[i]	<- score(align_normal) 
		alignOfnormal <- paste(seqTemp)
		
		align_complement <- pairwiseAlignment(refSeq, paste(complement(DNAString(seqTemp))), type="local-global", substitutionMatrix=mat)
		score_complement[i] <- score(align_complement)
		alignOfcomplement <- paste(complement(DNAString(seqTemp)))

		align_reverse <- pairwiseAlignment(refSeq, paste(reverse(DNAString(seqTemp))), type="local-global", substitutionMatrix=mat)
		score_reverse[i] <- score(align_reverse)
		alignOfreverse <- paste(reverse(DNAString(seqTemp)))

		align_reverseComplement	<- pairwiseAlignment(refSeq, paste(reverseComplement(DNAString(seqTemp))), type="local-global", substitutionMatrix=mat)
		score_reverseComplement[i] <- score(align_reverseComplement)
		alignOfreverseComplement <- paste(reverseComplement(DNAString(seqTemp)))

		alignment <- list(align_normal, align_complement, align_reverse, align_reverseComplement)
		scoreOfAligns <- c(score_normal[i],score_complement[i],score_reverse[i],score_reverseComplement[i])
		alignOfAlign <- c(alignOfnormal,alignOfcomplement,alignOfreverse,alignOfreverseComplement)

		allOfMessages <- c(": no changes are needed" , ": check if complemented" , ": chack if reversed" ,  ": check if reversed&complemented")
		maxAlign <- which (scoreOfAligns==max(scoreOfAligns))[1] 

			if (max(scoreOfAligns)- score_normal[i] > score_normal[i]/10) { 
			
				message[i]<- allOfMessages[maxAlign]
			} 
			
			else { 

			message[i] <-": no alignments changes are needed"
			} 

		changeAlign <- paste(alignOfAlign[maxAlign])  # tried also with paste

		if (missing(makeChange) || makeChange==TRUE) { 
			write(changeAlign,paste(paths[i],"/","QC_",files[i],sep=""))
		}

#################################
## first QC procedure: conversion
#################################

lengthRef <- nchar(refSeq)
align <- alignment[[maxAlign]]
startRef <- start(pattern(align))
endRef <- end(pattern(align))
newRef 	<- paste(substr(refSeq,1,startRef-1),pattern(align), substr(refSeq,endRef+1,lengthRef))
newSeq 	<- paste(paste(rep('-',startRef-1),collapse=""), subject(align), paste(rep('-',lengthRef-endRef)))	
c_ref_inf <- gregexpr("C",newRef,perl=TRUE)
c_seq_inf <- gregexpr("C",newSeq,perl=TRUE)
c_ref <- c_ref_inf [[1]]  [1: length(c_ref_inf[[1]])] 
c_seq <- c_seq_inf [[1]]  [1: length(c_seq_inf[[1]])]
n_c_ref <- length(c_ref)
n_c_seq	<- length(c_seq)
cg_ref_inf <- gregexpr("C-{0,}G",newRef,perl=TRUE)	
cg_seq_inf <- gregexpr("C-{0,}G",newSeq,perl=TRUE)			
cg_ref <- cg_ref_inf [[1]]  [1: length(cg_ref_inf[[1]])] 
cg_seq <- cg_seq_inf [[1]]  [1: length(cg_seq_inf[[1]])]
n_cg_ref <- length (cg_ref)
n_cg_seq <- length (cg_seq)	
n_c_nonCG_ref <- n_c_ref - n_cg_ref
n_c_nonCG_seq <- n_c_seq - n_cg_seq

if(n_c_nonCG_ref > 0) {

	conversion_ratio[i] <- ((n_c_nonCG_ref - n_c_nonCG_seq) / n_c_nonCG_ref)*100
}

else {

	conversion_ratio[i] <- 0
}

if (conversion_ratio[i] < conversion) {

	cat (seqName [i]," BAD CONVERSION RATE: ", conversion_ratio[i],"% !!!",'\n' )
}

###########################################
# seconde QC procedure: sequence identity 
###########################################

genomicAlign <- strsplit(coversionGenom(paste(pattern(align))),NULL)[[1]]
seqAlign <- strsplit( (paste(subject(align))),NULL ) [[1]]
newSeq <- strsplit(newSeq,NULL)[[1]]
lengthAlign <- length(seqAlign)
sumIdent <- sum(genomicAlign==seqAlign)
seq_identity[i]	<- sumIdent/lengthAlign*100

if (seq_identity[i] < identity) {

	cat (seqName [i]," BAD SEQUENCE IDENTITY: ", seq_identity[i],"% !!!",'\n' )
}

if(missing(makeChange) || makeChange==TRUE) {

	if (seq_identity[i] <identity || conversion_ratio[i] <conversion  ) {
		
		pf <-pf-1
	}

	else {
		FILE[pf] <- paste("QC_",files[i],sep="") 
		PATH[pf] <- paths[i]
	}
}

}  

if (missing(makeChange) || makeChange==TRUE) {

	newFilePath <- data.frame(FILE,PATH)
	newFilePath$FILE <- FILE
	newFilePath$PATH <- PATH
	write.table(newFilePath,file=paste(paths[i],"/","QC_methylData",Sys.Date(),sep=""),row.names = FALSE,quote=FALSE,sep='\t')
}

else { 

	cat ('you choosed not to use the autochange function, please be sure that the sequences are right aligned and in a good quality!','\n')
}

QCINFO	<- data.frame(seqName, message, score_normal, score_complement, score_reverse, score_reverseComplement, conversion_ratio, seq_identity)
QCINFO$seqName <- seqName
QCINFO$message <- message
QCINFO$score_normal <- score_normal
QCINFO$score_complement <- score_complement
QCINFO$score_reverse <- score_reverse
QCINFO$score_reverseComplement <- score_reverseComplement
QCINFO$conversion_ratio	<- conversion_ratio
QCINFO$seq_identity <- seq_identity
	
save(QCINFO, file = paste(paths[i], "/QCINFO.Rdata", sep = ""))
write.table(QCINFO, file = paste(paths[i], "/QCINFO.txt",sep = "" ),sep = "\t")

if (missing(makeChange) || makeChange == TRUE) {

	cat(paste("The new experiment data files after QC (QC_*) and the QCINFO.Rdata are under",paths[i]), "\n","please use ",paste ("load(",'"',paths[i],"/QCINFO.Rdata",'"',")",sep=""), " in order to view QC information","\n" )
}

else { 
	newFilePath <- methFileDataFrame
}

return(newFilePath)

}


###################################################
###################################################
## conversion of genomic sequence (C non CpG into T)
###################################################
###################################################

coversionGenom <- function(genomicSeq) {

cg_find <- gregexpr("C-{0,}G",genomicSeq)
c_find <- gregexpr("C",genomicSeq)
genomicSeq_bisulfite <- strsplit(genomicSeq,NULL)[[1]]
index_cg <- cg_find[[1]][1:length(cg_find[[1]])]
index_c <- c_find[[1]][1:length(c_find[[1]])]
is_cg <- as.numeric(index_c  %in% index_cg)

for (i in 1:length(is_cg)){

	if (is_cg[i] == 1) {

		genomicSeq_bisulfite[index_c[i]]<- 'C'
	}

	else {

		genomicSeq_bisulfite[index_c[i]]<- 'T'
	}

}		

genomicSeq_bisulfite <- paste(genomicSeq_bisulfite,collapse = "")

return(genomicSeq_bisulfite)

}

