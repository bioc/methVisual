####################################################
####################################################
## Analysis and Visualisation of Methylation Data ##
####################################################
####################################################


###########################
###########################
## absolute Methyl densenty
########################### 
###########################

plotAbsMethyl <- function(methData,real) {

absVec <- c()
methMatrix <- methData$methPos
methMatrix <- replace(methMatrix,which(methMatrix==2),0)
n <- dim(methMatrix)[2]

for(i in 1:n) {
	absVec[i] <- sum(methMatrix[,i])
}

if (missing(real) || real == FALSE  ) {
	x <- c(1:n)
}

else {
	x <- methData$positionCGIRef
}

if (missing(real) || real == TRUE  ) {
	plot(x,absVec,type="h", main="Absolute number of methylated CpGs", xlab= "position of CpG methylation on genomic sequence", ylab="absolute number of methylated CpG")
}

else {
	plot(x,absVec,type="h", main="Absolute number of methylated CpGs", xlab= "index of CpG methylation", ylab="absolute number of methylated CpG")
}

return(absVec)
gc(FALSE)

}


#################################
#################################
## Lollipops view of Methylation
#################################
#################################

MethLollipops <- function (methData) {

lRef <- methData$lengthRef
nrow <- lRef
ncol <- lRef
methMatrix <- findNonAligned (methData)
cgPos <- methData$positionCGIRef
ncgPos <- length(cgPos)
lengthData <- length(methData$seqName)
methyl <- c()
summeryMatrix_temp <- matrix(0, nrow = ncgPos, ncol = ncgPos)
summeryMatrix <- matrix(0, nrow = ncgPos, ncol = ncgPos)
cgPos <- c(1:ncgPos)

plot(c(0, ncgPos+1), c(0, lengthData + 1), yaxt = "n", type = "n", xlab = "index of CpG methylation", ylab = "index of clone sequences")

grid((ncgPos-1),NA,lwd=1)

methInAllSeq <- seq(0, 0, length.out = lRef)
vec1 <- seq(1, 1, length.out = ncgPos)
methInAllSeq[cgPos] <- vec1
y_methyl <- seq(lengthData + 1, lengthData + 1, length.out = ncgPos)
    
lines(cgPos, y_methyl, type = "p", col = "black", lwd = "2")

methInAllSeq <- seq(3, 3, length.out = lRef)

for (i in 1:lengthData) {

        seq_temp <- as.array(methMatrix[i, ])
        methInAllSeq_temp <- methInAllSeq
        methInAllSeq_temp[cgPos] <- seq_temp
        non_align <- which(methInAllSeq_temp == 2)
        methyl <- which(methInAllSeq_temp == 1)
        non_methyl <- which(methInAllSeq_temp == 0)
        y_non_align <- seq(i, i, length.out = length(non_align))
        y_methyl <- seq(i, i, length.out = length(methyl))
        y_non_methyl <- seq(i, i, length.out = length(non_methyl))
        axis(side = 2, at = c(1:(lengthData+1) ), labels = c(1:lengthData,"refSeq"), tick = TRUE, par(las = 2))
        lines(non_align, y_non_align, yaxt = "n", type = "p", col = NULL, lwd = "2")
        lines(methyl, y_methyl, type = "p", pch = 19, col = "black", lwd = "5")
        lines(non_methyl, y_non_methyl, type = "p", col = "black", lwd = "2")
}


LABEL_Y_AXIS <- c(1:(lengthData),"refSeq")
Experiment <- c(methData$seqName,"refernceSequence")
numberExp <- data.frame(LABEL_Y_AXIS = LABEL_Y_AXIS, Experiment= Experiment)

return(numberExp)
gc(FALSE)

}



###########################################
###########################################
## Cooccurence on Lollipops and Pie Diagram 
###########################################
###########################################


Cooccurrence <- function (methData,file,real,lolli) {

methPos <- methData$methPos
cgPos <- methData$positionCGIRef
n <- length(cgPos)
startRef <- min(methData$startEnd)
endRef <- max(methData$startEnd)
lengthData <- length(methData$seqName)
cgInAlign <- cgInAlign(methData)
whichCgExist <- which(cgInAlign > 0)

lengthOfX <- max(methData$positionCGIRef[  whichCgExist ]) -min(methData$positionCGIRef[  whichCgExist ])

if (lengthOfX >300) {

	real <- FALSE
	cat("span of aligned CpG is to big. The Cooccurrence image will be drown if parameter real=FALSE ",'\n')
}

CgExist <- c(1:length(whichCgExist))

if (missing(real) || real == FALSE  ) {

	methMatrix <- methData$methPos[,whichCgExist[1]:whichCgExist[length(CgExist)]]
}

else {

	methMatrix <- findNonAligned(methData)
}

if (missing(real) || real == FALSE  ) {

	cgPos <- CgExist
}

else {
	
	cgPos <- methData$positionCGIRef
}

if (missing(real) || real == FALSE  ) {

	lRef <-length(CgExist)
}

else {

	lRef <- methData$lengthRef
}

nrow <- lRef
ncol <- lRef
startMeth <- which(cgPos > startRef)[1]
length <- length(which(cgPos < endRef))
endMeth <- which(cgPos < endRef)[length]

if (missing(real)==TRUE || real == FALSE  ) {

	indexVec <- CgExist
}

else {

indexVec <- cgPos[startMeth:endMeth]
}

diff <- diff(indexVec)
minDiff <- min(diff)
summery <- matrixSNP(methData)
coo <- c()
relcoo <- c()
spanMeth <- indexVec[length(indexVec)] - indexVec[1]


if(missing(file)==FALSE) {
pdf(file, width = 20 * minDiff, height = lengthData, paper = "special", onefile = FALSE)
}

if(missing(file)==TRUE){
windows(n, lengthData,rescale="fixed")
}

plot(main ="Coocurrence Plot",c(indexVec[1] - spanMeth/50, indexVec[length(indexVec)] + spanMeth/50), c(0, lengthData + 2), yaxt = "n",xaxt="n", type = "n", xlab = "genomic position of CpG site" , ylab = "")

axis(1,at=cgPos)
coo <- c()

for (i in 1:(dim(summery)[1]) - 1) {

	coo[i] <- summery[i + 1, i]
        relcoo[i] <- coo[i]
}

k <- 0

for (i in 1:length(indexVec)) {

	k <- k + 1
        xlin <- c(indexVec[i], indexVec[i + 1])
        ylin <- c(lengthData + 1, lengthData + 1)
        xpostext <- (indexVec[i] + indexVec[i + 1])/2
        ypostext <- lengthData + 1.5
        lines(xlin, ylin, lwd = relcoo[i] * 10)
        result <- format(round(relcoo[i], 3))
        text(xpostext, ypostext, paste(result), col = "black", cex=20/n)
}

methInAllSeq <- seq(0, 0, length.out = lRef)

if(missing(real) || real==FALSE) {

vec1<- seq(1,1, length.out=length(CgExist))
}

else {
	vec1 <- seq(1, 1, length.out = n)
}

methInAllSeq[cgPos] <- vec1
methyl <- which(methInAllSeq == 1)
y_methyl <- seq(lengthData + 1, lengthData + 1, length.out = length(methyl))
k <- 1
lines(methyl, y_methyl, type = "p", col = "black", lwd = "2")
methInAllSeq <- seq(3, 3, length.out = lRef)


if (missing(lolli)==TRUE) {

	lolli <- 3
}



for (i in 1:lengthData) {

	seq_temp <- as.array(methMatrix[i, ])
        methInAllSeq_temp <- methInAllSeq
        methInAllSeq_temp[cgPos] <- seq_temp
        non_align <- which(methInAllSeq_temp == 2)
        methyl <- which(methInAllSeq_temp == 1)
        non_methyl <- which(methInAllSeq_temp == 0)
        y_non_align <- seq(k, k, length.out = length(non_align))
        y_methyl <- seq(k, k, length.out = length(methyl))
        y_non_methyl <- seq(k, k, length.out = length(non_methyl))
        axis(side = 2, at = 1:(lengthData + 1), labels = c(1:lengthData, ""), tick = TRUE, par(las = 2))
        lines(non_align, y_non_align, yaxt = "n", type = "p", col = NULL, lwd = "2")
        lines(methyl, y_methyl, type = "p", pch = 19, col = "black",cex=lolli)
        lines(non_methyl, y_non_methyl,pch=21, type = "p", col = "black",cex=lolli)
        k <- k + 1
}

absMeth <- c()
absNonMeth <- c()
methMatrixTemp <- replace(methData$methPos,which(methData$methPos==2),0)
relativeAmount <- c()
lengthData <- length(methData$seqName)
position <- which(cgPos %in% indexVec)
x <- indexVec
y <- rep(lengthData + 1, length(x))
n <- dim(methMatrix)[2]
nseq <- length(methData$seqName)

if(missing(real) || real==FALSE) {

	vecCgInAlign <- cgInAlign(methData)[whichCgExist]
}

else {

	vecCgInAlign <- cgInAlign(methData)
}


for (i in 1:length(x)) {

	if(missing(real) || real==FALSE) {

	        absMeth[i] <- sum(methMatrixTemp[, whichCgExist[i]])
	}

	else {

        	absMeth[i] <- sum(methMatrixTemp[, position[i]])
	}


	if(is.na(absMeth[i])==TRUE) { 

		absMeth[i]<-0
	}

	absNonMeth[i] <- nseq - absMeth[i]

	if(is.na(absNonMeth[i])==TRUE) {

		absNonMeth[i]<-0
	}

	relativeAmount[i] <- absMeth[i]/vecCgInAlign[position[i]]
       
	if (missing(real) || real == FALSE  ) {

		text(CgExist[i], lengthData + 2, paste(format(round(relativeAmount[i], 2))), col = "blue", cex=20/n)
	}

	else {

		text(indexVec[i], lengthData + 2, paste(format(round(relativeAmount[i], 2))), col = "blue", cex=20/n)
	}
 
}


vps <- baseViewports()

pushViewport(vps$inner, vps$figure, vps$plot, recording=FALSE)

grid.segments(x0 = unit(c(rep(0, length(x)), x), rep(c("npc", "native"), each = length(x))), x1 = unit(c(x, x), rep("native", 2 * length(x))), y0 = unit(c(y, rep(0, length(x))), rep(c("native", "npc"), each = length(x))), y1 = unit(c(y, y), rep("native", 2 * length(x))), gp = gpar(lty = "dashed", col = "grey"), draw = TRUE)

maxpiesize <- unit(8/length(x), "inches")
totals <- nseq
sizemult <- totals/max(totals)

for (i in 1:length(x)) {

	pushViewport(viewport(x = unit(x[i], "native"), y = unit(y[i], "native"), width = sizemult * maxpiesize, height = sizemult * maxpiesize), recording=FALSE)
        par(plt = gridPLT(), new = TRUE)

	if(absNonMeth[i]==0 && absMeth[i]==0) {

		pie(c(1,0), radius = 1,labels = rep("", 2))
popViewport()

	}

	else {

		pie(c(absNonMeth[i], absMeth[i]), radius = 1, labels = rep("",2))
	}

popViewport()


}

if(missing(file)==FALSE) {
dev.off()
}

LABEL_Y_AXIS <- c(1:(lengthData), "coocureceLolipop")
Experiment <- c(methData$seqName, "coocureceLolipop")
numberExp <- data.frame(LABEL_Y_AXIS = LABEL_Y_AXIS, Experiment = Experiment)

return(numberExp)
gc(FALSE)

}



###########
###########
## Heat-Map
###########
###########

heatMapMeth <- function (methData,file) {

cgPos <- methData$positionCGIRef
min <- min(methData$startEnd)
max <- max(methData$startEnd)
start <- which(cgPos >= min)[1]
length <- length(which(cgPos <= max))
end <- which(cgPos <= max)[length]
methMatrix <- replace(methData$methPos,which(methData$methPos==2),0)
methMatrix <- methMatrix[, start:end]
nrow <- dim(methMatrix)[1]
files <- paste(methData$seqName)

if(missing(file)==FALSE) {
pdf(file,onefile=FALSE)
}

heatmap(methMatrix,main="",scale="none", distfun=function(c){dist(c,method="binary")}) 

if(missing(file)==FALSE) {												
dev.off()
}
On_Plot_Nr <- c(1:nrow)
Experiment <- methData$seqName

numberExp <- data.frame(On_Plot_Nr= On_Plot_Nr, Experiment= Experiment)

return(numberExp)

}


###################################
###################################
## Summery Matrix for plotMatrixSNP
###################################
###################################

matrixSNP <-function (methData,correlation) {

cgPos <- methData$positionCGIRef
min <- min(methData$startEnd)
max <- max(methData$startEnd)
start <- which(cgPos >= min)[1]
length <- length(which(cgPos <= max))
end <- which(cgPos <= max)[length]
methMatrix <- replace(methData$methPos,which(methData$methPos==2),0)
methMatrix <- methMatrix[, start:end]
indexVec <-cgPos[start:end]
n <- length(indexVec)
summeryMatrix_temp <- matrix(0, nrow = n, ncol = n)
summeryMatrix <- matrix(0, nrow = n, ncol = n)
colnames(summeryMatrix) <- indexVec
rownames(summeryMatrix) <- indexVec
ncol <- length(methData$seqName)
cgInAlign <- cgInAlign(methData)
whichCgExist <- which(cgInAlign > 0)

if(missing(correlation)==TRUE || correlation==TRUE) {

	ncol <- dim(methMatrix)[2]
	summeryMatrix <- matrix(0, nrow = ncol, ncol = ncol)

	for (i in 1:ncol) {

		vec_temp <- methMatrix[,i]

	        for (y in 1:ncol) {

			if(sum(vec_temp)==0 || sum (methMatrix[,y])==0) {

				summeryMatrix[i,y] <- NA
			}

			else if(sum(vec_temp)==length(vec_temp) || sum (methMatrix[,y])== length(methMatrix[,y])){ 

				summeryMatrix[i,y] <- NA
			}

			else {

				summeryMatrix[i,y] <- round(cor(vec_temp,methMatrix[,y],use=,method="spearman" ),2)
			}


		}
	}

}   

else{

	for (i in 1:ncol) {

	        vec_temp <- methMatrix[i, ]

		for (y in 1:length(vec_temp)) {

			if (vec_temp[y] == 1) {

				summeryMatrix_temp[y, ] <- as.numeric(match(vec_temp, vec_temp[y], nomatch = 0))
			}
		
			else {

				summeryMatrix_temp[y, ] <- seq(0, 0, n)
			}
		}

		summeryMatrix <- summeryMatrix + summeryMatrix_temp
	}
}

return(summeryMatrix)

}


#############################################
#############################################
## drowing SNPPlot due to cooccurrence Matrix
#############################################
#############################################


plotMatrixSNP <- function (summeryMatrix,methData,file){

cgPos <- methData$positionCGIRef
min <- min(methData$startEnd)
max <- max(methData$startEnd)
start <- which(cgPos>=min)[1] 
length <- length (which(cgPos<=max))
end <- which(cgPos<=max)[length]
indexVec <-cgPos[start:end]
n <- dim(summeryMatrix)[1]
x <- c(1:nrow(summeryMatrix))
y <- c(1:nrow(summeryMatrix))	
col <- rev(gray(seq(0,1,1 / length(unique(as.vector(summeryMatrix))))))		
col.labels <-c( round( min(summeryMatrix,na.rm=TRUE),2), paste("cor",round(max(summeryMatrix, na.rm = TRUE),2)))
cex <- 0

for (i in 2:n){

	summeryMatrix[1:i-1,i]<-min(summeryMatrix)	
}	

if(missing(file)==FALSE) {
pdf(file,onefile=FALSE)
}

image(x=x,y=y,z=summeryMatrix, col = col,main="Distant cooccurrence plot", ylab="", xlab="",axes=FALSE, frame.plot=FALSE)

color.legend(0.5,length(x)/2,1,(length(x)/2)+(length(x)/4),col.labels,col ,gradient="y")

for(i in 1:n) {

	text(i,i,paste(indexVec[i]),col="red",cex=20/n)
}

for(i in 1:n) {

		for(x in i:n) {

			text(x,i,paste( round( summeryMatrix[x,i],3 )),col="blue",cex=10/n)
		}
}

if(missing(file)==FALSE) {
dev.off()
}

}

############
############
## using CA
############
############

methCA <- function (methData,file){


cgPos <- methData$positionCGIRef
min <- min(methData$startEnd)
max <- max(methData$startEnd)
start <- which(cgPos >= min)[1]
length <- length(which(cgPos <= max))
end <- which(cgPos <= max)[length]
methMatrix <- replace(methData$methPos,which(methData$methPos==2),0)
methMatrix <- methMatrix[, start:end]
nrow <- dim(methMatrix)[1]
ncol <- dim(methMatrix)[2]
cloneName <- c()
CpGi <- c()

for(i in 1:nrow) {

	cloneName[i] <- paste("clone",i,sep="")
}

for(i in 1:ncol) {

	CpGi[i] <- paste("CpG",i,sep="")
}

rownames(methMatrix) <- cloneName
colnames(methMatrix) <- CpGi
methMatrix <- methMatrix+0.1

if(missing(file)==FALSE) {
pdf(file,onefile=FALSE)
}

caMatrix <- ca(methMatrix)

plot(caMatrix)

if(missing(file)==FALSE) {
dev.off()
}

On_Plot_Nr <- c(1:nrow)
Experiment <- methData$seqName
numberExp <- data.frame(cloneName= cloneName, Experiment= Experiment)

return(numberExp)

}


###################################################################
###################################################################
## Fisher Test on two subsets of expariments over matched CpG sites
###################################################################
###################################################################

methFisherTest <- function(methData,set1,set2) {

cgPos <- methData$positionCGIRef
min <- min(methData$startEnd)
max <- max(methData$startEnd)
start <- which(cgPos >= min)[1]
length <- length(which(cgPos <= max))
end <- which(cgPos <= max)[length]
methMatrix <- methData$methPos
methMatrix <- replace(methMatrix,which(methMatrix==2),0)
methMatrix <- methMatrix[, start:end]


lengthCG <- dim(methMatrix)[2]
pVec <- c()
matSet1 <- methMatrix[set1,]
S1L <- dim(matSet1)[1]
matSet2 <- methMatrix[set2,]
S2L <- dim(matSet2)[1]

for(i in 1:lengthCG) {

	nMethS1 <- sum(matSet1[,i])
	nNonMethS1 <- S1L - nMethS1
	nMethS2 <- sum(matSet2[,i])
	nNonMethS2 <- S2L - nMethS2
	fisherMatTemp <- matrix(c(nMethS1,nMethS2,nNonMethS1,nNonMethS2),nrow=2,ncol=2)
	fisher <- fisher.test(fisherMatTemp)
	pVec[i] <- round(fisher$p.value,4)  
}

plot(c(1:lengthCG),pVec,main="Fisher's exact test",xlab="CpG index", ylab="p-value")

return(pVec)

}


######################
######################
## Mann Whitney U Test
######################
######################

methWhitneyUTest <- function (methData,set1,set2) {

cgPos <- methData$positionCGIRef
min <- min(methData$startEnd)
max <- max(methData$startEnd)
start <- which(cgPos >= min)[1]
length <- length(which(cgPos <= max))
end <- which(cgPos <= max)[length]
methMatrix <- methData$methPos
methMatrix <- replace(methMatrix,which(methMatrix==2),0)
methMatrix <- methMatrix[, start:end]


lengthCG <- dim(methMatrix)[2]
lset1 <- length(set1)
lset2 <- length(set2)
matSet1 <- methMatrix[set1,]
matSet2 <- methMatrix[set2,]
procentS1 <- c()
procentS2 <- c()

for(i in 1:lengthCG) {

	nMethS1 <- sum(matSet1[,i])
	procentS1[i] <- nMethS1/lset1
	nMethS2 <- sum(matSet2[,i])
	procentS2[i] <- nMethS2/lset2
}

Utest <- wilcox.test (procentS1,procentS2, paired=FALSE,exact=FALSE)

return(round(Utest$p.value,4))

}



#########################################################
#########################################################
## find non aligned start and end parts (used in the Lol) 
#########################################################
#########################################################

findNonAligned <- function (methData) {

startEnd <- methData$startEnd
positionCGIRef <- methData$positionCGIRef
methPos	<- replace(methData$methPos,which(methData$methPos==2),0)
vec3 <- c()

for (i in 1:dim(startEnd)[1]) {

start <- startEnd[i,1]
end <- startEnd[i,2]
vec1 <- which(positionCGIRef<start)
vec2 <- which(positionCGIRef>end)
vec3 <- c(vec1,vec2)
methPos[i,][vec3] <- seq(2,2,by=length(vec3))	
}
	
return(methPos)
gc(FALSE)

}


#####################################################
#####################################################
## finds sum of cg positions exists in ref and seq[i]
#####################################################
#####################################################

cgInAlign <- function(methData) {

lengthData <- length(methData$seqName)
posVec <- methData$positionCGIRef
startEnd <- methData$startEnd
sumWhichIn <- rep(0,0,length(posVec))

for (i in 1:lengthData) {

	whichIn	<- rep(0,0,length(posVec))
	start_one <- startEnd[i,1]
	end_one	<- startEnd[i,2]
	wichIn_Temp <- which(posVec>=start_one & posVec<=end_one)
	whichIn[wichIn_Temp] <- 1
	sumWhichIn <- sumWhichIn + whichIn
}

return(sumWhichIn)

}


