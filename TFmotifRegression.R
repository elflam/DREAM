###ENCODE - Dream
###Continuous probabilities of TF binding for 50 bp sections across genome (divided into chromosomes)
###Emily Flam
###September 2016
### Current errors: division of chromosomes does not produce correct object type. Must subset manually
          # e.g. chr1 <- genome$chr1

setwd("~/Desktop/DREAM")

#load BSgenome.Hsapiens.UCSC.hg19
library("BSgenome.Hsapiens.UCSC.hg19")

#Load transcription factor weight matrix
pwm <- read.table("pwm-dream/ATF2_HUMAN.H10MO.B.pwm",header = T,sep="",
                  quote="", na.strings="", comment.char="", stringsAsFactors=F, col.names = c(1,2,3,4))
colnames(pwm) <- c("A","C", "G", "T")

#load probability chart
probs <- read.table("pwm-dream/ATF2_HUMAN.H10MO.B.thr", header=T, sep="", quote="", na.strings="", comment.char="",
                    stringsAsFactors=F)
colnames(probs) <- c("Weight", "Probability")

#Subset human genome (UCSC hg19) into chromosome objects
#not currently functional
genome <- BSgenome.Hsapiens.UCSC.hg19
for (i in 1:22){
  assign(paste("chr", i, sep=""), parse(text=paste("genome$chr", i, sep="")))
}

#manually set starting coordinates for sliding frame 
#to be manually increased by one each later 
bpStart <- 1
bpEnd <- 9

#initiallize large vector of weights
#n is used to fill SummaryVec
#m is used to fill tmpSummary and is constantly set back to 100, should never exceed 150
seqSummary <- vector()
n<-1
tmpSummary <- vector()
m <- 1
  
#begin loop for chr1
#collect weight and corresponding probability for each 9 bp frame 
for (i in 1:length(chr1) - 8){
  
  #isolate 9bp sequence 
  seq <- chr1[bpStart:bpEnd]
  seqVec <- as.vector(seq)
  names(seqVec) <- c(1:9)
  

  #loop through sequence of 9 and grab weights for each position
  #add the weights together to get one weight for the whole 9bp sequence
  #Disregard Ns in sequence 
  weight = 0
  for (i in 1:length(seq)){
    if (!(seqVec[i]) == "N"){
      weight <- weight + pwm[names(seqVec)[i],seqVec[i]]
    }
  }
  
  #Get corresponding p-value for sequence weight
  #set seqP to 1 if weight is 0
  if (weight == 0){
    seqP <- 1
  }else{
    seqP <- probs[which(abs(probs$Weight - weight) == min(abs(probs$Weight - weight))), 2]
  }
  
  #Get -log of p.val
  seqLogP <- (-log10(seqP))
  
  #Add p-value to temporary summary vector and label with numeral 1:42 to help divide later
  #There will be 50-length(sequence) sequences in each 50 bp chunk
  tmpSummary[m] <- seqLogP
  names(tmpSummary)[m] <- m
  
  #Combine chunks of 150 bp and empty temporary summary vector
  if (bpEnd %% 150 == 0){
    seqSummary[n] <- max(tmpSummary)
    names(tmpSummary)[n] <- n
    n <- n+1
    tmpSummary <- tmpSummary()
    m<-100
  }
  
  #Increase bpStart, End, and m by 1 to "slide" sequence frame by one bp
  bpStart <- bpStart + 1 
  bpEnd <- bpEnd + 1
  m <- m+1
  
}  

 

