# Sleeve paried models for 2nd Time Point Only 
# input: "sleeve_seqtab_taxa_names.txt" - from DADA2


rm(list=ls())

setwd("/Users/whitneyware/IER_Sleeve/DADA2_analysis")

taxa <- c("phyla","genus" )
for( t in taxa )
{
  inFile <- ("sleeve_seqtab_taxa_names.txt")	
  myT <- read.table(inFile, sep="\t", header=TRUE, row.names = 1)	
  
  numCols <- ncol(myT)
  myColClasses <- c("character","character","character","numeric","numeric","numeric","character", rep("numeric", numCols-7))
  myT <-read.table(inFile,header=TRUE,sep="\t",colClasses=myColClasses)
  
  names <- colnames(myT)
  check <- cbind(names,myColClasses)
  # don't care about missing tumor for microbiome 
  #myT <- myT[!is.na(myT$Tumor.volume), ]
  myT = myT[ myT$Time.point..weeks. == 25,]
  
  groups <- sort(unique( myT$treatmentGroup))
  pValues <- vector()
  log10PValues <- vector()
  bugNames <- vector()
  xGroup <- vector()
  yGroup <- vector()
  meanXGroup <- vector()
  meanYGroup <- vector()
  xSampleSize <- vector()
  ySampleSize <- vector()
  
  index <- 1
  
  volume <- myT$Tumor.volume
  treatment <- myT$treatmentGroup
  
  endCol =which(names(myT) == "MDS1") -1
  for( i in 8:endCol)
  {
    bug <- myT[,i]
    
    if ( substr(names(myT)[i],1,3) != "MDS")
    {
      bug <- log( myT[,i] * 100000 + 1)	
    }
    
    if( sum( bug != 0) > nrow(myT) /4 ) 
    {
      for( x in 1:3 ) 
      {
        for( y in (x+1):4) 
        {
          bugNames[index] = names(myT)[i]
          
          xVals <- bug[myT$treatmentGroup==groups[x]];
          yVals <- bug[myT$treatmentGroup==groups[y]] 
          
          pValues[index] = t.test( xVals,yVals)$p.value
          xGroup[index] <- groups[x]
          yGroup[index] <- groups[y]
          meanXGroup[index] <- mean(xVals)
          meanYGroup[index] <- mean(yVals) 
          
          xSampleSize[index] <- length(xVals)
          ySampleSize[index] <- length(yVals)
          
          log10PValues[index] <- -log10(pValues[index])
          if( meanXGroup[index] < 	meanYGroup[index] ) 
            log10PValues[index] <- -log10PValues[index]
          
          index = index + 1				
        }
      }
    }
  }
  
  dFrame <- data.frame(bugNames,pValues,log10PValues,xGroup,yGroup,
                       xSampleSize, ySampleSize, meanXGroup,meanYGroup)
  
  dFrame <- dFrame[order(dFrame$pValues),]
  dFrame$pValuesAdjusted<- p.adjust( dFrame$pValues, method = "BH" )
  
  write.table(dFrame, file=paste0("pariedModelsBugs2ndTimePointOnly",t,".txt"), row.names=FALSE, sep="\t")	
}

