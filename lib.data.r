#####################################################################
# library for manipulating data
#####################################################################

# check two data frames with a given index whether two data frames have common rows
checkCommonRows <- function(dataA, dataB, index) {
  if (length(intersect(dataA[,index], dataB[,index])) > 0) {
    return (TRUE)
  }
  else return (FALSE)
}

# create complements from 2 data frames with a given index row (= dataA - dataB). Pay attention to the order of two sets.
createRowsComplementsData <- function(dataA, dataB, index) {
  newData = dataA[!(dataA[,index] %in% dataB[,index]),]
  return (newData)
}

# read validation data from cvs file at outputDir, if another file exists then row binding two datas 
readAndCombiValiData <- function(outputDir = "/home/dat/WORK/output/", dataFile1, dataFile2="") {
  evalData 	= read.csv(paste(outputDir,dataFile1,".csv", sep=""))
  if (dataFile2 != "") { # if another dataFile2 exists, then combine rows of the 2 datasets
    evalData2 	= read.csv(paste(outputDir,dataFile2,".csv", sep=""))
    evalData	= rbind(evalData, evalData2)
    # \TODO: (fix) quick and dirty to remove FXa from CSAR set
    FXa = read.table("/home/dat/WORK/output/FXa.txt")
    FXa = FXa[,1]
    evalData	= evalData[!evalData[, 1] %in% FXa,]
  }
  return (evalData)
}

# merge validation data after column
mergeScoresData <- function(outputDir = "/home/dat/WORK/output/", dataSet1, dataFile2) {	
  for (method in METHODS) {
    dataSet2 	= read.csv(paste(outputDir,"ML-Scores/",dataFile2,method,".csv", sep=""))
    mergeData	= merge(dataSet1, dataSet2[,1:2], by.x=1 , by.y=1)
    colnames(mergeData)[length(mergeData[1,])] = method
    dataSet1 	= mergeData
  }
  names(mergeData)[1] = "PDB"
  write.csv(mergeData, file=paste(outputDir,dataFile2,"allScores.csv",sep=""),row.names = FALSE)
  return (mergeData)
}

# concatenate path for reading descriptors
concatPath <- function(path, prefix, cutoff, binsize, desc) {
  if (desc == "credo" && binsize == 12) {
    fullPath = paste(path,prefix,desc,"_c",cutoff,".csv", sep="")
  }
  else {
    fullPath = paste(path,prefix,desc,"_c",cutoff,"b",binsize,".csv", sep="")
  }
  return (fullPath)
}

# create combinations data from merging 3 data sets
#
# USAGE: 
#source("/Users/knight/dev/classification_R/lib.data.r")
#createMergeData(path = "/Users/knight/dev/classification_R/data/", prefix = "CASF13_training_")
#
createMergeData <- function(path, prefix, cutoff = 12, binsize = 12, descName = c("elements", "sybyl", "credo")) {
  data1 = read.csv(concatPath(path, prefix, cutoff, binsize, descName[1]), na.strings=c(".", "NA", "", "?"))
  data2 = read.csv(concatPath(path, prefix, cutoff, binsize, descName[2]), na.strings=c(".", "NA", "", "?"))
  data3 = read.csv(concatPath(path, prefix, cutoff, binsize, descName[3]), na.strings=c(".", "NA", "", "?"))
  # \IMPORTANT remote the col which contains pKd information 
  data2_removePkD = data2[,-1]
  data3_removePkD = data3[,-1]
  nameIndex1 = length(data1[1,])
  nameIndex2 = length(data2_removePkD[1,])
  nameIndex3 = length(data3_removePkD[1,])
  mergeData12   = merge(data1, data2_removePkD, by.x = nameIndex1, by.y = nameIndex2)
  mergeData13   = merge(data1, data3_removePkD, by.x = nameIndex1, by.y = nameIndex3)
  mergeData23   = merge(data2, data3_removePkD, by.x = (nameIndex2+1), by.y = nameIndex3)
  mergeData123  = merge(mergeData12, data3_removePkD, by.x = 1, by.y = nameIndex3)
  # write to file
  write.table(mergeData12, file = paste(path,prefix,descName[1],"-",descName[2],"_c",cutoff,"b",binsize, sep=""), sep = ",", row.names = FALSE)
  write.table(mergeData13, file = paste(path,prefix,descName[1],"-",descName[3],"_c",cutoff,"b",binsize, sep=""), sep = ",", row.names = FALSE)
  write.table(mergeData23, file = paste(path,prefix,descName[2],"-",descName[3],"_c",cutoff,"b",binsize, sep=""), sep = ",", row.names = FALSE)
  write.table(mergeData123, file = paste(path,prefix,descName[1],"-",descName[2],"-",descName[3],"_c",cutoff,"b",binsize, sep=""), sep = ",", row.names = FALSE)
}