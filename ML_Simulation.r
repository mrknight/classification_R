#################################################################
# Machine Learning Based Regression (all rights reserved)	#
# 	Author: Quoc Dat Nguyen					#
#	Usage:  /TODO  						#
#################################################################

library(randomForest)
library(e1071)
library(mboost)

##################################################
# predict the scoring with choosen ML method     #
##################################################
predictScoreSimu <- function(trainData, testData, method = "BRT", targetIndex = 2, nameIndex = 1) {
  # check if the training data intersect with the test data
  if (length(intersect(testData[,nameIndex], trainData[,nameIndex])) > 0) {
    print("##################################################")
    print("## WARNING: intersect between training and test ##")
    print("##################################################")  
    print(intersect(testData[,nameIndex], trainData[,nameIndex]))
  }
  # targetIndex: column of target in table
  nTrainData 	= nrow(trainData)		# number of pdb complexes for training
  nTestData 	= nrow(testData) 		# number of pdb complexes for testing
  seed		= 1

  iTrain		= seq(1,nTrainData, 1)
  nSample		= nTrainData			
  set.seed(seed)		
  selectDesc	= setdiff(c(1:ncol(trainData)), c(nameIndex, targetIndex))

  iTrain		= sample(iTrain)[1:nSample]		# shuffle selected complexes
  trainTarget 	= trainData[ iTrain, targetIndex]
  trainDesc  	= trainData[ iTrain, selectDesc]
  rownames(trainDesc)[1:nrow(trainDesc)] = as.vector(trainData[ iTrain, 1])

  iTest		= seq(1,nTestData, 1)
  testTarget  	= testData[ iTest, targetIndex] 
  testDesc   	= testData[ iTest, selectDesc]
  rownames(testDesc)[1:nrow(testDesc)]   = as.vector(testData[ iTest, 1])

  ##################################################
  # SELECTING ML WITH BEST INTERNAL VALIDATION     #
  ##################################################
  if (method == "SVM") {
    ML_Score 	= svm(trainTarget ~ ., data=trainDesc, type = 'nu-regression', na.action=na.omit) 
  } 
  else if (method == "RF") {
    ML_Score = randomForest(trainTarget ~ ., data=trainDesc, ntree=500, mtry=round(length(selectDesc)/3), importance=TRUE, na.action=na.omit) #includes calculation of variable importance
  }
  else if (method == "BRT") {
    ML_Score 	= blackboost(trainTarget ~ ., data=trainDesc)
  }
  #######################################################
  # ML-SCORE PREDICTION OF TEST DATASET	       		#
  #######################################################

  testPred 	= predict(ML_Score, newdata = testDesc) 
  rmse 		= format(sqrt(mean(( testTarget - testPred)^2)), digits=3) 
  sdev 		= format(sd( testTarget - testPred ), digits=3)
  fitpoints 	= na.omit(cbind(testTarget, testPred))
  fitcorr  	= format(round(cor(fitpoints[,1], fitpoints[,2]), digits=3))   	#Pearson correlation coefficient (R)
  sprcorr     	= format(cor(rank(fitpoints[,1]), rank(fitpoints[,2])), digits=3) #Spearman rank correlation 

  #######################################################
  # WRITING OUT RF-SCORE PREDICTION OF TEST DATASET	#
  #######################################################
  return (fitcorr)
}


RFtrainData_07 		= read.csv("/home/dat/WORK/output/descriptors/RF/PDBbind_training07.csv", na.strings=c(".", "NA", "", "?"))
RFtestData_07		= read.csv("/home/dat/WORK/output/descriptors/RF/PDBbind_test07.csv", na.strings=c(".", "NA", "", "?"))

RFtrainData_12 		= read.csv("/home/dat/WORK/output/descriptors/RF/PDBbind_training12.csv", na.strings=c(".", "NA", "", "?"))
RFtestData_12		= read.csv("/home/dat/WORK/output/descriptors/RF/PDBbind_test12.csv", na.strings=c(".", "NA", "", "?"))

##################################################
# DATA PRE-PROCESSING				 #
##################################################

RFdesc	 	= c(82, 0, 1:3,6, 10:12, 15, 19:21, 24, 28:30, 33, 37:39, 42, 46:48, 51, 55:57, 60, 64:66, 69, 73:75, 78) 
RFdesc		= RFdesc + 1# move the index to 1 pos
RFnameIndex	= 1 # index of the PDB name

RFtrainData_12	= RFtrainData_12[,RFdesc]
RFtestData_12	= RFtestData_12[,RFdesc]

RFtrainData_07	= RFtrainData_07[,RFdesc]
RFtestData_07	= RFtestData_07[,RFdesc]

percent		= seq(10, 90, by=10)

trainData 	= RFtrainData_07
testData	= RFtestData_07


predictPartData <- function(percentSize, method, trainData, testData) {    
  nTrainData 	= nrow(trainData)
  siz		= round((percentSize*nTrainData)/100)
  # set seed to "randomness"
  as.numeric(Sys.time())-> t; set.seed((t - floor(t)) * 1e8 -> seed);
  choosenPDB	= sample(nTrainData, size=siz)
  #print(paste(mean(choosenPDB), length(choosenPDB)))
  newTrainData	= trainData[ choosenPDB,]  
  return (predictScoreSimu(newTrainData, testData, method=method))  
}

predictWholeData <- function(dummy, method) {
  sapply(percent, predictPartData, method, trainData, testData)
}

resultBRT 	= sapply(1:100, predictWholeData, method="BRT")
resultSVM 	= sapply(1:100, predictWholeData, method="SVM")
resultRF 	= sapply(1:100, predictWholeData, method="RF")

resultBRT 	= apply(resultBRT, 1, as.numeric)
resultSVM 	= apply(resultSVM, 1, as.numeric)
resultRF 	= apply(resultRF, 1, as.numeric)

plotBRT		= apply(resultBRT, 2, mean, na.rm=TRUE)
plotSVM		= apply(resultSVM, 2, mean, na.rm=TRUE)
plotRF		= apply(resultRF, 2, mean, na.rm=TRUE)
plotX		= round((percent*nrow(trainData))/100)
ymin	= min(c(plotBRT, plotSVM, plotRF))
ymax	= max(c(plotBRT, plotSVM, plotRF))

plot(plotX, plotBRT ,type = "b", col="red", pch=0, lwd=2, ylim=c(ymin, ymax), ylab="Pearson correlation coefficient R", xlab="Number of training complexes", main="PDBbind_v2007")
lines(plotX, plotSVM, type = "b", col="blue", pch=1, lwd=2)
lines(plotX, plotRF, type = "b", col="green", pch=2, lwd=2)
legend("bottomright", c("BRT","SVM","RF"), lty=c(1,1,1), pch=c(0,1,2), lwd=c(2,2,2), col=c("red","blue","green"))
grid()
