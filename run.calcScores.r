#################################################################
# Machine Learning Based Regression (all rights reserved)	#
# 	Author: Quoc Dat Nguyen					#
#	Usage:  \TODO  						#
#################################################################

library(randomForest)
library(e1071)
library(mboost)

source("init.readData.r")
source("lib.prediction.r")

#RFallData_12nowater = rbind(RFtrainData_12, RFtestData_12)
METHODS		= c("RF")
for (method in METHODS) 
{
	predictMLscore(trainData=RFtrainData_12, testData=RFtestData_12, method=method, outFile=paste("CASF12_", method, ".csv", sep=""))

#	RF_CSAR = RFtrainData_12nowater[!RFtrainData_12nowater[,1] %in% CSARset1[,1],]
#	predictMLscore(RF_CSAR, CSARset1, method=method, outFile=paste("CSAR_PDBbind12nowater-set1_", method, ".csv", sep=""))
	# give warning because of overlapping in training and test
#	predictMLscore(RFtrainData_12nowater, CSARset1, method=method, outFile=paste("CSAR_RF-Score12_set1_test-", method, ".csv", sep="")) 
	
#	RF_CSAR = RFtrainData_12nowater[!RFtrainData_12nowater[,1] %in% CSARset2[,1],]
#	predictMLscore(RF_CSAR, CSARset2, method=method, outFile=paste("CSAR_PDBbind12nowater-set2_", method, ".csv", sep=""))
	# give warning because of overlapping in training and test
#	predictMLscore(RFtrainData_12nowater, CSARset2, method=method, outFile=paste("CSAR_RF-Score12_set2_test-", method, ".csv", sep=""))
	
#	RF_CSAR = RFtrainData_12nowater[!RFtrainData_12nowater[,1] %in% CSARset_all[,1],]
#	predictMLscore(RF_CSAR, CSARset_all, method=method, outFile=paste("CSAR_PDBbind12nowater_", method, ".csv", sep=""))	

#	RF_CSAR = RFallData_12nowater[!RFallData_12nowater[,1] %in% CSARset_all[,1],]
#	predictMLscore(RF_CSAR, CSARset_all, method=method, outFile=paste("CSAR_PDBall12nowater_", method, ".csv", sep=""))	
}


