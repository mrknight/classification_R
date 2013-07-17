#################################################################
# Machine Learning Based Regression (all rights reserved)	#
# 	Author: Quoc Dat Nguyen					#
#	Usage:  \TODO  						#
#################################################################

library(randomForest)
library(e1071)
library(mboost)

source("ML_ReadData.r")
##################################################
# predict the scoring with choosen ML method     #
##################################################
predictMLscore <- function(trainData, testData, method = "BRT", targetIndex = 2, nameIndex = 1, outFile = "") {
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
	trainTarget = trainData[ iTrain, targetIndex]
	trainDesc  	= trainData[ iTrain, selectDesc]
	rownames(trainDesc)[1:nrow(trainDesc)] = as.vector(trainData[ iTrain, 1])

	iTest		= seq(1,nTestData, 1)
	testTarget  = testData[ iTest, targetIndex] 
	testDesc   	= testData[ iTest, selectDesc]
	rownames(testDesc)[1:nrow(testDesc)]   = as.vector(testData[ iTest, 1])

	##################################################
	# SELECTING ML WITH BEST INTERNAL VALIDATION     #
	##################################################
	if (method == "SVM") {
		nuValues	= seq(0.05, 1, by = 0.05)
		rmse_test_best	= 1e10
		for (nu in nuValues) {
			ML_Score 	= svm(trainTarget ~ ., data=trainDesc, nu = nu, type = 'nu-regression', na.action=na.omit) 
			testPred 	= predict(ML_Score, newdata = testDesc)   
			rmse_test = sqrt(mean(( testTarget - testPred)^2))
			if (rmse_test < rmse_test_best) {
				nu_best	= nu
				rmse_test_best	= rmse_test
				print(paste("nu_best=",nu_best, "rmse_test_best=",rmse_test_best))	
			}
			print(paste("nu=",nu))
		}
    epsValues	= seq(0, 1, by = 0.05)
    rmse_test_best_eps	= 1e10
    for (eps in epsValues) {
      ML_Score 	= svm(trainTarget ~ ., data=trainDesc, epsilon = eps, type = 'eps-regression', na.action=na.omit) 
      testPred 	= predict(ML_Score, newdata = testDesc)   
      rmse_test = sqrt(mean(( testTarget - testPred)^2))
      if (rmse_test < rmse_test_best_eps) {
	eps_best	= eps
	rmse_test_best_eps	= rmse_test
	print(paste("eps_best=",eps_best, "rmse_test_best=",rmse_test_best_eps))	
      }
      print(paste("eps=",eps))
    }
    if (rmse_test_best_eps < rmse_test_best) {
      ML_Score 	= svm(trainTarget ~ ., data=trainDesc, epsilon = eps_best, type = 'eps-regression', na.action=na.omit) 
    } else {
      ML_Score 	= svm(trainTarget ~ ., data=trainDesc, nu = nu_best, type = 'nu-regression', na.action=na.omit)     
    }
    print(ML_Score)
  } 
  else if (method == "RF") {
    set.seed(seed)
    rmse_OOB_best = 1e10		#dummy high value
    for (mtry in 2:ncol(trainDesc)) {
      RF_mtry = randomForest(trainTarget ~ ., data=trainDesc, ntree=500, mtry=mtry, na.action=na.omit)
      rmse_OOB = sqrt(mean(( RF_mtry$y - RF_mtry$predicted)^2)) 
      if (rmse_OOB < rmse_OOB_best) { 
	mbest 		= mtry
	rmse_OOB_best 	= rmse_OOB
	print(paste("mbest=",mbest, "rmse_OOB=",round(rmse_OOB,3)))
      }
      print(paste("mtry=",mtry))
    }
    ML_Score = randomForest(trainTarget ~ ., data=trainDesc, ntree=500, mtry=mbest, importance=TRUE, na.action=na.omit) #includes calculation of variable importance

    #####################################################
    # VARIABLE IMPORTANCE BY RF-SCORE		 	#
    #####################################################

    ML_Score$importance = importance(ML_Score, type=1) #only %incMSE
    varImpPlot(ML_Score, n.var=nrow(ML_Score$importance), main="Variable Importance")  #default is top 30 features
  }
  else if (method == "BRT") {
    ML_Score 	= blackboost(trainTarget ~ ., data=trainDesc)
  }
  ########################################################
  # ML-SCORE FIT OF TRAINING DATASET			 #
  ########################################################

  trainPred 	= predict(ML_Score, newdata = trainDesc) 
  rmse 		= format(sqrt(mean(( trainTarget - trainPred)^2)), digits=3)
  #print(rmse)  
  sdev 		= format(sd( trainTarget - trainPred ), digits=3)
  #print(sdev)
  fitpoints 	= na.omit(cbind(trainTarget, trainPred))
  fitcorr  	= format(round(cor(fitpoints[,1], fitpoints[,2]), digits=3))      #Pearson correlation coefficient 
  sprcorr   	= format(cor(rank(fitpoints[,1]), rank(fitpoints[,2])), digits=3) #Spearman rank correlation 

  plot(fitpoints[,1], fitpoints[,2], asp=1, xlab="Measured binding affinity (PDBbind DB)", ylab=paste("Predicted binding affinity (ML-Score) with", method))
  prline 	= lm(fitpoints[,2] ~ fitpoints[,1])						#Linear fit between predicted and measured binding affinity
  abline(prline) 
  title(main=paste("R=",fitcorr, "on training set (",nTrainData," complexes)"))
  grid()  

  #######################################################
  # ML-SCORE PREDICTION OF TEST DATASET	       		#
  #######################################################

  testPred 	= predict(ML_Score, newdata = testDesc) 
  rmse 		= format(sqrt(mean(( testTarget - testPred)^2)), digits=3) 
  sdev 		= format(sd( testTarget - testPred ), digits=3)
  fitpoints 	= na.omit(cbind(testTarget, testPred))
  fitcorr  	= format(round(cor(fitpoints[,1], fitpoints[,2]), digits=3))   	#Pearson correlation coefficient (R)
  sprcorr     	= format(cor(rank(fitpoints[,1]), rank(fitpoints[,2])), digits=3) #Spearman rank correlation 

  plot(fitpoints[,1], fitpoints[,2], asp=1, xlab="Measured binding affinity (PDBbind DB)", ylab=paste("Predicted binding affinity (ML-Score) with", method))
  prline 	= lm(fitpoints[,2] ~ fitpoints[,1])						#Linear fit between predicted and measured binding affinity
  abline(prline) 
  title(main=paste("R=",fitcorr, "on independent test set (",nTestData," complexes)"))
  grid()  

  #######################################################
  # WRITING OUT RF-SCORE PREDICTION OF TEST DATASET	#
  #######################################################
  if (outFile != "") {
    write.csv(cbind(testPred, testTarget), row.names=testData[seq(1,nTestData,1), nameIndex], file=outFile)
    write.csv(cbind(trainPred, trainTarget), row.names=trainData[seq(1,nTrainData,1), nameIndex], file=paste("train-",outFile,sep=""))
  }
  else {
#    fitpoints 	= na.omit(cbind(testTarget, testPred))
#    fitcorr  	= format(round(cor(fitpoints[,1], fitpoints[,2]), digits=3))   	#Pearson correlation coefficient (R)
#    sprcorr     	= format(cor(rank(fitpoints[,1]), rank(fitpoints[,2])), digits=3) #Spearman rank correlation 

  }
  return (fitcorr)
}


for (method in METHODS)
{
#  method="BRT"
  predictMLscore(RFtrainData_07, RFtestData_07, method=method, outFile=paste("RF-Score07-", method, ".csv", sep=""))
  predictMLscore(RFtrainData_12, RFtestData_12, method=method, outFile=paste("RF-Score12-", method, ".csv", sep=""))
  predictMLscore(RFtrainData_07proton, RFtestData_07proton, method=method, outFile=paste("RF-Score07proton-", method, ".csv", sep=""))
  predictMLscore(RFtrainData_12proton, RFtestData_12proton, method=method, outFile=paste("RF-Score12proton-", method, ".csv", sep=""))

  predictMLscore(SFCtrainData, SFCtestData, method=method, outFile=paste("SFC-Score07-", method, ".csv", sep=""))
  
  predictMLscore(trainData, testData, method=method, outFile=paste("RF_SFC-Score07-",method,".csv", sep=""))
  predictMLscore(trainData_proton, testData_proton, method=method, outFile=paste("RF_SFC-Score07proton-",method,".csv", sep=""))


  RF_CSAR = RFtrainData_07[!RFtrainData_07[,1] %in% CSARset1[,1],]
#  predictMLscore(RF_CSAR, CSARset1, method=method, outFile=paste("CSAR_RF-Score07_set1-", method, ".csv", sep=""))
  #predictMLscore(RFtrainData_07proton, CSARset1, method=method, outFile=paste("CSAR_RF-Score07_set1proton-", method, ".csv", sep=""))
  
  RF_CSAR = RFtrainData_07[!RFtrainData_07[,1] %in% CSARset2[,1],]
#  predictMLscore(RF_CSAR, CSARset2, method=method, outFile=paste("CSAR_RF-Score07_set2-", method, ".csv", sep=""))  
  #predictMLscore(RFtrainData_07proton, CSARset2, method=method, outFile=paste("CSAR_RF-Score07_set2proton-", method, ".csv", sep=""))
  
  RF_CSAR = RFtrainData_12[!RFtrainData_12[,1] %in% CSARset1[,1],]
#  predictMLscore(RF_CSAR, CSARset1, method=method, outFile=paste("CSAR_RF-Score12_set1-", method, ".csv", sep=""))
  #predictMLscore(RFtrainData_12proton, CSARset1, method=method, outFile=paste("CSAR_RF-Score12_set1proton-", method, ".csv", sep=""))
  
  RF_CSAR = RFtrainData_12[!RFtrainData_12[,1] %in% CSARset2[,1],]
#  predictMLscore(RF_CSAR, CSARset2, method=method, outFile=paste("CSAR_RF-Score12_set2-", method, ".csv", sep=""))
  #predictMLscore(RFtrainData_12proton, CSARset2, method=method, outFile=paste("CSAR_RF-Score12_set2proton-", method, ".csv", sep=""))

}


