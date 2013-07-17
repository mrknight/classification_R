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


##############################################################
# Analysis 
##############################################################
readEvalData <- function(outputDir = "/home/dat/WORK/output/", dataFile1, dataFile2="") {
  evalData 	= read.csv(paste(outputDir,dataFile1,".csv", sep=""))
  if (dataFile2 != "") { # if another dataFile2 exists, then combine the 2 datasets
    evalData2 	= read.csv(paste(outputDir,dataFile2,".csv", sep=""))
    evalData	= rbind(evalData, evalData2)
    # TODO: quick and dirty to remove FXa from CSAR set
    FXa = read.table("/home/dat/WORK/output/FXa.txt")
    FXa = FXa[,1]
    evalData	= evalData[!rownames(evalData) %in% FXa,]    
  }
  return (evalData)
}

calcCor_CON <- function(evalData) {
# read the result of all CONventional scoring functions and calculate the correlation coefficients 
  scores 	= evalData
  scoreList	= colnames(scores)[-1]
  numScores	= length(colnames(scores)[-1])
  calcCorIndex <- function(iCol, data) {
    score	= abs(data[,iCol+1])
    scaleScores	= (score - min(score))/(max(score) - min(score))
    correl 	= c(cor(data[,1], scaleScores, method="pearson"), cor(data[,1], scaleScores, method="spearman"))
    #correl 	= c(abs(cor(data[,1], data[,iCol+1], method="pearson")), abs(cor(data[,1], data[,iCol+1], method="spearman")))
    names(correl) = c("pearson", "spearman")
    return (correl)
  }
  correl = t(sapply(c(1:numScores), calcCorIndex, scores))
  rownames(correl) = scoreList
  return (correl)
}

calcCor_CON(readEvalData(dataFile1="PDBbind_v2007_unprepared"))
calcCor_CON(readEvalData(dataFile1="PDBbind_v2007_prepared"))
calcCor_CON(readEvalData(dataFile1="PDBbind_v2012_unprepared"))
calcCor_CON(readEvalData(dataFile1="PDBbind_v2012_prepared"))
calcCor_CON(readEvalData(dataFile1="CSAR_set1_unprepared"))
calcCor_CON(readEvalData(dataFile1="CSAR_set2_unprepared"))

calcCor_CON(readEvalData(dataFile1="CSAR_set1_unprepared",dataFile2="CSAR_set2_unprepared"))


readEvalDataML <- function(outputDir = "/home/dat/WORK/output/ML-Scores/", method, dataFile1, dataFile2="") {
  evalData 	= read.csv(paste(outputDir,dataFile1,method,".csv", sep=""))
  if (dataFile2 != "") { # if another dataFile2 exists, then combine the 2 datasets
    evalData2 	= read.csv(paste(outputDir,dataFile2,method,".csv", sep=""))
    evalData	= rbind(evalData, evalData2)
    # TODO: quick and dirty to remove FXa from CSAR set
    FXa = read.table("/home/dat/WORK/output/FXa.txt")
    FXa = FXa[,1]
    evalData	= evalData[!evalData[,1] %in% FXa,]    
  }
  return (evalData)
}

calcCorrel <- function(method, dataFile1, dataFile2="", withPlot=FALSE, withTrain=FALSE) {
  if (withTrain) {
    dataFile1 = paste("train",dataFile1,sep="")
    dataFile2 = paste("train",dataFile2,sep="")
  }
  MLscores 	= readEvalDataML(method=method, dataFile1=dataFile1, dataFile2=dataFile2)
  scores	= abs(MLscores[,2])
  scaleScores	= (scores - min(scores))/(max(scores) - min(scores))
  #correl 	= c(cor(MLscores[,2], MLscores[,3], method="pearson"), cor(MLscores[,2], MLscores[,3], method="spearman"))
  correl 	= c(cor(scaleScores, MLscores[,3], method="pearson"), cor(scaleScores, MLscores[,3], method="spearman"))
  names(correl) = c("pearson", "spearman")
  if (withPlot) {
    x = MLscores[,3]
    y = scores
    plot(x, y, col = "blue", xlab="Measured binding affinity (PDBbind DB)", ylab=paste("Predicted binding affinity (ML-Score) with", method))
    prline = lm(y ~ x)	#Linear fit between predicted and measured binding affinity
    abline(prline)     
    new.x = data.frame( x = seq(from = range(x)[1]-2, to = range(x)[2]+2) )    
    #lines(x = new.x$x, y = predict( prline, new.x, interval = "confidence" )[ , "fit" ], col = "red" )
    lines(x = new.x$x, y = predict( prline, new.x, interval = "prediction" )[ ,"upr" ], col = "violet")
    lines(x = new.x$x, y = predict( prline, new.x, interval = "prediction" )[ ,"lwr" ], col = "violet")
    title(main=paste("R=",round(correl[1], digits=3), "on independent test set (",nrow(MLscores)," complexes)"))
  }
  return (correl)
}

for (i in c(1:11)) {
#  print(DATA_FNAMES[i])
#  cat(t(sapply(METHODS, calcCorrel, DATA_FNAMES[i])))
}

t(sapply(METHODS, calcCorrel, DATA_FNAMES[8], DATA_FNAMES[9]))
t(sapply(METHODS, calcCorrel, DATA_FNAMES[10], DATA_FNAMES[11]))

t(sapply(METHODS, calcCorrel, DATA_FNAMES[12]))
t(sapply(METHODS, calcCorrel, DATA_FNAMES[13]))
t(sapply(METHODS, calcCorrel, DATA_FNAMES[1], "", withPlot=TRUE))

calcCorrel(method="RF", dataFile1=DATA_FNAMES[2], withPlot=TRUE)
calcCorrel(method="RF", dataFile1=DATA_FNAMES[6], withPlot=TRUE)

mergeScoresData <- function(outputDir = "/home/dat/WORK/output/", dataFile1, dataFile2) {
  dataSet1 	= read.csv(paste(outputDir,dataFile1,".csv", sep=""))
  for (method in METHODS) {
    dataSet2 	= read.csv(paste(outputDir,"ML-Scores/",dataFile2,method,".csv", sep=""))
    mergeData	= merge(dataSet1, dataSet2[,1:2], by.x=1 , by.y=1)
    colnames(mergeData)[length(mergeData[1,])] = method
    dataSet1 	= mergeData
  }
  names(mergeData)[1] = "PDB"
  write.csv(mergeData, file=paste(outputDir,dataFile1,"_all.csv",sep=""),row.names = FALSE)
}

#mergeScoresData(dataFile1="PDBbind_v2007_unprepared", dataFile2=DATA_FNAMES[1])
#mergeScoresData(dataFile1="PDBbind_v2007_prepared", dataFile2=DATA_FNAMES[3])
#mergeScoresData(dataFile1="PDBbind_v2012_unprepared", dataFile2=DATA_FNAMES[2])
#mergeScoresData(dataFile1="PDBbind_v2012_prepared", dataFile2=DATA_FNAMES[4])

#mergeScoresData(dataFile1="PDBbind_v2007_unprepared_SFC", dataFile2=DATA_FNAMES[5])
#mergeScoresData(dataFile1="PDBbind_v2007_unprepared_SFC-RF", dataFile2=DATA_FNAMES[6])

#dataFile1 	= "PDBbind_v2007_unprepared"
#dataFile2 	= DATA_FNAMES[5]

for (method in METHODS) {
  RF_CSAR = SFCtrainData[!SFCtrainData[,1] %in% CSAR_SFC[,1],]
  predictMLscore(RF_CSAR, CSAR_SFC, method=method, outFile=paste("CSAR_SFC-", method, ".csv", sep=""))

  predictMLscore(trainData, CSARtestData, method=method, outFile=paste("CSAR_RF_SFC-", method, ".csv", sep=""))
}
