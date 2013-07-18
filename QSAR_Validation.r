#####################################################################
# calculate all external validation metrics
#####################################################################

# check NA value in each column of a given data frame, remove the column if NA exists
checkNA <- function(data) {
	colNA	= apply(data, 1, existNA)
	return (data[!colNA, ])
}
# check if a NA value exists in a given vector, return TRUE if exists else FALSE
existNA <- function(vec) {
	return (sum(is.na(vec)) > 0)
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
# calculate the squared correlation cofficient, depends on yMean and yiCalc
calcXSquare <- function(yiCalc, yi, yMean){
	return( 1 - sum( (yiCalc - yi) * (yiCalc - yi) ) /
				sum( (yi - yMean) * (yi - yMean) )
	)
}
# calculate the root mean squared error (RMSE) / or residual standard error (RSE)
calcRMSE <- function(yiCalc, yi, RSE=FALSE){
	n 	= NROW(yiCalc)
	if (RSE)
		n = n - 2
	return( sqrt( sum( (yiCalc - yi) * (yiCalc - yi) ) / n ) )
}
# calculate the r²_m metric and reverse r²_m (r'²_m) metric by interchange the axes
calcR2m <- function(yObs, yPred, REVERSE=FALSE) {	
	if (REVERSE) { # to calculate the reverse need to interchange the yObs and yPred value
		tmp 	= yObs
		yObs	= yPred
		yPred	= tmp
	}
	k 	= sum(yObs*yPred)/(sum(yPred)^2)
	r2_0	= 1 - ( sum((yObs - k*yPred)^2) )/( (yObs - mean(yObs))^2 )
	r2		= cor(yObs, yPred)
	r2m		= r2 * ( 1 - sqrt(r2 - r2_0) )
	return (r2m)
}
# \TODO: merge validation data after column
mergeScoresData <- function(outputDir = "/home/dat/WORK/output/", dataSet1, dataFile2) {	
	for (method in METHODS) {
		dataSet2 	= read.csv(paste(outputDir,"ML-Scores/",dataFile2,method,".csv", sep=""))
		mergeData	= merge(dataSet1, dataSet2[,1:2], by.x=1 , by.y=1)
		colnames(mergeData)[length(mergeData[1,])] = method
		dataSet1 	= mergeData
	}
	names(mergeData)[1] = "PDB"
	write.csv(mergeData, file=paste(outputDir,dataFile2,"all.csv",sep=""),row.names = FALSE)
}


# \TODO: explain
calcValidationMetric <- function(yObs, yPred,  method = "lm") {
	modelData 	= data.frame(cbind(yObs, yPred))
	# calculate the fitting from linear model
	if (method == "lm")
		model = lm(data = modelData, formula = yObs ~ yPred)
	else if (method == "glm")
		model = glm(data = modelData, formula = yPred ~ yObs)
		
	yFit 		= predict(model)
	# we could calculate the r2 yourself, or we could get the same metrics from summary of model in R
	sm = summary(model)
	obsMean 	= mean(yObs)	
	r2			= sm$r.squared #calcXSquare(yFit, yObs, obsMean)
	r2.nofit	= calcXSquare(yPred, yObs, obsMean)
	r2.pearson	= cor(yPred, yObs)^2
	# n 			= NROW(yFit)
	r2.adj		= sm$adj.r.squared #(r2*(n-1) - 1) / (n-2)
	rmse		= calcRMSE(yFit,  yObs)
	rse			= sm$sigma #calcRMSE(yFit,  yObs, RSE=TRUE)
	r2m			= calcR2m(yObs = yObs, yPred = yPred)
	r2m.reverse	= calcR2m(yObs = yObs, yPred = yPred, REVERSE = TRUE)
	
	return (list (	"r2" 	= r2,
					"r2.adj"= r2.adj, 
					"r2.nofit"		= r2.nofit, 
					"r2.pearson"	= r2.pearson, 
					"rmse" 	= rmse,
					"r2m"	= r2m, 
					"r2m.reverse"	= r2m.reverse, 
					"rse"	= rse))
}

outputDir	= "/home/dat/WORK/output/unprepared/"
set1File	= "CSAR_set1_unprepared"
set2File	= "CSAR_set2_unprepared"
CSAR		= checkNA(readAndCombiValiData(outputDir, set1File, set2File))

yObs		= CSAR[, 2] # experimental
yPred		= CSAR[, 3]

test 		= calcValidationMetric(yObs, yPred, method = "lm")
#testdata 	= data.frame(cbind(yObs, yPred))
#cvq2(testdata, yObs ~ yPred)

CSAR_all	= mergeScoresData(dataSet1 = CSAR, dataFile2 = 'CSAR_PDBbind12nowater_')
CSAR_set1 	= mergeScoresData(dataSet1 = read.csv(paste(outputDir,set1File,".csv", sep="")), dataFile2 = 'CSAR_PDBbind12nowater-set1_')
CSAR_set2 	= mergeScoresData(dataSet1 = read.csv(paste(outputDir,set2File,".csv", sep="")), dataFile2 = 'CSAR_PDBbind12nowater-set2_')	

