#####################################################################
# calculate all external validation metrics
#####################################################################

METHODS		= c("BRT", "SVM-nu", "SVM-eps", "RF")
source('lib.validation.extern.r')

outputDir	= "/home/dat/WORK/output/unprepared/"
set1File	= "CSAR_set1_unprepared"
set2File	= "CSAR_set2_unprepared"
CSAR		= checkNA(readAndCombiValiData(outputDir, set1File, set2File))

#CSAR_set1 	= mergeScoresData(dataSet1 = read.csv(paste(outputDir,set1File,".csv", sep="")), dataFile2 = 'CSAR_PDBbind12nowater-set1_')
#CSAR_set2 	= mergeScoresData(dataSet1 = read.csv(paste(outputDir,set2File,".csv", sep="")), dataFile2 = 'CSAR_PDBbind12nowater-set2_')	

CSAR_all	= mergeScoresData(dataSet1 = CSAR, dataFile2 = 'CSAR_PDBbind12nowater_')

yObs		= CSAR_all[, 2] # experimental
result.RF   = t(apply(CSAR_all[, 3:NCOL(CSAR_all)], 2, getAllMetrics, yObs))

#yPred		= CSAR_set1[, 13] 
#result.RF   = calcValidationMetric(yObs, yPred, method = "lm")

CSAR_all	= mergeScoresData(dataSet1 = CSAR, dataFile2 = 'CSAR_PDBall12nowater_')

yObs			= CSAR_all[, 2] # experimental
result.RF_all   = t(apply(CSAR_all[, 3:NCOL(CSAR_all)], 2, getAllMetrics, yObs))

#testdata 	= data.frame(cbind(yObs, yPred))
#cvq2(testdata, yObs ~ yPred)
