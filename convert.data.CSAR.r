core1 = read.csv("/home/dat/WORK/dev/weka-3-6-11/CSARset1_elements_c12b12.csv", na.strings=c(".", "NA", "", "?"))
core2 = read.csv("/home/dat/WORK/dev/weka-3-6-11/CSARset2_elements_c12b12.csv", na.strings=c(".", "NA", "", "?"))
refine = read.csv("/home/dat/WORK/dev/weka-3-6-11/CASF13_refined_elements_c12b12.csv", na.strings=c(".", "NA", "", "?"))
core  = rbind(core1, core2)
# \TODO: quick and dirty to remove FXa from CSAR set
FXa 		= read.table("/home/dat/WORK/output/FXa.txt")
FXa 		= FXa[,1]
core		= core[!core[,1] %in% FXa, ]

testData  = core
trainData = refine
nameIndex = length(testData[1,])
source("lib.validation.extern.r")

print(checkCommonRows(trainData, testData, nameIndex))
newTrainData = createRowsComplementsData(trainData, testData, nameIndex)
print(checkCommonRows(testData, newTrainData, nameIndex))

#newTrainData = trainData[!(trainData[,nameIndex] %in% testData[,nameIndex]),]
#length(intersect(testData[,nameIndex], newTrainData[,nameIndex]))

write.table(newTrainData, file = "/home/dat/WORK/dev/weka-3-6-11/CSAR13_training_elements_c12b12.csv", sep = ",", row.names = FALSE)
write.table(testData, file = "/home/dat/WORK/dev/weka-3-6-11/CSAR13_test_elements_c12b12.csv", sep = ",", row.names = FALSE)