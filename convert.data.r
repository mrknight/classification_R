core = read.csv("/Users/knight/dev/classification_R/weka-3-6-11/CASF13_core_elements_c12b12.csv", na.strings=c(".", "NA", "", "?"))
refine = read.csv("/Users/knight/dev/classification_R/weka-3-6-11/CASF13_refined_elements_c12b12.csv", na.strings=c(".", "NA", "", "?"))
testData  = core
trainData = refine
nameIndex = length(testData[1,])
source("/Users/knight/dev/classification_R/lib.validation.extern.r")

print(checkCommonRows(testData, trainData, nameIndex))
newTrainData = createRowsComplementsData(trainData, testData, nameIndex)
print(checkCommonRows(testData, newTrainData, nameIndex))

#newTrainData = trainData[!(trainData[,nameIndex] %in% testData[,nameIndex]),]
#length(intersect(testData[,nameIndex], newTrainData[,nameIndex]))

write.table(newTrainData, file = "/Users/knight/dev/classification_R/weka-3-6-11/CASF13_training_elements_c12b12.csv", sep = ",", row.names = FALSE)
write.table(testData, file = "/Users/knight/dev/classification_R/weka-3-6-11/CASF13_test_elements_c12b12.csv", sep = ",", row.names = FALSE)