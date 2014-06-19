core = read.csv("/home/dat/WORK/dev/rfscore/bin/CASF13_test_credo.csv", na.strings=c(".", "NA", "", "?"))
refine = read.csv("/home/dat/WORK/dev/rfscore/bin/CASF13_refine_credo.csv", na.strings=c(".", "NA", "", "?"))
testData  = core
trainData = refine
nameIndex = length(testData[1,])
source("lib.validation.extern.r")

print(checkCommonRows(testData, trainData, nameIndex))
newTrainData = createRowsComplementsData(trainData, testData, nameIndex)
print(checkCommonRows(testData, newTrainData, nameIndex))

#newTrainData = trainData[!(trainData[,nameIndex] %in% testData[,nameIndex]),]
#length(intersect(testData[,nameIndex], newTrainData[,nameIndex]))

write.table(newTrainData, file = "/home/dat/WORK/dev/weka-3-6-11/CASF13_training_credo.csv", sep = ",", row.names = FALSE)
write.table(testData, file = "/home/dat/WORK/dev/weka-3-6-11/CASF13_test_credo.csv", sep = ",", row.names = FALSE)