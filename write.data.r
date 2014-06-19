#####################################################
# MERGE 2 data sets				 					#
#####################################################
mergeDesc2 <- function(data1, data2, nameIndex1, nameIndex2, selectCol1, selectCol2) {
	selectData1 	= data1[,selectCol1]
	selectData2 	= data2[,selectCol2]
	mergeData 		= merge(selectData1, selectData2, by.x = nameIndex1, by.y = nameIndex2)
	# \TODO: check descriptors and remove ZERO descriptors
	#for (i in 1:ncol(mergeData)) {
	#  print(sum(mergeData[,i]==0))
	#}
	return (mergeData)
}

# READ the csv data
RFtrainData_12 		= read.csv("/home/dat/WORK/output/descriptors/RF/PDBbind_training12.csv", na.strings=c(".", "NA", "", "?"))
RFtestData_12		= read.csv("/home/dat/WORK/output/descriptors/RF/PDBbind_test12.csv", na.strings=c(".", "NA", "", "?"))

#####################################################
# DATA PRE-PROCESSING				 				#
#####################################################

RFdesc	 	= c(82, 0, 1:3,6, 10:12, 15, 19:21, 24, 28:30, 33, 37:39, 42, 46:48, 51, 55:57, 60, 64:66, 69, 73:75, 78) 
RFdesc		= RFdesc + 1 # move the index to next 1 pos
RFnameIndex	= 1 # index of the PDB name

RFtrainData_12	= RFtrainData_12[,RFdesc]
RFtestData_12	= RFtestData_12[,RFdesc]
##### DEFINE all the ML methods #####

write.table(RFtrainData_12, "/home/dat/WORK/dev/weka/CASF12_training.csv", sep = ",", row.names = FALSE)
write.table(RFtestData_12, "/home/dat/WORK/dev/weka/CASF12_test.csv", sep = ",", row.names = FALSE)