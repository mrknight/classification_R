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
SFCtrainData 		= read.csv("/home/dat/WORK/output/descriptors/SFC/Training_set_fix.txt", na.strings=c(".", "NA", "", "?"))
SFCtestData 		= read.csv("/home/dat/WORK/output/descriptors/SFC/PDBbind_test_set.csv", na.strings=c(".", "NA", "", "?"))

RFtrainData_07 		= read.csv("/home/dat/WORK/output/descriptors/RF/PDBbind_training07.csv", na.strings=c(".", "NA", "", "?"))
RFtestData_07		= read.csv("/home/dat/WORK/output/descriptors/RF/PDBbind_test07.csv", na.strings=c(".", "NA", "", "?"))

RFtrainData_07proton 	= read.csv("/home/dat/WORK/output/descriptors/RF/PDBbind_training07_proton.csv", na.strings=c(".", "NA", "", "?"))
RFtestData_07proton		= read.csv("/home/dat/WORK/output/descriptors/RF/PDBbind_test07_proton.csv", na.strings=c(".", "NA", "", "?"))

RFtrainData_12proton 	= read.csv("/home/dat/WORK/output/descriptors/RF/PDBbind_training12_proton.csv", na.strings=c(".", "NA", "", "?"))
RFtestData_12proton		= read.csv("/home/dat/WORK/output/descriptors/RF/PDBbind_test12_proton.csv", na.strings=c(".", "NA", "", "?"))

RFtrainData_12nowater 	= read.csv("/home/dat/WORK/output/descriptors/RF/PDBbind_training12_nowater.csv", na.strings=c(".", "NA", "", "?"))
RFtestData_12nowater	= read.csv("/home/dat/WORK/output/descriptors/RF/PDBbind_test12_nowater.csv", na.strings=c(".", "NA", "", "?"))

RFtrainData_12 		= read.csv("/home/dat/WORK/output/descriptors/RF/PDBbind_training12.csv", na.strings=c(".", "NA", "", "?"))
RFtestData_12		= read.csv("/home/dat/WORK/output/descriptors/RF/PDBbind_test12.csv", na.strings=c(".", "NA", "", "?"))

CSARset1		= read.csv("/home/dat/WORK/output/descriptors/RF/CSARset1.csv", na.strings=c(".", "NA", "", "?"))
CSARset2		= read.csv("/home/dat/WORK/output/descriptors/RF/CSARset2.csv", na.strings=c(".", "NA", "", "?"))
CSARset_all		= rbind(CSARset1, CSARset2)
CSAR_SFC		= read.csv("/home/dat/WORK/output/descriptors/SFC/CSAR2010_test_set.csv", na.strings=c(".", "NA", "", "?"))
#####################################################
# DATA PRE-PROCESSING				 				#
#####################################################

SFCdesc		= c(69, 3:68)
SFCnameIndex= 1 # index of the PDB name
RFdesc	 	= c(82, 0, 1:3,6, 10:12, 15, 19:21, 24, 28:30, 33, 37:39, 42, 46:48, 51, 55:57, 60, 64:66, 69, 73:75, 78) 
RFdesc		= RFdesc + 1 # move the index to 1 pos
RFnameIndex	= 1 # index of the PDB name

# test data for PDBbind which merges RF and SFC descriptors together
trainData 	= mergeDesc2(RFtrainData_07, SFCtrainData, RFnameIndex, SFCnameIndex, RFdesc, SFCdesc)
testData 	= mergeDesc2(RFtestData_07, SFCtestData, RFnameIndex, SFCnameIndex, RFdesc, SFCdesc)

trainData_proton 	= mergeDesc2(RFtrainData_07proton, SFCtrainData, RFnameIndex, SFCnameIndex, RFdesc, SFCdesc)
testData_proton 	= mergeDesc2(RFtestData_07proton, SFCtestData, RFnameIndex, SFCnameIndex, RFdesc, SFCdesc)

# test data for CSAR which merges RF and SFC descriptors together
CSARtestData 	= mergeDesc2(CSARset_all, CSAR_SFC, RFnameIndex, SFCnameIndex, RFdesc, SFCdesc)

iSFCdesc		= c(69, 70, 3:68)

SFCtrainData	= SFCtrainData[,iSFCdesc]
SFCtestData		= SFCtestData[,iSFCdesc]

RFtrainData_07	= RFtrainData_07[,RFdesc]
RFtestData_07	= RFtestData_07[,RFdesc]

RFtrainData_07proton	= RFtrainData_07proton[,RFdesc]
RFtestData_07proton		= RFtestData_07proton[,RFdesc]

RFtrainData_12	= RFtrainData_12[,RFdesc]
RFtestData_12	= RFtestData_12[,RFdesc]

RFtrainData_12proton	= RFtrainData_12proton[,RFdesc]
RFtestData_12proton		= RFtestData_12proton[,RFdesc]

RFtrainData_12nowater	= RFtrainData_12nowater[,RFdesc]
RFtestData_12nowater	= RFtestData_12nowater[,RFdesc]

CSARset1	= CSARset1[,RFdesc]
CSARset2	= CSARset2[,RFdesc]
CSARset_all	= CSARset_all[,RFdesc]
# \TODO: quick and dirty to remove FXa from CSAR set
FXa 		= read.table("/home/dat/WORK/output/FXa.txt")
FXa 		= FXa[,1]
CSARset_all	= CSARset_all[!CSARset_all[,1] %in% FXa, ]

CSAR_SFC	= CSAR_SFC[,iSFCdesc]

##### DEFINE all the ML methods #####
METHODS		= c("BRT", "SVM-nu", "SVM-eps", "RF")
DATA_FNAMES	= c("RF-Score07-", "RF-Score12-", "RF-Score07proton-", "RF-Score12proton-", 
		    "SFC-Score07-",
		    "RF_SFC-Score07-", "RF_SFC-Score07proton-",
		    "CSAR_RF-Score07_set1-",
		    "CSAR_RF-Score07_set2-",
		    "CSAR_RF-Score12_set1-",
		    "CSAR_RF-Score12_set2-",
		    "CSAR_SFC-",
		    "CSAR_RF_SFC-")