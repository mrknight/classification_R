#####################################################################
# library for extern validation
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
	k 		= sum(yObs*yPred)/sum(yPred^2)
	r2_0	= 1 - sum( (yObs - k*yPred)^2 )/sum( (yObs - mean(yObs))^2 )
	r2		= cor(yObs, yPred)^2
	r2m		= r2 * ( 1 - sqrt(r2 - r2_0) )
	return (r2m)
}

# calc some validation metrics from a given observer and predicted values
# \param: yMean if given then is a vector of 2 values
calcValidationMetric <- function(yObs, yPred,  method = "lm", yMean = NULL) {
	yPred = abs(yPred)
	modelData 	= data.frame(cbind(yObs, yPred))
	# calculate the fitting from linear model
	if (method == "lm")
		model = lm(data = modelData, formula = yObs ~ yPred)
	else if (method == "glm")
		model = glm(data = modelData, formula = yPred ~ yObs)
		
	yFit 		= predict(model)
	# we could calculate the r2 yourself, or we could get the same metrics from summary of model in R
	sm 			= summary(model)
	obsMean 	= mean(yObs)
	r2.pred		= calcXSquare(yFit, yObs, obsMean) # sm$r.squared
#	r2.nofit	= calcXSquare(yPred, yObs, obsMean)
	r2.pearson	= cor(yPred, yObs)^2
# 	n 			= NROW(yFit)
	r2.adj		= sm$adj.r.squared #(r2*(n-1) - 1) / (n-2)
#	rmse		= calcRMSE(yPred, yObs)
	rse			= sm$sigma #calcRMSE(yFit, yObs, RSE=TRUE)
	r2m			= calcR2m(yObs = yObs, yPred = yPred)
#	r2m.reverse	= calcR2m(yObs = yObs, yPred = yPred, REVERSE = TRUE)
#	r2m.delta	= abs(r2m - r2m.reverse)
#	r2m.average	= (r2m + r2m.reverse) / 2
	
	r2.train	= NA
	r2.overall	= NA
	if (!is.null(yMean)) { # if yMean is given
		if (!is.na(yMean[1])) r2.train		= calcXSquare(yFit, yObs, yMean[1])
		if (!is.na(yMean[2])) r2.overall	= calcXSquare(yFit, yObs, yMean[2])
	}
	return (list (	"r2.pearson"	= r2.pearson, 
					"r2.pred" 		= r2.pred,
					"r2.train"		= r2.train, 
					"r2.overall"	= r2.overall, 
					"r2.adj"		= r2.adj, 
#					"r2.nofit"		= r2.nofit, 
					"r2m"			= r2m, 
#					"r2m.reverse"	= r2m.reverse, 
#					"r2m.delta"		= r2m.delta, 
#					"r2m.average"	= r2m.average, 
#					"rmse" 	= rmse,
					"rse"	= rse))
}

# same as calcValidationMetric but with other param order and return as a vector
getAllMetrics <- function(yPred, yObs,  yMean = NULL) {
	return ( unlist(calcValidationMetric(yObs = yObs, yPred = yPred, yMean = yMean, method = "lm")) )
}
