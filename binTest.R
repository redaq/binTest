####################################################################
### BinTest.R
### Test whether there is clustering of causal variants in the region
### Author: Dandi Qiao
### 1/6/2014
######################################################################

library(MASS)


binTest <- function(snpData, subjInfo, geno,  covName, family, quantilePCutoff=0.005, quantileRCutoff=0.001, permBinom=FALSE, numPerm=2000, threshold=0.05)

## snpData: Information of the variants. A data frame containing at least four columns of a map file, including the snp ID (column 1), chromosome number (column 2), genetic distance (column 3), and physical distance (column 4).

## subjInfo: Infomration of the subjects. A data frame containing the family ID (column 1 with column name FID), subject ID (column 2 with column name IID), phenotype (column 3 with specified column name), other covariates (column 4 and above, with specified column names).

## geno: genotypes of the subjects. A matrix containing the genotypes. The number of rows is the number of variants in the region, which should match the number of rows of snpData. The number of columns is the number of subjects in the dataset, which should match the number of rows of subjInfo.

## covName: The names of the covariates used in the association tests. They should correspond to the columns names in subjInfo
## link: The link function used in the glm function to calculate single-marker association p-values

## quantilePCutoff: The quantile cutoff of the single-marker association p-values used to select a subset of all the variants in the region. The default value is 0.5%.

## quantileRCutoff: The quantile cutoff of the number of neighboring variants used to select a subset of the neighboring variants around each variant for calculating the distances between pairs of variants. The default value is 0.1%.


## permBinom: Whether use binomial test to determine whether the clustering p-value is converged, then the parameter numPerm will be the number of permutations used in each permutation set. If permBinom is FALSE, the numPerm should be used to sepcify the total number of permutation to be used in the test. The default value is FALSE here.

## numPerm: see permBinom description.

## threshold: threshold used to check whether the clustering p-value has converged. The default is 0.05.


{	## check the input data
	nSNP <- nrow(snpData)
	nInd <- nrow(subjInfo)
	pheno <- subjInfo[,3]
	
	if(nrow(geno)!=nSNP | ncol(geno)!=nInd) { stop("The dimension of the genotype matrix does not match the snpData or subjInfo!")}
	
	## calculate the observed single-marker association p-values
	pvalue <- apply(geno, 1, getPvalue, subjInfo, covName, family)
	
	## cutoff
	sigP <- quantile(pvalue, quantilePCutoff)
	pvalueSig <- pvalue[pvalue < sigP]
	indexSig <- (1:nSNP)[pvalue < sigP]
	
	## permutation
	permMat <- list()
	indexSigPerm <- list()
	cat("Start permutation...\n")
	
	if(permBinom)
	{	countLoop <- 0
		finalPvalue <- -1
		ptest <- 0
		
		while(ptest < threshold)
		{	cat("Loop ", countLoop, "\n")
			for( B in 1:numPerm)
			{	cat(B, " ")
				phenoSim <- sample(pheno)
				pvaluePerm <- apply(geno, 1, getPvalue, subjInfo, covName, family)
				sigTPerm <- quantile(pvaluePerm, quantilePCutoff)
				indexSigPerm[[countLoop*numPerm+B]] <- (1:nSNP)[pvaluePerm < sigTPerm]
				permMat[[countLoop*numPerm+B]] <- pvaluePerm[pvaluePerm < sigTPerm]
			}
			
			## Test
			pvalue1 <- cal_bin(snpData, pvalueSig, indexSig, permMat, indexSigPerm, numPerm, quantileRCutoff)
			totalB <- numPerm + countLoop*numPerm*2
			cat("P value is ", pvalue1, " with ", totalB, " number of permutations. \n")
			
			## another test set
			testPermMat <- list()
			testIndexSigPerm <- list()
			cat("Test set ", countLoop, " \n")
			for( B in 1:numPerm)
			{
				cat(B, " ")
				phenoSim <- sample(pheno)
				pvaluePerm <- apply(geno, 1, getPvalue, subjInfo, covName, family)
				sigTPerm <- quantile(pvaluePerm, quantilePCutoff)
				testIndexSigPerm[[B]] <- (1:nSNP)[pvaluePerm < sigTPerm]
				testPermMat[[B]] <- pvaluePerm[pvaluePerm < sigTPerm] 
				
			} 
			pvalue2 <- cal_bin(snpData, pvalueSig, indexSig, testPermMat, testIndexSigPerm, numPerm, quantileRCutoff)
			tempX <- numPerm*pvalue2
			cat("Number of perm larger than observed in test set is ", tempX, "\n")
			ptest <- binom.test(tempX, numPerm, p=pvalue1, alternative=c("two.sided"))$p.value
			cat("Binomial p-value is ", ptest, " so we stop getting more permutations if p > ", threshold , ".\n")
			
			for( B in 1:numPerm)
			{
				indexSigPerm[[totalB+B]] <- testIndexSigPerm[[B]]
				permMat[[totalB+B]] <- testPermMat[[B]]
				
			}
			## Test
			finalP <- cal_bin(snpData, pvalueSig, indexSig, permMat, indexSigPerm, totalB+numPerm, quantileRCutoff)
			cat("The final p value is ", finalP, "\n")
			countLoop <- countLoop + 1
			
		}
		return(list(binTest_p = finalP, numberPerm=totalB+numPerm))
	}else
	{	for(B in 1:numPerm)
		{	cat(B, " ")
			phenoSim <- sample(pheno)
			pvaluePerm <- apply(geno, 1, getPvalue, subjInfo, covName, family)
			sigTPerm <- quantile(pvaluePerm, quantilePCutoff)
			indexSigPerm[[B]] <- (1:nSNP)[pvaluePerm < sigTPerm]
			permMat[[B]] <- pvaluePerm[pvaluePerm < sigTPerm]
		}
	
		## Test
		finalP <- cal_bin(snpData, pvalueSig, indexSig, permMat, indexSigPerm, numPerm, quantileRCutoff)
		return(list(binTest_p=finalP))
		
	}
	

}


getPvalue <- function(genotypes, subjInfo, covName, family)
## Calculating the single-marker association p-values
## genotypes: the genotypes of this marker, a vector of length nInd, which also match the number of rows in subjInfo.

## subjInfo: the information of the all the subjects.  A data frame containing the family ID (column 1 with column name FID), subject ID (column 2 with column name IID), phenotype (column 3 with specified column name), other covariates (column 4 and above, with specified column names).

## covName: the names of the covariants used in the test, they should be the column names of the corresponding columns in subjInfo

## family: the family used in glm function, a character string naming a family function, a family function or the result of a cal to a family function. 
 
{	genotypes <- as.numeric(genotypes)
	result <- glm(as.formula(paste("subjInfo[,3] ~ genotypes +", paste(covName, collapse="+"), sep="")), data=subjInfo, family=family)
	return (summary(result)$coef[2,4])
}

cal_bin <- function(snpData, pvalueSig, indexSig, permMat, indexSigPerm, numPerm, quantileRCutoff)
### This function is used to calculate the p-value in the clustering test
## snpData
## pvalueSig: the observed single-marker p-values that are smaller than the cutoff
## indexSig: the index of the variants with the single-marker p-values less than the cutoff
## permMat: the list of the p-values of the selected single-marker association tests, obtained under permutations, the length of the list equals the number of permutations
## indexSigPerm:: the list of the corresponding varint index of the selected single-marker p-values, obtained under permutation, the length of the list equals the number of permutations
## numPerm
## quantileRCutoff
{	numSnps <- ceiling(nrow(snpData)*quantileRCutoff)
	pos <- snpData[,4]
	Dnear <- calDistSig(pos[indexSig], pvalueSig, numSnps)
	
	chrDistPerm <- list()
	for( B in 1:numPerm)
	{
		chrDistPerm[[B]] <- calDistSig(pos[indexSigPerm[[B]]], permMat[[B]], numSnps)
	}
	chrDistPermUnlist <- unlist(chrDistPerm)
	
	quan <- 0.1*(1:9)
	binEdgeChr <- quantile(chrDistPermUnlist, quan, na.rm=TRUE)
	
	## under the null distribution
	percBMatChr <- matrix(0, 10, numPerm)
	
	for( B in 1:numPerm)
	{
		if(sum(is.na(chrDistPerm[[B]]))==0)
		{
			temp <- ecdf(chrDistPerm[[B]])
			cutoffv <- c(0, temp(binEdgeChr), 1)
			percBMatChr[,B] <- cutoffv[2:11]-cutoffv[1:10]
		}
	}
	Schr <- cov(t(percBMatChr))
	SChrinv <- ginv(Schr)
	statPermChr <- rep(0, numPerm)
	exp <- rep(0.1, 10)
	for( B in 1:numPerm)
	{
		if(sum(is.na(chrDistPerm[[B]]))==0)
		{
			statPermChr[B] <- t(percBMatChr[,B]-exp) %*% SChrinv %*% (percBMatChr[,B]-exp)
		}
	}
	
	DnearTempChr <- ecdf(Dnear)
	cutoff3 <- c(0, DnearTempChr(binEdgeChr),1)
	obsMatChr <- cutoff3[2:11]-cutoff3[1:10]
	
	statObsChr <- t(obsMatChr - exp) %*% SChrinv %*% (obsMatChr - exp)
	FnChr <- ecdf(statPermChr)
	binChr <- 1-FnChr(statObsChr)
	return(binChr)
}

calDistSig <- function(Dvector, Pvector, numAway)
## Function used to obtain all the distances
## Dvector: the physical distance in bps of the selected variants with small p-values than the cutoff value
## Pvector: the vector of all the selected p-values 
## numAway: The number of neighboring variants used in the calculation of distances
{
	numSigChr <- length(Dvector)
	numAway <- min(numSigChr, numAway)
	Dtemp <- rep(0, numAway*(numAway-1)/2+(numSigChr-numAway)*numAway)
	if(length(Dtemp)!=0)
	{
		for(snpindex in 2:numAway)
		{
			chrDistTemp <- (Dvector[snpindex]-Dvector[1:(snpindex-1)])
			chrPTemp <- sqrt(Pvector[snpindex]*Pvector[1:(snpindex-1)])
			Dtemp[((snpindex-1)*(snpindex-2)/2+1):(snpindex*(snpindex-1)/2)] <- chrDistTemp*chrPTemp
		}
		if(numAway!=numSigChr)
		{
			for(snpindex in (numAway+1):numSigChr)
			{
				chrDistTemp <- (Dvector[snpindex] - Dvector[(snpindex-numAway):(snpindex-1)])
				chrPTemp <- sqrt(Pvector[snpindex]*Pvector[(snpindex-numAway):(snpindex-1)])
				Dtemp[(numAway*(numAway-1)/2+numAway*(snpindex-numAway-1)+1):(numAway*(numAway-1)/2+numAway*(snpindex-numAway))] <- chrDistTemp*chrPTemp
			}
		}
		return (Dtemp)
	}else
	{
		stop("Error: length of Dtemp is 0")
	}
}