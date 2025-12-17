#### test binTest.R

geno <- read.table("sample.geno")
snpData <- read.table("sample.snp")
subjInfo <- (read.table("sample.ind"))
covPC <- read.table("samplePC.txt")

subjInfo <- cbind(subjInfo[,1], subjInfo[,1], subjInfo[,3], subjInfo[,2], covPC[,1:5])
colnames(subjInfo) <- c("FID", "IID", "Affect", "Sex", "PC1", "PC2", "PC3", "PC4", "PC5")


covName <- c("Sex", "PC1", "PC2")
family=binomial(link="logit")
quantilePCutoff=0.005
quantileRCutoff=0.001
permBinom=FALSE
numPerm=3
threshold=0.05
binTest(snpData, subjInfo, geno, covName, family, quantilePCutoff=0.1, quantileRCutoff=0.1, permBinom=F, numPerm=5, threshold=0.05)
binTest(snpData, subjInfo, geno, covName, family, quantilePCutoff=0.1, quantileRCutoff=0.1, permBinom=T, numPerm=5, threshold=0.05)