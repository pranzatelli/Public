# It's important to load in required libraries. These libraries
# are the most widely used libraries for interpreting microarray
# data in R.
suppressMessages(library('limma'))
suppressMessages(library('sva'))
suppressMessages(library('Biobase'))

# The first thing to do is load the files, skipping the first
# nine lines (all metadata) and populating a matrix of probe
# intensities per sample.
files1 <- c('JPNSS/SS1.txt','JPNSS/SS4.txt','JPNSS/SS5.txt',
	'JPNSS/SS7.txt','JPNSS/SS8.txt','JPNSS/SS9.txt',
	'JPNSS/SS11.txt','JPNSS/SS12.txt')
skip <- 9

ExprsMatrixJPN <- matrix(ncol=45220,nrow=8)
idx <- vector("list",45220)
i <- 1
for(file in files1) {
	DataFrame <- read.delim(file,skip=skip)
	ExprsMatrixJPN[i,] <- DataFrame$gProcessedSignal
	idx <- DataFrame$ProbeName
	i <- i + 1
}

# And, of course, we have to make sure that the columns are
# appropriately labeled by probe and the rows labeled by
# sample.
rownames(ExprsMatrixJPN) <- files
colnames(ExprsMatrixJPN) <- idx


# Almost identical operation as above.
files2 <- c('USSS/37.txt','USSS/94.txt','USSS/101.txt',
	'USSS/105.txt','USSS/106.txt')
skip <- 9

ExprsMatrixUS <- matrix(ncol=45225,nrow=5)
idx2 <- vector("list",45225)
i <- 1
for(file in files2) {
	DataFrame <- read.delim(file,skip=skip)
	ExprsMatrixUS[i,] <- DataFrame$gProcessedSignal
	idx2 <- DataFrame$ProbeName
	i <- i + 1
}

rownames(ExprsMatrixUS) <- files
colnames(ExprsMatrixUS) <- idx2

# We turn the matrices to dataframes, which makes it easy
# to merge into a single frame of all samples. The five
# extra files are deleted using na.omit, which removes lines
# containing NA values (because they do not appear in all 
# samples).
dataframe1 <- data.frame(t(ExprsMatrixUS))
dataframe2 <- data.frame(t(ExprsMatrixJPN))

ExprsMatrix <- as.matrix(merge(dataframe2,dataframe1,by='row.names')[-1])
ExprsMatrix <- na.omit(ExprsMatrix)

# Here we need to estimate the gene values, and remove genes
# with too low expression. The standard cutoff is 20th percentile.
EM.rows.means <- rowMeans(ExprsMatrix)
twenty_cutoff <- quantile(EM.rows.means,c(.2))
ExprsMatrix <- ExprsMatrix[apply(as.array(EM.rows.means), 
	1, function(x) x>twenty_cutoff),]

# We read in a short file describing some facts about our samples.
pheno <- read.delim('metadata/pheno.delim',header=TRUE)
rownames(pheno) <- pheno$ID
pheno <- pheno[-1]

# Here we perform the surrogate variable analysis. The surrogate
# variables are eigenvectors that stand in for the complicated
# processes of environmental impact on any sample - like reducing
# a high-dimensional batch effect to its principal components.
mod <- model.matrix(~as.factor(Country),data=pheno)
mod0 <- model.matrix(~1,data=pheno)
svobj <- sva(ExprsMatrix,mod,mod0)


# Extra stuff, checking FDR-corrected q-values to make sure there
# are still usable genes and we haven't wiped out all information.
pValues <- f.pvalue(ExprsMatrix,mod,mod0)
qValues <- p.adjust(pValues,method="BH")

modSv <- cbind(mod,svobj$sv)
mod0Sv <- cbind(mod0,svobj$sv)
pValuesSv <- f.pvalue(ExprsMatrix,modSv,mod0Sv)
qValuesSv <- p.adjust(pValuesSv,method="BH")





