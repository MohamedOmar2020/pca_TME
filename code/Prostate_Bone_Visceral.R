###################################################################################
### Mohamed Omar
### Jul 12 2022
### ### Goal: test whether the PRN stromal signature can classify prostate cancer metastasis into Bone // Visceral metastasis

#################################################################################

###### 
# Clean Work space
rm(list = ls())

############################################################################
### Load library
require(Biobase)
require(limma)
require(pROC)
require(caret)
require(RColorBrewer)
require(ggplot2)
require(reshape)
require(plotROC)
library(GEOquery)

Sys.setenv(VROOM_CONNECTION_SIZE = 5000072)

############################################################################
## Download the data sets

#Dataset1 <- getGEO("GSE101607", GSEMatrix = TRUE, AnnotGPL = TRUE)
#Dataset1 <- Dataset1$GSE101607_series_matrix.txt.gz

#Dataset2 <- getGEO("GSE74685", GSEMatrix = TRUE, AnnotGPL = TRUE)
#Dataset2 <- Dataset2$GSE74685_series_matrix.txt.gz

save(Dataset1, Dataset2, file = './objs/bone_visceral_datasets.rda')

load('./objs/bone_visceral_datasets.rda')

############################################################

## Get the expression matrices

expr1 <- exprs(Dataset1)
expr2 <- exprs(Dataset2)

## Get feature data
FeatData1 <- fData(Dataset1)
FeatData2 <- fData(Dataset2)

## Get phenotype
pheno1 <- pData(Dataset1)
pheno2 <- pData(Dataset2)

#############################################################

## Annotate exprs1
rownames(expr1) <- FeatData1$`Gene symbol`
summary(is.na(rownames(expr1)))
rownames(expr1) <- gsub("-", "", rownames(expr1))
rownames(expr1) <- gsub("_", "", rownames(expr1))
sel <- which(apply(expr1, 1, function(x) all(is.finite(x)) ))
expr1 <- expr1[sel, ]
expr1 <- expr1[!(rownames(expr1) == ""), ]
# log transform
expr1 <- log2(expr1)


## Annotate expr2
rownames(expr2) <- FeatData2$GENE_SYMBOL
summary(is.na(rownames(expr2)))
rownames(expr1) <- gsub("-", "", rownames(expr1))
rownames(expr1) <- gsub("_", "", rownames(expr1))
sel <- which(apply(expr2, 1, function(x) all(is.finite(x)) ))
expr2 <- expr2[sel, ]
expr2 <- expr2[!(rownames(expr2) == ""), ]
dim(expr2)

#########################################

## Modify pheno1 

# Keep only prostate samples
pheno1 <- pheno1[grep("prostate", pheno1$`disease state:ch1`), ]
# Remove untreated prostate cancer (keep only castration-resistant)
pheno1 <- pheno1[grep("Castration", pheno1$`disease state:ch1`), ]

pheno1$Mets <- as.factor(pheno1$`tissue:ch1`)
levels(pheno1$Mets) <- c("Bone_Mets")

## Modify expr1
expr1 <- expr1[,colnames(expr1) %in% rownames(pheno1)]

############################################################
## Modify pheno2

# # Remove xenograft samples
# pheno2 <- pheno2[-c(1:20), ]
# 
# # Remove Primary tumors
# pheno2 <- pheno2[-c(46:56), ]
# 
# pheno2$Mets <- pheno2$title
# 
# pheno2$Mets <- gsub(".+_", "", pheno2$Mets)
# 
# pheno2$Mets[pheno2$Mets == "LIVER CRPC metastasis"] <- "Visceral_Mets"
# pheno2$Mets[pheno2$Mets == "LN CRPC metastasis"] <- "Visceral_Mets"
# pheno2$Mets[pheno2$Mets == "RETROPERITONEAL CRPC metastasis"] <- "Visceral_Mets"
# pheno2$Mets[pheno2$Mets == "LUNG CRPC metastasis"] <- "Visceral_Mets"
# pheno2$Mets[pheno2$Mets == "KIDNEY CRPC metastasis"] <- "Visceral_Mets"
# pheno2$Mets[pheno2$Mets == "APPENDIX CRPC metastasis"] <- "Visceral_Mets"
# pheno2$Mets[pheno2$Mets == "PERITONEUM CRPC metastasis"] <- "Visceral_Mets"
# 
# pheno2$Mets[pheno2$Mets == "BONE CRPC metastasis"] <- "Bone_Mets"
# 
# pheno2$Mets <- as.factor(pheno2$Mets)
# levels(pheno2$Mets)
# table(pheno2$Mets)
# 
# ## Modify expr2
# expr2 <- expr2[, colnames(expr2) %in% rownames(pheno2)]
# 

###################################################################
## Modify pheno3
pheno2$Mets <- pheno2$title
pheno2$Mets <- gsub(".+_", "", pheno2$Mets)

pheno2$Mets[pheno2$Mets == "LIVER CRPC metastasis"] <- "Visceral_Mets"
pheno2$Mets[pheno2$Mets == "LN CRPC metastasis"] <- "Visceral_Mets"
pheno2$Mets[pheno2$Mets == "LUNG CRPC metastasis"] <- "Visceral_Mets"
pheno2$Mets[pheno2$Mets == "RETROPERITONEAL CRPC metastasis"] <- "Visceral_Mets"
pheno2$Mets[pheno2$Mets == "RENAL CRPC metastasis"] <- "Visceral_Mets"
pheno2$Mets[pheno2$Mets == "ADRENAL CRPC metastasis"] <- "Visceral_Mets"
pheno2$Mets[pheno2$Mets == "KIDNEY CRPC metastasis"] <- "Visceral_Mets"
pheno2$Mets[pheno2$Mets == "APPENDIX CRPC metastasis"] <- "Visceral_Mets"
pheno2$Mets[pheno2$Mets == "PERITONEUM CRPC metastasis"] <- "Visceral_Mets"
pheno2$Mets[pheno2$Mets == "PERITONEAL CRPC metastasis"] <- "Visceral_Mets"
pheno2$Mets[pheno2$Mets == "SCROTUM CRPC metastasis"] <- "Visceral_Mets"
pheno2$Mets[pheno2$Mets == "SKIN CRPC metastasis"] <- "Visceral_Mets"
pheno2$Mets[pheno2$Mets == "SPLEEN CRPC metastasis"] <- "Visceral_Mets"

pheno2$Mets[pheno2$Mets == "BONE CRPC metastasis"] <- "Bone_Mets"

pheno2$Mets <- as.factor(pheno2$Mets)

#############################################################
#############################################################

## Combine the expression and phenotype

AllExprs <- list(expr1, expr2)
names(AllExprs) <- c("GSE101607", "GSE74685")

AllPheno <- list(pheno1, pheno2)
names(AllPheno) <- c("GSE101607", "GSE74685")

GroupMets <- c(pheno1$BoneMets, pheno2$Mets)
GroupMets <- as.factor(GroupMets)
levels(GroupMets) <- c("Bone_Mets", "Visceral_Mets")
table(GroupMets)

### Find commom subset of genes
commonGenes <- Reduce("intersect", lapply(AllExprs, rownames))

### Filter expression for the common genes
exprsMetastasis <- mapply(x=AllExprs, FUN=function(x, gns) {
  x <- x[ gns ,]
}, MoreArgs=list(gns=commonGenes))

##########
## Assemble in one data frame
AllMat <- do.call("cbind", exprsMetastasis)
## Check if sample names are identical
all(colnames(AllMat) == names(GroupMets))

## Divide data into training and testing data
# ind <- createDataPartition(y = GroupMets, p = 0.7, list = FALSE)
# 
# TrainGroup <- GroupMets[ind]
# TestGroup <- GroupMets[-ind]
# 
# trainMat <- AllMat[, ind]
# testMat <- AllMat[, -ind]

###################################################################
## Load the PRN stromal signature
load("./objs/PRN_stromal_signature.rda")

### filter out genes not present in the matrix
keepGns <- intersect(as.vector(ktspPredictorRes$TSPs), rownames(AllMat))

keep <- ktspPredictorRes$TSPs[,1] %in% rownames(AllMat) & ktspPredictorRes$TSPs[,2] %in% rownames(AllMat)
table(keep)

## Subset
ktspPredictor_metOrigin <- ktspPredictorRes
ktspPredictor_metOrigin$score <- ktspPredictor_metOrigin$score[keep]
ktspPredictor_metOrigin$TSPs <- ktspPredictor_metOrigin$TSPs[keep, ]
ktspPredictor_metOrigin$tieVote <- droplevels(ktspPredictor_metOrigin$tieVote[keep])
ktspPredictor_metOrigin$name <- paste0(nrow(ktspPredictor_metOrigin$TSPs), 'TSPs')


## Normalization between arrays
usedMat <- normalizeBetweenArrays(AllMat, method = "quantile")
boxplot(usedMat)

usedGroup <- GroupMets

############################################################################
## Predict the met origin
ktspPredictor_metOrigin$labels <- c('Visceral_Mets', 'Bone_Mets')
usedGroup <- factor(usedGroup, levels = c('Visceral_Mets', 'Bone_Mets'))

### Compute the sum and find the best threshold
ktspStats_metOrigin <- SWAP.KTSP.Statistics(inputMat = usedMat, classifier = ktspPredictor_metOrigin, CombineFunc = sum)
summary(ktspStats_metOrigin$statistics)

### Threshold
thr <- coords(roc(usedGroup, ktspStats_metOrigin$statistics, levels = c("Visceral_Mets", "Bone_Mets"), direction = "<",), transpose = TRUE, "best")["threshold"]
thr

### Print ROC curve local maximas
coords(roc(usedTrainGroup, ktspStats_metOrigin$statistics, levels = c("Bone_Mets", "Visceral_Mets"), direction = "<", ), transpose = TRUE ,"local maximas")

### Plot Curve: note that you must reorder the levels!!!
### ("good" goes first, "bad" goes second, the opposite of confusionMatrix)
roc(usedTrainGroup, ktspStats_metOrigin$statistics, plot = TRUE, print.thres=thr, print.thres.adj=c(0.01,1.25), print.auc=TRUE, print.auc.col="black", levels = c("Bone_Mets", "Visceral_Mets"), direction = "<", col="blue", lwd=2, grid=TRUE, main="Mechanistic KTSP performance in the training data")

### Get predictions based on best threshold from ROC curve
usedTrainPredictionRes <- SWAP.KTSP.Classify(usedTrainMat, ktspPredictorRes, DecisionFunc = function(x) sum(x) > thr)

### Resubstitution performance in the TRAINING set
confusionMatrix(usedTrainPredictionRes, usedTrainGroup, positive = "Bone_Mets")

#########################################################################
#########################################################################