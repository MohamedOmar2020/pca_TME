###################################################################################
### Mohamed Omar
### ### Goal: Creating the restricted ktsp classifier.
#################################################################################

###### 
# Clean Work space
rm(list = ls())

# Set work directory
setwd('/Users/mohamedomar/Documents/Research/Projects/Pca_TME')
############################################################################
### Load library
require(switchBox)
require(Biobase)
require(limma)
require(pROC)
require(caret)
require(RColorBrewer)
require(ggplot2)
require(reshape)
require(plotROC)
library(enrichR)
library(mltools)
library(xtable)
library(precrec)
library(tidyverse)
library(data.table)

# Increase Java Heap Space in R:
options(java.parameters = "-Xmx8000m")
library(xlsx)

###########################################################################
### Load expression and phenotype data
load("./data/bulk/MetastasisDataGood.rda")

############################################################################
# load the pairs

# load mechanistic pairs
load('./data/PRN_pairs_top.rda')


myTSPs <- as.matrix(PRN_pairs_top)


### Quantile normalize
usedTrainMat <- normalizeBetweenArrays(mixTrainMat, method = "quantile")
usedTestMat <- normalizeBetweenArrays(mixTestMat, method = "quantile")

### Common genes
keepGns <- intersect(as.vector(myTSPs), rownames(usedTrainMat))

### Associated groups
usedTrainGroup <- mixTrainGroup
usedTestGroup <- mixTestGroup

### For the TSP
myTSPs <- myTSPs[myTSPs[,1] %in% keepGns & myTSPs[,2] %in% keepGns , ]

########################
# save for source data
write.xlsx(myTSPs, file = 'tables/Source_data.xlsx', append = TRUE, row.names = F, sheetName = '6a - PRN gene pairs')

###########################################################################
### TRAINING using restricted pairs
###########################################################################

### Set Feature number and max k
ktsp <- c(3:25)
featNo <- nrow(usedTrainMat)

### Train a classifier using default filtering function based on Wilcoxon
set.seed(333)

ktspPredictorRes <- SWAP.Train.KTSP(
  usedTrainMat, usedTrainGroup, krange = ktsp,
  FilterFunc = SWAP.Filter.Wilcoxon, featureNo=featNo, RestrictedPairs = myTSPs, disjoint = F)

ktspPredictorRes

### Check consistency with biology
keepTest <- ktspPredictorRes$TSPs[,1] %in% myTSPs[,"up"] & ktspPredictorRes$TSPs[,2] %in% myTSPs[,"down"]

table(keepTest)

## Subset
ktspPredictorRes$score <- ktspPredictorRes$score[keepTest]
ktspPredictorRes$TSPs <- ktspPredictorRes$TSPs[keepTest, ]
ktspPredictorRes$tieVote <- droplevels(ktspPredictorRes$tieVote[keepTest])
ktspPredictorRes$name <- paste0(nrow(ktspPredictorRes$TSPs), 'TSPs')

save(ktspPredictorRes, file = './objs/PRN_stromal_signature.rda')

############################################################################
### Compute the sum and find the best threshold: All training samples
ktspStatsTrainRes <- SWAP.KTSP.Statistics(inputMat = usedTrainMat, classifier = ktspPredictorRes, CombineFunc = sum)
summary(ktspStatsTrainRes$statistics)

######
# save for source data
train_stats <- as.data.frame(ktspStatsTrainRes$statistics)
colnames(train_stats) <- 'votes'
write.xlsx(train_stats, file = 'tables/Source_data.xlsx', append = TRUE, row.names = T, sheetName = '6a - training votes')

### Threshold
thr <- coords(roc(usedTrainGroup, ktspStatsTrainRes$statistics, levels = c("No_Mets", "Mets"), direction = "<"), "best")["threshold"]
thr 

### Print ROC curve local maximas
coords(roc(usedTrainGroup, ktspStatsTrainRes$statistics, levels = c("No_Mets", "Mets"), direction = "<"), "local maximas")

### Plot Curve: note that you must reorder the levels!!!
### ("good" goes first, "bad" goes second, the opposite of confusionMatrix)
roc_train <- roc(usedTrainGroup, ktspStatsTrainRes$statistics, plot = F, print.thres=thr, print.auc=TRUE, print.auc.col="black", levels = c("No_Mets", "Mets"), direction = "<", col="blue", lwd=2, grid=TRUE, main="Mechanistic KTSP performance in the training data")

### Get predictions based on best threshold from ROC curve
usedTrainPredictionRes <- SWAP.KTSP.Classify(usedTrainMat, ktspPredictorRes, DecisionFunc = function(x) sum(x) > thr)

### Resubstitution performance in the TRAINING set
confusionMatrix(usedTrainPredictionRes, usedTrainGroup, positive = "Mets", mode = "everything")

MCC_Mechanistic_Train <- mltools::mcc(pred = usedTrainPredictionRes, actuals = usedTrainGroup)
MCC_Mechanistic_Train

#########################################################################
#########################################################################
### Testing

## Compute the sum and find the best threshold
ktspStatsTestRes <- SWAP.KTSP.Statistics(inputMat = usedTestMat, classifier = ktspPredictorRes, CombineFunc = sum)
summary(ktspStatsTestRes$statistics)

######
# save for source data
test_stats <- as.data.frame(ktspStatsTestRes$statistics)
colnames(test_stats) <- 'votes'
write.xlsx(test_stats, file = 'tables/Source_data.xlsx', append = TRUE, row.names = T, sheetName = '6a - testing votes')

## Plot curve
roc_test <- roc(usedTestGroup, ktspStatsTestRes$statistics, plot = F, print.auc=TRUE, print.auc.col="black", levels = c("No_Mets", "Mets"), direction = "<", col="blue", lwd=2, grid=TRUE, main= "Mechanistic KTSP using TF_MiR Gns")

### Get predictions based on best threshold from ROC curve
usedTestPredictionRes <- SWAP.KTSP.Classify(usedTestMat, ktspPredictorRes, DecisionFunc = function(x) sum(x) > thr)

### Resubstitution performance in the Test set
confusionMatrix(usedTestPredictionRes, usedTestGroup, positive = "Mets", mode = "everything")

MCC_Mechanistic_Test <- mltools::mcc(pred = usedTestPredictionRes, actuals = usedTestGroup)
MCC_Mechanistic_Test

############################################################################
## ROC curve for paper

scores <- join_scores(ktspStatsTrainRes$statistics, ktspStatsTestRes$statistics, chklen = F)

labels <- join_labels(usedTrainGroup, usedTestGroup, chklen = F)

ROC_data <- mmdata(scores = scores, labels = labels, dsids = c(1,2), modnames = c('Training','Testing'), expd_first = "dsids")

sscurves <- evalmod(ROC_data, modnames = c('1','2'), raw_curves = T, expd_first = "mids")
sscurves


pdf("./figures/PRN_signature.pdf", width=8, height=8)
autoplot(sscurves, curvetype = c("ROC"), raw_curves = TRUE, show_legend = TRUE, show_cb = FALSE, ret_grob = TRUE) + 
  labs(title = "") + 
  annotate("text", x = .65, y = .25, label = "Training AUC = 0.69", size = 3) + 
  annotate("text", x = .65, y = .18, label = "Testing AUC = 0.70", size = 3) +
  theme(axis.title.y = element_text(family = 'Arial', face = 'bold.italic', size = 7, angle = 90), 
        axis.title.x = element_text(family = 'Arial', face = 'bold.italic', size = 7, angle = 0),
        axis.text.x = element_text(family = 'Arial', size = 6), 
        axis.text.y = element_text(family = 'Arial', size = 6), 
        legend.text = element_text(family = 'Arial', face = 'bold', size = 7), 
        )
dev.off()
