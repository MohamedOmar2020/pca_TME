###################################################################################
### Mohamed Omar
### 5/5/2019
### ### Goal: Creating the restricted ktsp classifier.
### Combined Adhesion and activation genes
#################################################################################

###### 
# Clean Work space
rm(list = ls())

# Set work directory
#setwd("/Volumes/Macintosh/Dropbox (MechPred)/MechPred/User/Mohamed/MechanisticModels/Prostate")

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

###########################################################################
### Load expression and phenotype data
load("./data/bulk/MetastasisDataGood.rda")

############################################################################
# load the pairs

# load mechanistic pairs
load('./data/MYC_pairs.rda')
load('./data/MYC_pairs_top.rda')


myTSPs <- as.matrix(MYC_pairs_top)

### Quantile normalize
usedTrainMat <- normalizeBetweenArrays(mixTrainMat, method = "quantile")
usedTestMat <- normalizeBetweenArrays(mixTestMat, method = "quantile")

### Common genes
keepGns <- intersect(as.vector(myTSPs), rownames(usedTrainMat))
#keepGns_TF_MiR <- keepGns
#save(keepGns_TF_MiR, file = "./Objs/KTSP/KeepGns_TF_MiR.rda")

#usedTrainMat <- usedTrainMat[keepGns, ]
#usedTestMat <- usedTestMat[keepGns, ]

### Associated groups
usedTrainGroup <- mixTrainGroup
usedTestGroup <- mixTestGroup

### For the TSP
myTSPs <- myTSPs[myTSPs[,1] %in% keepGns & myTSPs[,2] %in% keepGns , ]

#print(xtable(myTSPs, type = "latex"), file = "./Objs/KTSP/Restricted_Pairs.tex")
#write.csv(myTSPs, file = "./Objs/KTSP/Restricted_Pairs.csv")

###########################################################################
### TRAINING using restricted pairs
###########################################################################

### Set Feature number and max k
ktsp <- c(3:25)
featNo <- nrow(usedTrainMat)

### Train a classifier using default filtering function based on Wilcoxon
set.seed(333)

ktspPredictorRes <- SWAP.Train.KTSP(
  usedTrainMat, usedTrainGroup, krange = 25,
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

### Threshold
thr <- coords(roc(usedTrainGroup, ktspStatsTrainRes$statistics, levels = c("No_Mets", "Mets"), direction = "<"), "best")["threshold"]
thr 

### Print ROC curve local maximas
coords(roc(usedTrainGroup, ktspStatsTrainRes$statistics, levels = c("No_Mets", "Mets"), direction = "<"), "local maximas")

### Plot Curve: note that you must reorder the levels!!!
### ("good" goes first, "bad" goes second, the opposite of confusionMatrix)
roc(usedTrainGroup, ktspStatsTrainRes$statistics, plot = F, print.thres=thr, print.auc=TRUE, print.auc.col="black", levels = c("No_Mets", "Mets"), direction = "<", col="blue", lwd=2, grid=TRUE, main="Mechanistic KTSP performance in the training data")

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

## Plot curve
roc(usedTestGroup, ktspStatsTestRes$statistics, plot = F, print.auc=TRUE, print.auc.col="black", levels = c("No_Mets", "Mets"), direction = "<", col="blue", lwd=2, grid=TRUE, main= "Mechanistic KTSP using TF_MiR Gns")

### Get predictions based on best threshold from ROC curve
usedTestPredictionRes <- SWAP.KTSP.Classify(usedTestMat, ktspPredictorRes, DecisionFunc = function(x) sum(x) > thr)

### Resubstitution performance in the Test set
confusionMatrix(usedTestPredictionRes, usedTestGroup, positive = "Mets", mode = "everything")

MCC_Mechanistic_Test <- mltools::mcc(pred = usedTestPredictionRes, actuals = usedTestGroup)
MCC_Mechanistic_Test

############################################################################
###########################################################################
### Plot genes in the training set
## Which TSPs
i <- 1:nrow(ktspPredictorRes$TSPs)

## Assemble in a data frame
tsp <- lapply(i , function(i, x) as.character(unlist(x[i, 1:2])), x=ktspPredictorRes$TSPs)

## Assemble
dfTspTrain <- lapply(tsp, function(i,x,g){
  out <- data.frame(t(x[i, ]), Group=g)
  out <- cbind(pair= paste("TSP:", paste(i, collapse = "-")), out)
}, x=usedTrainMat, g=usedTrainGroup)

names(dfTspTrain) <- rownames(ktspPredictorRes$TSPs)

# Change the names of elements inside each element in dfTspTrain (For Plotting)  
for(i in seq_along(dfTspTrain)) names(dfTspTrain[[i]]) <-  c("pair", "Gene1", "Gene2", "Group")



## Reduce
datTrain <- Reduce("rbind", dfTspTrain)

## Rename columns
#colnames(datTrain)[colnames(datTrain) %in% c("variable", "value")] <- c("Gene", "Expression")

##############
## Scatter plots
png(filename = "./figures/bulk/MechanisticKTSP_Train_ScatterPlot_Combined.png", width=3000, height=1500, res=300)
### Prepare ggplot
sctplt <- ggplot(na.omit(datTrain), aes(x=Gene1, y=Gene2, color=Group)) +
  geom_point(size=0.5) + geom_abline(intercept=0, slope=1, alpha=0.5) +
  facet_wrap(~ pair, scales="free", nrow=2) +
  theme(axis.text=element_text(face="bold", size = 6),
        axis.title=element_text(face="bold", size = 8),
        legend.position="bottom",
        legend.text=element_text(face="bold", size=10),
        legend.title = element_text(face="bold", size=10),
        strip.text.x = element_text(face="bold", size=6))
sctplt
dev.off()


########################################################################
#######################################################################
### Plot genes in the testing set

## Which TSPs
i <- 1:nrow(ktspPredictorRes$TSPs)

## Assemble in a data frame
tsp <- lapply(i, function(i, x) as.character(unlist(x[i,1:2])), x=ktspPredictorRes$TSPs)

## Assemble
dfTspTest <- lapply(tsp, function(i, x, g){
  out <- data.frame(t(x[i, ]), Group=g)
  out <- cbind(pair=paste("TSP:", paste(i, collapse = "-")), out)
}, x=usedTestMat, g=usedTestGroup)


names(dfTspTest) <- rownames(ktspPredictorRes$TSPs)

# Change the names of elements inside each element in dfTspTrain 
for(i in seq_along(dfTspTest)) names(dfTspTest[[i]]) <-  c("pair", "Gene1", "Gene2", "Group")


## Reduce
datTest <- Reduce("rbind", dfTspTest)

## Rename columns
#colnames(datTest)[colnames(datTest) %in% c("variable", "value")] <- c("Gene", "Expression")

########################
## Make paired boxplot
# png("./Figs/KTSP/mechanistic.testKTSPexprs.png", width = 3000, height = 1500, res = 200)
# bxplt <- ggplot(na.omit(datTest), aes(x=Gene, y=Expression, fill=Group)) +
#   geom_boxplot(outlier.shape = NA) + 
#   geom_point(aes(group=Group), position = position_jitterdodge(), size=0.75, color=rgb(0.2,0.2,0.2,0.5)) + 
#   facet_wrap(~pair, scales = "free", nrow = 2) +
#   theme(axis.text = element_text(face = "bold", size = 12), axis.title = element_text(face = "bold", size = 12), legend.position = "bottom", legend.text = element_text(face = "bold", size = 17.5), legend.title = element_text(face = "bold", size= 17.5), strip.text.x = element_text(face = "bold", size = 11))
# bxplt
# dev.off()

#####################
## Scatter plots
png(filename = "./figures/bulk/MechanisticKTSP_Test_ScatterPlot_Combined.png", width=3000, height=1500, res=300)
### Prepare ggplot
sctplt <- ggplot(na.omit(datTest), aes(x=Gene1, y=Gene2, color=Group)) +
  geom_point(size=0.5) + geom_abline(intercept=0, slope=1, alpha=0.5) +
  facet_wrap(~ pair, scales="free", nrow=2) +
  theme(axis.text=element_text(face="bold", size = 6),
        axis.title=element_text(face="bold", size = 8),
        legend.position="bottom",
        legend.text=element_text(face="bold", size=10),
        legend.title = element_text(face="bold", size=10),
        strip.text.x = element_text(face="bold", size=6))
sctplt
dev.off()

#########################################################################
######################################################

## GGPLOT to compare the two classifiers
### Prepare the legend
forLegend_KTSP <- apply(rbind(
  ci(roc(usedTrainGroup, ktspStatsTrainRes$statistics, levels = c("No_Mets", "Mets"), direction = "<")),
  ci(roc(usedTestGroup, ktspStatsTestRes$statistics, levels = c("No_Mets", "Mets"), direction = "<"))
  #ci(roc(usedTrainGroup, ktspStatsTrainUnRes$statistics, levels = c("No_Mets", "Mets"), direction = "<")),
  #ci(roc(usedTestGroup, ktspStatsTestUnRes$statistics, levels = c("No_Mets", "Mets"), direction = "<"))
),  1, function(x) {
  x <- format(round(x, digits=2), nsmall=2)
  paste("AUC: ", x[[2]], ";", "95% CI: ", x[[1]], "-", x[[3]])
})


#################################################################
### ROC curves Using ggplot2

### Training
datTrn_KTSP <- melt(data.frame(
  ## Training Group
  Training= usedTrainGroup,
  ## Agnostic KTSP SUM: the lowest mus be for disease status
  #Agnostic.Training = ktspStatsTrainUnRes$statistics,
  ## Mechanistic KTSP SUM training
  PRN_signature.Training= ktspStatsTrainRes$statistics))
### Change Colnames
colnames(datTrn_KTSP) <- c("Status", "KTSP_type", "KTSP_sum")


### Testing
datTst_KTSP <- melt(data.frame(
  ## Testing group
  Testing= usedTestGroup,
  ## Agnostic KTSP SUM: the lowest mus be for disease status
  #Agnostic.Testing=ktspStatsTestUnRes$statistics,
  ## Mechanistic KTSP SUM training
  PRN_signature.Testing=ktspStatsTestRes$statistics))
### Change Colnames
colnames(datTst_KTSP) <- c("Status", "KTSP_type", "KTSP_sum")

### Combine
dat_KTSP <- rbind(datTrn_KTSP, datTst_KTSP)
dat_KTSP$Status <- as.numeric(dat_KTSP$Status)-1
####
### Replace levels
levels(dat_KTSP$KTSP_type) <- gsub("\\.", "-", levels(dat_KTSP$KTSP_type))
levels(dat_KTSP$KTSP_type) <- paste(levels(dat_KTSP$KTSP_type), forLegend_KTSP[c(1,2)])

#################################################################
### Plot Curve
png("./figures/bulk/CompareAUCggplot_Combined.png",
    width=3000, height=3000, res=360)
### Color
myCol <- brewer.pal(3, "Dark2")[c(2,1)]
### Plot and legend titles
plotTitle <- "Performance of the PRN stromal signature at predicting metastasis"
#legendTitle <- paste("Mechanistic (", nrow(ktspPredictorRes$TSPs), " pairs)",
#                     " Agnostic (", nrow(ktspPredictorUnRes$TSPs), " pairs)",  sep="")
### Plot
basicplot_KTSP_Combined <- ggplot(dat_KTSP, aes(d=Status, m=KTSP_sum, color=KTSP_type,
                                                linetype = KTSP_type)) +
  geom_roc(cutoffs.at = seq(1,20,2)) +
  style_roc(theme = theme_grey) + ggtitle(plotTitle) +
  theme(plot.title = element_text(face="bold", size=16, hjust = 0.5),
        axis.text=element_text(face="plain", size = 11),
        axis.title=element_text(face="bold", size = 13),
        legend.justification=c(1,0),  legend.position=c(1,0),
        legend.background=element_rect(fill="lightblue1"),
        legend.text=element_text(face="plain", size = 10),
        legend.title = element_text(face="bold", size=12)) +
  #scale_color_manual(legendTitle, values=rep(myCol, 2)) +
  scale_linetype_manual(values=rep(c("solid", "dotted"), each=2)) +
  guides(colour = guide_legend(override.aes = list(size=3)))
### Plot
basicplot_KTSP_Combined
### Close device
dev.off()
