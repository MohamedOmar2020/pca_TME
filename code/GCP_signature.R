###################################################################################
### Mohamed Omar
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
library(glmnet)

###########################################################################
### Load expression and phenotype data
load("./data/bulk/MetastasisDataGood.rda")

### Quantile normalize
usedTrainMat <- normalizeBetweenArrays(mixTrainMat, method = "quantile")
usedTestMat <- normalizeBetweenArrays(mixTestMat, method = "quantile")

### Associated groups
usedTrainGroup <- mixTrainGroup
usedTestGroup <- mixTestGroup

# Prepare the data
Data_train <- as.data.frame(cbind(t(usedTrainMat), usedTrainGroup))
Data_train$usedTrainGroup <- as.factor(Data_train$usedTrainGroup)
levels(Data_train$usedTrainGroup) <- c('No_Mets', 'Mets')
colnames(Data_train)[colnames(Data_train) == 'usedTrainGroup'] <- 'group'
table(Data_train$group)

Data_test <- as.data.frame(cbind(t(usedTestMat), usedTestGroup))
Data_test$usedTestGroup <- as.factor(Data_test$usedTestGroup)
levels(Data_test$usedTestGroup) <- c('No_Mets', 'Mets')
colnames(Data_test)[colnames(Data_test) == 'usedTestGroup'] <- 'group'
table(Data_test$group)


############################################################################
# load PRN signature
load('./objs/PRN_stromal_signature.rda')

PRN_genes <- as.vector(ktspPredictorRes$TSPs)

# define GCP signature
GCP_genes <- c("FOXM1", "ASPM", "TK1", "PRC1",
                 "CDC20", "BUB1B", "PBK", "DTL",
                 "CDKN3", "RRM2", "ASF1B", "CEP55",
                 "CDC2", "DLGAP5", "C18orf24", "RAD51",
                 "KIF11", "BIRC5", "RAD54L", "CENPM",
                 "KIAA0101", "KIF20A", "PTTG1", "CDCA8",
                 "NUSAP1", "PLK1", "CDCA3", "ORC6L",
                 "CENPF", "TOP2A", "MCM10")

summary(GCP_genes %in% rownames(mixTestMat))
intersect(PRN_genes, GCP_genes)


### Common genes
GCP_genes_fil <- intersect(GCP_genes, rownames(usedTrainMat))

################################################################
# fit the models using logistic regression
################################################################
PRN_model <- glm(as.formula((paste("group ~", paste(PRN_genes, collapse = "+")))), data = Data_train, family = "binomial")
summary(PRN_model)

GCP_model <- glm(as.formula((paste("group ~", paste(GCP_genes_fil, collapse = "+")))), data = Data_train, family = "binomial")
summary(GCP_model)
save(GCP_model, file = './objs/GCP_model.rda')

############################################################################
# Make predictions
############################################################################

# Training data
prob_train_PRN <- PRN_model %>% predict(Data_train , type = "response")
prob_train_GCP <- GCP_model %>% predict(Data_train , type = "response")

### Threshold
thr_PRN <- coords(roc(Data_train$group, prob_train_PRN, levels = c('No_Mets', 'Mets'), direction = "<"), "best")["threshold"]
thr_GCP <- coords(roc(Data_train$group, prob_train_GCP, levels = c('No_Mets', 'Mets'), direction = "<"), "best")["threshold"]


### ROC Curve
ROC_train_PRN <- roc(Data_train$group, prob_train_PRN, plot = F, print.auc=TRUE, print.auc.col="black", ci = T, levels = c('No_Mets', 'Mets'), direction = "<", col="blue", lwd=2, grid=TRUE)
ROC_train_PRN

ROC_train_GCP <- roc(Data_train$group, prob_train_GCP, plot = F, print.auc=TRUE, print.auc.col="black", ci = T, levels = c('No_Mets', 'Mets'), direction = "<", col="blue", lwd=2, grid=TRUE)
ROC_train_GCP

### Get predictions based on best threshold from ROC curve
predClasses_train_PRN <- ifelse(prob_train_PRN >= thr_PRN$threshold, 'No_Mets', 'Mets')
table(predClasses_PRN)
predClasses_train_PRN <- factor(predClasses_PRN, levels = c('No_Mets', 'Mets'))

predClasses_train_GCP <- ifelse(prob_train_GCP >= thr_GCP$threshold, 'No_Mets', 'Mets')
table(predClasses_train_GCP)
predClasses_train_GCP <- factor(predClasses_train_GCP, levels = c('No_Mets', 'Mets'))

##################################
# Testing data
##################################
prob_test_PRN <- PRN_model %>% predict(Data_test , type = "response")
prob_test_GCP <- GCP_model %>% predict(Data_test , type = "response")

### ROC Curve
ROC_test_PRN <- roc(Data_test$group, prob_test_PRN, plot = F, print.auc=TRUE, print.auc.col="black", ci = T, levels = c('No_Mets', 'Mets'), direction = "<", col="blue", lwd=2, grid=TRUE)
ROC_test_PRN

ROC_test_GCP <- roc(Data_test$group, prob_test_GCP, plot = F, print.auc=TRUE, print.auc.col="black", ci = T, levels = c('No_Mets', 'Mets'), direction = "<", col="blue", lwd=2, grid=TRUE)
ROC_test_GCP

### Get predictions based on best threshold from ROC curve
predClasses_test_PRN <- ifelse(prob_test_PRN >= thr_PRN$threshold, 'No_Mets', 'Mets')
table(predClasses_test_PRN)
predClasses_test_PRN <- factor(predClasses_test_PRN, levels = c('No_Mets', 'Mets'))

predClasses_test_GCP <- ifelse(prob_test_GCP >= thr_GCP$threshold, 'No_Mets', 'Mets')
table(predClasses_test_GCP)
predClasses_test_GCP <- factor(predClasses_test_GCP, levels = c('No_Mets', 'Mets'))

