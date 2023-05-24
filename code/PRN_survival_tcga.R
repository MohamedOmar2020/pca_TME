############################################################################


# Clean the working directory
rm(list = ls())

# Load necessary packages
require(Biobase)
require(limma)
require(pROC)
require(caret)
require(RColorBrewer)
require(ggplot2)
require(reshape)
require(plotROC)
library(mltools)
library(xtable)
library(dplyr)
library(precrec)
library(patchwork)
library(survminer)
library(survival)
library(tidyverse)
library(dvmisc)
library(pheatmap)

#################
PRN_signature <- load('./objs/PRN_stromal_signature.rda')
PRN_signature 

################
# Load the TCGA expression and pheno data
load("./objs/for_survival.rda")

# load the original training and testing data
load("./data/bulk/MetastasisDataGood.rda")
### Quantile normalize
usedTrainMat <- normalizeBetweenArrays(mixTrainMat, method = "quantile")
usedTestMat <- normalizeBetweenArrays(mixTestMat, method = "quantile")

train_groups <- mixTrainGroup
test_groups <- mixTestGroup

###############
### combine in 1 dataset: Training
Data_train <- as.data.frame(cbind(t(usedTrainMat), train_groups))
Data_train$train_groups <- as.factor(Data_train$train_groups)
levels(Data_train$train_groups) <- c(0, 1)
colnames(Data_train)[colnames(Data_train) %in% c('train_groups')] <- c('label')

##########
### combine in 1 dataset: testing
Data_test <- as.data.frame(cbind(t(usedTestMat), test_groups))
Data_test$test_groups <- as.factor(Data_test$test_groups)
levels(Data_test$test_groups) <- c(0, 1)
colnames(Data_test)[colnames(Data_test) %in% c('test_groups')] <- c('label')

##########
### combine in 1 dataset: tcga
Data_tcga <- as.data.frame(cbind(t(Expr_tcga), group_tcga))
Data_tcga$group_tcga <- as.factor(Data_tcga$group_tcga)
levels(Data_tcga$group_tcga) <- c(0, 1)
colnames(Data_tcga)[colnames(Data_tcga) %in% c('group_tcga')] <- c('label')

###########################################################################
### TRAINING using logistic regression
###########################################################################

# get the genes in common
PRN_signature <- as.vector(ktspPredictorRes$TSPs)

PRN_signature_fil <- PRN_signature[PRN_signature %in% rownames(usedTrainMat) & PRN_signature %in% rownames(usedTestMat) & PRN_signature %in% rownames(Expr_tcga)]

#############################################################################################################
##############################################################################################################
# the model

logReg_model <- glm(as.formula((paste("label ~", paste(PRN_signature_fil, collapse = "+")))), data = Data_train, family = "binomial")
summary(logReg_model)

save(logReg_model, file = "./objs/logReg_model.rda")

###########################################################################
############################################################################
### predict in the training dataset
# Make predictions

Train_prob_logReg <- logReg_model %>% predict(Data_train , type = "response")

### Threshold
thr <- coords(roc(Data_train$label, Train_prob_logReg, levels = c(0, 1), direction = "<"), "best")["threshold"]
thr

### ROC Curve
ROCTrain_logReg <- roc(Data_train$label, Train_prob_logReg, plot = F, print.thres=thr$threshold, print.auc=TRUE, print.auc.col="black", ci = T, levels = c("0", "1"), direction = "<", col="blue", lwd=2, grid=TRUE)
ROCTrain_logReg

### Get predictions based on best threshold from ROC curve
Train_predClasses_logReg <- ifelse(Train_prob_logReg >= thr$threshold, "1", "0")
table(Train_predClasses_logReg)
Train_predClasses_logReg <- factor(Train_predClasses_logReg, levels = c('0', '1'))


### Resubstitution performance in the TRAINING set
Confusion_train <- confusionMatrix(Train_predClasses_logReg, Data_train$label, positive = "1", mode = "everything")
Confusion_train

## MCC
MCC_train <- mltools::mcc(pred = Train_predClasses_logReg, actuals = Data_train$label)
MCC_train


#########################################################################
#########################################################################
### Testing

Test_prob_logReg <- logReg_model %>% predict(Data_test , type = "response")

### ROC
ROC_test <- roc(Data_test$label, Test_prob_logReg, plot = F, print.thres=thr$threshold, print.auc=TRUE, print.auc.col="black", ci = T, levels = c(0, 1), direction = "<", col="blue", lwd=2, grid=TRUE)
ROC_test

#################
### Get predictions based on best threshold from ROC curve
Test_predClasses_logReg <- ifelse(Test_prob_logReg >= thr$threshold, "1", "0")
table(Test_predClasses_logReg)
Test_predClasses_logReg <- factor(Test_predClasses_logReg, levels = c('0', '1'))

####################
### CI in testing 1
Confusion_test <- confusionMatrix(Test_predClasses_logReg, Data_test$label, positive = "1", mode = "everything")
Confusion_test

################
## MCC
MCC_test <- mltools::mcc(pred = Test_predClasses_logReg, actuals = Data_test$label)
MCC_test

#########################################################################
#########################################################################
### TCGA

tcga_prob_logReg <- logReg_model %>% predict(Data_tcga , type = "response")

### ROC
ROC_tcga <- roc(Data_tcga$label, tcga_prob_logReg, plot = F, print.thres=thr$threshold, print.auc=TRUE, print.auc.col="black", ci = T, levels = c(0, 1), direction = "<", col="blue", lwd=2, grid=TRUE)
ROC_tcga

#################
### Get predictions based on best threshold from ROC curve
tcga_predClasses_logReg <- ifelse(tcga_prob_logReg >= thr$threshold, "1", "0")
table(tcga_predClasses_logReg)
tcga_predClasses_logReg <- factor(tcga_predClasses_logReg, levels = c('0', '1'))

####################
### CI in tcgaing 1
Confusion_tcga <- confusionMatrix(tcga_predClasses_logReg, Data_tcga$label, positive = "1", mode = "everything")
Confusion_tcga

################
## MCC
MCC_tcga <- mltools::mcc(pred = tcga_predClasses_logReg, actuals = Data_tcga$label)
MCC_tcga

##########################
# load another TCGA metadata file to get the tumor purity score
Pheno_tcga2 <- read.delim2('./data/bulk/TCGA/tcga_meta_data_all.csv', sep = ',')
Pheno_tcga2 <- Pheno_tcga2[!duplicated(Pheno_tcga2$patient_id), ]
rownames(Pheno_tcga2) <- Pheno_tcga2$patient_id
CommonSamples <- intersect(rownames(Pheno_tcga), rownames(Pheno_tcga2))
Pheno_tcga2 <- Pheno_tcga2[CommonSamples, ]
all(rownames(Pheno_tcga2) == rownames(Pheno_tcga))

# add the gleason score to CoxData_tcga
Pheno_tcga$purity <- Pheno_tcga2$purity
Pheno_tcga$purity <- as.numeric(Pheno_tcga$purity)

## Keep only the relevant information (Metastasis Event and Time)
Phenotype_tcga <- cbind(Pheno_tcga[, c("Overall.Survival.Status", "Overall.Survival..Months.", "Progression.Free.Status", "Progress.Free.Survival..Months.", "purity")], 
                        tcga_prob_logReg, tcga_predClasses_logReg)

#Expr_metabric <- Expr_metabric[ClassifierGenes, ]
#Expr_tcga <- Expr_tcga[ClassifierGenes, ]

# create a merged pdata and Z-scores object
CoxData_tcga <- data.frame(Phenotype_tcga)

# divide the PRN probability into quartiles
CoxData_tcga <- CoxData_tcga %>%
  mutate(PRN_quartiles = ntile(tcga_prob_logReg, 4))

# Divide purity into quartiles
CoxData_tcga$purity_q <- quant_groups(CoxData_tcga$purity, 4) 
CoxData_tcga_purity_q1 <- CoxData_tcga[CoxData_tcga$purity_q == '[0.16,0.44]', ]
CoxData_tcga_purity_q2 <- CoxData_tcga[CoxData_tcga$purity_q == '(0.44,0.61]', ]
CoxData_tcga_purity_q3 <- CoxData_tcga[CoxData_tcga$purity_q == '(0.61,0.78]', ]
CoxData_tcga_purity_q4 <- CoxData_tcga[CoxData_tcga$purity_q == '(0.78,1]', ]

########################################################################  
## Fit survival curves

# OS
## TCGA all genes
Fit_sig_tcga_os_logReg <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ tcga_predClasses_logReg, data = CoxData_tcga)


## by PRN quartiles
Fit_sig_tcga_os_logReg_quartiles <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ PRN_quartiles, data = CoxData_tcga)

################
# PFS
## metabric all genes
Fit_sig_tcga_pfs_logReg <- survfit(Surv(Progress.Free.Survival..Months., Progression.Free.Status) ~ tcga_predClasses_logReg, data = CoxData_tcga)

## by PRN quartiles
Fit_sig_tcga_pfs_logReg_quartiles <- survfit(Surv(Progress.Free.Survival..Months., Progression.Free.Status) ~ PRN_quartiles, data = CoxData_tcga)

## in samples with highest stromal content/least purity/Q1 
Fit_sig_tcga_pfs_q1 <- survfit(Surv(Progress.Free.Survival..Months., Progression.Free.Status) ~ tcga_predClasses_logReg, data = CoxData_tcga_purity_q1)

## in samples with lowest stromal content/highest purity/Q4 
Fit_sig_tcga_pfs_q4 <- survfit(Surv(Progress.Free.Survival..Months., Progression.Free.Status) ~ tcga_predClasses_logReg, data = CoxData_tcga_purity_q4)

############################################################################
############################################################################
# plot OS

pdf("./figures/survival/LogReg_tcga_os_all.pdf", width = 8, height = 8, onefile = F)
ggsurvplot(Fit_sig_tcga_os_logReg,
           risk.table = FALSE,
           pval = TRUE,
           pval.size = 10,
           legend.labs = c('prediction: 0', 'prediction: 1'),
           ggtheme = theme_survminer(base_size = 30, font.x = c(30, 'bold.italic', 'black'), font.y = c(30, 'bold.italic', 'black'), font.tickslab = c(30, 'plain', 'black'), font.legend = c(30, 'bold', 'black')),
           palette = 'jco',
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = 'PRN signature and TCGA OS')
dev.off()

########
# by quartiles
pdf("./figures/survival/logReg_tcga_os_quartiles.pdf", width = 8, height = 8, onefile = F)
ggsurvplot(Fit_sig_tcga_os_logReg_quartiles,
           risk.table = FALSE,
           pval = TRUE,
           pval.size = 10,
           legend.labs = c('Q1', 'Q2', 'Q3', 'Q4'),
           ggtheme = theme_survminer(base_size = 30, font.x = c(30, 'bold.italic', 'black'), font.y = c(30, 'bold.italic', 'black'), font.tickslab = c(30, 'plain', 'black'), font.legend = c(30, 'bold', 'black')),
           palette = 'jco',
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = 'PRN signature and TCGA OS: quartiles')
dev.off()

######################################
# plot PFS

tiff("./figures/survival/logReg_tcga_pfs_all.tiff", width = 2000, height = 2000, res = 400)
ggsurvplot(Fit_sig_tcga_pfs_logReg,
           risk.table = FALSE,
           pval = TRUE,
           pval.size = 8,
           legend.labs = c('predicted PFS: 0', 'predicted PFS: 1'),
           legend.title	= '',
           ggtheme = theme_survminer(base_size = 12, font.x = c(12, 'bold.italic', 'black'), font.y = c(12, 'bold.italic', 'black'), font.tickslab = c(12, 'plain', 'black'), font.legend = c(12, 'bold', 'black')),
           palette = 'jco',
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           #title = 'PRN signature and TCGA PFS'
           )
dev.off()

########
# by quartiles
tiff("./figures/survival/logreg_tcga_PFS_quartiles.tiff", width = 2000, height = 2000, res = 400)
ggsurvplot(Fit_sig_tcga_pfs_logReg_quartiles,
           risk.table = FALSE,
           pval = TRUE,
           pval.size = 10,
           legend.labs = c('Q1', 'Q2', 'Q3', 'Q4'),
           legend.title	= '',
           ggtheme = theme_survminer(base_size = 12, font.x = c(12), font.y = c(12, 'bold.italic', 'black'), font.tickslab = c(12, 'plain', 'black'), font.legend = c(12, 'bold', 'black')),
           palette = 'jco',
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           #title = 'PRN signature and TCGA PFS: quartiles'
           )
dev.off()

########
# PFS in least pure samples/q1
tiff("./figures/survival/logreg_tcga_PFS_leastPureQ1.tiff", width = 2000, height = 2000, res = 400)
ggsurvplot(Fit_sig_tcga_pfs_q1,
           risk.table = FALSE,
           pval = TRUE,
           pval.size = 10,
           legend.labs = c('predicted PFS: 0', 'predicted PFS: 1'),
           legend.title	= '',
           ggtheme = theme_survminer(base_size = 12, font.x = c(12), font.y = c(12, 'bold.italic', 'black'), font.tickslab = c(12, 'plain', 'black'), font.legend = c(12, 'bold', 'black')),
           palette = 'jco',
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           title = 'PRN signature and TCGA PFS in samples with lowest purity: 145 samples'
)
dev.off()

########
# PFS in most pure samples/q4
tiff("./figures/survival/logreg_tcga_PFS_MostPureQ4.tiff", width = 2000, height = 2000, res = 400)
ggsurvplot(Fit_sig_tcga_pfs_q4,
           risk.table = FALSE,
           pval = TRUE,
           pval.size = 10,
           legend.labs = c('predicted PFS: 0', 'predicted PFS: 1'),
           legend.title	= '',
           ggtheme = theme_survminer(base_size = 12, font.x = c(12), font.y = c(12, 'bold.italic', 'black'), font.tickslab = c(12, 'plain', 'black'), font.legend = c(12, 'bold', 'black')),
           palette = 'jco',
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           title = 'PRN signature and TCGA PFS in samples with highest purity: 139 samples'
)
dev.off()
##############################################################################
## fit coxph model:

########
## PFS

# by probability

Fit_sig_tcga_PFS_coxph_logReg <- coxph(Surv(Progress.Free.Survival..Months., Progression.Free.Status) ~ tcga_prob_logReg, data = CoxData_tcga)
summary(Fit_sig_tcga_PFS_coxph_logReg)

########
## by quartiles

# make a factor with Q1 (lowest risk) being the reference
CoxData_tcga$quartiles <- factor(CoxData_tcga$quartiles, levels = c('1', '2', '3', '4'))
levels(CoxData_tcga$quartiles) <- paste0('Q', levels(CoxData_tcga$quartiles))

# fit
Fit_sig_tcga_pfs_coxph_logReg_quartiles <- coxph(Surv(Progress.Free.Survival..Months., Progression.Free.Status) ~ quartiles, data = CoxData_tcga)
summary_tcga_pfs_coxph_logReg_quartiles <- summary(Fit_sig_tcga_pfs_coxph_logReg_quartiles)

##########################
## multivariate cox with gleason
# read another version of the clinical data (Firehose legacy) which contains gleason score
pheno2 <- read.delim2('data/bulk/TCGA/prad_tcga_clinical_data.tsv')
pheno2$Patient.ID <- gsub("\\-", "\\.", pheno2$Patient.ID)
pheno2 <- pheno2[!duplicated(pheno2$Patient.ID), ]
rownames(pheno2) <- pheno2$Patient.ID
CommonSamples <- intersect(rownames(Pheno_tcga), rownames(pheno2))
pheno2 <- pheno2[CommonSamples, ]
all(rownames(pheno2) == rownames(CoxData_tcga))

# add the gleason score to CoxData_tcga
CoxData_tcga$gleason <- pheno2$Radical.Prostatectomy.Gleason.Score.for.Prostate.Cancer
CoxData_tcga$gleason <- as.factor(CoxData_tcga$gleason)
levels(CoxData_tcga$gleason)

colnames(CoxData_tcga)[colnames(CoxData_tcga) == 'tcga_prob_logReg'] <- 'PRN signature'

# fit the multivariate COX with gleason as cofactor 
Fit_sig_tcga_PFS_coxph_logReg_withGS <- coxph(Surv(Progress.Free.Survival..Months., Progression.Free.Status) ~ `PRN signature` + gleason, data = CoxData_tcga)
summary(Fit_sig_tcga_PFS_coxph_logReg_withGS)

tiff('./figures/survival/multivariateCox.tiff', width = 2500, height = 3000, res = 400)
ggforest(Fit_sig_tcga_PFS_coxph_logReg_withGS, fontsize = 1, main = 'HR')
dev.off()

###################
# fit multivariate COX with gleason and tumor purity as cofactor 
Fit_sig_tcga_PFS_coxph_logReg_withGS_purity <- coxph(Surv(Progress.Free.Survival..Months., Progression.Free.Status) ~ `PRN signature` + gleason + purity, data = CoxData_tcga)
summary(Fit_sig_tcga_PFS_coxph_logReg_withGS_purity)

tiff('./figures/survival/multivariateCox_gleason_purity.tiff', width = 2500, height = 3000, res = 400)
ggforest(Fit_sig_tcga_PFS_coxph_logReg_withGS_purity, fontsize = 1, main = 'HR')
dev.off()

