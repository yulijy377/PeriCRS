# Loading the Working Environment ------------------------------------------------------------------
library(plyr)
library(GSVA)
library(readxl)
library(reshape2)
library(tidyverse)
library(ggpubr)
library(ggthemes)
library(ggsci)
library(ggplot2)
library(pheatmap) 
library(RColorBrewer)
library(ConsensusClusterPlus)
library(survival)
library(survminer)
library(export)
library(tidyverse)
library(limma)
library(clusterProfiler)
library(org.Hs.eg.db)
library(glmnet)
library(survival)
library(dplyr)
library(data.table)
library(ggplot2) 
library(ggprism) 
library(ggrepel)
# Loading Data --------------------------------------------------------------------
circa_gene <- read_xlsx(path = "cicra_gene.xlsx")
circa_gene2 <- as.data.frame(circa_gene$`gene symbol`)
load(file="Transcriptome_cohort.Rdata")
# Quantification of Periodontitis Circadian Rhythm Score (PeriCRS) -------------------------------------------------------------------
clock_control_gene <- read_xlsx(path = "clock_control.xlsx")
clock_core_gene <-  read_xlsx(path = "clock_core.xlsx")
Hgeneset <- GSEABase::getGmt("h.all.v2023.1.Hs.symbols.gmt.txt")
exp_control <- exp_train_anno_final[rownames(exp_train_anno_final)%in%clock_control_gene$ID,]
exp_core <- exp_train_anno_final[rownames(exp_train_anno_final)%in%clock_core_gene$ID,]
clock_control_gsva <- gsva(as.matrix(exp_control),Hgeneset,method = "ssgsea",kcdf="Gaussian")
clock_core_gsva <- gsva(as.matrix(exp_core),Hgeneset,method = "ssgsea",kcdf="Gaussian",min.sz=1)
samepath <- intersect(rownames(clock_control_gsva),rownames(clock_core_gsva))
identical(colnames(clock_control_gsva),colnames(clock_core_gsva))
CRS1 <- clock_control_gsva[samepath,]
CRS2 <- clock_core_gsva[samepath,]
identical(colnames(CRS1),colnames(CRS2))
identical(rownames(CRS1),rownames(CRS2))
CRS <- CRS1-CRS2
circadian_rhythm_score <- as.data.frame(t(as.data.frame(colSums(CRS))))
rownames(circadian_rhythm_score) <- "circadian_rhythm_score"

# Constructing circadian hierarchy framework in transcriptomic cohorts ---------------------------------------------------------------------
circadian_rhythm_score2 <- sweep(circadian_rhythm_score,1,apply(circadian_rhythm_score,1,median))
input_path <- "003.ConsensusClusterPlus/"
fig_path <- "003.ConsensusClusterPlus/"
circadian_rhythm_score3 <- circadian_rhythm_score2
apply(circadian_rhythm_score3, 2, class)
circadian_rhythm_score3 %>% data.frame() %>% mutate(across(where(is.character), as.numeric))  %>% as.matrix() -> circadian_rhythm_score4
circadian_rhythm_score2 <- as.matrix(circadian_rhythm_score2)
circadian_rhythm_model <- ConsensusClusterPlus(circadian_rhythm_score2, 
                                               maxK = 10, 
                                               reps = 1000, 
                                               pItem = 0.8, 
                                               pFeature = 1, 
                                               clusterAlg = "hc", 
                                               corUse = "pairwise.complete.obs",
                                               distance='euclidean',
                                               title="circadian_rhythm_score",
                                               seed=123, 
                                               plot="png",
                                               writeTable=T)
# Analyzing Clinical Characteristics of Disrupted and Homeostasis Circadian Rhythm subgroups----------------------------------------------------------
kgroup <- read.csv(file = "circadian_rhythm_score/circadian_rhythm_score.k=3.consensusClass.csv",header = F)
table(kgroup$V2)
PCAPlot(kgroup)
colnames(kgroup) <- c("GSM","subgroup")
finalgroup <- merge(kgroup,pd,by="GSM")
table(finalgroup[,2:3])

#Differentially Expressed Genes (DEGs) Screening
subperio <- dplyr::filter(finalgroup,group== "Periodontitis" & !subgroup == 2)
subperio$subgroup <- factor(subperio$subgroup,
                            levels = c(3,1))
exp_subperio <- exp_train_anno_final[,colnames(exp_train_anno_final)%in%subperio$GSM]
identical(colnames(exp_subperio), subperio$GSM)
design <- model.matrix(~ subgroup, data = subperio)
fit=lmFit(exp_subperio,design)
fit=eBayes(fit)
res <- decideTests(fit, p.value=0.05)
summary(res)
DEG=topTable(fit,coef=2,n=Inf)
DEG=na.omit(DEG)
DEG_subperio <- rownames_to_column(DEG)
# Functional Enrichment Analysis of DEGs------------------------------------------------------------------
upgene <- filter(DEG_subperio,logFC > 0.5 & adj.P.Val < 0.01)
up_entrezid <- bitr(geneID = upgene$rowname
                    , fromType = "SYMBOL" 
                    , toType = "ENTREZID"
                    , OrgDb = "org.Hs.eg.db"
)
dngene <- filter(DEG_subperio,logFC < -0.5 & adj.P.Val < 0.01)
dn_entrezid <- bitr(geneID = dngene$rowname
                    , fromType = "SYMBOL" 
                    , toType = "ENTREZID" 
                    , OrgDb = "org.Hs.eg.db"
)
kkup <- enrichKEGG(gene = up_entrezid$ENTREZID, organism = "hsa", pvalueCutoff =0.05, qvalueCutoff =0.05)
kkdn <- enrichKEGG(gene = dn_entrezid$ENTREZID, organism = "hsa", pvalueCutoff =0.05, qvalueCutoff =0.05)
barplot(kkup, drop = TRUE, showCategory = 20)
barplot(kkdn, drop = TRUE, showCategory = 20)


# Identification of Key Molecular Clock Gene Using Machine Learning -------------------------------------------------------------------
CRSrisk <- finalgroup
CRSrisk <- filter(CRSrisk,!subgroup==2)
Lasso_sample <- CRSrisk
Lasso_sample <- filter(CRSrisk,group=="Periodontitis")
Lasso_sample <- Lasso_sample[,1:2]
colnames(Lasso_sample) <- c("sample","CRS_Risk")
clock_control_gene <-  read_xlsx(path = "clock_control.xlsx")
clock_core <- read_xlsx(path="clock_core.xlsx")
clock_all <- rbind(clock_control_gene,clock_core)
Lasso_exp <- exp_train_anno_final[rownames(exp_train_anno_final)%in%clock_all$ID,colnames(exp_train_anno_final)%in%Lasso_sample$sample]

Lasso_data <- as.data.frame(t(Lasso_exp))
Lasso_data <- rownames_to_column(Lasso_data,var = "sample")
Lasso_data <- merge(Lasso_sample,Lasso_data)
rownames(Lasso_data) <- Lasso_data$sample
Lasso_data <- Lasso_data[,-1]
lasso_data <- Lasso_data
set.seed(1234)
x <- as.matrix(lasso_data[ ,c(2:ncol(lasso_data))])
y <- lasso_data$CRS_Risk
y <- as.factor(lasso_data$CRS_Risk)
alpha1_fit <- glmnet(x, y, alpha = 1, family = "binomial", nlambda = 100)
plot(alpha1_fit, xvar = "lambda", label = TRUE)
set.seed(1234)
alpha1.fit.cv <- cv.glmnet(x, y, type.measure = "deviance", alpha = 1, family = "binomial")
plot(alpha1.fit.cv)
coef(alpha1.fit.cv, s = alpha1.fit.cv$lambda.1se)
feature_all <- as.data.frame(as.matrix(coef(alpha1.fit.cv, s = alpha1.fit.cv$lambda.1se)))
colnames(feature_all) <- "coff"
feature_opt <-  feature_all %>% filter(abs(coff) > 0)
rownames(feature_opt)
DEG_lasso <- DEG_subperio[DEG_subperio$rowname %in% rownames(feature_opt),]
Target_gene <- filter(DEG_lasso,logFC > 0.5 & adj.P.Val < 0.01|logFC < -0.5 & adj.P.Val < 0.01)
# Construction of Target_Genes'PPI Network  -------------------------------





