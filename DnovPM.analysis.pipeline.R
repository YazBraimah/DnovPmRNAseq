# ____                   ____            ____  _   _    _                   
#|  _ \ _ __   _____   _|  _ \ _ __ ___ |  _ \| \ | |  / \   ___  ___  __ _ 
#| | | | '_ \ / _ \ \ / / |_) | '_ ` _ \| |_) |  \| | / _ \ / __|/ _ \/ _` |
#| |_| | | | | (_) \ V /|  __/| | | | | |  _ <| |\  |/ ___ \\__ \  __/ (_| |
#|____/|_| |_|\___/ \_/ |_|   |_| |_| |_|_| \_\_| \_/_/   \_\___/\___|\__, |
#                                                                        |_|

## Load packages:
req_packages = c("Biobase", "cluster", "cowplot", "cummeRbund", 
                 "data.table", "DESeq", "edgeR", "ggplot2", "ggrepel", 
                 "ggthemes", "GO.db", "goseq", "grid", "gridExtra", 
                 "plotly", "qvalue", "reshape", "Rmisc", "splitstackshape", 
                 "statmod", "VennDiagram", "ggthemr")

lapply(req_packages, require, character.only = TRUE)

## The Cowplot package changes the default themes of ggplot2. Set the default theme like so:
theme_set(theme_gray())

## Load functions:
source("Functions2.R")

## Load annotations:
grpTrinotate = read.csv("Annotations/Trinotate_report_dvir1.06_subset.txt", header = T, sep = "\t", na.strings = ".", stringsAsFactors=FALSE)
GO_info = read.table("Annotations/Trinotate_report_dvir1.06_gene_ontology.txt", header=F, row.names=1,stringsAsFactors=F)

gffRecord = read.table("Annotations/FBgn_ID_name_coordinates.txt", header = T)

melOrths = read.table(file = "Annotations/mel_orths.txt", header = T)
melOrthsAll = aggregate(mel_GeneSymbol~FBgn_ID, data = melOrths, toString)

Annots = merge(merge(melOrthsAll, grpTrinotate, all=TRUE), gffRecord, all=TRUE)

## Load PAML data

tmp.FB.names = unique(subset(Annots, select=c("FBgn_ID", "FBtr_ID")))
paml.data = read.csv(file = "Annotations/PAML.branchSite.ALL.results.txt", header = T, sep = "\t")
paml.data = merge(tmp.FB.names, paml.data, all=T)
paml.data = merge(gffRecord, paml.data, all=T)
KaKs.data = read.csv(file = "Annotations/KaKs.ALL.results.txt", header = T, sep = "\t", check.names = F)
KaKs.data = merge(tmp.FB.names, KaKs.data, all=T)
KaKs.data = merge(gffRecord, KaKs.data, all=T)

## Read in expression data
DnovPM.dvir1.06.CountsMatrix = read.table("ExpressionData/genes_DnovPM_dvi1.06.counts.matrix", header=T, row.names=1, com='', check.names=F)
DnovPM.dvir1.06.TpmMatrix.cbmt = read.table("ExpressionData/genes_DnovPM_dvi1.06.TPM.not_cross_norm.counts_by_min_TPM", header = T)
DnovPM.dvir1.06.TmmMatrix = read.table("ExpressionData/genes_DnovPM_dvi1.06.TMM.EXPR.matrix", header=T, row.names=1, com='', check.names=F)

# Read in sample info
DnovPM.Samples_data = read.table("ExpressionData/samples.txt", header=F, check.names=F, fill=T)
DnovPM.Samples_data = DnovPM.Samples_data[DnovPM.Samples_data[,2] != '',]

## Barplot of gene counts by sample
libSizes <- as.data.frame(colSums(DnovPM.dvir1.06.CountsMatrix))
libSizes <- cbind(sample = row.names(libSizes), libSizes)
row.names(libSizes)<- NULL
colnames(libSizes) = c("sample", "Total_reads")
ggplot(libSizes, aes(sample, Total_reads)) + 
    geom_bar(stat="identity") + 
    theme(axis.text.x = element_text(angle = -90, hjust = 0)) + 
    geom_hline(yintercept = 20000000)


## Boxplot of log10(TPM) across all samples
m.expData<-melt(as.matrix(DnovPM.dvir1.06.TmmMatrix))
colnames(m.expData) <- c("gene_id", "replicate", "TPM")
m.expData.exp<- within(m.expData, replicate<-data.frame(do.call('rbind', strsplit(as.character(replicate),'_',fixed=TRUE))))
m.expData<-data.frame(m.expData, m.expData.exp$replicate$X1, m.expData.exp$replicate$X2, m.expData.exp$replicate$X3)
colnames(m.expData) <- c("gene_id", "replicate", "TPM", "sample", "tissue", "rep_num")
m.expData$TPM <- m.expData$TPM + 1
ggplot(m.expData) + 
    geom_boxplot(aes(x = replicate, y = log10(TPM), fill = sample), size = 0.3, alpha = I(1/3)) + 
    facet_wrap(~tissue, scales = "free_x") +
    theme(axis.text.x = element_text(angle = -90, hjust = 0)) + 
    scale_fill_hue(l = 50, h.start = 200)

## Estimate of the number of expressed genes (Brian Haas' method)
# extract the data between 10 TPM and 100 TPM
DnovPM_filt_data = DnovPM.dvir1.06.TpmMatrix.cbmt[DnovPM.dvir1.06.TpmMatrix.cbmt[,1] > -100 & DnovPM.dvir1.06.TpmMatrix.cbmt[,1] < -10,]
# perform a linear regression on this filtered subset of the data
DnovPM_fit = lm(DnovPM_filt_data[,2] ~ DnovPM_filt_data[,1])
print(DnovPM_fit)
# plot it
ggplot(DnovPM.dvir1.06.TpmMatrix.cbmt, aes(neg_min_tpm,num_features)) + 
    geom_point() + 
    scale_x_continuous(limits=c(-100,0)) + 
    scale_y_continuous(limits=c(0,20000)) + 
    geom_smooth(data=DnovPM_filt_data, method = "lm") + 
    geom_hline(yintercept = 9185, colour = "green") + ggtitle("expressed genes", subtitle = "something")


## Using the minimum CPM method to filter-out low expression genes
all.CPM <- cpm(DnovPM.dvir1.06.CountsMatrix)
thresh <- all.CPM > 1
table(rowSums(thresh))
keep <- rowSums(thresh) >= 2
counts.keep <- DnovPM.dvir1.06.CountsMatrix[keep,]
dim(counts.keep)

DnovPM_CountsMatrix_RT = subset(DnovPM.dvir1.06.CountsMatrix, select=grepl("RT", colnames(DnovPM.dvir1.06.CountsMatrix)))
head(DnovPM_CountsMatrix_RT)
RT.CPM <- cpm(DnovPM_CountsMatrix_RT)
thresh <- RT.CPM > 0.5
table(rowSums(thresh))
keep <- rowSums(thresh) >= 3
RT.counts.keep <- DnovPM_CountsMatrix_RT[keep,]
dim(RT.counts.keep)

RT.group <- factor(c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6,7,7,7))
RT.design <- model.matrix(~0+RT.group)
colnames(RT.design)<-c("C12_RT", "C3_RT", "C6_RT", "H12_RT", "H3_RT", "H6_RT", "V_RT")
DnovPM_DGElist_RT<-DGEList(counts = RT.counts.keep, group = RT.group)
DnovPM_DGElist_RT<-calcNormFactors(DnovPM_DGElist_RT)
DnovPM_DGElist_RT<-estimateDisp(DnovPM_DGElist_RT, RT.design, robust = T)
DnovPM_RT_fit <- glmFit(DnovPM_DGElist_RT, RT.design)

plotBCV(DnovPM_DGElist_RT)
g <- gof(DnovPM_RT_fit)
z <- zscoreGamma(g$gof.statistics,shape=g$df/2,scale=2)
qqnorm(z); qqline(z, col = 4,lwd=4,lty=4)

head(g$outlier)


## calculate dispersion
d <- DGEList(counts = counts.keep, group = DnovPM.Samples_data$V1)
d <- calcNormFactors(d)
d <- estimateCommonDisp(d)
d <- estimateTagwiseDisp(d)
summary(d$tagwise.dispersion)
## Plot biological coefficient of variation
plotBCV(d)
## Plot grouping of samples
plotMDS(d, method = "bcv", col=as.numeric(d$samples$group))

## Plot sample correlation
data = log2(DnovPM.dvir1.06.CountsMatrix+1)
data = as.matrix(data)
sample_cor = cor(data, method='pearson', use='pairwise.complete.obs')
sample_dist = as.dist(1-sample_cor)
hc_samples = hclust(sample_dist, method='complete')
heatmap.3(sample_cor, dendrogram='both', Rowv=as.dendrogram(hc_samples), Colv=as.dendrogram(hc_samples), col = colorpanel(75, '#dd70cd','black','#afc64f'), scale='none', symm=TRUE, key=TRUE,density.info='none', trace='none', symkey=FALSE, symbreaks=F, margins=c(10,10), cexCol=1, cexRow=1, cex.main=0.75, main=paste("sample correlation matrix"))

# normalize by DESeq method:
meta <- data.frame(row.names=colnames(DnovPM.dvir1.06.CountsMatrix), condition=DnovPM.Samples_data$V1)
DnovPMCountData<-round(DnovPM.dvir1.06.CountsMatrix)
DnovPMCountData_normByDESeq = newCountDataSet(DnovPMCountData, meta)
DnovPMCountData_normByDESeq = estimateSizeFactors(DnovPMCountData_normByDESeq)
DnovPMCountData_normByDESeq = data.frame(counts(DnovPMCountData_normByDESeq, normalized=T))

MA_BPlot(DnovPMCountData_normByDESeq, "C6_RT_3", "C6_RT_1")

## Filter count data by minimum count across ANY sample
DnovPM_max_gene_expr_per_row = apply(DnovPM.dvir1.06.CountsMatrix, 1, max)
DnovPM.dvir1.06.CountsMatrix.min400count = DnovPM.dvir1.06.CountsMatrix[DnovPM_max_gene_expr_per_row >= 400,,drop=F ]

##########################################################################################
################# Remove "bad" replicates  for DE analysis ###############################
### Based on the QC analysis above, some replicates show inconsistencies that are likely due to cross tissue contamination
### during dissections. The DE analysis will exclude these replicates

## Define good replicates (propper replicate grouping and correlation)
DnovPM.GoodReps = as.character(subset(DnovPM.Samples_data, V2 != "H3_RT_1" & V2 != "H6_RT_1" & V2 != "C6_RT_1")$V2)
DnovPM.GoodSamples = subset(DnovPM.Samples_data, V2 != "H3_RT_1" & V2 != "H6_RT_1" & V2 != "C6_RT_1")

## Create counts matrix with good replicates only
DnovPM.dvir1.06.CountsMatrix.BRR=subset(DnovPM.dvir1.06.CountsMatrix, select=DnovPM.GoodReps)

## Create normalized TPM matrix with good replicates only
DnovPM.dvir1.06.TmmMatrix.BRR=subset(DnovPM.dvir1.06.TmmMatrix, select=DnovPM.GoodReps)

## Rename columns to keep replicate order
# count matrices
colnames(DnovPM.dvir1.06.CountsMatrix.BRR) = DnovPM.GoodReps

# TPM matrices
colnames(DnovPM.dvir1.06.TmmMatrix.BRR) = colnames(DnovPM.dvir1.06.CountsMatrix.BRR)


#########################################################################################
### Summary TPM table and matrix for gene level plots (includes good replicates only) ######### 

DnovPM.TPM.tmp<-DnovPM.dvir1.06.TmmMatrix.BRR
colnames(DnovPM.TPM.tmp) <- DnovPM.GoodSamples$V1
m.DnovPM.TPM.tmp <- as.data.frame(melt(as.matrix(DnovPM.TPM.tmp)))
m.DnovPM.TPM.tmp <- within(m.DnovPM.TPM.tmp, X2<-data.frame(do.call('rbind', strsplit(as.character(X2),'_',fixed=TRUE))))
m.DnovPM.TPM.tmp<-data.frame(m.DnovPM.TPM.tmp$X1, m.DnovPM.TPM.tmp$X2$X1, m.DnovPM.TPM.tmp$X2$X2, m.DnovPM.TPM.tmp$value)
colnames(m.DnovPM.TPM.tmp) <- c("FBgn_ID", "sample", "tissue", "TPM")
m.DnovPM.TPM.tmp$condition <- ifelse(grepl("C", m.DnovPM.TPM.tmp$sample, ignore.case = F), "conspecific", ifelse(grepl("H", m.DnovPM.TPM.tmp$sample, ignore.case = F), "heterospecific", "virgin"))
m.DnovPM.TPM.tmp$time <- ifelse(grepl("3", m.DnovPM.TPM.tmp$sample), "3hpm", ifelse(grepl("6", m.DnovPM.TPM.tmp$sample), "6hpm", ifelse(grepl("12", m.DnovPM.TPM.tmp$sample), "12hpm","virgin")))
m.DnovPM.TPM.tmp.c = summarySE(m.DnovPM.TPM.tmp, measurevar = "TPM", groupvars = c("FBgn_ID", "sample", "tissue", "condition", "time"))
fbgn_to_geneName<-subset(gffRecord, select=c(FBgn_ID, gene_name))
TPMse_DnovPM <- merge(fbgn_to_geneName, m.DnovPM.TPM.tmp.c, all=TRUE)
DnovPM_MeanTPMmatrix<-cast(m.DnovPM.TPM.tmp.c, FBgn_ID~sample+tissue, value ="TPM")
TPMse_DnovPM$condition = factor (TPMse_DnovPM$condition, levels = c("virgin", "conspecific", "heterospecific"))
TPMse_DnovPM$time = factor (TPMse_DnovPM$time, levels = c("virgin", "3hpm", "6hpm", "12hpm"))

## plot a gene's expression like this:
plotGenePM(TPMse_DnovPM, "snRNA:U6:3")
plotGeneG(TPMse, "snRNA:U6:3")
## Generate tissue specificity matrix:
Dnov_virgin_tissue_MeanTPMmatrix <- subset(DnovPM_MeanTPMmatrix, select=c(FBgn_ID, V_CR, V_H, V_OV, V_RT))
rownames(Dnov_virgin_tissue_MeanTPMmatrix) <- Dnov_virgin_tissue_MeanTPMmatrix[,1]
Dnov_virgin_tissue_MeanTPMmatrix[,1] <- NULL
Dnov_virgin_Specificity_table = calcSpecificity(Dnov_virgin_tissue_MeanTPMmatrix)
Dnov_virgin_Specificity_table = as.data.frame(Dnov_virgin_Specificity_table)

########################################## edgeR analysis ###################################################
## Define groups for DE analysis.
# simple way (for pair-wise comparisons):

## calculate dispersion
DnovPM.group <- factor(c(1,1,1,2,2,2,3,3,4,4,5,5,6,6,6,7,7,8,8,9,9,10,10,11,11,12,12,13,13,14,14,14))
DnovPM.design <- model.matrix(~0+DnovPM.group)
colnames(DnovPM.design)<-levels(DnovPM.GoodSamples$V1)

## Filter count data by minimum count across ANY sample
DnovPM_max_gene_expr_per_row = apply(DnovPM.dvir1.06.CountsMatrix.BRR, 1, max)
DnovPM.dvir1.06.CountsMatrix.BRR.min400count = DnovPM.dvir1.06.CountsMatrix.BRR[DnovPM_max_gene_expr_per_row >= 400,,drop=F ]

DnovPM_ExpStd<-DGEList(counts = DnovPM.dvir1.06.CountsMatrix.BRR.min400count, group = DnovPM.group)
DnovPM_ExpStd<-calcNormFactors(DnovPM_ExpStd)
DnovPM_ExpStd<-estimateDisp(DnovPM_ExpStd, DnovPM.design, robust = T)

DnovPM_fit <- glmFit(DnovPM_ExpStd, DnovPM.design)

virgin_RT_contrasts<- makeContrasts(V_RT.vs.V_CR=V_RT-V_CR, 
                                        V_RT.vs.V_H=V_RT-V_H,
                                        V_RT.vs.V_OV=V_RT-V_OV,
                                        levels=DnovPM.design)

virgin_OV_contrasts<- makeContrasts(V_OV.vs.V_CR=V_OV-V_CR, 
                                    V_OV.vs.V_H=V_OV-V_H,
                                    V_OV.vs.V_RT=V_OV-V_RT,
                                    levels=DnovPM.design)

virgin_H_contrasts<- makeContrasts(V_H.vs.V_CR=V_H-V_CR, 
                                    V_H.vs.V_OV=V_H-V_OV,
                                    V_H.vs.V_RT=V_H-V_RT,
                                    levels=DnovPM.design)

## Find RT-biased genes
lrt.RT.v.rest <- glmLRT(DnovPM_fit, contrast = virgin_RT_contrasts)
lrt.RT.v.rest.tTags <- topTags(lrt.RT.v.rest, n = NULL)
lrt.RT.v.rest.tTags.table <- lrt.RT.v.rest.tTags$table
Dnov.dvir1.06.RT.list<-rownames(subset(lrt.RT.v.rest.tTags.table, logFC.V_RT.vs.V_CR > 2 & logFC.V_RT.vs.V_H > 2 & logFC.V_RT.vs.V_OV > 2 & FDR<0.001))

## Find OV-biased genes
lrt.OV.v.rest <- glmLRT(DnovPM_fit, contrast = virgin_OV_contrasts)
lrt.OV.v.rest.tTags <- topTags(lrt.OV.v.rest, n = NULL)
lrt.OV.v.rest.tTags.table <- lrt.OV.v.rest.tTags$table
Dnov.dvir1.06.OV.list<-rownames(subset(lrt.OV.v.rest.tTags.table, logFC.V_OV.vs.V_CR > 2 & logFC.V_OV.vs.V_H > 2 & logFC.V_OV.vs.V_RT > 2 & FDR<0.001))

## Find H-biased genes
lrt.H.v.rest <- glmLRT(DnovPM_fit, contrast = virgin_H_contrasts)
lrt.H.v.rest.tTags <- topTags(lrt.H.v.rest, n = NULL)
lrt.H.v.rest.tTags.table <- lrt.H.v.rest.tTags$table
Dnov.dvir1.06.H.list<-rownames(subset(lrt.H.v.rest.tTags.table, logFC.V_H.vs.V_CR > 2 & logFC.V_H.vs.V_OV > 2 & logFC.V_H.vs.V_RT > 2 & FDR<0.001))


Dnov.dvir1.06.RT.matrix <- subset(Dnov_virgin_tissue_MeanTPMmatrix, rownames(Dnov_virgin_tissue_MeanTPMmatrix) %in% Dnov.dvir1.06.RT.list)
plotHeatmap(Dnov.dvir1.06.RT.matrix, clustering = "both", labRow = T)

# pdf("Plots/Dnov_RT-biased_barPlots.pdf", height = 4)
# lapply(Dnov.dvir1.06.RT.list, plotGenePM, object = TPMse_DnovPM)
# dev.off()

ggplot(subset(paml.data, FBgn_ID %in% Dnov.dvir1.06.RT.list & omega < 800 & grepl("Chr", chromosome)), aes(max, omega)) + 
    geom_point(size=2, alpha=0.5, colour = "#7d49c3") + 
    geom_point(data=subset(paml.data, FBgn_ID %in% SFP_elements$`D.ame,D.lum,D.nov,D.vir` ), aes(max, omega), inherit.aes = F, size=2, alpha=0.5, colour = "#4f922a") + 
    geom_hline(yintercept = 0.15, linetype="dashed", colour = "yellow") + 
    geom_hline(yintercept = 1, linetype="dashed", colour = "gray")  + 
    facet_grid(~chromosome, scales = "free_x") + 
    scale_colour_manual(name = "", values =c("#7aa457"="#7aa457","#9e6ebd"="#9e6ebd"), labels = c("SFPs","EB biased")) + 
    scale_x_continuous(breaks=c(5000000, 10000000, 15000000, 20000000, 25000000, 30000000), labels=expression("5", "10", "15", "20", "25", "30")) + 
    xlab ("Chromosome coordinates (Mb)") + 
    labs(y=expression(K[a]/K[s])) + 
    geom_text_repel(data=subset(paml.data, FBgn_ID %in% Dnov.dvir1.06.RT.list & omega > 0.95 & omega < 800), aes(label = gene_name), size =3, force = 30, colour = "#7d49c3") +
    geom_text_repel(data=subset(paml.data, gene_name == "GJ21515"), aes(label = gene_name), size =3, force = 30, colour = "#4f922a") +
    geom_text_repel(data=subset(paml.data, FBgn_ID %in% SFP_elements$`D.ame,D.lum,D.nov,D.vir` & omega > 0.8), aes(label = gene_name), size =3, force = 4, colour = "#4f922a") +
    theme(axis.title.x = element_text(face = "bold", size = 10, vjust=0.1), axis.text.x=element_text(face = "bold", size = 12),axis.text.y = element_text(face = "bold", size = 12), axis.title.y = element_text(face = "bold.italic", size = 12, vjust=0.1), strip.text=element_text(face="bold", size = 12))


########### RT contrasts
# create RT-specific count matrix and define groups/design
DnovPM_CountsMatrix_RT = subset(DnovPM.dvir1.06.CountsMatrix.BRR.min400count, select=grepl("RT", colnames(DnovPM.dvir1.06.CountsMatrix.BRR.min400count)))

## filter this matrix by minimum read count
DnovPM_max_gene_expr_per_row = apply(DnovPM_CountsMatrix_RT, 1, max)
DnovPM_CountsMatrix_RT.min400 = DnovPM_CountsMatrix_RT[DnovPM_max_gene_expr_per_row >= 200,,drop=F ]

######
RT.group <- factor(c(1,1,1,2,2,2,3,3,4,4,4,5,5,6,6,7,7,7))
RT.design <- model.matrix(~0+RT.group)
colnames(RT.design)<-c("C12_RT", "C3_RT", "C6_RT", "H12_RT", "H3_RT", "H6_RT", "V_RT")

# create edgeR DE object and run glm
DnovPM_DGElist_RT<-DGEList(counts = DnovPM_CountsMatrix_RT.min400, group = RT.group)
DnovPM_DGElist_RT<-calcNormFactors(DnovPM_DGElist_RT)
DnovPM_DGElist_RT<-estimateDisp(DnovPM_DGElist_RT, RT.design, robust = T)
DnovPM_RT_fit <- glmFit(DnovPM_DGElist_RT, RT.design)

# define contrasts
con_virgin_contrasts <- makeContrasts(C3.vs.V=C3_RT-V_RT, C6.vs.V=C6_RT-V_RT, C12.vs.V=C12_RT-V_RT, levels = RT.design)
het_virgin_contrasts <- makeContrasts(H3.vs.V=H3_RT-V_RT, H6.vs.V=H6_RT-V_RT, H12.vs.V=H12_RT-V_RT, levels = RT.design)

# identify overall con_vs_het
RT_con.vs.virgin <- glmLRT(DnovPM_RT_fit, contrast = con_virgin_contrasts)
RT_con.vs.virgin.tTags <- topTags(RT_con.vs.virgin, n = NULL)
RT_con.vs.virgin.tTags.table <- RT_con.vs.virgin.tTags$table
RT_con.vs.virgin.tTags.table$FBgn_ID <- rownames(RT_con.vs.virgin.tTags.table)

RT_het.vs.virgin <- glmLRT(DnovPM_RT_fit, contrast = het_virgin_contrasts)
RT_het.vs.virgin.tTags <- topTags(RT_het.vs.virgin, n = NULL)
RT_het.vs.virgin.tTags.table <- RT_het.vs.virgin.tTags$table
RT_het.vs.virgin.tTags.table$FBgn_ID <- rownames(RT_het.vs.virgin.tTags.table)

# Pairwise contrasts
# identify DE genes: C3_RT vs V_RT
RT_con.3hrs.vs.virgin <- glmLRT(DnovPM_RT_fit, contrast = con_virgin_contrasts[,"C3.vs.V"])
RT_con.3hrs.vs.virgin.tTags <- topTags(RT_con.3hrs.vs.virgin, n = NULL)
RT_con.3hrs.vs.virgin.tTags.table <- RT_con.3hrs.vs.virgin.tTags$table
RT_con.3hrs.vs.virgin.Up.list <- rownames(subset(RT_con.3hrs.vs.virgin.tTags.table, logFC > 1 & FDR < 0.001))
RT_con.3hrs.vs.virgin.Down.list <- rownames(subset(RT_con.3hrs.vs.virgin.tTags.table, logFC < -1 & FDR < 0.001))
# identify DE genes: C6_RT vs V_RT
RT_con.6hrs.vs.virgin <- glmLRT(DnovPM_RT_fit, contrast = con_virgin_contrasts[,"C6.vs.V"])
RT_con.6hrs.vs.virgin.tTags <- topTags(RT_con.6hrs.vs.virgin, n = NULL)
RT_con.6hrs.vs.virgin.tTags.table <- RT_con.6hrs.vs.virgin.tTags$table
RT_con.6hrs.vs.virgin.Up.list <- rownames(subset(RT_con.6hrs.vs.virgin.tTags.table, logFC > 1 & FDR < 0.001))
RT_con.6hrs.vs.virgin.Down.list <- rownames(subset(RT_con.6hrs.vs.virgin.tTags.table, logFC < -1 & FDR < 0.001))
# identify DE genes: C12_RT vs V_RT
RT_con.12hrs.vs.virgin <- glmLRT(DnovPM_RT_fit, contrast = con_virgin_contrasts[,"C12.vs.V"])
RT_con.12hrs.vs.virgin.tTags <- topTags(RT_con.12hrs.vs.virgin, n = NULL)
RT_con.12hrs.vs.virgin.tTags.table <- RT_con.12hrs.vs.virgin.tTags$table
RT_con.12hrs.vs.virgin.Up.list <- rownames(subset(RT_con.12hrs.vs.virgin.tTags.table, logFC > 1 & FDR < 0.001))
RT_con.12hrs.vs.virgin.Down.list <- rownames(subset(RT_con.12hrs.vs.virgin.tTags.table, logFC < -1 & FDR < 0.001))
# identify DE genes: H3_RT vs V_RT
RT_het.3hrs.vs.virgin <- glmLRT(DnovPM_RT_fit, contrast = het_virgin_contrasts[,"H3.vs.V"])
RT_het.3hrs.vs.virgin.tTags <- topTags(RT_het.3hrs.vs.virgin, n = NULL)
RT_het.3hrs.vs.virgin.tTags.table <- RT_het.3hrs.vs.virgin.tTags$table
RT_het.3hrs.vs.virgin.Up.list <- rownames(subset(RT_het.3hrs.vs.virgin.tTags.table, logFC > 1 & FDR < 0.001))
RT_het.3hrs.vs.virgin.Down.list <- rownames(subset(RT_het.3hrs.vs.virgin.tTags.table, logFC < -1 & FDR < 0.001))
# identify DE genes: H6_RT vs V_RT
RT_het.6hrs.vs.virgin <- glmLRT(DnovPM_RT_fit, contrast = het_virgin_contrasts[,"H6.vs.V"])
RT_het.6hrs.vs.virgin.tTags <- topTags(RT_het.6hrs.vs.virgin, n = NULL)
RT_het.6hrs.vs.virgin.tTags.table <- RT_het.6hrs.vs.virgin.tTags$table
RT_het.6hrs.vs.virgin.Up.list <- rownames(subset(RT_het.6hrs.vs.virgin.tTags.table, logFC > 1 & FDR < 0.001))
RT_het.6hrs.vs.virgin.Down.list <- rownames(subset(RT_het.6hrs.vs.virgin.tTags.table, logFC < -1 & FDR < 0.001))
# identify DE genes: H12_RT vs V_RT
RT_het.12hrs.vs.virgin <- glmLRT(DnovPM_RT_fit, contrast = het_virgin_contrasts[,"H12.vs.V"])
RT_het.12hrs.vs.virgin.tTags <- topTags(RT_het.12hrs.vs.virgin, n = NULL)
RT_het.12hrs.vs.virgin.tTags.table <- RT_het.12hrs.vs.virgin.tTags$table
RT_het.12hrs.vs.virgin.Up.list <- rownames(subset(RT_het.12hrs.vs.virgin.tTags.table, logFC > 1 & FDR < 0.001))
RT_het.12hrs.vs.virgin.Down.list <- rownames(subset(RT_het.12hrs.vs.virgin.tTags.table, logFC < -1 & FDR < 0.001))


RT_UP_3hrs_candidates <- list(Conspecific = RT_con.3hrs.vs.virgin.Up.list, Heterospecific = RT_het.3hrs.vs.virgin.Up.list)
RT_UP_6hrs_candidates <- list(Conspecific = RT_con.6hrs.vs.virgin.Up.list, Heterospecific = RT_het.6hrs.vs.virgin.Up.list)
RT_UP_12hrs_candidates <- list(Conspecific = RT_con.12hrs.vs.virgin.Up.list, Heterospecific = RT_het.12hrs.vs.virgin.Up.list)
RT_Down_3hrs_candidates <- list(Conspecific = RT_con.3hrs.vs.virgin.Down.list, Heterospecific = RT_het.3hrs.vs.virgin.Down.list)
RT_Down_6hrs_candidates <- list(Conspecific = RT_con.6hrs.vs.virgin.Down.list, Heterospecific = RT_het.6hrs.vs.virgin.Down.list)
RT_Down_12hrs_candidates <- list(Conspecific = RT_con.12hrs.vs.virgin.Down.list, Heterospecific = RT_het.12hrs.vs.virgin.Down.list)

# Rearrange into lists of lists, and partition by species
RT_UP_3hrs_combs <- unlist(lapply(1:length(RT_UP_3hrs_candidates), function(j) combn(names(RT_UP_3hrs_candidates), j, simplify = FALSE)), recursive = FALSE)
RT_UP_6hrs_combs <- unlist(lapply(1:length(RT_UP_6hrs_candidates), function(j) combn(names(RT_UP_6hrs_candidates), j, simplify = FALSE)), recursive = FALSE)
RT_UP_12hrs_combs <- unlist(lapply(1:length(RT_UP_12hrs_candidates), function(j) combn(names(RT_UP_12hrs_candidates), j, simplify = FALSE)), recursive = FALSE)
RT_Down_3hrs_combs <- unlist(lapply(1:length(RT_Down_3hrs_candidates), function(j) combn(names(RT_Down_3hrs_candidates), j, simplify = FALSE)), recursive = FALSE)
RT_Down_6hrs_combs <- unlist(lapply(1:length(RT_Down_6hrs_candidates), function(j) combn(names(RT_Down_6hrs_candidates), j, simplify = FALSE)), recursive = FALSE)
RT_Down_12hrs_combs <- unlist(lapply(1:length(RT_Down_12hrs_candidates), function(j) combn(names(RT_Down_12hrs_candidates), j, simplify = FALSE)), recursive = FALSE)

names(RT_UP_3hrs_combs) <- sapply(RT_UP_3hrs_combs, function(i) paste0(i, collapse = ","))
names(RT_UP_6hrs_combs) <- sapply(RT_UP_6hrs_combs, function(i) paste0(i, collapse = ","))
names(RT_UP_12hrs_combs) <- sapply(RT_UP_12hrs_combs, function(i) paste0(i, collapse = ","))
names(RT_Down_3hrs_combs) <- sapply(RT_Down_3hrs_combs, function(i) paste0(i, collapse = ","))
names(RT_Down_6hrs_combs) <- sapply(RT_Down_6hrs_combs, function(i) paste0(i, collapse = ","))
names(RT_Down_12hrs_combs) <- sapply(RT_Down_12hrs_combs, function(i) paste0(i, collapse = ","))

RT_UP_3hrs_elements <- lapply(RT_UP_3hrs_combs, function(i) Setdiff(RT_UP_3hrs_candidates[i], RT_UP_3hrs_candidates[setdiff(names(RT_UP_3hrs_candidates), i)]))
RT_UP_6hrs_elements <- lapply(RT_UP_6hrs_combs, function(i) Setdiff(RT_UP_6hrs_candidates[i], RT_UP_6hrs_candidates[setdiff(names(RT_UP_6hrs_candidates), i)]))
RT_UP_12hrs_elements <- lapply(RT_UP_12hrs_combs, function(i) Setdiff(RT_UP_12hrs_candidates[i], RT_UP_12hrs_candidates[setdiff(names(RT_UP_12hrs_candidates), i)]))
RT_Down_3hrs_elements <- lapply(RT_Down_3hrs_combs, function(i) Setdiff(RT_Down_3hrs_candidates[i], RT_Down_3hrs_candidates[setdiff(names(RT_Down_3hrs_candidates), i)]))
RT_Down_6hrs_elements <- lapply(RT_Down_6hrs_combs, function(i) Setdiff(RT_Down_6hrs_candidates[i], RT_Down_6hrs_candidates[setdiff(names(RT_Down_6hrs_candidates), i)]))
RT_Down_12hrs_elements <- lapply(RT_Down_12hrs_combs, function(i) Setdiff(RT_Down_12hrs_candidates[i], RT_Down_12hrs_candidates[setdiff(names(RT_Down_12hrs_candidates), i)]))

### Draw a VennDiagram of each element
RT_UP_3hrs_Vdiag<-venn.diagram(RT_UP_3hrs_candidates, NULL, fill=c("#b067a3", "#9c954d"), alpha=c(0.75,0.75), cex = 1.5, cat.fontface= 2, cat.cex = 0, resolution = 1000, main = "3hpm")
RT_UP_6hrs_Vdiag<-venn.diagram(RT_UP_6hrs_candidates, NULL, fill=c("#b067a3", "#9c954d"), alpha=c(0.75,0.75), cex = 1.5, cat.fontface= 2, cat.cex = 0, resolution = 1000, main = "6hpm")
RT_UP_12hrs_Vdiag<-venn.diagram(RT_UP_12hrs_candidates, NULL, fill=c("#b067a3", "#9c954d"), alpha=c(0.75,0.75), cex = 1.5, cat.fontface= 2, cat.cex = 0, resolution = 1000, main = "12hpm")
grid.arrange(gTree(children=RT_UP_3hrs_Vdiag), gTree(children=RT_UP_6hrs_Vdiag), gTree(children=RT_UP_12hrs_Vdiag))

RT_Down_3hrs_Vdiag<-venn.diagram(RT_Down_3hrs_candidates, NULL, fill=c("#b067a3", "#9c954d"), alpha=c(0.75,0.75), cex = 1.5, cat.fontface= 2, cat.cex = 0, resolution = 1000, main = "3hpm")
RT_Down_6hrs_Vdiag<-venn.diagram(RT_Down_6hrs_candidates, NULL, fill=c("#b067a3", "#9c954d"), alpha=c(0.75,0.75), cex = 1.5, cat.fontface= 2, cat.cex = 0, resolution = 1000, main = "6hpm")
RT_Down_12hrs_Vdiag<-venn.diagram(RT_Down_12hrs_candidates, NULL, fill=c("#b067a3", "#9c954d"), alpha=c(0.75,0.75), cex = 1.5, cat.fontface= 2, cat.cex = 0, resolution = 1000, main = "12hpm")
grid.arrange(gTree(children=RT_Down_3hrs_Vdiag), gTree(children=RT_Down_6hrs_Vdiag), gTree(children=RT_Down_12hrs_Vdiag))


# Test difference between con- and heterospecific
condition.group = subset(m.DnovPM.TPM.tmp, FBgn_ID == "FBgn0202928" & tissue == "RT")
condition.group$replicates = colnames(DnovPM_CountsMatrix_RT.min400)
rownames(condition.group) <- condition.group$replicates
condition.group = subset(condition.group, select = c(sample, condition, time))

condition.design <- model.matrix(~0+condition, data = condition.group)
colnames(condition.design) <- unique(condition.group$condition)

condition.contrasts <- makeContrasts(con.vs.het = conspecific-heterospecific,
                                     con.vs.vir = conspecific - virgin,
                                     het.vs.vir = heterospecific - virgin,
                                     levels = condition.design)

###### KEY ANALYSIS FOR DIFFERENCE BETWEEN CON and HET ###############
lrt <- glmFit(DnovPM_DGElist_RT, condition.design)
con.vs.het.RT.all.contrast <- glmLRT(lrt, contrast = condition.contrasts[,"con.vs.het"])
con.vs.het.RT.all.tTags <- topTags(con.vs.het.RT.all.contrast, n = NULL)
con.vs.het.RT.all.tTags.table <- con.vs.het.RT.all.tTags$table
con.vs.het.RT.het.Up.list <- rownames(subset(con.vs.het.RT.all.tTags.table, logFC < -1 & FDR < 0.001))
con.vs.het.RT.con.Up.list <- rownames(subset(con.vs.het.RT.all.tTags.table, logFC > 1 & FDR < 0.001))

con.vs.vir.RT.all.contrast <- glmLRT(lrt, contrast = condition.contrasts[,"con.vs.vir"])
con.vs.vir.RT.all.tTags <- topTags(con.vs.vir.RT.all.contrast, n = NULL)
con.vs.vir.RT.all.tTags.table <- con.vs.vir.RT.all.tTags$table
con.vs.vir.RT.con.Up.list <- rownames(subset(con.vs.vir.RT.all.tTags.table, logFC > 1 & FDR < 0.001))
con.vs.vir.RT.con.Down.list <- rownames(subset(con.vs.vir.RT.all.tTags.table, logFC < -1 & FDR < 0.001))

het.vs.vir.RT.all.contrast <- glmLRT(lrt, contrast = condition.contrasts[,"het.vs.vir"])
het.vs.vir.RT.all.tTags <- topTags(het.vs.vir.RT.all.contrast, n = NULL)
het.vs.vir.RT.all.tTags.table <- het.vs.vir.RT.all.tTags$table
het.vs.vir.RT.het.Up.list <- rownames(subset(het.vs.vir.RT.all.tTags.table, logFC > 1 & FDR < 0.001))
het.vs.vir.RT.het.Down.list <- rownames(subset(het.vs.vir.RT.all.tTags.table, logFC < -1 & FDR < 0.001))


PM.vs.vir_Up_candidates <- list(conspecific = con.vs.vir.RT.con.Up.list, heterospecific = het.vs.vir.RT.het.Up.list)
PM.vs.vir_Down_candidates <- list(conspecific = con.vs.vir.RT.con.Down.list, heterospecific = het.vs.vir.RT.het.Down.list)


PM.vs.vir_Up_candidates_Vdiag<-venn.diagram(PM.vs.vir_Up_candidates, NULL, fill=c("#b067a3", "#9c954d"), alpha=c(0.75,0.75), cex = 1.5, cat.fontface= 4, cat.cex = 0, resolution = 1000)
PM.vs.vir_Down_candidates_Vdiag<-venn.diagram(PM.vs.vir_Down_candidates, NULL, fill=c("#b067a3", "#9c954d"), alpha=c(0.75,0.75), cex = 1.5, cat.fontface= 4, cat.cex = 0, resolution = 1000)
grid.arrange(gTree(children=PM.vs.vir_Up_candidates_Vdiag))
grid.arrange(gTree(children=PM.vs.vir_Down_candidates_Vdiag))

pdf("con.vs.het.RT.het.Up.list.pdf", height = 4)
lapply(con.vs.het.RT.het.Up.list, plotGenePM_RT, object = TPMse_DnovPM)
dev.off()
pdf("con.vs.vir.RT.con.Up.list.pdf", height = 4)
lapply(con.vs.vir.RT.con.Up.list, plotGenePM_RT, object = TPMse_DnovPM)
dev.off()
pdf("het.vs.vir.RT.het.Up.list.pdf", height = 4)
lapply(het.vs.vir.RT.het.Up.list, plotGenePM_RT, object = TPMse_DnovPM)
dev.off()

C2inv.qtl = data.frame(xmin=17747413.5, xmax=34500000, ymin=0, ymax = 1.5, chromosome = "Chr_2")
C5.qtl = data.frame(xmin=c(13800000, 16300000, 22800000), xmax=c(14750000, 21700000, 25000000), ymin=0, ymax = 1.5, chromosome = "Chr_5")

Dnov.dvir1.06.RT.list.sigP = unique(subset(Annots, FBgn_ID %in% Dnov.dvir1.06.RT.list & SignalP == "YES"))$FBgn_ID

ggplot(subset(paml.data, FBgn_ID %in% Dnov.dvir1.06.RT.list & omega < 800 & grepl("Chr", chromosome)), aes(max, omega)) + 
    geom_point(size=2, alpha=0.75, aes(colour = "#7d49c3")) + 
    geom_point(data=subset(paml.data, FBgn_ID %in% SFP_elements$`D.ame,D.lum,D.nov,D.vir` ), aes(max, omega, colour = "#4f922a"), size=2, alpha=0.75) + 
    geom_hline(yintercept = 0.15, linetype="dashed", colour = "red") + 
    geom_hline(yintercept = 1, linetype="dashed", colour = "black") + 
    scale_colour_manual(name = "", values =c("#7d49c3"="#7d49c3","#4f922a"="#4f922a"), labels = c("SFPs","Female RT")) + 
    facet_grid(~chromosome, scales = "free_x") + 
    xlab ("Chromosome coordinates (Mb)") + 
    labs(y=expression(K[a]/K[s])) + 
    theme(axis.title.x = element_text(face = "bold", size = 10, vjust=0.1), axis.text.x=element_text(face = "bold", size = 12),axis.text.y = element_text(face = "bold", size = 12), axis.title.y = element_text(face = "bold.italic", size = 12, vjust=0.1), strip.text=element_text(face="bold", size = 12), legend.text = element_text(size = 12, face="bold")) + 
    geom_rect(data=C2inv.qtl, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill="red"), alpha=0.1, inherit.aes = FALSE) + 
    geom_rect(data=C5.qtl, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="red", alpha=0.1, inherit.aes = FALSE) + 
    scale_fill_manual(name = '', values = "red",labels = "PMPZ\nQTL region")+ scale_x_continuous(breaks=c(5000000, 10000000, 15000000, 20000000, 25000000, 30000000), labels=expression("5", "10", "15", "20", "25", "30")) + theme_bw()



## Only chromosome 2
ggplot(subset(paml.data, FBgn_ID %in% Dnov.dvir1.06.RT.list & omega < 800 & chromosome == "Chr_2"), aes(max, omega)) + 
    geom_point(size=2, alpha=0.75, aes(colour = "#7d49c3")) + 
    geom_point(data=subset(paml.data, FBgn_ID %in% SFP_elements$`D.ame,D.lum,D.nov,D.vir` ), aes(max, omega, colour = "#4f922a"), size=2, alpha=0.75) + 
    geom_hline(yintercept = 0.15, linetype="dashed", colour = "black") + 
    geom_hline(yintercept = 1, linetype="dashed", colour = "red") + 
    scale_colour_manual(name = "", values =c("#7d49c3"="#7d49c3","#4f922a"="#4f922a"), labels = c("SFPs","Female RT")) + 
    scale_x_continuous(breaks=c(5000000, 10000000, 15000000, 20000000, 25000000, 30000000), labels=expression("5", "10", "15", "20", "25", "30")) + 
    xlab ("Chromosome coordinates (Mb)") + 
    labs(y=expression(K[a]/K[s])) + 
    geom_text(data=subset(paml.data, FBgn_ID %in% Dnov.dvir1.06.RT.list & omega > 0.95 & omega < 800 & chromosome == "Chr_2"), aes(label = gene_name), size =4, colour = "#7d49c3", fontface = "bold.italic", vjust = 1.4, hjust = 0.1) + 
    geom_text(data=subset(paml.data, FBgn_ID %in% SFP_elements$`D.ame,D.lum,D.nov,D.vir` & omega > 0.8), aes(label = gene_name), size =4, colour = "#4f922a", fontface = "bold.italic", hjust = 1.1) + 
    theme(axis.title.x = element_text(face = "bold", size = 10, vjust=0.1), axis.text.x=element_text(face = "bold", size = 12),axis.text.y = element_text(face = "bold", size = 12), axis.title.y = element_text(face = "bold.italic", size = 12, vjust=0.1), strip.text=element_text(face="bold", size = 12), plot.title=element_text(face="bold"), legend.text = element_text(size = 12, face="bold")) + 
    ggtitle("Chromosome 2 (Muller element E)") + 
    geom_rect(data=C2inv.qtl, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="red", alpha=0.1, inherit.aes = FALSE) + scale_fill_manual(name = '', values = "red",labels = "PMPZ\nQTL region")


###########
###########

# create gene lists and factor labeling
RT_factors = as.data.frame(Dnov.dvir1.06.RT.list)
RT_factors$V1 = "RT-biased"
rownames(RT_factors) = Dnov.dvir1.06.RT.list
RT_factors = subset(RT_factors, select = "V1")

OV_factors = as.data.frame(Dnov.dvir1.06.OV.list)
OV_factors$V1 = "OV-biased"
rownames(OV_factors) = Dnov.dvir1.06.OV.list
OV_factors = subset(OV_factors, select = "V1")

factor_labeling = rbind(RT_factors, OV_factors)
colnames(factor_labeling) = c('type')
factor_list = unique(factor_labeling[,1])

cat_genes_vec = as.integer(features_with_GO %in% rownames(factor_labeling))
pwf=nullp(cat_genes_vec,bias.data=lengths_features_with_GO)
rownames(pwf) = names(GO_info_listed)

GO_enriched_list = list()

for (feature_cat in factor_list) {
    message('Processing category: ', feature_cat)
    cat_genes_vec = as.integer(features_with_GO %in% rownames(factor_labeling)[factor_labeling$type == feature_cat])
    pwf$DEgenes = cat_genes_vec
    res = goseq(pwf,gene2cat=GO_info_listed)
    ## over-represented categories:
    pvals = res$over_represented_pvalue
    pvals[pvals > 1 -1e-10] = 1-1e-10
    q = qvalue(pvals)
    res$over_represented_FDR = q$qvalues
    go_enrich_factor = feature_cat
    enrich_result_table = res[res$over_represented_pvalue<=0.05,]
    descr = unlist(lapply(enrich_result_table$category, get_GO_term_descr))
    enrich_result_table$go_term = descr
    enrich_result_table$factor = go_enrich_factor
    GO_enriched_list[[feature_cat]] = enrich_result_table
}

GO_enrichment_data = rbindlist(GO_enriched_list)


ggplot(subset(GO_enrichment_data, over_represented_FDR < 0.01 & factor == "RT-biased"), 
       aes(category, -log10(over_represented_pvalue), size = numDEInCat, colour = ontology)) + 
    geom_point()  + 
    xlab(NULL) + 
    geom_text_repel(data = subset(GO_enrichment_data, over_represented_FDR < 0.01 & factor == "RT-biased" & numDEInCat > 20), 
                    aes(category, -log10(over_represented_pvalue),label=term), 
                    force = 20, 
                    inherit.aes = F, 
                    box.padding = unit(0.35, "lines"), 
                    point.padding = unit(0.5, "lines"), 
                    fontface = "bold", 
                    size = 3) + 
    theme(axis.text.x = element_text(angle = 45, face = "bold", vjust = 1, hjust = 1)) + 
    scale_size(range = c(0,12)) + 
    scale_colour_manual(values=c("#c27d92", "#a8a34b", "#8e61bf")) + 
    scale_y_continuous(limits=c(3, 10))



###########

## Get DE genes that are RT-biased:
## Here they are: FBgn0200179, FBgn0197453, 
at3hrs=unique(subset(TPMse_DnovPM, FBgn_ID %in% unlist(RT_UP_3hrs_elements) & FBgn_ID %in% Dnov.dvir1.06.RT.list)$FBgn_ID)
at6hrs=unique(subset(TPMse_DnovPM, FBgn_ID %in% unlist(RT_UP_6hrs_elements) & FBgn_ID %in% Dnov.dvir1.06.RT.list)$FBgn_ID)
at12hrs=unique(subset(TPMse_DnovPM, FBgn_ID %in% unlist(RT_UP_12hrs_elements) & FBgn_ID %in% Dnov.dvir1.06.RT.list)$FBgn_ID)
TheGenes = union(at3hrs, union(at6hrs, at12hrs))

pdf(file = "RT-biased-sig.ofRT.pdf",width = 9, height = 3)
lapply(TheGenes, plotBoth)
dev.off()

