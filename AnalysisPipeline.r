
req_packages = c("Biobase", "cluster", "cowplot", "cummeRbund", 
                 "data.table", "DESeq", "edgeR", "ggplot2", 
                 "ggrepel", "ggthemes", "ggthemr", "Glimma", 
                 "GO.db", "goseq", "gplots", "grid", "gridExtra", 
                 "heatmap3", "imager", "plotly", "qvalue", "RColorBrewer", 
                 "reshape", "Rmisc", "splitstackshape", "statmod", 
                 "VennDiagram")

lapply(req_packages, require, character.only = TRUE)

# The Cowplot package changes the default themes of ggplot2. Set the default theme like so:
theme_set(theme_gray())

source("Functions2.R")

grpTrinotate = read.csv("Annotations/Trinotate_report_dvir1.06_subset.txt", header = T, sep = "\t", na.strings = "", stringsAsFactors=FALSE)
GO_info = read.table("Annotations/Trinotate_report_dvir1.06_gene_ontology.txt", header=F, row.names=1,stringsAsFactors=F)

gffRecord = read.table("Annotations/FBgn_ID_name_coordinates.txt", header = T)

melOrths = read.table(file = "Annotations/mel_orths.txt", header = T)
melOrthsAll = aggregate(mel_GeneSymbol~FBgn_ID, data = melOrths, toString)

Annots = merge(merge(melOrthsAll, grpTrinotate, all=TRUE), gffRecord, all=TRUE)

tmp.FB.names = unique(subset(Annots, select=c("FBgn_ID", "FBtr_ID")))
paml.data = read.csv(file = "Annotations/PAML.branchSite.ALL.results.txt", header = T, sep = "\t")
paml.data = merge(tmp.FB.names, paml.data, all=T)
paml.data = merge(gffRecord, paml.data, all=T)
KaKs.data = read.csv(file = "Annotations/KaKs.ALL.results.txt", header = T, sep = "\t", check.names = F)
KaKs.data = merge(tmp.FB.names, KaKs.data, all=T)
KaKs.data = merge(gffRecord, KaKs.data, all=T)

annot.sum = unique(subset(Annots, select=c("FBgn_ID", "gene_name", "SwissProt_BlastX_Description")))
annot.sum = aggregate(SwissProt_BlastX_Description ~ FBgn_ID + gene_name, data = annot.sum, toString)

countsMatrix = read.table("ExpressionData/genes_DnovPM_dvi1.06.counts.matrix", header=T, row.names=1, com='', check.names=F)

tpmMatrix = read.table("ExpressionData/genes_DnovPM_dvi1.06.TMM.EXPR.matrix", header=T, row.names=1, com='', check.names=F)

cbmtMatrix = read.table("ExpressionData/genes_DnovPM_dvi1.06.TPM.not_cross_norm.counts_by_min_TPM", header = T)

sampleData = read.table("ExpressionData/samples.txt", header=F, check.names=F, fill=T)
sampleData = sampleData[sampleData[,2] != '',]
sampleInfo = read.table("ExpressionData/SampleInfo.txt", header=T, check.names=F, fill=T)

## TPM average plots
tmp.tpmMatrix<-tpmMatrix
colnames(tmp.tpmMatrix) <- sampleData$V1
tmp.tpmMatrix.m <- as.data.frame(melt(as.matrix(tmp.tpmMatrix)))
tmp.tpmMatrix.m <- within(tmp.tpmMatrix.m, X2<-data.frame(do.call('rbind', strsplit(as.character(X2),'_',fixed=TRUE))))
tmp.tpmMatrix.m <- data.frame(tmp.tpmMatrix.m$X1, tmp.tpmMatrix.m$X2$X1, tmp.tpmMatrix.m$X2$X2, tmp.tpmMatrix.m$value)
colnames(tmp.tpmMatrix.m) <- c("FBgn_ID", "sample", "tissue", "TPM")
tmp.tpmMatrix.m$condition <- ifelse(grepl("C", tmp.tpmMatrix.m$sample, ignore.case = F), "conspecific", ifelse(grepl("H", tmp.tpmMatrix.m$sample, ignore.case = F), "heterospecific", "virgin"))
tmp.tpmMatrix.m$time <- ifelse(grepl("3", tmp.tpmMatrix.m$sample), "3hpm", ifelse(grepl("6", tmp.tpmMatrix.m$sample), "6hpm", ifelse(grepl("12", tmp.tpmMatrix.m$sample), "12hpm","virgin")))
tmp.tpmMatrix.m.c = summarySE(tmp.tpmMatrix.m, measurevar = "TPM", groupvars = c("FBgn_ID", "sample", "tissue", "condition", "time"))
fbgn_to_geneName<-subset(gffRecord, select=c(FBgn_ID, gene_name))
TPMse <- merge(fbgn_to_geneName, tmp.tpmMatrix.m.c, all=TRUE)
tmpMat<-cast(tmp.tpmMatrix.m.c, FBgn_ID~sample+tissue, value ="TPM")
meanTPMmatrix <- tmpMat[,-1]
rownames(meanTPMmatrix) <- tmpMat[,1]
TPMse$condition = factor (TPMse$condition, levels = c("virgin", "conspecific", "heterospecific"))
TPMse$time = factor (TPMse$time, levels = c("virgin", "3hpm", "6hpm", "12hpm"))


libSizes <- as.data.frame(colSums(countsMatrix))
libSizes <- cbind(sample = row.names(libSizes), libSizes)
row.names(libSizes)<- NULL
colnames(libSizes) = c("sample", "Total_reads")
options(repr.plot.width = 6, repr.plot.height = 3)
ggplot(libSizes, aes(sample, Total_reads)) + 
    geom_bar(stat="identity") + 
    theme(axis.text.x = element_text(angle = -90, hjust = 0)) + 
    geom_hline(yintercept = 20000000)

m.expData<-melt(as.matrix(tpmMatrix))
colnames(m.expData) <- c("gene_id", "replicate", "TPM")
m.expData.exp<- within(m.expData, replicate<-data.frame(do.call('rbind', strsplit(as.character(replicate),'_',fixed=TRUE))))
m.expData<-data.frame(m.expData, m.expData.exp$replicate$X1, m.expData.exp$replicate$X2, m.expData.exp$replicate$X3)
colnames(m.expData) <- c("gene_id", "replicate", "TPM", "sample", "tissue", "rep_num")
m.expData$TPM <- m.expData$TPM + 1
options(repr.plot.width = 9, repr.plot.height = 5)
ggplot(m.expData) + 
    geom_boxplot(aes(x = replicate, y = log10(TPM), fill = sample), size = 0.3, alpha = I(1/3)) + 
    facet_wrap(~tissue, scales = "free_x") +
    theme(axis.text.x = element_text(angle = -90, hjust = 0)) + 
    scale_fill_hue(l = 50, h.start = 200)

cbmt.sub = cbmtMatrix[cbmtMatrix[,1] > -100 & cbmtMatrix[,1] < -10,]

cbmt.sub_fit = lm(cbmt.sub[,2] ~ cbmt.sub[,1])
print(cbmt.sub_fit)

options(repr.plot.width = 4, repr.plot.height = 3)
ggplot(cbmtMatrix, aes(neg_min_tpm,num_features)) + 
    geom_point() +  
    scale_x_continuous(limits=c(-100,10)) + 
    scale_y_continuous(limits=c(0,20000)) + 
    geom_smooth(data=cbmt.sub, method = "lm") + 
    geom_hline(yintercept = 9185, colour = "green") + ggtitle("expressed genes", subtitle = "something")

all.CPM <- cpm(countsMatrix)

thresh <- all.CPM > 1

keep <- rowSums(thresh) >= 2

counts.keep <- countsMatrix[keep,]
dim(counts.keep)

countsMatrix.virgin = subset(countsMatrix, select=grepl("^V", colnames(countsMatrix)))
head(countsMatrix.virgin)

cpmMatrix.virgin <- cpm(countsMatrix.virgin)
v.thresh <- cpmMatrix.virgin > 5
## look at number of columns with above requirements
table(rowSums(v.thresh))

v.keep <- rowSums(v.thresh) >= 3
countsMatrix.virgin.filt <- countsMatrix.virgin[v.keep,]
## Check how many genes remain in the matrix.
dim(countsMatrix.virgin.filt)

# Let's look at the first column
plot(cpmMatrix.virgin[,1],countsMatrix.virgin[,1],ylim=c(0,100),xlim=c(0,5))
# Add a vertical line at 5 CPM
abline(v=5)

sampleInfo.v = subset(sampleInfo, Status == "virgin")
sampleInfo.v

groups.v = factor(sampleInfo.v$Tissue)
design.v = model.matrix( ~ 0 + groups.v)
colnames(design.v) <- levels(groups.v)
rownames(design.v) <- sampleInfo.v$Replicate
design.v

dgeList.v <- DGEList(counts = countsMatrix.virgin.filt, group = groups.v)
dgeList.v <- calcNormFactors(dgeList.v)
dgeList.v <- estimateCommonDisp(dgeList.v)
dgeList.v <- estimateTagwiseDisp(dgeList.v)
dgeList.v_fit <- glmFit(dgeList.v, design.v)

# Extract annotation for genes in the fit object
ann.v = subset(annot.sum, FBgn_ID %in% rownames(dgeList.v_fit))
# convert factors to characters
ann.v = data.frame(lapply(ann.v, as.character), stringsAsFactors=FALSE)
# align the fit object's rownames with gene ID's
ann.v = ann.v[match(rownames(dgeList.v_fit), ann.v$FBgn_ID),]
# convert factors to characters, again
ann.v <- data.frame(lapply(ann.v, as.character), stringsAsFactors=FALSE)
# Rename "FBgn_ID" as "GeneID"
colnames(ann.v) = c ("GeneID", "gene_name", "SwissProt_BlastX_Description")
# Check that the fit rownames match the annotation file's gene ID's
table(ann.v$GeneID==rownames(dgeList.v_fit))
# Add the annotations to the fit object in the "genes" slot
dgeList.v_fit$genes = ann.v

dgeList.v_fit

summary(dgeList.v$tagwise.dispersion)

options(repr.plot.width = 9, repr.plot.height = 6)
par(mfrow=c(2,2))
# Biological coefficient of variation
plotBCV(dgeList.v)
# mean-variance trend
virgin.voom = voom(dgeList.v, design.v, plot=TRUE)
# QQ-plot
g.v <- gof(dgeList.v_fit)
z.v <- zscoreGamma(g.v$gof.statistics,shape=g.v$df/2,scale=2)
qqnorm(z.v); qqline(z.v, col = 4,lwd=1,lty=1)
# log2 transformed and normalize boxplot of counts across samples
boxplot(virgin.voom$E, xlab="", ylab="Log2 counts per million",las=2,main="Voom transformed logCPM")
abline(h=median(virgin.voom$E),col="blue")

## colour samples by tissue-type
col.tissue <- c("#8d75ca","#78a450","#c16786","#c7733b")[sampleInfo.v$Tissue]
options(repr.plot.width = 7, repr.plot.height = 4)
plotMDS(dgeList.v, col=col.tissue, pch= 17, cex = 1)
legend("top",fill=c("#8d75ca","#78a450","#c16786","#c7733b"),legend=levels(sampleInfo.v$Tissue), cex = 0.85)
# Add a title
title("Virgin tissue")

glMDSPlot(dgeList.v, groups = dgeList.v$samples$group, labels = sampleInfo.v$Replicate)

## Plot sample correlation
data = log2(countsMatrix.virgin.filt+1)
data = as.matrix(data)
sample_cor = cor(data, method='pearson', use='pairwise.complete.obs')
sample_dist = as.dist(1-sample_cor)
hc_samples = hclust(sample_dist, method='complete')

options(repr.plot.width = 4, repr.plot.height = 4)
heatmap.3(sample_cor, dendrogram='both', Rowv=as.dendrogram(hc_samples), Colv=as.dendrogram(hc_samples), col = colorpanel(75, '#dd70cd','black','#afc64f'), scale='none', symm=TRUE, key=TRUE,density.info='none', trace='none', symkey=FALSE, symbreaks=F, cexCol=1, cexRow=1, cex.main=0.75, main=paste("sample correlation matrix"))

# normalize counts by the DESeq method (this is only used for these plots:
meta.v <- data.frame(row.names=colnames(countsMatrix.virgin.filt), condition=sampleInfo.v$Tissue)
countData.v<-round(countsMatrix.virgin.filt)
countData.v_normByDESeq = newCountDataSet(countData.v, meta.v)
countData.v_normByDESeq = estimateSizeFactors(countData.v_normByDESeq)
countData.v_normByDESeq = data.frame(counts(countData.v_normByDESeq, normalized=T))

options(repr.plot.width = 9, repr.plot.height = 5)
MA_BPlot(countData.v_normByDESeq, "V_CR_1", "V_CR_2")
MA_BPlot(countData.v_normByDESeq, "V_H_1", "V_H_2")
MA_BPlot(countData.v_normByDESeq, "V_OV_1", "V_OV_2")
MA_BPlot(countData.v_normByDESeq, "V_RT_1", "V_RT_2")
MA_BPlot(countData.v_normByDESeq, "V_RT_1", "V_RT_3")
MA_BPlot(countData.v_normByDESeq, "V_RT_2", "V_RT_3")

cont.v.repTract <- makeContrasts(V_RT.vs.V_CR=repTract-carcass,
                                 V_RT.vs.V_HD=repTract-head,
                                 V_RT.vs.V_OV=repTract-ovaries,
                                 levels=design.v)
cont.v.ovaries <- makeContrasts(V_OV.vs.V_CR=ovaries-carcass,
                                V_OV.vs.V_HD=ovaries-head,
                                V_OV.vs.V_RT=ovaries-repTract,
                                levels=design.v)
cont.v.head <- makeContrasts(V_H.vs.V_CR=head-carcass,
                             V_H.vs.V_OV=head-ovaries,
                             V_H.vs.V_RT=head-repTract,
                             levels=design.v)

lrt.v.repTract <- glmLRT(dgeList.v_fit, contrast = cont.v.repTract)
lrt.v.repTract.tTags <- topTags(lrt.v.repTract, n = NULL)
lrt.v.repTract.tTags.table <- lrt.v.repTract.tTags$table
repTract.list<-subset(lrt.v.repTract.tTags.table, logFC.V_RT.vs.V_CR > 2 & logFC.V_RT.vs.V_HD > 2 & logFC.V_RT.vs.V_OV > 2 & FDR<0.001)$GeneID
length(repTract.list)
head(lrt.v.repTract.tTags.table)

lrt.v.ovaries <- glmLRT(dgeList.v_fit, contrast = cont.v.ovaries)
lrt.v.ovaries.tTags <- topTags(lrt.v.ovaries, n = NULL)
lrt.v.ovaries.tTags.table <- lrt.v.ovaries.tTags$table
ovaries.list<-subset(lrt.v.ovaries.tTags.table, logFC.V_OV.vs.V_CR > 2 & logFC.V_OV.vs.V_HD > 2 & logFC.V_OV.vs.V_RT > 2 & FDR<0.001)$GeneID
length(ovaries.list)
head(lrt.v.ovaries.tTags.table)

lrt.v.head <- glmLRT(dgeList.v_fit, contrast = cont.v.head)
lrt.v.head.tTags <- topTags(lrt.v.head, n = NULL)
lrt.v.head.tTags.table <- lrt.v.head.tTags$table
head.list<-subset(lrt.v.head.tTags.table, logFC.V_H.vs.V_CR > 2 & logFC.V_H.vs.V_OV > 2 & logFC.V_H.vs.V_RT > 2 & FDR<0.001)$GeneID
length(head.list)
head(lrt.v.head.tTags.table)

############## Proceed from here:

## Generate a factor labeling matrix to annotate tissue-biased genes
RT_factors = as.data.frame(repTract.list)
RT_factors$V1 = "RT-biased"
rownames(RT_factors) = repTract.list
RT_factors = subset(RT_factors, select = "V1")

OV_factors = as.data.frame(ovaries.list)
OV_factors$V1 = "OV-biased"
rownames(OV_factors) = ovaries.list
OV_factors = subset(OV_factors, select = "V1")

H_factors = as.data.frame(head.list)
H_factors$V1 = "H-biased"
rownames(H_factors) = head.list
H_factors = subset(H_factors, select = "V1")

virgin.factor.labeling = rbind(RT_factors, OV_factors, H_factors)
colnames(virgin.factor.labeling) = c('tissue_bias')
virgin.factor_list = unique(virgin.factor.labeling[,1])

## For a heatmap of tissue-biased genes:
tissueBiased.meanTPM = subset(meanTPMmatrix, rownames(meanTPMmatrix) %in% rownames(virgin.factor.labeling))
tissueBiased.meanTPM.v = subset(tissueBiased.meanTPM, select=grepl("^V", colnames(tissueBiased.meanTPM)))
data = tissueBiased.meanTPM.v
gene_factors = unique(virgin.factor.labeling[,1])
gene_factor_colors = rainbow(length(gene_factors))
names(gene_factor_colors) = gene_factors
data = log2(data+1)
data = as.matrix(data) # convert to matrix
sample_cor = cor(data, method='pearson', use='pairwise.complete.obs')
sample_dist = as.dist(1-sample_cor)
hc_samples = hclust(sample_dist, method='complete')
gene_cor = NULL
if (is.null(gene_cor)) { gene_cor = cor(t(data), method='pearson', use='pairwise.complete.obs') }
gene_dist = as.dist(1-gene_cor)
hc_genes = hclust(gene_dist, method='complete')
myheatcol = colorpanel(75, '#f8e95a','black','#e1526d')
data = t(scale(t(data), scale=F)) # center rows, mean substracted
heatmap_data = data
gene_factor_row_vals = as.factor(virgin.factor.labeling[rownames(heatmap_data),])
names(gene_factor_row_vals) = rownames(heatmap_data)
gene_factors_here = unique(gene_factor_row_vals)
names(gene_factors_here) = gene_factors_here
num_gene_factors_here = length(gene_factors_here)
geneFactorColors = c("#e1526d", "#b85516", "#647700")
geneFactorAnnotations = matrix(nrow=nrow(heatmap_data), ncol=num_gene_factors_here)
for (i in 1:num_gene_factors_here) {
    geneFactorAnnotations[,i] = gene_factor_row_vals %in% gene_factors_here[i]
}
geneFactorAnnotations = apply(geneFactorAnnotations, 1:2, function(x) as.logical(x))
geneFactorAnnotations = t(sample_matrix_to_color_assignments(t(geneFactorAnnotations), col=geneFactorColors))
rownames(geneFactorAnnotations) = rownames(heatmap_data)
colnames(geneFactorAnnotations) = gene_factors_here
heatmap_data[heatmap_data < -2] = -2
heatmap_data[heatmap_data > 2] = 2

## plot it
heatmap3(heatmap_data, 
         col=myheatcol, 
         cexCol = 2, 
         labRow = "", 
         legendfun=function() showLegend(legend=c("OV-biased", "H-biased", "RT-biased"), col=c("#e1526d", "#b85516", "#647700")), 
         cex.main=0.75, 
         RowSideColors = geneFactorAnnotations)

heatmap3(x = tissueBiased.meanTPM.v, balanceColor = T, ColSideLabs = "something", RowSideColors = geneFactorAnnotations)


## Calculate tissue specificity for all genes:
virginSpecificity = calcSpecificity(tissueBiased.meanTPM.v)
head(virginSpecificity)
##############################  GAP HERE ##############################

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
options(repr.plot.width = 7, repr.plot.height = 3)

plotGenePM(TPMse_DnovPM, "GJ22636")

subset(paml.data, gene_name == "GJ19434")

### this is incomplete
Dnov_virgin_tissue_MeanTPMmatrix <- subset(DnovPM_MeanTPMmatrix, select=c(FBgn_ID, V_CR, V_H, V_OV, V_RT))
rownames(Dnov_virgin_tissue_MeanTPMmatrix) <- Dnov_virgin_tissue_MeanTPMmatrix[,1]
Dnov_virgin_tissue_MeanTPMmatrix[,1] <- NULL

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

## RT-biased genes
lrt.RT.v.rest <- glmLRT(DnovPM_fit, contrast = virgin_RT_contrasts)
lrt.RT.v.rest.tTags <- topTags(lrt.RT.v.rest, n = NULL)
lrt.RT.v.rest.tTags.table <- lrt.RT.v.rest.tTags$table
Dnov.dvir1.06.RT.list<-rownames(subset(lrt.RT.v.rest.tTags.table, logFC.V_RT.vs.V_CR > 2 & logFC.V_RT.vs.V_H > 2 & logFC.V_RT.vs.V_OV > 2 & FDR<0.001))

## OV-biased genes
lrt.OV.v.rest <- glmLRT(DnovPM_fit, contrast = virgin_OV_contrasts)
lrt.OV.v.rest.tTags <- topTags(lrt.OV.v.rest, n = NULL)
lrt.OV.v.rest.tTags.table <- lrt.OV.v.rest.tTags$table
Dnov.dvir1.06.OV.list<-rownames(subset(lrt.OV.v.rest.tTags.table, logFC.V_OV.vs.V_CR > 2 & logFC.V_OV.vs.V_H > 2 & logFC.V_OV.vs.V_RT > 2 & FDR<0.001))

## H-biased genes
lrt.H.v.rest <- glmLRT(DnovPM_fit, contrast = virgin_H_contrasts)
lrt.H.v.rest.tTags <- topTags(lrt.H.v.rest, n = NULL)
lrt.H.v.rest.tTags.table <- lrt.H.v.rest.tTags$table
Dnov.dvir1.06.H.list<-rownames(subset(lrt.H.v.rest.tTags.table, logFC.V_H.vs.V_CR > 2 & logFC.V_H.vs.V_OV > 2 & logFC.V_H.vs.V_RT > 2 & FDR<0.001))

Dnov.dvir1.06.RT.matrix <- subset(Dnov_virgin_tissue_MeanTPMmatrix, rownames(Dnov_virgin_tissue_MeanTPMmatrix) %in% Dnov.dvir1.06.RT.list)

length(Dnov.dvir1.06.RT.list) 

options(repr.plot.width = 4, repr.plot.height = 4)
plotHeatmap(Dnov.dvir1.06.RT.matrix, clustering = "both", labRow = F)

options(repr.plot.width = 9, repr.plot.height = 2)
ggplot(subset(paml.data, FBgn_ID %in% Dnov.dvir1.06.RT.list & omega < 800 & grepl("Chr", chromosome)), aes(max, omega)) + 
    geom_point(size=2, alpha=0.5, colour = "#7d49c3") + 
  #  geom_point(data=subset(paml.data, FBgn_ID %in% SFP_elements$`D.ame,D.lum,D.nov,D.vir` ), aes(max, omega), inherit.aes = F, size=2, alpha=0.5, colour = "#4f922a") + 
    geom_hline(yintercept = 0.15, linetype="dashed", colour = "yellow") + 
    geom_hline(yintercept = 1, linetype="dashed", colour = "gray")  + 
    facet_grid(~chromosome, scales = "free_x") + 
#    scale_colour_manual(name = "", values =c("#7aa457"="#7aa457","#9e6ebd"="#9e6ebd"), labels = c("SFPs","EB biased")) + 
    scale_x_continuous(breaks=seq(5000000,30000000,5000000), labels=expression("5", "10", "15", "20", "25", "30")) + 
    xlab ("Chromosome coordinates (Mb)") + 
    labs(y=expression(K[a]/K[s])) + 
#    geom_text_repel(data=subset(paml.data, FBgn_ID %in% Dnov.dvir1.06.RT.list & omega > 0.95 & omega < 800), aes(label = gene_name), size =3, force = 30, colour = "#7d49c3") +
  #  geom_text_repel(data=subset(paml.data, gene_name == "GJ21515"), aes(label = gene_name), size =3, force = 30, colour = "#4f922a") +
  #  geom_text_repel(data=subset(paml.data, FBgn_ID %in% SFP_elements$`D.ame,D.lum,D.nov,D.vir` & omega > 0.8), aes(label = gene_name), size =3, force = 4, colour = "#4f922a") +
    theme(axis.title.x = element_text(face = "bold", size = 10, vjust=0.1), axis.text.x=element_text(face = "bold", size = 12),axis.text.y = element_text(face = "bold", size = 12), axis.title.y = element_text(face = "bold.italic", size = 12, vjust=0.1), strip.text=element_text(face="bold", size = 12))



gene_lengths = read.table("GO.analysis/FBgn_lengths.txt", header=T, row.names=1)
gene_lengths = as.matrix(gene_lengths[,1,drop=F])
GO_info = read.table("GO.analysis/Trinotate_report_dvir1.06_gene_ontology.txt", header=F, row.names=1,stringsAsFactors=F)
GO_info_listed = apply(GO_info, 1, function(x) unlist(strsplit(x,',')))
names(GO_info_listed) = rownames(GO_info)
features_with_GO = rownames(GO_info)
lengths_features_with_GO = gene_lengths[features_with_GO,]
get_GO_term_descr =  function(x) {
    d = 'none';
    go_info = GOTERM[[x]];
    if (length(go_info) >0) { d = paste(Ontology(go_info), Term(go_info), sep=' ');}
    return(d);
}

RT_factors = as.data.frame(Dnov.dvir1.06.RT.list)
RT_factors$V1 = "RT-biased"
rownames(RT_factors) = Dnov.dvir1.06.RT.list
RT_factors = subset(RT_factors, select = "V1")

OV_factors = as.data.frame(Dnov.dvir1.06.OV.list)
OV_factors$V1 = "OV-biased"
rownames(OV_factors) = Dnov.dvir1.06.OV.list
OV_factors = subset(OV_factors, select = "V1")

H_factors = as.data.frame(Dnov.dvir1.06.H.list)
H_factors$V1 = "H-biased"
rownames(H_factors) = Dnov.dvir1.06.H.list
H_factors = subset(H_factors, select = "V1")

factor_labeling = rbind(RT_factors, OV_factors, H_factors)
colnames(factor_labeling) = c('type')
factor_list = unique(factor_labeling[,1])

# build pwf based on ALL DE features
options(repr.plot.width = 6, repr.plot.height = 3)

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

### RT plot
options(repr.plot.width = 7, repr.plot.height = 4)
ggplot(subset(GO_enrichment_data, over_represented_FDR < 0.05 & factor == "RT-biased"), 
       aes(category, -log10(over_represented_pvalue), size = numDEInCat, colour = ontology)) + 
    geom_point()  + 
    xlab(NULL) + 
    geom_text_repel(data = subset(GO_enrichment_data, factor == "RT-biased" & over_represented_FDR < 0.05 & numDEInCat > 20), 
                    aes(category, -log10(over_represented_pvalue),label=term), 
                    force = 8, 
                    inherit.aes = F, 
                    box.padding = unit(0.35, "lines"), 
                    point.padding = unit(0.5, "lines"), 
                    fontface = "bold", 
                    size = 4) + 
    theme(axis.text.x = element_text(angle = 45, face = "bold", vjust = 1, hjust = 1)) + 
    scale_size(range = c(0,12)) + 
    scale_colour_manual(values=c("#e00068",
                                 "#01777b",
                                 "#9d31e4")) + 
    scale_y_continuous(limits=c(3, 9))

DnovPM_CountsMatrix_RT = subset(DnovPM.dvir1.06.CountsMatrix.BRR.min400count, select=grepl("RT", colnames(DnovPM.dvir1.06.CountsMatrix.BRR.min400count)))

DnovPM_max_gene_expr_per_row = apply(DnovPM_CountsMatrix_RT, 1, max)
DnovPM_CountsMatrix_RT.min400 = DnovPM_CountsMatrix_RT[DnovPM_max_gene_expr_per_row >= 400,,drop=F ]

head(DnovPM_CountsMatrix_RT.min400)

RT.group <- factor(c(1,1,1,2,2,2,3,3,4,4,4,5,5,6,6,7,7,7))
RT.design <- model.matrix(~0+RT.group)
colnames(RT.design)<-c("C12_RT", "C3_RT", "C6_RT", "H12_RT", "H3_RT", "H6_RT", "V_RT")

DnovPM_DGElist_RT<-DGEList(counts = DnovPM_CountsMatrix_RT.min400, group = RT.group)
DnovPM_DGElist_RT<-calcNormFactors(DnovPM_DGElist_RT)
DnovPM_DGElist_RT<-estimateDisp(DnovPM_DGElist_RT, RT.design, robust = T)
DnovPM_RT_fit <- glmFit(DnovPM_DGElist_RT, RT.design)

g <- gof(DnovPM_RT_fit)
z <- zscoreGamma(g$gof.statistics,shape=gof$df / 2,scale=2)

con_virgin_contrasts <- makeContrasts(C3.vs.V=C3_RT-V_RT, C6.vs.V=C6_RT-V_RT, C12.vs.V=C12_RT-V_RT, levels = RT.design)
het_virgin_contrasts <- makeContrasts(H3.vs.V=H3_RT-V_RT, H6.vs.V=H6_RT-V_RT, H12.vs.V=H12_RT-V_RT, levels = RT.design)

RT_con.vs.virgin <- glmLRT(DnovPM_RT_fit, contrast = con_virgin_contrasts)
RT_con.vs.virgin.tTags <- topTags(RT_con.vs.virgin, n = NULL)
RT_con.vs.virgin.tTags.table <- RT_con.vs.virgin.tTags$table
RT_con.vs.virgin.tTags.table$FBgn_ID <- rownames(RT_con.vs.virgin.tTags.table)

RT_het.vs.virgin <- glmLRT(DnovPM_RT_fit, contrast = het_virgin_contrasts)
RT_het.vs.virgin.tTags <- topTags(RT_het.vs.virgin, n = NULL)
RT_het.vs.virgin.tTags.table <- RT_het.vs.virgin.tTags$table
RT_het.vs.virgin.tTags.table$FBgn_ID <- rownames(RT_het.vs.virgin.tTags.table)

colnames(RT_con.vs.virgin.tTags.table) = c("3hpm", "6hpm", "12hpm", "logCPM", "LR", "pvalue", "FDR", "FBgn_ID")
colnames(RT_het.vs.virgin.tTags.table) = c("3hpm", "6hpm", "12hpm", "logCPM", "LR", "pvalue", "FDR", "FBgn_ID")

head(RT_con.vs.virgin.tTags.table)
head(RT_het.vs.virgin.tTags.table)

PM.RT.con_vs_virgin.m <-melt(RT_con.vs.virgin.tTags.table, id.vars = c("pvalue", "FDR", "FBgn_ID", "LR", "logCPM"))
PM.RT.con_vs_virgin.m$cross = "Conspecific"
PM.RT.het_vs_virgin.m <-melt(RT_het.vs.virgin.tTags.table, id.vars = c("pvalue", "FDR", "FBgn_ID", "LR", "logCPM"))
PM.RT.het_vs_virgin.m$cross = "Heterospecific"

RT.PM_vs_virgin.m = rbind(PM.RT.con_vs_virgin.m, PM.RT.het_vs_virgin.m)
RT.PM_vs_virgin.m$sig = ifelse(RT.PM_vs_virgin.m$FDR < 0.00001, "YES", "NO")
head(RT.PM_vs_virgin.m)

y <- RT.PM_vs_virgin.m$pvalue
qqnorm(y); qqline(y, col = 2)
qqplot(y, rt(300, df = 5))

ggplot(data, aes(sample = RT.PM_vs_virgin.m$pvalue)) + stat_qq(color="firebrick2", alpha=1) + geom_abline(intercept = mean(RT.PM_vs_virgin.m$pvalue), slope = sd(RT.PM_vs_virgin.m$pvalue))

options(repr.plot.width = 9, repr.plot.height = 4)
ggplot(RT.PM_vs_virgin.m, aes(value, -log10(FDR), colour = sig)) + 
    geom_point(alpha=I(1/2)) + 
    facet_grid(cross~variable, scales = "free")

length(unique(subset(RT.PM_vs_virgin.m, sig == "YES")$FBgn_ID))

# C3_RT vs V_RT
RT_con.3hrs.vs.virgin <- glmLRT(DnovPM_RT_fit, contrast = con_virgin_contrasts[,"C3.vs.V"])
RT_con.3hrs.vs.virgin.tTags <- topTags(RT_con.3hrs.vs.virgin, n = NULL)
RT_con.3hrs.vs.virgin.tTags.table <- RT_con.3hrs.vs.virgin.tTags$table
RT_con.3hrs.vs.virgin.Up.list <- rownames(subset(RT_con.3hrs.vs.virgin.tTags.table, logFC > 1 & FDR < 0.001))
RT_con.3hrs.vs.virgin.Down.list <- rownames(subset(RT_con.3hrs.vs.virgin.tTags.table, logFC < -1 & FDR < 0.001))
# C6_RT vs V_RT
RT_con.6hrs.vs.virgin <- glmLRT(DnovPM_RT_fit, contrast = con_virgin_contrasts[,"C6.vs.V"])
RT_con.6hrs.vs.virgin.tTags <- topTags(RT_con.6hrs.vs.virgin, n = NULL)
RT_con.6hrs.vs.virgin.tTags.table <- RT_con.6hrs.vs.virgin.tTags$table
RT_con.6hrs.vs.virgin.Up.list <- rownames(subset(RT_con.6hrs.vs.virgin.tTags.table, logFC > 1 & FDR < 0.001))
RT_con.6hrs.vs.virgin.Down.list <- rownames(subset(RT_con.6hrs.vs.virgin.tTags.table, logFC < -1 & FDR < 0.001))
# C12_RT vs V_RT
RT_con.12hrs.vs.virgin <- glmLRT(DnovPM_RT_fit, contrast = con_virgin_contrasts[,"C12.vs.V"])
RT_con.12hrs.vs.virgin.tTags <- topTags(RT_con.12hrs.vs.virgin, n = NULL)
RT_con.12hrs.vs.virgin.tTags.table <- RT_con.12hrs.vs.virgin.tTags$table
RT_con.12hrs.vs.virgin.Up.list <- rownames(subset(RT_con.12hrs.vs.virgin.tTags.table, logFC > 1 & FDR < 0.001))
RT_con.12hrs.vs.virgin.Down.list <- rownames(subset(RT_con.12hrs.vs.virgin.tTags.table, logFC < -1 & FDR < 0.001))

# H3_RT vs V_RT
RT_het.3hrs.vs.virgin <- glmLRT(DnovPM_RT_fit, contrast = het_virgin_contrasts[,"H3.vs.V"])
RT_het.3hrs.vs.virgin.tTags <- topTags(RT_het.3hrs.vs.virgin, n = NULL)
RT_het.3hrs.vs.virgin.tTags.table <- RT_het.3hrs.vs.virgin.tTags$table
RT_het.3hrs.vs.virgin.Up.list <- rownames(subset(RT_het.3hrs.vs.virgin.tTags.table, logFC > 1 & FDR < 0.001))
RT_het.3hrs.vs.virgin.Down.list <- rownames(subset(RT_het.3hrs.vs.virgin.tTags.table, logFC < -1 & FDR < 0.001))
# H6_RT vs V_RT
RT_het.6hrs.vs.virgin <- glmLRT(DnovPM_RT_fit, contrast = het_virgin_contrasts[,"H6.vs.V"])
RT_het.6hrs.vs.virgin.tTags <- topTags(RT_het.6hrs.vs.virgin, n = NULL)
RT_het.6hrs.vs.virgin.tTags.table <- RT_het.6hrs.vs.virgin.tTags$table
RT_het.6hrs.vs.virgin.Up.list <- rownames(subset(RT_het.6hrs.vs.virgin.tTags.table, logFC > 1 & FDR < 0.001))
RT_het.6hrs.vs.virgin.Down.list <- rownames(subset(RT_het.6hrs.vs.virgin.tTags.table, logFC < -1 & FDR < 0.001))
# H12_RT vs V_RT
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
RT_UP_3hrs_Vdiag<-venn.diagram(RT_UP_3hrs_candidates, NULL, fill=c("#b067a3", "#9c954d"), alpha=c(0.75,0.75), cex = 1, cat.fontface= 2, cat.cex = 0, resolution = 1000, main = "3hpm")
RT_UP_6hrs_Vdiag<-venn.diagram(RT_UP_6hrs_candidates, NULL, fill=c("#b067a3", "#9c954d"), alpha=c(0.75,0.75), cex = 1, cat.fontface= 2, cat.cex = 0, resolution = 1000, main = "6hpm")
RT_UP_12hrs_Vdiag<-venn.diagram(RT_UP_12hrs_candidates, NULL, fill=c("#b067a3", "#9c954d"), alpha=c(0.75,0.75), cex = 1, cat.fontface= 2, cat.cex = 0, resolution = 1000, main = "12hpm")

options(repr.plot.width = 4, repr.plot.height = 3)
grid.arrange(gTree(children=RT_UP_3hrs_Vdiag), gTree(children=RT_UP_6hrs_Vdiag), gTree(children=RT_UP_12hrs_Vdiag))

RT_Down_3hrs_Vdiag<-venn.diagram(RT_Down_3hrs_candidates, NULL, fill=c("#b067a3", "#9c954d"), alpha=c(0.75,0.75), cex = 1.5, cat.fontface= 2, cat.cex = 0, resolution = 1000, main = "3hpm")
RT_Down_6hrs_Vdiag<-venn.diagram(RT_Down_6hrs_candidates, NULL, fill=c("#b067a3", "#9c954d"), alpha=c(0.75,0.75), cex = 1.5, cat.fontface= 2, cat.cex = 0, resolution = 1000, main = "6hpm")
RT_Down_12hrs_Vdiag<-venn.diagram(RT_Down_12hrs_candidates, NULL, fill=c("#b067a3", "#9c954d"), alpha=c(0.75,0.75), cex = 1.5, cat.fontface= 2, cat.cex = 0, resolution = 1000, main = "12hpm")

options(repr.plot.width = 4, repr.plot.height = 3)
grid.arrange(gTree(children=RT_Down_3hrs_Vdiag), gTree(children=RT_Down_6hrs_Vdiag), gTree(children=RT_Down_12hrs_Vdiag))

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


PM.vs.vir_Up_candidates_Vdiag<-venn.diagram(PM.vs.vir_Up_candidates, NULL, fill=c("#b067a3", "#9c954d"), alpha=c(0.75,0.75), cex = 1.5, cat.fontface= 4, cat.cex = 0, resolution = 1000, main = "Upregulated")
PM.vs.vir_Down_candidates_Vdiag<-venn.diagram(PM.vs.vir_Down_candidates, NULL, fill=c("#b067a3", "#9c954d"), alpha=c(0.75,0.75), cex = 1.5, cat.fontface= 4, cat.cex = 0, resolution = 1000, main = "Downregulated")
options(repr.plot.width = 4, repr.plot.height = 1)
grid.arrange(gTree(children=PM.vs.vir_Up_candidates_Vdiag))
grid.arrange(gTree(children=PM.vs.vir_Down_candidates_Vdiag))


RT_UP_3hrs_candidates

options(repr.plot.width = 6, repr.plot.height = 2)
lapply(RT_UP_6hrs_candidates$Heterospecific, plotGenePM, object = TPMse_DnovPM)

########### RT contrasts
# create RT-specific count matrix and define groups/design
DnovPM_CountsMatrix_OV = subset(DnovPM.dvir1.06.CountsMatrix.BRR.min400count, select=grepl("OV", colnames(DnovPM.dvir1.06.CountsMatrix.BRR.min400count)))

## filter this matrix by minimum read count
DnovPM_max_gene_expr_per_row = apply(DnovPM_CountsMatrix_OV, 1, max)
DnovPM_CountsMatrix_OV.min200 = DnovPM_CountsMatrix_OV[DnovPM_max_gene_expr_per_row >= 200,,drop=F ]
head(DnovPM_CountsMatrix_OV.min200)

######
OV.group <- factor(c(1,1,2,2,3,3))
OV.design <- model.matrix(~0+OV.group)
colnames(OV.design)<-c("C6_OV", "H6_OV", "V_OV")
OV.group
OV.design

# create edgeR DE object and run glm
DnovPM_DGElist_OV<-DGEList(counts = DnovPM_CountsMatrix_OV.min200, group = OV.group)
DnovPM_DGElist_OV<-calcNormFactors(DnovPM_DGElist_OV)
DnovPM_DGElist_OV<-estimateDisp(DnovPM_DGElist_OV, OV.design, robust = T)
DnovPM_OV_fit <- glmFit(DnovPM_DGElist_OV, OV.design)

# define contrasts
OV_contrasts <- makeContrasts(C6.vs.V=C6_OV-V_OV, H6.vs.V=H6_OV-V_OV, C6.vs.H6=C6_OV-H6_OV, levels = OV.design)
OV_contrasts

# identify overall con_vs_het
OV_all_comaprisons <- glmLRT(DnovPM_OV_fit, contrast = OV_contrasts)
OV_all_comaprisons.tTags <- topTags(OV_all_comaprisons, n = NULL)
OV_all_comaprisons.tTags.table <- OV_all_comaprisons.tTags$table
OV_all_comaprisons.tTags.table$FBgn_ID <- rownames(OV_all_comaprisons.tTags.table)
head(OV_all_comaprisons.tTags.table)

colnames(OV_all_comaprisons.tTags.table) = c("conspceific", "heterospecific", "logCPM", "LR", "p-value", "FDR", "FBgn_ID")
OV_all_comaprisons.tTags.table.m <-melt(OV_all_comaprisons.tTags.table, id.vars = c("p-value", "FDR", "FBgn_ID", "LR", "logCPM"))
OV_all_comaprisons.tTags.table.m$sig = ifelse(OV_all_comaprisons.tTags.table.m$FDR < 0.001 & OV_all_comaprisons.tTags.table.m$value > 1, "YES", "NO")
head(OV_all_comaprisons.tTags.table.m)

options(repr.plot.width = 6, repr.plot.height = 2)
ggplot(OV_all_comaprisons.tTags.table.m, aes(value, -log10(FDR), colour = sig)) + 
    geom_point(alpha=I(1/2)) + 
    facet_wrap(~variable, scales = "free") +
    geom_text_repel(data=subset(OV_all_comaprisons.tTags.table.m, -log10(FDR) > 10), aes(label = FBgn_ID), size =2.5, colour = "#4f922a")

OV_DE.list = unique(subset(OV_all_comaprisons.tTags.table.m, sig == "YES")$FBgn_ID)
lapply(OV_DE.list, plotGenePM, object = TPMse_DnovPM)

########### RT contrasts
# create RT-specific count matrix and define groups/design
DnovPM_CountsMatrix_H = subset(DnovPM.dvir1.06.CountsMatrix.BRR.min400count, select=grepl("_H", colnames(DnovPM.dvir1.06.CountsMatrix.BRR.min400count)))

## filter this matrix by minimum read count
DnovPM_max_gene_expr_per_row = apply(DnovPM_CountsMatrix_H, 1, max)
DnovPM_CountsMatrix_H.min200 = DnovPM_CountsMatrix_H[DnovPM_max_gene_expr_per_row >= 200,,drop=F ]
head(DnovPM_CountsMatrix_H.min200)

######
H.group <- factor(c(1,1,2,2,3,3))
H.design <- model.matrix(~0+H.group)
colnames(H.design)<-c("C6_H", "H6_H", "V_H")
H.design

# create edgeR DE object and run glm
DnovPM_DGElist_H<-DGEList(counts = DnovPM_CountsMatrix_H.min200, group = H.group)
DnovPM_DGElist_H<-calcNormFactors(DnovPM_DGElist_H)
DnovPM_DGElist_H<-estimateDisp(DnovPM_DGElist_H, H.design, robust = T)
DnovPM_H_fit <- glmFit(DnovPM_DGElist_H, H.design)

# define contrasts
H_contrasts <- makeContrasts(C6.vs.V=C6_H-V_H, H6.vs.V=H6_H-V_H, C6.vs.H6=C6_H-H6_H, levels = H.design)
H_contrasts

# identify overall con_vs_het
H_all_comaprisons <- glmLRT(DnovPM_H_fit, contrast = H_contrasts)
H_all_comaprisons.tTags <- topTags(H_all_comaprisons, n = NULL)
H_all_comaprisons.tTags.table <- H_all_comaprisons.tTags$table
H_all_comaprisons.tTags.table$FBgn_ID <- rownames(H_all_comaprisons.tTags.table)
head(H_all_comaprisons.tTags.table)

colnames(H_all_comaprisons.tTags.table) = c("conspceific", "heterospecific", "logCPM", "LR", "p-value", "FDR", "FBgn_ID")
H_all_comaprisons.tTags.table.m <-melt(H_all_comaprisons.tTags.table, id.vars = c("p-value", "FDR", "FBgn_ID", "LR", "logCPM"))
H_all_comaprisons.tTags.table.m$sig = ifelse(H_all_comaprisons.tTags.table.m$FDR < 0.001, "YES", "NO")
head(H_all_comaprisons.tTags.table.m)

options(repr.plot.width = 6, repr.plot.height = 2)
ggplot(H_all_comaprisons.tTags.table.m, aes(value, -log10(FDR), colour = sig)) + 
    geom_point(alpha=I(1/2)) + 
    facet_wrap(~variable, scales = "free")
#    geom_text_repel(data=subset(H_all_comaprisons.tTags.table.m, -log10(FDR) > 10), aes(label = FBgn_ID), size =2.5, colour = "#4f922a")

H_DE.list = unique(subset(H_all_comaprisons.tTags.table.m, sig == "YES")$FBgn_ID)
lapply(H_DE.list, plotGenePM, object = TPMse_DnovPM)

subset(melOrths, mel_GeneSymbol == "hdly")

plotGenePM(TPMse_DnovPM, "GJ10165")
