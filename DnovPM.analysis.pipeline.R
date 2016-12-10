
############# dvir1.06-based analysis
grpTrinotate<-read.csv(file = "~/Dropbox/RNAseq/Male_RNAseq/eXpress.analyses/MaleRNAseq.R.files/Annotations/Trinotate_report_dvir1.06_subset.txt", sep = "\t", header = T, na.strings = ".")
gffRecord <- read.table(file = "~/Dropbox/RNAseq/Male_RNAseq/eXpress.analyses/MaleRNAseq.R.files/Annotations/FBgn_ID_name_coordinates.txt", header = T)
melOrths <- read.table(file = "~/Dropbox/RNAseq/Male_RNAseq/eXpress.analyses/MaleRNAseq.R.files/Annotations/mel_orths.txt", header = T)
melOrthsAll<-aggregate(mel_GeneSymbol~FBgn_ID, data = melOrths, toString)
tmp.merged <- merge(melOrthsAll, grpTrinotate, all=TRUE)
Annots <- merge(tmp.merged, gffRecord, all=TRUE)
rm(grpTrinotate, melOrthsAll, tmp.merged)

DnovPM.dvir1.06.CountsMatrix = read.table("ExpressionData/genes_DnovPM_dvi1.06.counts.matrix", header=T, row.names=1, com='', check.names=F)
DnovPM.dvir1.06.TpmMatrix.cbmt = read.table("ExpressionData/genes_DnovPM_dvi1.06.TPM.not_cross_norm.counts_by_min_TPM", header = T)
DnovPM.dvir1.06.TmmMatrix = read.table("ExpressionData/genes_DnovPM_dvi1.06.TMM.EXPR.matrix", header=T, row.names=1, com='', check.names=F)

DnovPM.Samples_data = read.table("ExpressionData/samples.txt", header=F, check.names=F, fill=T)
DnovPM.Samples_data = DnovPM.Samples_data[DnovPM.Samples_data[,2] != '',]

## Barplot of gene counts by sample
libSizes <- as.data.frame(colSums(DnovPM.dvir1.06.CountsMatrix))
libSizes <- cbind(sample = row.names(libSizes), libSizes)
row.names(libSizes)<- NULL
colnames(libSizes) = c("sample", "Total_reads")
ggplot(libSizes, aes(sample, Total_reads)) + geom_bar(stat="identity") + theme(axis.text.x = element_text(angle = -90, hjust = 0)) + geom_hline(yintercept = 20000000)


## Boxplot of log10(TPM) across all samples
m.expData<-melt(as.matrix(DnovPM.dvir1.06.TmmMatrix))
colnames(m.expData) <- c("gene_id", "replicate", "TPM")
m.expData.exp<- within(m.expData, replicate<-data.frame(do.call('rbind', strsplit(as.character(replicate),'_',fixed=TRUE))))
m.expData<-data.frame(m.expData, m.expData.exp$replicate$X1, m.expData.exp$replicate$X2, m.expData.exp$replicate$X3)
colnames(m.expData) <- c("gene_id", "replicate", "TPM", "sample", "tissue", "rep_num")
m.expData$TPM <- m.expData$TPM + 1
p <- ggplot(m.expData)
p <- p + geom_boxplot(aes(x = replicate, y = log10(TPM), fill = tissue, colour = sample), size = 0.3, alpha = I(1/3))
p <- p + theme(axis.text.x = element_text(angle = -90, hjust = 0))
p <- p + scale_fill_hue(l = 50, h.start = 200)
p


## Estimate of the number of expressed genes (Brian Haas' method)
# extract the data between 10 TPM and 100 TPM
DnovPM_filt_data = DnovPM.dvir1.06.TpmMatrix.cbmt[DnovPM.dvir1.06.TpmMatrix.cbmt[,1] > -100 & DnovPM.dvir1.06.TpmMatrix.cbmt[,1] < -10,]
# perform a linear regression on this filtered subset of the data
DnovPM_fit = lm(DnovPM_filt_data[,2] ~ DnovPM_filt_data[,1])
print(DnovPM_fit)
# plot it
ggplot(DnovPM.dvir1.06.TpmMatrix.cbmt, aes(neg_min_tpm,num_features)) + geom_point() + scale_x_continuous(limits=c(-100,0)) + scale_y_continuous(limits=c(0,20000)) + geom_smooth(data=DnovPM_filt_data, method = "lm") + geom_hline(yintercept = 9185, colour = "green") + ggtitle("dvir1.06")

## calculate dispersion
d <- DGEList(counts = DnovPM.dvir1.06.CountsMatrix, group = DnovPM.Samples_data$V1)
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
heatmap.3(sample_cor, dendrogram='both', Rowv=as.dendrogram(hc_samples), Colv=as.dendrogram(hc_samples), col = greenred(75), scale='none', symm=TRUE, key=TRUE,density.info='none', trace='none', symkey=FALSE, symbreaks=F, margins=c(10,10), cexCol=1, cexRow=1, cex.main=0.75, main=paste("sample correlation matrix"))

# normalize by DESeq method:
meta <- data.frame(row.names=colnames(DnovPM.dvir1.06.CountsMatrix), condition=DnovPM.Samples_data$V1)
DnovPMCountData<-round(DnovPM.dvir1.06.CountsMatrix)
DnovPMCountData_normByDESeq = newCountDataSet(DnovPMCountData, meta)
DnovPMCountData_normByDESeq = estimateSizeFactors(DnovPMCountData_normByDESeq)
DnovPMCountData_normByDESeq = data.frame(counts(DnovPMCountData_normByDESeq, normalized=T))

#MA_BPlot(DnovPMCountData_normByDESeq, "V_RT_2", "V_RT_1")

## Filter count data by minimum count across ANY sample
DnovPM_max_gene_expr_per_row = apply(DnovPM.dvir1.06.CountsMatrix, 1, max)
DnovPM.dvir1.06.CountsMatrix.min400count = DnovPM.dvir1.06.CountsMatrix[DnovPM_max_gene_expr_per_row >= 400,,drop=F ]

#########################################################################################
### Summary TPM table and matrix for gene level plots (includes all replicates) ######### 

DnovPM.TPM.tmp<-DnovPM.dvir1.06.TmmMatrix
colnames(DnovPM.TPM.tmp) <- DnovPM.Samples_data$V1
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
plotGenePM(TPMse_DnovPM, "FBgn0283024")

## Generate tissue specificity matrix:
Dnov_virgin_tissue_MeanTPMmatrix <- subset(DnovPM_MeanTPMmatrix, select=c(FBgn_ID, V_CR, V_H, V_OV, V_RT))
rownames(Dnov_virgin_tissue_MeanTPMmatrix) <- Dnov_virgin_tissue_MeanTPMmatrix[,1]
Dnov_virgin_tissue_MeanTPMmatrix[,1] <- NULL
Dnov_virgin_Specificity_table = YazSpecificity(Dnov_virgin_tissue_MeanTPMmatrix)
Dnov_virgin_Specificity_table = as.data.frame(Dnov_virgin_Specificity_table)

########################################## edgeR analysis ###################################################
## Define groups for DE analysis.
# simple way (for pair-wise comparisons):

DnovPM.group <- factor(c(1,1,1,2,2,2,3,3,4,4,5,5,5,6,6,6,7,7,7,8,8,9,9,10,10,10,11,11,12,12,13,13,14,14,14))
DnovPM.design <- model.matrix(~0+DnovPM.group)
colnames(DnovPM.design)<-levels(d$samples$group)
DnovPM_ExpStd<-DGEList(counts = DnovPM.dvir1.06.CountsMatrix.min400count, group = DnovPM.group)
DnovPM_ExpStd<-calcNormFactors(DnovPM_ExpStd)
DnovPM_ExpStd<-estimateDisp(DnovPM_ExpStd, DnovPM.design, robust = T)
DnovPM_fit <- glmFit(DnovPM_ExpStd, DnovPM.design)

virgin_RT_contrasts<- makeContrasts(V_RT.vs.V_CR=V_RT-V_CR, 
                                        V_RT.vs.V_H=V_RT-V_H,
                                        V_RT.vs.V_OV=V_RT-V_OV,
                                        levels=DnovPM.design)

lrt.RT.v.rest <- glmLRT(DnovPM_fit, contrast = virgin_RT_contrasts)
lrt.RT.v.rest.tTags <- topTags(lrt.RT.v.rest, n = NULL)
lrt.RT.v.rest.tTags.table <- lrt.RT.v.rest.tTags$table
Dnov.dvir1.06.RT.list<-rownames(subset(lrt.RT.v.rest.tTags.table, logFC.V_RT.vs.V_CR > 2 & logFC.V_RT.vs.V_H > 2 & logFC.V_RT.vs.V_OV > 2 & FDR<0.001))

Dnov.dvir1.06.RT.matrix <- subset(Dnov_virgin_tissue_MeanTPMmatrix, rownames(Dnov_virgin_tissue_MeanTPMmatrix) %in% Dnov.dvir1.06.RT.list)
YazHeatmap(Dnov.dvir1.06.RT.matrix, clustering = "both", labRow = T)

# pdf("Plots/Dnov_RT-biased_barPlots.pdf", height = 4)
# lapply(Dnov.dvir1.06.RT.list, plotGenePM, object = TPMse_DnovPM)
# dev.off()

ggplot(subset(paml.data, FBgn_ID %in% Dnov.dvir1.06.RT.list & omega < 800 & grepl("Chr", chromosome)), aes(max, omega)) + geom_point(size=2, alpha=0.5, colour = "#7d49c3") + geom_point(data=subset(paml.data, FBgn_ID %in% SFP_elements$`D.ame,D.lum,D.nov,D.vir` ), aes(max, omega), inherit.aes = F, size=2, alpha=0.5, colour = "#4f922a") + geom_hline(yintercept = 0.15, linetype="dashed", colour = "yellow") + geom_hline(yintercept = 1, linetype="dashed", colour = "gray")  + facet_grid(~chromosome, scales = "free_x") + scale_colour_manual(name = "", values =c("#7aa457"="#7aa457","#9e6ebd"="#9e6ebd"), labels = c("SFPs","EB biased")) + scale_x_continuous(breaks=c(5000000, 10000000, 15000000, 20000000, 25000000, 30000000), labels=expression("5", "10", "15", "20", "25", "30")) + xlab ("Chromosome coordinates (Mb)") + labs(y=expression(K[a]/K[s])) + geom_text_repel(data=subset(paml.data, FBgn_ID %in% Dnov.dvir1.06.RT.list & omega > 0.95 & omega < 800 & chromosome == "Chr_2"), aes(label = gene_name), size =3, force = 30, colour = "#7d49c3") + geom_text_repel(data=subset(paml.data, FBgn_ID %in% SFP_elements$`D.ame,D.lum,D.nov,D.vir` & omega > 0.8), aes(label = gene_name), size =3, force = 4, colour = "#4f922a") + theme(axis.title.x = element_text(face = "bold", size = 10, vjust=0.1), axis.text.x=element_text(face = "bold", size = 12),axis.text.y = element_text(face = "bold", size = 12), axis.title.y = element_text(face = "bold.italic", size = 12, vjust=0.1), strip.text=element_text(face="bold", size = 12))

# better way
df1 = subset(m.DnovPM.TPM.tmp, FBgn_ID == "FBgn0202928")
df1$replicates = DnovPM.Samples_data$V2
rownames(df1) <- df1$replicates
df1 = subset(df1, select = c(tissue, condition, time))


cbind(df1,Group=RT.Group)
design <- model.matrix(~0+RT.Group)
colnames(design) <- levels(RT.Group)

########### RT contrasts
# create RT-specific count matrix and define groups/design
DnovPM_CountsMatrix_min400_RT = subset(DnovPM.dvir1.06.CountsMatrix.min400count, select=grepl("RT", colnames(DnovPM.dvir1.06.CountsMatrix.min400count)))
RT.group <- factor(c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6,7,7,7))
RT.design <- model.matrix(~0+RT.group)
colnames(RT.design)<-c("C12_RT", "C3_RT", "C6_RT", "H12_RT", "H3_RT", "H6_RT", "V_RT")

# create edgeR DE object and run glm
DnovPM_DGElist_RT<-DGEList(counts = DnovPM_CountsMatrix_min400_RT, group = RT.group)
DnovPM_DGElist_RT<-calcNormFactors(DnovPM_DGElist_RT)
DnovPM_DGElist_RT<-estimateDisp(DnovPM_DGElist_RT, RT.design, robust = T)
DnovPM_RT_fit <- glmFit(DnovPM_DGElist_RT, RT.design)

# define contrasts
con_virgin_contrasts <- makeContrasts(C3.vs.V=C3_RT-V_RT, C6.vs.V=C6_RT-V_RT, C12.vs.V=C12_RT-V_RT, levels = RT.design)
het_virgin_contrasts <- makeContrasts(H3.vs.V=H3_RT-V_RT, H6.vs.V=H6_RT-V_RT, H12.vs.V=H12_RT-V_RT, levels = RT.design)

# identify DE genes: C3_RT vs V_RT
RT_con.3hrs.vs.virgin <- glmLRT(DnovPM_RT_fit, contrast = con_virgin_contrasts[,"C3.vs.V"])
RT_con.3hrs.vs.virgin.tTags <- topTags(RT_con.3hrs.vs.virgin, n = NULL)
RT_con.3hrs.vs.virgin.tTags.table <- RT_con.3hrs.vs.virgin.tTags$table
RT_con.3hrs.vs.virgin.list <- rownames(subset(RT_con.3hrs.vs.virgin.tTags.table, logFC > 1 & FDR < 0.05))
# identify DE genes: C6_RT vs V_RT
RT_con.6hrs.vs.virgin <- glmLRT(DnovPM_RT_fit, contrast = con_virgin_contrasts[,"C6.vs.V"])
RT_con.6hrs.vs.virgin.tTags <- topTags(RT_con.6hrs.vs.virgin, n = NULL)
RT_con.6hrs.vs.virgin.tTags.table <- RT_con.6hrs.vs.virgin.tTags$table
RT_con.6hrs.vs.virgin.list <- rownames(subset(RT_con.6hrs.vs.virgin.tTags.table, logFC > 1 & FDR < 0.05))
# identify DE genes: C12_RT vs V_RT
RT_con.12hrs.vs.virgin <- glmLRT(DnovPM_RT_fit, contrast = con_virgin_contrasts[,"C12.vs.V"])
RT_con.12hrs.vs.virgin.tTags <- topTags(RT_con.12hrs.vs.virgin, n = NULL)
RT_con.12hrs.vs.virgin.tTags.table <- RT_con.12hrs.vs.virgin.tTags$table
RT_con.12hrs.vs.virgin.list <- rownames(subset(RT_con.12hrs.vs.virgin.tTags.table, logFC > 1 & FDR < 0.05))
# identify DE genes: H3_RT vs V_RT
RT_het.3hrs.vs.virgin <- glmLRT(DnovPM_RT_fit, contrast = het_virgin_contrasts[,"H3.vs.V"])
RT_het.3hrs.vs.virgin.tTags <- topTags(RT_het.3hrs.vs.virgin, n = NULL)
RT_het.3hrs.vs.virgin.tTags.table <- RT_het.3hrs.vs.virgin.tTags$table
RT_het.3hrs.vs.virgin.list <- rownames(subset(RT_het.3hrs.vs.virgin.tTags.table, logFC > 1 & FDR < 0.05))
# identify DE genes: H6_RT vs V_RT
RT_het.6hrs.vs.virgin <- glmLRT(DnovPM_RT_fit, contrast = het_virgin_contrasts[,"H6.vs.V"])
RT_het.6hrs.vs.virgin.tTags <- topTags(RT_het.6hrs.vs.virgin, n = NULL)
RT_het.6hrs.vs.virgin.tTags.table <- RT_het.6hrs.vs.virgin.tTags$table
RT_het.6hrs.vs.virgin.list <- rownames(subset(RT_het.6hrs.vs.virgin.tTags.table, logFC > 1 & FDR < 0.05))
# identify DE genes: H12_RT vs V_RT
RT_het.12hrs.vs.virgin <- glmLRT(DnovPM_RT_fit, contrast = het_virgin_contrasts[,"H12.vs.V"])
RT_het.12hrs.vs.virgin.tTags <- topTags(RT_het.12hrs.vs.virgin, n = NULL)
RT_het.12hrs.vs.virgin.tTags.table <- RT_het.12hrs.vs.virgin.tTags$table
RT_het.12hrs.vs.virgin.list <- rownames(subset(RT_het.12hrs.vs.virgin.tTags.table, logFC > 1 & FDR < 0.05))

RT_UP_candidates <- list(Up.at.C3.RT = RT_con.3hrs.vs.virgin.list, 
                          Up.at.C6.RT = RT_con.6hrs.vs.virgin.list, 
                          Up.at.C12.RT = RT_con.12hrs.vs.virgin.list,
                         Up.at.H3.RT = RT_het.3hrs.vs.virgin.list, 
                         Up.at.H6.RT = RT_het.6hrs.vs.virgin.list, 
                         Up.at.H12.RT = RT_het.12hrs.vs.virgin.list)

RT_UP_combs <- unlist(lapply(1:length(RT_UP_candidates), function(j) combn(names(RT_UP_candidates), j, simplify = FALSE)), recursive = FALSE)
names(RT_UP_combs) <- sapply(RT_UP_combs, function(i) paste0(i, collapse = ","))
RT_UP_elements <- lapply(RT_UP_combs, function(i) Setdiff(RT_UP_candidates[i], RT_UP_candidates[setdiff(names(RT_UP_candidates), i)]))

RT_UP_candidates_Vdiag<-venn.diagram(RT_UP_candidates, NULL, fill=c("#7d35ae", "#c33d92", "#c471d1", "#c7415f", "#d46c41", "#c8392a"), alpha=c(0.5,0.5,0.5,0.5,0.5,0.5), cex = 1.5, cat.fontface= 4, cat.cex = 1.25, resolution = 1000, main = "RT Upregulated")
grid.arrange(gTree(children=RT_con_candidates_Vdiag))


# Test difference between con- and heterospecific
condition.group = subset(m.DnovPM.TPM.tmp, FBgn_ID == "FBgn0202928" & tissue == "RT")
condition.group$replicates = colnames(DnovPM_CountsMatrix_min400_RT)
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
con.vs.het.RT.all.list <- rownames(subset(con.vs.het.RT.all.tTags.table, logFC < -1 & FDR < 0.05))

het.vs.vir.RT.all.contrast <- glmLRT(lrt, contrast = condition.contrasts[,"het.vs.vir"])
het.vs.vir.RT.all.tTags <- topTags(het.vs.vir.RT.all.contrast, n = NULL)
het.vs.vir.RT.all.tTags.table <- het.vs.vir.RT.all.tTags$table
het.vs.vir.RT.all.list <- rownames(subset(het.vs.vir.RT.all.tTags.table, logFC > 1 & FDR < 0.05))

con.vs.vir.RT.all.contrast <- glmLRT(lrt, contrast = condition.contrasts[,"con.vs.vir"])
con.vs.vir.RT.all.tTags <- topTags(con.vs.vir.RT.all.contrast, n = NULL)
con.vs.vir.RT.all.tTags.table <- con.vs.vir.RT.all.tTags$table
con.vs.vir.RT.all.list <- rownames(subset(con.vs.vir.RT.all.tTags.table, logFC > 1 & FDR < 0.05))

PM.vs.vir_candidates <- list(conspecific = con.vs.vir.RT.all.list, heterospecific = het.vs.vir.RT.all.list, between = con.vs.het.RT.all.list)

PM.vs.vir_candidates_Vdiag<-venn.diagram(PM.vs.vir_candidates, NULL, fill=c("#7d35ae", "#c7415f", "#d46c41"), alpha=c(0.5,0.5, 0.5), cex = 1.5, cat.fontface= 4, cat.cex = 1.25, resolution = 1000, main = "RT Upregulated")
grid.arrange(gTree(children=PM.vs.vir_candidates_Vdiag))

pdf("Plots/con.vs.het.RT.all.list.pdf", height = 4)
lapply(con.vs.het.RT.all.list, plotGenePM_RT, object = TPMse_DnovPM)
dev.off()

C2inv.qtl = data.frame(xmin=17747413.5, xmax=34500000, ymin=0, ymax = 1.5, chromosome = "Chr_2")
C5.qtl = data.frame(xmin=c(13800000, 16300000, 22800000), xmax=c(14750000, 21700000, 25000000), ymin=0, ymax = 1.5, chromosome = "Chr_5")

ggplot(subset(paml.data, FBgn_ID %in% Dnov.dvir1.06.RT.list & omega < 800 & grepl("Chr", chromosome)), aes(max, omega)) + geom_point(size=2, alpha=0.75, aes(colour = "#7d49c3")) + geom_point(data=subset(paml.data, FBgn_ID %in% SFP_elements$`D.ame,D.lum,D.nov,D.vir` ), aes(max, omega, colour = "#4f922a"), size=2, alpha=0.75) + geom_hline(yintercept = 0.15, linetype="dashed", colour = "red") + geom_hline(yintercept = 1, linetype="dashed", colour = "black") + scale_colour_manual(name = "", values =c("#7d49c3"="#7d49c3","#4f922a"="#4f922a"), labels = c("SFPs","Female RT")) + facet_grid(~chromosome, scales = "free_x") + xlab ("Chromosome coordinates (Mb)") + labs(y=expression(K[a]/K[s])) + theme(axis.title.x = element_text(face = "bold", size = 10, vjust=0.1), axis.text.x=element_text(face = "bold", size = 12),axis.text.y = element_text(face = "bold", size = 12), axis.title.y = element_text(face = "bold.italic", size = 12, vjust=0.1), strip.text=element_text(face="bold", size = 12), legend.text = element_text(size = 12, face="bold")) + geom_rect(data=C2inv.qtl, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill="red"), alpha=0.1, inherit.aes = FALSE) + geom_rect(data=C5.qtl, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="red", alpha=0.1, inherit.aes = FALSE) + scale_fill_manual(name = '', values = "red",labels = "PMPZ\nQTL region")+ scale_x_continuous(breaks=c(5000000, 10000000, 15000000, 20000000, 25000000, 30000000), labels=expression("5", "10", "15", "20", "25", "30"))

ggplot(subset(paml.data, FBgn_ID %in% Dnov.dvir1.06.RT.list & omega < 800 & chromosome == "Chr_2"), aes(max, omega)) + geom_point(size=2, alpha=0.75, aes(colour = "#7d49c3")) + geom_point(data=subset(paml.data, FBgn_ID %in% SFP_elements$`D.ame,D.lum,D.nov,D.vir` ), aes(max, omega, colour = "#4f922a"), size=2, alpha=0.75) + geom_hline(yintercept = 0.15, linetype="dashed", colour = "black") + geom_hline(yintercept = 1, linetype="dashed", colour = "red") + scale_colour_manual(name = "", values =c("#7d49c3"="#7d49c3","#4f922a"="#4f922a"), labels = c("SFPs","Female RT")) + scale_x_continuous(breaks=c(5000000, 10000000, 15000000, 20000000, 25000000, 30000000), labels=expression("5", "10", "15", "20", "25", "30")) + xlab ("Chromosome coordinates (Mb)") + labs(y=expression(K[a]/K[s])) + geom_text(data=subset(paml.data, FBgn_ID %in% Dnov.dvir1.06.RT.list & omega > 0.95 & omega < 800 & chromosome == "Chr_2"), aes(label = gene_name), size =4, colour = "#7d49c3", fontface = "bold.italic", vjust = 1.4, hjust = 0.1) + geom_text(data=subset(paml.data, FBgn_ID %in% SFP_elements$`D.ame,D.lum,D.nov,D.vir` & omega > 0.8), aes(label = gene_name), size =4, colour = "#4f922a", fontface = "bold.italic", hjust = 1.1) + theme(axis.title.x = element_text(face = "bold", size = 10, vjust=0.1), axis.text.x=element_text(face = "bold", size = 12),axis.text.y = element_text(face = "bold", size = 12), axis.title.y = element_text(face = "bold.italic", size = 12, vjust=0.1), strip.text=element_text(face="bold", size = 12), plot.title=element_text(face="bold"), legend.text = element_text(size = 12, face="bold")) + ggtitle("Chromosome 2 (Muller element E)") + geom_rect(data=C2inv.qtl, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="red", alpha=0.1, inherit.aes = FALSE) + scale_fill_manual(name = '', values = "red",labels = "PMPZ\nQTL region")

