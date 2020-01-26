
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

## Load trinotate annotation info
trnTrinotate = read.csv("Annotations/DeNovo-based/Trinotate_report_novTrin4_subset.xls", header = T, sep = "\t", na.strings = "", stringsAsFactors=FALSE)
trnGO_info = read.table("Annotations/DeNovo-based/Trinotate_report.xls.gene_ontology", header=F, row.names=1,stringsAsFactors=F)

## subset annotation for fit function
annot.sum = unique(subset(trnTrinotate, select=c("Trin_ID", "sprot_Top_BLASTX_hit", "dvir1.06_pep_BLASTX")))
annot.sum = aggregate(cbind(sprot_Top_BLASTX_hit, dvir1.06_pep_BLASTX) ~ Trin_ID, data = annot.sum, toString)

## Load counts matrix for male and female data
countsMatrix = read.table("ExpressionData/DeNovo-based/gene.novTrin4.maleANDfemale.counts.matrix", header = T, row.names=1, com='', check.names = F)

## Load normalized TPM matrix
tpmMatrix = read.table("ExpressionData/DeNovo-based/gene.novTrin4.maleANDfemale.TMM.EXPR.matrix", header=T, row.names=1, com='', check.names=F)

## Load sample info (use combined male and female info)
sampleData = read.table("ExpressionData/DeNovo-based/samples.MF.txt", header=F, check.names=F, fill=T)
sampleData = sampleData[sampleData[,2] != '',]
sampleInfo = read.table("ExpressionData/DeNovo-based/SampleInfo.MF.txt", header=T, check.names=F, fill=T)

## Create mean TPM data tables for plotting and heatmap functions
tmp.tpmMatrix<-tpmMatrix
colnames(tmp.tpmMatrix) <- sampleData$V1
tmp.tpmMatrix.m <- as.data.frame(melt(as.matrix(tmp.tpmMatrix)))
tmp.tpmMatrix.m <- within(tmp.tpmMatrix.m, X2<-data.frame(do.call('rbind', strsplit(as.character(X2),'_',fixed=TRUE))))
tmp.tpmMatrix.m <- data.frame(tmp.tpmMatrix.m$X1, tmp.tpmMatrix.m$X2$X1, tmp.tpmMatrix.m$X2$X2, tmp.tpmMatrix.m$value)
colnames(tmp.tpmMatrix.m) <- c("Trin_ID", "sample", "tissue", "TPM")
tmp.tpmMatrix.m$condition <- ifelse(grepl("C", tmp.tpmMatrix.m$sample, ignore.case = F), 
                                    "conspecific", 
                                    ifelse(grepl("H", tmp.tpmMatrix.m$sample, ignore.case = F), 
                                           "heterospecific",
                                           ifelse(grepl("nov", tmp.tpmMatrix.m$sample, ignore.case = F),
                                                  "conspecific",
                                                  "virgin")))
tmp.tpmMatrix.m$sex <- ifelse(grepl("C24", tmp.tpmMatrix.m$sample, ignore.case = F), "male", "female")
tmp.tpmMatrix.m$time <- ifelse(grepl("3", tmp.tpmMatrix.m$sample), "3hpm", ifelse(grepl("6", tmp.tpmMatrix.m$sample), "6hpm", ifelse(grepl("12", tmp.tpmMatrix.m$sample), "12hpm", ifelse(grepl("24", tmp.tpmMatrix.m$sample), "24hpm","virgin"))))
tmp.tpmMatrix.m.c = summarySE(tmp.tpmMatrix.m, measurevar = "TPM", groupvars = c("Trin_ID", "sample", "tissue", "condition", "sex", "time"))
TPMse = tmp.tpmMatrix.m.c

## Use above mean calculation for matrix 
tmpMat<-cast(tmp.tpmMatrix.m.c, Trin_ID~sample+tissue, value ="TPM")
meanTPMmatrix <- tmpMat[,-1]
rownames(meanTPMmatrix) <- tmpMat[,1]

## Set the factor levels to order condition and time variables
TPMse$condition = factor (TPMse$condition, levels = c("virgin", "conspecific", "heterospecific"))
TPMse$time = factor (TPMse$time, levels = c("virgin", "3hpm", "6hpm", "12hpm", "24hpm"))


### Quality Assessment of sequence libraries

## Plot barplots of library sizes 
libSizes <- as.data.frame(colSums(countsMatrix))
libSizes <- cbind(sample = row.names(libSizes), libSizes)
row.names(libSizes)<- NULL
colnames(libSizes) = c("sample", "Total_reads")
options(repr.plot.width = 6, repr.plot.height = 3)
ggplot(libSizes, aes(sample, Total_reads)) + 
    geom_bar(stat="identity") + 
    theme(axis.text.x = element_text(angle = -90, hjust = 0)) + 
    geom_hline(yintercept = 20000000)


## Filter data based on minimum CPM (>1)
all.CPM <- cpm(countsMatrix)
thresh <- all.CPM > 1
keep <- rowSums(thresh) >= 2
counts.keep <- countsMatrix[keep,]
dim(counts.keep)


## Generate a counts matrix of virgin female samples only
countsMatrix.virgin = subset(countsMatrix, select=grepl("^V", colnames(countsMatrix)))
head(countsMatrix.virgin)

## Re-filter this matrix with minimum CPM (>1) 
cpmMatrix.virgin <- cpm(countsMatrix.virgin)
v.thresh <- cpmMatrix.virgin > 1
## look at number of columns with above requirements
table(rowSums(v.thresh))

v.keep <- rowSums(v.thresh) >= 2
countsMatrix.virgin.filt <- countsMatrix.virgin[v.keep,]
## Check how many genes remain in the matrix.
dim(countsMatrix.virgin.filt)

# Let's look at the first column
plot(cpmMatrix.virgin[,2],countsMatrix.virgin[,2],ylim=c(0,200),xlim=c(0,20))
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
ann.v = subset(annot.sum, Trin_ID %in% rownames(dgeList.v_fit))
# convert factors to characters
ann.v = data.frame(lapply(ann.v, as.character), stringsAsFactors=FALSE)
# align the fit object's rownames with gene ID's
ann.v = ann.v[match(rownames(dgeList.v_fit), ann.v$Trin_ID),]
# convert factors to characters, again
ann.v <- data.frame(lapply(ann.v, as.character), stringsAsFactors=FALSE)

# Check that the fit rownames match the annotation file's gene ID's
table(ann.v$Trin_ID==rownames(dgeList.v_fit))
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
repTract.list<-subset(lrt.v.repTract.tTags.table, logFC.V_RT.vs.V_CR > 2 & logFC.V_RT.vs.V_HD > 2 & logFC.V_RT.vs.V_OV > 2 & FDR<0.001)$Trin_ID
length(repTract.list)
head(lrt.v.repTract.tTags.table)

lrt.v.ovaries <- glmLRT(dgeList.v_fit, contrast = cont.v.ovaries)
lrt.v.ovaries.tTags <- topTags(lrt.v.ovaries, n = NULL)
lrt.v.ovaries.tTags.table <- lrt.v.ovaries.tTags$table
ovaries.list<-subset(lrt.v.ovaries.tTags.table, logFC.V_OV.vs.V_CR > 2 & logFC.V_OV.vs.V_HD > 2 & logFC.V_OV.vs.V_RT > 2 & FDR<0.001)$Trin_ID
length(ovaries.list)
head(lrt.v.ovaries.tTags.table)

lrt.v.head <- glmLRT(dgeList.v_fit, contrast = cont.v.head)
lrt.v.head.tTags <- topTags(lrt.v.head, n = NULL)
lrt.v.head.tTags.table <- lrt.v.head.tTags$table
head.list<-subset(lrt.v.head.tTags.table, logFC.V_H.vs.V_CR > 2 & logFC.V_H.vs.V_OV > 2 & logFC.V_H.vs.V_RT > 2 & FDR<0.001)$Trin_ID
length(head.list)
head(lrt.v.head.tTags.table)

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

options(repr.plot.width = 7, repr.plot.height = 5)
heatmap3(heatmap_data, 
         col=myheatcol, 
         cexCol = 2, 
         labRow = "", 
         legendfun=function() showLegend(legend=c("OV-biased", "H-biased", "RT-biased"), col=c("#e1526d", "#b85516", "#647700")), 
        cex.main=0.75, 
             RowSideColors = geneFactorAnnotations)

virginSpecificity = as.data.frame(calcSpecificity(tissueBiased.meanTPM.v))
head(virginSpecificity)

options(repr.plot.width = 5, repr.plot.height = 3)
plot(x=log10(tissueBiased.meanTPM.v$V_RT), y=virginSpecificity$V_RT)


