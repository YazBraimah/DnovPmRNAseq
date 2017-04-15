# ____                   ____            ____  _   _    _                   
#|  _ \ _ __   _____   _|  _ \ _ __ ___ |  _ \| \ | |  / \   ___  ___  __ _ 
#| | | | '_ \ / _ \ \ / / |_) | '_ ` _ \| |_) |  \| | / _ \ / __|/ _ \/ _` |
#| |_| | | | | (_) \ V /|  __/| | | | | |  _ <| |\  |/ ___ \\__ \  __/ (_| |
#|____/|_| |_|\___/ \_/ |_|   |_| |_| |_|_| \_\_| \_/_/   \_\___/\___|\__, |
#                                                                        |_|

############# dvir1.06-based analysis

# Read in expression data
DnovPM.novTrin4.CountsMatrix = read.table("ExpressionData/gene.novTrin4.female.counts.matrix", header=T, row.names=1, com='', check.names=F)
DnovPM.novTrin4.TpmMatrix.cbmt = read.table("ExpressionData/gene.novTrin4.female.TPM.not_cross_norm.counts_by_min_TPM", header = T)
DnovPM.novTrin4.TmmMatrix = read.table("ExpressionData/gene.novTrin4.female.TMM.EXPR.matrix", header=T, row.names=1, com='', check.names=F)

# Read in sample info
DnovPM.Samples_data = read.table("ExpressionData/samples.txt", header=F, check.names=F, fill=T)
DnovPM.Samples_data = DnovPM.Samples_data[DnovPM.Samples_data[,2] != '',]

## Estimate of the number of expressed genes (Brian Haas' method)
# extract the data between 10 TPM and 100 TPM
DnovPM_filt_data = DnovPM.novTrin4.TpmMatrix.cbmt[DnovPM.novTrin4.TpmMatrix.cbmt[,1] > -100 & DnovPM.novTrin4.TpmMatrix.cbmt[,1] < -10,]
# perform a linear regression on this filtered subset of the data
DnovPM_fit = lm(DnovPM_filt_data[,2] ~ DnovPM_filt_data[,1])
print(DnovPM_fit)
# plot it
ggplot(DnovPM.dvir1.06.TpmMatrix.cbmt, aes(neg_min_tpm,num_features)) + 
    geom_point() + 
    geom_smooth(data=DnovPM_filt_data, method = "lm") + 
    scale_x_continuous(limits=c(-100,0)) + 
    scale_y_continuous(limits=c(0,20000)) + 
    geom_hline(yintercept = 12062, colour = "green") + ggtitle("expressed genes", subtitle = "something")

## Filter count data by minimum count across ANY sample
DnovPM_max_gene_expr_per_row = apply(DnovPM.novTrin4.CountsMatrix, 1, max)
DnovPM.novTrin4.CountsMatrix.min240count = DnovPM.novTrin4.CountsMatrix[DnovPM_max_gene_expr_per_row >= 240,,drop=F ]

## Barplot of gene counts by sample
libSizes <- as.data.frame(colSums(DnovPM.novTrin4.CountsMatrix.min240count))
libSizes <- cbind(sample = row.names(libSizes), libSizes)
row.names(libSizes)<- NULL
colnames(libSizes) = c("sample", "Total_reads")
ggplot(libSizes, aes(sample, Total_reads)) + 
    geom_bar(stat="identity") + 
    theme(axis.text.x = element_text(angle = -90, hjust = 0)) + 
    geom_hline(yintercept = 20000000)

## Boxplot of log10(TPM) across all samples
m.expData<-melt(as.matrix(DnovPM.novTrin4.CountsMatrix.min240count))
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

## calculate dispersion
d <- DGEList(counts = DnovPM.novTrin4.CountsMatrix.min240count, group = DnovPM.Samples_data$V1)
d <- calcNormFactors(d)
d <- estimateCommonDisp(d)
d <- estimateTagwiseDisp(d)
summary(d$tagwise.dispersion)
## Plot biological coefficient of variation
plotBCV(d)
## Plot grouping of samples
plotMDS(d, method = "bcv", col=as.numeric(d$samples$group))

## Plot sample correlation
data = log2(DnovPM.novTrin4.CountsMatrix.min240count+1)
data = as.matrix(data)
sample_cor = cor(data, method='pearson', use='pairwise.complete.obs')
sample_dist = as.dist(1-sample_cor)
hc_samples = hclust(sample_dist, method='complete')
heatmap.3(sample_cor, dendrogram='both', Rowv=as.dendrogram(hc_samples), Colv=as.dendrogram(hc_samples), col = colorpanel(75, '#dd70cd','black','#afc64f'), scale='none', symm=TRUE, key=TRUE,density.info='none', trace='none', symkey=FALSE, symbreaks=F, margins=c(10,10), cexCol=1, cexRow=1, cex.main=0.75, main=paste("sample correlation matrix"))

# normalize by DESeq method:
meta <- data.frame(row.names=colnames(DnovPM.novTrin4.CountsMatrix.min240count), condition=DnovPM.Samples_data$V1)
DnovPMCountData<-round(DnovPM.novTrin4.CountsMatrix.min240count)
DnovPMCountData_normByDESeq = newCountDataSet(DnovPMCountData, meta)
DnovPMCountData_normByDESeq = estimateSizeFactors(DnovPMCountData_normByDESeq)
DnovPMCountData_normByDESeq = data.frame(counts(DnovPMCountData_normByDESeq, normalized=T))

MA_BPlot(DnovPMCountData_normByDESeq, "H3_RT_2", "H3_RT_1")


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
