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
                 "qvalue", "reshape", "Rmisc", "splitstackshape", 
                 "statmod", "VennDiagram", "ggthemr", "Glimma")

lapply(req_packages, require, character.only = TRUE)

## The Cowplot package changes the default themes of ggplot2. Set the default theme like so:
theme_set(theme_gray())

## Load functions:
source("Functions2.R")

## Load annotations:
grpTrinotate = read.csv("Annotations/Trinotate_report_dvir1.06_subset.txt", header = T, sep = "\t", na.strings = "unknown", stringsAsFactors=FALSE)
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

#### QC and Differential expression

## Virgin Tissue comparisons

# Extract counts matrix for virgin tissue samples
counts.virgin.samples = subset(DnovPM.dvir1.06.CountsMatrix, select=grepl("^V", colnames(DnovPM.dvir1.06.CountsMatrix)))
head(counts.virgin.samples)

# Using the minimum CPM method to filter-out low expression genes
cpm.virgin.sample <- cpm(counts.virgin.samples)
v.thresh <- cpm.virgin.sample > 5
table(rowSums(v.thresh))
v.keep <- rowSums(v.thresh) >= 2
counts.virgin.keep <- counts.virgin.samples[v.keep,]
dim(counts.virgin.keep)

# Create DGE object and look at variance statistics
virgin.group <- factor(c(1,1,2,2,3,3,4,4,4))
virgin.design <- model.matrix(~0+virgin.group)
colnames(virgin.design)<-c("V_CR", "V_H", "V_OV", "V_RT")
virgin.DGElist<-DGEList(counts = counts.virgin.keep, group = virgin.group)
virgin.DGElist<-calcNormFactors(virgin.DGElist)
virgin.DGElist<-estimateDisp(virgin.DGElist, virgin.design, robust = T)
summary(virgin.DGElist$tagwise.dispersion)

virgin.fit <- glmFit(virgin.DGElist, virgin.design)

################################################################
################################################################
################################################################
################################################################

## Here's the voom method for DE
virgin.voom = voom(virgin.DGElist, virgin.design, plot=TRUE)
names(virgin.voom)

boxplot(virgin.voom$E, xlab="", ylab="Log2 counts per million",las=2,main="Voom transformed logCPM")
## Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(virgin.voom$E),col="blue")

# Fit the linear model
fit <- lmFit(virgin.voom)
names(fit)

## set-up contrasts:
cont.matrix <- makeContrasts(V_RT.vs.V_CR=V_RT-V_CR, 
                             V_RT.vs.V_H=V_RT-V_H,
                             V_RT.vs.V_OV=V_RT-V_OV, levels=virgin.design)
cont.matrix

fit.cont <- contrasts.fit(fit, cont.matrix)

fit.cont <- eBayes(fit.cont)
dim(fit.cont)

summa.fit <- decideTests(fit.cont)
summary(summa.fit)

vennDiagram(summa.fit)

fit.treat = treat(fit.cont, lfc = 1)
res.treat = decideTests(fit.treat)
summary(res.treat)

topTable(fit.treat, coef = 1, sort.by = "p")

# create an interactive plot:
glMDPlot(fit.treat, coef=1, counts=virgin.DGElist$counts, groups=virgin.group,
         status=res.treat, id.column="ENTREZID", main="B.PregVsLac",
         folder="md")

################################################################
################################################################
################################################################
################################################################

## Plot the biological coefficient of variation as a function of log CPM
plotBCV(virgin.DGElist)
g.virgin <- gof(virgin.fit)
z.virgin <- zscoreGamma(g.virgin$gof.statistics,shape=g.virgin$df/2,scale=2)
qqnorm(z.virgin); qqline(z.virgin, col = 4,lwd=4,lty=4)

## Plot two dimensional grouping of samples
plotMDS(virgin.DGElist, col=as.numeric(virgin.DGElist$samples$group))

## Plot sample correlation
data = log2(counts.virgin.samples+1)
data = as.matrix(data)
sample_cor = cor(data, method='pearson', use='pairwise.complete.obs')
sample_dist = as.dist(1-sample_cor)
hc_samples = hclust(sample_dist, method='complete')

heatmap.3(sample_cor, dendrogram='both', Rowv=as.dendrogram(hc_samples), Colv=as.dendrogram(hc_samples), col = colorpanel(75, '#dd70cd','black','#afc64f'), scale='none', symm=TRUE, key=TRUE, density.info='density', trace='none', symkey=FALSE, symbreaks=F, cexCol=1, cexRow=1, cex.main=0.75, main=paste("virgin samples correlation matrix"))

## Normalize by DESeq and generate MA plots to examine differences between replicates
meta <- data.frame(row.names=colnames(counts.virgin.keep), condition=substr(colnames(counts.virgin.keep), start=1, stop = nchar(colnames(counts.virgin.keep))-2))
virgin.counts.rounded<-round(counts.virgin.keep)
virgin.counts.rounded_normByDESeq = newCountDataSet(virgin.counts.rounded, meta)
virgin.counts.rounded_normByDESeq = estimateSizeFactors(virgin.counts.rounded_normByDESeq)
virgin.counts.rounded_normByDESeq = data.frame(counts(virgin.counts.rounded_normByDESeq, normalized=T))

MA_BPlot(virgin.counts.rounded_normByDESeq, "V_CR_1", "V_CR_2")
MA_BPlot(virgin.counts.rounded_normByDESeq, "V_H_1", "V_H_2")
MA_BPlot(virgin.counts.rounded_normByDESeq, "V_OV_1", "V_OV_2")
MA_BPlot(virgin.counts.rounded_normByDESeq, "V_RT_1", "V_RT_2")
MA_BPlot(virgin.counts.rounded_normByDESeq, "V_RT_1", "V_RT_3")
MA_BPlot(virgin.counts.rounded_normByDESeq, "V_RT_2", "V_RT_3")

## Generate tissue specificity matrix:
Dnov_virgin_tissue_MeanTPMmatrix <- subset(DnovPM_MeanTPMmatrix, select=c(FBgn_ID, V_CR, V_H, V_OV, V_RT))
rownames(Dnov_virgin_tissue_MeanTPMmatrix) <- Dnov_virgin_tissue_MeanTPMmatrix[,1]
Dnov_virgin_tissue_MeanTPMmatrix[,1] <- NULL
Dnov_virgin_Specificity_table = calcSpecificity(Dnov_virgin_tissue_MeanTPMmatrix)
Dnov_virgin_Specificity_table = as.data.frame(Dnov_virgin_Specificity_table)

## Identify tissue-biased transcripts
# first, define contrasts:
virgin_RT_contrasts<- makeContrasts(V_RT.vs.V_CR=V_RT-V_CR, 
                                    V_RT.vs.V_H=V_RT-V_H,
                                    V_RT.vs.V_OV=V_RT-V_OV,
                                    levels=virgin.design)

virgin_OV_contrasts<- makeContrasts(V_OV.vs.V_CR=V_OV-V_CR, 
                                    V_OV.vs.V_H=V_OV-V_H,
                                    V_OV.vs.V_RT=V_OV-V_RT,
                                    levels=virgin.design)

virgin_H_contrasts<- makeContrasts(V_H.vs.V_CR=V_H-V_CR, 
                                   V_H.vs.V_OV=V_H-V_OV,
                                   V_H.vs.V_RT=V_H-V_RT,
                                   levels=virgin.design)
## RT-biased genes
lrt.RT.v.rest <- glmLRT(virgin.fit, contrast = virgin_RT_contrasts)
lrt.RT.v.rest.tTags <- topTags(lrt.RT.v.rest, n = NULL)
lrt.RT.v.rest.tTags.table <- lrt.RT.v.rest.tTags$table
Dnov.dvir1.06.RT.list<-rownames(subset(lrt.RT.v.rest.tTags.table, logFC.V_RT.vs.V_CR > 2 & logFC.V_RT.vs.V_H > 2 & logFC.V_RT.vs.V_OV > 2 & FDR<0.001))

## OV-biased genes
lrt.OV.v.rest <- glmLRT(virgin.fit, contrast = virgin_OV_contrasts)
lrt.OV.v.rest.tTags <- topTags(lrt.OV.v.rest, n = NULL)
lrt.OV.v.rest.tTags.table <- lrt.OV.v.rest.tTags$table
Dnov.dvir1.06.OV.list<-rownames(subset(lrt.OV.v.rest.tTags.table, logFC.V_OV.vs.V_CR > 2 & logFC.V_OV.vs.V_H > 2 & logFC.V_OV.vs.V_RT > 2 & FDR<0.001))

## H-biased genes
lrt.H.v.rest <- glmLRT(virgin.fit, contrast = virgin_H_contrasts)
lrt.H.v.rest.tTags <- topTags(lrt.H.v.rest, n = NULL)
lrt.H.v.rest.tTags.table <- lrt.H.v.rest.tTags$table
Dnov.dvir1.06.H.list<-rownames(subset(lrt.H.v.rest.tTags.table, logFC.V_H.vs.V_CR > 2 & logFC.V_H.vs.V_OV > 2 & logFC.V_H.vs.V_RT > 2 & FDR<0.001))

## Generate a factor labeling matrix to annotate tissue-biased genes
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

virgin.factor.labeling = rbind(RT_factors, OV_factors, H_factors)
virgin.factor.labeling$FBgn_ID = row.names(virgin.factor.labeling)
colnames(virgin.factor.labeling) = c('FBgn_ID', 'tissue_bias')
virgin.factor_list = unique(factor_labeling[,1])

# obtain mean TPM natrix of tissue-biased transcripts
tissue.biased.matrix = subset(Dnov_virgin_tissue_MeanTPMmatrix, row.names(Dnov_virgin_tissue_MeanTPMmatrix) %in% row.names(virgin.factor.labeling))
tissue.biased.matrix$FBgn_ID = row.names(tissue.biased.matrix)

## combine mean TPM matrix of tissue biased transcripts with PAML data and plot along chromosome coordinates
trail = merge(tissue.biased.matrix, virgin.factor.labeling)
trail2 = subset(paml.data, FBgn_ID %in% trail$FBgn_ID & omega < 800 & grepl("Chr", chromosome))
trail3 = merge (trail, trail2)

ggplot(trail3, aes(max, omega, colour = tissue_bias)) + 
    geom_point(size=2, alpha=0.5) + 
    geom_hline(yintercept = 0.15, linetype="dashed", colour = "yellow") + 
    geom_hline(yintercept = 1, linetype="dashed", colour = "gray")  + 
    facet_grid(tissue_bias~chromosome, scales = "free") +
    scale_x_continuous(breaks=seq(5000000,30000000,5000000), labels=expression("5", "10", "15", "20", "25", "30")) + 
    xlab ("Chromosome coordinates (Mb)") + 
    labs(y=expression(K[a]/K[s])) + 
    theme(axis.title.x = element_text(face = "bold", size = 10, vjust=0.1), axis.text.x=element_text(face = "bold", size = 12),axis.text.y = element_text(face = "bold", size = 12), axis.title.y = element_text(face = "bold.italic", size = 12, vjust=0.1), strip.text=element_text(face="bold", size = 12))

## Gene Ontology analysis
# get gene lengths and GO info
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

# build pwf based on ALL DE features
options(repr.plot.width = 6, repr.plot.height = 3)

cat_genes_vec = as.integer(features_with_GO %in% rownames(factor_labeling))
pwf=nullp(cat_genes_vec,bias.data=lengths_features_with_GO)
rownames(pwf) = names(GO_info_listed)

# Perform GO analysis
GO_enriched_list = list()

for (feature_cat in factor_list) {
    message('Processing category: ', feature_cat)
    cat_genes_vec = as.integer(features_with_GO %in% rownames(virgin.factor.labeling)[virgin.factor.labeling$tissue_bias == feature_cat])
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

