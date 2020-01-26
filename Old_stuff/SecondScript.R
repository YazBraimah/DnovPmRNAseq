# ____                   ____            ____  _   _    _                   
#|  _ \ _ __   _____   _|  _ \ _ __ ___ |  _ \| \ | |  / \   ___  ___  __ _ 
#| | | | '_ \ / _ \ \ / / |_) | '_ ` _ \| |_) |  \| | / _ \ / __|/ _ \/ _` |
#| |_| | | | | (_) \ V /|  __/| | | | | |  _ <| |\  |/ ___ \\__ \  __/ (_| |
#|____/|_| |_|\___/ \_/ |_|   |_| |_| |_|_| \_\_| \_/_/   \_\___/\___|\__, |
#                                                                        |_|

############# dvir1.06-based analysis


# Read in expression data
DnovPM.Trin4.CountsMatrix = read.table("ExpressionData/gene.novTrin4.female.counts.matrix", header=T, row.names=1, com='', check.names=F)
DnovPM.Trin4.TpmMatrix.cbmt = read.table("ExpressionData/genes_DnovPM_dvi1.06.TPM.not_cross_norm.counts_by_min_TPM", header = T)
DnovPM.Trin4.TmmMatrix = read.table("ExpressionData/genes_DnovPM_dvi1.06.TMM.EXPR.matrix", header=T, row.names=1, com='', check.names=F)


hatchData = read.csv("~/Dropbox/RNAseq/PostMating.RNAseq/Classical_genetics/Hatch_rate_data.csv", header = T, sep = ",")

myHueSwatch = c("#c05243",
                "#8dd055",
                "#784ac0",
                "#cda754",
                "#cc53a1",
                "#87cbaf",
                "#4e304b",
                "#56633c",
                "#9b98bf")
somTheme = define_palette(swatch = myHueSwatch, gradient = c("#a58d48", "#a670b0"), background = "#1f1317", text = c("#fff1f4", "#c6e395"), line = c("#fff1f4", "#ffd7c5"), gridline = "#6d5e63")
ggthemr(somTheme, spacing = 1, type = "outer", layout = "scientific")

ggplot(hatchData, aes(Female, Hatch_rate, fill = Cross_Type)) + geom_bar(stat = "identity") + geom_errorbar(aes(ymin=Hatch_rate-se, ymax=Hatch_rate+se), width=.2, position=position_dodge(.9)) + facet_grid(~Male) + theme(axis.text.x = element_text(angle=45, vjust = 0.5)) 

ggthemr_reset()


het3hrs = RT_het.3hrs.vs.virgin.tTags.table
het3hrs$time = "3hpm"
het3hrs$cross = "heterospecific"
het3hrs$FBgn_ID = rownames(het3hrs)

het6hrs = RT_het.6hrs.vs.virgin.tTags.table
het6hrs$time = "6hpm"
het6hrs$cross = "heterospecific"
het6hrs$FBgn_ID = rownames(het6hrs)

het12hrs = RT_het.12hrs.vs.virgin.tTags.table
het12hrs$time = "12hpm"
het12hrs$cross = "heterospecific"
het12hrs$FBgn_ID = rownames(het12hrs)

con3hrs = RT_con.3hrs.vs.virgin.tTags.table
con3hrs$time = "3hpm"
con3hrs$cross = "conspecific"
con3hrs$FBgn_ID = rownames(con3hrs)

con6hrs = RT_con.6hrs.vs.virgin.tTags.table
con6hrs$time = "6hpm"
con6hrs$cross = "conspecific"
con6hrs$FBgn_ID = rownames(con6hrs)

con12hrs = RT_con.12hrs.vs.virgin.tTags.table
con12hrs$time = "12hpm"
con12hrs$cross = "conspecific"
con12hrs$FBgn_ID = rownames(con12hrs)



pairWise.DE.data = rbind(con3hrs, con6hrs, con12hrs, het3hrs, het6hrs, het12hrs)
pairWise.DE.data$time = factor(pairWise.DE.data$time, levels = c("3hpm", "6hpm", "12hpm"))
pairWise.DE.data$sig = ifelse(pairWise.DE.data$FDR < 0.001 & pairWise.DE.data$logFC > 1 | pairWise.DE.data$FDR < 0.001 & pairWise.DE.data$logFC < -1 , "YES", "NO")
pairWise.DE.data$direction = ifelse(pairWise.DE.data$logFC > 0, "Upregulated", "Downregulated")
ggplot(pairWise.DE.data, aes(logFC, -log10(PValue), colour = sig)) + geom_point(alpha = 0.75) + facet_grid(time~cross, scales = "free") + theme_solarized() +scale_colour_manual(values=c("#94008a","#009b79"))


ggplot(subset(pairWise.DE.data, sig == "YES"), aes(logFC, fill = cross)) + geom_histogram(alpha = 0.5, binwidth = 0.1) + facet_grid(time~direction, scales = "free") + 

ggplot(con.vs.het.RT.all.tTags.table, aes(logFC, -log10(PValue))) + geom_point(alpha = 0.75) + theme_solarized() +scale_colour_manual(values=c("#94008a","#009b79"))


#####################
con.vs.vir=con.vs.vir.RT.all.tTags.table
con.vs.vir$comparison = "con vs. vir"
con.vs.vir$FBgn_ID = rownames(con.vs.vir)

het.vs.vir=het.vs.vir.RT.all.tTags.table
het.vs.vir$comparison = "het vs. vir"
het.vs.vir$FBgn_ID = rownames(het.vs.vir)

condition.DE.data = rbind(con.vs.vir, het.vs.vir)
condition.DE.data$sig = ifelse(condition.DE.data$FDR < 0.001 & condition.DE.data$logFC > 1 | condition.DE.data$FDR < 0.001 & condition.DE.data$logFC < -1 , "YES", "NO")
ggplot(condition.DE.data, aes(logFC, -log10(PValue), colour = sig)) + geom_point(alpha = 0.75) + facet_grid(~comparison, scales = "free") + theme_solarized_2() +scale_colour_manual(values=c("#94008a","#009b79"))

#########################

het.v.con = con.vs.het.RT.all.tTags.table
het.v.con$FBgn_ID = rownames(het.v.con)
het.v.con$sig = ifelse(het.v.con$FDR < 0.001 & het.v.con$logFC > 1 | het.v.con$FDR < 0.001 & het.v.con$logFC < -1 , "YES", "NO")
ggplot(het.v.con, aes(logFC, -log10(PValue), colour = sig)) + geom_point(alpha = 0.75) + theme_solarized_2() +scale_colour_manual(values=c("#94008a","#009b79"))

#############
#############
genome.omega = subset(paml.data, omega != "NA" & omega < 100)
genome.omega = subset(genome.omega, select=c(FBtr_ID, gene_name, chromosome, omega))
genome.omega$Class = "All genes"

RT.omega = subset(paml.data, FBgn_ID %in% Dnov.dvir1.06.RT.list & omega != "NA" & omega < 100)
RT.omega = subset(RT.omega, select=c(FBtr_ID, gene_name, chromosome, omega))
RT.omega$Class = "RT biased"

SFP.omega = subset(paml.data, FBgn_ID %in% SFP_elements$`D.ame,D.lum,D.nov,D.vir` & omega != "NA" & omega < 100)
SFP.omega = subset(SFP.omega, select=c(FBtr_ID, gene_name, chromosome, omega))
SFP.omega$Class = "SFPs"

omegaData.df = rbind(genome.omega, RT.omega, SFP.omega)
omegaData.df$Class = factor (omegaData.df$Class, levels = c("All genes", "RT biased", "SFPs"))

genome.avg.omega = mean(subset(omegaData.df, Class == "All genes" & omega != "NA" & omega < 100)$omega)
genome.se.omega = sd(subset(omegaData.df, Class == "All genes" & omega != "NA" & omega < 100)$omega)/sqrt(length(subset(omegaData.df, Class == "All genes" & omega != "NA" & omega < 100)$omega))

RT.avg.omega = mean(subset(omegaData.df, Class == "RT biased" & omega != "NA" & omega < 100)$omega)
RT.se.omega = sd(subset(omegaData.df, Class == "RT biased" & omega != "NA" & omega < 100)$omega)/sqrt(length(subset(omegaData.df, Class == "RT biased" & omega != "NA" & omega < 100)$omega))

SFP.avg.omega = mean(subset(omegaData.df, Class == "SFPs" & omega != "NA" & omega < 100)$omega)
SFP.se.omega = sd(subset(omegaData.df, Class == "SFPs" & omega != "NA" & omega < 100)$omega)/sqrt(length(subset(omegaData.df, Class == "SFPs" & omega != "NA" & omega < 100)$omega))

meanOmega.df=data.frame(Class = c("All", "RT biased", "SFPs"), omega=c(genome.avg.omega, RT.avg.omega, SFP.avg.omega), se = c(genome.se.omega, RT.se.omega, SFP.se.omega))

meanOmega.df$Class = factor (meanOmega.df$Class, levels = c("All", "RT biased", "SFPs"))

ggplot(meanOmega.df, aes(Class, omega, colour=Class)) + 
    geom_point(size = 2) + 
    geom_errorbar(aes(ymin=omega-se, ymax=omega+se), width=.2, position=position_dodge(.9)) +
    theme_bw() + 
    theme(legend.position="none") +
    ylab(expression(omega)) + 
    xlab("Gene class") +
    scale_colour_manual(values = c("black", "#7d49c3", "#4f922a"))

RT.paml.data = subset(paml.data, FBgn_ID %in% Dnov.dvir1.06.RT.list)
RT.paml.data.sig = subset(RT.paml.data, Damr_FDR < 0.05 | Dlum_FDR < 0.05 | Dnov_FDR < 0.05 | Dvir_FDR < 0.05)
RT.paml.data.sig.list = as.character(unique(RT.paml.data.sig$FBgn_ID))

#################################################
#################################################
#################################################

# create gene lists and factor labeling
con.vs.vir.Up = as.data.frame(con.vs.vir.RT.con.Up.list)
con.vs.vir.Up$V1 = "con.vs.vir.Up"
rownames(con.vs.vir.Up) = con.vs.vir.RT.con.Up.list
con.vs.vir.Up = subset(con.vs.vir.Up, select = "V1")

het.vs.vir.Up = as.data.frame(het.vs.vir.RT.het.Up.list)
het.vs.vir.Up$V1 = "het.vs.vir.Up"
rownames(het.vs.vir.Up) = het.vs.vir.RT.het.Up.list
het.vs.vir.Up = subset(het.vs.vir.Up, select = "V1")

con.vs.vir.Down = as.data.frame(con.vs.vir.RT.con.Down.list)
con.vs.vir.Down$V1 = "con.vs.vir.Down"
rownames(con.vs.vir.Down) = con.vs.vir.RT.con.Down.list
con.vs.vir.Down = subset(con.vs.vir.Down, select = "V1")

het.vs.vir.Down = as.data.frame(con.vs.het.RT.het.Up.list)
het.vs.vir.Down$V1 = "het.vs.vir.Down"
rownames(het.vs.vir.Down) = het.vs.vir.RT.het.Down.list
het.vs.vir.Down = subset(het.vs.vir.Down, select = "V1")

factor_labeling = rbind(het.vs.vir.Down)
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


ggplot(subset(GO_enrichment_data, over_represented_FDR < 0.05), 
       aes(category, -log10(over_represented_pvalue), size = numDEInCat, colour = ontology)) + 
    geom_point()  + 
    facet_grid(~factor) +
    xlab(NULL) + 
    geom_text_repel(data = subset(GO_enrichment_data, over_represented_FDR < 0.05 & numDEInCat > 5), 
                    aes(category, -log10(over_represented_pvalue),label=term), 
                    force = 8, 
                    inherit.aes = F, 
                    box.padding = unit(0.35, "lines"), 
                    point.padding = unit(0.5, "lines"), 
                    fontface = "bold", 
                    size = 3) + 
    theme(axis.text.x = element_text(angle = 45, face = "bold", vjust = 1, hjust = 1)) + 
    scale_size(range = c(0,12)) + 
    scale_colour_manual(values=c("#c27d92", "#a8a34b", "#8e61bf")) + 
    scale_y_continuous(limits=c(3, 10))





plotBoth<-function(gene_id){
    male = plotGeneG(TPMse, gene_id) + ggtitle("Male")
    female = plotGenePM(TPMse_DnovPM, gene_id) + ggtitle("Female (D. nov)") + theme(axis.text.x = element_text(angle=45))
    both = plot_grid(female, male, ncol = 2)
    return(both)
}


###############


TSmaleDerived.hetSig = intersect(het.vs.vir.RT.het.Up.list, unlist(TS_elements))
TSmaleDerived.conSig = intersect(con.vs.vir.RT.con.Up.list, unlist(TS_elements))
TS.DE.male.Derived = union(TSmaleDerived.conSig, TSmaleDerived.hetSig)

AGmaleDerived.conSig = intersect(con.vs.vir.RT.con.Up.list, unlist(AG_elements))
AGmaleDerived.hetSig = intersect(het.vs.vir.RT.het.Up.list, unlist(AG_elements))
AG.DE.male.Derived = union(AGmaleDerived.conSig, AGmaleDerived.hetSig)

pdf("TS.DE.male.Derived.pdf", height = 3, width = 9)
lapply(TS.DE.male.Derived, plotBoth)
dev.off()

pdf("AG.DE.male.Derived.pdf", height = 3, width = 9)
lapply(AG.DE.male.Derived, plotBoth)
dev.off()