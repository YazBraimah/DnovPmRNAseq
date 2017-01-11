## This script loads all the required packages and several custom plotting functions

## Define required packages
req_packages = c("Biobase","cluster","cowplot","cummeRbund","data.table","DESeq","edgeR","ggplot2","ggrepel","ggthemes","GO.db","goseq","grid","gridExtra","plotly","qvalue","reshape", "Rmisc","splitstackshape","statmod","VennDiagram", "ggthemr")

## Load them
lapply(req_packages, require, character.only = TRUE)

## Set ggplot's theme back to defaults (cowplot changes it)
theme_set(theme_gray())

## Load heatmap script and miscellaneous R functions (From the Trinity package)
source("~/Dropbox/RNAseq/Male_RNAseq/VirilisGroupMaleRNAseq.Rproject/R.functions/heatmap.3.R")
source("~/Dropbox/RNAseq/Male_RNAseq/VirilisGroupMaleRNAseq.Rproject/R.functions/misc_rnaseq_funcs.R")

### Here're several convenient plotting functions:

## MA plots between two columns of a matrix. Also calculates the proportion of genes that are >2-fold against the logCount.
MA_BPlot <- function(data, col1, col2) {
    arguments <- as.list(match.call())
    y = eval(as.name(arguments$col2), data)
    x = eval(as.name(arguments$col1), data)
    M = log2(x/y)
    A = 0.5*log2(x*y);
    res = data.frame(row.names(data),x=M, y=A, row.names = 1)
    res$bins <- cut(x=res$y, breaks=seq(from=0, to=20, by = 0.5), labels=as.character(seq(0.5,20,0.5)))
    res$count = 1
    bad.res = subset(res, x >= 1 | x <= -1)
    bad.res = subset(bad.res, x != Inf & x != (-Inf))
    badBygroup = tapply(bad.res$count, bad.res$bins, sum)
    allBygroup = tapply(res$count, res$bins, sum)
    off.2fold = badBygroup/allBygroup
    off.2fold = data.frame(off.2fold)
    off.2fold <- cbind(log2Bin = row.names(off.2fold), off.2fold)
    rownames(off.2fold) <- NULL
    off.2fold$log2Bin <- factor(off.2fold$log2Bin, levels=as.character(seq(0.5,20,0.5)))
    b.plot <- ggplot(off.2fold, aes(log2Bin, off.2fold)) + geom_bar(stat = "identity", fill="#F8766D", colour = "#00BFC4") + labs(y="prop. >2-fold", x="log2 read count bin") + theme(legend.position="none") + labs(title = paste(col1, "vs", col2, sep = " "), face = "bold")
    res$FC = NA
    res$FC[res$x > 1 ] = "above" 
    res$FC[res$x < -1 ] = "above"
    ma.plot = ggplot(res, aes(y, x, colour = FC)) + geom_point() + labs(y = "M [log2(x/y)]", x = "A [0.5xlog2(xy)]")+ scale_x_continuous(limits=c(0,20)) + theme(legend.position="none")
    #ma.plot = qplot(y, x, data = res, colour = FC, ylab = "M [log2(x/y)]", xlab = "A [0.5xlog2(xy)]") + scale_x_continuous(limits=c(0,20)) + theme(legend.position="none")
    plots = plot_grid(b.plot, ma.plot, ncol = 1)
    return(plots)
}
#
plotGenePM<-function(object, gene_id, logMode=FALSE){
  if (grepl("FBgn", gene_id)){
    geneName<-subset(Annots, FBgn_ID == gene_id)$gene_name
  } else {geneName<-subset(Annots, gene_name == gene_id)$FBgn_ID}
  swisprotName<-subset(Annots, FBgn_ID == gene_id | gene_name == gene_id)$SwissProt_BlastX_Description
  melOrth<-subset(Annots, FBgn_ID == gene_id | gene_name == gene_id)$mel_GeneSymbol
  omega<-subset(paml.data, FBgn_ID == gene_id | gene_name == gene_id)$omega
  coords.tmp<-subset(gffRecord, FBgn_ID == gene_id | gene_name == gene_id)
  coords<-paste(coords.tmp$chromosome, ":", coords.tmp$min,"-",coords.tmp$max, sep = "")
  p <- ggplot(subset(object, FBgn_ID == gene_id | gene_name == gene_id), aes(time, TPM, fill = condition)) + geom_bar(position=position_dodge(), stat="identity") + geom_errorbar(aes(ymin=TPM-se, ymax=TPM+se), width=.2, position=position_dodge(.9)) + facet_grid(~tissue, scales="free_x", space = "free_x") + labs(title = paste(gene_id," (", geneName,"), ",coords,"\n","Ka/Ks = ", omega,"        mel. orth.: ",melOrth,"\n",swisprotName, sep = "")) + scale_fill_manual(values = c("#77756b", "#b067a3", "#9c954d"))
  if (logMode)
  {
    p <- p + scale_y_log10()
  }
  if (logMode)
  {
    p <- p + ylab("log10 TPM")
  } else {
    p <- p + ylab("TPM")
  }
  return(p)
}
#
plotGenePM_RT<-function(object, gene_id, logMode=FALSE){
  if (grepl("FBgn", gene_id)){
    geneName<-subset(Annots, FBgn_ID == gene_id)$gene_name
  } else {geneName<-subset(Annots, gene_name == gene_id)$FBgn_ID}
  swisprotName<-subset(Annots, FBgn_ID == gene_id | gene_name == gene_id)$SwissProt_BlastX_Description
  melOrth<-subset(Annots, FBgn_ID == gene_id | gene_name == gene_id)$mel_GeneSymbol
  omega<-subset(paml.data, FBgn_ID == gene_id | gene_name == gene_id)$omega
  coords.tmp<-subset(gffRecord, FBgn_ID == gene_id | gene_name == gene_id)
  coords<-paste(coords.tmp$chromosome, ":", coords.tmp$min,"-",coords.tmp$max, sep = "")
  p <- ggplot(subset(object, FBgn_ID == gene_id & tissue == "RT"| gene_name == gene_id & tissue == "RT"), aes(time, TPM, fill = condition)) + geom_bar(position=position_dodge(), stat="identity") + geom_errorbar(aes(ymin=TPM-se, ymax=TPM+se), width=.2, position=position_dodge(.9)) + labs(title = paste(gene_id," (", geneName,"), ",coords,"\n","Ka/Ks = ", omega,"        mel. orth.: ",melOrth,"\n",swisprotName, sep = "")) + scale_fill_manual(values = c("#77756b", "#b067a3", "#9c954d"))
  if (logMode)
  {
    p <- p + scale_y_log10()
  }
  if (logMode)
  {
    p <- p + ylab("log10 TPM")
  } else {
    p <- p + ylab("TPM")
  }
  return(p)
}
# Plot mel gene
plotGeneMel<-function(object, gene_id, logMode=FALSE){
    if (grepl("FBgn", gene_id)){
        geneName<-subset(object, mel_FBgn_ID == gene_id)$mel_GeneSymbol
    } else {geneName<-subset(object, mel_GeneSymbol == gene_id)$mel_FBgn_ID}
    p <- ggplot(subset(object, mel_FBgn_ID == gene_id | mel_GeneSymbol == gene_id), aes(tissue, as.numeric(RPKM), fill = status)) + 
        geom_bar(stat = "identity", position = "dodge") + 
        facet_grid(age~sex, scales = "free") + 
        theme(axis.text.x=element_text(angle=45, vjust = 0.1)) + 
        labs(title = paste(gene_id,", ", geneName, sep = ""))
    if (logMode)
    {
        p <- p + scale_y_log10()
    }
    if (logMode)
    {
        p <- p + ylab("log10 RPKM")
    } else {
        p <- p + ylab("RPKM")
    }
    return(p)
}

## Heatmap plotter (based on cummerBund's heatmap function)
plotHeatmap <- function (object, rescaling = "row", clustering = "row", labCol = T, labRow = F, logMode = T, pseudocount = 1, border = FALSE, heatscale = c(low = "blue", mid = "black", high = "yellow"), heatMidpoint = NULL, method = "none", fullnames = T, replicates = FALSE, ...) {  
    
    m = object
    m = m[!apply(m, 1, sum) == 0, ]
    if (logMode) {
        m = log2(m + pseudocount)
    }
    
    if (!is.function(method)) {
        method = function(mat) {
            JSdist(makeprobs(t(mat)))
        }
    }
    if (clustering == "row") {
        m = m[hclust(method(m))$order, ]
    }
    if (clustering == "column") {
        m = m[, hclust(method(t(m)))$order]
    }
    if (clustering == "both") {
        m = m[hclust(method(m))$order, hclust(method(t(m)))$order]
    }
    if (is.function(rescaling)) {
        m = rescaling(m)
    } else {
        if (rescaling == "column") {
            m = scale(m, center = T)
            m[is.nan(m)] = 0
        }
        if (rescaling == "row") {
            m = t(scale(t(m), center = T))
            m[is.nan(m)] = 0
        }
    }
    rows = dim(m)[1]
    cols = dim(m)[2]
    melt.m = cbind(rowInd = rep(1:rows, times = cols), colInd = rep(1:cols, each = rows), melt(m))
    g = ggplot(data = melt.m)
    if (border == TRUE) 
        g2 = g + geom_rect(aes(xmin = colInd - 1, xmax = colInd, 
                               ymin = rowInd - 1, ymax = rowInd, fill = value), 
                           colour = "grey")
    if (border == FALSE) 
        g2 = g + geom_rect(aes(xmin = colInd - 1, xmax = colInd, 
                               ymin = rowInd - 1, ymax = rowInd, fill = value))
    if (labCol == T) {
        g2 = g2 + scale_x_continuous(breaks = (1:cols) - 
                                         0.5, labels = colnames(m))
    }
    if (labCol == F) {
        g2 = g2 + scale_x_continuous(breaks = (1:cols) - 
                                         0.5, labels = rep("", cols))
    }
    if (labRow == T) {
        g2 = g2 + scale_y_continuous(breaks = (1:rows) - 
                                         0.5, labels = rownames(m))
    }
    if (labRow == F) {
        g2 = g2 + scale_y_continuous(breaks = (1:rows) - 
                                         0.5, labels = rep("", rows))
    }
    g2 <- g2 + theme(axis.ticks = element_blank())
    g2 = g2 + theme(panel.grid.minor = element_line(colour = NA), 
                    panel.grid.major = element_line(colour = NA), panel.background = element_rect(fill = NA, 
                                                                                                  colour = NA))
    g2 = g2 + theme(axis.text.x = element_text(angle = -90, 
                                               hjust = 0))
    if (logMode & rescaling == "row") {
        legendTitle <- bquote(paste("centered ", log[2], " (TPM + ", 
                                    .(pseudocount), ")", sep = ""))
    }
    else if (logMode & rescaling == "none") {
        legendTitle <- bquote(paste(log[2], " (TPM + ", 
                                    .(pseudocount), ")",sep = ""))
    }
    else {
        legendTitle <- "TPM"
    }
    if (length(heatscale) == 2) {
        g2 <- g2 + scale_fill_gradient(low = heatscale[1], 
                                       high = heatscale[2], name = legendTitle)
    }
    else if (length(heatscale) == 3) {
        if (is.null(heatMidpoint)) {
            heatMidpoint = (max(m) + min(m))/2
        }
        g2 <- g2 + scale_fill_gradient2(low = heatscale[1], 
                                        mid = heatscale[2], high = heatscale[3], midpoint = heatMidpoint, 
                                        name = legendTitle)
    }
    return(g2)
}  

## A function to calculate the tissue specificity index (based on CummerBund's S function)
calcSpecificity<-function(matrix,logMode=T,pseudocount=1,relative=FALSE){
    tpms<-matrix
    if(logMode){
        tpms<-log10(tpms+pseudocount)
    }
    tpms<-t(makeprobs(t(tpms)))
    d<-diag(ncol(tpms))
    res<-apply(d,MARGIN=1,function(q){
        JSdistFromP(tpms,q)
    })
    colnames(res)<-paste(colnames(tpms))
    
    if(relative){
        res<-res/max(res)
    }
    1-res
}

## Lookup annotation information of a given gene (dvir1.06 data)
geneLookupG <- function(gene, complete=F) {
    if (complete){
        result <- noquote(t(subset(Annots, FBgn_ID == gene | gene_name == gene)))
        return (result)
    } else {
        result <- noquote(t(subset(Annots, FBgn_ID == gene | gene_name == gene, select = c(gene_name, FBgn_ID, FBtr_ID, chromosome, SwissProt_BlastX_Description, mel_GeneSymbol))))
    }
    return (result) 
}

## Modifications of functions to compare groups of lists 
## (from http://stackoverflow.com/questions/23559371/how-to-get-the-list-of-items-in-venn-diagram-in-r)
Intersect <- function (x) {  
    # Multiple set version of intersect
    # x is a list
    if (length(x) == 1) {
        unlist(x)
    } else if (length(x) == 2) {
        intersect(x[[1]], x[[2]])
    } else if (length(x) > 2){
        intersect(x[[1]], Intersect(x[-1]))
    }
}
#
Union <- function (x) {  
    # Multiple set version of union
    # x is a list
    if (length(x) == 1) {
        unlist(x)
    } else if (length(x) == 2) {
        union(x[[1]], x[[2]])
    } else if (length(x) > 2) {
        union(x[[1]], Union(x[-1]))
    }
}
#
Setdiff <- function (x, y) {
    # Remove the union of the y's from the common x's. 
    # x and y are lists of characters.
    xx <- Intersect(x)
    yy <- Union(y)
    setdiff(xx, yy)
}

## Ouput the color IDs used by ggplot
gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
}

## Miscellaneous operators
'%!in%' <- function(x,y)!('%in%'(x,y))
