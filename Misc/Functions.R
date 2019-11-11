####
####
#### A few useful functions used in the analysis



##------------------------------------------------------------------------##
##------------------------------------------------------------------------##

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

##------------------------------------------------------------------------##
##------------------------------------------------------------------------##

setClusters <- function(matrix, k){
    df = matrix
    df = log2(df+1)
    df = as.matrix(df) # convert to matrix
    df = t(scale(t(df), scale=F))
    sample_cor = cor(df, method='pearson', use='pairwise.complete.obs')
    sample_dist = as.dist(1-sample_cor)
    hc_samples = hclust(sample_dist, method='complete')
    gene_cor = NULL
    if (is.null(gene_cor)) { gene_cor = cor(t(df), method='pearson', use='pairwise.complete.obs') }
    gene_dist = as.dist(1-gene_cor)
    hc_genes = hclust(gene_dist, method='complete')
    gene_partition_assignments <- cutree(as.hclust(hc_genes), k=k)
    
    cluster_list = list()
    max_cluster_count = max(gene_partition_assignments)
    for (i in 1:max_cluster_count) {
        partition_i = (gene_partition_assignments == i)
        partition_data = df[partition_i,,drop=F]
        cluster_list[[i]] = as.data.frame(partition_data)
        cluster_list[[i]]$Gene_ID = rownames(cluster_list[[i]])
        cluster_list[[i]]$cluster = paste("cluster_", i, sep = "")
    }

    #cluster_list
    cluster_df = rbindlist(cluster_list)
    return(cluster_df)
}

##------------------------------------------------------------------------##
##------------------------------------------------------------------------##

geneBoxPlot <- function (object, gene) {
        if (grepl("FBgn", gene)) {
            swisprotName <- subset(annot.sum, FBgn_ID == gene)$SwissProt_BlastX_Description
            melOrth <- subset(melOrthsAll, FBgn_ID == gene)$mel_GeneSymbol
            geneName <- subset(annot.sum, FBgn_ID == gene)$gene_name
            
            p <- ggplot(subset(object, FBgn_ID == gene), aes(Time, TPM, colour = Status)) + 
            geom_boxplot(position = "dodge", outlier.size = 0, width = 0.3) + 
            geom_point(pch = 21, position = position_jitterdodge(), aes(colour = Status)) + 
            facet_grid(~Tissue, scales = "free_x", space = "free_x") + 
            labs(title = paste(gene, " (", geneName, ")", sep = ""), subtitle = paste("mel. orth.: ", melOrth, "\n", swisprotName, sep = "")) + 
            theme_bw() + 
            scale_colour_manual(values = c("#3f5a2a","#ffb200","#00b5d4"))
            
        } else if (grepl("MSTRG", gene)) {
            swisprotName <- subset(Annots.g, gene_id == gene)$sprot_Top_BLASTX_hit_description
            melOrth <- subset(Annots.g, gene_id == gene)$mel_GeneSymbol
            geneName <- subset(Annots.g, gene_id == gene)$gene_name
            
            p <- ggplot(subset(object, gene_id == gene), aes(Time, TPM, colour = Status)) + 
            geom_boxplot(position = "dodge", outlier.size = 0, width = 0.3) + 
            geom_point(pch = 21, position = position_jitterdodge(), aes(colour = Status)) + 
            facet_grid(~Tissue, scales = "free_x", space = "free_x") + 
            labs(title = paste(gene, " (", geneName, ")", sep = ""), subtitle = paste("mel. orth.: ", melOrth, "\n", swisprotName, sep = "")) + 
            theme_bw() + 
            scale_colour_manual(values = c("#3f5a2a","#ffb200","#00b5d4"))
        } else if (grepl("gene", gene) | grepl("TRINITY", gene)){
            swisprotName <- subset(tTrinotate.sub.some, gene_id.y == gene)$sprot_Top_BLASTX_hit_description
            geneName <- subset(tTrinotate.sub.some, gene_id.y == gene)$gene_name
            
            p <- ggplot(subset(object, gene_id == gene), aes(Time, TPM, colour = Status)) + 
            geom_boxplot(position = "dodge", outlier.size = 0, width = 0.3) + 
            geom_point(pch = 21, position = position_jitterdodge(), aes(colour = Status)) + 
            facet_grid(~Tissue, scales = "free_x", space = "free_x") + 
            labs(title = paste(gene, " (", geneName, ")", sep = ""), subtitle = paste(swisprotName, sep = "")) + 
            theme_bw() + 
            scale_colour_manual(values = c("#3f5a2a","#ffb200","#00b5d4"))
            
        } else {
            fbgn_id <- subset(annot.sum, gene_name == gene)$FBgn_ID
            melOrth <- subset(Annots, gene_name == gene)$mel_GeneSymbol
            swisprotName <- subset(Annots, FBgn_ID == fbgn_id)$sprot_Top_BLASTX_hit_description
            
            p <- ggplot(subset(object, gene_name == gene), aes(Time, TPM, colour = Status)) + 
            geom_boxplot(position = "dodge", outlier.size = 0, width = 0.3) + 
            geom_point(pch = 21, position = position_jitterdodge(), aes(colour = Status)) + 
            facet_grid(~Tissue, scales = "free_x", space = "free_x") + 
            labs(title = paste(gene, " (", fbgn_id, ")", sep = ""), subtitle = paste("mel. orth.: ", melOrth, "\n", swisprotName, sep = "")) + 
            theme_bw() + 
            scale_colour_manual(values = c("#3f5a2a","#ffb200","#00b5d4"))
        }
        
    return(p)
}

##------------------------------------------------------------------------##
##------------------------------------------------------------------------##

heatmap_mean <- function(tpmMatrix, gene_list, title, x = 2, melNames = FALSE, show_names = FALSE) {

    col_annot = unique(dplyr::select(sampleInfo, SampleName, Sex, Tissue, Time))
    rownames(col_annot) = col_annot$SampleName
    col_annot = subset(col_annot, select = c("Sex", "Tissue"))

    data <- subset(tpmMatrix, rownames(tpmMatrix) %in% gene_list)
    
## process tmm matrix
data = log2(data+1)
data = as.data.frame(t(scale(t(data), scale=F)))
data[data < -x] = -x
data[data > x] = x
    if(melNames) {
        init_cols = colnames(data)
        data$FBgn_ID = rownames(data)
        data = merge(data, melOrthsAll_dupsMarked, by.x = "FBgn_ID", by.y = "FBgn_ID", all.x = T)
        data = merge(data, fbgn_to_geneName, by.x = "FBgn_ID", by.y = "FBgn_ID", all.x = T)
        data$mel_GeneSymbol = ifelse(is.na(data$mel_GeneSymbol), paste("unk. (", data$gene_name, ")", sep = ""), data$mel_GeneSymbol)
        rownames(data) = data$mel_GeneSymbol
        data = subset(data, select = init_cols)
    } else {
        init_cols = colnames(data)
        data$FBgn_ID = rownames(data)
        data = merge(data, fbgn_to_geneName, by.x = "FBgn_ID", by.y = "FBgn_ID", all.x = T)
        rownames(data) = data$gene_name
        data = subset(data, select = init_cols)
    }

    if(show_names){
        p <- pheatmap(mat = data, main = title, color = inferno(100), border_color = NA, show_colnames = TRUE, show_rownames = TRUE, annotation_col = col_annot, drop_levels = TRUE,annotation_names_row = F, fontsize = 8)  

    } else {
        p <- pheatmap(mat = data, main = title, color = inferno(100), border_color = NA, show_colnames = TRUE, show_rownames = FALSE, annotation_col = col_annot, drop_levels = TRUE,annotation_names_row = F, fontsize = 8)
    }

    return(p)
    
}

##------------------------------------------------------------------------##
##------------------------------------------------------------------------##

heatmap_mean_ra <- function (tpmMatrix, gene_list, title, x = 2, 
    show_names = FALSE, row_annots) 
{
    col_annot = unique(dplyr::select(sampleInfo, SampleName, Sex, Tissue, 
        Time))
    rownames(col_annot) = col_annot$SampleName
    col_annot = subset(col_annot, select = c("Sex", "Tissue"))
    data <- subset(tpmMatrix, rownames(tpmMatrix) %in% gene_list)
    data = log2(data + 1)
    data = as.data.frame(t(scale(t(data), scale = F)))
    data[data < -x] = -x
    data[data > x] = x
#     init_cols = colnames(data)
#         data$FBgn_ID = rownames(data)
#         data = merge(data, fbgn_to_geneName, by.x = "FBgn_ID", 
#             by.y = "FBgn_ID", all.x = T)
#         rownames(data) = data$gene_name
#         data = subset(data, select = init_cols)
    if (show_names) {
        p <- pheatmap(mat = data, main = title, color = viridis(100), 
            border_color = NA, show_colnames = TRUE, show_rownames = TRUE, 
            annotation_col = col_annot, annotation_row = row_annots, drop_levels = TRUE, annotation_names_row = T, 
            fontsize = 8)
    }
    else {
        p <- pheatmap(mat = data, main = title, color = viridis(100), 
            border_color = NA, show_colnames = TRUE, show_rownames = FALSE, 
            annotation_col = col_annot, annotation_row = row_annots, drop_levels = TRUE, annotation_names_row = T, 
            fontsize = 8)
    }
    return(p)
}

##------------------------------------------------------------------------##
##------------------------------------------------------------------------##

## extract gene IDs based on GO term and enriched gene list:
extract_GO_genes = function(go_term, gene_set, trinity = FALSE){
    if (trinity){
        rownames(subset(GOinfo_pasa, row.names(GOinfo_pasa) %in% gene_set & grepl(go_term, GOinfo_pasa$V2)))
    } else {
    rownames(subset(GOinfo_annotated, row.names(GOinfo_annotated) %in% gene_set & grepl(go_term, GOinfo_annotated$V2)))
    }
}

##------------------------------------------------------------------------##
##------------------------------------------------------------------------##

corr_eqn <- function(x,y, digits = 2) {
  corr_coef <- round(cor(x, y, use = "pairwise.complete.obs"), digits = digits)
  paste("italic(R)^2 == ", corr_coef)
}

##------------------------------------------------------------------------##
##------------------------------------------------------------------------##


RT.geneBoxPlot <- function(object, gene) {
    swisprotName <- subset(annot.sum, FBgn_ID == gene)$SwissProt_BlastX_Description
        melOrth <- subset(melOrthsAll, FBgn_ID == gene)$mel_GeneSymbol
        geneName <- subset(annot.sum, FBgn_ID == gene)$gene_name
        p <- ggplot(subset(object, FBgn_ID == gene & Tissue == "repTract"), aes(Time, TPM, colour = Status)) + 
                geom_boxplot(position = "dodge", outlier.size = 0, width = 0.3) + 
                geom_point(pch = 21, position = position_jitterdodge(), aes(colour = Status)) +
#                 facet_grid(~Tissue, scales = "free_x", space = "free_x") +
                labs(title = paste(geneName, " (", melOrth, ")", sep = ""), subtitle = paste(swisprotName,  sep = "")) + 
                theme_minimal() + 
                theme(legend.position="none", axis.text.x = element_text(angle = 45, hjust = 1, size = 12), axis.text.y = element_text(size = 12), axis.title.x = element_blank()) +
                scale_colour_manual(values = c("#3f5a2a", "#ffb200", "#00b5d4"))
    return(p)
}

##----------------------------------------##
##----------------------------------------##

OV.geneBoxPlot <- function(object, gene) {
    swisprotName <- subset(annot.sum, FBgn_ID == gene)$SwissProt_BlastX_Description
        melOrth <- subset(melOrthsAll, FBgn_ID == gene)$mel_GeneSymbol
        geneName <- subset(annot.sum, FBgn_ID == gene)$gene_name
        p <- ggplot(subset(object, FBgn_ID == gene & Tissue == "ovaries"), aes(Time, TPM, colour = Status)) + 
                geom_boxplot(position = "dodge", outlier.size = 0, width = 0.3) + 
                geom_point(pch = 21, position = position_jitterdodge(), aes(colour = Status)) +
#                 facet_grid(~Tissue, scales = "free_x", space = "free_x") +
                labs(title = paste(geneName, " (", melOrth, ")", sep = ""), subtitle = paste(swisprotName,  sep = "")) + 
                theme_minimal() + 
                theme(legend.position="none", axis.text.x = element_text(angle = 45, hjust = 1, size = 12), axis.text.y = element_text(size = 12), axis.title.x = element_blank()) +
                scale_colour_manual(values = c("#3f5a2a", "#ffb200", "#00b5d4"))
    return(p)
}

##----------------------------------------##
##----------------------------------------##

H.geneBoxPlot <- function(object, gene) {
    swisprotName <- subset(annot.sum, FBgn_ID == gene)$SwissProt_BlastX_Description
        melOrth <- subset(melOrthsAll, FBgn_ID == gene)$mel_GeneSymbol
        geneName <- subset(annot.sum, FBgn_ID == gene)$gene_name
        p <- ggplot(subset(object, FBgn_ID == gene & Tissue == "head"), aes(Time, TPM, colour = Status)) + 
                geom_boxplot(position = "dodge", outlier.size = 0, width = 0.3) + 
                geom_point(pch = 21, position = position_jitterdodge(), aes(colour = Status)) +
#                 facet_grid(~Tissue, scales = "free_x", space = "free_x") +
                labs(title = paste(geneName, " (", melOrth, ")", sep = ""), subtitle = paste(swisprotName,  sep = "")) + 
                theme_minimal() + 
                theme(legend.position="none", axis.text.x = element_text(angle = 45, hjust = 1, size = 12), axis.text.y = element_text(size = 12), axis.title.x = element_blank()) +
                scale_colour_manual(values = c("#3f5a2a", "#ffb200", "#00b5d4"))
    return(p)
}
##------------------------------------------------------------------------##
##------------------------------------------------------------------------##


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

##------------------------------------------------------------------------##
##------------------------------------------------------------------------##

## estimate number of cluster for K-means
findK<-function(object, k.range=c(2:20), logMode=T, pseudocount=1,...){
    require(cluster)
    m<-as.data.frame(object)
    m<-m[rowSums(m)>0,]
    if(logMode){
        m<-log10(m+pseudocount)
    }
    n<-JSdist(makeprobs(t(m)))
    myWidths<-c()
    for (k in k.range){
        #print(k)
        myWidths<-c(myWidths,pam(n,k,...)$silinfo$avg.width)
    }
    plot(k.range,myWidths)
}

##------------------------------------------------------------------------##
##------------------------------------------------------------------------##

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

##------------------------------------------------------------------------##
##------------------------------------------------------------------------##

## Ouput the color IDs used by ggplot
gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
}

##------------------------------------------------------------------------##
##------------------------------------------------------------------------##

# TPM function
tpm <- function(counts, lengths) {
    rate <- counts / lengths
    rate / sum(rate) * 1e6
}

##------------------------------------------------------------------------##
##------------------------------------------------------------------------##

## Miscellaneous operators
'%!in%' <- function(x,y)!('%in%'(x,y))

##------------------------------------------------------------------------##
##------------------------------------------------------------------------##

## pulled from here, and then tweaked slightly: http://www.biostars.org/p/18211/
 
# CODE
 
heatmap.3 <- function(x,
                      Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE,
                      distfun = dist,
                      hclustfun = hclust,
                      dendrogram = c("both","row", "column", "none"),
                      symm = FALSE,
                      scale = c("none","row", "column"),
                      na.rm = TRUE,
                      revC = identical(Colv,"Rowv"),
                      add.expr,
                      breaks,
                      symbreaks = max(x < 0, na.rm = TRUE) || scale != "none",
                      col = "heat.colors",
                      colsep,
                      rowsep,
                      sepcolor = "white",
                      sepwidth = c(0.05, 0.05),
                      cellnote,
                      notecex = 1,
                      notecol = "cyan",
                      na.color = par("bg"),
                      trace = c("none", "column","row", "both"),
                      tracecol = "cyan",
                      hline = median(breaks),
                      vline = median(breaks),
                      linecol = tracecol,
                      margins = c(5,5),
                      ColSideColors,
                      RowSideColors,
                      side.height.fraction=0.1,
                      #cexRow = 0.2 + 1/log10(max(nr,2)),
                      #cexCol = 0.2 + 1/log10(max(nc,2)),
        cexRow = 0.2,
        cexCol = 0.2,                  

        scaleRangeMin,
        scaleRangeMax,


    cex.main = 1,
                      labRow = NULL,
                      labCol = NULL,
                      key = TRUE,
                      keysize = 1.5,
                      density.info = c("none", "histogram", "density"),
                      denscol = tracecol,
                      symkey = max(x < 0, na.rm = TRUE) || symbreaks,
                      densadj = 0.25,
                      main = NULL,
                      xlab = NULL,
                      ylab = NULL,
                      lmat = NULL,
                      lhei = NULL,
                      lwid = NULL,
                      NumColSideColors = 1,
                      NumRowSideColors = 1,
                      KeyValueName="Value",...){
 
    invalid <- function (x) {
      if (missing(x) || is.null(x) || length(x) == 0)
          return(TRUE)
      if (is.list(x))
          return(all(sapply(x, invalid)))
      else if (is.vector(x))
          return(all(is.na(x)))
      else return(FALSE)
    }



    x <- as.matrix(x)
    scale01 <- function(x, low = min(x), high = max(x)) {
        x <- (x - low)/(high - low)
        x
    }

    retval <- list()


    scale <- if (symm && missing(scale))
        "none"
    else match.arg(scale)

    dendrogram <- match.arg(dendrogram)

    trace <- match.arg(trace)

    density.info <- match.arg(density.info)

    if (length(col) == 1 && is.character(col))
        col <- get(col, mode = "function")

    if (!missing(breaks) && (scale != "none"))
        warning("Using scale=\"row\" or scale=\"column\" when breaks are",
            "specified can produce unpredictable results.", "Please consider using only one or the other.")

    if (is.null(Rowv) || is.na(Rowv))
        Rowv <- FALSE

    if (is.null(Colv) || is.na(Colv))
        Colv <- FALSE
    else if (Colv == "Rowv" && !isTRUE(Rowv))
        Colv <- FALSE

    if (length(di <- dim(x)) != 2 || !is.numeric(x))
        stop("`x' must be a numeric matrix")

    nr <- di[1]
    nc <- di[2]

    if (nr <= 1 || nc <= 1)
        stop("`x' must have at least 2 rows and 2 columns")
    #print(paste("nr:", nr, "nc:", nc, "cexCol:", cexCol, "cexRow:", cexRow))
    #stop("debug")



    if (!is.numeric(margins) || length(margins) != 2)
        stop("`margins' must be a numeric vector of length 2")

    if (missing(cellnote))
        cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))

    if (!inherits(Rowv, "dendrogram")) {
        if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in% c("both", "row"))) {
            if (is.logical(Colv) && (Colv))
                dendrogram <- "column"
            else dedrogram <- "none"

            warning("Discrepancy: Rowv is FALSE, while dendrogram is `",
                dendrogram, "'. Omitting row dendogram.")
        }
    }

    if (!inherits(Colv, "dendrogram")) {
        if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in% c("both", "column"))) {
            if (is.logical(Rowv) && (Rowv))
                dendrogram <- "row"
            else dendrogram <- "none"

            warning("Discrepancy: Colv is FALSE, while dendrogram is `",
                dendrogram, "'. Omitting column dendogram.")
        }
    }
 
   if (inherits(Rowv, "dendrogram")) {
        ddr <- Rowv
        rowInd <- order.dendrogram(ddr)
    }
    else if (is.integer(Rowv)) {
        hcr <- hclustfun(distfun(x))
        ddr <- as.dendrogram(hcr)
        ddr <- reorder(ddr, Rowv)
        rowInd <- order.dendrogram(ddr)
        if (nr != length(rowInd))
            stop("row dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Rowv)) {
        Rowv <- rowMeans(x, na.rm = na.rm)
        hcr <- hclustfun(distfun(x))
        ddr <- as.dendrogram(hcr)
        ddr <- reorder(ddr, Rowv)
        rowInd <- order.dendrogram(ddr)
        if (nr != length(rowInd))
            stop("row dendrogram ordering gave index of wrong length")
    }
    else {
        rowInd <- nr:1
    }
 
   if (inherits(Colv, "dendrogram")) {
        ddc <- Colv
        colInd <- order.dendrogram(ddc)
    }
    else if (identical(Colv, "Rowv")) {
        if (nr != nc)
            stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
        if (exists("ddr")) {
            ddc <- ddr
            colInd <- order.dendrogram(ddc)
        }
        else colInd <- rowInd
    }
    else if (is.integer(Colv)) {
        hcc <- hclustfun(distfun(if (symm)
            x
        else t(x)))
        ddc <- as.dendrogram(hcc)
        ddc <- reorder(ddc, Colv)
        colInd <- order.dendrogram(ddc)
        if (nc != length(colInd))
            stop("column dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Colv)) {
        Colv <- colMeans(x, na.rm = na.rm)
        hcc <- hclustfun(distfun(if (symm)
            x
        else t(x)))
        ddc <- as.dendrogram(hcc)
        ddc <- reorder(ddc, Colv)
        colInd <- order.dendrogram(ddc)
        if (nc != length(colInd))
            stop("column dendrogram ordering gave index of wrong length")
    }
    else {
        colInd <- 1:nc
    }

    retval$rowInd <- rowInd
    retval$colInd <- colInd
    retval$call <- match.call()

    x <- x[rowInd, colInd]  # rearrange matrix according to dendrograms
    x.unscaled <- x

    cellnote <- cellnote[rowInd, colInd]  # also rearrange the cellnotes

    # get labels 
    if (is.null(labRow))
        labRow <- if (is.null(rownames(x)))
            (1:nr)[rowInd]
        else rownames(x)
    else labRow <- labRow[rowInd]
    if (is.null(labCol))
        labCol <- if (is.null(colnames(x)))
            (1:nc)[colInd]
        else colnames(x)
    else labCol <- labCol[colInd]


    ## do scaling of matrix according to Z-scores
    if (scale == "row") {
        retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
        x <- sweep(x, 1, rm)
        retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
        x <- sweep(x, 1, sx, "/")
    }
    else if (scale == "column") {
        retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
        x <- sweep(x, 2, rm)
        retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
        x <- sweep(x, 2, sx, "/")
    }

    # number of breaks
    if (missing(breaks) || is.null(breaks) || length(breaks) < 1) {
        if (missing(col) || is.function(col))
            breaks <- 16
        else breaks <- length(col) + 1
    }

    # set breakpoints
    if (length(breaks) == 1) {
        if (missing(scaleRangeMin))
            scaleRangeMin = min(x, na.rm=na.rm)

        if (missing(scaleRangeMax))
            scaleRangeMax = max(x, na.rm=na.rm)


        if (!symbreaks) {
            breaks <- seq(scaleRangeMin, scaleRangeMax, length=breaks);
        } else {
            #extreme <- max(abs(x), na.rm = TRUE)
            extreme = max(abs(c(scaleRangeMin,scaleRangeMax)), na.rm=na.rm)
            breaks <- seq(-extreme, extreme, length = breaks)
        }
    }

    nbr <- length(breaks)
    ncol <- length(breaks) - 1

    if (class(col) == "function")
        col <- col(ncol)

    min.breaks <- min(breaks)
    max.breaks <- max(breaks)

    # adjust for out-of-range given break settings
    x[x < min.breaks] <- min.breaks
    x[x > max.breaks] <- max.breaks

    # layout height
    if (missing(lhei) || is.null(lhei))
        lhei <- c(keysize, 4)

    # layout width
    if (missing(lwid) || is.null(lwid))
        lwid <- c(keysize, 4)

    # define the layout
    if (missing(lmat) || is.null(lmat)) {
        lmat <- rbind(4:3, 2:1)
 
        if (!missing(ColSideColors)) {
            if (!is.character(ColSideColors) || ncol(ColSideColors) != nc)
                stop("'ColSideColors' must be a matrix of ncol(x) ", nc, " columns")
            lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
            #lhei=c(lhei[1], side.height.fraction*NumColSideColors, lhei[2])
            side_height = min(side.height.fraction*nrow(ColSideColors), 1);
            lhei=c(lhei[1], side_height, lhei[2])
        }
 
        if (!missing(RowSideColors)) {
            if (!is.character(RowSideColors) || nrow(RowSideColors) != nr)
                stop("'RowSideColors' must be a matrix of nrow(x) ", nr, " rows.  It currently has ", nrow(RowSideColors), " rows.")
            lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 1), lmat[,2] + 1)
            #lwid <- c(lwid[1], side.height.fraction*NumRowSideColors, lwid[2])
            side_width = min(side.height.fraction*ncol(RowSideColors), 1);
            lwid <- c(lwid[1], side_width, lwid[2])
        }
        lmat[is.na(lmat)] <- 0
    }
 
    if (length(lhei) != nrow(lmat))
        stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
    if (length(lwid) != ncol(lmat))
        stop("lwid must have length = ncol(lmat) =", ncol(lmat))
    

    op <- par(no.readonly = TRUE)
    on.exit(par(op))
 
    layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
 
    ###########################################
    ## Draw the colorbars for the annotations:
    ########################################### 

    if (!missing(RowSideColors)) {
        if (!is.matrix(RowSideColors)){
                par(mar = c(margins[1], 0, 0, 0.5))
                image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
        } else {
            par(mar = c(margins[1], 0, 0, 0.5))
            rsc = t(RowSideColors[rowInd, , drop=F])
            rsc.colors = matrix()
            rsc.names = names(table(rsc))
            rsc.i = 1
            for (rsc.name in rsc.names) {
                rsc.colors[rsc.i] = rsc.name
                rsc[rsc == rsc.name] = rsc.i
                rsc.i = rsc.i + 1
            }
            # print(rsc)
            rsc = matrix(as.numeric(rsc), nrow = dim(rsc)[1])
            #print("RSC: ", rsc)
            #print(rsc.colors)    
            image(1:nrow(rsc), 1:ncol(rsc), rsc, col = as.vector(rsc.colors), axes = FALSE, xlab="", ylab="")
        
            # add labels
            if (length(colnames(RowSideColors)) > 0) {  
                #axis(1, 0:(dim(rsc)[2] - 1)/(dim(rsc)[2] - 1), rownames(RowSideColors), las = 2, tick = FALSE)
                #axis(1, 0:(nrow(rsc)-1), colnames(RowSideColors), las = 2, tick = T) # ncol because transposed
                axis(1, 1:ncol(RowSideColors), labels=colnames(RowSideColors), las=2, cex.axis=0.5, tick=F, xlab="", ylab="")

            }
        }
    }
    


    if (!missing(ColSideColors)) {
 
        if (!is.matrix(ColSideColors)){
            par(mar = c(0.5, 0, 0, margins[2]))
            image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
        } else {
            par(mar = c(0.5, 0, 0, margins[2]))
            csc = ColSideColors[, colInd, drop=F]
            csc.colors = matrix()
            csc.names = names(table(csc))
            csc.i = 1
            for (csc.name in csc.names) {
                csc.colors[csc.i] = csc.name
                csc[csc == csc.name] = csc.i
                csc.i = csc.i + 1
            }
            csc = matrix(as.numeric(csc), nrow = dim(csc)[1])
            #print(csc)
            image(1:nrow(t(csc)), 1:ncol(t(csc)), t(csc), col = as.vector(csc.colors), axes = FALSE, xlab="", ylab="")

            # add labels
            if (length(rownames(ColSideColors)) > 0) {
                #axis(2, 0:(dim(csc)[2] - 1)/max(1,(dim(csc)[2] - 1)), colnames(ColSideColors), las = 2, tick = FALSE)
                axis(2, 1:(nrow(ColSideColors)), labels=rownames(ColSideColors), las = 2, tick = FALSE, cex.axis=0.5)
            }
        }
    }
 


    par(mar = c(margins[1], 0, 0, margins[2]))
    x <- t(x)
    cellnote <- t(cellnote)
    if (revC) {
        iy <- nr:1
        if (exists("ddr"))
            ddr <- rev(ddr)
        x <- x[, iy]
        cellnote <- cellnote[, iy]
    }
    else iy <- 1:nr
    
    # draw the central heatmap
    image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, breaks = breaks, ...)
    
    # store the matrix drawn
    retval$carpet <- x
    
    # store the dendrograms
    if (exists("ddr"))
        retval$rowDendrogram <- ddr
    if (exists("ddc"))
        retval$colDendrogram <- ddc
    
    # store the breaks
    retval$breaks <- breaks
    
    # store the colormap used
    retval$col <- col
    
    # specially color in the na values  
    if (!invalid(na.color) & any(is.na(x))) { # load library(gplots)
        mmat <- ifelse(is.na(x), 1, NA)
        image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "", col = na.color, add = TRUE)
    }

    # X-axis column labels
    axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0, cex.axis = cexCol)

    # X-axis title
    if (!is.null(xlab))
        mtext(xlab, side = 1, line = margins[1] - 1.25)

    # Y-axis row labeling
    axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0,
        cex.axis = cexRow)

    # Y-axis title
    if (!is.null(ylab))
        mtext(ylab, side = 4, line = margins[2] - 1.25)
 
    if (!missing(add.expr))
        eval(substitute(add.expr))
    if (!missing(colsep))
        for (csep in colsep) rect(xleft = csep + 0.5, ybottom = rep(0, length(csep)), xright = csep + 0.5 + sepwidth[1], ytop = rep(ncol(x) + 1, csep), lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
    if (!missing(rowsep))
        for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) + 1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) + 1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
    

    min.scale <- min(breaks)
    max.scale <- max(breaks)
    x.scaled <- scale01(t(x), min.scale, max.scale)
    
    # column trace
    if (trace %in% c("both", "column")) {
        retval$vline <- vline
        vline.vals <- scale01(vline, min.scale, max.scale)
        for (i in colInd) {
            if (!is.null(vline)) {
                abline(v = i - 0.5 + vline.vals, col = linecol, lty = 2)
            }
            xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
            xv <- c(xv[1], xv)
            yv <- 1:length(xv) - 0.5
            lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
        }
    }

    # row trace
    if (trace %in% c("both", "row")) {
        retval$hline <- hline
        hline.vals <- scale01(hline, min.scale, max.scale)
        for (i in rowInd) {
            if (!is.null(hline)) {
                abline(h = i + hline, col = linecol, lty = 2)
            }
            yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
            yv <- rev(c(yv[1], yv))
            xv <- length(yv):1 - 0.5
            lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
        }
    }

    # add cell labels
    if (!missing(cellnote))
        text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote), col = notecol, cex = notecex)

    ###########################
    ## Plot the row dendrogram
    ###########################

    par(mar = c(margins[1], 0, 0, 0))
    if (dendrogram %in% c("both", "row")) {
        plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
    }
    else plot.new()

    #############################
    ## Plot the column dendrogram
    #############################

    par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
    if (dendrogram %in% c("both", "column")) {
        plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
    }
    else plot.new()

    if (!is.null(main))
        title(main, cex.main=cex.main) #cex.main = 1.5 * op[["cex.main"]])


    ############################
    ## Add the Color Chart
    ############################

    if (key) {
        par(mar = c(5, 4, 2, 1), cex = 0.75)
        tmpbreaks <- breaks
        if (symkey) {
            max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
            min.raw <- -max.raw
            tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
            tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)
        }
        else {
            min.raw <- min(c(x,breaks), na.rm = TRUE)
            max.raw <- max(c(x,breaks), na.rm = TRUE)
        }
 
        message('for plotting:: min.raw: ', min.raw, ' max.raw: ', max.raw);
        
        z <- seq(min.raw, max.raw, length = length(col))
        image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks,
            xaxt = "n", yaxt = "n")
        par(usr = c(0, 1, 0, 1))
        lv <- pretty(breaks)
        xv <- scale01(as.numeric(lv), min.raw, max.raw)
        axis(1, at = xv, labels = lv)
        if (scale == "row")
            mtext(side = 1, "Row Z-Score", line = 2)
        else if (scale == "column")
            mtext(side = 1, "Column Z-Score", line = 2)
        else mtext(side = 1, KeyValueName, line = 2)
        if (density.info == "density") {
            dens <- density(x, adjust = densadj, na.rm = TRUE)
            omit <- dens$x < min(breaks) | dens$x > max(breaks)
            dens$x <- dens$x[-omit]
            dens$y <- dens$y[-omit]
            dens$x <- scale01(dens$x, min.raw, max.raw)
            lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol,
                lwd = 1)
            axis(2, at = pretty(dens$y)/max(dens$y) * 0.95, pretty(dens$y))
            title("Color Key\nand Density Plot")
            par(cex = 0.5)
            mtext(side = 2, "Density", line = 2)
        }
        else if (density.info == "histogram") {
            h <- hist(x, plot = FALSE, breaks = breaks)
            hx <- scale01(breaks, min.raw, max.raw)
            hy <- c(h$counts, h$counts[length(h$counts)])
            lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s",
                col = denscol)
            axis(2, at = pretty(hy)/max(hy) * 0.95, pretty(hy))
            title("Color Key\nand Histogram")
            par(cex = 0.5)
            mtext(side = 2, "Count", line = 2)
        }
        else title("Color Key")
    }
    else plot.new()

    retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)], high = retval$breaks[-1], color = retval$col)

    invisible(retval)
}



# EXAMPLE USAGE
 
# example of colsidecolors rowsidecolors (single column, single row)
#mat <- matrix(1:100, byrow=T, nrow=10)
#column_annotation <- sample(c("red", "blue", "green"), 10, replace=T)
#column_annotation <- as.matrix(column_annotation)
#colnames(column_annotation) <- c("Variable X")
 
#row_annotation <- sample(c("red", "blue", "green"), 10, replace=T)
#row_annotation <- as.matrix(t(row_annotation))
#rownames(row_annotation) <- c("Variable Y")
 
#heatmap.3(mat, RowSideColors=row_annotation, ColSideColors=column_annotation)
 
# multiple column and row
#mat <- matrix(1:100, byrow=T, nrow=10)
#column_annotation <- matrix(sample(c("red", "blue", "green"), 20, replace=T), ncol=2)
#colnames(column_annotation) <- c("Variable X1", "Variable X2")
 
#row_annotation <- matrix(sample(c("red", "blue", "green"), 20, replace=T), nrow=2)
#rownames(row_annotation) <- c("Variable Y1", "Variable Y2")
 
#heatmap.3(mat, RowSideColors=row_annotation, ColSideColors=column_annotation)
 


#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------

# Borrowing the following from gplots  (gplots isn't compatible with R 3.0 (yet), and so bypassing it for now).

colorpanel = function (n, low, mid, high) 
{
    if (missing(mid) || missing(high)) {
        low <- col2rgb(low)
        if (missing(high)) 
            high <- col2rgb(mid)
        else high <- col2rgb(high)
        red <- seq(low[1, 1], high[1, 1], length = n)/255
        green <- seq(low[3, 1], high[3, 1], length = n)/255
        blue <- seq(low[2, 1], high[2, 1], length = n)/255
    }
    else {
        isodd <- odd(n)
        if (isodd) {
            n <- n + 1
        }
        low <- col2rgb(low)
        mid <- col2rgb(mid)
        high <- col2rgb(high)
        lower <- floor(n/2)
        upper <- n - lower
        red <- c(seq(low[1, 1], mid[1, 1], length = lower), seq(mid[1, 
            1], high[1, 1], length = upper))/255
        green <- c(seq(low[3, 1], mid[3, 1], length = lower), 
            seq(mid[3, 1], high[3, 1], length = upper))/255
        blue <- c(seq(low[2, 1], mid[2, 1], length = lower), 
            seq(mid[2, 1], high[2, 1], length = upper))/255
        if (isodd) {
            red <- red[-(lower + 1)]
            green <- green[-(lower + 1)]
            blue <- blue[-(lower + 1)]
        }
    }
    rgb(red, blue, green)
}


greenred = function (n)  {
    colorpanel(n, "green", "black", "red")
}

odd = function (x) {
    x%%2 == 1
}

even = function (x) {
    x%%2 == 0
}


plot_counts_matrix_log2_dist = function(matrix_file) {

    
    data = read.table(file=matrix_file, com='', row.names=1, header=T)

    conditions = colnames(data)
    colors = rainbow(length(conditions))


    plot(density(log2(data[,1])), col=colors[1], main=matrix_file, xlab='log2(frag_counts)', ylab='density')

    for (i in 2:length(data[1,])) {

        points(density(log2(data[,i])), type='l', col=colors[i])

    }

    legend('topright', conditions, col=colors, pch=15)

}


matrix_to_color_assignments = function(matrix_m, col=NULL, by=c("matrix", "row", "col")) {

    if (! is.matrix(matrix_m))
        stop("Error, matrix_to_color_assignments() requires a matrix as parameter.")
    num_colors = 0
    
    if (is.null(col)) {
        num_colors = min(nrow(matrix_m), ncol(matrix_m))
        col = rainbow(num_colors)
    }
    else {
        num_colors = length(col)
    }
    
    by = match.arg(by)
    
    if (by == "matrix") {

        min_val = min(matrix_m)
        matrix_m = matrix_m - min_val
        max_val = max(matrix_m)
        matrix_m = matrix_m / max_val * num_colors
        #print(matrix_m)
        matrix_m = apply(matrix_m, 1:2, function(x) ifelse (x<1, as.character(col[1]), as.character(col[x])));
        
        matrix_m = matrix(as.character(matrix_m), nrow=dim(matrix_m)[1])
    }
    else {

        row_or_col_only_color_selector_func = function(x) { 
                a = min(x); 
                b = max(x); 
                c = (x-a)/(b-a) * num_colors;
                c = round(c);
                c = ifelse (c<1, 1, c); 
                #print(paste(c("color selection: (a)", a, " (b)", b, " (c)", paste(c, sep=',')))); 
                colors = as.character(col[c]);
                return(colors);
        }
    
        if (by == "row") {
            matrix_m = t(apply(matrix_m, 1, row_or_col_only_color_selector_func));
        }
        else {
            # by column
            matrix_m = apply(matrix_m, 2, row_or_col_only_color_selector_func);
        }
    }
    
    #print(matrix_m)
    return(matrix_m)
}

sample_matrix_to_color_assignments = function(sampleAnnotationsMatrix, colors) {

    if (missing(colors))
        colors = rainbow(nrow(sampleAnnotationsMatrix))

    nsamples = nrow(sampleAnnotationsMatrix);

    if (length(colors) < nrow(sampleAnnotationsMatrix))
        stop("Error, only ", length(colors), " colors specified, but have ", nsamples, " samples");

    for (i in 1:nrow(sampleAnnotationsMatrix)) {
        c = colors[i]
        sampleAnnotationsMatrix[i,] = sapply(sampleAnnotationsMatrix[i,], function(x) ifelse( x, as.character(c), 'white'))
    }

    return(sampleAnnotationsMatrix);

}

##########
qqunif.plot<-function(pvalues, 
  should.thin=T, thin.obs.places=2, thin.exp.places=2, 
  xlab=expression(paste("Expected (",-log[10], " p-value)")),
  ylab=expression(paste("Observed (",-log[10], " p-value)")), 
  draw.conf=TRUE, conf.points=1000, conf.col="lightgray", conf.alpha=.05,
  already.transformed=FALSE, pch=20, aspect="iso", prepanel=prepanel.qqunif,
  par.settings=list(superpose.symbol=list(pch=pch)), ...) {
  
  
  #error checking
  if (length(pvalues)==0) stop("pvalue vector is empty, can't draw plot")
  if(!(class(pvalues)=="numeric" || 
    (class(pvalues)=="list" && all(sapply(pvalues, class)=="numeric"))))
    stop("pvalue vector is not numeric, can't draw plot")
  if (any(is.na(unlist(pvalues)))) stop("pvalue vector contains NA values, can't draw plot")
  if (already.transformed==FALSE) {
    if (any(unlist(pvalues)==0)) stop("pvalue vector contains zeros, can't draw plot")
  } else {
    if (any(unlist(pvalues)<0)) stop("-log10 pvalue vector contains negative values, can't draw plot")
  }
  
  
  grp<-NULL
  n<-1
  exp.x<-c()
  if(is.list(pvalues)) {
    nn<-sapply(pvalues, length)
    rs<-cumsum(nn)
    re<-rs-nn+1
    n<-min(nn)
    if (!is.null(names(pvalues))) {
      grp=factor(rep(names(pvalues), nn), levels=names(pvalues))
      names(pvalues)<-NULL
    } else {
      grp=factor(rep(1:length(pvalues), nn))
    }
    pvo<-pvalues
    pvalues<-numeric(sum(nn))
    exp.x<-numeric(sum(nn))
    for(i in 1:length(pvo)) {
      if (!already.transformed) {
        pvalues[rs[i]:re[i]] <- -log10(pvo[[i]])
        exp.x[rs[i]:re[i]] <- -log10((rank(pvo[[i]], ties.method="first")-.5)/nn[i])
      } else {
        pvalues[rs[i]:re[i]] <- pvo[[i]]
        exp.x[rs[i]:re[i]] <- -log10((nn[i]+1-rank(pvo[[i]], ties.method="first")-.5)/(nn[i]+1))
      }
    }
  } else {
    n <- length(pvalues)+1
    if (!already.transformed) {
      exp.x <- -log10((rank(pvalues, ties.method="first")-.5)/n)
      pvalues <- -log10(pvalues)
    } else {
      exp.x <- -log10((n-rank(pvalues, ties.method="first")-.5)/n)
    }
  }


  #this is a helper function to draw the confidence interval
  panel.qqconf<-function(n, conf.points=1000, conf.col="gray", conf.alpha=.05, ...) {
    require(grid)
    conf.points = min(conf.points, n-1);
    mpts<-matrix(nrow=conf.points*2, ncol=2)
          for(i in seq(from=1, to=conf.points)) {
                mpts[i,1]<- -log10((i-.5)/n)
                mpts[i,2]<- -log10(qbeta(1-conf.alpha/2, i, n-i))
                mpts[conf.points*2+1-i,1]<- -log10((i-.5)/n)
                mpts[conf.points*2+1-i,2]<- -log10(qbeta(conf.alpha/2, i, n-i))
          }
          grid.polygon(x=mpts[,1],y=mpts[,2], gp=gpar(fill=conf.col, lty=0), default.units="native")
      }

  #reduce number of points to plot
  if (should.thin==T) {
    if (!is.null(grp)) {
      thin <- unique(data.frame(pvalues = round(pvalues, thin.obs.places),
        exp.x = round(exp.x, thin.exp.places),
        grp=grp))
      grp = thin$grp
    } else {
      thin <- unique(data.frame(pvalues = round(pvalues, thin.obs.places),
        exp.x = round(exp.x, thin.exp.places)))
    }
    pvalues <- thin$pvalues
    exp.x <- thin$exp.x
  }
  gc()
  
  prepanel.qqunif= function(x,y,...) {
    A = list()
    A$xlim = range(x, y)*1.02
    A$xlim[1]=0
    A$ylim = A$xlim
    return(A)
  }

  #draw the plot
  xyplot(pvalues~exp.x, groups=grp, xlab=xlab, ylab=ylab, aspect=aspect,
    prepanel=prepanel, scales=list(axs="i"), pch=pch,
    panel = function(x, y, ...) {
      if (draw.conf) {
        panel.qqconf(n, conf.points=conf.points, 
          conf.col=conf.col, conf.alpha=conf.alpha)
      };
      panel.xyplot(x,y, ...);
      panel.abline(0,1);
    }, par.settings=par.settings, ...
  )
}

##----------------------------------------##
##----------------------------------------##

theme_black = function(base_size = 12, base_family = "") {
 
  theme_grey(base_size = base_size, base_family = base_family) %+replace%
 
    theme(
      # Specify axis options
      axis.line = element_blank(),  
      axis.text.x = element_text(size = base_size*0.8, color = "white", lineheight = 0.9),  
      axis.text.y = element_text(size = base_size*0.8, color = "white", lineheight = 0.9),  
      axis.ticks = element_line(color = "white", size  =  0.2),  
      axis.title.x = element_text(size = base_size, color = "white", margin = margin(0, 10, 0, 0)),  
      axis.title.y = element_text(size = base_size, color = "white", angle = 90, margin = margin(0, 10, 0, 0)),  
      axis.ticks.length = unit(0.3, "lines"),   
      # Specify legend options
      legend.background = element_rect(color = NA, fill = "black"),  
      legend.key = element_rect(color = "white",  fill = "black"),  
      legend.key.size = unit(1.2, "lines"),  
      legend.key.height = NULL,  
      legend.key.width = NULL,      
      legend.text = element_text(size = base_size*0.8, color = "white"),  
      legend.title = element_text(size = base_size*0.8, face = "bold", hjust = 0, color = "white"),  
      legend.position = "right",  
      legend.text.align = NULL,  
      legend.title.align = NULL,  
      legend.direction = "vertical",  
      legend.box = NULL, 
      # Specify panel options
      panel.background = element_rect(fill = "black", color  =  NA),  
      panel.border = element_rect(fill = NA, color = "white"),  
      panel.grid.major = element_line(color = "grey35"),  
      panel.grid.minor = element_line(color = "grey20"),  
      panel.margin = unit(0.5, "lines"),   
      # Specify facetting options
      strip.background = element_rect(fill = "grey30", color = "grey10"),  
      strip.text.x = element_text(size = base_size*0.8, color = "white"),  
      strip.text.y = element_text(size = base_size*0.8, color = "white",angle = -90),  
      # Specify plot options
      plot.background = element_rect(color = "black", fill = "black"),  
      plot.title = element_text(size = base_size*1.2, color = "white"),  
      plot.margin = unit(rep(1, 4), "lines")
 
    )
 
}

##----------------------------------------##
##----------------------------------------##

theme_monokai_full <- function(base_size = 14, base_family = ""){
  color.background = "#232323"
  color.grid.major = "#232323"
  color.text = "#ffffff"
  color.axis = "#ffffff"

  theme_bw(base_size=base_size) +
    theme(

      panel.background=element_rect(fill=color.background, color=NA),
      plot.background=element_rect(fill=color.background, color=color.background),
      panel.border=element_rect(color=color.axis),

      panel.grid.major=element_line(color="grey95",size=.1, linetype=3),
      panel.grid.minor=element_blank(),
#       axis.line.x=element_line(color=color.grid.major, size=1),
#       axis.line.y=element_line(color=color.grid.major, size=1),
#       axis.ticks=element_line(color=NA),

      legend.background = element_rect(fill=color.background),
      legend.key = element_rect(fill=color.background, color=NA),
      legend.text = element_text(size=rel(.8),color=color.text),#color.axis.title),
      legend.title = element_text(color=color.text),

      plot.title=element_text(color=color.text, size=rel(1.2)),
      plot.subtitle=element_text(color=color.text, size=rel(0.8)),
      axis.text.x=element_text(size=rel(.95),color=color.text,angle = 45, hjust = 1),
      axis.text.y=element_text(size=rel(.95),color=color.text),
      axis.title.x=element_text(size=rel(1),color=color.text, vjust=0),
      axis.title.y=element_text(size=rel(1),color=color.text, vjust=1.25),
      strip.background = element_rect(fill = color.background, colour = color.axis),
      strip.text = element_text(colour = color.text)
    )

}

##----------------------------------------##
##----------------------------------------##

theme_black_full <- function(base_size = 14, base_family = ""){
  color.background = "#000000"
  color.grid.major = "#000000"
  color.text = "#ffffff"
  color.axis = "#ffffff"

  theme_bw(base_size=base_size) +
    theme(

      panel.background=element_rect(fill=color.background, color=NA),
      plot.background=element_rect(fill=color.background, color=color.background),
      panel.border=element_rect(color=color.axis),

      panel.grid.major=element_line(color="grey95",size=.1, linetype=3),
      panel.grid.minor=element_blank(),
#       axis.line.x=element_line(color=color.grid.major, size=1),
#       axis.line.y=element_line(color=color.grid.major, size=1),
#       axis.ticks=element_line(color=NA),

      legend.background = element_rect(fill=color.background),
      legend.key = element_rect(fill=color.background, color=NA),
      legend.text = element_text(size=rel(.8),color=color.text),#color.axis.title),
      legend.title = element_text(color=color.text),

      plot.title=element_text(color=color.text, size=rel(1.2)),
      plot.subtitle=element_text(color=color.text, size=rel(0.8)),
      axis.text.x=element_text(size=rel(.95),color=color.text),
      axis.text.y=element_text(size=rel(.95),color=color.text),
      axis.title.x=element_text(size=rel(1),color=color.text, vjust=0),
      axis.title.y=element_text(size=rel(1),color=color.text, vjust=1.25),
      strip.background = element_rect(fill = color.background, colour = color.axis),
      strip.text = element_text(colour = color.text)
    )

}

##----------------------------------------##
##----------------------------------------##

geneBoxPlot_fa2 <- function (tpmTable, gene) {

        if (grepl("FBgn", gene)) {
            description <- subset(snapshots, FBgn_ID == gene)$GeneName
            geneName <- subset(tpmTable, ref_gene_id == gene)$gene_name
            detail <- subset(snapshots, FBgn_ID == gene)$gene_snapshot_text
            
            p <- ggplot(subset(tpmTable, ref_gene_id == gene), aes(Library_Name, TPM, fill = Sex, colour = Sex)) + 
                geom_boxplot(lwd = 0.3) + 
#                 geom_jitter() +
                facet_grid(.~dev_stage, scale = "free_x", space = "free_x") + 
                labs(title = paste(gene, " (", geneName, "): ", description, sep = ""), subtitle = paste(str_wrap(detail, width = 110), sep = "")) + 
                theme_bw() +
                theme(axis.text.x=element_text(size=rel(.95),angle = 45, hjust = 1)) +
                scale_fill_brewer(palette="Dark2") +
                scale_colour_brewer(palette="Set1")
            } else {
            description <- subset(snapshots, GeneSymbol == gene)$GeneName
            geneName <- subset(tpmTable, gene_name == gene)$ref_gene_id
            detail <- subset(snapshots, GeneSymbol == gene)$gene_snapshot_text
            
            p <- ggplot(subset(tpmTable, gene_name == gene), aes(Library_Name, TPM, fill = Sex, colour = Sex)) + 
                geom_boxplot(lwd = 0.3) + 
#                 geom_jitter() +
                facet_grid(.~dev_stage, scale = "free_x", space = "free_x") + 
                labs(title = paste(gene, " (", geneName, "): ", description, sep = ""), subtitle = paste(str_wrap(detail, width = 110), sep = "")) + 
                theme_bw() +
                theme(axis.text.x=element_text(size=rel(.95),angle = 45, hjust = 1)) +
                scale_fill_brewer(palette="Dark2") +
                scale_colour_brewer(palette="Set1")
        }
    return(p)
}

##----------------------------------------##
##----------------------------------------##

heatmap_fa2 <- function(matrix, gene_list, fly_atlas = F, title) {

if (fly_atlas) {
    ## make an annotation bar for "sex"
    col_annot = unique(dplyr::select(sample.info, Sample_Name, Sex, dev_stage))
    rownames(col_annot) = col_annot$Sample_Name
    rownames(col_annot) = gsub("_", " ", rownames(col_annot))
    col_annot = subset(col_annot, select = c("Sex", "dev_stage"))

    # # set colors
    mat_colors = list(Sex = brewer.pal(3, "Set1"), dev_stage = brewer.pal(2, "Set2"))
    names(mat_colors$Sex) = unique(sample.info$Sex)
    names(mat_colors$dev_stage) = unique(sample.info$dev_stage)
} else {
    col_annot = unique(dplyr::select(sampleInfo, Replicate, Status, Female, Male, Handler))
    rownames(col_annot) = col_annot$Replicate
    rownames(col_annot) = gsub("_", " ", rownames(col_annot))
    col_annot = subset(col_annot, select = c("Status", "Female", "Male", "Handler"))

    # # set colors
    # mat_colors <- list(group = brewer.pal(3, "Set1"))
    # names(mat_colors$group) <- unique(col_groups)
}
    

## process tmm matrix
data = subset(matrix, rownames(matrix) %in% gene_list)
data = log2(data+1)
data = as.data.frame(t(scale(t(data), scale=F)))
data[data < -2] = -2
data[data > 2] = 2
init_cols = colnames(data)
data$FBgn_ID = rownames(data)
data = merge(data, FBgn_to_symbol, by.x = "FBgn_ID", by.y = "primary_FBgn", all.x = T)
rownames(data) = data$gene_symbol
data = subset(data, select = init_cols)
colnames(data) = gsub("_", " ", colnames(data))

# plotter
p <- pheatmap(
  mat               = data,
  main              = title,
  color             = inferno(100),
  border_color      = NA,
  show_colnames     = TRUE,
  show_rownames     = TRUE,
  annotation_col    = col_annot,
  drop_levels       = TRUE,
#   cluster_col    = FALSE,
  annotation_names_row = F,
  fontsize          = 8    
)
    return(p)
    }

plot.qq <- function(vec, title.str="qqplot", hit.idx=NULL) {
  col.vec <- rep(1, length(vec))
  if (length(hit.idx) > 0) {
    col.vec[1:length(hit.idx)] <- 2
    col.vec <- rev(col.vec)
  } 
  qqplot(-log10(ppoints(length(vec))), -log10(vec), ylab="-log10(observed)", xlab="-log10(expected)", main=title.str, col=col.vec)
  abline(a=0, b=1)
}

##----------------------------------------##
##----------------------------------------##

geneBarPlot <- function (object, gene, show_reps = F) 
{
    swisprotName <- subset(annot.sum, FBgn_ID == gene)$SwissProt_BlastX_Description
    melOrth <- subset(melOrthsAll, FBgn_ID == gene)$mel_GeneSymbol
    geneName <- subset(annot.sum, FBgn_ID == gene)$gene_name
    
    tmpDF = subset(object, FBgn_ID == gene)
    tmpDF.se = summarySE(tmpDF, measurevar = "TPM", groupvars = c("FBgn_ID", "gene_name", "Sample", "Sex", "Tissue", "Status", "Time"))
    
    p <- ggplot(subset(object, FBgn_ID == gene), 
        aes(Time, TPM, fill = Status)) + 
        facet_grid(.~Tissue, scale = "free_x", space = "free_x") +
        geom_bar(data = tmpDF.se, mapping = aes(Time, fill = Status), stat = "identity", position = position_dodge(.6), size = 0.7, width = 0.5) +
                geom_errorbar(data = tmpDF.se, aes(Time, colour = Status, ymin=TPM-se, ymax=TPM+se), width=.2, position=position_dodge(.6), size = 0.5) + 
        labs(title = paste(melOrth, sep = ""), 
            subtitle = paste(swisprotName, sep = "")) + theme_bw() + 
        theme(axis.text.x = element_text(angle = 45, 
            hjust = 1, size = 12), axis.text.y = element_text(size = 12), 
            axis.title.x = element_blank()) + 
            scale_fill_manual(values = c("#3f5a2a", 
        "#ffb200", "#00b5d4")) + scale_colour_manual(values = c("black","black","black"))
    if (show_reps){
            p <- p + geom_point(data = tmpDF, mapping = aes(Time, TPM, colour = Status), position = position_jitterdodge(jitter.width = 0,dodge.width = 0.6), size = 0.8)
            }
    return(p)
}

##----------------------------------------##
##----------------------------------------##

RT.geneBarPlot <- function (object, gene, show_reps = F) 
{
    swisprotName <- subset(annot.sum, FBgn_ID == gene)$SwissProt_BlastX_Description
    melOrth <- subset(melOrthsAll, FBgn_ID == gene)$mel_GeneSymbol
    geneName <- subset(annot.sum, FBgn_ID == gene)$gene_name
    
    tmpDF = subset(object, FBgn_ID == gene & Tissue == "repTract")
    tmpDF.se = summarySE(tmpDF, measurevar = "TPM", groupvars = c("FBgn_ID", "gene_name", "sample", "Sex", "Tissue", "Status", "Time"))
    
    p <- ggplot(subset(object, FBgn_ID == gene & Tissue == "repTract"), 
        aes(Time, TPM, fill = Status)) + 
        geom_bar(data = tmpDF.se, mapping = aes(Time, fill = Status), stat = "identity", position = position_dodge(.6), size = 0.7, width = 0.5) +
        geom_errorbar(data = tmpDF.se, aes(Time, ymin=TPM-se, ymax=TPM+se), width=.2, position=position_dodge(.6), size = 1) + 
        labs(title = paste(melOrth, sep = ""), 
            subtitle = paste(swisprotName, sep = "")) + 
        theme(legend.position = "none", axis.text.x = element_text(angle = 45, 
            hjust = 1, size = 12), axis.text.y = element_text(size = 12), 
            axis.title.x = element_blank()) + scale_fill_manual(values = c("#3f5a2a", 
        "#ffb200", "#00b5d4"))
    if (show_reps){
            p <- p + geom_point(data = tmpDF, mapping = aes(Time, TPM, colour = Status), position = position_jitterdodge(jitter.width = 0,dodge.width = 0.6), size = 0.8)
            }
    return(p)
}

##----------------------------------------##
##----------------------------------------##


RT.genePointPlot <- function (object, gene, show_reps = F, nonRuv = F) 
{
    swisprotName <- subset(annot.sum, FBgn_ID == gene)$SwissProt_BlastX_Description
    melOrth <- subset(melOrthsAll, FBgn_ID == gene)$mel_GeneSymbol
    geneName <- subset(annot.sum, FBgn_ID == gene)$gene_name
    tmpDF = subset(object, FBgn_ID == gene)
    if (nonRuv) {
        tmpDF = subset(tmpDF, Tissue == "repTract")
    }
    tmpDF.se = summarySE(tmpDF, measurevar = "TPM", groupvars = c("FBgn_ID", 
        "sample", "Status", "Time", "gene_name"))
    tmpDF.se$Time <- factor(tmpDF.se$Time, levels = c("virgin", 
        "3hpm", "6hpm", "12hpm"))
    tmpDF.se$Status <- factor(tmpDF.se$Status, levels = c("virgin", 
        "conspecific", "heterospecific"))
    p <- ggplot() + 
        geom_point(data = tmpDF.se, mapping = aes(Time, TPM, colour = Status), stat = "identity", position = position_dodge(0.6), size = 2) + 
        geom_errorbar(data = tmpDF.se, aes(Time, ymin = TPM - se, ymax = TPM + se, colour = Status), width = 0.2, position = position_dodge(0.6), size = 0.5) + 
        labs(title = paste(gene, " (", geneName, ")", sep = ""), subtitle = paste(melOrth, ": ", swisprotName, sep = "")) + 
        theme_bw() + 
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 13), axis.text.y = element_text(size = 12), 
        axis.title.x = element_blank()) + scale_colour_manual(values = c("#3f5a2a", "#ffb200", "#00b5d4"))
    if (show_reps) {
        p <- p + geom_point(data = tmpDF, mapping = aes(y = TPM, x = Time, colour = Status), position = position_jitter(w = 0.1, h = 0), size = 0.8)
    }
    return(p)
}