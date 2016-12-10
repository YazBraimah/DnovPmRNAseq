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
  p <- ggplot(subset(object, FBgn_ID == gene_id | gene_name == gene_id), aes(time, TPM, fill = condition)) + geom_bar(position=position_dodge(), stat="identity") + geom_errorbar(aes(ymin=TPM-se, ymax=TPM+se), width=.2, position=position_dodge(.9)) + facet_grid(~tissue, scales="free_x", space = "free_x") + labs(title = paste(gene_id," (", geneName,"), ",coords,"\n","Ka/Ks = ", omega,"        mel. orth.: ",melOrth,"\n",swisprotName, sep = "")) + scale_fill_manual(values = c("#84a955", "#965da7", "#bc5d41"))
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
  p <- ggplot(subset(object, FBgn_ID == gene_id & tissue == "RT"| gene_name == gene_id & tissue == "RT"), aes(time, TPM, fill = condition)) + geom_bar(position=position_dodge(), stat="identity") + geom_errorbar(aes(ymin=TPM-se, ymax=TPM+se), width=.2, position=position_dodge(.9)) + labs(title = paste(gene_id," (", geneName,"), ",coords,"\n","Ka/Ks = ", omega,"        mel. orth.: ",melOrth,"\n",swisprotName, sep = "")) + scale_fill_manual(values = c("#84a955", "#965da7", "#bc5d41"))
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
