library(pacman)
#p_load(cmapR)
p_load(magrittr)
p_load(limma)
p_load(pheatmap)

#setwd('h:/LabMembers/Karsten/Projects/CPTAC3/PTRC/Westbrook/20180817_limma/')

#source('c:/Users/karsten/Dropbox/Devel/R-code/my_io.r')


##gct.in <- 'proteome-ratio-norm-NArm-UP2_PTPN12status.gct'
#gct.in <- c( 
# 
#  Proteome='//flynn-cifs/prot_proteomics/LabMembers/Karsten/Projects/CPTAC3/PTRC/Westbrook/data/v2/Karsten_MultiOmics_analysis/proteome-ratio-norm-NArm-UP2_PTPN12status.gct'#,
  #RNA='//flynn-cifs/prot_proteomics/LabMembers/Karsten/Projects/CPTAC3cmap/PTRC/Westbrook/data/v2/Karsten_MultiOmics_analysis/mRNA_TPM_gct_row-median_2comp.gct',
  #pSTY='//flynn-cifs/prot_proteomics/LabMembers/Karsten/Projects/CPTAC3/PTRC/Westbrook/data/v2/Karsten_MultiOmics_analysis/phosphoproteome-ratio-norm-NArm_PTPN12status.gct'
#)



## ###################################################################
## linear model for data with repeated measurements
lm.repeats <- function(data, design, contrasts, repeats='pdx'){
  
  ## correlation of repeated measurements
  dupfit <- duplicateCorrelation(data, design, block=design[ , repeats])
  
  ## linear model
  fit1 <- lmFit(data, design, block=design[, repeats], correlation = dupfit$consensus.correlation)
  fit2 <- eBayes(fit1, robust=T)
  
  ## contrasts
  fit3 <- contrasts.fit(fit1, contrasts)
  FC <- fit3$coefficients
  
  ## add p-values etc
  tt <- topTable(fit2, number = nrow(data), sort.by = 'none', adjust.method = 'fdr' )
  
  ## assemble results
  results <- data.frame(FC, tt, data)
  
  return(invisible(results))
}


## ##########################################################
##
## can be applied to a single gct file
##
run.lm.repeats <- function(gct,                     ## path to gct file
                   group = 'PTPN12Class_IHCbased',      ## group to compare; track in gct cdesc
                   repeats='Participant
',              ## repeated measurements; track in gct cdesc
                   group.levels=c('high', 'low'), ## group levels of interest
                   plot=F,                     ## heatmap of significant features
                   plot.fdr=0.05,              ## FDR cutoff for heatmap
                   plot.cdesc=rev(c('PTPN12Class', 'PTPN12IHCscores', 'CSresponseRank')), ## heatmap annotation tracks
                   plot.main='',                ## heatmap title
                   show_rownames = F, 
                   scale = 'row', 
                   cluster_rows = F,
                   pdf=F,
                   #gct=F,
                   ...                          ## passed to pheatmap
                   ){
  require(pacman)
  #p_load(cmapR)
  p_load(glue)
  p_load(limma)
  p_load(magrittr)
  
  
  ## parse gct
  #gct <- parse.gctx2(gct.in)
  mat <- gct@mat
  cdesc <- gct@cdesc
  rdesc <- gct@rdesc
  rid <- gct@rid
  
  ## exlcude samples without annotation
  idx <- which(cdesc[, group] %in% group.levels)
  mat <- mat[, idx]
  cdesc <- cdesc[idx, ]

  ## ###############################################
  ## design matrix
  group <- cdesc[, group] %>% as.factor
  rep <- cdesc[, repeats] %>% as.factor %>% as.numeric

  #design <- model.matrix(~ 0 + group * rep)
  design <- model.matrix(~ 0 + group + rep)
  
  colnames(design) <- make.names(colnames(design))
  

  ## ###############################################
  ## contrast matrix
  contrasts <- makeContrasts( contrasts=glue('group{group.levels[1]}-group{group.levels[2]}') , levels=design)
  
  ## ################################################
  ## run the model
  lm <- lm.repeats(mat, design, contrasts, repeats='rep')

  ## ################################################
  ## append to rdesc
  rdesc <- data.frame(lm[, 1:grep('^adj.P.Val$', colnames(lm))], rdesc)
  
  ## GCT
  out <- new('GCT') 
  out@mat <- mat
  out@cdesc <- cdesc
  out@rdesc <- rdesc
  out@cid <- colnames(mat)
  out@rid <- rownames(mat)
    
  
  if(plot){
    p_load(pheatmap)
    sig.idx <- which(lm$adj.P.Val < plot.fdr)
    filename = NA
    if(pdf)
      filename = glue('pheatmap_limma_{plot.main}_nsig_{length(sig.idx)}.pdf')
    
    p <- try(pheatmap( mat[ sig.idx, ], annotation_col = cdesc[,  plot.cdesc], 
              show_rownames = show_rownames, scale = scale, cluster_rows = cluster_rows,
              main=glue('{plot.main} N={length(sig.idx)} @ FDR {plot.fdr*100}%'),
              filename = filename, 
              ...))
    if(class(p) == 'try-error'){
      p <- try(pheatmap( mat[ sig.idx, ], annotation_col = cdesc[,  plot.cdesc], 
                         show_rownames = show_rownames, scale = scale, cluster_rows = F,
                         main=glue('{plot.main} N={length(sig.idx)} @ FDR {plot.fdr*100}%'),
                         filename = filename, 
                         ...))
    }
  }
  ## export GCT
  fn <-  sub('\\.gct','-limma', gct.in) %>% sub('.*/' ,'', .)

  return(out) 
}


