#################################################
#                                               #
#        KEGG MAPPER UTILITY FUNCTIONS          #
#                                               #
#################################################

### MAPPING DATA LOCATIONS (PART OF REPO, SHOULD UPDATE REGULARLY OR CREATE PIPELINE TO GENERATE AUTOMATICALLY)
KEGGGENE_TO_KEGGPATHWAY_MAP='../etc/kegg_data/kegg_human_genes_to_pathways.csv'
GENESYMBOL_TO_KEGGGENE_MAP='../ncbi_data/gene_to_entrez_and_kegg.csv'

map_kegg_pathways <- function() {
  #################################################
  #                                               #
  #          MAIN KEGG MAPPING WORKFLOW           #
  #                                               #
  #################################################
  
  # Gather potential comparisons
  
  # Map gene symbols to KEGG ids
  
  # Main loop for each comparison in Protigy results
  
    # Mark genes up/down/both based on statistics (set the filtering level!)
  
    # Associate pathways to genes
  
    # Compute pathway enrichments
  
  # Organize conditions and data
  
  # Generate HTML output document
  
}


prepareKEGGhelperFunction <- function(report, fcfield = 'logFC', sigfield = 'adj.P.Val', sig_threshold = 0.1) {
  
  report$direction_site <-apply(
    report,
    1,
    function (x) { 
      if (as.numeric(x[sigfield]) > sig_threshold 
          | 
          is.na(as.numeric(x[sigfield])) 
          | 
          is.na(as.numeric(x[fcfield]))
          ) 
      {
        return("unch") 
      } else { 
        return(if (as.numeric(x[fcfield]) < 0) "down" else "up") 
      }
    }
  )
  
  combinations<-(unique(report[,c('KEGG','direction_site')]))
  collapsedDirections<-"unch"
  if (all(c('up','down') %in% colnames(table(combinations)))) {
    collapsedDirections<-apply(table(combinations),1,function (x) { if (x['up']==1 & x['down']==1) { return("both") } 
      else { if (x['up']==1) {return ("up")} else if (x['down']==1) {return ("down")} else {return("unch")} }})
  } else {
    if (any(c('up','down') %in% colnames(table(combinations)))) { warning("Only one type of directionality (up or down) is present. Color plots are wrong.")}
  }
  combinations.annotated<-cbind(table(combinations),direction=collapsedDirections)
  keggprep.final<-merge(report,combinations.annotated,by.x='KEGG',by.y=0,all.x=TRUE,sort=FALSE)
  keggprep.final$UNIPROTKB<-keggprep.final$best_accession_number
  
  return(keggprep.final)
  
}
