#################################################
#                                               #
#        KEGG MAPPER UTILITY FUNCTIONS          #
#                                               #
#################################################

### MAPPING DATA LOCATIONS (PART OF REPO, SHOULD UPDATE REGULARLY OR CREATE PIPELINE TO GENERATE AUTOMATICALLY)
KEGGGENE_TO_KEGGPATHWAY_MAP='etc/kegg_data/kegg_human_genes_to_pathways.csv'
GENESYMBOL_TO_KEGGGENE_MAP='etc/ncbi_data/gene_to_entrez_and_kegg.csv'

### OTHER PROTIGY NAMES THAT SHOULD BE GLOBALIZED
COMPARISON_LIST_FIELD_NAME = 'grp.comp.all'
FOLD_CHANGE_NAME_BASE = 'logFC'
SIGNIFICANCE_NAME_BASE = 'adj.P.Val'
DATA_OBJECT_NAME = 'data'
DATA_OUTPUT_OBJECT_NAME = 'output'
GENE_SYMBOL_NAME = 'id.mapped'  # POTENTIALLY TO BE REPLACED BY ENTREZ GENE ID LATER
ID_CONCAT_NAME = 'id.concat'

### OTHER MISC GLOBALS
JOINER_CHARACTER = '.'
KEGG_GENE_ACCESSION_FIELDNAME = 'keggGene'


map_kegg_pathways <- function(protigy_parameters, protigy_results) {
  #######################################################################################
  #                                                                                     #
  #          MAIN KEGG MAPPING WORKFLOW                                                 #
  #                                                                                     #
  #                                                                                     #
  # param: protigy_parameters - global.param.imp from protigy data structure            #
  #                           - must contain 'grp.comp.all' list                        #
  #                                                                                     #
  # param: protigy_results    - global.results.imp from protigy data structure          #
  #                           - must contain 'data' list with 'output' dataframe        #
  #                                                                                     #
  #                                                                                     #
  #######################################################################################

  # Gather potential comparisons
  
  # Map gene symbols to KEGG ids
  
  # Main loop for each comparison in Protigy results
  
    # Mark genes up/down/both based on statistics (set the filtering level!)
  
    # Associate pathways to genes
  
    # Compute pathway enrichments
  
  # Organize conditions and data
  
  # Generate HTML output document
  
}


gather_comparisons_and_make_column_names <- function(protigy_parameters) {
  #######################################################################################
  #                                                                                     #
  # function: gather_comparisons_and_make_column_names                                  #
  #   makes a dataframe containing all relevant comparisons and names the columns for   #
  #   use in determining significance and directionality                                #
  #                                                                                     #
  # param: protigy_parameters - global.param.imp from protigy data structure            #
  #                           - must contain 'grp.comp.all' list                        #
  #                                                                                     #
  # returns: dataframe with the following fields:                                       #
  #    comparison_name                                                                  #
  #    significance_name                                                                #
  #    foldchange_name                                                                  #
  #                                                                                     #
  #######################################################################################
  
  comparison_names = protigy_parameters[[COMPARISON_LIST_FIELD_NAME]]
  
  comparison_nomenclature_df = data.frame(
    comparison_name = comparison_names,
    significance_name = paste(SIGNIFICANCE_NAME_BASE,comparison_names, sep=JOINER_CHARACTER),
    foldchange_name = paste(SIGNIFICANCE_NAME_BASE,comparison_names, sep=JOINER_CHARACTER),
    row.names = comparison_names
  )

  return(comparison_nomenclature_df)    
}


annotate_genes_with_kegg_ids <- function(protigy_results){

  
  # Load the gene symbol to kegg accession mapping
  gene_to_kegg_table = read.csv(GENESYMBOL_TO_KEGGGENE_MAP, stringsAsFactors = FALSE, row.names = 1)
  
  # Extract the whole protigy result table
  protigy_result_table = protigy_results[[DATA_OBJECT_NAME]][[DATA_OUTPUT_OBJECT_NAME]]

  # Extract just the gene symbols and row ids from the protigy results, note uppercasing of gene symbol names
  temporary_data_frame_for_merging = data.frame(
    id = rownames(protigy_result_table),
    id.concat = protigy_result_table[[ID_CONCAT_NAME]],
    id.mapped = toupper(protigy_result_table[[GENE_SYMBOL_NAME]]),
    row.names = rownames(protigy_result_table),
    stringsAsFactors = FALSE 
  )
  
  # Perform the merge
  temp_merge_result = merge(
    x = temporary_data_frame_for_merging,
    y = gene_to_kegg_table,
    by.x = 'id.mapped',
    by.y = 'Approved.symbol',
    all.x = TRUE,
    sort = FALSE
  )

  # Reorder to match original order and select the kegg accessions
  rownames(temp_merge_result) <- temp_merge_result[['id']]    
  ordered_kegg_accessions <- temp_merge_result[rownames(protigy_result_table),KEGG_GENE_ACCESSION_FIELDNAME]
  
  # Rebind to original result table
  protigy_result_table[[KEGG_GENE_ACCESSION_FIELDNAME]] <- ordered_kegg_accessions
  
  return(protigy_result_table)
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
