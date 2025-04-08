################################################################################################################
## Filename: global.R
## Created: October 09, 2015
## Author(s): Karsten Krug, Ozan Aygun
##
## Purpose: Shiny-app to perform differential expression analysis, primarily on proteomics data, to perform
##          simple data QC, to interactively browse through the results and to download high-quality result
##          figures.
##
## This file defines global parameters, loads all required R-packages and defines the actual functions to perform data filtering,
## data normalization, the moderated test statistics, and visualization. The code for moderated t-tests,
## two-component normalization and the reproducibility filter has been written by D. R. Mani and
## adopted by me for integration into a Shiny-Server environment.
##
##
## Last updated April 20, 2022 by Natalie Clark (nclark@broadinstitute.org) - v 1.0
##################################################################################################################
#options( stringsAsFactors = F )

## R package managing tool
if (!require("pacman")) install.packages ("pacman")
require('pacman')

if (!require("BiocManager")) install.packages ("BiocManager")
require(BiocManager)
options(repos = BiocManager::repositories())

#################################################################
## global parameters
#################################################################
## use PACMAN R package manager?
## set to FALSE if deployed to RStudio Connect 
PACMAN <- TRUE
## version number
VER <- "1.1.8"
## maximal file size for upload
MAXSIZEMB <<- 1024
## list of strings indicating missing data
NASTRINGS <<- c("NA", "<NA>", "#N/A", "#NUM!", "#DIV/0!", "#NA", "#NAME?", "na", "#VALUE!","Na","nA","NaN")
## separator tested in the uploaded file
SEPARATOR <<- c('\t', ',', ';')
## Colors used throughout the app to color the defined groups
GRPCOLORS <<- c(RColorBrewer::brewer.pal(8, "Set1"), RColorBrewer::brewer.pal(8, "Dark2"), RColorBrewer::brewer.pal(8, "Set2"), terrain.colors(20), cm.colors(20), topo.colors(20))
## number of characters to display in plots/tables for column names
STRLENGTH <<- 25
## operating system
OS <<- Sys.info()['sysname']
## temp directory to write the Excel file
TMPDIR <<- ifelse(OS=='Windows', "./", "/tmp/")
## app name
APPNAME <<- 'Protigy'
## app folder
APPDIR <<- getwd()
## directory to store data files
DATADIR <<- ifelse(OS=='Linux', "/shiny-data/", tempdir())
## User database for shared sessions 
USERDB <<- file.path(APPDIR, 'conf/user-roles.txt')
## R version
RVERSION <- as.numeric(paste(R.version$major, sub('\\..*','',R.version$minor), sep='.')) 
## email for trouble shooting
MAIL <<- 'nclark@broadinstitute.org'
## URL to configuration app (SSP/RSC only)
CONFAPP <<- 'https://rstudio-connect.broadapps.org/protigy-sessions/'
## PIWIK location
PIWIKURL <<- ''

#################################################################
## load required packages
#################################################################
if(PACMAN){
  p_load(RColorBrewer)
  p_load(shiny)
  p_load(shinydashboard)
  p_load(shinyjs)
  p_load(shinyjqui)
  p_load(shinyalert)
  p_load(markdown)
  p_load(tippy)

  p_load(magrittr)
  p_load(tibble)
  p_load(dplyr)
  p_load(stringr)
  
  ## heatmap
  p_load(pheatmap)
  p_load(heatmaply)

  ## clustering  
  p_load(ape)
  p_load(dendextend)
  p_load(scales)
  p_load(gtable)
  p_load(fastcluster)
  
  ## rmarkdown
  p_load(rmarkdown)
  p_load(knitr)
  
  ## moderated tests
  p_load(limma)
  p_load(statmod)
  
  ## multiscatter
  p_load(hexbin)
  p_load(Hmisc)
  p_load(grid)

  p_load(scatterplot3d)
  p_load(plotly)
  
  ## Excel
  p_load(readxl)
  p_load(WriteXLS)
  
  ## reproducibility filter
  p_load(reshape)
  p_load(nlme)
  p_load(BlandAltmanLeh)
  
  ## normalization
  p_load(preprocessCore)
  p_load(vsn)
  
  ## normalization 2-component
  p_load (mice)
  p_load (mixtools)
  p_load (mclust)
  
  ## table preview
  p_load(DT)
  
  ## label placements without overlap
  p_load(maptools)
  p_load(ggrepel)
  
  p_load(UpSetR)
  p_load(gridExtra)
  p_load(vioplot)

  ## id mapping
  p_load(RSQLite)
  p_load(org.Hs.eg.db)
  p_load(org.Mm.eg.db)
  p_load(org.Rn.eg.db)
  p_load(org.Dr.eg.db)
   
  ## enrichment analysis
  p_load(gprofiler2)
  
  #other needed packages
  p_load(seriation)
  p_load(preprocessCore)
  p_load(car)

} else { ## RStudio Connect
  
  library(RColorBrewer)
  library(shiny)
  library(shinydashboard)
  library(shinyjs)
  library(shinyjqui)
  library(shinyalert)
  library(tippy)

  library(magrittr)
  library(tibble)
  library(dplyr)
  library(stringr)

  ## heatmap
  library(pheatmap)
  library(heatmaply)
  
  ## clustering
  library(ape)
  library(dendextend)
  library(scales)
  library(gtable)
  library(fastcluster)
  
  ## rmarkdown
  library(rmarkdown)
  library(markdown)
  library(knitr)
  
  ## moderated tests
  library(limma)
  library(statmod)
  
  ## multiscatter
  library(hexbin)
  library(Hmisc)
  library(grid)
  
  library(scatterplot3d)
  library(plotly)
  
  ## Excel
  library(WriteXLS)
  library(readxl)
  
  ## reproducibility filter
  library(reshape)
  library(nlme)
  library(BlandAltmanLeh)
  
  ## normalization
  library(preprocessCore)
  library(vsn)
  
  ## normalization 2-component
  library (mice)
  library (mixtools)
  library (mclust)
  
  ## table preview
  library(DT)
  
  ## label placements without overlap
  library(maptools)
  library(ggrepel)
  
  library(UpSetR)
  library(gridExtra)
  library(vioplot)

  ## id mapping
  library(RSQLite)
  library(org.Hs.eg.db)
  library(org.Mm.eg.db)
  library(org.Rn.eg.db)
  library(org.Dr.eg.db)
   
  ## enrichment analysis
  library(gprofiler2)
}

## required to install limma in R >= 3.5
if(RVERSION >= 3.5){
  if(PACMAN)
    p_load(BiocManager)
  else
    library(BiocManager)
  options(repos = BiocManager::repositories())
}

## pca
# ChemometricsWithR - removed from CRAN
# use GH version instead
if(!require("ChemometricsWithR")){
  install.packages("remotes")
  library(remotes)
  install_github("rwehrens/ChemometricsWithR")
  library(ChemometricsWithR)
}

## morpheus
#p_load_gh('morpheus')
#if(!require(morpheus))
#  devtools::install_github('cmap/morpheus.R')
#p_load(morpheus)

source('src/modT.R')
source('src/helptext.R')
source('src/my_io.R')
source('src/gct-io.R')
source('src/plots.R')
source('src/manage_sessions.R')
#source('src/linear_model.R')


##################################################
##    import header descriptions
headerDesc <- read_excel('docs/description-column-headers.xlsx')


## #####################################
## CSS for loading animation
appCSS <- "
#loading-content {
  position: absolute;
  background: #000000;
  opacity: 0.9;
  z-index: 100;
  left: 0;
  right: 0;
  height: 100%;
  text-align: center;
  color: #FFFFFF;
}
"
## #################################################################
##
##              color scheme for PPI databases
##
## #################################################################
ppi.db.col <- c('InWeb'='deepskyblue',
                'BioGRID'='yellow',
                'Reactome'='magenta',
                'Shared'='orange'
                )
## ##################################################################
##
##                    import PPi databases
##
## ##################################################################
import.ppi.db <- function(){
  
  if( file.exists('ppi/ppi.RData')){
    load('ppi/ppi.RData')
  } else {
    ppi <- list()
    
    ## ###################################
    ##             InWeb
    ppi$iw <- readRDS('ppi/core.psimitab.rds')
    ## uniprot
    ppi$iw$V1 <- toupper(sub('uniprotkb\\:', '', ppi$iw$V1))
    ppi$iw$V2 <- toupper(sub('uniprotkb\\:', '', ppi$iw$V2))
    ## gene names
    ppi$iw$V5 <- toupper(sub('^uniprotkb\\:(.*?)\\(gene name\\).*', '\\1', ppi$iw$V5))
    ppi$iw$V6 <- toupper(sub('^uniprotkb\\:(.*?)\\(gene name\\).*', '\\1', ppi$iw$V6))
    
    ## ###################################
    ##            biogrid
    ppi$bg <- readRDS('ppi/BIOGRID-HUMAN-3.4.147.mitab.rds')
    ppi$bg$Alt.IDs.Interactor.A <- sapply(strsplit(ppi$bg$Alt.IDs.Interactor.A, '\\|'), function(x) sub('.*\\:','',x[2]) )
    ppi$bg$Alt.IDs.Interactor.B <- sapply(strsplit(ppi$bg$Alt.IDs.Interactor.B, '\\|'), function(x) sub('.*\\:','',x[2]) )
    
    ## ###################################
    ##           Reactome
    ppi$react <- readRDS('ppi/homo_sapiens.mitab.interactions.rds')
    ppi$react$alternative.id.A <- sub('^.*_(.*?)\\(.*$','\\1', ppi$react$alternative.id.A)
    ppi$react$alternative.id.B <- sub('^.*_(.*?)\\(.*$','\\1', ppi$react$alternative.id.B)
    
    # export 
    save(ppi, file='ppi/ppi.RData')
  }
  return(ppi)
}

# run the function to import 
ppi <- import.ppi.db()


## ###############################################
##
##        map uniprot/refseq to gene names
## n.try = number of ids taken from 'ids' to try to
##         determine organism
## ###############################################
mapIDs <- function(ids,rdesc=NULL,
                   n.try=100
                   ){
    withProgress(message='Mapping gene names...', {

            ## ###################################
            ##           id type
            ## ###################################
            keytype <- 'UNKNOWN'
            ## Uniprot or RefSeq or Ensembl?
            if(length(grep('^(Q|P|O|A|E|H|F)', ids)) > 0)
                keytype='UNIPROT'
            if(length(grep('^(NP_|XP_|YP_)', ids)) > 0)
                keytype='REFSEQ'
            if(length(grep('ENSP', ids)) > 0)
              keytype='ENSEMBLPROT'

            ## ###################################
            ##        extract query strings
            ## ###################################
            if(keytype == 'UNIPROT' ){
              id.query <- sub('(-|;|\\.|_|\\|).*', '', ids) ## first id
            } else if(keytype == 'REFSEQ') {
              id.query <- sub('(\\.|;).*', '', ids) ## first id
            } else if(keytype=='ENSEMBLPROT') {
              id.query <- sub('(\\.|;).*', '', ids) ## first id
            }else {
              id.query <- ids
            }
            names(id.query) <- ids

            ## ###################################
            ##          determine organism
            ## ###################################
            orgtype <- 'UNKNOWN'
            if(keytype != 'UNKNOWN'){
              
              # try human 
              if(orgtype == 'UNKNOWN'){
                cat('trying human...\n')
                id.map.tmp <- try( mapIds(org.Hs.eg.db, keys=id.query[ sample(1:length(id.query), n.try)] , column=c('SYMBOL'), keytype=keytype, multiVals='first') )
                if(class(id.map.tmp) != 'try-error'){
                  orgtype='HSA'
                }
              }
              # try mouse
              if(orgtype == 'UNKNOWN'){ 
                cat('trying mouse...\n')
                id.map.tmp <- try( mapIds(org.Mm.eg.db, keys=id.query[ sample(1:length(id.query), n.try)] , column=c('SYMBOL'), keytype=keytype, multiVals='first') )
                  if(class(id.map.tmp) != 'try-error'){
                    orgtype='MMU'
                  }
              }
              # try rat
              if(orgtype == 'UNKNOWN'){
                cat('trying rat...\n')
                id.map.tmp <- try( mapIds(org.Rn.eg.db, keys=id.query[ sample(1:length(id.query), n.try)] , column=c('SYMBOL'), keytype=keytype, multiVals='first') )
                if(class(id.map.tmp) != 'try-error'){
                  orgtype='RNO'
                }
              }
              # try zebrafish
              if(orgtype == 'UNKNOWN'){
                cat('trying zebrafish...\n')
                id.map.tmp <- try( mapIds(org.Dr.eg.db, keys=id.query[ sample(1:length(id.query), n.try)] , column=c('SYMBOL'), keytype=keytype, multiVals='first') )
                if(class(id.map.tmp) != 'try-error'){
                  orgtype='DRE'
                }
              }
            }
  
            ## ##################################
            ## map
            
            #if geneSymbol column is included in rdesc, use that
            if(!is.null(rdesc)&("geneSymbol"%in%colnames(rdesc))){
              id.map.tmp <- sub('(-|;|\\.|_|\\|).*', '', rdesc$geneSymbol) #take first if there is a list
            }else{
              if(keytype != 'UNKNOWN' & orgtype != 'UNKNOWN'){
                if(orgtype == 'HSA')
                  id.map.tmp <- try(mapIds(org.Hs.eg.db, keys=id.query , column=c('SYMBOL'), keytype=keytype, multiVals='first'))
                if(orgtype == 'MMU')
                  id.map.tmp <- try(mapIds(org.Mm.eg.db, keys=id.query , column=c('SYMBOL'), keytype=keytype, multiVals='first'))
                if(orgtype == 'RNO')
                  id.map.tmp <- try(mapIds(org.Rn.eg.db, keys=id.query , column=c('SYMBOL'), keytype=keytype, multiVals='first'))
                if(orgtype == 'DRE')
                  id.map.tmp <- try(mapIds(org.Dr.eg.db, keys=id.query , column=c('SYMBOL'), keytype=keytype, multiVals='first'))
              } else {
                id.map.tmp <- c()
              }
            }

            if(class(id.map.tmp) == 'try-error' | is.null( class(id.map.tmp) ) | class(id.map.tmp) == 'NULL' ){

              id.map <- data.frame(id=names(id.query), id.query=id.query, id.mapped=names(id.query), id.concat=ids, stringsAsFactors=F)
              #cat('test2\n')
              keytype <- 'UNKNOWN'
              
            } else {
    
              ## if successful
              id.mapped <- id.map.tmp
              id.mapped[which(is.na(id.mapped) | id.mapped=="")] <- ids[which(is.na(id.mapped) | id.mapped=="")]
              id.map.tmp[which(is.na(id.map.tmp) | id.map.tmp=="")] <- 'NotFound'
              
              id.map <- data.frame(id=names(id.query), id.query=id.query, id.mapped=as.character(id.mapped), id.concat=paste(ids, id.map.tmp, sep='_'), stringsAsFactors=F)
            
            }

            ## results
            res <- list()
            res[[1]] <- keytype
            res[[2]] <- id.map
            res[[3]] <- orgtype
            names(res) <- c('keytype', 'id.map', 'orgtype')
    })

    return(res)
}

############################################################
## given a vector of group names, return all distinct
## pairwise comparison
## mode = 'sort' : group names will be sorted alphanumerically
############################################################
get_pairwise_group_comparisons <- function(groups.unique, 
                                           mode=c('sort') ## sort: group names are ordered alphanumerically
                                                          ##       e.g. "Exp_B" and "Exp_A" -> "Exp_A.vs.Exp_B"
                                                          ##            which translates to Exp_B over Exp_A             
                                           
                                           ){
  
  mode <- match.arg(mode)
  
  groups.comp <- c()
  count=1
  
  for(i in 1:(length(groups.unique)-1))
    for(j in (i+1):length(groups.unique)){
      
      ## order alphabetically
      if(mode == 'sort')
        groups.tmp <- sort( groups.unique[c(i,j)] )
      else
        groups.tmp <- groups.unique[c(i,j)]
      
      groups.comp[count] <- paste(groups.tmp[1], groups.tmp[2], sep='.vs.')
      count <- count+1
  }
  return(groups.comp)
  
}


#################################################################################
##
##   calculate average abundance per group groups 
##
#################################################################################
calculate_fc <- function(tab, grp.vec, groups.comp, test, 
                         mode='sort',  ## for assigning pairwise comparisons
                                      ## sort: group names are ordered alphanumerically
                                      ##       e.g. "Exp_B" and "Exp_A" -> "Exp_A.vs.Exp_B",
                         center=F     ## if TRUE median centering will be applied
                         ){
  #browser()
  groups <- unique(grp.vec)
  
  ## ##########################################
  ## average groups
  group_avg <- lapply(groups, function(gg){
    gg.idx = names(grp.vec)[ which(grp.vec == gg) ]
    apply(tab[, gg.idx], 1, mean, na.rm=T)
  })
  group_avg <- data.frame(Reduce('cbind', group_avg))
  colnames(group_avg) <- groups
    
  ## ##########################################
  ## one sample & mod F: just average the log values
  if(test %in% c("One-sample mod T", "none", "mod F")){
      group_fc <- group_avg
  }

  ## ##########################################
  ## moderated F: average log values and subtract
  ## groups as specified by the user
  #if(test == "mod F"){
  #  groups.comp <- get_pairwise_group_comparisons(groups, mode=mode)
  #}
  
  ## ##########################################
  ## two sample: subtract averaged groups
  if(test %in% c("Two-sample mod T")){
    
    groups.comp.split <- strsplit(groups.comp, split = '\\.vs\\.')
    names(groups.comp.split) <- groups.comp
    
    ## calculate log FCs
    group_fc <- lapply(groups.comp.split, function(x){
      g1 <- x[1] 
      g2 <- x[2]
      
      ## subtract logs
      g2.over.g1 <- group_avg[, g2] - group_avg[, g1]
      g2.over.g1
    })
    group_fc <- data.frame(Reduce('cbind', group_fc))
    colnames(group_fc) <- groups.comp
  }
  
  ## median center data
  if(center)
    group_fc <- apply(group_fc, 2, function(column) column - median(column, na.rm=T))
  
  ## for one-sample/mod F/none, it is average, not FC
  if(test %in% c("One-sample mod T", "none", "mod F")){
    colnames(group_fc) <- paste0('RawAveExpr.', colnames(group_fc))
  }else{
    colnames(group_fc) <- paste0('RawlogFC.', colnames(group_fc))
  }
  return(group_fc)
}

#################################################################################
##
##   - deduplicate ids and show warning in a modal shiny window
##
#################################################################################
de_duplicate_ids <- function(ids, global.param=NULL, show_modal = TRUE){
  
  if(is.null(global.param)) global.param$id.col.value <- 'id'
  
  ## dups
  dup_idx <- which(duplicated(ids))
  dup_n <- length(dup_idx)
  dup_ids <- ids[ dup_idx ]
  
  ## make dups unique
  ids <- make.unique(ids , sep='_')
  
  if(show_modal){  
    
      require(shiny)
    
      ## dups after making unique
      dup_ids_fix <- ids[dup_idx]
      
      
      ## show the first 3 dups in a table
      html_tab <- paste0('<table border="1">
                                       <th style="padding:10px">duplciated</th><th style="padding:10px">deduplicated</th>')
      for(tr in 1:min(c(3, dup_n))){
        html_tab <- paste0(html_tab, 
                           '<tr><td style="padding:10px">', dup_ids[tr], '</td><td style="padding:10px">', dup_ids_fix[tr], '</td></tr>\n')  
      }
      html_tab <- paste0(html_tab, '</table>')
      
      #######################
      ## give a warning about duplicated ids
      showModal(modalDialog(
        size='m',
        title = paste("Warning: ids are not unique!"),
        HTML(paste0('Found <b>', dup_n, '</b> duplicated entries in column <b>', 
                    global.param$id.col.value,
                    '</b> which will be made unique as shown below (first ', min(c(3, dup_n)), ' duplicated ids shown).<br><br>')),
        HTML(html_tab),
        HTML('<br><br>Click "Dismiss" to continue.')
      
        ))
  } ## end if(show.modal)
  
  return(ids)
}


## ##############################################################################
##
##            generate links to external databases
##
## ##############################################################################
link.db <- function(id, # vetcor of ids
                    keytype=c('UNKNOWN', 'UNIPROT', 'REFSEQ','ENSEMBLPROT'),
                    db=c('GENECARDS', 'UNIPROT')){
  
  keytype <- match.arg(keytype)
  db <- match.arg(db)
  
  if(keytype == 'UNIPROT'){
    up.link <- paste("<a href='https://www.uniprot.org/uniprot/", sub('(_|,|;|\\.).*', '', id),"' target='_blank'>", id, "</a>", sep='')
  }
  if(keytype %in% c('REFSEQ', 'ENSEMBLPROT','UNKNOWN')){
    up.link <- paste("<a href='http://www.genecards.org/Search/Keyword?queryString=", sub('^(NP_|NM_|NR_.*?)(_|,|;|\\.).*', '\\1', id),"' target='_blank'>", id, "</a>", sep='')
  }
  return(up.link)
}

#############################################################################################
normalize.data  <- function(data, id.col,
                         method=c('Median',
                                  'Median (non-zero)',
                                  'Quantile', 
                                  'VSN', 
                                  'Median-MAD',
                                  'Median-MAD (non-zero)',
                                  '2-component', 
                                  'Upper-quartile'),
                         grp.vec=NULL  ## if NULL apply global normalization strategy
                                       ## if vector of group assingments, normalization will be applied to each group separately  
                         ){
  
  #################################
  ## global normalization
  if(is.null(grp.vec)){
    
    data.norm <- normalize.data.helper(data, id.col, method=match.arg(method))
  
  ################################  
  ## group-specific normalization  
  } else { 
    
    ## extract groups
    groups <- unique(grp.vec)
    
    ############################################
    ## loop over replicate groups
    for(gg in groups){
      
      gg.idx = names(grp.vec)[ which(grp.vec == gg) ]
    
      data.group <- data[, c(id.col, gg.idx)]
      
      data.group.norm <- normalize.data.helper(data.group, id.col, method=match.arg(method), per_group = TRUE)
     
      if(gg == groups[1]){
        data.norm <- data.group.norm  
      } else {
        data.norm <- cbind(data.norm, data.group.norm[, gg.idx])
      }
    } ## end for
    
  } ## end else
  
  return(data.norm)
}
#############################################################################################
##
##              different normalization methods for expression data
##
## 20160235
#############################################################################################
normalize.data.helper <- function(data, id.col, 
                           method=c('Median',
                                    'Median (non-zero)',
                                    'Quantile', 
                                    'VSN', 
                                    'Median-MAD',
                                    'Median-MAD (non-zero)',
                                    '2-component', 
                                    'Upper-quartile'),
                           per_group=FALSE ## for Median & Median-MAD
                           ){
    cat('\n\n-- normalize data --\n\n')

    method = match.arg(method)

    cat('   normalization method: ', method, '\n')
    
    ids = data[, id.col]
    data = data[ , -grep(paste('^', id.col, '$', sep=''), colnames(data))]
    data <- data.matrix(data)

    ## quantile
    if(method == 'Quantile'){
        p_load("preprocessCore")
        data.norm <- normalize.quantiles(data)
        rownames(data.norm) <- rownames(data)
        colnames(data.norm) <- paste( colnames(data))

        ## shift median to zero
        ## data.norm <- apply(data.norm, 2, function(x) x - median(x, na.rm=T))
    }
    ## median only
    if(method == 'Median'){
      
      data.norm <- apply(data, 2, function(x) x - median(x, na.rm=T))
      colnames(data.norm) <- paste( colnames(data), sep='.')
      
      if(per_group){
        all_medians <- apply(data, 2, median, na.rm=T)
        data.norm <- data.norm + median( all_medians, na.rm=T)
      }
    }
    ## median plus shifting by medians of medians
    if(method == 'Median (non-zero)'){
      
      all_medians <- apply(data, 2, median, na.rm=T)
      data.norm <- apply(data, 2, function(x) x - median(x, na.rm=T))
      
      data.norm <- data.norm + median( all_medians, na.rm=T )
      colnames(data.norm) <- paste( colnames(data), sep='.')
    }
    ## median & MAD
    if(method == 'Median-MAD'){
        data.norm <- apply(data, 2, function(x) (x - median(x, na.rm=T))/mad(x, na.rm=T) )
        colnames(data.norm) <- paste( colnames(data), sep='.')
        
        if(per_group){
          all_medians <-  apply(data, 2, median, na.rm=T)
          data.norm <- data.norm + median( all_medians, na.rm=T)
        }
    }
    ## median & MAD plus shifting by medians of medians
    if(method == 'Median-MAD (non-zero)'){
  
      all_medians <- apply(data, 2, median, na.rm=T)
      data.norm <- apply(data, 2, function(x) (x - median(x, na.rm=T))/mad(x, na.rm=T) )
      
      data.norm <- data.norm + median( all_medians, na.rm=T )
      colnames(data.norm) <- paste( colnames(data), sep='.')

    }
    
    ## 2-component normalization
    if(method == '2-component'){
      
      data.norm.list <- vector('list', ncol(data))
      names(data.norm.list) <- colnames(data)
      
      for(x in colnames(data)){  
          res <- try(two.comp.normalize(data[, x], type="unimodal"))
          data.norm.list[[x]] <- res
          if(class(res) == 'try-error') break;
       }
          
        ## check if all runs were successful
        ## return the 'try-error' object to 
        ## catch the error in server.R
        for(i in 1:length(data.norm.list)){
            if(class(data.norm.list[[i]]) == 'try-error'){
                    msg <- data.norm.list[[i]]
                    return(msg)
            }
        }

      ## if 2-comp was successful on all data column convert list to matrix 
      data.norm = matrix( unlist(lapply(data.norm.list, function(x)x$norm.sample)), ncol=length(data.norm.list), dimnames=list(rownames(data), names(data.norm.list)) )
    }
    ## Upper quartile
    if(method == 'Upper-quartile'){
      data.norm <- apply(data, 2, function(x) x - quantile(x, c(0.75),na.rm=T))
      colnames(data.norm) <- paste( colnames(data), sep='.')
    }
    
    ## VSN - variance stabilizing normalization
    if(method == 'VSN'){
      p_load(vsn)
      data.norm <- justvsn(data)
    }
    
    ## add id column
    data.norm <- data.frame(ids, data.norm)
    colnames(data.norm)[1] <- id.col
    
    cat('\n\n-- normalize data exit--\n\n')
    
    return(data.norm)
}


##############################################################################
##
##  - perform principle component analysis
##  - calculate variances explained by components
##  - plot the results
##
## changelog: 20131001 implementation
##            20131007 3D plot now plot pc1 vs. pc3 vs. pc2
##                     instead of pc1 vs. pc2 vs. pc3
##            20151208 legend
##            20161103 check number of rows (N), 2D plot only
##############################################################################
my.prcomp.static <- function(x, pca.x, pca.y, pca.z, col=NULL, cor=T, plot=T, rgl=F, scale=T, pch=20, cex.points=3, rgl.point.size=30, main="PCA", leg.vec=NULL, leg.col=NULL, ...){

    cex.font = 1.8

    ## number of data columns, N=2 -> 2D plot only
    N <- nrow(x)

    ## color
    if( is.null(col) ) col="black"

    ## perform pca
    pca <- prcomp(x, scale=scale)

    ## calculate variance
    comp.var <- eigen(cov(pca$x))$values

    ## extract the principle components
    pc1=pca$x[,1]
    pc2=pca$x[,2]
    if(N>2)
        pc3=pca$x[,3]

    ##############
    # rgl plot
    ##############
    if(rgl & N > 2){
        p_load(rgl)
        plot3d(pc1, pc2, pc3, xlab=paste("PC 1 (", round(100*comp.var[1]/sum(comp.var),1),"%)", sep=""), ylab=paste("PC 2 (",round(100*comp.var[2]/sum(comp.var),1),"%)", sep=""), zlab=paste("PC 3 (", round(100*comp.var[3]/sum(comp.var),1),"%)", sep=""), type="s", col=col, expand=1.2, size=rgl.point.size)
    }

    ########################################
    # scatterplot 2D/3D
    ########################################
    if( plot){

         p_load(scatterplot3d)

        if(N > 2)
            par(mfrow=c(1,3), mar=c(7,7,3,1))
        if(N <= 2)
            par(mfrow=c(1,2), mar=c(7,7,3,1))

         ## PC 1-2
         plot(pc1, pc2, xlab=paste("PC 1 (", round(100*comp.var[1]/sum(comp.var),1),"%)", sep=""), ylab=paste("PC 2 (",round(100*comp.var[2]/sum(comp.var),1),"%)", sep=""), pch=pch, main=main, col=col, sub=paste("Cumulative variance = ", round(100*sum(comp.var[1:2]/sum(comp.var)),1),"%", sep=""), cex=cex.points, ylim=c( min(pc2),  max(pc2)+.15*max(pc2)), cex.axis=cex.font, cex.lab=cex.font, cex.sub=cex.font )


        if(N > 2) {
            ## PC 1-3
            scatterplot3d( pca$x[,1], pca$x[,3], pca$x[,2], xlab=paste("PC 1 (", round(100*comp.var[1]/sum(comp.var),1),"%)", sep=""), zlab=paste("PC 2 (",round(100*comp.var[2]/sum(comp.var),1),"%)", sep=""), ylab=paste("PC 3 (", round(100*comp.var[3]/sum(comp.var),1),"%)", sep=""), color=col,  cex.symbols=cex.points, pch=pch, main=main, sub=paste("Cumulative variance = ", round(100*sum(comp.var[1:3]/sum(comp.var)),1),"%", sep=""), type="h" )

        }
        ## legend
         plot.new()
         plot.window(xlim=c(0,1), ylim=c(0, 1))
         if(!is.null(leg.vec) & !is.null(leg.col))
             legend('topleft', legend=leg.vec, col=leg.col, pch=pch, pt.cex=max(1, cex.points-1.5), ncol=ifelse( length(leg.vec)> 10, 2, 1), bty='n', cex=2 )
       par(mfrow=c(1,1))
    }

    return(pca)
}

#####################################################
##
## calculate PCs and their variance variance
##
## res - results after testing
## grp - class vector, names are column names of input matrix
#####################################################
my.prcomp2 <- function(res, grp){

  ## remove missing values
  rm.idx <- apply(res, 1, function(x) sum(is.na(x)) + sum(is.infinite(x)))
  rm.idx <- which(rm.idx > 0)
  if(length(rm.idx)>0) res <- res[-rm.idx, ]

  ## extract expression data
  res = res[, names(grp)]

    ## perform pca
  pca <- PCA(scale(t(res)))
 
  return(pca)
}



# 
# #####################################################
# ##
# ##
# ##
# #####################################################
# plotPCAloadings <- function(pca, topn, pca.x, pca.y, pca.z){
# 
#     ##load = loadings(pca)[, c(pc1, pc2, pc3)]
# 
#     ## ###################
#     ## extract loadings
#     load.pca.x <- pca$loadings[, pca.x]
#     load.pca.y <- pca$loadings[, pca.y]
#     load.pca.z <- pca$loadings[, pca.z]
# 
#     n=length(load.pca.x)
# 
#     ## ###################
#     ## choose top N
#     x <- rev(sort(abs( load.pca.x ), decreasing=T )[1:min(topn, n)])
#     y <- rev(sort(abs( load.pca.y ), decreasing=T )[1:min(topn, n)])
#     z <- rev(sort(abs( load.pca.z ), decreasing=T )[1:min(topn, n)])
# 
# 
#     ## ###################################################
#     ## base plotting system
#     par(mfrow=c(1,3))
#     barplot(x, horiz=T, main=paste('PC', pca.x), las=2, border='blue', space=0, col='grey95', ylab='Features', xlab='Absolute coefficient', names.arg=rev(1:length(x)))
#     text(rep(0, length(x)), 1:length(x)-.5, labels=names(x), pos=4)
# 
#     barplot(y, horiz=T, main=paste('PC', pca.y), las=2, border='blue', space=0, col='grey95', axisnames=T, xlab='Absolute coefficient',names.arg=rev(1:length(x)))
#     text(rep(0, length(y)), 1:(length(y))-.5, labels=names(y), pos=4)
# 
#     barplot(z, horiz=T, main=paste('PC', pca.z), las=2, border='blue', space=0, col='grey95', axisnames=T, xlab='Absolute coefficient', names.arg=rev(1:length(x)))
#     text(rep(0, length(z)), 1:(length(z))-.5, labels=names(z), pos=4)
# 
# }
# 
# #####################################################
# ##
# ##            scatterPlotPCAloadings
# ##
# #####################################################
# 
# scatterPlotPCAloadings <- function(pca, topn, pca.x, pca.y, pca.z){
# 
#         ## extract loadings
#         load.pca.x <- pca$loadings[, pca.x]
#         load.pca.y <- pca$loadings[, pca.y]
#         load.pca.z <- pca$loadings[, pca.z]
# 
#         n=length(load.pca.x)
# 
#         ## ###################
#         ## choose top N
#         x <- rev(sort(abs( load.pca.x ), decreasing=T )[1:min(topn, n)])
#         y <- rev(sort(abs( load.pca.y ), decreasing=T )[1:min(topn, n)])
#         z <- rev(sort(abs( load.pca.z ), decreasing=T )[1:min(topn, n)])
# 
# 
#         tmp.pcaloadings <- as.data.frame(pca$loadings[,c(pca.x,pca.y,pca.z)])
#         names(tmp.pcaloadings) <- make.names(names(tmp.pcaloadings))
# 
#         PC1 <- names(tmp.pcaloadings)[pca.x]
#         PC2 <- names(tmp.pcaloadings)[pca.y]
#         PC3 <- names(tmp.pcaloadings)[pca.z]
# 
#         my.scatter <- function(datafr,xa,ya,topx,topy){
# 
#                 #scatterplot dimensions pairwise (pairs passed to function as x,y)
#                 x.axis <- xa ; y.axis <- ya ;topxy <- c(names(topx),names(topy))
#                 # This is the data.frame of the TopN loadings to be marked
#                 mark.frame <- tmp.pcaloadings[which(row.names(tmp.pcaloadings) %in% topxy),]
# 
#                 xmin = min(datafr[,x.axis]);xmax = max(datafr[,x.axis])
#                 ymin = min(datafr[,y.axis]);ymax = max(datafr[,y.axis])
# 
#                 #make the scatterplot
#                 ggplot(data = datafr,aes_string(x = x.axis,y = y.axis))+
#                         geom_hline(yintercept = 0, linetype = "dashed",size =0.7)+
#                         geom_vline(xintercept = 0, linetype = "dashed",size = 0.7)+
#                         geom_point(color = "navy",alpha = 0.15,size=0.1, show.legend = FALSE)+
#                         geom_label_repel(data = mark.frame, size = 3, label.r = unit(0.45,"lines"),
#                                          color = "deeppink", bg = "plum1",
#                                          segment.size = 0.1, box.padding = unit(1,"lines"),
#                                          aes(label = rownames(mark.frame)))+
#                         geom_point(data = mark.frame,
#                                    color = "deeppink",
#                                    size = 2,show.legend = FALSE)+
#                         xlim(c(xmin-0.01,xmax+0.01))+ylim(c(ymin-0.01,ymax+0.01))+
#                         theme_bw()+
#                         theme(panel.grid.major = element_blank(),
#                               panel.grid.minor = element_blank()
#                         )
#                 #mark the TopN ids
# 
#         }
# 
#         # call my.scatter 3 times to make 3 ggplots
#         g1 <- my.scatter(tmp.pcaloadings,PC1,PC2,x,y)
#         g2 <- my.scatter(tmp.pcaloadings,PC2,PC3,y,z)
#         g3 <- my.scatter(tmp.pcaloadings,PC1,PC3,x,z)
#         #combine them in a row by using multiplot() function
#         multiplot(g1,g2,g3,cols = 3)
# }


#################################################
##   Given a string and a number of characters
##   the function chops the string to the
##   specified number of characters and adds
##   '...' to the end.
## parameter
##   string     - character
##   nChar      - numeric
## value
##   string of 'nChar' characters followed
##     by '...'
##################################################
chopString <- function(string, nChar=10, add.dots=T)
{

    string.trim <- strtrim(string, nChar)

    if(add.dots)
        string.trim[ which(nchar(string) > nChar) ] <-  paste(string.trim[which(nchar(string) > nChar) ], '...')
    if(!add.dots)
        string.trim[ which(nchar(string) > nChar) ] <-  paste(string.trim[which(nchar(string) > nChar) ])

    return(string.trim)

}

## ######################################################################
## 20170223
##
##                  Standard deviation filter
##
## ######################################################################
sd.filter <- function(tab, grp.vec, id.col, sd.perc){

    perc <- as.numeric(sd.perc)

    ## extract groups
    groups <- unique(grp.vec)

    ## list to store index of filtered values per group
    ##values.filt <- vector('list', length(groups))
    ##names(values.filt) <- groups

    ## ##########################################
    ## get expression data
    ids=tab[, id.col]
    tab=data.matrix(tab[, names(grp.vec)])

    ## #########################################
    ## calculate sd across all measurements
    sd.tab <- apply(tab, 1, sd, na.rm=T)

    ## #########################################
    ## determine percentile value used to filter
    sd.perc.val <- quantile(sd.tab, sd.perc/100, na.rm=T)

    ## #########################################
    ## index of values to filter
    filt.idx <- which(sd.tab < sd.perc.val)
    not.filt.idx <- which(sd.tab >= sd.perc.val)
    
    ## set filtered values to NA
    tab[filt.idx, ] <- NA

    tab <- data.frame(ids, tab)
    colnames(tab)[1] <- id.col

    values.filt <- lapply(groups, function(x) filt.idx)
    names(values.filt) <- groups
    
    return( 
      list(
        table=tab, 
        values.filtered=values.filt, 
        sd.perc.val=sd.perc.val
      )
    )
}
########################################################################
## 20210429
##                missing data filter
na.filter <- function(tab, id.col, grp.vec, na.filt.val){
  
  ## ##########################################
  ## get expression data
  ids=tab[, id.col]
  tab=data.matrix(tab[, names(grp.vec)])
  
  ## #########################################
  ## calculate sd across all measurements
  na.tab <- apply(tab, 1, function(x) sum(is.na(x))/length(x) )
  
  ## #########################################
  ## index of values to keep
  keep.idx <- which(na.tab <= (na.filt.val/100))
  
  ## filter
  tab <- tab[ keep.idx, ]
  ids <- ids[ keep.idx ]
  
  tab <- data.frame(ids, tab)
  colnames(tab)[1] <- id.col
  
  return(list(
    table=tab, 
    ids=ids
  ))
}

########################################################################
## 20160224
##                   reproducibility filter
##
## n=2: Bland-Altman
## n>2: lmm-model written by DR Mani
##
## - replaces not reproducibly measuered values in 'tab' with 'NA'
## - done separately for each group in 'grp.vec'
##
########################################################################
my.reproducibility.filter <- function(tab, grp.vec, id.col='id', alpha=0.05){

    alpha <- as.numeric(alpha)

    ## extract groups
    groups <- unique(grp.vec)

    ## list to store index of filtered values per group
    values.filt <- vector('list', length(groups))
    names(values.filt) <- groups

    ## add rownames to tab
    ##rownames(tab) <- tab[, id.col]

    ##tab.repro.filter <- tab
    ## View(tab.repro.filter)

    ############################################
    ## loop over replicate groups
    for(gg in groups){

        gg.idx = names(grp.vec)[ which(grp.vec == gg) ]

        ########################################
        ## if there are more than 2 replicates
        ## use the Mani's lmm model
        if( length(gg.idx) > 2 ){
            repro.idx <- reproducibility.filter( tab[, c(id.col, gg.idx)], id.col=id.col, alpha=alpha)

            if(length(repro.idx) != nrow(tab)) stop('Reproducibility vector not of same length as matrix!\n')

            not.repro.idx <- which(!repro.idx)

            if(length(not.repro.idx) > 0)
                tab[not.repro.idx, gg.idx] <- NA

            values.filt[[gg]] <- rownames(tab)[ not.repro.idx ]#not.repro.idx
        }
        ########################################
        ## if there are two replicates use
        ## Blandt-Altmann filter
        ## R-package 'BlandAltmanLeh'
        if( length(gg.idx) == 2 ){

            ## Bland-Altman
            ##ba <-  bland.altman.stats(as.numeric( as.character( tab[, gg.idx[1] ]) ), as.numeric( as.character( tab[,  gg.idx[2] ] )), two=3.290527 )
            ba <-  bland.altman.stats(as.numeric( as.character( tab[, gg.idx[1] ]) ), as.numeric( as.character( tab[,  gg.idx[2] ] )), two=qnorm(1-alpha/2) )
            ## calculate diffs on my own..
            my.diffs <- tab[, gg.idx[1]] - tab[, gg.idx[2]]
            ## index of outliers
            ##not.repro.idx <- which( ba$diffs < ba$lower.limit | ba$diffs > ba$upper.limit)
            not.repro.idx <- which( my.diffs < ba$lower.limit | my.diffs > ba$upper.limit)

            ## set values of outliers to NA
            if(length(not.repro.idx) > 0)
                tab[not.repro.idx, gg.idx] <- NA

            ## store the results
            values.filt[[gg]] <- rownames(tab)[ not.repro.idx ]
            rm(not.repro.idx)
        }

    }
    return(list(table=tab, values.filtered=values.filt))
}


##################################################################
## function to dynamically determine the height (in px) of the heatmap
## depending on the number of genes
dynamicHeightHM <- function(n, ## number of genes to plot in the heatmap 
                            unit=c('px', 'in')){

    unit=match.arg(unit)

    if(is.null(n))
        return(0)

    ## pixel
    if( n < 50)
        height=500
    if( n >= 50 & n <= 100)
        height=800
    if(n >=100)
        height=800+n

    ## inches
    if(unit == 'in'){
        height=height*11/800
    }
    return(height)
}
##################################################################
## determine cell width for heatmap
cwHM <- function(n){
    if(is.null(n))
        return(0)

    cw=55
    if(n < 6) cw=60
    if(n > 10) cw=50
    if(n > 15) cw=30
    if(n > 20) cw=25
    if(n > 30) cw=20
    if(n > 40) cw=15
    if(n > 60) cw=12
    if(n > 80) cw=9
    if(n > 100) cw=6

    if(n > 120) cw=3

    return(cw)
}

##################################################################
## function to dynamically determine the width of the heatmap
## depending on the number of data columns
dynamicWidthHM <- function(n, style=c('none', 'One-sample mod T', 'Two-sample mod T', 'mod F'), unit=c('px', 'in')){

    style=match.arg(style)
    unit=match.arg(unit)

    cw <- cwHM(n)

    width=max(cw * n, 1000)

    if(style == 'One-sample mod T')
        width=width+200

    if(unit == 'in')
        width=width*11/800

    return(width)
}


###################################################################
##
##       generate the profile plots under the 'QC' tab
##
###################################################################
makeProfileplot <- function(tab, id.col, grp, grp.col, grp.col.leg, 
                            legend=T, cex.lab=1.5, mar=c(5,5,3,1), 
                            xlim.mode=c('symmetric', 'as-is', 'as-is (80%)'),
                            plotly=F,
                            main='',
                            ... ){

    cat('\n-- makeProfileplot --\n')

    ## table
    tab <- tab[, setdiff(colnames(tab), id.col)]

    ## calculate densities
    dens <- apply(tab, 2, density, na.rm=T)
   
    ## ylim
    ylim <- max(unlist(lapply(dens, function(x) max(x$y))))

    ## xlim
    xlim.mode <- match.arg(xlim.mode)
    if(xlim.mode == 'symmetric'){
      xlim=max(abs(tab), na.rm=T)
      xlim <- c(-xlim, xlim)
    }
    if(xlim.mode == 'as-is'){
      xlim_min <- min(unlist(lapply(dens, function(x) min(x$x))))
      xlim_max <- max(unlist(lapply(dens, function(x) max(x$x))))
      xlim <- c(xlim_min, xlim_max)
    }
    if(xlim.mode == 'as-is (80%)'){
      x_all <- unlist(lapply(dens, function(x) x$x))
      x_robust_quant <- quantile(x_all, c(0.1, 0.9))
      x_robust <- x_all[x_all > x_robust_quant[1] & x_all < x_robust_quant[2] ]
      
      xlim_min <- min(x_robust)
      xlim_max <- max(x_robust)
      xlim <- c(xlim_min, xlim_max)
    }
    ##########################################
    ## plot
    if(!plotly){
      p <- NULL
      par(mar=mar)
      for(i in 1:ncol(tab)){
          if(i == 1)
              plot(dens[[i]], xlab='expression', xlim=xlim, ylim=c(0, ylim), col=my.col2rgb(grp.col[i], alpha=100), lwd=3, cex.axis=2, cex.lab=2, cex.main=1.5, main=main, ...)
          else
              lines(dens[[i]], col=my.col2rgb(grp.col[i], alpha=100), lwd=3)
  
          ## divide legend if there are too many experiments
          N.exp <- length(names(grp.col.leg))
          if( N.exp > 15){
              legend('topright', legend=names(grp.col.leg)[1:floor(N.exp/2)], col=grp.col.leg[1:floor(N.exp/2)], lty='solid', bty='n', cex=1.5, lwd=3)
              legend('topleft', legend=names(grp.col.leg)[ceiling(N.exp/2):N.exp], col=grp.col.leg[ceiling(N.exp/2):N.exp], lty='solid', bty='n', cex=1.5, lwd=3)
  
          } else
              legend('topright', legend=names(grp.col.leg), col=grp.col.leg, lty='solid', bty='n', cex=1.5, lwd=3)
      }
    } else { ## plotly
     
      grp.unique <- names(grp.col.leg)
      #save(dens, grp, grp.unique, grp.col.leg, file='debug.RData')
      for(g in grp.unique) {
        
        grp_tmp <- names(grp)[grp == g]
        
        x_tmp <- lapply(dens[grp_tmp], function(d) d$x) %>% unlist
        y_tmp <- lapply(dens[grp_tmp], function(d) d$y) %>% unlist
        label_tmp <- lapply(grp_tmp, function(l) rep(l, length(dens[[l]]$x)) ) %>% unlist
        
        if(g == grp.unique[1]) {
          p <- plot_ly(x=x_tmp, y=y_tmp, mode='lines', name=g, type='scatter', 
                     color = I(grp.col.leg[g]), text=label_tmp ) 
        } else {
          p <- p %>% add_trace(x=x_tmp, y=y_tmp, mode='lines', name=g,
                               color = I(grp.col.leg[g]), text=label_tmp)
        }
      }
      p <- p %>% layout(title = main, 
                        xaxis = list(title='expression', range=xlim),
                        yaxis = list(title='Density'))
    }
    cat('\n-- makeProfileplot exit --\n')
    return(p)
}


##############################################
# Multiple plot function
##############################################
###################################################################################
# Multiple plot function
# Source: R-cookbook
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
###################################################################################
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
        p_load(grid)

        # Make a list from the ... arguments and plotlist
        plots <- c(list(...), plotlist)

        numPlots = length(plots)

        # If layout is NULL, then use 'cols' to determine layout
        if (is.null(layout)) {
                # Make the panel
                # ncol: Number of columns of plots
                # nrow: Number of rows needed, calculated from # of cols
                layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                                 ncol = cols, nrow = ceiling(numPlots/cols))
        }

        if (numPlots==1) {
                print(plots[[1]])

        } else {
                # Set up the page
                grid.newpage()
                pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

                # Make each plot, in the correct location
                for (i in 1:numPlots) {
                        # Get the i,j matrix positions of the regions that contain this subplot
                        matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

                        print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                                        layout.pos.col = matchidx$col))
                }
        }


}


## #################################################################################
## ppi.bait - character
## db   - character, ppi databases
## IDs  - vector of ids in dataset
## sig.idx - index if significnat features in the dataset
##
## ppi.db      - list of ppi databases, defined in 'global.r'
## ppi.db.col  - colors for different ppi databases, defined in 'global.r'
get.interactors <- function(ppi.bait, IDs, sig.idx, db=c('iw', 'bg', 'react'), ppi.db, ppi.db.col){
    
    # ######################################
    # extract gene symbols
    # last string after underscore
    IDs <- toupper(sub('.*_(.*)$', '\\1', IDs) )
    names(IDs) <- IDs
    
    # gene symbol of bait protein
    ppi.bait <-  toupper(sub('.*_(.*)$', '\\1', ppi.bait))
    
    ## ###################################################
    ## index of bait protein
    ppi.bait.idx <- which( toupper(IDs) == toupper(ppi.bait) ) ## bait in data set
    
  
    ## ##############################################
    ## return if no database has been selected
    if(length(db) == 0){
        ppi.int.idx <- NULL
        ppi.col <- rep('', length(IDs))
        leg=''
        ppi.col.leg.sig.tmp='white'
        ppi.db.col.tmp='white'
    
        ## ##########################
        ## assemble output
        out <- list(
          ppi.int.idx=ppi.int.idx,
          ppi.bait.idx=ppi.bait.idx,
          leg=leg,
          leg.col=ppi.db.col.tmp,
          ppi.col=ppi.col
        )
        return(out)
    }

    ## #########################################################
    ## list to store ALL interactors of 'ppi.bait' found in the
    ## selected PPi databases
    ppi.int.all.l <- vector('list', length(db))
    names(ppi.int.all.l) <- db

    ## #########################################################
    ## list to store DETECTED interactors
    ppi.int.detect.l <- ppi.int.all.l

    ## #########################################################
    ## list to store DETECTED and SIGNIFICANT interactors
    ppi.int.signif.l <- ppi.int.all.l

    ## #########################################
    ##             InWeb
    ## #########################################
    if( 'iw' %in% db ){

        iw <- ppi$iw

        ## ##################################
        ## try gene names
        ppi.bait.gn <- ppi.bait

        i1.gn <- iw$V5
        i2.gn <- iw$V6
        
        iw.int <- i2.gn[ which(i1.gn == ppi.bait.gn) ]
        iw.int <- c(iw.int, i1.gn[which(i2.gn == ppi.bait.gn)])

        iw.int <- setdiff( unique(iw.int), ppi.bait )

        ## store all interactors
        ppi.int.all.l[['iw']] <- iw.int
        names(ppi.int.all.l)[which(names(ppi.int.all.l) == 'iw')] <- 'InWeb'

        ## interactors found in dataset
        ppi.int.detect.l[['iw']] <- intersect( IDs, iw.int)
        names(ppi.int.detect.l)[which(names(ppi.int.detect.l) == 'iw')] <- 'InWeb'

        ## significant interactors
        ppi.int.signif.l[['iw']] <- intersect( IDs[ sig.idx ], iw.int)

    }
    ## #################################################
    ##                BioGRID
    ## #################################################
    if( 'bg' %in% db ){

        bg <- ppi$bg

        ##ppi.bait.gn <- toupper(sub('.*_', '', ppi.bait))
        ppi.bait.gn <- ppi.bait

        i1.gn <- bg$Alt.IDs.Interactor.A
        i2.gn <- bg$Alt.IDs.Interactor.B

        bg.int <- i2.gn[which(i1.gn == ppi.bait.gn)]
        bg.int <- c(bg.int, i1.gn[which(i2.gn == ppi.bait.gn)])

        bg.int <- setdiff( unique(bg.int), ppi.bait)

        ## all interactors
        ppi.int.all.l[['bg']] <- bg.int
        names(ppi.int.all.l)[ which(names(ppi.int.all.l) == 'bg')] <- 'BioGRID'

        ## interactors found in dataset
        ppi.int.detect.l[['bg']] <- intersect( IDs, bg.int)
        names(ppi.int.detect.l)[which(names(ppi.int.detect.l) == 'bg')] <- 'BioGRID'

        ## significant interactors
        ppi.int.signif.l[['bg']] <- intersect( IDs[ sig.idx ], bg.int)

    }
    ## #################################################
    ##              Reactome
    ## #################################################
    if( 'react' %in% db ){
        react <- ppi$react


        ppi.bait.gn <- ppi.bait

        i1.gn <- react$alternative.id.A
        i2.gn <- react$alternative.id.B

        react.int <- i2.gn[which(i1.gn == ppi.bait.gn)]
        react.int <- c(react.int, i1.gn[which(i2.gn == ppi.bait.gn)])

        ## make unique
        react.int <-setdiff(  unique(react.int), ppi.bait )

        ## all interactors
        ppi.int.all.l[['react']] <- react.int
        names(ppi.int.all.l)[which(names(ppi.int.all.l) == 'react')] <- 'Reactome'

        ## interactors found in dataset
        ppi.int.detect.l[['react']] <- intersect( IDs, react.int)
        names(ppi.int.detect.l)[which(names(ppi.int.detect.l) == 'react')] <- 'Reactome'

        ## significant interactors
        ppi.int.signif.l[['react']] <- intersect( IDs[ sig.idx ], react.int)

    }

    ## #################################################
    ##
    ## - Interactions found in multiple ppi databases
    ## #################################################
    if(length(db) > 1){
        ## all interactors
        ppi.int.all.l.ol <- table(unlist(ppi.int.all.l))
        ppi.int.all.l.ol <- names(ppi.int.all.l.ol[ ppi.int.all.l.ol > 1 ])

        ppi.int.all.l <- append(ppi.int.all.l, list(ppi.int.all.l.ol))
        names(ppi.int.all.l)[length(ppi.int.all.l)] <- 'Shared'

        ## detected interactors
        ppi.int.detect.l.ol <- table(unlist(ppi.int.detect.l))
        ppi.int.detect.l.ol <- names(ppi.int.detect.l.ol[ppi.int.detect.l.ol > 1])

        ppi.int.detect.l <- append(ppi.int.detect.l, list(ppi.int.detect.l.ol))
        names(ppi.int.detect.l)[length(ppi.int.detect.l)] <- 'Shared'

        ## significant interactors
        ppi.int.signif.l.ol <- table(unlist(ppi.int.signif.l))
        ppi.int.signif.l.ol <- names(ppi.int.signif.l.ol[ppi.int.signif.l.ol > 1])
        ppi.int.signif.l <- append(ppi.int.signif.l, list(ppi.int.signif.l.ol))
        names(ppi.int.signif.l)[length(ppi.int.signif.l)] <- 'Shared'
    }
    ##save(ppi.int.all.l, ppi.int.detect.l, file='tmp.RData')

    ## ##################################################
    ##            combine
    ## ##################################################
    ## all interactions
    ppi.int <- unique(unlist(ppi.int.all.l))
    #ppi.int <- unlist(ppi.int.all.l)
    
    ppi.int.idx <- NULL
    leg <- NULL
    leg.col <- NULL

    ## ##################################################
    ## if there are interactors
    if(length(ppi.int) > 0){

        ## ###############################
        ## index of interactors in dataset
        ppi.int.idx <- which( IDs %in% unlist(ppi.int.detect.l) )

        ## ######################################################################
        ## check ppi source: different colors
        ## - consider overlap as separate class
        ## - counts based on this vector are mutually axclusive, i.e. to count
        ##   number of interactions in InWeb one has to sum #InWeb + #Overlap
        ppi.col <- rep('', length(IDs))
        names(ppi.col) <- names(IDs)
        
        for( d in 1:length(ppi.int.detect.l))
          ppi.col[ which( IDs %in% ppi.int.detect.l[[d]] )] <- ppi.db.col[ names(ppi.int.detect.l)[d]  ]


        ## #########################################
        ## colors, defined at the beginning of
        ## this file
        ## #########################################
        ppi.db.col.tmp <- ppi.db.col[ names(ppi.int.all.l) ]


        ## #########################################
        ## numbers for the legend
        leg.all <- sapply(ppi.int.all.l, length)
        leg.detect <- sapply(ppi.int.detect.l, length)
        leg.signif <- sapply(ppi.int.signif.l, length)

        leg <- paste(names(ppi.int.all.l), ' ',leg.signif, '/' ,leg.detect, '/' ,leg.all, sep='')
        #leg <- c(paste('bait:',ppi.bait), leg)
        ##leg <- paste(names(ppi.col.leg), ' (',ppi.col.leg.sig,'/' ,ppi.col.leg,')', sep='')
        ##legend('topleft', legend=leg, col=ppi.db.col.tmp, pch=16, bty='n', cex=1.5, title=paste('Known interactors (sig/tot)'))

    } else {  ## end if there are interactors
        ppi.col <- rep('', length(IDs))
        leg=''
        ppi.col.leg.sig.tmp='white'
        ppi.db.col.tmp='white'
    }

    ## ##############################################
    ## return if no database has been selected
    if(length(db) == 0 | (length(ppi.int.idx) == 0 &  length(ppi.bait.idx) == 0)){
      ppi.int.idx <- NULL
      ppi.col <- rep('', length(IDs))
      leg=''
      ppi.col.leg.sig.tmp='white'
      ppi.db.col.tmp='white'
      
      ## assemble output
      out <- list(
        ppi.int.idx=ppi.int.idx,
        ppi.bait.idx=ppi.bait.idx,
        leg=leg,
        leg.col=ppi.db.col.tmp,
        ppi.col=ppi.col
      )
    } else if(length(ppi.int.idx) == 0 & length(ppi.bait.idx) == 0){
      
      ## ##########################
      ## if neither bait nor interactors were found ...
      ppi.int.idx <- NULL
      ppi.bait.idx <- NULL
      ppi.col <- rep('', length(IDs))
      leg=''
      ppi.col.leg.sig.tmp='white'
      ppi.db.col.tmp='white'
      
      ## assemble output
      out <- list(
        ppi.int.idx=ppi.int.idx,
        ppi.bait.idx=ppi.bait.idx,
        leg=leg,
        leg.col=ppi.db.col.tmp,
        ppi.col=ppi.col
      )
    } else { 
      
      ## ##########################
      ## all normal
      out <- list(
        ppi.int.idx=ppi.int.idx,
        ppi.bait.idx=ppi.bait.idx,
        leg=leg,
        leg.col=ppi.db.col.tmp,
        ppi.col=ppi.col
      )
      
    }
    return(out)
}

# ##########################################################
# Function to plot color bar
# https://stackoverflow.com/questions/9314658/colorbar-from-custom-colorramppalette
# example: 
# color.bar(colorRampPalette(c("light green", "yellow", "orange", "red"))(100), -1)

color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
  scale = (length(lut)-1)/(max-min)
  
  dev.new(width=1.75, height=5)
  plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
  axis(2, ticks, las=1)
  for (i in 1:(length(lut)-1)) {
    y = (i-1)/scale + min
    rect(0,y,10,y+1/scale, col=lut[i], border=NA)
  }
}

#############################################################
## create and export xlsx file
##
export2xlsx <- function(res.comb, grp, grp.comp, rdesc, which.test, headerDesc=NULL, fn, session.dir){

  grp.sorted <- sort(grp)
  
  ## append annotation columns
  if(!is.null(rdesc)){
    res.comb <- left_join(res.comb,  rdesc, 'id')
  }
  
  ## #########################################################
  ## Two sample moderated T-test:
  ## adjust labels in the header of the Excel sheet
  ## WT.vs.KO -> KO.over.WT
  if(which.test == 'Two-sample mod T'){
    
    colnames.tmp <- colnames(res.comb)
    grp.comp <- unique(grp.comp)
    
    for(g in grp.comp){
      g.new <- strsplit(g, '\\.vs\\.') %>% unlist
      g.new <- paste(g.new[2], '.over.', g.new[1], sep='')
      colnames.tmp <- gsub(paste0(g,'$'), g.new, colnames.tmp)
    }
    colnames(res.comb) <- colnames.tmp
  }
  
  ############################################################
  ## experimental design
  expDesign <- data.frame(Column=names(grp.sorted), Experiment=grp.sorted)
  
  ##########################################################
  ## add column descriptions
  ## 
  if(!is.null(headerDesc )){
    
    ## descriptions available in the in the "library" file
    ColumnHeaderLibrary <- headerDesc$ColumnHeader
    DescriptionLibrary <- headerDesc$Description
    SourceLibrary <- headerDesc$Source
    names(ColumnHeaderLibrary) <- names(DescriptionLibrary) <- names(SourceLibrary) <- headerDesc$ColumnHeader
    
    ## result table to annotate
    ColumnHeader <- colnames(res.comb)
    Description <- Source <- rep(NA, length(ColumnHeader))
    names(Description) <- names(Source) <- names(ColumnHeader) <- ColumnHeader
    
    ## to store index of already matched column names
    matched <- c()
    
    ## expression columns used
    for(cc in 1:nrow(expDesign)) {
      match.idx <- which(ColumnHeader == expDesign[cc, 'Column'])
      
      Description[ match.idx ] <- paste0('Abundance in group "', expDesign[cc, 'Experiment'],'"') 
      Source[ match.idx ] <- 'Protigy'
      
      matched <- c( matched, match.idx)
    }
    
    ############################    
    ## loop over terms in the 'library' file
    for(cc in  ColumnHeaderLibrary){
      #if(cc == 'logFC') browser()
      ## exact match
      match.idx <- which( ColumnHeader == cc)
      
      if(length(match.idx) == 1){

        ## if hasn't been matched ...
        if( !(match.idx %in% matched) ){
          if(SourceLibrary[cc] == 'limma'){
            Description[ match.idx ] <- paste0(DescriptionLibrary[cc], ' "', which.test, '"')
          } else {
            Description[ match.idx ] <- paste0(DescriptionLibrary[cc])
          }
        
          Source[ match.idx ] <- SourceLibrary[ cc ]
        
          matched <- c( matched, match.idx)
        } ## end if already matched
        
      } ##else { 
      
      ## no exact match
     else if(length(match.idx) == 0) {

        ## test prefix: protigy/limma columns
        match.idx <- grep( paste0('^',cc,'\\..*'), ColumnHeader)
        
        match.idx <- match.idx[ which(!match.idx %in% matched) ]
        
        
        if(length(match.idx) > 0 ){
        
          
          ## if hasn't been matched ...
         # if( !(match.idx %in% matched) ) {
              if(SourceLibrary[cc] == 'limma'){
                Description[ match.idx ] <- paste0(DescriptionLibrary[cc], ' "', which.test, '" in group "',   sub(paste0('^', cc, '.'), '', ColumnHeader[match.idx]), '"')
             } else {
                Description[ match.idx ] <- paste0(DescriptionLibrary[cc], ' "',sub(paste0('^', cc, '.'), '', ColumnHeader[match.idx]), '"')
             }
              Source[ match.idx ] <- SourceLibrary[ cc ]
              matched <- c( matched, match.idx)
          #}  ## end if already matched
          
        } else { ## end test prefix
          
          ## text suffix
          match.idx <- grep( paste0('.*\\.', cc,'$'), ColumnHeader)
          if(length(match.idx) > 0 ){
            
            ## if hasn't been matched ...
            if( !(match.idx %in% matched) ){
                if(SourceLibrary[cc] == 'SM'){
                  Description[ match.idx ] <- paste0(DescriptionLibrary[cc], ' in directory "',   sub(paste0('.', cc, '$'), '', ColumnHeader[match.idx]), '"')
                } else {
                  Description[ match.idx ] <- paste0(DescriptionLibrary[cc], ' "',sub(paste0('^', cc, '.'), '', ColumnHeader[match.idx]), '"')
                }
                Source[ match.idx ] <- SourceLibrary[ cc ]
                matched <- c( matched, match.idx)
            }  ## end if already matched
          }
        } ## end test suffix 
      } ## end no extact match
    } ## end loop over ColumnHeaderLibrary
    
    allColumnsInTable <- data.frame(ColumnHeader, Description, Source)
  } 
  
  #############################
  ## export
  render.xlsx <- try(
     WriteXLS(c('res.comb', 'expDesign', 'allColumnsInTable'), ExcelFileName=paste(session.dir, fn, sep='/'), 
              FreezeRow=1, FreezeCol=1, 
              SheetNames=c(which.test, 'Class vector', 'Description of table header'), 
              row.names=F, BoldHeaderRow=T, AutoFilter=T)
  )
  #browser()
  ## error checking
  if(class(render.xlsx) == 'try-error' | !render.xlsx){
    showModal(modalDialog(
      size='m',
      title = "Problem generating Excel sheet",
      HTML(paste('Error msg:', render.xlsx[[1]], '<br><br>')),
      HTML('If you are a Broadie and encounter this problem using Protigy on RStudio Connect please reach out to your administrator.<br><br>'),
      HTML('If you downloaded <a href="https://github.com/broadinstitute/protigy" target="_blank_">Protigy from GitHub</a> and running it on your local computer, one possible reason could be that Perl is not installed on your computer.<br><br>'), 
      HTML('For Windows OS you can install Strawberry Perl (<a href="http://strawberryperl.com/" target="_blank_">http://strawberryperl.com/</a>).<br>Please note that R needs to be restarted after installing Perl.)')
      
         ))
    }
    return(render.xlsx)
}


###############################################################
## filename for result tables
## .xlsx
## .gct
##
create.fn <- function(label, which.test, log.transform, repro.filt, filt.data, suffix=c('xlsx', 'gct')){
  suffix <- match.arg(suffix)
  fn <- sub(' ', '_',
                paste0(
                  label, '_',
                  sub(' ', '_', which.test),
                  ifelse(log.transform != 'none', paste( '_', log.transform, '_', sep=''), '_'),
                  ifelse(repro.filt=='yes', paste(filt.data, sep=''), '_'),
                  sub(' .*', '', Sys.time()), 
                ".", suffix) 
  )
  return(fn)
}


#######################################################
##
## returns the index of significant features in "group" in the 
## UNFILTERED dataset
##
get_sig_features <- function(global.results, global.param, group){
  
  ## unfiltered data set
  res = as.data.frame( global.results$data$output )
  
  ## #############################
  ## one-sample T or two sample T
  if(global.param$which.test %in% c('One-sample mod T', 'Two-sample mod T')){
    
    if(global.param$filter.type == 'adj.p')
      pval <- res[, paste('adj.P.Val.', group, sep='')]
    else if(global.param$filter.type == 'nom.p')
      pval <- res[, paste('P.Value.', group, sep='')]
    else
      pval <- rep(1, nrow(res))
  }
  ## ###############################
  ## p-values mod F
  if(global.param$which.test == 'mod F') {
    
    if(global.param$filter.type == 'adj.p')
      pval <- res[, paste('adj.P.Val', sep='')]
    else if(global.param$filter.type == 'nom.p')
      pval <- res[, paste('P.Value', sep='')]
    else
      pval <- rep(1, nrow(res))
  }
 
  ## ###############################
  ## p-values no test
  if(global.param$which.test == 'none') {
    pval <- rep(1, nrow(res))
  }
  
  ##################################
  ## indices of significant features
  sig.idx <- which(pval < as.numeric(global.param$filter.value ))
  
  return(sig.idx)
}



