################################################################################################################
## Filename: global.r
## Created: October 09, 2015
## Author(s): Karsten Krug, Ozan Aygun
##
## Purpose: Shiny-app to perform differential expression analysis, primarily on proteomics data, to perform
##          simple data QC, to interactively browse through the results and to download high-quality result
##          figures.
##
## This file defines global parameters, loads all required R-packages and defines the actual functions to perform data filtering,
## data normalization, the moderated test statistics, and visualization. The code for moderated t-tests,
## two-component normalization and the reproducibility filter has been written by Mani DR and
## adopted by me for intergration into a Shiny-Server environment.
##
##
## required packages:
##
## cran.pckg <- c('pheatmap', 'RColorBrewer', 'hexbin', 'Hmisc', 'grid', 'scatterplot3d', 'plotly', 'WriteXLS', 'reshape','nlme', 'BlandAltmanLeh', 'mice','mixtools', 'mclust')
## bioc.pgkg <- c( 'preprocessCore', 'limma')
##
## changelog: 20160614 - included 'na' to indicate missing values
##                     - outsourced Mani's code to a separate file 'modT.r'
################################################################################################################

## R package managing tool
if (!require("pacman")) install.packages ("pacman")
require('pacman')
p_load (RColorBrewer)

#################################################################
## global parameters
#################################################################
## version number
VER="0.8.1"
## maximal filesize for upload
MAXSIZEMB <<- 500
## list of strings indicating missing data
NASTRINGS <<- c("NA", "<NA>", "#N/A", "#NUM!", "#DIV/0!", "#NA", "#NAME?", "na", "#VALUE!")
## speparator tested in the uploaded file
SEPARATOR <<- c('\t', ',', ';')
## Colors used throughout the app to color the defined groups
GRPCOLORS <<- c(RColorBrewer::brewer.pal(8, "Set1"), RColorBrewer::brewer.pal(8, "Dark2"), RColorBrewer::brewer.pal(8, "Set2"), terrain.colors(20), cm.colors(20), topo.colors(20))
## number of characters to display in plots/tables for column names
STRLENGTH <<- 20
## operating system
OS <<- Sys.info()['sysname']
## temp directory to write the Excel file
TMPDIR <<- ifelse(OS=='Windows', "./", "/tmp/")
## app name
##APPNAME <<- sub('.*/','',getwd())
APPNAME <<- 'Protigy'
## app folder
APPDIR <<- getwd()
## directory to store data files
##DATADIR <<- ifelse(OS=='Windows', ".", "/local/shiny-data/")
DATADIR <<- ifelse(OS=='Linux', "/local/shiny-data/", '.')

## email for trouble shooting
MAIL <<- 'karsten@broadinstitute.org'
## URL to configuration app (SSP only)
CONFAPP <<- 'http://shiny-proteomics.broadinstitute.org:3838/modTconf/'
## PIWIK location
PIWIKURL <<- '//shiny-proteomics.broadinstitute.org/piwik/'

#################################################################
## load required packages
#################################################################

p_load(shiny)
p_load(shinydashboard)
p_load(shinyjs)
## colors

## heatmap
##p_load(pheatmap)
p_load(heatmaply)

# clustering
p_load(ape)
p_load(dendextend)
p_load(scales)
p_load(gtable)
p_load(fastcluster)

# correlations
#p_load(GO.db) ## reguired by WGCNA
#p_load(impute) ## reguired by WGCNA
#p_load(WGCNA)

# markdown reports
p_load(rmarkdown)
p_load(knitr)

## moderated tests
p_load(limma)
p_load(statmod)

## multiscatter
p_load(hexbin)
p_load(Hmisc)
p_load(grid)
## pca
p_load(ChemometricsWithR)
p_load(scatterplot3d)
p_load(plotly)
## export
p_load(WriteXLS)
## reproducibility filter
p_load(reshape)
p_load(nlme)
p_load(BlandAltmanLeh)
## normalization Quantile
p_load(preprocessCore)
## normalization 2-component
p_load (mice)
p_load (mixtools)
p_load (mclust)
## table preview
p_load(DT)
## label placements without overlap
p_load(maptools)
p_load(ggrepel)
p_load(dplyr)

## id mapping
p_load(RSQLite)
p_load(org.Hs.eg.db)
p_load(org.Mm.eg.db)
p_load(org.Rn.eg.db)
p_load(org.Dr.eg.db)


## morpheus
#p_load_gh('morpheus')
#if(!require(morpheus))
#  devtools::install_github('cmap/morpheus.R')
#p_load(morpheus)

# Required for cmpaR gctx file format. Fails to install on shiny-proteomics, but not reuqired as of now.
#p_load (rhdf5) 



source('src/modT.r')
source('src/pheatmap.r')
source('src/helptext.r')
source('src/gct-io.r')
source('src/plots.r')


## #####################################
## CSS for loading animantion
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
mapIDs <- function(ids,
                   n.try=10
                   ){
    withProgress(message='Mapping gene names...', {

            ## ###################################
            ##           id type
            ## ###################################
            keytype <- 'UNKNOWN'
            ## Uniprot or RefSeq?
            if(length(grep('^(Q|P|O|A|E|H|F)', ids)) > 0)
                keytype='UNIPROT'
            if(length(grep('^(NP_|XP_|YP_)', ids)) > 0)
                keytype='REFSEQ'

            ## ###################################
            ##        extract query strings
            ## ###################################
            if(keytype == 'UNIPROT' ){
              id.query <- sub('(-|;|\\.|_|\\|).*', '', ids) ## first id
            } else if(keytype == 'REFSEQ') {
              id.query <- sub('(\\.|;).*', '', ids) ## first id
            } else {
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
                id.map.tmp <- try( mapIds(org.Hs.eg.db, keys=id.query[ sample(1:length(id.query), n.try)] , column=c('SYMBOL'), keytype=keytype, multiVals='first') )
                if(class(id.map.tmp) != 'try-error'){
                  orgtype='HSA'
                }
              }
              # try mouse
              if(orgtype == 'UNKNOWN'){ 
                id.map.tmp <- try( mapIds(org.Mm.eg.db, keys=id.query[ sample(1:length(id.query), n.try)] , column=c('SYMBOL'), keytype=keytype, multiVals='first') )
                  if(class(id.map.tmp) != 'try-error'){
                    orgtype='MMU'
                  }
              }
              # try rat
              if(orgtype == 'UNKNOWN'){ 
                id.map.tmp <- try( mapIds(org.Rn.eg.db, keys=id.query[ sample(1:length(id.query), n.try)] , column=c('SYMBOL'), keytype=keytype, multiVals='first') )
                if(class(id.map.tmp) != 'try-error'){
                  orgtype='RNO'
                }
              }
              # try zebrafish
              if(orgtype == 'UNKNOWN'){ 
                id.map.tmp <- try( mapIds(org.Dr.eg.db, keys=id.query[ sample(1:length(id.query), n.try)] , column=c('SYMBOL'), keytype=keytype, multiVals='first') )
                if(class(id.map.tmp) != 'try-error'){
                  orgtype='DRE'
                }
              }
            }
  
            ## ##################################
            ## map
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

            if(class(id.map.tmp) == 'try-error' | is.null( class(id.map.tmp) ) | class(id.map.tmp) == 'NULL' ){

              id.map <- data.frame(id=names(id.query), id.query=id.query, id.mapped=names(id.query), id.concat=ids, stringsAsFactors=F)
              #cat('test2\n')
              keytype <- 'UNKNOWN'
              
            } else {
    
              id.map.tmp[which(is.na(id.map.tmp))] <- 'NotFound'
              id.map <- data.frame(id=names(id.query), id.query=id.query, id.mapped=id.map.tmp, id.concat=paste(ids, id.map.tmp, sep='_'), stringsAsFactors=F)
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

## ##############################################################################
##
##            generate links to external databases
##
## ##############################################################################
link.db <- function(id, # vetcor of ids
                    keytype=c('UNKNOWN', 'UNIPROT', 'REFSEQ'),
                    db=c('GENECARDS', 'UNIPROT')){
  
  keytype <- match.arg(keytype)
  db <- match.arg(db)
  
  if(keytype == 'UNIPROT'){
    up.link <- paste("<a href='https://www.uniprot.org/uniprot/", sub('(_|,|;|\\.).*', '', id),"' target='_blank'>", id, "</a>", sep='')
  }
  if(keytype %in% c('REFSEQ', 'UNKNOWN')){
    up.link <- paste("<a href='http://www.genecards.org/Search/Keyword?queryString=", sub('^(NP_|NM_|NR_.*?)(_|,|;|\\.).*', '\\1', id),"' target='_blank'>", id, "</a>", sep='')
  }
  return(up.link)
}





#############################################################################################
##
##              different normalization methods for expression data
##
## 20160235
#############################################################################################
normalize.data <- function(data, id.col, method=c('Median', 'Quantile', 'Median-MAD', '2-component')){
    cat('\n\n-- normalize data --\n\n')

    method = match.arg(method)

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
        ##rownames(data.norm) <- rownames(data)
        colnames(data.norm) <- paste( colnames(data), sep='.')
    }
    ## median & MAD
    if(method == 'Median-MAD'){
        data.norm <- apply(data, 2, function(x) (x - median(x, na.rm=T))/mad(x, na.rm=T) )
        ##rownames(data.norm) <- rownames(data)
        colnames(data.norm) <- paste( colnames(data), sep='.')
    }
    ## 2-component normalization
    if(method == '2-component'){
        data.norm.list = apply(data, 2, two.comp.normalize, type="unimodal")
        ##cat('\n\n here 1\n\nlength list:', length(data.norm.list), '\n')
        ##save(data.norm.list, file='test.RData')
        ## check if successful
        for(i in 1:length(data.norm.list)){
            ##cat('\ni=', i, '\n')
            if(length(data.norm.list[[i]]) == 1){
                if(data.norm.list[[i]] == 'No_success'){
                    ##cat('\n\nno success\n\n')
                    return(paste( colnames(data)[i] ))
                }
            }
            ##cat('\n length:', length(data.norm.list[[i]]), '\n')
        }
        ##cat('\n\n here \n\n')
        data.norm = matrix( unlist(lapply(data.norm.list, function(x)x$norm.sample)), ncol=length(data.norm.list), dimnames=list(rownames(data), names(data.norm.list)) )
    }
    ## add id column
    data.norm <- data.frame(ids, data.norm)
    colnames(data.norm)[1] <- id.col
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
my.prcomp <- function(x, pca.x, pca.y, pca.z, col=NULL, cor=T, plot=T, rgl=F, scale=T, pch=20, cex.points=3, rgl.point.size=30, main="PCA", leg.vec=NULL, leg.col=NULL, ...){

    cex.font = 1.8

    ##View(x)

    ## number of data columns, N=2 -> 2D plot only
    N <- nrow(x)

    ## color
    if( is.null(col) ) col="black"

    ## perform pca
    pca <- prcomp(x, scale=scale)

    ##View(pca$x)

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
    ##res <- data.matrix(res)
    rm.idx <- apply(res, 1, function(x) sum(is.na(x)) + sum(is.infinite(x)))
    rm.idx <- which(rm.idx > 0)
    if(length(rm.idx)>0) res <- res[-rm.idx, ]

  ## extract expression data
  res = res[, names(grp)]

   ## View(res)

  ## perform pca
  ##pca <- prcomp(x, scale=scale)
  pca <- PCA(scale(t(res)))


  return(pca)
}

#####################################################
##
##
##
#####################################################
plotPCAloadings <- function(pca, topn, pca.x, pca.y, pca.z){

    ##load = loadings(pca)[, c(pc1, pc2, pc3)]

    ## ###################
    ## extract loadings
    load.pca.x <- pca$loadings[, pca.x]
    load.pca.y <- pca$loadings[, pca.y]
    load.pca.z <- pca$loadings[, pca.z]

    n=length(load.pca.x)

    ## ###################
    ## choose top N
    x <- rev(sort(abs( load.pca.x ), decreasing=T )[1:min(topn, n)])
    y <- rev(sort(abs( load.pca.y ), decreasing=T )[1:min(topn, n)])
    z <- rev(sort(abs( load.pca.z ), decreasing=T )[1:min(topn, n)])


    ## ###################################################
    ## base plotting system
    par(mfrow=c(1,3))
    barplot(x, horiz=T, main=paste('PC', pca.x), las=2, border='blue', space=0, col='grey95', ylab='Features', xlab='Absolute coefficient', names.arg=rev(1:length(x)))
    text(rep(0, length(x)), 1:length(x)-.5, labels=names(x), pos=4)

    barplot(y, horiz=T, main=paste('PC', pca.y), las=2, border='blue', space=0, col='grey95', axisnames=T, xlab='Absolute coefficient',names.arg=rev(1:length(x)))
    text(rep(0, length(y)), 1:(length(y))-.5, labels=names(y), pos=4)

    barplot(z, horiz=T, main=paste('PC', pca.z), las=2, border='blue', space=0, col='grey95', axisnames=T, xlab='Absolute coefficient', names.arg=rev(1:length(x)))
    text(rep(0, length(z)), 1:(length(z))-.5, labels=names(z), pos=4)

}

#####################################################
##
##            scatterPlotPCAloadings
##
#####################################################

scatterPlotPCAloadings <- function(pca, topn, pca.x, pca.y, pca.z){

        ## extract loadings
        load.pca.x <- pca$loadings[, pca.x]
        load.pca.y <- pca$loadings[, pca.y]
        load.pca.z <- pca$loadings[, pca.z]

        n=length(load.pca.x)

        ## ###################
        ## choose top N
        x <- rev(sort(abs( load.pca.x ), decreasing=T )[1:min(topn, n)])
        y <- rev(sort(abs( load.pca.y ), decreasing=T )[1:min(topn, n)])
        z <- rev(sort(abs( load.pca.z ), decreasing=T )[1:min(topn, n)])


        tmp.pcaloadings <- as.data.frame(pca$loadings[,c(pca.x,pca.y,pca.z)])
        names(tmp.pcaloadings) <- make.names(names(tmp.pcaloadings))

        PC1 <- names(tmp.pcaloadings)[pca.x]
        PC2 <- names(tmp.pcaloadings)[pca.y]
        PC3 <- names(tmp.pcaloadings)[pca.z]

        my.scatter <- function(datafr,xa,ya,topx,topy){

                #scatterplot dimensions pairwise (pairs passed to function as x,y)
                x.axis <- xa ; y.axis <- ya ;topxy <- c(names(topx),names(topy))
                # This is the data.frame of the TopN loadings to be marked
                mark.frame <- tmp.pcaloadings[which(row.names(tmp.pcaloadings) %in% topxy),]

                xmin = min(datafr[,x.axis]);xmax = max(datafr[,x.axis])
                ymin = min(datafr[,y.axis]);ymax = max(datafr[,y.axis])

                #make the scatterplot
                ggplot(data = datafr,aes_string(x = x.axis,y = y.axis))+
                        geom_hline(yintercept = 0, linetype = "dashed",size =0.7)+
                        geom_vline(xintercept = 0, linetype = "dashed",size = 0.7)+
                        geom_point(color = "navy",alpha = 0.15,size=0.1, show.legend = FALSE)+
                        geom_label_repel(data = mark.frame, size = 3, label.r = unit(0.45,"lines"),
                                         color = "deeppink", bg = "plum1",
                                         segment.size = 0.1, box.padding = unit(1,"lines"),
                                         aes(label = rownames(mark.frame)))+
                        geom_point(data = mark.frame,
                                   color = "deeppink",
                                   size = 2,show.legend = FALSE)+
                        xlim(c(xmin-0.01,xmax+0.01))+ylim(c(ymin-0.01,ymax+0.01))+
                        theme_bw()+
                        theme(panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank()
                        )
                #mark the TopN ids

        }

        # call my.scatter 3 times to make 3 ggplots
        g1 <- my.scatter(tmp.pcaloadings,PC1,PC2,x,y)
        g2 <- my.scatter(tmp.pcaloadings,PC2,PC3,y,z)
        g3 <- my.scatter(tmp.pcaloadings,PC1,PC3,x,z)
        #combine them in a row by using multiplot() function
        multiplot(g1,g2,g3,cols = 3)
}


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
    ##tab=tab[, names(grp.vec)]

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

    tab[filt.idx, ] <- NA

    tab <- data.frame(ids, tab)
    colnames(tab)[1] <- id.col

    ##View(tab)
    values.filt <- lapply(groups, function(x) filt.idx)

    return(list(table=tab, values.filtered=values.filt))
}



########################################################################
## 20160224
##                   reproducibility filter
##
## n=2: Bland-Altman
## n>2: lmm-model written by Mani DR
##
## - replaces not reprodicibly measuered values in 'tab' with 'NA'
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
dynamicHeightHM <- function(n, unit=c('px', 'in')){

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
makeProfileplot <- function(tab, id.col, grp, grp.col, grp.col.leg, legend=T, cex.lab=1.5, mar=c(5,5,3,1), ... ){

    cat('\n-- makeProfileplot --\n')

    ## table
    tab <- tab[, setdiff(colnames(tab), id.col)]

    xlim=max(abs(tab), na.rm=T)

    ## caclulate densities
    dens <- apply(tab, 2, density, na.rm=T)

    ## ylim
    ylim <- max(unlist(lapply(dens, function(x) max(x$y))))

    ##########################################
    ## plot
    par(mar=mar)
    for(i in 1:ncol(tab)){
        if(i == 1)
            plot(dens[[i]], xlab='expression', xlim=c(-xlim, xlim), ylim=c(0, ylim), col=my.col2rgb(grp.col[i], alpha=100), lwd=3, cex.axis=2, cex.lab=2, cex.main=1.5, ...)
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

    cat('\n-- makeProfileplot exit --\n')
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
        ##
        ppi.col <- rep('', length(IDs))
        #names(ppi.col) <- IDs
        names(ppi.col) <- names(IDs)
        
        for( d in 1:length(ppi.int.detect.l))
          ppi.col[ which( IDs %in% ppi.int.detect.l[[d]] )] <- ppi.db.col[ names(ppi.int.detect.l)[d]  ]
        
          #ppi.col[ which(names(ppi.col) %in% ppi.int.detect.l[[d]] )] <- ppi.db.col[ names(ppi.int.detect.l)[d]  ]
        
          #ppi.col[ ppi.int.detect.l[[d]] ] <- ppi.db.col[ names(ppi.int.detect.l)[d]  ]


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
