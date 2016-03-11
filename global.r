################################################################################################################
## Filename: global.r
## Created: October 09, 2015
## Author(s): Karsten Krug, Mani DR
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
################################################################################################################

#################################################################
## global parameters
#################################################################
## version number
VER=0.3
## maximal filesize for upload
MAXSIZEMB <<- 400
## list of strings indicating missing data
NASTRINGS <<- c("NA", "<NA>", "#NUM!", "#DIV/0!", "#NA", "#NAME?")
## speparator tested in the uplosded file
SEPARATOR <<- c('\t', ',', ';')
## Colors used throughout the app to color the defined groups
GRPCOLORS <<- c(RColorBrewer::brewer.pal(9, "Set1"), RColorBrewer::brewer.pal(8, "Dark2"), RColorBrewer::brewer.pal(8, "Set2"))
## number of characters to display in plots/tables for column names
STRLENGTH <<- 20
## temp directory to write the Excel file
TMPDIR <<- ifelse(Sys.info()['sysname']=='Windows', "./", "/tmp/")
## app name
APPNAME <<- sub('.*/','',getwd())

#################################################################
## load required packages
#################################################################
library(shiny)
## heatmap
library(pheatmap)
## moderated t-test
library(limma)
## colors
library (RColorBrewer)
## multiscatter
library(hexbin)
library(Hmisc)
library(grid)
## pca
library(scatterplot3d)
library(plotly)
## export
library(WriteXLS)
## reproducibility filter
library(reshape)
library(nlme)
library(BlandAltmanLeh)
## normalization Quantile
library(preprocessCore)
## normalization 2-component
library (mice)
library (mixtools)
library (mclust)


#################################################################################################
##                     multiscatterplot using hexagonal binning
## - mat    numerical matrix of expression values, rows are features, columns are samples
##
## changelog: 2015116 implementation
#################################################################################################
my.multiscatter <- function(mat, hexbin=30, hexcut=5, cor=c('pearson', 'spearman', 'kendall'), repro.filt=NULL, grp, grp.col.legend, define.max=F, max.val=3, min.val=-3){

    ## sort table according to group vector
    ##grp <- sort(grp)
    ##mat <- mat[, names(grp)]

    ## cor method
    corm = match.arg(cor)
    ## correlation
    cm = cor(mat, use='pairwise.complete', method=corm)

    ## number of samples to compare
    N = ncol(mat)

    ## define limits
    if(define.max){
        lim=c(min.val, max.val)

    } else{
        ## determine x and y limits
        lim=max( abs( mat ), na.rm=T )
        lim=c(-lim, lim)
    }

    ##cat(grp.col.legend, "\n", names(grp.col.legend), '\n\n')

    ###########################################################################
    ## help function to set up the viewports
    ## original code from:  http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
    multiplot <- function(plots, cols=1) {
        ## Make a list from the ... arguments and plotlist
        ##plots <- c(list(...))
        ## number of plots
        numPlots = length(plots)
        ## layout matrix
        layout <- matrix(seq(1, cols * ceiling(numPlots/cols)), ncol = cols, nrow = ceiling(numPlots/cols))
        ## Set up the page
        grid.newpage()
        ## grid layout
        la <-  grid.layout(nrow(layout), ncol(layout))
        pushViewport(viewport(layout = la))
        ## Make each plot, in the correct location
        for (i in numPlots:1) {
            ## Get the i,j matrix positions of the regions that contain this subplot
            matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
            vp = viewport(layout.pos.row = matchidx$row, layout.pos.col = matchidx$col)

            ## textplot: correlation coefficient
            if(matchidx$row < matchidx$col){
                numb = plots[[i]]
                col='black'
                size = min(max(abs(90*as.numeric(numb)), 25), 50)
                grid.rect(width=unit(.85, 'npc'), height=unit(.85, 'npc'), vp=vp, gp=gpar(fill='grey95', col='transparent'))
                grid.text(numb, vp=vp, gp=gpar(fontsize=size, col=col))
            } else {
                print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                                layout.pos.col = matchidx$col))
            }
        }
    } ## end function 'multiplot'
    #########################################################################
    ##
    ##                     do the actual plotting
    ##
    #########################################################################

    ## list to store the plots
    plotList=vector('list', N*N)
    count=1
    for(i in 1:N)
        for(j in 1:N){

            ## extract pairwise data
            dat <- data.frame(x=mat[,i], y=mat[,j])
            rownames(dat) <- rownames(mat)

            ## extract groups
            current.group <- unique(grp[names(grp)[c(i,j)]])

            ##cat(current.group, '\n')
            ##str(add.points)

            ###########################
            ## lower triangle
            if(i < j){

                ## hexbin
                hex <- hexbin(dat$x, dat$y, hexbin)
                gghex <- data.frame(hcell2xy(hex), c = cut2(hex@count, g = hexcut))
                p <- ggplot(gghex) + geom_hex(aes(x = x, y = y, fill = c) ,stat = "identity") + guides(fill=FALSE) + theme( plot.margin=unit(rep(0, 4), 'cm')) + xlab('') + ylab('') + xlim(lim[1], lim[2]) + ylim(lim[1], lim[2])

                ##if(length(current.group) == 1)
                ##    p <- p + scale_fill_manual( values=paste(rep( grp.col.legend[current.group], hexcut) ))
                ##else

                p <- p + scale_fill_manual( values=paste('grey', ceiling(seq(70, 20, length.out=hexcut)), sep=''))

                ## add filtered values
                if(!is.null(repro.filt) & length(current.group) == 1){
                    not.valid.idx <- repro.filt[[current.group]]
                    dat.repro <- dat[not.valid.idx, ]
                    ##cat(not.valid.idx)
                    ##cat(dim(dat.repro))
                    p = p + geom_point( aes(x=x, y=y ), data=dat.repro, colour=my.col2rgb('red', 100), size=.5)
                }
            }
            ###########################
            ## diagonal
            if(i == j){
                p = ggplot(dat, aes(x=x)) + geom_histogram(fill=grp.col.legend[current.group], colour=grp.col.legend[current.group], binwidth=sum(abs(range(dat$x, na.rm=T)))/50) + ggtitle(colnames(mat)[i]) + theme(plot.title=element_text(size=9)) + theme( panel.background = element_blank(), plot.margin=unit(rep(0, 4), 'cm')) + xlab(paste('N',sum(!is.na(dat$x)), sep='=')) + ylab('') + xlim(lim[1], lim[2]) ##+ annotate('text', label=sum(!is.na(dat$x)), x=unit(0, 'npc'), y=unit(0, 'npc'))


            }
            ###########################
            ## upper triangle
            if(i > j){
                cortmp = cm[i,j]

                p=paste(round(cortmp,2))

            }

            plotList[[count]] <- p
            count=count+1
        }

    multiplot( plotList, cols=N)
}


##########################################################################################################
##                     translate a color name into rgb space
##
## changelog:  20100929 implementation
##########################################################################################################
my.col2rgb <- function(color, alpha=80, maxColorValue=255){

    out <- vector( "character", length(color) )

    for(col in 1:length(color)){

        col.rgb <- col2rgb(color[col])

        out[col] <- rgb(col.rgb[1], col.rgb[2], col.rgb[3], alpha=alpha, maxColorValue=maxColorValue)

    }
    return(out)
}

#####################################################################################
##
##                        two sample moderated t-test
##
## - code written by mani dr
## - code modified by karsten krug
##   20151211 'label'
#####################################################################################
modT.test.2class <- function (d, output.prefix, groups, id.col=NULL, data.col=NULL,
                              group.na.rm=FALSE, nastrings=c("NA", "<NA>", "#NUM!", "#DIV/0!", "#NA", "#NAME?"), label=NULL) {

    ## store group names
    groups.org <- groups
    groups <- as.numeric(as.factor(groups))

    id <- d[ , id.col]

    ## extract data columns
    if (is.null (data.col)) data <- d [, setdiff (colnames (d), id.col)]
    else data <- d [, make.names (data.col)]


    ## moderated t test for 2 classes
    design.mat <- cbind (ref=1, comparison=groups)
    mod.t.result <- moderated.t (data, design.mat)

    ## 20151211 kk
    mod.t.result <- data.frame( mod.t.result, Log.P.Value=-10*log(mod.t.result$P.Value,10))

    ## add label
    if(!is.null(label))
        colnames(mod.t.result) <- paste(colnames(mod.t.result), label, sep='.')

    mod.t <- data.frame ( cbind (data.frame (id), data, mod.t.result) )
    ##mod.t <- data.frame ( cbind (data.frame (id),  mod.t.result, data) )
    rownames(mod.t) <- id ##make.unique(as.character(mod.t[, 1]), sep='_')
    colnames(mod.t)[1] <- 'id'

    ##write.csv (final.results, paste (output.prefix, ".csv", sep=''), row.names=FALSE)

    ## write out / return results
    final.results <- mod.t

    ##invisible (final.results)
    return( list(input=d, output=final.results, groups=groups.org) )
}


######################################################################################################
##                               One-sample moderated t-test
##
## run moderated t-test, and plot results
## mainly for iTRAQ, but can be used of other data
##
## code written by mani dr
## code modified by karsten krug
## 20151210 'label'
######################################################################################################
modT.test <- function (d, output.prefix, id.col=NULL, data.col=NULL, fix.id=FALSE,
                       p.value.alpha=0.05, use.adj.pvalue=TRUE, apply.log=FALSE,
                       na.rm=FALSE, nastrings=c("NA", "<NA>", "#NUM!", "#DIV/0!", "#NA", "#NAME?"),
                       plot=TRUE, pairs.plot.2rep=FALSE, limits=NULL, xlab="", ylab="", label='', ...) {
  #
  # data.file should contain one peptide in each row.
  # The columns contain the normalized log-ratio from each replicate
  # (technical or biological). The ratio is based on classes of interest
  # that need to be distinguished: i.e., ratio = intensity_A / intensity_B;
  # this test calculates the p-value for determining if peptide p is
  # differentially regulated between classes A and B (i.e., if the log
  # ratio is different from 0).
  # While the standard scatter plot routine can only handle 2 replicates,
  # a pairs plot is created when there are more than 2 replicates.
  # The moderated t-test can be applied to any number of replicates.
  # An id column can be optionally included in the data.file to track
  # peptides (row numbers are used as id if a column is not specified).
  #
  # graphics can be controlled using ...
  #  when using scatterhist (for 2 replicates), this can be arguments to points
  #  when > 2 replicates are present, ... can include arguments to points in
  #   addition to: plot.col, subset.col, hist.col, hist.breaks,
  #                prefix (for correlation), cex.cor

  id <- d[,id.col]
  ##if ( any (duplicated (id)) ) stop ('IDs are not unique. Use fix.id=TRUE option')

  # extract data columns
  if (is.null (data.col)) data <- d [, setdiff (colnames (d), id.col)]
  else data <- d [, make.names (data.col)]


  # log transform is required
  if (apply.log) data <- log2 (data)

  # moderated t test
  mod.t.result <- moderated.t (data)
  if (use.adj.pvalue) mod.sig <- mod.t.result [,'adj.P.Val'] <= p.value.alpha
  else  mod.sig <- mod.t.result [,'P.Value'] <= p.value.alpha
  change <- apply (data, 1,
                   function (x) {
                     x <- x [is.finite (x)]
                     ret.value <- '?'
                     if ( all (x < 0) ) ret.value <- 'down'
                     else if ( all (x > 0)) ret.value <- 'up'
                     return (ret.value)
                   })
    ## 20151210 kk
    mod.t.result <- data.frame( mod.t.result, change=change, significant=mod.sig, Log.P.Value=-10*log(mod.t.result$P.Value,10))

    ## add label
    if(!is.null(label))
        colnames(mod.t.result) <- paste(colnames(mod.t.result), label, sep='.')

    ##mod.t <- data.frame ( cbind (data.frame (id), data, mod.t.result, change=change, significant=mod.sig) )
    mod.t <- data.frame ( cbind (data.frame (id), data, mod.t.result) )
    colnames (mod.t)[1] <- id.col   # retain id.col (if provided)
    rownames(mod.t) <- make.unique( as.character(mod.t[,1]), sep='_' )
    colnames(mod.t)[1] <- 'id'

    final.results <- mod.t
    return( list(input=d, output=final.results) )
}

#############################################################################################
##
##              different normalization methods for expression data
##
## 20160235
#############################################################################################
normalize.data <- function(data, id.col, method=c('Median', 'Quantile', 'Median-MAD', '2-component')){
    ##cat('norm start\n')

    method = match.arg(method)

    ids = data[, id.col]
    data = data[ , -grep(paste('^', id.col, '$', sep=''), colnames(data))]

    data <- data.matrix(data)

    ## quantile
    if(method == 'Quantile'){
        require("preprocessCore")
        data.norm <- normalize.quantiles(data)
        rownames(data.norm) <- rownames(data)
        colnames(data.norm) <- paste( colnames(data))

        ## shift median to zero
        data.norm <- apply(data.norm, 2, function(x) x - median(x, na.rm=T))
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
        ## check if successful
        for(i in 1:length(data.norm.list)){
            if(length(data.norm.list[[i]]) == 1){
                if(data.norm.list[[i]] == 'No_success')
                    return('No_success')
            }
        }
        data.norm = matrix( unlist(lapply(data.norm.list, function(x)x$norm.sample)), ncol=length(data.norm.list), dimnames=list(rownames(data), names(data.norm.list)) )
    }
    ## add id column
    data.norm <- data.frame(ids, data.norm)
    colnames(data.norm)[1] <- id.col
    return(data.norm)
}

###########################################################################
##
##           Moderated t-test for significance testing
##
## code written by Mani DR
## code modified by karsten krug
##    - single function for one/two sample test
##
##########################################################################
moderated.t <- function (data, design=NULL) {
    ## data is a table with rows representing peptides/proteins/genes
    ## and columns representing replicates

    data.matrix <- data.frame (data)
    ## the design matrix is expected to be:
    ##    ref    comparison
    ##     1         0
    ##    1         0
    ##         ...
    ##     1         1
    ##  where comparison has 0's and 1's corresponding
    ##  to the column position of the 2 classes in data
    ##  (see limma user manual section 13)


    #############################################
    ## two sample test
    if(!is.null(design)){
        m <- lmFit (data.matrix, design)
        m <- eBayes (m)
        sig <- topTable (m, coef=colnames (design)[2], number=nrow(data), sort.by='none')
    } else {
    #############################################
    ## one sample test
        m <- lmFit (data.matrix, method='robust')
        m <- eBayes (m)
        sig <- topTable (m, number=nrow(data), sort.by='none')
    }
  return (sig)
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
##############################################################################
my.prcomp <- function(x, col=NULL, cor=T, plot=T, rgl=F, scale=T, pch=20, cex.points=3, rgl.point.size=30, main="PCA", leg.vec=NULL, leg.col=NULL, ...){

    cex.font = 1.8

    # color
    if( is.null(col) ) col="black"

    # perform pca
    pca <- prcomp(x, scale=scale)

    # calculate variance
    comp.var <- eigen(cov(pca$x))$values

    ## extract the principle components
    pc1=pca$x[,1]
    pc2=pca$x[,2]
    pc3=pca$x[,3]

    ##############
    # rgl plot
    ##############
    if(rgl){
        require(rgl)
        plot3d(pc1, pc2, pc3, xlab=paste("PC 1 (", round(100*comp.var[1]/sum(comp.var),1),"%)", sep=""), ylab=paste("PC 2 (",round(100*comp.var[2]/sum(comp.var),1),"%)", sep=""), zlab=paste("PC 3 (", round(100*comp.var[3]/sum(comp.var),1),"%)", sep=""), type="s", col=col, expand=1.2, size=rgl.point.size)
    }

    ########################################
    # scatterplot 2D/3D
    ########################################
    if( plot){
         require(scatterplot3d)

         par(mfrow=c(1,3), mar=c(7,7,3,1))

         ## PC 1-2
         plot(pc1, pc2, xlab=paste("PC 1 (", round(100*comp.var[1]/sum(comp.var),1),"%)", sep=""), ylab=paste("PC 2 (",round(100*comp.var[2]/sum(comp.var),1),"%)", sep=""), pch=pch, main=main, col=col, sub=paste("Cumulative variance = ", round(100*sum(comp.var[1:2]/sum(comp.var)),1),"%", sep=""), cex=cex.points, ylim=c( min(pc2),  max(pc2)+.15*max(pc2)), cex.axis=cex.font, cex.lab=cex.font, cex.sub=cex.font )

         ##if(!is.null(leg.vec) & !is.null(leg.col))
         ##    legend('top', legend=leg.vec, col=leg.col, pch=pch, pt.cex=max(1, cex.points-1.5), ncol=length(leg.vec))
             ##legend('top', legend=leg.vec, col=leg.col, pch=pch, pt.cex=cex.points, ncol=length(leg.vec))


         ## PC 1-3
         scatterplot3d( pca$x[,1], pca$x[,3], pca$x[,2], xlab=paste("PC 1 (", round(100*comp.var[1]/sum(comp.var),1),"%)", sep=""), zlab=paste("PC 2 (",round(100*comp.var[2]/sum(comp.var),1),"%)", sep=""), ylab=paste("PC 3 (", round(100*comp.var[3]/sum(comp.var),1),"%)", sep=""), color=col,  cex.symbols=cex.points, pch=pch, main=main, sub=paste("Cumulative variance = ", round(100*sum(comp.var[1:3]/sum(comp.var)),1),"%", sep=""), type="h" )


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
## color ramp
##
## ToDo: opacity!
####################################################
myColorRamp <- function(colors, values, range=NULL) {

    if(is.null(range))
        v <- (values - min(values))/diff(range(values))
    else
        v <- (values - min(values, na.rm=T))/diff( range )

    x <- colorRamp(colors)(v)
    rgb(x[,1], x[,2], x[,3], maxColorValue = 255)
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
chopString <- function(string, nChar=10)
{

    string.trim <- strtrim(string, nChar)
    string.trim[ which(nchar(string) > nChar) ] <-  paste(string.trim[which(nchar(string) > nChar) ], '...')

    return(string.trim)

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

    ## extract groups
    groups = unique(grp.vec)

    ## list to store index of filtered values per group
    values.filt <- vector('list', length(groups))
    names(values.filt) <- groups

    ## add rownames to tab
    ##rownames(tab) <- tab[, id.col]

    tab.repro.filter <- tab
    ##View(tab.repro.filter)

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
            values.filt[[gg]] <- not.repro.idx
        }
        ########################################
        ## if there are two replicates use
        ## Blandt-Altmann filter
        ## R-package 'BlandAltmanLeh'
        if( length(gg.idx) == 2 ){

            ## Bland-Altman
            ba <-  bland.altman.stats(as.numeric( as.character( tab[, gg.idx[1] ]) ), as.numeric( as.character( tab[,  gg.idx[2] ] )), two=3.290527 )
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

############################################################################################
##
##              Generalized reprodicibility filter for > 2 replicates
##
## written by Mani DR
############################################################################################
reproducibility.filter <- function (data, id.col='id', alpha=0.05) {
  ##
  ## Reproducibility Filter using lme4
  ## Theory: MethComp book (pp 58-61). Comparing Clinical Measurement Methods by Bendix Carstensen
  ## Implementation: MethComp book pg 142, but don't include item (=id) in the fixed effects
  ## -- this is unnecessary for the application and makes computations very time consuming;
  ## all we really need to assess reproducibility are the variances
  ##
  ## NB: using library (nlme) and lmer is much more convoluted since incorporating the
  #      var-cov matrix stratified by replicate (=method) is not easy
  #      (see https://stat.ethz.ch/pipermail/r-sig-mixed-models/2007q3/000248.html)
  #      something like:
  #        model <- lmer (y ~ rep + (rep1|id)+(rep2|id)+..., data=data.long)
  #      is possible, but the above does not work for 2 replicates (gives same results for >2 reps)
  #

  d <- data [, setdiff (colnames (data), id.col)]     # data part of input
  data.long <- melt (data.frame (data), id=id.col)    # convert to long format
  colnames (data.long) <- c ('id', 'rep', 'y')
  # keep column order in data so that (i,j) below correctly corresponds to columns
  data.long [,'rep'] <- factor (data.long[,'rep'], levels=colnames (d))

  # exclude missing data points (only missing measurements are removed instead of entire rows)
  data.long <- data.long [ !is.na (data.long[,'y']), ]

  # Model: y_mi = a_m + c_mi + e_mi,  c_mi ~ N(0,tau_m^2), e_mi ~ N(0, sigma_m^2)
  # where m = method and i = item (=id)
  # [Eq 5.2, pg 58, MethComp book]. Also see interpretation of effect on pg 59-61
  model <- lme (y ~ rep,
                random=list (id=pdIdent(~rep)),
                weights=varIdent(form=~1|rep),
                data=data.long)
  n <- nlevels (data.long[,'rep'])
  p <- length (unique (data.long[,'id']))
  df <- p - 1    # approx df for confidence interval (p=# of independent items)

  rep.all <- rep (TRUE, nrow (d))  # vector summarizing reproducibility of each input id
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      # variance of method_i - method_j: pg 58
      # var (y_i0 - y_j0) = tau_i^2 + tau_i^2 + sigma_i^2 + sigma_j^2
      # where tau is the sd of the between-item variation for a method
      # and sigma is the sd of the within-item variation for the method
      tau <- as.numeric (unlist (VarCorr(model)[c(i,j),'StdDev']))  # returns tau_i and tau_j
      # VarCorr(model)[n+1,'StdDev'] is sigma_1 (ie sigma for method 1)
      # for methods 2-n, coef (model$modelStruct$varStruct, uncons=F, allcoef=T) has the
      #  multiplying factor to obtain sigma_m using sigma_1
      sigma <- as.numeric (unlist (VarCorr(model)[n+1,'StdDev']))   #
      sigma1 <- sigma * ifelse (i==1, 1, coef (model$modelStruct$varStruct, uncons=F, allcoef=T)[i-1])
      sigma2 <- sigma * coef (model$modelStruct$varStruct, uncons=F, allcoef=T)[j-1]
      total.sd <- sqrt (tau[1]^2 + tau[2]^2 + sigma1^2 + sigma2^2)

      # bias of method_i - method_j: alpha_i - alpha_j
      alpha1 <- ifelse (i==1, 0, fixef(model)[i])
      alpha2 <- fixef (model)[j]
      bias <- alpha1 - alpha2

      # limits of agreement (assuming approx df = p-1)
      t.crit <- qt ( (1-alpha/2), df ) * sqrt ( (p+1) / p )
      ci.low <- bias - t.crit * total.sd
      ci.high <- bias + t.crit * total.sd

      # record reproducibility for method_i - method_j
      rep.ij <- (d[,i] - d[,j]) >= ci.low & (d[,i] - d[,j]) <= ci.high
      # if data is missing, assume that data is reproducible
      rep.ij [ is.na (rep.ij) ] <- TRUE
      rep.all <- rep.all & rep.ij

      # print bias and LoA for sanity check
      ##cat ('rep', i, ' - rep', j, ': bias=', bias, ' ci=(', ci.low, ',', ci.high, ')\n', sep='')
    }
  }
    ## karsten krug 20160301
    ## return rownames of data matrix
    return(rep.all)
    ##return( rownames(data)[ which(!rep.all) ] )
}


##################################################################
## function to dynamically determine the height (in px) of the heatmap
## depending on the number of genes
dynamicHeightHM <- function(n){
    if( n < 50)
        height=500
    if( n >= 50 & n <= 100)
        height=800
    if(n >=100)
        height=800+n
    return(height)
}

###################################################################
##
##       generate the boxplots under the 'QC' tab
##
###################################################################
makeBoxplot <- function(tab, id.col, grp, grp.col, grp.col.leg, legend=T, cex.lab=1.5, mar=c(4,15,2,6)){

    ## table
    tab <- tab[, setdiff(colnames(tab), id.col)]

    ##########################################
    ## order after groups
    ord.idx <- order(grp)
    grp <- grp[ord.idx]
    tab <- tab[, ord.idx]
    grp.col <- grp.col[ ord.idx]


    at.vec=1:ncol(tab)
    ##########################################
    ## plot
    par(mar=mar)
    boxplot(tab, pch=20, col='white', outline=T, horizontal=T, las=2, xlab=expression(log[2](ratio)), border=grp.col, at=at.vec, axes=F, main='', cex=2, xlim=c(0, ifelse(legend, ncol(tab)+2, ncol(tab)) ))
    ##legend('top', legend=names(grp.col.leg), ncol=2, bty='n', border = names(grp.col.leg), fill='white', cex=1.5)
    if(legend)
        legend('top', legend=names(grp.col.leg), ncol=length(grp.col.leg), bty='n', border = grp.col.leg, fill=grp.col.leg, cex=cex.lab)
    ##legend('top', legend=c(input$label.g1, input$label.g2), ncol=2, bty='n', border = c('grey10', 'darkblue'), fill='white', cex=1.5, lwd=3)
    mtext( paste('N=',unlist(apply(tab,2, function(x)sum(!is.na(x)))), sep=''), at=at.vec, side=4, las=2, adj=0, cex.lab=cex.lab)
    axis(1)
    axis(2, at=at.vec, labels=chopString(colnames(tab), STRLENGTH), las=2, cex=cex.lab)


}
###########################################################################################
##
##                    two-component mixture model normalization
## written by Mani DR
##
##########################################################################################
two.comp.normalize <- function (sample, type) {
  #   1. For all sample types, fit a 2-component gaussian mixture model using normalmixEM.
  #   2. For the bimodal samples, find the major mode M1 by kernel density estimation
  #     2a. Fit the model with one component mean constrained to equal M1
  #     2b. Normalize (standardize) samples using mean (M1) and resulting std. dev.
  #   3. For unimodal samples, find the mode M using kernel density estimation
  #     3a. Fit the model with mean for both components constrained to be equal to M
  #     3b. Normalize (standardize) samples using mean M and smaller std. dev. from model fit

  # WARNING:
  # This code has a lot of hacks to fix the flakiness of normalmixEM, and the idiosyncracies
  # of the actual data. Carefully re-examine code for new or altered input data

  data <- sample [ !is.na (sample) ]
  dens <- density (data, kernel='gaussian', bw='SJ')     # gaussian kernel with S-J bandwidth
                                                         # (see Venalbles & Ripley, 2002, pg, 129)
  # find major (highest) mode > -3 (to avoid problems with lower mode having higher density than higher mode)
  x.range <- dens$x > -3
  dens.x <- dens$x [x.range];  dens.y <- dens$y [x.range]
  mode <- dens.x[which.max(dens.y)]
  if (type=='bimodal') mean.constr <- c (NA, mode) else mean.constr <- c (mode, mode)
  model <- normalmixEM (data, k=2, mean.constr=mean.constr)
  model.rep <- normalmixEM (data, k=2, mean.constr=mean.constr)
  model.alt <- Mclust (data, G=2, modelNames="V")
  alt.mu <- model.alt$parameters$mean
  alt.sd <- sqrt (model.alt$parameters$variance$sigmasq)
  # find reproducible model fit that is close to Mclust fit
  # if not, re-fit model -- without this condition
  # normalmixEM produces one-off model fits
  n.try <- 1
  if (type=='unimodal') model.mode <- which(model$mu==mode)[which.min (model$sigma)]
  else model.mode <- which(model$mu==mode)
  model.other <- model.mode %% 2 + 1
  alt.mode <- which.min(abs(model.alt$par$mean-mode)); alt.other <- alt.mode %% 2 + 1
  while ( abs (model$mu[model.mode] - alt.mu[alt.mode]) > 3e-1 || abs (model$sigma[model.mode]-alt.sd[alt.mode]) > 3e-1 ||
          model$sigma[model.mode] < 0.1 ||
          (type=='bimodal' && (abs (model$mu[model.other] - alt.mu[alt.other]) > 1e1)) ||
          !all (c (model$mu, model$sigma) - c (model.rep$mu, model.rep$sigma) < 1e-3) ) {
    # if major mode (and SD of mode) is not within 0.3, or if the other mean (for bimodals only)
    # is not within 1 of the Mclust result, try again
    model <- normalmixEM (data, k=2, mean.constr=mean.constr)
    model.rep <- normalmixEM (data, k=2, mean.constr=mean.constr)

      if (n.try > 50){
          return("No_success")
          ##stop (paste ("Can't fit mixture model ... giving up\n"))
      }
    n.try <- n.try + 1
  }


  if (type=='bimodal') {
    # sometimes (esp. in phosphoproteome) the minor (lower) mode can be larger than the major (higher) mode
    # this situation is not possible in the unimodal samples
    corrected.mode <- model$mu [which.max(model$mu)]
    if (corrected.mode != mode) {
      cat ('  Lower mode larger than higher mode\n')
      mode <- corrected.mode
    }
  }
  norm.mean <- mode
  norm.sd <- ifelse (type=='bimodal', model$sigma[which(model$mu==mode)], min (model$sigma))

  # normalize by standardizing
  data <- data - norm.mean
  data <- data / norm.sd

  # return normalized data reorganized to original order
  sample [ !is.na (sample) ] <- data
  return ( list (norm.sample=sample, norm.mean=norm.mean, norm.sd=norm.sd, fit=unlist (c(model$mu, model$sigma))) )
}
