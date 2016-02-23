
## modified by Karsten Krug, 10/09/2015

# $Id: two-class.r 60 2015-08-24 19:04:15Z manidr $
##print ('$Id: two-class.r 60 2015-08-24 19:04:15Z manidr $')

##library(limma)
##library (gplots)
##library (RColorBrewer)

##library(hexbin)
##library(ggplot2)
##library(Hmisc)
##library(grid)

##library(scatterplot3d)

#################################################################################################
##            multiscatterplot using hexagonal binning
## - mat    numerical matrix of expression values, rows are features, columns are samples
##
## changelog: 2015116 implementation
#################################################################################################
my.multiscatter <- function(mat, hexbin=30, hexcut=5, cor=c('pearson', 'spearman', 'kendall')){


    ###########################################################################
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
                ##col=
                size = max(abs(90*as.numeric(numb)), 25)
                grid.rect(width=unit(.85, 'npc'), height=unit(.85, 'npc'), vp=vp, gp=gpar(fill='grey95', col='transparent'))
                grid.text(numb, vp=vp, gp=gpar(fontsize=size))
            } else {
                print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                                layout.pos.col = matchidx$col))
            }
        }
    } ## end function 'multiplot'
    #########################################################################

    ## cor method
    corm = match.arg(cor)
    ## correlation
    cm = cor(mat, use='pairwise.complete', method=corm)

    ## number of samples to compare
    N = ncol(mat)

    ## list to store the plots
    plotList=vector('list', N*N)
    count=1
    for(i in 1:N)
        for(j in 1:N){

            dat <- data.frame(x=mat[,i], y=mat[,j])

            ###########################
            ## lower triangle
            if(i < j){

                ## hexbin
                hex <- hexbin(dat$x, dat$y, hexbin)
                gghex <- data.frame(hcell2xy(hex), c = cut2(hex@count, g = hexcut))
                p <- ggplot(gghex) + geom_hex(aes(x = x, y = y, fill = c) ,stat = "identity") + guides(fill=FALSE) + theme( plot.margin=unit(rep(0, 4), 'cm')) + scale_fill_manual( values=paste('grey', ceiling(seq(70, 20, length.out=hexcut)), sep='')) + xlab('') + ylab('')
            }
            ###########################
            ## diagonal
            if(i == j){
                p = ggplot(dat, aes(x=x)) + geom_histogram(fill='grey70', colour='black', binwidth=sum(abs(range(dat$x, na.rm=T)))/50) + ggtitle(colnames(mat)[i]) + theme(plot.title=element_text(size=20)) + theme( panel.background = element_blank(), plot.margin=unit(rep(0, 4), 'cm')) + xlab(paste('N',sum(!is.na(dat$x)), sep='=')) + ylab('') ##+ annotate('text', label=sum(!is.na(dat$x)), x=unit(0, 'npc'), y=unit(0, 'npc'))

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
## translate a color name into rgb space
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

    ## read data file
    ## d <- read.csv (data.file, na.strings=nastrings)


    ## extract id column
    ## if id column is not specified, use row numbers
    row.num <- 1:nrow (d)
    if (is.null (id.col)) {
        d <- cbind (id=row.num, d)
        id.col <- c ('id')
    }

    id.col <- make.names (id.col)
    if (group.na.rm) {
        ## ... and remove rows with all missing values for a group
        d.gr <- d [, setdiff (colnames (d), id.col)]
        nas.gr <- apply (d.gr, 1, function (x) { all (is.na (x[groups==0])) || all (is.na (x[groups==1])) })
        no.na <- !nas.gr
        d <- d [no.na, ]
    }
    id <- d[,id.col]

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
    rownames(mod.t) <- make.unique(as.character(mod.t[,1]), sep='_')
    colnames(mod.t)[1] <- 'id'

    ##write.csv (final.results, paste (output.prefix, ".csv", sep=''), row.names=FALSE)

    ## write out / return results
    final.results <- mod.t

    ##invisible (final.results)
    return( list(input=d, output=final.results, groups=groups.org) )
}


######################################################################################################
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

  # read data file
  ##d <- read.csv (data.file, na.strings=nastrings)
  if (na.rm) {
    # ... and remove rows with any missing values
    no.na <- apply (d, 1, function (x) { all (!is.na (x))} )
    d <- d [no.na,]
  }


  # extract id column
  # if id column is not specified, use row numbers
  row.num <- 1:nrow (d)
  if (is.null (id.col)) {
    d <- cbind (id=row.num, d)
    id.col <- c ('id')
  }

  id.col <- make.names (id.col)
  if (fix.id) d[,id.col] <- make.unique (d[,id.col])
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

  # plot and/or determine prediction CI's, etc.
  # the output/behavior of scatterhist.ci can be controlled using ... for specifying args
  results <- NULL
  if (plot) {
    if (ncol (data)==2 && !pairs.plot.2rep) {
      if (!exists ('libdir')) source ('scatterhist.r')    # source only when not in GenePattern
      # when 2 replicates are present, create a scatterplot with marginal histograms
      pdf (paste (output.prefix, ".pdf", sep=''), width=6, height=6, pointsize=12)
      results <- scatterhist.ci (data[,1], data[,2], id=id, ci.level=ci.level,
                                 special.subsets=list (mod.t[,'significant']),
                                 subset.colors=c('red'), limits=limits, xlab=xlab, ylab=ylab, ...)
      if (length (setdiff (colnames (results), c ('x','y','id'))) != 0) {
        results <- results [, setdiff (colnames (results), c ('x','y'))]  # remove x and y -- data already has this
        colnames (results)[1] <- id.col   # retain id.col (if provided)
      } else results <- NULL
      dev.off ()
    } else if (ncol (data) >= 2) {
      # create a pairs plot when more than 2 replicates are present
      if (!exists ('libdir')) source ('pairs-plot.r')     # source only when not in GenePattern
      keep <- !is.na (mod.sig)
      plot.data <- data [keep,]
      significant <<- mod.sig [keep]

      if (is.null (limits)) {
        # use full range of x, y if limits not specified
        r <- apply (plot.data, 2, range, na.rm=TRUE)
        limits <- c ( min(r), max(r) )
      }

      size.in <- ncol (data) + 1
      pdf (paste (output.prefix, ".pdf", sep=''), height=size.in, width=size.in, pointsize=8)
      pairs (plot.data, diag.panel=panel.hist, lower.panel=panel.scatterplot, upper.panel=panel.cor,
             subset=significant, xlim=limits, ylim=limits, ...)
      dev.off()
    } else warning ("No plots generated.")
  }



  # write out / return results
  if (is.null (results)) final.results <- mod.t
  else final.results <- merge (results, mod.t, by=id.col)
  ##write.csv (final.results, paste (output.prefix, ".csv", sep=''), row.names=FALSE)



    ##invisible (final.results)
    return( list(input=d, output=final.results) )
##    list(input=d, output=final.results)
}



##############################################################################################
##
## Support functions
##


###########################################################################
##
##           Moderated t-test for significance testing
##
## code written by mani dr
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

    ###############
    # scatterplot 2D/3D
    ###############
    if( plot){
         require(scatterplot3d)

         par(mfrow=c(1,2))

         ## PC 1-2
         plot(pc1, pc2, xlab=paste("PC 1 (", round(100*comp.var[1]/sum(comp.var),1),"%)", sep=""), ylab=paste("PC 2 (",round(100*comp.var[2]/sum(comp.var),1),"%)", sep=""), pch=pch, main=main, col=col, sub=paste("Cumulative variance = ", round(100*sum(comp.var[1:2]/sum(comp.var)),1),"%", sep=""), cex=cex.points, ylim=c( min(pc2),  max(pc2)+.15*max(pc2)) )

         if(!is.null(leg.vec) & !is.null(leg.col))
             legend('top', legend=leg.vec, col=leg.col, pch=pch, pt.cex=max(1, cex.points-1.5), ncol=length(leg.vec))
             ##legend('top', legend=leg.vec, col=leg.col, pch=pch, pt.cex=cex.points, ncol=length(leg.vec))


         ## PC 1-3
         scatterplot3d( pca$x[,1], pca$x[,3], pca$x[,2], xlab=paste("PC 1 (", round(100*comp.var[1]/sum(comp.var),1),"%)", sep=""), zlab=paste("PC 2 (",round(100*comp.var[2]/sum(comp.var),1),"%)", sep=""), ylab=paste("PC 3 (", round(100*comp.var[3]/sum(comp.var),1),"%)", sep=""), color=col,  cex.symbols=cex.points, pch=pch, main=main, sub=paste("Cumulative variance = ", round(100*sum(comp.var[1:3]/sum(comp.var)),1),"%", sep=""), type="h" )

       par(mfrow=c(1,1))
    }

    return(pca)
}

##################################################
## color ramp
myColorRamp <- function(colors, values) {
    v <- (values - min(values))/diff(range(values))
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
