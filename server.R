library(shiny)


library(pheatmap)
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

##library (gplots)

## maxumum file size
options(shiny.maxRequestSize = 200*1024^2)

#################################################################
## global parameters
NASTRINGS <<- c("NA", "<NA>", "#NUM!", "#DIV/0!", "#NA", "#NAME?")
SEPARATOR <<- c(',', '\t', ';')
GRPCOLORS <<- RColorBrewer::brewer.pal(9, "Set1")
STRLENGTH <<- 15 ## number of characters to display in plots/tables for column names


###########################################################################################################
##                         Define server logic
###########################################################################################################
shinyServer(

    function(input, output, session) {

        #####################################
        ## reactive variables to store
        ## data accross a session
        #####################################

        ## test results
        global.results <-  reactiveValues()
        ## input data
        global.input <- reactiveValues()
        ## parameters
        global.param <-  reactiveValues(
            grp=NULL,           ## the actual group assignment
            N.grp=NULL,         ## number of defined groups
            ##grp.labels=NULL,    ## group label assignment
            grp.colors=NULL,    ## group color assignment
            grp.colors.legend=NULL, ## group colors, names are group names
            grp.done=F          ## group assignment finished?
        )
        ## coordinates in volcano plot
        volc <- reactiveValues()


        ################################################################################
        ##
        ##                                instructions
        ##
        ################################################################################

        #############################
        ## getting started
        output$help.start <- renderText({

            if(!is.null(input$file)) return()
            if(!is.null( global.input$file)) return()

            HTML( paste('<br><br><p><font size=\"5\">This Shiny App allows you to perform <b>one-sample</b> and <b>two-sample moderated T-tests</b> and to interactively explore the results of the analysis.
<br><br>The input table that you have to specify on the left has to be in the same format that is required to run Mani\'s R-scripts, i.e. they should contain only numeric expression values and a single id column and should be saved as text file (cvs or tab-delimited), e.g.</font></p><br><br>
<table style=\"width:50%\" align=\"center\" border=\"2\">
<tr><td>Unique ID</td><td>Exprs 1</td><td>Exprs 2</td><td>...</td><td>Exprs M</td></tr>
<tr><td>id 1</td><td>0.91</td><td>2.11</td><td>...</td><td>-1.53</td></tr>
<tr><td>id 2</td><td>1.77</td><td>0.31</td><td>...</td><td>3.13</td></tr>
<tr><td>id 3</td><td>-0.42</td><td>1.45</td><td>...</td><td>-0.61</td></tr>
<tr><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td></tr>
<tr><td>id N</td><td>-4.12</td><td>2.05</td><td>...</td><td>0.51</td></tr>
</table>
<br><br><p><font size=\"5\">Please note that right now the expression values in your data matrix should be <b>log-transformed</b>, especially in case of ratio data!</font></p>
<br><p><font size=\"5\">You can also load a <b>sample data set</b> by pressing the button on the left.</font></p>'

) )
        })
        ##############################
        ## help upload
        output$help.id.column <- renderText({

            if(is.null(input$id.col)) return()
            if(input$id.col > 0 && !is.null(input$id.col.value)) return()

            HTML( paste('<br><br><p><font size=\"5\">The table that you just uploaded contains the columns that are listed on the left side. Please specify the column containing unique IDs.</font></p>' ) )
        })

        ##############################
        ## help experimental design
        output$help.exp.design <- renderText({

            if(is.null(input$id.col)) return()
            if(input$id.col ==0) return()
            if(global.param$grp.done == T) return()

            HTML( paste('<br><br><p><font size=\"5\">Now you can specify your experimental design, i.e. assign columns of your table to groups that you want to compare. To do so please enter a name in the text field and select all columns that belong to that group. After pressing the \"Next\" button all selected columns will be removed from the list and you can continue to assign the remaining columns to another group. This process will go on until there are no unassigned columns left.<br><br>Please note that there is no limit in the number of groups you can define!</font></p>

<br><p align=\"center\"><font size=\"5\"><mark>Don\'t use any special characters (e.g. spaces, colons, etc.) in the group names!</mark></font></p>' ) )

        })

        ################################
        ## help test
        output$help.test <- renderText({

            if(global.param$grp.done == F) return()
            if(input$run.test > 0) return()

            HTML( paste('<br><br><p><font size=\"5\">You can choose between a one-sample and two-sample moderate T-test. If you select the one-sample test each group will be tested separately whether the group mean is significantly  different from zero. In case of a two-sample test each possible pairwise comparison will be performed and tested whether the group means are significantly different from each other.
<br><br>The moderate F-test is not supported yet (but will come soon...)!
</font></p>' ) )
        })

        ################################
        ## help results
        output$help.results <- renderText({

            if(global.param$grp.done == F) return()
            if(input$run.test == 0) return()

            HTML( paste('<p><font size=\"4\">This page allows you to interactively explore the results of you analyis. On the left you can choose between different filters, the results will be updated immediately. You can change the appearance of the heatmap by modifying the parameters below, you can select points shown in the Volcano plots and browse through the result table. Just feel free to play around and let me know if something is not working.
</font></p><br>') )
        })

        output$help <- renderText({

            HTML('<p align=\"center\"><font size=\"5\"><mark>To analyze another data set or to start over hit the F5 button.</mark></font></p>' )

        })

        ################################################################################
        ##
        ##                                navigation bar
        ##
        ################################################################################
        output$navbar <- renderUI({

            ##if(is.null(input$file)) return()
            if(is.null(input$file) && is.null( global.input$file)) return()
            ##if(is.null( global.input$file)) return()
            if(is.null(input$id.col)) return()
            if( !is.null(input$id.col))
                if( input$id.col == 0 ) return()
            if(is.null(global.results$data)) return()


            ## determine the number of group comparisons, e.g. for the number of volcano plots to draw
            groups.comp <- unique(global.param$grp.comp)

            volc.tabs <- list()
            volc.tabs[[1]] <- 'Volcanos'
            for(i in 1:length(unique(groups.comp))){
                volc.tabs[[i+1]]=tabPanel(paste0( groups.comp[i] ),

                                          fluidPage(
                                              fluidRow(
                                                  column(1, numericInput( paste("cex.volcano",groups.comp[i],sep='.'), "Point size", value=3, min=1, step=1)),
                                                  column(1, numericInput( paste("opac.volcano",groups.comp[i],sep='.'), "Opacity %", value=50, min=0, max=100, step=10)),
                                                  column(1, numericInput( paste("cex.volcano.lab",groups.comp[i],sep='.'), "Label size", value=1, min=.1, step=.1)),
                                                  column(1, selectInput( paste("grid.volcano",groups.comp[i],sep='.'), "Grid", c(T, F), selected=T)),

                                                  column(7),
                                                  column(1, h4('Export'))
                                              ),
                                              fluidRow(
                                                  column(11),
                                                  column(1, downloadButton(paste('downloadVolcano', groups.comp[i],sep='.'), 'Download this plot'))
                                              ),
                                              tags$br(),
                                              tags$hr(),
                                              tags$br(),
                                              fluidRow(
                                                  column(12, align='center', tableOutput(paste('info',groups.comp[i], sep='.')))
                                              ),
                                              fluidRow(
                                                  column(9, align='center',
                                                         plotOutput( paste("volcano",groups.comp[i], sep='.'), width=800, height=800, click=paste('plot_click', groups.comp[i], sep='.'), hover=paste('plot_hover', groups.comp[i], sep='.'))),
                                                  column(3, tableOutput(paste('volc.tab.selected', groups.comp[i], sep='.')))
                                              )
                                          )
                                        ) ## end tabPanel
            } ## end for i



            navbarPage('',

                           #######################################
                           ## HM
                           tabPanel('Heatmap',

                                    fluidPage(


                                        fluidRow(
                                            column(1, h4('Column labels')),
                                            column(1, h4('Row labels')),
                                            ##column(1, h4('Color key')),
                                            column(3),
                                            column(1, h4('Dimensions')),
                                            column(5),
                                            column(1, h4('Export'))
                                        ),
                                        fluidRow(
                                            column(1, numericInput( "cexCol", "Size", value=8, min=1, step=1)),
                                            column(1, numericInput( "cexRow", "Size", value=8, min=1, step=1)),
                                            column(1, selectInput( "scale", "Scale", c("row","column","none"), selected="none")),
                                            ##column(2, selectInput( "hm.trace", "Trace", c("column","row","both","none"), selected="none")),
                                            column(2, selectInput( "hm.clust", "Cluster", c("column","row","both","none"), selected="none")),
                                            column(1, numericInput( "pixRow", "Pixel rows", value=800, step=50)),
                                            column(1, numericInput( "pixCol", "Pixel columns", value=800, step=50)),
                                            column(4),
                                            column(1, downloadButton('downloadHM', 'Download this plot'))
                                            ##fluidRow(column(1, downloadButton('downloadHM', 'Download this plot')), column(11))
                                        ),
                                        ##fluidRow(
                                            ##column(1, numericInput( "srtCol", "Rotation", value=45, step=5)),
                                            ##column(1, numericInput( "srtRow", "Rotation", value=0, step=5)),
                                        ##    column(5),
                                        ##    column(1, numericInput( "pixCol", "Pixel columns", value=800, step=50)),
                                        ##    column(6)
                                        ##),
                                        tags$br(),
                                        tags$hr(),
                                        tags$br(),
                                        fluidRow(
                                            column(12, plotOutput("HM") )
                                        ),
                                        tags$br(),
                                        tags$hr(),
                                        tags$br()

                                    )

                                    ),
                       ###################################################################################
                       ## insert volcanos
                       do.call(navbarMenu, volc.tabs),


                           #######################################
                           ## PCA
                           tabPanel('PCA',
                                     fluidPage(
                                         fluidRow(
                                             column(12, tags$h3('Principle component analysis:'))
                                         ),
                                         fluidRow(
                                             column(11),
                                             column(1, h4('Export'))
                                         ),
                                         fluidRow(
                                             column(11),
                                             column(1, downloadButton('downloadPCA', 'Download this plot'))
                                         ),
                                         ##tags$br(),
                                         tags$hr(),
                                        ## tags$br(),
                                         fluidRow(
                                             column(12, plotOutput("pca", width=800, height=400) )
                                         ),
                                         fluidRow(
                                             ##column(12, plotOutput("pca", width=800, height=400) )
                                             column(12, plotlyOutput("pca.plotly", width=400, height=400) )
                                         )
                                    )
                                   ),
                           #######################################
                           ## table
                           tabPanel('Table',
                                    fluidPage(
                                        fluidRow(column(8, tags$h3('Result table (filtered):')), column(4, downloadButton('downloadTable', 'Download entire table')) ),
                                        fluidRow( column(12, tags$br())),
                                        fluidRow(column(12, dataTableOutput("tableprev")))
                                    )
                                    ),


                           #######################################
                           ## QC
                           navbarMenu( "QC",
                                      ## boxplots
                                      tabPanel('Boxplots',
                                               fluidPage(
                                                   fluidRow(
                                                       column(12, plotOutput("expr.boxplot"))
                                                       ##column(12, plotlyOutput("expr.boxplot.plotly"))
                                                   )
                                               )
                                               ),

                                      ## P-value distribution
                                      tabPanel('P-values',
                                               fluidPage(
                                                   fluidRow(
                                                        column(12, plotOutput("pval.hist"))
                                                        ##column(12, plotlyOutput("pval.hist.plotly"))
                                                   )
                                               )
                                               ),
                                      #########################################################
                                      ## correlation
                                      tabPanel('Multi scatter',
                                               fluidPage(
                                                   column(1, h3('Dimensions')),
                                                   column(10),
                                                   column(1, h3('Export'))

                                               ),
                                               fluidPage(
                                                   column(1, numericInput( "pixRowMS", "Pixel rows", value=1000, step=50)),
                                                   column(10),
                                                   column(1, downloadButton('downloadMS', 'Download this plot'))
                                               ),
                                               fluidPage(
                                                   column(1, numericInput( "pixColMS", "Pixel columns", value=1000, step=50)),
                                                   column(11)
                                               ),
                                               tags$br(),
                                               tags$hr(),
                                               tags$br(),
                                               fluidPage(
                                                   fluidRow(
                                                       ##column(12, plotOutput("expr.cor",width=1000, height=1000))
                                                       column(12, plotOutput("expr.cor"))
                                                   )
                                               )
                                               )

                                      )
                        ) ## end navbarpage
        })

        #########################################################################################
        ##
        ##                                  user input
        ##
        #########################################################################################

        ######################################
        ## file upload
        output$file.upload <- renderUI({

            if(!is.null(input$file)) return()
            if(!is.null( global.input$file)) return()
            list(
                HTML('<font size=\"3\"><b>Import data:</b></font>'),
                fileInput("file", "", accept=c('text/csv',
                       'text/comma-separated-values,text/plain',
                       '.csv', '.txt', '.tsv')),
                HTML('<hr border-width:\"10px\">')
            )

        })

        ######################################
        ## file testdata
        output$file.testdata <- renderUI({

            if(!is.null(input$file)) return()
            if(!is.null( global.input$file)) return()

            list(
                HTML('<br><br><p><font size=\"3\"><b>Use sample data:</b></font></p>'),
                actionButton("testdata", "Load test data set")
)

        })

        ######################################
        ## pick id column
        output$choose.id.column <- renderUI({

            ##if(is.null(input$file)) return()

            if(is.null( global.input$file) && is.null( input$file)) return()

            if(!is.null(input$id.col))
                if(input$id.col > 0 && !is.null(input$id.col.value)) return()

            tab <- global.input$table
            tab.colnames <- global.input$table.colnames


            ## radio button to pick id column
            list(
                radioButtons( "id.col.value", "Choose ID column", colnames(tab)),
                actionButton("id.col", 'OK'),
                ## experimental design
                HTML('<br><hr size=\"5\">'),
                downloadButton("exportTemplate", 'Export experimental design template')
            )

        })

        ######################################
        ## define groups
        output$define.groups <- renderUI({

            if(is.null(input$id.col)) return()
            if( !is.null(input$id.col))
                if( input$id.col == 0 ) return()
            if(!is.null(global.results$data)) return()
            if(!is.null(global.param$grp))
                if(sum(is.na(global.param$grp)) == 0) return()

            ## get groups
            tabhead <- global.param$grp
            tabhead <- names(tabhead[which(is.na(tabhead))])

            ## number of assigned groups
            N.grp = global.param$N.grp + 1

            list(
                ## upload template
                fileInput("exp.file", "Upload experimental design file", accept=c('text/plain','.txt')), HTML('<p align=\"center\"><font size=\"3\" color=\"blue\"><b>or define below</b></font></p>'),
                textInput(paste('label.grp', sep=''), paste('group', N.grp), value=paste('g', N.grp, sep='')),
                checkboxGroupInput( "groups", paste("Select group", N.grp), tabhead),
                actionButton( 'update.grp', 'Next')
            )
        })

        ######################################
        ## filter type
        output$filter.type <- renderUI({
            if( (is.null(input$file) && (is.null(global.input$file )))| is.null(input$run.test)) return()
            if(!is.null(input$run.test)) if(input$run.test == 0) return()

            list(selectInput('filter.type', 'Filter based on:', c('nom.p', 'adj.p', 'top.n', 'none'), selected='nom'))
        })

        #####################################
        ## select test
        output$list.groups <- renderUI({

            if( !global.param$grp.done ) return()
                list(
                    radioButtons('which.test', 'Select test', choices=c('One-sample mod T', 'Two-sample mod T', 'mod F'), selected='One-sample mod T'),
                    actionButton( 'run.test', 'Run test!')
            )
        })

        ################################################################################################
        ##
        ##                              calculations
        ##
        ################################################################################################

        ###############################################
        ## initialize group assignemnt
        observeEvent( input$id.col ,{

            if( is.null( global.input$table) | is.null(input$id.col.value) ) return()

            groups <- rep(NA, ncol(global.input$table))
            names(groups) <- colnames(global.input$table)

            ## remove id column
            groups <- groups[-c( which(input$id.col.value == names(groups))) ]

            ## set group assingment
            global.param$grp <- groups
            ## set number of assinged groups
            global.param$N.grp <- 0
        })
        ########################################
        ## update group assignment
        observeEvent( input$update.grp ,{

            label.grp.tmp = input$label.grp
            grp.tmp = input$groups

            ## current group assignment
            grp = global.param$grp
            ## update
            grp[ grp.tmp ] = label.grp.tmp
            global.param$grp <- grp

            ## check if done
            if(sum(is.na(grp)) == 0){
                global.param$grp.done = T
                ##global.param$grp.label = unique(grp)

                ## group colors
                grp.col <- rep(GRPCOLORS[1], length(grp))
                ##names(grp.col) <- names(grp)
                ##for(i in 2:global.param$N.grp) grp.col[ which(grp == unique(grp)[i]) ] <- GRPCOLORS[i]
                for(i in 2:length(unique(grp))) grp.col[ which(grp == unique(grp)[i]) ] <- GRPCOLORS[i]
                global.param$grp.colors <- grp.col
                ## group colors unique, e.g. to plot in a legend
                idx <- !duplicated(grp)
                grp.col.legend = grp.col[idx]
                names(grp.col.legend) <- grp[idx]
                global.param$grp.colors.legend <- grp.col.legend
            }
            ## update number of groups
            global.param$N.grp = length(unique( na.omit(grp)) )

        })
        ###############################################
        ## test data set
        observeEvent( input$testdata, {

            ##cat( getwd())
            ## cat('jaaaa')
            tab <- read.table( '/local/shiny-server/modT/testdata.txt', sep='\t', header=T, stringsAsFactors=F, na.strings=NASTRINGS  )

            global.input$table <- tab
            rm(tab)
            global.input$file  <- TRUE
        })


        ################################################
        ## upload file
        observeEvent( input$file, {

                ###########################################################
                ## determine the separator
                tab.sep=NULL
                ## try to figure out the separator
                for(s in SEPARATOR){
                    tab <- read.table(input$file$datapath, sep=s, header=T, stringsAsFactors=F, nrows=1)
                    if(length(tab) > 1){
                        global.param$tabsep <- s
                        break;
                    }
                }
                ###########################################################
                ## import the table
                tab <- read.table( input$file$datapath, sep=global.param$tabsep, header=T, stringsAsFactors=F, na.strings=NASTRINGS  )
                ids <- make.unique(as.character(tab[, input$id.col.value]), sep='_')
                tab[, input$id.col.value] <- ids

                ## store table
                global.input$table <- tab

                ## shorten column names and store together with the original names
                colnames.tmp <- chopString(colnames(tab), STRLENGTH)
                names(colnames.tmp) <- colnames(tab)
                global.input$table.colnames <- colnames.tmp
                rm(tab, colnames.tmp)
        })

        ##################################################
        ## upload experimental design file
        observeEvent( input$exp.file, {

            ## read the file
            grp.file <- read.delim(input$exp.file$datapath, header=T, stringsAsFactors=F)

            ## extract all non-empty cells in the 'Experiment' column
            grp.tmp <- grp.file[which(nchar(grp.file$Experiment) > 0 ), ]

            ## class vector
            grp=grp.tmp$Experiment
            names(grp)=grp.tmp$Column.Name

            ## update number of groups
            global.param$N.grp <- length(unique( na.omit(grp)) )
            ## store group assignment
            global.param$grp <- grp
            ## group colors
            grp.col <- rep(GRPCOLORS[1], length(grp))
            for(i in 2:length(unique(grp))) grp.col[ which(grp == unique(grp)[i]) ] <- GRPCOLORS[i]
            global.param$grp.colors <- grp.col
            ## group colors for figure legend
            idx <- !duplicated(grp)
            grp.col.legend = grp.col[idx]
            names(grp.col.legend) <- grp[idx]
            global.param$grp.colors.legend <- grp.col.legend

            ## all done
            global.param$grp.done = T

        })

        ################################################
        ## once the 'run test' button was pressed...
        observeEvent(input$run.test, {

            ## determine which test should be performed
            test = input$which.test

            ## specify which comparisons should be performed
            if(test == 'One-sample mod T'){
                ## each group separetely
                groups.comp <- global.param$grp
                global.param$grp.comp <- groups.comp
            }
            if(test == 'Two-sample mod T'){
                ## all pairwise combinations
                groups.unique <- unique(global.param$grp)

                groups.comp <- c()
                count=1
                for(i in 1:(length(groups.unique)-1))
                    for(j in (i+1):length(groups.unique)){
                        groups.comp[count] <- paste(groups.unique[i], groups.unique[j], sep='.vs.')
                        count <- count+1
                    }
                global.param$grp.comp <- groups.comp
            }

            #############################################
            ##   prepare the data matrix

            ## table
            tab <- global.input$table
            ## id column
            id.col = input$id.col.value
            ## all group labels
            groups=global.param$grp



            ##################################
            ## two sample
            if(test == 'Two-sample mod T'){

                withProgress(message='Two-sample test', value=0, {

                    count=0
                    ## loop over groups
                    for(g in unique(groups.comp)){

                        ## extract current groups
                        groups.tmp <- groups[grep(paste( unlist( strsplit(g, '\\.vs\\.')), collapse='|'  ) , groups) ]

                        ## progress bar
                        incProgress(1/length(unique(groups.comp)), detail=g)

                        ## extract table of current group
                        tab.group <- cbind(tab[, id.col], tab[, names(groups.tmp)])
                        colnames(tab.group)[1] <- id.col
                        res.tmp <-  modT.test.2class( tab.group, groups=groups.tmp, id.col=id.col, label=g )$output

                        if(count == 0)
                            res.comb <- res.tmp
                        else
                            ##res.comb <- merge(res.comb, res.tmp, by='id')
                            res.comb <- cbind(res.comb, res.tmp)

                        count=count + 1

                    }
                    global.results$data$output <- res.comb
                })
            }
            ##################################
            ## one sample
            if(test == 'One-sample mod T'){

                withProgress(message='One-sample test', value=0, {

                    count=0
                    ## loop over groups
                    for(g in unique(groups.comp)){

                        ## progress bar
                        incProgress(1/length(unique(groups.comp)), detail=g)

                        ## extract table of current group
                        tab.group <- cbind(tab[, id.col], tab[, names(groups)[which(groups == g)]])
                        colnames(tab.group)[1] <- id.col
                        res.tmp <- modT.test( tab.group, id.col=id.col, plot=F, nastrings=NASTRINGS, label=g)$output

                        if(count == 0)
                            res.comb <- res.tmp
                        else
                            res.comb <- merge(res.comb, res.tmp, by='id')

                        count=count + 1

                    }
                    ## add ids as rownames
                    rownames(res.comb) <- res.comb$id
                })
                global.results$data$output <- res.comb
            }
            ###################################################################
            ##            insert the panels for the volcanos
            ###################################################################
           ins.volc()
        })

        ############################################################################
        ##
        ##     filter the test results accross multiple groups
        ##
        ############################################################################
        filter.res  <-  reactive({

            groups.comp=unique(global.param$grp.comp)

            ## test results
            res <- global.results$data$output

            #################################
            ## top N
            if(input$filter.type=='top.n'){

                if(length(groups.comp) > 1){
                    ## order according to p-value
                    res <- res[ order( unlist(apply( res[, grep('^P.Value', colnames(res) )], 1, min))  ), ]
                    res <- res[ 1:input$top.n, ]

                    ## order according to FC
                    res <- res[order( unlist(apply( res[, grep('^logFC', colnames(res) )], 1, min))  ), ]

                    ## now separate for each group comparison
                    res.groups <- lapply(groups.comp, function(g){
                        res.filt=res[ which(!is.na(paste('P.Value.', g, sep='') )), ]
                        res.filt=res.filt[ order( res.filt[, paste('P.Value.', g, sep='') ], decreasing=F) , ]
                        res.filt[1:input$top.n, ]
                    })
                    names(res.groups) <- groups.comp

                } else {
                     res <- res[order(res[, paste('P.Value.', groups.comp, sep='')]) , ]
                     res <- res[1:input$top.n, ]
                     res <- res[order(res[, paste('logFC.', groups.comp, sep='')]) , ]

                     res.groups <- list(res)
                     names(res.groups) <- groups.comp
                }
                global.results$filter.cutoff <- input$top.n
            }
            #################################
            ## nominal p-value
            if(input$filter.type=='nom.p'){
                if(length(groups.comp) > 1){
                    res <- res[ which( unlist(apply( res[, grep('^P.Value', colnames(res) )], 1, function(x) sum(x < input$p.val))) > 0), ]
                    res <- res[order( unlist(apply( res[, grep('^logFC', colnames(res) )], 1, min))  ), ]
                    ## now separate for each group comparison
                    res.groups <- lapply(groups.comp, function(g){
                        res[ which( res[, paste('P.Value.', g, sep='')] < input$p.val) , ]
                    })
                    names(res.groups) <- groups.comp
                } else {
                    res <- res[which(res[, paste('P.Value.', groups.comp, sep="")] < input$p.val), ]
                    res <-  res[order(res[, paste('logFC.', groups.comp, sep="")]), ]
                    res.groups <- list(res)
                    names(res.groups) <- groups.comp
                }
                global.results$filter.cutoff <- input$p.val
            }

            #################################
            ## adjusted p-value
            if(input$filter.type=='adj.p'){
                if(length(groups.comp) > 1){
                    res <- res[ which( unlist(apply( res[, grep('^adj.P.Val', colnames(res) )], 1, function(x) sum(as.numeric(x) < input$adj.p ))) > 0), ]
                    res <- res[order( unlist(apply( res[, grep('^logFC', colnames(res) )], 1, min))  ), ]
                    ## now separate for each group comparison
                    res.groups <- lapply(groups.comp, function(g){
                         res[ which( res[, paste('adj.P.Val.', g, sep='')] < input$adj.p) , ]
                    })
                    names(res.groups) <- groups.comp
                } else {
                    res <- res[which(res[, paste('adj.P.Val.', groups.comp,sep='')] < input$adj.p), ]
                    res <-  res[order(res[, paste('logFC.', groups.comp, sep="")]), ]
                    res.groups <- list(res)
                    names(res.groups) <- groups.comp
                }
                global.results$filter.cutoff <- input$adj.p
            }

            ###################################
            ## global FDR

            #################################
            ## no filter
            if(input$filter.type=='none'){
                global.results$filter.cutoff <- 'none'
                res.groups <- lapply(groups.comp, function(g) res)
                names(res.groups) <- groups.comp
            }

            ###################################################
            ## global filter accross all experiments
            global.results$filtered <- res
            global.results$filter.type <- input$filter.type
            global.results$filtered.groups <- res.groups

        })

        ###################################################################################################
        ##
        ##                                     output
        ##
        ###################################################################################################

        ##############################
        ## download table
        output$downloadTable <- downloadHandler(
            filename = function(){ paste("results_", sub(' ', '_',input$which.test), '_', sub(' .*', '', Sys.time()),".txt", sep='') },
            content = function(file){
                tab <- global.results$data$output
                colnames(tab) <- sub('^X','',colnames(tab))
                write.table(  tab, file, sep='\t', quote=F, na='', row.names=F  )
            }
        )

        #########################################
        ## download experimental design template
        output$exportTemplate <- downloadHandler(
            filename = function(){ 'experimentalDesign.txt' },
            content = function(file){
                tab <- global.input$table
                exp.design <- cbind(colnames(tab), rep('', ncol(tab)))
                colnames(exp.design) <- c('Column.Name', 'Experiment')
                write.table(  exp.design, file, sep='\t', quote=F, na='', row.names=F  )
            }
        )


        #####################################################################################
        ##
        ##             display the corresponding part of the table
        ##
        #####################################################################################
        output$tableprev <- renderDataTable({

            if(is.null(global.results$data)) return()

            filter.res()

            tab <- global.results$filtered
            colnames(tab) <- sub('^X','',colnames(tab))

            if(nrow(tab) > 0){
                ## add links to uniprot
                up.id <- tab[, 'id']
                up.link <- paste("<a href='http://www.uniprot.org/uniprot/", sub('(_|,|;).*', '', up.id),"' target='_blank'>", up.id, "</a>", sep='')
                tab[, 'id'] <- up.link
            }
            tab

        }, options = list( pageLength = 50), escape=F)


        #####################################################################################
        ##
        ##                             Volcano plot
        ##
        #####################################################################################

        ###################################################################
        ## function to generate the panels for the volcanos
        ## insert the plots into the webpage
        ###################################################################
        ins.volc <- reactive({

            grp.comp <- unique( global.param$grp.comp )

            for(i in 1:length(grp.comp)){
              local({
                  my_i <- i
                  ##########################
                  ## the actual plots
                  output[[paste("volcano", grp.comp[my_i], sep='.')]] <- renderPlot({
                      plotVolcano( grp.comp[my_i] )
                  })

                  ##########################
                  ## download button
                  output[[paste("downloadVolcano", grp.comp[my_i], sep='.')]] <- downloadHandler(
                      ##filename = paste('volcano_', global.results$filter.type, '_', global.results$filter.cutoff, '.pdf', sep=''),
                      filename =  paste( 'volcano_',grp.comp[my_i],'.pdf'),
                      content = function(file){

                          pdf(file, height=11, width=11)

                          for(j in 1:length(grp.comp)){
                              local({
                                  my_j=j
                                  plotVolcano(grp.comp[my_j])
                              })
                          }
                          dev.off()
                      }
                  ) ## end download handler

                  ##################################
                  ## info table
                  output[[paste('info', grp.comp[my_i], sep='.')]] <-  renderTable({
                      if(is.null(global.results$data)) return()
                      res <- as.data.frame( global.results$data$output )
                      text.tmp <- nearPoints(res, input[[paste('plot_hover', grp.comp[my_i], sep='.')]], threshold=10, maxpoints =  1, xvar=paste('logFC', grp.comp[my_i], sep='.'), yvar=paste('Log.P.Value', grp.comp[my_i], sep='.'))
                      text.tmp <- text.tmp[, c('id', paste('logFC', grp.comp[my_i], sep='.'), paste('P.Value',grp.comp[my_i], sep='.'), paste('adj.P.Val', grp.comp[my_i], sep='.'))]
                      rownames(text.tmp) <- NULL
                      if( nrow(text.tmp) == 1 )
                          text.tmp
                      else{
                          text.tmp[1, ] <- rep(' ', ncol(text.tmp))
                          text.tmp
                      }
                  })

                  ##################################
                  ## observe clicks
                  observeEvent(input[[paste('plot_click', grp.comp[my_i], sep='.')]], {
                      res <- as.data.frame( global.results$data$output )
                      ##group.comp <- unique(global.param$grp.comp)

                      text.tmp <- nearPoints(res, input[[paste('plot_click', grp.comp[my_i], sep='.')]], threshold=10, maxpoints =  1, xvar=paste('logFC', grp.comp[my_i], sep='.'), yvar=paste('Log.P.Value', grp.comp[my_i], sep='.'))

                      if(nrow(text.tmp) == 1){
                          ################################################
                          ## first click
                          if(is.null(volc[[ paste('x', grp.comp[my_i], sep='.')]] )){
                              volc[[paste('x', grp.comp[my_i], sep='.')]] = text.tmp[paste('logFC', grp.comp[my_i], sep='.')]
                              volc[[paste('y', grp.comp[my_i], sep='.')]] = text.tmp[paste('Log.P.Value', grp.comp[my_i], sep='.')]
                              volc[[paste('text', grp.comp[my_i], sep='.')]] = text.tmp['id']
                              volc[[paste('xy', grp.comp[my_i], sep='.')]] = paste(text.tmp[paste('logFC', grp.comp[my_i], sep='.')], text.tmp[paste('Log.P.Value', grp.comp[my_i], sep='.')])
                              volc[[paste('P.Value', grp.comp[my_i], sep='.')]] <- text.tmp[paste('P.Value', grp.comp[my_i], sep='.')]
                              volc[[paste('adj.P.Val', grp.comp[my_i], sep='.')]] <- text.tmp[paste('adj.P.Val', grp.comp[my_i], sep='.')]
                          } else {
                              ######################################################
                              ## REMOVE: check if point is present already
                              if( paste(text.tmp[paste('logFC', grp.comp[my_i], sep='.')], text.tmp[paste('Log.P.Value', grp.comp[my_i], sep='.')]) %in% volc[[paste('xy', grp.comp[my_i], sep='.')]]){

                                  ## if so remove from the list
                                  idx = which( volc[[paste('xy', grp.comp[my_i], sep='.')]] == paste(text.tmp[paste('logFC', grp.comp[my_i], sep='.')], text.tmp[paste('Log.P.Value', grp.comp[my_i], sep='.')]))

                                  if(length(volc[[paste('xy', grp.comp[my_i], sep='.')]] > 1)){
                                      volc[[paste('x', grp.comp[my_i], sep='.')]] <- volc[[paste('x', grp.comp[my_i], sep='.')]][-idx]
                                      volc[[paste('y', grp.comp[my_i], sep='.')]] <- volc[[paste('y', grp.comp[my_i], sep='.')]][-idx]
                                      volc[[paste('text', grp.comp[my_i], sep='.')]] <- volc[[paste('text', grp.comp[my_i], sep='.')]][-idx]
                                      volc[[paste('xy', grp.comp[my_i], sep='.')]] <- volc[[paste('xy', grp.comp[my_i], sep='.')]][-idx]
                                      volc[[paste('P.Value', grp.comp[my_i], sep='.')]] <- volc[[paste('P.Value', grp.comp[my_i], sep='.')]][-idx]
                                      volc[[paste('adj.P.Val', grp.comp[my_i], sep='.')]] <- volc[[paste('adj.P.Val', grp.comp[my_i], sep='.')]][-idx]
                                  } else {
                                      volc[[paste('text', grp.comp[my_i], sep='.')]] <- volc[[ paste('y', grp.comp[my_i], sep='.') ]] <- volc[[paste('x', grp.comp[my_i], sep='.')]] <- volc[[paste('xy', grp.comp[my_i], sep='.')]] <- volc[[paste('adj.P.Val', grp.comp[my_i], sep='.')]] <- volc[[paste('P.Value', grp.comp[my_i], sep='.')]]<- NULL
                                                    }
                                  ################################################
                                  ## ADD: if selected point is not present add it to the list
                              } else{
                                  volc[[paste('x', grp.comp[my_i], sep='.')]]=c( volc[[paste('x', grp.comp[my_i], sep='.')]],
                                                                      text.tmp[paste('logFC', grp.comp[my_i], sep='.')])
                                  volc[[paste('y', grp.comp[my_i], sep='.')]]=c(volc[[paste('y', grp.comp[my_i], sep='.')]], text.tmp[paste('Log.P.Value', grp.comp[my_i], sep='.')])
                                  volc[[paste('text', grp.comp[my_i], sep='.')]]=c(volc[[paste('text', grp.comp[my_i], sep='.')]],  text.tmp[ 'id'] )
                                  volc[[paste('xy', grp.comp[my_i], sep='.')]] = c(volc[[paste('xy', grp.comp[my_i], sep='.')]], paste(text.tmp[paste('logFC', grp.comp[my_i], sep='.')], text.tmp[paste('Log.P.Value', grp.comp[my_i], sep='.')]) )

                                  volc[[paste('P.Value', grp.comp[my_i], sep='.')]]=c(volc[[paste('P.Value', grp.comp[my_i], sep='.')]],  text.tmp[paste('P.Value', grp.comp[my_i], sep='.')] )
                                  volc[[paste('adj.P.Val', grp.comp[my_i], sep='.')]]=c(volc[[paste('adj.P.Val', grp.comp[my_i], sep='.')]],  text.tmp[paste('adj.P.Val', grp.comp[my_i], sep='.')] )
                              }
                          }
                      }
                  }) ## end observe clicks

                  ######################################################
                  ## table of selected features
                  output[[paste('volc.tab.selected', grp.comp[my_i], sep='.')]] <- renderTable({

                      if(is.null(volc[[paste('x', grp.comp[my_i], sep='.')]])) return()
                      if(length(volc[[paste('x', grp.comp[my_i], sep='.')]]) == 0) return()
                      tags$h4('Selection:')

                      id.tmp <- volc[[paste('text', grp.comp[my_i], sep='.')]]
                      ## dat.select = data.frame(id=unlist(volc[[paste('text', grp.comp[my_i], sep='.')]]), logFC=unlist(volc[[paste('x', grp.comp[my_i], sep='.')]]), xy=unlist(volc[[paste('xy', grp.comp[my_i], sep='.')]]))
                       dat.select = data.frame(id=unlist(volc[[paste('text', grp.comp[my_i], sep='.')]]), logFC=unlist(volc[[paste('x', grp.comp[my_i], sep='.')]]), P.Value=unlist(volc[[paste('P.Value', grp.comp[my_i], sep='.')]]), adj.P.Value=unlist(volc[[paste('adj.P.Val', grp.comp[my_i], sep='.')]]) )
                      up.id <- dat.select[, 'id']
                      up.link <- paste("<a href='http://www.uniprot.org/uniprot/", sub('(_|,|;).*', '', up.id),"' target='_blank'>", up.id, "</a>", sep='')
                      dat.select[, 'id'] <- up.link

                      dat.select
                  }, sanitize.text.function = function(x) x)



              }) ## end local

            } ## end for loop

        })

        ############################################
        ## actual plot
        plotVolcano <- function(group){

            ##grp.comp <- unique( global.param$grp.comp )
            ##group <-  grp.comp[i]

            filter.res()

            ## pch for significant points
            sig.pch=18

            ## unfiltered
            res = as.data.frame( global.results$data$output )
            ##rownames(res) <- res[, input$id.col.value]

            ## filtered
            res.filt = as.data.frame(global.results$filtered.groups[[group]])
            ##rownames(res.filt) <- res.filt$id

            ## extract fc and p
            logFC <- res[, paste('logFC.', group, sep='')]
            logPVal <- res[, paste('Log.P.Value.', group, sep='')]

            ## index of missing values
            rm.idx <- union( which(is.na(logFC)), which(is.na(logPVal)) )
            if(length(rm.idx) > 0){
                res <- res[-rm.idx, ]
                logFC <- logFC[-rm.idx]
                logPVal <- logPVal[-rm.idx]
            }

            ## which filter?
            filter.str <- paste('filter:', global.results$filter.type, '\ncutoff:', global.results$filter.cutoff)

            ## filter
            if(global.results$filter.type == 'top.n'){
                PVal <- res[, paste('P.Value.', group, sep='')]
                sig.idx = order(PVal, decreasing=F)[1:global.results$filter.cutoff]
            }
            if(global.results$filter.type == 'nom.p'){
                PVal <- res[, paste('P.Value.', group, sep='')]
                sig.idx = which(PVal <= global.results$filter.cutoff)
            }
            if(global.results$filter.type == 'adj.p'){
                adjPVal <- res[, paste('adj.P.Val.', group, sep='')]
                sig.idx = which(adjPVal <= global.results$filter.cutoff)
            }
            if(global.results$filter.type == 'none')
                sig.idx = 1:length(logFC)

            ## colors
            ##col <- rep(my.col2rgb('black', (input[[paste('opac.volcano',  group, sep='.')]]*255)/100), nrow(res))
            ##names(col) = rownames(res)

            pch.vec=rep(19, nrow(res))
            cex.vec=rep( input[[paste('cex.volcano', group, sep='.')]], nrow(res))

            if(length(sig.idx) > 0){
                pch.vec[sig.idx] <- sig.pch
                cex.vec[sig.idx] <- cex.vec[1]+1
            }

            ##    col[sig.idx] <- my.col2rgb('darkred', (input[[paste('opac.volcano', group, sep='.')]]*255)/100 )

            ## remove NA



            col=myColorRamp(c('black', 'grey20', 'darkred', 'red', 'deeppink'), na.omit(logPVal))


            ## limits
            xlim = max(abs(logFC), na.rm=T)
            xlim = xlim + xlim*.1
            xlim = c(-xlim, xlim)

            ylim = max(logPVal, na.rm=T)
            ylim = c(0, ylim+.2*ylim)

            ## plot
            par(mar=c(4,5,5,2))
            plot.new()
            plot.window( xlim=xlim, ylim=ylim, cex.axis=1.8, cex.lab=1.8, main=group, cex.main=2)
            mtext(expression(log(FC)), side=1, cex=1.8, line=3)
            mtext(expression(-10*log[10](P-value)), side=2, cex=1.8, line=3)
            axis(1, cex.axis=1.8)
            axis(2, las=2, cex.axis=1.8)

            if( input[[paste('grid.volcano', group, sep='.')]] )
                grid()
            ##abline(h=-10*log(global.results$filter.cutoff, 10), col='grey', lwd=2, lty='dashed')
            ##plot(logFC, logPVal, col=col, xlim=xlim, ylim=ylim, pch=pch.vec, cex=cex.vec, xlab=expression(log(FC)), ylab=expression(-10*log[10](P-value)), cex.axis=1.8, cex.lab=1.8, main=group, cex.main=2)
            points(logFC, logPVal, col=col, pch=pch.vec, cex=cex.vec)
            ## filter
            abline(h=min(logPVal[sig.idx], na.rm=T), col='grey30', lwd=2, lty='dashed')
            text( xlim[2]-(xlim[2]*.05), min(logPVal[sig.idx], na.rm=T), paste(global.results$filter.type, global.results$filter.cutoff, sep='='), pos=3, col='grey30')
            ##text( xlim[2], min(logPVal[sig.idx], na.rm=T), paste(global.results$filter.type, global.results$filter.cutoff, sep='='), adj=c(-1,-1), col='grey')

            legend('top', bty='n', legend=paste(filter.str, '\nN=', sum(!is.na(logFC) & !is.na(logPVal)), sep=''), cex=1.5)
            ## up/down
            ##legend('topleft', legend=paste('down:', sum()))

            if(!is.null( volc[[paste('x', group, sep='.')]] ) & length(volc[[paste('x', group, sep='.')]]) ){
                for(i2 in 1:length(unlist(volc[[paste('x', group, sep='.')]])))
                   text(unlist(volc[[paste('x', group, sep='.')]][i2]), unlist(volc[[paste('y', group, sep='.')]][i2]), unlist(volc[[paste('text', group, sep='.')]][i2]),pos=ifelse(volc[[paste('x', group, sep='.')]][i2] < 0, 2, 4), cex=input[[paste('cex.volcano.lab', group, sep='.')]])
            }
            ##if( input[[paste('grid.volcano', group, sep='.')]] )
            ##    grid()
        }

        #######################################################################################
        ##
        ##                                    boxplots
        ##
        #######################################################################################
        output$expr.boxplot <- renderPlot({

            if(is.null(global.results$data)) return()

            ## dataset
            tab <- global.input$table
            ## id column
            id.col <- input$id.col.value
            ## table
            tab <- tab[, setdiff(colnames(tab), id.col)]
            ## get groups
            ## grp <- global.results$data$groups
            grp <- global.param$grp
            ##grp.lab <- global.param$grp.label

            grp.col <- global.param$grp.colors
            grp.col.leg <- global.param$grp.colors.legend
            ## mapping to colors
            ##grp.col <- rep('grey10', length(grp))
            ##grp.col[which(grp == input$label.g2)] <- 'darkblue'

            ##########################################
            ## order after groups
            ord.idx <- order(grp)
            grp <- grp[ord.idx]
            tab <- tab[, ord.idx]
            grp.col <- grp.col[ ord.idx]


            at.vec=1:ncol(tab)
            ##########################################
            ## plot
            par(mar=c(4,15,2,6))
            boxplot(tab, pch=20, col='white', outline=T, horizontal=T, las=2, xlab=expression(log[2](ratio)), border=grp.col, at=at.vec, axes=F, main='', cex=2, xlim=c(0, ncol(tab)+2))
            ##legend('top', legend=names(grp.col.leg), ncol=2, bty='n', border = names(grp.col.leg), fill='white', cex=1.5)
            legend('top', legend=names(grp.col.leg), ncol=length(grp.col.leg), bty='n', border = grp.col.leg, fill=grp.col.leg, cex=1.5)
            ##legend('top', legend=c(input$label.g1, input$label.g2), ncol=2, bty='n', border = c('grey10', 'darkblue'), fill='white', cex=1.5, lwd=3)
            mtext( paste('N=',unlist(apply(tab,2, function(x)sum(!is.na(x)))), sep=''), at=at.vec, side=4, las=2, adj=0, cex=.8)
            axis(1)
            axis(2, at=at.vec, labels=chopString(colnames(tab), STRLENGTH), las=2)

        },
        width = function(){ width=1000},
        height= function(){ height=20*ncol( global.input$table) } )

        ####################################################################
        ## plotly boxplots
        output$expr.boxplot.plotly <- renderPlotly({

            if(is.null(global.results$data)) return()

            ## dataset
            tab <- global.input$table
            ## id column
            id.col <- input$id.col.value
            ## table
            tab <- tab[, setdiff(colnames(tab), id.col)]
            ## get groups
            ## grp <- global.results$data$groups
            grp <- global.param$grp
            ##grp.lab <- global.param$grp.label

            grp.col <- global.param$grp.colors
            grp.col.leg <- global.param$grp.colors.legend

            ##########################################
            ## convert to a data frame
            exprs = unlist(tab)
            column = group = rep('', length(exprs))
            for(i in 1:ncol(tab)){
                column[ max(1+((i-1)*nrow(tab))) : (((i)*nrow(tab))) ] = colnames(tab)[i]
            }
            for(i in 1:length(grp)){
                group[which(column == names(grp)[i])] <- grp[i]
            }

            tab2 = data.frame(exprs, column, group)
            rm(exprs, column, group)

            ##########################################
            ## plot
            ##plot_ly(tab, x=exprs, group=column, type='box')
            ##bp <- plot_ly(tab, x=exprs, color=group, type='box')
            ## bp.l <- layout(bp, boxmode=column)

            ##p=plot_ly(tab2, x=exprs, group=column, type='box', colors=grp.col) %>% layout(boxmode='group')
            plot_ly(tab2, x=exprs, type='box', group=column, marker=list( color=grp.col) )
            ##subplot(p)
        })

        ######################################################################################
        ##
        ##                                   Correlation
        ##
        ######################################################################################
        output$expr.cor <- renderPlot({
            if(is.null(global.results$data)) return()
            ##plotMultiScatter()
            plotCorrMat(lower='pearson', upper='spearman')
        },
        width = function(){ input$pixColMS},
        height= function(){ input$pixRowMS}
        )

        ###############################
        ## actual plot
        plotMultiScatter <- function(){
            ## dataset
            tab <- global.input$table
            ## id column
            id.col <- input$id.col.value
            ## table
            tab <- tab[, setdiff(colnames(tab), id.col)]
            ## get groups
            grp <- global.results$data$groups
            ## mapping to colors
            grp.col <- rep('grey10', length(grp))
            grp.col[which(grp == input$label.g2)] <- 'darkblue'

            colnames(tab) <- chopString(colnames(tab), STRLENGTH)
            ###############################
            ## plot
            withProgress({
                setProgress(message = 'Processing...', detail= 'Calculation correlations')
                my.multiscatter(tab)
            })
        }
        ## download image, Multiscatter
        output$downloadMS <- downloadHandler(
            filename =  paste( 'multiscatter.pdf'),
            content = function(file){
                pdf(file, height=input$pixRowMS*(11/800), width=input$pixColMS*(11/800))
                plotMultiScatter()
                dev.off()
            }
        )

        ###################################################
        ## correlation matrix
        ###################################################
        plotCorrMat <- function(filename=NA, lower=c('pearson', 'spearman', 'kendall', 'pcor'), upper=c('pearson', 'spearman', 'kendall', 'pcor')){

            ## dataset
            tab <- global.input$table
            ## id column
            id.col <- input$id.col.value
            ## class vector
            grp <- sort(global.param$grp)
            grp.col.legend <- global.param$grp.colors.legend

            ## table
            tab <- tab[, setdiff(colnames(tab), id.col)]
            tab <- tab[, names(grp)]

            ###########################
            ## calculate correlations
            cm.upper <- cor(tab, use='pairwise', method=match.arg(upper))
            cm.lower <- cor(tab, use='pairwise', method=match.arg(lower))

            ###########################
            ## initialize correlation matrix
            cm <- matrix(NA, ncol=ncol(cm.upper),nrow=nrow(cm.upper), dimnames=dimnames(cm.upper))
            cm[ lower.tri(cm, diag=F) ] <- cm.lower[lower.tri(cm.lower, diag=F)]
            cm[ upper.tri(cm, diag=F) ] <- cm.upper[upper.tri(cm.upper, diag=F)]


            col=colorRampPalette(c('blue', 'grey80' , 'red'))

            ## annotation of rows/columns
            anno=data.frame(Group=grp)
            anno.color=list(Group=grp.col.legend)

            Rowv=F
            Colv=F

            ##setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=0.9, name="vp", just=c("right","top"))), action="prepend")
            setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=0.9, name="vp", just=c("right","top"))))
            pheatmap(cm, fontsize_row=10, fontsize_col=10,
                     cluster_rows=Rowv, cluster_cols=Colv, border_col=NA, col=col(100), filename=filename, labels_col=chopString(colnames(cm), STRLENGTH), labels_row=chopString(rownames(cm), STRLENGTH), main='', annotation_col=anno, annotation_colors=anno.color,  annotation_row=anno, display_numbers=T, fontsize_number=100/ncol(cm)+10)
            setHook("grid.newpage", NULL, "replace")
            grid.text(paste(match.arg(lower)), y=.995, x=.4, gp=gpar(fontsize=25))
            grid.text(paste(match.arg(upper)), x=-0.01, rot=90, gp=gpar(fontsize=25))
            ##popViewport()

        }

        ####################################################################################
        ##
        ##                                  Heatmap
        ##
        ####################################################################################
        output$HM <- renderPlot({
            if(is.null(global.results$data)) return()
            plotHM()
        },
        width = function(){ input$pixCol},
        height= function(){ input$pixRow}
        )
        #######################################
        ## actual plot
        plotHM <- function(filename=NA){

            ## groups
            ## grp = global.results$data$groups
            grp=global.param$grp
            grp.col <- global.param$grp.colors
            grp.col.legend <- global.param$grp.colors.legend

            filter.res()
            res = global.results$filtered
            ##rownames(res)

            ## extract expression values
            res = res[, names(grp)]

            ## which filter has been used
            filter.str <- paste('filter:', global.results$filter.type, '\ncutoff:', global.results$filter.cutoff)

            if(nrow(res) < 3) return()
            colnames(res) <- sub('^X','', colnames(res) )

            res <- data.matrix(res)

            ## reordering of columns/rows
            if(input$hm.clust == 'column'){
                Rowv=FALSE
                ##colv.dist = dist(res, method='euclidean')
                ##colv.clust = hclust( colv.dist, method='complete')
                ##Colv=as.dendrogram( colv.clust  )
                Colv=TRUE
            } else if(input$hm.clust == 'row'){
                Rowv=TRUE
                ##Rowv=as.dendrogram( hclust( dist(res, method='euclidean'), 'complete') )
                Colv=FALSE
            } else if(input$hm.clust == 'both'){
                Rowv=Colv=TRUE
                ##Rowv=Colv=as.dendrogram( hclust( dist(res, method='euclidean'), 'complete') )
            } else {
                Rowv=Colv=FALSE
            }

           ############################################
           ## heatmap.2
           ## heatmap.2(data.matrix(res), Rowv=Rowv, Colv=Colv, dendrogram=input$hm.clust, scale='row', col=rev(brewer.pal (11, "RdBu")), margins=c(15,15), density.info='density', keysize=input$keysize, cexRow=input$cexRow, cexCol=input$cexCol, trace=input$hm.trace, main=paste(filter.str, '\nN=', nrow(res),sep=''), key.title='', srtCol=input$srtCol, srtRow=input$srtRow, ColSideColors=grp.col, key=F)


           ##############################################
           ## annotation of columns
           anno.col=data.frame(Group=grp)
           anno.col.color=list(Group=grp.col.legend)


           ############################################
           ## heatmap
           pheatmap(res, fontsize_row=input$cexRow, fontsize_col=input$cexCol,
                     cluster_rows=Rowv, cluster_cols=Colv, border_col=NA, col=rev(brewer.pal (11, "RdBu")), filename=filename, main=paste(filter.str, '\nN=', nrow(res),sep=''), annotation_col=anno.col, annotation_colors=anno.col.color, labels_col=chopString(colnames(res), STRLENGTH))
        }
        ## download image
        output$downloadHM <- downloadHandler(
            filename =  paste( 'heatmap.pdf'),
            content = function(file){
                ##pdf(file, height=input$pixRow*(11/800), width=input$pixCol*(11/800))
                plotHM(filename=file)
                ##dev.off()
            }
        )

        #####################################################################################
        ## histogram of p-values
        #####################################################################################
        output$pval.hist <- renderPlot({

            if(is.null(global.results$data)) return()
            groups.comp <- unique(global.param$grp.comp)

            res = global.results$data$output

            par(mfrow=c(length(groups.comp),1))
            for(g in groups.comp){

                pval <- res[, paste('P.Value', g, sep='.')]
                hist(pval, breaks=50, main=paste('Histogram of P-values (N=', sum(!is.na(pval)), ')',sep=''), xlab='P-value', cex.main=3, cex.axis=2, cex.lab=2, col='darkblue', border=NA)
                legend('top', legend=g, cex=2)
            }
        },
        width = function(){ width=1000},
        height= function(){ height=500*length(unique(global.param$grp.comp))} )



        ################################################################################
        ## Histogram PLOTLY
        output$pval.hist.plotly <- renderPlotly({

            if(is.null(global.results$data)) return()
            groups.comp <- unique(global.param$grp.comp)

            res = global.results$data$output

            ##par(mfrow=c(length(groups.comp),1))
            for(g in groups.comp){

                pval <- res[, paste('P.Value', g, sep='.')]


                ##hist(pval, breaks=50, col='grey', main=paste('Histogram of P-values\nN=', sum(!is.na(pval)), sep=''), xlab='P-value')
                ##legend('top', legend=g, cex=2)
            }
            plot_ly(x=pval, type = 'histogram')

        })


        ######################################################################################
        ##
        ##                                 PCA
        ##
        ######################################################################################
        output$pca <- renderPlot({
            if(is.null(global.results$data) | is.na(input$top.n)) return()
            pca=plotPCA()

        })

        ##################################
        ## plot_ly
        output$pca.plotly <- renderPlotly({
            if(is.null(global.results$data) | is.na(input$top.n)) return()
            pca=plotPCA(plot=F)

            pca.mat = data.frame(
                PC1=pca$x[,1],
                PC2=pca$x[,2]
            )
            rownames(pca.mat) <- rownames(pca$x)

            plot_ly( pca.mat, x=PC1, y=PC2, type='scatter', mode='markers', marker=list(size=20, color=global.param$grp.colors), text=chopString(rownames(pca.mat), STRLENGTH) )

        })

        #########################
        ## actual plot
        plotPCA <- function(plot=T){

            filter.res()

            res <- global.results$filtered

            if(nrow(res) < 3) return()

            ## get groups
            grp <- global.param$grp
            ## mapping to colors
            grp.col <- global.param$grp.colors
            grp.col.leg <- global.param$grp.colors.legend

            ## remove missing values
            rm.idx <- apply(res, 1, function(x) sum(is.na(x)) + sum(is.infinite(x)))
            rm.idx <- which(rm.idx > 0)
            if(length(rm.idx)>0) res <- res[-rm.idx, ]
            if(nrow(res) < 3) return()

            ## plot
            pca <- my.prcomp(t(res[, names(grp)]), col=grp.col, plot=plot, rgl=F, main='', cex.points=5, leg.vec=names(grp.col.leg), leg.col=grp.col.leg)
            return(pca)
        }
        ## download PCA
        output$downloadPCA <- downloadHandler(
            filename = paste('pca_', global.results$filter.type, '_', global.results$filter.cutoff, '.pdf', sep=''),
            ##filename =  paste( 'pca.pdf'),
            content = function(file){
                pdf(file, height=11, width=22)
                pca=plotPCA()
                dev.off()
            }
        )
})


