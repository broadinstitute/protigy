################################################################################################################
## Filename: server.r
## Created: October 09, 2015
## Author(s): Karsten Krug, Ozan Aygun
##
## Purpose: Shiny-app to perform differential expression analysis, primarily on proteomics data, to perform
##          simple data QC, to interactively browse through the results and to download high-quality result
##          figures.
##
## This file defines the server logical of the app. It also takes care of most of the user interface.
##
################################################################################################################
p_load(shiny)

############################################
## set maximum file size for upload
############################################
options(shiny.maxRequestSize = MAXSIZEMB*1024^2, stringsAsFactors=F)



###########################################################################################################
##                         Define server logic
############################################################################################
shinyServer(

    function(input, output, session) {

        hide(id = "loading-content", anim = TRUE, animType = "fade")

        ##show("app-content")

        #####################################
        ## reactive variables to store
        ## data accross a session
        #####################################

        ## error messages
        error <- reactiveValues()
        
        # flags
        
        
        ## test results
        global.results <-  reactiveValues(

            ## data tables
            data=NULL,
            table.norm=NULL,
            table.log=NULL,
            table.repro.filt=NULL,
            filtered=NULL,

            export.results=F,
            pca=NULL,
            repro.filt=NULL,

            ## ids and gene names
            ids=NULL,
            gene.names=NULL
        )
        
        ## input data
        global.input <- reactiveValues(
          
          ## meta data (gct 1.3)
          rdesc=NULL,
          cdesc=NULL
          
        )

        ##################################################
        ## parameters
        global.param <-  reactiveValues(
          
            file.done=F,
            id.done=F,
            grp.done=F,                  ## group assignment finished?
            session.imported=F,          ## flag whether this is an imported session
            session.import.init=F,       ## flag for initiate the session (used for to switch between default values for 'filter.value' and user=sepcified, i.e. after importing a session)
            analysis.run=F,              ## flag whether the analysis has been run
            
            
            file.gct3=F,                 ## special case: GCT v1.3 file
            
            session=NULL,                ## session id
            user=NULL,                   ## user
            grp=NULL,                    ## the actual group assignment
            N.grp=NULL,                  ## number of defined groups
            grp.colors=NULL,             ## group color assignment
            grp.colors.legend=NULL,      ## group colors, names are group names
            
            which.test='One-sample mod T', ## specify test
            
            log.transform='none',        ## log transformation
            norm.data='none',            ## data normalization
            filt.data='none',            ## data filtering

            ##repro.filt='no',           ## reproducibility filter
            repro.filt.val=0.001,
            ##sd.filt='no',              ## sd filter
            sd.filt.val=10,              ## remove lower 10 percent of features with lowest sd

            
            filter.type='adj.p',         ## default filter
            filter.value=0.05,           ## default filter value

            run.test=0,                   ## number of times the 'Run analysis' button has been pressed during a user session

            update.ppi.select= FALSE,     ## trigger selectize, volcano
            update.ppi.select.scat=FALSE, ## trigger selectize, scatterplot
            collapse.ppi=TRUE             ## should PPI query panel be collapsed?

        )

        #####################################################
        ## plotting parameter
        global.plotparam <- reactiveValues(

            ## multiscatter
            ms.max=FALSE,
            ms.min.val=-4,
            ms.max.val=4,

            ## volcano
            volc.ps=1,         ## point size
            volc.ls=1,         ## label size
            volc.grid=T,       ## grid
            volc.maxp=100,     ## max. -log10 p-value
            volc.hyper.fc=1,    ## min. FC for hyperbolic curve
            volc.hyper.curv=3,  ## curvation parameter for hyperbol. curve

            ## heatmap
            hm.cexCol=8,
            hm.cexRow=3,
            hm.scale="none",
            hm.max=FALSE,
            hm.max.val=4,

            ## PCA
            pca.x='PC 1',
            pca.y='PC 2',
            pca.z='PC 3',
            pca.load.topn=20,

            ## correlation matrix
            cm.upper='pearson',
            cm.lower='spearman',
            cm.numb=TRUE
        )

        ## coordinates in volcano plot
        volc <- reactiveValues()
        volc.brush <- reactiveValues()


        ## ###################################################
        ##
        ##               RAM usage indicator
        ##
        ## ###################################################
        ## if(OS != 'Windows'){

        ##     getFreeMem <- function(){
        ##          return( as.numeric(system("awk '/MemFree/ {print $2}' /proc/meminfo",  intern=TRUE)) )
        ##     }
        ##     getUsedMemPerc <- function(){
        ##         free= as.numeric(system("awk '/MemFree/ {print $2}' /proc/meminfo",  intern=TRUE))
        ##         tot= as.numeric(system("awk '/MemTotal/ {print $2}' /proc/meminfo",  intern=TRUE))

        ##         return( (tot-free)/tot)
        ##     }
        ##    ## getCPUutil <- function()
        ##    ##     ##return(as.numeric(system("mpstat -P ALL | awk '/all/ {print $4}'", intern=T)))
        ##    ##     return(as.numeric(system( "mpstat 1 1 | grep '[A|P]M.*all' | awk '{print $4}'", intern=T )))


        ##     ## update every 1 seconds
        ##     RAMused <- reactivePoll(10000, session, getFreeMem, getUsedMemPerc)
        ##     ##CPUused <- reactivePoll(2000, session, getCPUutil, getCPUutil)

        ## }

        ################################################################################
        ##
        ##                                instructions / help pages
        ##
        ################################################################################
        callModule(printHTML, id='getting.started', what='gs', global.input = global.input)
        callModule(printHTML, id='change.log', what='cl', global.input = global.input)
        callModule(printHTML, id='id.column', what='id', global.param=global.param, global.input=global.input)
        callModule(printHTML, id='gct3.file', what='gct3',global.param=global.param, global.input=global.input)
        callModule(printHTML, id='exp.design', what='ed', global.param=global.param, global.input=global.input)
        callModule(printHTML, id='analysis', what='ana', global.param=global.param, global.input=global.input)
        
        #################################
        ## F5 hint
        output$F5hint <- renderText({

            HTML('<p align=\"center\"><font size=\"5\" color=\"red\">To analyze another data set or to start over hit the F5 button.</font></p>' )

        })

        #####################################
        ## Error messages
        output$error <- renderText({
            if( is.null(error$msg) ) return()
            HTML(paste('<p align=\"center\"><font size=\"5\" color=\"red\">', error$msg,'</font></p>'))
        })

        ## #########################################################
        ## toggle all options for export
        ## #########################################################
        observeEvent(input$export.toggle.all, {

            updateCheckboxInput(session, "export.hm", "Heatmap", value=!input$export.toggle.all)
            updateCheckboxInput(session, 'export.box', 'Boxplots', value=!input$export.toggle.all)
            updateCheckboxInput(session, 'export.volc', 'Volcano plot', value=!input$export.toggle.all)
            updateCheckboxInput(session, 'export.phist', 'P-value histogram', value=!input$export.toggle.all)
            updateCheckboxInput(session, 'export.pca', 'PCA', value=!input$export.toggle.all)
            updateCheckboxInput(session, 'export.pca.loadings', 'PCA loadings (xls)', value=!input$export.toggle.all)
            updateCheckboxInput(session, 'export.ms', 'Multiscatter', value=!input$export.toggle.all)
            updateCheckboxInput(session, 'export.excel', 'Excel sheet', value=!input$export.toggle.all)
            updateCheckboxInput(session, 'export.cm', 'Correlation matrix', value=!input$export.toggle.all)
            updateCheckboxInput(session, 'export.profile', 'Profile plot', value=!input$export.toggle.all)

        })

        ## ########################################################
        ##  Memory usage
        ## ########################################################
        ## output$memfree <- renderMenu({

        ##     if(is.null(session$user)) return()
        ##     if(OS == 'Windows') return()

        ##     RAM <- round(RAMused()*100,1)
        ##     status=ifelse(RAM < 70, 'success', 'danger' )

        ##     notificationItem( paste( RAM,'% RAM'), shiny::icon("alert", "fa-1x", lib='glyphicon'), status=status)

        ## })

        ## ########################################################
        ## logged user
        output$logged.user <-renderMenu({

            if(is.null(session$user)) return()

            user <- session$user
            ##user='test'

            notificationItem(user, shiny::icon("user", "fa-1x", lib='glyphicon'), status='success') ##, ic='users', status='info', text='test'),
          ##  notificationItem(user, shiny::icon("users"), status='info') ##, icon='users', status='info', text='test'),


        })

        ## #########################################################
        ## logout user
        output$logout <- renderMenu({
            if(is.null(session$user)) return()
            notificationItem('Logout', icon=shiny::icon("sign-out", "fa-1x"), status='success', href="__logout__")
        })

        ## #########################################################
       ## session name
       ## output$session.name <- renderMenu({
       ##     if(is.null( global.param$label)) return()
       ##     ##session <- global.param$session
       ##
       ##     notificationItem(global.param$label, shiny::icon("star", "fa-1x", lib='glyphicon'), status='success')
       ## })



        ################################################################################
        ##
        ##                      navbar -   render UI
        ##
        ################################################################################
        output$navbar <- renderUI({

            if(!global.param$analysis.run) return()

            ##############################################
            ## determine the number of group comparisons,
            ## e.g. for the number of volcano plots to draw
            ##############################################
            groups.comp <- unique(global.param$grp.comp)
            ##groups <- global.param$grp.comp

            ## ##########################################
            ## class vector
            grp <- global.param$grp
            grp.unique <- unique(grp)

            ################################################################################################################
            ##
            ##                           define the different tabs shown in the navigation bar
            ##
            ################################################################################################################

            ###############################################################
            ## EXPORT tab
            ## - export all figures/tables at once and generate a zip file
            ## - download the zip file
            ###############################################################
            export.tab <- tabPanel('Export',
                                   if(!(global.results$export.results)){
                                      fluidRow(
                                         column(width=6,
                                                box(title="Session name:",
                                                    textInput( 'label', '', value=global.param$label, width=200),
                                                    checkboxInput('export.save.session', 'Save as session', value=F),
                                                    status = "primary",
                                                    solidHeader = T,
                                                    width=NULL
                                                    ),

                                                box(title="Export results", actionButton('export.results', 'Export/Save session'),
                                                    status = "primary",
                                                    solidHeader = T,
                                                    width=NULL)
                                                ),
                                         column(width=6,
                                                box(title="Specify what to export:",
                                                    checkboxInput('export.toggle.all', 'Toggle all', value=F),
                                                    tags$hr(),
                                                    checkboxInput('export.hm', 'Heatmap',value=T),
                                                    checkboxInput('export.box', 'Boxplots',value=T),
                                                    checkboxInput('export.volc', 'Volcano plot',value=T),
                                                    checkboxInput('export.phist', 'P-value histogram',value=T),
                                                    checkboxInput('export.pca', 'PCA',value=T),
                                                    checkboxInput('export.pca.loadings', "PCA loadings (xls)", value = T),
                                                    checkboxInput('export.ms', 'Multiscatter',value=T),
                                                    checkboxInput('export.excel', 'Excel sheet',value=T),
                                                    checkboxInput('export.cm', 'Correlation matrix',value=T),
                                                    checkboxInput('export.profile', 'Profile plot',value=T),

                                                    status = "primary",
                                                    solidHeader = T,
                                                    width=NULL

                                                    )
                                                )
                                       )

                                   } else {
                                       fluidRow(
                                         column(width=6,
                                                box(title="Session name:",
                                                    textInput( 'label', '', value=global.param$label, width=200),
                                                    checkboxInput('export.save.session', 'Save as session', value=F),
                                                    status = "primary",
                                                    solidHeader = T,
                                                    width=NULL
                                                    ),

                                                box(title="Export results", actionButton('export.results', 'Export (.zip)'),
                                                    status = "primary",
                                                    solidHeader = T,
                                                    width=NULL),

                                                box(title='Download results', downloadButton('download.results', 'Download (.zip)'),
                                                    status = "primary",
                                                    solidHeader = T,
                                                    width=NULL)
                                                ),
                                         column(width=6,
                                                box(title="Specify what to export:",
                                                    checkboxInput('export.toggle.all', 'Toggle all', value=F),
                                                    tags$hr(),
                                                    checkboxInput('export.hm', 'Heatmap',value=T),
                                                    checkboxInput('export.box', 'Boxplots',value=T),
                                                    checkboxInput('export.volc', 'Volcano plot',value=T),
                                                    checkboxInput('export.phist', 'P-value histogram',value=T),
                                                    checkboxInput('export.pca', 'PCA',value=T),
                                                    checkboxInput('export.pca.loadings', "PCA loadings (xls)", value = T),
                                                    checkboxInput('export.ms', 'Multiscatter',value=T),
                                                    checkboxInput('export.excel', 'Excel sheet',value=T),
                                                    checkboxInput('export.cm', 'Correlation matrix',value=T),
                                                    checkboxInput('export.profile', 'Profile plot',value=T),

                                                    status = "primary",
                                                    solidHeader = T,
                                                    width=NULL

                                                    )
                                                )
                                       )
                                   } ## end else
                                   )

            ## ##########################################
            ## SUMMARY
            ##    some general numbers on the
            ##    uploaded data
            ##
            ## ##########################################
            summary.tab <-  tabPanel('Summary',


                                     fluidRow(
                                        box(title="Dataset:", solidHeader = T, status = "primary", width = 4,
                                         tableOutput('summary.data')),

                                        box(title="Workflow:", solidHeader = T, status = "primary",width = 4,
                                          tableOutput('summary.workflow')),

                                        box(title="Test results:", solidHeader = T, status = "primary",width = 4,
                                          tableOutput('summary.test'))
                                      ),
                                     fluidRow(
                                      box(title="Quantified features", solidHeader = T, status = "primary",width = 12,
                                          plotlyOutput('summary.nonmissing.data'))
                                    ),
                                    fluidRow(
                                        box(title="Quantified features (cumulative)", solidHeader = T, status = "primary",width = 12,
                                            plotlyOutput('summary.missing.data.row'))
                                    )
                          ) ## end tab panel


            ## ##########################################
            ## SCATTERPLOTS plotly
            ## - for each group comparison
            ##
            ## ##########################################
            scat.tabs <- list()
            scat.tabs[[1]] <- 'Scatterplots'

            for(i in 1:length( grp.unique )){

                ## extract data columns of current experiment
                scat.x <- names(grp)[which(grp == grp.unique[i])]
                scat.y <- names(grp)[which(grp == grp.unique[i])]

                scat.tabs[[i+1]]=tabPanel(paste0( grp.unique[ i ] ),
                                          fluidPage(

                                              ## ###############################
                                              ## select data columns to plot
                                              box( title='Choose expression columns', status = 'primary', solidHeader = T, width=12,

                                                  fluidRow(
                                                      column(1),
                                                      column(5, selectInput( paste0('scat.x.', grp.unique[i]), 'x-axis', scat.x, selected=scat.x[1], multiple=FALSE, selectize=TRUE)),
                                                      column(5, selectInput( paste0('scat.y.', grp.unique[i]), 'y-axis', scat.y, selected=scat.y[2], multiple=FALSE, selectize=TRUE)),
                                                      column(1)
                                                  )

                                                  ), # end box

                                              ## ###########################################################
                                              ## PPI stuff
                                              box( title='Protein-protein interactions', status = 'primary', width=12, collapsible = TRUE, solidHeader=T, collapsed=global.param$collapse.ppi,
                                                  fluidRow(
                                                      column(3,
                                                             selectizeInput( inputId=gsub('\\.','', gsub('\\.', '', paste0('ppi.bait.scat.', grp.unique[i])) ), label=NULL,

                                                                            choices=NULL,
                                                                            selected=NULL,
                                                                            options = list(
                                                                                maxOptions=5,
                                                                                placeholder='Bait protein',
                                                                                onInitialize = I('function() { this.setValue(""); }')
                                                                            )

                                                                            )
                                                             ),
                                                      column(3, checkboxGroupInput(paste('ppi.db.scat', grp.unique[i], sep='.' ), 'Source data', choices=c('BioGRID (human)' = 'bg', 'InWeb' = 'iw', 'Reactome (human)' = 'react'), selected=c('bg', 'iw', 'react') ))#,
                                                      #column(3, checkboxInput( paste('ppi.show.labels.scat', grp.unique[i], sep='.' ), 'Show labels', value=F) )

                                                      ##column(1)
                                                  )),


                                              ## ###############################
                                              ## plot
                                              box( title=grp.unique[i], status = 'primary', solidHeader = T, width=12,
                                                  fluidRow(
                                                       column(12, align='center', plotlyOutput(  paste0('scatterplot.', grp.unique[i]) , width=600, height=600) )
                                                  ) ## end fluiRow

                                                  ) ## end box


                                          ) ## end fluidPage

                                          ) ## end tabPanel
            } ## end for


            ## ##########################################
            ## VOLCANO
            ##      tabs for the volcano plots
            ## NOT for F test
            ############################################
            if( !(global.param$which.test %in% c('mod F', 'none'))) {
              
                volc.tabs <- list()
                volc.tabs[[1]] <- 'Volcanos'
    
                for(i in 1:length(unique(groups.comp))){
                    volc.tabs[[i+1]]=tabPanel(paste0( groups.comp[i] ),
    
                                              fluidPage(
    
                                                  ## ###########################################################
                                                  ## volcano plotting parameters
                                                  box( title='Plotting parameters', status = 'primary', solidHeader = T, width=12,
                                                     fluidRow(
    
    
                                                              column(3, numericInput( paste("cex.volcano",groups.comp[i], sep='.'), "Point size", value=global.plotparam$volc.ps, min=1, step=1, width='100px')),
                                                              ##column(1, numericInput( paste("opac.volcano",groups.comp[i],sep='.'), "Opacity %", value=50, min=0, max=100, step=10)),
                                                              column(3, numericInput( paste("cex.volcano.lab",groups.comp[i],sep='.'), "Label size", value=global.plotparam$volc.ls, min=.1, step=.1, width='100px')),
                                                              column(3, selectInput( paste("grid.volcano",groups.comp[i],sep='.'), "Grid", c(T, F), selected=global.plotparam$volc.grid, width='100px')),
                                                              column(3, numericInput( paste( "max.logP", groups.comp[i], sep='.'), "Max. Log10(P-value)", value=global.plotparam$volc.maxp, min=20, max=300, step=10, width='100px') )
    
                                                              ##column(1, downloadButton(paste('downloadVolcano', groups.comp[i],sep='.'), 'Download (pdf)'))
                                                          )
                                                      ),
    
                                                  ## ###########################################################
                                                  ## PPI stuff
                                                  box( title='Protein-protein interactions', status = 'primary', width=12, collapsible = TRUE, solidHeader=T, collapsed=global.param$collapse.ppi,
                                                      fluidRow(
                                                          column(3,
                                                                  selectizeInput( inputId=gsub('\\.','', gsub('\\.', '', paste0('ppi.bait.', groups.comp[i])) ), label=NULL,
                                                                               ## choices=global.results$id.map$id.concat,
                                                                                choices=NULL,
                                                                                selected=NULL,
                                                                                options = list(
                                                                                    maxOptions=5,
                                                                                    placeholder='Bait protein',
                                                                                    onInitialize = I('function() { this.setValue(""); }')
                                                                                )
    
                                                                                )
                                                                 ),
                                                          column(3, checkboxGroupInput(paste('ppi.db', groups.comp[i], sep='.' ), 'Source data', choices=c('BioGRID (human)' = 'bg', 'InWeb' = 'iw', 'Reactome (human)' = 'react'), selected=c('bg', 'iw', 'react') )),
                                                          column(3, checkboxInput( paste('ppi.show.labels', groups.comp[i], sep='.' ), 'Show labels', value=F) ),
    
                                                          column(1, checkboxInput( paste('ppi.hyper.filt', groups.comp[i], sep='.' ), 'Hyperbolic curve', value=F)),
                                                          column(1, numericInput(paste( "ppi.min.fc", groups.comp[i], sep='.'), 'Min. FC', value=global.plotparam$volc.hyper.fc, min=0, max=100, step=0.1)),
                                                          column(1, numericInput(paste( "ppi.curve", groups.comp[i], sep='.'), 'Curvation', value=global.plotparam$volc.hyper.fc, min=0.1, max=100, step=0.1) )
                                                          ##column(1)
                                                      )),
    
    
                                                  ## ############################################################
                                                  ## the actual plot  plus table
                                                  fluidRow(
                                                      ## plot
                                                      column(width=8,
                                                             box( width=NULL,  title='Volcano plot', status = 'primary', solidHeader = T,
                                                                 plotOutput( paste("volcano", groups.comp[i], sep='.'), height=600, click=paste('plot_click', groups.comp[i], sep='.'), hover=hoverOpts(id=paste('plot_hover', groups.comp[i], sep='.'), delay=10), brush=brushOpts(id=paste('plot_brush', groups.comp[i], sep='.'), resetOnNew=T, delayType='debounce', delay='1000' ), dblclick=paste('plot_dblclick', groups.comp[i], sep='.'))
                                                                 )),
                                                      ## table
                                                      column(width=4,
                                                             box(width=NULL,  title='Selected points', status = 'primary', solidHeader = T,
                                                                 actionButton(inputId=paste('volc.tab.reset', groups.comp[i], sep='.'), label='Reset'),
                                                                 tags$hr(),
                                                                 tableOutput(paste('volc.tab.selected', groups.comp[i], sep='.'))
                                                                 )))
                                              ) ## end fluidPage
    
                                              ) ## end tabPanel
                } ## end for i
            } # end if 

            ################################################
            ##
            ##                 clustering
            ##
            ## #############################################
            
            ############################################
            ## HEATMAP
            ############################################
            hm.tab <-  tabPanel('Static heatmap',
                                      box(title='Static Heatmap', status = 'primary', solidHeader = T, width="100%", height="100%",
                                        fluidRow(
                                        column(2, numericInput( "cexCol", "Font size column", value=ifelse( !is.null(global.plotparam$hm.cexCol), global.plotparam$hm.cexCol, 12  ), min=1, step=1)),
                                        column(2, numericInput( "cexRow", "Font size row", value=ifelse( !is.null(global.plotparam$hm.cexRow), global.plotparam$hm.cexRow, 6), min=1, step=1)),
                                        column(2, selectInput( "hm.scale", "Scale", c("row","column","none"), selected=global.plotparam$hm.scale)),
                                        column(2, selectInput( "hm.clust", "Cluster", c("column","row","both","none"), selected=ifelse(global.param$which.test != "mod F", "none" ,"both"))),
                                        column(2, checkboxInput('hm.max', 'Cap values', value=global.plotparam$hm.max)),
                                        column(2, numericInput( "hm.max.val", "Max. value", value=global.plotparam$hm.max.val, step=1, min=2))
                                        ),
                                        fluidRow(
                                            column(12, align='center', plotOutput("HM", height=min( dynamicHeightHM( nrow(global.results$filtered)), 1200 ), width=dynamicWidthHM(length(global.param$grp))) )
                                        )
                                      )
                                    
            )
            ## ###########################################
            ##   HEATMAPPLY
            ## ###########################################
            hm.int.tab <-  tabPanel('Interactive heatmap',
                                box(title='Interactive Heatmap', status = 'primary', solidHeader = T, width="100%", height="100%",
                                    fluidRow(
                                      column(2, selectInput( "hm.int.scale", "Scale", c("row","column","none"), selected=global.plotparam$hm.scale)),
                                      column(2, selectInput( "hm.int.clust", "Cluster", c("column","row","both","none"), selected=ifelse(global.param$which.test != "mod F", "none" ,"both"))),
                                      column(2, checkboxInput('hm.int.max', 'Cap values', value=global.plotparam$hm.max)),
                                      column(2, numericInput( "hm.int.max.val", "Max. value", value=global.plotparam$hm.max.val, step=1, min=2))
                                    ),
                                    fluidRow(
                                      column(12, align='center', plotlyOutput("HM.int", height=1000, width=1000 ))
                                    )
                                )
            )
            
            ## #############################################
            ##   FANPLOT
            ## #############################################
            hc.fanplot <-  tabPanel('Fanplot',
                                    box(title='Fanplot', status = 'primary', solidHeader = T, width="100%", height="100%",
                                        fluidRow(
                                          column(12, align='center', plotOutput("HC.fan" ))
                                        )
                                    )
            )
            
            clust.tab <- vector('list', 3)
            clust.tab[[1]] <- hm.tab
            clust.tab[[2]] <- hm.int.tab
            clust.tab[[3]] <- hc.fanplot
            
            ## ##############################################
            ##  MORPHEUS widget
            ##
            ## ##############################################
            morph.tab <- tabPanel('Morpheus',
                                    fluidRow(
                                      column(12, morpheusOutput("HM.morpheus",height=min( dynamicHeightHM( nrow(global.results$filtered)), 1200 ), width=dynamicWidthHM(length(global.param$grp))))
                                    )
            )
            
            

            #############################################
            ##
            ##               PCA
            ##
            #############################################
            pca.tab2 <- vector('list', 3)
            ## ##################################################
            ## Run PCA
            ## ##################################################
            pca.tab2[[1]] <- tabPanel('Run me first',
             fluidRow(

                      box(title='Summary', status='primary', solidHeader = T, htmlOutput('run.pca')),

                      box(title='Variance', status='primary', solidHeader = T, plotOutput('pca.var'))
                    )
            )
            ## ##################################################
            ## PC plots
            ## ##################################################
            pca.tab2[[2]] <- tabPanel('PC plots',
                                      fluidPage(

                                          fluidRow(
                                            box( title = 'Select principle components', status='primary', solidHeader = T, align='center',
                                              column(width=4, selectInput('pca.x', 'x-axis', paste('PC', 1:10)), selected=global.plotparam$pca.x),
                                              column(width=4, selectInput('pca.y', 'y-axis', paste('PC', 1:10), selected=global.plotparam$pca.y)),
                                              column(width=4, selectInput('pca.z', 'z-axis', paste('PC', 1:10), selected=global.plotparam$pca.z))
                                            )
                                          ),
                                          fluidRow(
                                            box( title='2D', status = 'primary', solidHeader = T, width = 1000, height = 700,
                                                 column(12, align='center', plotlyOutput("pcaxy.plotly", width=600, height=600))
                                            )
                                          ),

                                          fluidRow(
                                            box( title='3D', status = 'primary', solidHeader = T, width = 800, height = 900,
                                                 column(12, align='center', plotlyOutput("pcaxyz.plotly", width=800, height=800))
                                            )
                                          )##,


                                          ##fluidRow(
                                          ##  box(title='Figure to download', status = 'primary', solidHeader = T, width = 1200, height = 400,
                                          ##      column(12, align='center', plotOutput("pca", width=900, height=300))
                                          ##  )
                                          ##)


                                      )
                                      )
            ## ####################################################################
            ## Loadings plot
            ## ####################################################################
            pca.tab2[[3]] <- tabPanel('Loadings',
                                      fluidPage(
                                        fluidRow(
                                            column(width=12,
                                            box(title='Loadings by Ozan Aygun', solidHeader=T, status='primary', width=1000,## height=min( nrow(global.results$filtered), global.plotparam$pca.load.topn )*20+50,
                                                sliderInput("pca.load.topn", "Choose number of loadings", 1, 100, 20),
                      ##                          plotOutput("pca.loadings")##, width=1000, height=min(nrow(global.results$filtered), global.plotparam$pca.load.topn )*20 )
                      ##                        ),
                      ##                      box(title = "PCA loadings scatterplots",solidHeader = T, status = 'primary',width = 1000,
                                                background = "navy",
                                                plotOutput("scatter.pca.loadings")
                                            ))

                                        )
                                      )
            )


            ## ###########################################
            ## TABLE: filtered result table
            ##
            ## ###########################################
            table.tab <- tabPanel('Table',
                     fluidPage(
                         fluidRow(column(12, tags$h3('Result table (filtered):'))),
                         fluidRow(column(12, tags$br())),
                         fluidRow(column(12, dataTableOutput("tableprev")))
                     )
                     )
            #############################################
            ## QC tabs
            ##
            #############################################
            qc.tabs <- vector('list', 5)
            ##names(qc.tabs) <- c('Boxplots', 'P-values', 'Multi scatter', 'Correlation matrix', 'Correlation matrix transposed')

            ###########################
            ## boxplots
            ##qc.tabs[['Boxplots']] <- tabPanel('Boxplots',
            qc.tabs[[1]] <- tabPanel('Boxplots',
                                              fluidPage(
                                                  fluidRow(

                                                      box(title="Before normalization", solidHeader=T, status="primary",
                                                          column(width=12, plotOutput("expr.boxplot", width="100%", height=max( 30*(ncol( global.input$table)+2), 500)))
                                                          ),
                                                      if(!is.null(global.results$table.norm)){
                                                          box(title="After normalization", solidHeader=T, status="primary",
                                                          column(width=12, plotOutput("expr.boxplot.norm", width="100%", height=max( 30*(ncol( global.input$table)+2), 500)))
                                                          )
                                                      }

                                                  )

                                              )
                                     )

            ###########################
            ## profile plots
            ##qc.tabs[['Profile plots']] <- tabPanel('Profile plots',
            qc.tabs[[2]] <- tabPanel('Profile plots',

                                                       if(is.null(global.results$table.norm)){
                                                           fluidPage(
                                                               fluidRow(
                                                                   box(title="Before normalization", solidHeader=T, status="primary",
                                                                       column(width=12, plotOutput("expr.profile"))
                                                               ))
                                                            )
                                                       } else {
                                                           fluidPage(

                                                               fluidRow(
                                                                   box(title="Before normalization", solidHeader=T, status="primary",
                                                                      column(width=12, plotOutput("expr.profile"))
                                                                      ),
                                                                   box(title="After normalization", solidHeader=T, status="primary",
                                                                       column(width=12,plotOutput("expr.profile.norm"))
                                                                       )
                                                               )
                                                           )
                                                       }
                                                   )

            ###########################
            ## P-value distribution
            ##qc.tabs[['P-values']] <- tabPanel('P-values',
            qc.tabs[[3]] <- tabPanel('P-values',
                                              fluidPage(
                                                  fluidRow(
                                                      box(title="Distribution of P-values", solidHeader=T, status="primary", width=1000, height=600*ifelse( global.param$which.test != 'mod F', length(unique(global.param$grp.comp)), 1 ),
                                                          column(12, plotOutput("pval.hist"))
                                                          )
                                                  )
                                              )
                                              )
            ############################
            ## correlation multiscatter
            ##qc.tabs[['Multi scatter']] <- tabPanel('Multi scatter',
            qc.tabs[[4]] <- tabPanel('Multi scatter',
                                                   fluidPage(
                                                       fluidRow(
                                                           box(title="Parameters", solidHeader=T, status="primary",
                                                               column(2, checkboxInput('ms.max', 'Define limits', value=FALSE)),
                                                               column(2, numericInput( "ms.min.val", "min.", value=-4, step=1)),
                                                               column(2, numericInput( "ms.max.val", "max.", value=4, step=1)),
                                                               column(6)
                                                               )
                                                       ),
                                                       fluidRow(
                                                           box(title='Multiscatter', solidHeader=T, status="primary", width=100*(ncol(data.frame(global.input$table))-1), height=130*(ncol(data.frame(global.input$table)) - 1),
                                                               column(12, plotOutput("multi.scatter"))
                                                               )
                                                       )
                                                   )

                                     )
             ###########################
             ## correlation matrix
             ## qc.tabs[['Correlation matrix']] <- tabPanel('Correlation matrix',
             qc.tabs[[5]] <- tabPanel('Correlation matrix',

                                                         fluidPage(
                                                             fluidRow(
                                                                 box(title='Parameters', solidHeader=T, status="primary",
                                                                     column(3,  selectInput( "cm.upper", "Upper triangle", c("pearson","spearman","kendall"), selected=global.plotparam$cm.upper)),
                                                                     column(3,  selectInput( "cm.lower", "Lower triangle", c("pearson","spearman","kendall"), selected=global.plotparam$cm.lower)),
                                                                     column(1,  checkboxInput('cm.numb', 'Show numbers', value=global.plotparam$cm.numb)),
                                                                     column(5))
                                                             ),
                                                             fluidRow(
                                                                 box(title='Correlation matrix', solidHeader=T, status="primary", width=dynamicWidthHM( length(global.param$grp) ), height=dynamicWidthHM( length(global.param$grp) ),
                                                                     column(12, plotOutput("correlation.matrix"))
                                                                 )
                                                             )
                                                         )
                                                         )
             ###########################
             ## correlation matrix
             ##qc.tabs[['Correlation matrix transposed']] <- tabPanel('Correlation matrix transposed',
             qc.tabs[[7]] <- tabPanel('Correlation matrix transposed',

                                                         fluidPage(
                                                             fluidRow(
                                                                 column(1,  selectInput( "cm.upper", "Upper triangle", c("pearson","spearman","kendall"), selected=global.plotparam$cm.upper)),
                                                                 column(1,  selectInput( "cm.lower", "Lower triangle", c("pearson","spearman","kendall"), selected=global.plotparam$cm.lower)),
                                                                 column(1,  checkboxInput('cm.numb', 'Show numbers', value=global.plotparam$cm.numb)),
                                                                 column(9)
                                                             ),
                                                             fluidRow(
                                                                 column(11),
                                                                 column(1, downloadButton('downloadCMtrans', 'Download (pdf)'))
                                                             )
                                                         ),
                                                         tags$br(),tags$hr(),tags$br(),
                                                         fluidPage(
                                                             fluidRow(
                                                                 column(12, plotOutput("correlation.matrix.trans"))
                                                             )
                                                         )
                                      )



            ## #################################################################################
            ##
            ##                         insert the tabs
            ##
            ## #################################################################################
            if( !(global.param$which.test %in% c('mod F', 'none'))){
            
              navbarPage(title='', id='mainPage',


                        #######################################
                        ##            insert summary tab
                        #######################################
                        summary.tab,
                        #######################################
                        ##              insert heatmap
                        #######################################
                        ##hm.tab,
                        navbarMenu('Clustering', clust.tab[[1]], clust.tab[[2]], clust.tab[[3]]),
                        
                        ##do.call(navbarMenu, hm.tab),
                        #morph.tab,
                        
                        #######################################
                        ##              insert volcanos
                        #######################################
                        #if( !(global.param$which.test %in% c('mod F', 'none'))){
                           do.call(navbarMenu, volc.tabs),
                      #  } else {
                      #       list()
                      #     },
                       
                       ## ####################################
                       ##          insert scatterplots
                       ## ####################################
                       do.call(navbarMenu, scat.tabs),

                       #######################################
                       ##              insert PCA
                       #######################################
                       navbarMenu('PCA', pca.tab2[[1]], pca.tab2[[2]], pca.tab2[[3]]),

                       #######################################
                       ##           insert table preview
                       #######################################
                       table.tab,

                       #######################################
                       ## QC
                       if(global.param$which.test != 'none') {
                           navbarMenu('QC', qc.tabs[[1]], qc.tabs[[2]], qc.tabs[[3]], qc.tabs[[4]], qc.tabs[[5]])
                       } else {
                           navbarMenu('QC', qc.tabs[[1]], qc.tabs[[2]], qc.tabs[[4]], qc.tabs[[5]])
                       },

                       ## #######################################
                       ## export
                       export.tab
                       ) ## end navbarpage
            } else {
              
              navbarPage(title='', id='mainPage',
                      #######################################
                      ##            insert summary tab
                      #######################################
                      summary.tab,
                      #######################################
                      ##              insert heatmap
                      #######################################
                      ##hm.tab,
                      navbarMenu('Clustering', clust.tab[[1]], clust.tab[[2]], clust.tab[[3]]),
                      ##do.call(navbarMenu, hm.tab),
                      #morph.tab,
                      
                      ## ####################################
                      ##          insert scatterplots
                      ## ####################################
                      do.call(navbarMenu, scat.tabs),
                      
                      #######################################
                      ##              insert PCA
                      #######################################
                      navbarMenu('PCA', pca.tab2[[1]], pca.tab2[[2]], pca.tab2[[3]]),
                      
                      #######################################
                      ##           insert table preview
                      #######################################
                      table.tab,
                      
                      #######################################
                      ## QC
                      if(global.param$which.test != 'none') {
                        navbarMenu('QC', qc.tabs[[1]], qc.tabs[[2]], qc.tabs[[3]], qc.tabs[[4]], qc.tabs[[5]])
                      } else {
                        navbarMenu('QC', qc.tabs[[1]], qc.tabs[[2]], qc.tabs[[4]], qc.tabs[[5]])
                      },
                      
                      ## #######################################
                      ## export
                      export.tab
                ) ## end navbarpage
            }
            

            ## global.param$ins.volc.n <- global.param$ins.volc.n + 1
        }) ## end renderUI


        #########################################################################################
        #  
        #                                  user input
        #
        #########################################################################################

        # #########################################################################
        #
        # 1) UI file upload
        # - entry point for the app
        # - session id is generated here
        # - search path to browse saved sessions
        #
        output$file.upload <- renderUI({

            #if(!is.null(input$file)) return()
            #if(!is.null( global.input$file)) return()
            if(global.param$file.done) return()
          
          
            ##View(ppi$iw[1:10, ])

            ## ######################################
            ## generate session ID and prepare data
            ## directory
            #########################################

            ## generate 'session id'
            if(is.null(global.param$session))
                global.param$session <- paste(paste(letters[sample(26, 5)], collapse=''), paste(sample(100,5), collapse=''), sep='')

            ## ##############################################
            ## authenticated session
            ##if(is.null(global.param$user)){
            if(!is.null(session$user)){

                ## ########################################
                ## user name
                global.param$user <- sub('@.*', '', session$user) ## to avoid problems...

                ## ########################################
                ## default search path
                search.path <- c()
                search.path[1] <- paste(DATADIR, global.param$user, sep='')

                ## #########################################
                ## parse 'user_roles.txt'
                ##
                ## - determine the folder on the server
                ##   the current user has access to
                ##
                ## #########################################
                user.roles <- read.delim(paste(APPDIR, '/conf/user-roles.txt', sep=''), stringsAsFactors=F)

                ## check whether the user appears as collaborator in a project
                idx <- grep( paste('(^|;)', session$user, '($|;)', sep=''), user.roles$collaborator)

                ## if so add the project path to the search path
                if(length(idx) > 0){
                    ##search.path <- c()
                    for(i in 1:length(idx)){

                        ## folder of the project OWNER to be parsed
                        dir.owner <- paste(DATADIR, sub('@.*', '', user.roles$owner[idx[i]]), sep='')

                        ## check if the folder exists (if not, 'user-roles.txt' has not been updated)
                        if(dir.exists(dir.owner)){
                            tmp <- grep( paste(user.roles$project[ idx[i] ], '_session.*RData$', sep='' ),
                                        dir( dir.owner, full.names=T, recursive=T), value=T)
                            ##search.path[i+1] <- sub('^(.*/).*' , '\\1', tmp)
                            search.path <- c( search.path, sub('^(.*/).*' , '\\1', tmp) )
                        }
                    }
                }
                ## store the search path
                global.param$search.path <- search.path
            }

            ##########################################
            ## upload form
            list(
              br(),
              HTML('<font size=\"3\"><b>Upload file (txt, csv, gct):</b></font>'),
              fileInput("file", "", accept=c('text/csv',
                       'text/comma-separated-values,text/plain',
                       '.csv', '.txt', '.tsv', '.gct')),
              HTML('<hr border-width:\"10px\">')
            )
        })

        
        ##@##################################
        ## 1b) UI browse saved sessions
        ## only available on the server
        output$browse.sessions <- renderUI({

            if(!is.null( global.input$file)) return()
            if(is.null(global.param$user)) return()
            if(length(global.param$user)==0) return()

            ## #########################################
            ## identify all sessions the user has
            ## access to
            search.path <-  global.param$search.path

            saved.sessions <- list()
            for(i in 1:length(search.path))
                saved.sessions[[i]] <- grep( '_session.*RData$', dir( search.path[i], full.names=T, recursive=T ), value=T )
            saved.sessions <- unlist(saved.sessions)


            ## subfolders in user directory
            ##saved.sessions <- grep( '_session.RData', dir( global.param$search.path, full.names=T, recursive=T ), value=T )

            ## don't show the panel if there is no saved session
            if(length(saved.sessions) == 0) return()

            ## get the time stamp of the files
            time.tmp <- file.info(saved.sessions)$ctime
            names(saved.sessions) <-  paste( sub('_.*','', sub('.*/','',saved.sessions)), time.tmp, sep='_' )

            ## order by time
            saved.sessions <- saved.sessions[ order(time.tmp) ]

            ## store saved sessions
            global.param$saved.sessions <- saved.sessions

            list(
                ##selectInput('session.browse', paste('My saved sessions',sep=''), choices=names(saved.sessions)),
                selectizeInput(inputId = 'session.browse', 
                               label = 'Saved sessions:',
                               choices=names(saved.sessions), 
                               options=list( 
                                 placeholder = 'Search sessions',
                                 onInitialize = I('function() { this.setValue(""); }')
                               )
                               ),
              
                ##selectInput('session.browse', paste('Saved sessions (', sub('_at_','@',global.param$user),')',sep=''), choices=sort(names(saved.sessions))),
                actionButton('session.browse.import', 'Import'),
                actionButton('session.manage', 'Manage sessions', onclick =paste("window.open('", CONFAPP,"', 'newwindow', 'width=500 height=600'); return false;", sep=''))
            )
        })

        
        # ###################################################
        #       group assignment for gct 1.3 files
        # - id column defined as first column in gct
        # - group assingment from column annotations
        #
        output$define.groups.gct3 <- renderUI({
          
          if(!global.param$file.gct3) return()
          if(global.param$grp.done) return()
          

          list(
            ## experimental design
            #HTML('<font size=\"3\"><b>Found GCT v1.3 file with', ncol(global.input$cdesc),'annotation columns. Choose one column as class vector for marker selection.</b></font><br>'),
            radioButtons("grp.gct3", "Choose column", colnames(global.input$cdesc)),
            actionButton("update.grp.gct3", 'OK')
          )
          
        })
        # preview levels current selection
        output$grp.gct3.prev <- renderText({
          if(is.null(input$grp.gct3)) return()
          if(global.param$grp.done) return()
          
          HTML(paste('<br><p><font size=\"4\"><b>Current selection:</b><b>', input$grp.gct3,'</b></p><br>'))
          
        })
        
        # preview levels current selection
        output$grp.gct3.prev.tab <- renderTable({
          if(is.null(input$grp.gct3)) return()
          if(global.param$grp.done) return()
          
          tab <- table(global.input$cdesc[, input$grp.gct3]) 
          tab <- data.frame(Level=names(tab), Freq=unlist(tab))
          
          list(
            tab
          )
        })
        
        # #####################################
        # define groups for GCTv3
        observeEvent(input$update.grp.gct3, {
          
          tab <- global.input$table
          cdesc <- global.input$cdesc
          
          #cat('\n\n', cdesc[, input$grp.gct3],'\n\n')
          
          # store grp coloum
          global.param$grp.gct3 <- input$grp.gct3
          
          # initialize grp file
          Column.Name <- colnames(tab)
          Experiment <- rep('', length(Column.Name))
          names(Experiment) <- Column.Name
          Experiment[ rownames(cdesc) ] <- cdesc[, input$grp.gct3]
          
          #grp.gct3 <- input$grp.gct3
          #save(tab, cdesc, grp.gct3, file='tmp.RData')
          
          
          global.param$cdesc.all <- global.param$cdesc.selection <- setdiff(colnames(cdesc),  input$grp.gct3)
          
          # data frame
          grp.file=data.frame(
            Column.Name,
            Experiment,
            stringsAsFactors = F
              )
          
          ## ################################
          ## ANNOTATION: extract empty cells
          ## - corresponding columns will be carried over as
          ##   annotation columns in the result file
          grp.anno <- grp.file[which(nchar( Experiment) == 0 ), ]
          grp.anno <- setdiff( grp.anno$Column.Name, global.param$id.col.value )
          

          if(length(grp.anno)>0)
            global.input$table.anno <- data.frame(id=global.results$id.map[, 'id'], global.input$table[ , grp.anno])
          
          ## ################################
          ## EXPRESSION
          ## - extract all non-empty cells in the 'Experiment' column
          exprs.idx <- rownames(cdesc)
          grp.exprs <- grp.file[exprs.idx, ]
          
          ## order alphabetically to make coloring consistent
          grp.exprs <- grp.exprs[order(grp.exprs$Experiment), ]
          
          ## class vector
          grp=grp.exprs$Experiment
          names(grp)=grp.exprs$Column.Name
          
          ## update input table, keep id and expression columns
          global.input$table <- global.input$table[ , c(global.param$id.col.value, names(grp))]
          
          ################################
          ## update number of groups
          global.param$N.grp <- length(unique( na.omit(grp)) )
          
          ## store group assignment
          global.param$grp <- global.param$grp.all <- grp
          global.param$grp.comp.all <- global.param$grp.comp <- unique(grp)
          
          ## group colors
          grp.col <- rep(GRPCOLORS[1], length(grp))
          names(grp.col) <- names(grp)
          
          for(i in 2:length(unique(grp))) grp.col[ which(grp == unique(grp)[i]) ] <- GRPCOLORS[i]
          global.param$grp.colors <- grp.col
          
          ## group colors for figure legend
          idx <- !duplicated(grp)
          grp.col.legend = grp.col[idx]
          names(grp.col.legend) <- grp[idx]
          
          global.param$grp.colors.legend.all <- global.param$grp.colors.legend <- grp.col.legend
          
          # ####################################################  
          # colors for other annotation tracks
          col.tmp <- cdesc.colors(cdesc, global.param$grp.gct3, grp.col.legend)
          global.param$anno.col.all <- global.param$anno.col <- col.tmp$anno.col
          global.param$anno.col.color.all <- global.param$anno.col.color <- col.tmp$anno.col.color
          

          ## all done
          global.param$grp.done = T
          
        })
        
        
        ## ##################################################
        ##
        ##           2) UI pick id column
        ##
        output$choose.id.column <- renderUI({

            if(global.param$analysis.run) return()
            if(!global.param$file.done) return()
            if(global.param$id.done) return()
            
            ## get uploaded table and column names
            tab <- global.input$table
            tab.colnames <- global.input$table.colnames

            ## try to find 'id' and move it to the first position
            id.idx <- grep('id|ID', tab.colnames)
            if(length(id.idx) > 0){
                tab.colnames <- c(tab.colnames[id.idx], tab.colnames[-id.idx])
            }

            ## radio button to pick id column
            list(
                ## experimental design
                HTML('<font size=\"3\"><b>Export experimental design file:</b></font><br>'),
                downloadButton("exportTemplate", 'Export'),
                HTML('<br><hr size=\"5\">'),
                radioButtons( "id.col.value", "Choose ID column", tab.colnames),
                actionButton("id.col", 'OK')
            )

        })

        ######################################
        ## 3) UI upload experimental design file
        output$define.groups <- renderUI({

            if(is.null(input$id.col)) return() ## no id colum selected
            if( !is.null(input$id.col))
                if( input$id.col == 0 ) return() ## not pressed yet
            if(!is.null(global.results$data)) return() ## test has been run
            if(!is.null(global.param$grp)){            ## group assignment has been RUN
                if(sum(is.na(global.param$grp)) == 0) return() ## group assignemnt has been DONE
            }

            list(
                ## upload template
                fileInput("exp.file", "Upload experimental design file", accept=c('text/plain','.txt')),
                actionButton( 'update.grp', 'Next')
            )
        })

        ######################################
        ## UI filter type
        output$filter.type <- renderUI({
            ## show after analysis had been run
            if(!global.param$analysis.run) return()

            ## if a test has been performed
            if(global.param$which.test != 'none')
                 res <- list( selectInput('filter.type', 'Filter based on:', c('nom.p', 'adj.p', 'top.n', 'none'), selected=global.param$filter.type) )
             else
                res <- list(selectInput('filter.type', 'Filter based on:', c('none'), selected='none'))
            res
        })

        ######################################
        ## UI conditional filter value
        ## - depends on 'input$filter.type'
        output$filter.value <- renderUI({

            ## show after analysis had been run
            if(!global.param$analysis.run) return()

            ## choose different default values depending on
            ## chosen filter type
            res <- list(
                  conditionalPanel(condition = "input['filter.type'] == 'top.n'", numericInput( "filter.value.top.n", "Top N features", value=50, min=2, step=1)),
                  conditionalPanel(condition = "input['filter.type'] == 'nom.p'", numericInput( "filter.value.nom.p", "P-value filter", value=0.01, min=0, max=1, step=1e-2)),
                  conditionalPanel(condition = "input['filter.type'] == 'adj.p'", numericInput( "filter.value.adj.p", "Corrected P-Value (FDR)", value=0.05, min=0, max=1, step=1e-2))
            )

            res
        })

        
        # ##########################################################
        #
        #     Modal window: group selection
        #
        # ###########################################################
        observeEvent(input$select.groups.button, {
          
        ##  if(global.param$run.test > 0){
              showModal(modalDialog(
                size='m',
                title = "Modify selection",
                footer = fluidRow(
                  column(6),
                  column(3, actionButton(inputId = 'update.groups.button.modal' , label='Update')),
                  column(3, modalButton(label='Close'))
                  ),
                #footer = actionButton(inputId = '.groups.button.modal' , label='OK'),
                fluidPage(
                  fluidRow(
                    column(6, HTML('<b>Select groups</b>')),
                    column(6, HTML('<b>Select annotation tracks</b>'))
                    
                  ),
                  fluidRow(
                    column( 6, checkboxGroupInput('select.groups', label=' ', choices = unique(global.param$grp.comp.all), selected = unique(global.param$grp.comp.selection))),
                    column( 6, checkboxGroupInput('select.anno', label=' ', choices = unique(global.param$cdesc.all), selected = unique(global.param$cdesc.selection)))
                  )
                ),
                easyClose = FALSE
              ))
          # } else {
          #   showModal(modalDialog(
          #     size='m',
          #     title = "Modify selection",
          #     footer = fluidRow(
          #       column(6),
          #       column(3, actionButton(inputId = 'update.groups.button.modal' , label='Update')),
          #       column(3, modalButton(label='Close'))
          #     ),
          #     #footer = actionButton(inputId = '.groups.button.modal' , label='OK'),
          #     fluidPage(
          #       fluidRow(
          #         column(6, HTML('<b>Select groups</b>')),
          #         column(6, HTML('<b>Select annotation tracks</b>'))
          #         
          #       ),
          #       fluidRow(
          #         column( 6, checkboxGroupInput('select.groups', label=' ', choices = unique(global.param$grp.comp.all), selected = unique(global.param$grp.comp.selection))),
          #         column( 6)
          #       )
          #     ),
          #     easyClose = FALSE
          #   ))
          #   
          # }
          # 
        })
        # ##################################################
        # update based on selection
        observeEvent(input$update.groups.button.modal,{  
          
          ## ####################
          ## update class vector
          
          # super important! to memorize previous selections in the modal window
          global.param$grp.comp.selection <- input$select.groups
          grp.selection <- global.param$grp.all
          
          ## extract groups selected in the modal window
          grp.unique <- unique( unlist( strsplit( sub('\\.vs\\.', ' ', input$select.groups), ' ')))
          
          ## update selection
          grp.selection <- grp.selection[ grep(paste('^', paste(grp.unique, collapse='|'), '$', sep=''), grp.selection) ]
          global.param$grp.selection <- grp.selection
          
          ## ANNOTATION TRACKS: GCT v1.3
          if(!is.null(global.param$anno.col)){
            
            cdesc.selection <- input$select.anno
            global.param$cdesc.selection <- cdesc.selection
            
            # update data tracks and colors
            anno.col <- global.param$anno.col.all
            anno.col.color <- global.param$anno.col.color.all
            
            # preserve last column (class vector)
            global.param$anno.col <- anno.col[names(global.param$grp.selection), c( cdesc.selection, global.param$grp.gct3) ]
            global.param$anno.col.color <- anno.col.color[ c(cdesc.selection, global.param$grp.gct3 ) ]
          }
          
          #save(cdesc.selection, file='tmp.RData')
          ## update class vector
          #cdesc.selection <- global.param$cdesc.all
          #grp.unique <- unique( unlist( strsplit( sub('\\.vs\\.', ' ', input$select.groups), ' ')))
          #grp.selection <- grp.selection[ grep(paste('^', paste(grp.unique, collapse='|'), '$', sep=''), grp.selection) ]
          #global.param$cdesc.selection <- cdesc.selection
          
    
        })
        
        # ##################################################
        # update ANNOTATION TRACKS based on selection
        #observeEvent(input$select.anno, {
        #  global.param$cdesc.selection <- input$select.anno
        #  
        #  ## update class vector
        #  cdesc.selection <- global.param$cdesc.all
        #  #grp.unique <- unique( unlist( strsplit( sub('\\.vs\\.', ' ', input$select.groups), ' ')))
        #  #grp.selection <- grp.selection[ grep(paste('^', paste(grp.unique, collapse='|'), '$', sep=''), grp.selection) ]
        #  global.param$cdesc.selection <- cdesc.selection
          
        #})
        
        
        ## #####################################################################
        ## UI: set up analysis
        ##
        ##
        ## #####################################################################
        output$list.groups <- renderUI({

            if( !global.param$grp.done ) return()

            ####################################
            ## initialize
            if( is.null(input$filt.data)){

                list(
                    
                     radioButtons('log.transform', 'Log-transformation', choices=c('none', 'log10', 'log2'), selected=global.param$log.transform),
                     radioButtons('norm.data', 'Data normalization', choices=c('Median', 'Median-MAD', '2-component', 'Quantile', 'none'), selected=global.param$norm.data),

                     radioButtons('filt.data', 'Filter data', choices=c('Reproducibility', 'StdDev', 'none'), selected=global.param$filt.data),

                     ##radioButtons('sd.filt', 'SD filter (NOT WORKING)', choices=c('yes', 'no'), selected=global.param$sd.filt),
                     radioButtons('which.test', 'Select test', choices=c('One-sample mod T', 'Two-sample mod T', 'mod F', 'none'), selected=global.param$which.test),

                     actionButton('run.test', 'Run analysis!'),
                     br(),
                     hr(),
                     actionButton('select.groups.button', 'Select Groups')
                     #fluidRow(
                    #  column(6, actionButton('select.groups.button', 'Select Groups')),
                    #  column(6, actionButton('select.anno.button', 'Select Annotation Tracks'))
                    # )
                )
            }
            ## ###################################################
            ## no filter
            else if(input$filt.data == 'none'){

                list(
                  
                     radioButtons('log.transform', 'Log-transformation', choices=c('none', 'log10', 'log2'), selected=global.param$log.transform),
                     radioButtons('norm.data', 'Data normalization', choices=c('Median', 'Median-MAD', '2-component', 'Quantile', 'none'), selected=global.param$norm.data),

                     radioButtons('filt.data', 'Filter data', choices=c('Reproducibility', 'StdDev', 'none'), selected='none'),

                     radioButtons('which.test', 'Select test', choices=c('One-sample mod T', 'Two-sample mod T', 'mod F', 'none'), selected=global.param$which.test),

                     actionButton('run.test', 'Run analysis!'),
                     br(),
                     hr(),
                    # fluidRow(
                    #   column(6, actionButton('select.groups.button', 'Select Groups')),
                    #   column(6, actionButton('select.anno.button', 'Select Annotation Tracks'))
                     #)
                     actionButton('select.groups.button', 'Select Groups')
                )
            }

            ## ###################################################
            ## Reproducibility filter
            else if(input$filt.data == 'Reproducibility'){

                list(
                  
                     radioButtons('log.transform', 'Log-transformation', choices=c('none', 'log10', 'log2'), selected=global.param$log.transform),
                     radioButtons('norm.data', 'Data normalization', choices=c('Median', 'Median-MAD', '2-component', 'Quantile', 'none'), selected=global.param$norm.data),

                     radioButtons('filt.data', 'Filter data', choices=c('Reproducibility', 'StdDev', 'none'), selected='Reproducibility'),
                     selectInput('repro.filt.val', 'alpha', choices=c(.1, .05, 0.01, 0.001 ), selected=global.param$repro.filt.val),

                     radioButtons('which.test', 'Select test', choices=c('One-sample mod T', 'none'), selected='One-sample mod T'),

                     actionButton('run.test', 'Run analysis!'),
                     br(),
                     hr(),
                     actionButton('select.groups.button', 'Select Groups')
                )
            }

            ## ###################################################
            ## StdDev filter
            else if(input$filt.data == 'StdDev'){

                list(
                  
                    radioButtons('log.transform', 'Log-transformation', choices=c('none', 'log10', 'log2'), selected=global.param$log.transform),
                    radioButtons('norm.data', 'Data normalization', choices=c('Median', 'Median-MAD', '2-component', 'Quantile', 'none'), selected=global.param$norm.data),

                    radioButtons('filt.data', 'Filter data', choices=c('Reproducibility', 'StdDev', 'none'), selected='StdDev'),
                    sliderInput('sd.filt.val', 'Percentile StdDev', min=10, max=90, value=global.param$sd.filt.val),

                    radioButtons('which.test', 'Select test', choices=c('One-sample mod T', 'Two-sample mod T', 'mod F', 'none'), selected=global.param$which.test),
                    actionButton('run.test', 'Run analysis!'),
                    br(),
                    hr(),
                    actionButton('select.groups.button', 'Select Groups')
                )
            }

        })


        #############################################################################
        ## observer
        ##                  export all analysis results
        ##
        ## - generate a zip-file
        #############################################################################
        observeEvent(input$export.results, {

            ##cat('User:', session$user, '\n')
            if(!is.null(error$msg)) return()

            ## update label
            global.param$label <- gsub('_| |,|;|\\:|\\+|\\*', '-', input$label)


            ###############################
            ## apply filter
            res = global.results$filtered

            #######################################
            ## extract expression values
            res = res[, names(global.param$grp)]

            ## groups to compare
            grp.comp <- unique( global.param$grp.comp )


            #############################################################
            ##                correlation matrix
            #############################################################
            if(input$export.cm){
                withProgress(message='Exporting', detail='correlation matrix',{
                fn.cm <- paste(global.param$session.dir, '/correlation_matrix_lo_',global.plotparam$cm.lower, '_up_',global.plotparam$cm.upper, '.pdf', sep='')
                ##                plotCorrMat(lower=input$cm.lower, upper=input$cm.upper, display_numbers=F, filename=fn.cm)
                ## cat('------------ ', dynamicWidthHM(length(global.param$grp), unit='in'), '\n')

                plotCorrMat(lower=input$cm.lower, upper=input$cm.upper, display_numbers=input$cm.numb, filename=fn.cm, width=dynamicWidthHM(length(global.param$grp), unit='in'), height=dynamicWidthHM(length(global.param$grp), unit='in') )


                })
            }
            ############################################################
            ##                   heatmap
            ## require at least three significant hits
            ############################################################
            if(input$export.hm){

                ##filter.res()
                ##res = global.results$filtered

                if(nrow(res) >= 3){
                    withProgress(message='Exporting', detail='heatmap',{
                    fn.hm <- paste(global.param$session.dir, 'heatmap.pdf', sep='/')
                    ## heatmap title
                    hm.title <- paste('filter:', global.param$filter.type, ' / cutoff:', global.param$filter.value, sep='')
                    hm.title <- paste(hm.title, '\nsig / total: ', nrow(res), ' / ', nrow( global.results$data$output ), sep='')

                    # column annotation
                    if(!is.null(global.input$cdesc))
                      hm.cdesc <- global.input$cdesc[global.param$cdesc.selection,  ]
                    else
                      hm.cdesc <- NULL
                    
                    if(input$hm.max){
                        ##plotHM(res=res, grp=global.param$grp, grp.col=global.param$grp.colors, grp.col.legend=global.param$grp.colors.legend,  hm.clust=input$hm.clust, hm.title=hm.title, hm.scale=input$hm.scale, cellwidth=ifelse(ncol(res)<40, 40, 20), fontsize_row=input$cexRow, fontsize_col=input$cexCol, max.val=input$hm.max.val, style=global.param$which.test, filename=fn.hm, cellheight=min( dynamicHeightHM( nrow(global.results$filtered)), 1500 )/nrow(global.results$filtered))
                        #plotHM(res=res, grp=global.param$grp, grp.col=global.param$grp.colors, grp.col.legend=global.param$grp.colors.legend,  hm.clust=input$hm.clust, hm.title=hm.title, hm.scale=input$hm.scale, fontsize_row=input$cexRow, fontsize_col=input$cexCol, max.val=input$hm.max.val, style=global.param$which.test, filename=fn.hm, width=dynamicWidthHM(length(global.param$grp), unit='in'), height=min( dynamicHeightHM( nrow(global.results$filtered), unit='in'), 20 ), cdesc=hm.cdesc, cdesc.grp=global.param$grp.gct3)
                        plotHM(res=res, grp=global.param$grp, grp.col=global.param$grp.colors, grp.col.legend=global.param$grp.colors.legend,  hm.clust=input$hm.clust, hm.title=hm.title, hm.scale=input$hm.scale, fontsize_row=input$cexRow, fontsize_col=input$cexCol, max.val=input$hm.max.val, style=global.param$which.test, anno.col=global.parma$anno.col, anno.col.color=global.param$anno.col.color)
                      
                    } else {
                      plotHM(res=res, grp=global.param$grp, grp.col=global.param$grp.colors, grp.col.legend=global.param$grp.colors.legend,  hm.clust=input$hm.clust, hm.title=hm.title, hm.scale=input$hm.scale, fontsize_row=input$cexRow, fontsize_col=input$cexCol, style=global.param$which.test, anno.col=global.param$anno.col, anno.col.color=global.param$anno.col.color)
                      
                      ##plotHM(res=res, grp=global.param$grp, grp.col=global.param$grp.colors, grp.col.legend=global.param$grp.colors.legend,  hm.clust=input$hm.clust, hm.title=hm.title, hm.scale=input$hm.scale, fontsize_row=input$cexRow, fontsize_col=input$cexCol, style=global.param$which.test, filename=fn.hm,  width=dynamicWidthHM(length(global.param$grp), unit='in'), height=min( dynamicHeightHM( nrow( global.results$filtered ), unit='in'), 20 ), cdesc=hm.cdesc, cdesc.grp=global.param$grp.gct3)
                    }
                    })
                } ## end if nrow(res)>3

            }


            ############################################################
            ##                   volcanos
            ############################################################
            if(input$export.volc & !(global.param$which.test %in% c('mod F', 'none'))){
                withProgress(message='Exporting', detail='volcano plot',{
                fn.volc <- paste(global.param$session.dir, 'volcano.pdf', sep='/')

                    pdf(fn.volc, height=11, width=11)
                    for(j in 1:length(grp.comp)){
                        local({
                            my_j=j
                            ##plotVolcano(grp.comp[my_j], max.logP=input$max.logP)
                            plotVolcano(grp.comp[my_j])
                        })
                    }
                dev.off()
                })
            }
            ############################################################
            ##                 PCA
            ############################################################
            if(input$export.pca){

                if(nrow(res) >  3 && ncol(res) > 2){
                    withProgress(message='Exporting', detail='pca',{
                        fn.volc <- paste(global.param$session.dir, 'pca.pdf', sep='/')
                        pdf(fn.volc, height=5, width=15)
                        pca=plotPCA()
                        dev.off()

                    })
                }
            }

            ############################################################
            ##                 PCA Loadings as excel sheet
            ############################################################

            if(input$export.pca.loadings){
                    withProgress(message='Exporting', detail='PCA loadings as Excel sheet',{
                            if(is.null(global.results$data) | is.na(global.param$filter.value) | is.null(global.results$pca)) return()
                            if(!is.null(error$msg)) return()
                            if(length(global.param$grp) < 3) return()

                            pca <- global.results$pca
                            pca.loadings <- as.data.frame(pca$loadings)

                            ## generate_filename
                            fn.tmp <- sub(' ','_', paste(global.param$session.dir, '/', 'pcaLoadings_', sub(' ', '_',global.param$which.test),  ifelse(global.param$log.transform != 'none', paste('_', global.param$log.transform, sep=''), '_'), ifelse(global.param$norm.data != 'none', paste('_', global.param$norm.data, sep=''), '_'), ifelse(input$repro.filt=='yes', paste('_reprofilt', sep=''), '_'), sub(' .*', '', Sys.time()),".xlsx", sep=''))
                            global.param$pcaloadings.ExcelFileName <- fn.tmp
                            ## Excel
                            WriteXLS("pca.loadings", ExcelFileName= fn.tmp, SheetNames= "pcaLoadings", row.names=T, BoldHeaderRow=T, AutoFilter=T)

                            cat("-- writing pca_loadings --")

                    })

            }
            ############################################################
            ##                    boxplots
            ############################################################
            if(input$export.box){
                withProgress(message='Exporting', detail='box plot',{
                    fn.box <- paste(global.param$session.dir, 'boxplots_unnormalized.pdf', sep='/')

                    ###############################################
                    ## unnormalized ratios
                    if(is.null(global.results$table.log))
                        tab <- data.frame(global.input$table)
                    else
                        tab <- data.frame(global.results$table.log)

                    ## id column
                    id.col.value <- global.param$id.col.value
                    ## group vector
                    grp <- global.param$grp
                    ## group colors
                    grp.col <- global.param$grp.colors
                    grp.col.leg <- global.param$grp.colors.legend

                    # update selection
                    tab <- tab[, c(id.col.value, names(grp))]
                    grp.col <- grp.col[names(grp)]
                    grp.col.leg <- grp.col.leg[unique(grp)]
                    
                    
                    ## pdf
                    #pdf(fn.box, 12, max(3, .6*ncol(global.input$table)))
                    pdf(fn.box, 12, max(3, .6*length(global.param$grp)))
                    makeBoxplot(tab, id.col.value, grp, grp.col, grp.col.leg)
                    dev.off()

                    ################################################
                    ## normalized
                    if(!is.null(global.results$table.norm) ){
                        fn.box <- paste(global.param$session.dir, paste('boxplots_', global.param$norm.data,'.pdf', sep=''), sep='/')
                        ## normalized ratios
                        tab <- data.frame(global.results$table.norm)
                        
                        ## pdf
                        pdf(fn.box, 12, max(3, .6*length(global.param$grp)))
                        makeBoxplot(tab, id.col.value, grp, grp.col, grp.col.leg)
                        dev.off()
                    }
                })
            }

            ############################################################
            ##                profile plots
            ############################################################
             if(input$export.profile){
                withProgress(message='Exporting', detail='profile plot',{
                    fn.profile <- paste(global.param$session.dir, 'profile_plot.pdf', sep='/')

                    ###############################################
                    ## unnormalized ratios
                    if(is.null(global.results$table.log))
                        tab <- data.frame(global.input$table)
                    else
                        tab <- data.frame(global.results$table.log)

                    ## id column
                    id.col.value <- global.param$id.col.value
                    ## group vector
                    grp <- global.param$grp
                    ## group colors
                    grp.col <- global.param$grp.colors
                    grp.col.leg <- global.param$grp.colors.legend

                    # update selection
                    tab <- tab[, c(id.col.value, names(grp))]
                    grp.col <- grp.col[names(grp)]
                    grp.col.leg <- grp.col.leg[unique(grp)]
                    
                    
                    
                    ################################################
                    ## normalized
                    if(!is.null(global.results$table.norm) ){
                        pdf(fn.profile, 14, 7)
                        par(mfrow=c(1,2))

                        ## normalized ratios
                        tab.norm <- data.frame(global.results$table.norm)
                        makeProfileplot(tab, id.col.value, grp, grp.col, grp.col.leg, main='Before normalization')
                        makeProfileplot(tab.norm, id.col.value, grp, grp.col, grp.col.leg, main=paste(global.param$norm.data, 'normalized'))
                        dev.off()

                    } else{
                        pdf(fn.profile, 7, 7)
                        par(mfrow=c(1,1))
                        makeProfileplot(tab, id.col.value, grp, grp.col, grp.col.leg, main='unnormalized')
                        dev.off()

                    }
                })
            }

            ############################################################
            ##                 multi scatter
            ############################################################
            if(input$export.ms){
                withProgress(message='Exporting', detail='multiscatter',{


                    fn.ms <- paste(global.param$session.dir, 'multiscatter.pdf', sep='/')
                    pdf(fn.ms, height=120*ncol(global.input$table)*(11/800), width=120*ncol(global.input$table)*(11/800))
                    plotMultiScatter( define.max=input$ms.max, min.val=input$ms.min.val, max.val=input$ms.max.val )
                    dev.off()
                })
            }

            ############################################################
            ##                   p-values
            ############################################################
            if(input$export.phist & global.param$which.test != 'none'){

                withProgress(message='Exporting', detail='p-value histogram',{
                fn.pval <- paste(global.param$session.dir, 'histogram_P-values.pdf', sep='/')

                ## unfiltered results
                res.all = global.results$data$output

                pdf(fn.pval, 10, 5*ifelse( global.param$which.test != 'mod F', length(grp.comp), 1 ))

                ############################################
                ## mod T
                if(!(global.param$which.test %in% c('mod F'))){
                    par(mfrow=c(length(grp.comp),1))
                    for(g in grp.comp){
                        pval <- res.all[, paste('P.Value', g, sep='.')]
                        hist(pval, breaks=50, main=paste('Histogram of P-values (N=', sum(!is.na(pval)), ')',sep=''), xlab='P-value', cex.main=2.2, cex.axis=2, cex.lab=2, col='darkblue', border=NA)
                        legend('top', legend=g, cex=2)
                    }
                 ############################################
                 ## mod F
                } else {

                    pval <- res.all[, paste('P.Value')]
                    hist(pval, breaks=50, main=paste('Histogram of P-values (N=', sum(!is.na(pval)), ')',sep=''), xlab='P-value', cex.main=2.2, cex.axis=2, cex.lab=2, col='darkblue', border=NA)
                }
                dev.off()
                })

            }

            #########################################################
            ##               Excel sheet
            #########################################################
            if(input$export.excel){
                withProgress(message='Exporting', detail='Excel sheet',{

                    res.comb <- global.results$data$output
                    tmp <- sort(global.param$grp)

                    ## append annotation columns
                    if(!is.null(global.input$table.anno))
                        res.comb <- left_join(res.comb,  global.input$table.anno, 'id')
                        ##res.comb <- cbind(res.comb, Annotation.starts.here=rep('', nrow(res.comb)), global.input$table.anno)

                    expDesign <- data.frame(Column=names(tmp), Experiment=tmp)

                    ## generate_filename
                    fn.tmp <- sub(' ','_',
                                  paste(global.param$session.dir, '/', 'results_',
                                        sub(' ', '_',global.param$which.test), '_',
                                        ifelse(global.param$log.transform != 'none', paste( global.param$log.transform, '_', sep=''), '_'),
                                        ifelse(global.param$norm.data != 'none', paste( global.param$norm.data, '_', sep=''), '_'),
                                        ifelse(input$repro.filt=='yes', paste(global.param$filt.data, sep=''), '_'),
                                        sub(' .*', '', Sys.time()),".xlsx", sep=''))

                    global.param$ExcelFileName <- fn.tmp
                    ## Excel
                    ## WriteXLS(c('res.comb', 'expDesign'), ExcelFileName=fn.tmp, FreezeRow=1, FreezeCol=1, SheetNames=c('modT', 'class vector'), row.names=F, BoldHeaderRow=T, AutoFilter=T)
                    WriteXLS(c('res.comb', 'expDesign'), ExcelFileName=fn.tmp, FreezeRow=1, FreezeCol=1, SheetNames=c(global.param$which.test, 'class vector'), row.names=F, BoldHeaderRow=T, AutoFilter=T)

                })
            }

            #########################################################
            ##          save session parameters
            #########################################################
            global.input.imp <- reactiveValuesToList(global.input)
            global.param.imp <- reactiveValuesToList(global.param)
            global.results.imp <- reactiveValuesToList(global.results)

            #########################################################
            ## save current plotting parameters
            #########################################################

            ## heatmap
            global.plotparam$hm.scale <- input$hm.scale
            global.plotparam$hm.clust <- input$hm.clust
            global.plotparam$hm.cexCol <- input$hm.cexCol
            global.plotparam$hm.cexRow <- input$hm.cexRow
            global.plotparam$hm.max <- input$hm.max
            global.plotparam$hm.max.val <- input$hm.max.val
            ## pca
            global.plotparam$pca.x <- input$pca.x
            global.plotparam$pca.y <- input$pca.y
            global.plotparam$pca.z <- input$pca.z
            ## multiscatter
            global.plotparam$ms.max <- input$ms.max
            global.plotparam$ms.min <- input$ms.min
            global.plotparam$ms.max <- input$ms.max
            ## correlation matrix
            global.plotparam$cm.upper <- input$cm.upper
            global.plotparam$cm.lower <- input$cm.lower
            global.plotparam$cm.numb <- input$cm.numb

            ## convert to list
            global.plotparam.imp <- reactiveValuesToList(global.plotparam)

            ## volcano coordinates
            volc.imp <-  reactiveValuesToList(volc)

            #################################
            ## save as R-object
            ## no label present
            if(is.null(global.param$label)){

                fn.tmp <- paste(global.param$session.dir, paste('session_', gsub('( |\\:)', '-', Sys.time()), '.RData', sep=''), sep='/')

                save(global.input.imp, global.param.imp, global.results.imp, global.plotparam.imp, volc.imp, file=fn.tmp)
                fn.zip <- paste( gsub('\\:', '', gsub(' ','-', gsub('-','',Sys.time()))),'.zip', sep='')

            }
            ## label present
            if(!is.null(global.param$label) | nchar(global.param$label) == 0){

                fn.tmp <- paste(global.param$session.dir, paste(global.param$label, paste('_session_', gsub('( |\\:)', '-', Sys.time()), '.RData', sep=''), sep=''), sep='/')

                save(global.input.imp, global.param.imp, global.results.imp, global.plotparam.imp, volc.imp, file=fn.tmp)
                fn.zip <- paste( global.param$label, '_', gsub('\\:', '', gsub(' ','-', gsub('-','',Sys.time()))),'.zip', sep='')

            }

            #########################################################
            ##   export parameters as small text file
            #########################################################
            params <- global.param.imp[c('log.transform', 'norm.data', 'filt.data',  'repro.filt.val', 'sd.filt.val', 'which.test', 'filter.type', 'filter.value')]

            params.txt <- unlist(lapply(params, paste, collapse=';'))
            params.txt <-  paste(names(params.txt),params.txt, sep='\t')

            ## hostname
            params.txt <- c(paste('## hostname: ', session$clientData$url_hostname, sep=''), '' , params.txt)

            ## user names
            params.txt <- c(paste('## user: ', session$user, sep=''), '' , params.txt)

            ## session id
            params.txt <- c(paste('## session id: ', global.param$session ), params.txt)

            ## version
            params.txt <- c(paste('## ', APPNAME, ' (v', VER, ')', sep='') , params.txt)

            ## add date and time
            params.txt <- c(paste('##', as.character(Sys.time())), params.txt)

            ## export
            writeLines(params.txt, con=paste(global.param$session.dir, 'params.txt', sep='/'))

            ## ############################################################
            ##            create an archive
            ## ############################################################
            fn.all <- grep('pdf$|xlsx$|txt$|gct$|RData$',  dir(global.param$session.dir) , value=T, ignore.case=T)
	          fn.all.abs <- grep('pdf$|xlsx$|txt$|gct$|RData$', dir(global.param$session.dir, full.names=T, ignore.case=T), value=T)

            ## #################################################
            ## handle special characters in file names
            ## #################################################

            ## for windows: "
            fn.all.abs <- paste('"',fn.all.abs,'"', sep='')
            ## unix/linux: '
            fn.all <- paste('\'',fn.all,'\'', sep='')

            ###################################################
            ## run command
	          if(OS == 'Windows')
                system( paste('zip -0 ', paste(global.param$session.dir, fn.zip, sep='/'), ' ', paste(fn.all.abs, collapse=' '),sep='') )
            else
                system( paste('cd ', global.param$session.dir,' && zip -0 ', fn.zip, ' ', paste(fn.all, collapse=' '), sep='') )

            ## store file.name
	          global.param$zip.name=fn.zip
            ## cat('test16\n')
	          
            #----------------------------------------------------------------------
            #             clean up
	          # remove archived files: all files except RData
	          fn.rm <- gsub('"|\'', '', fn.all.abs[-grep('\\.RData[\\"|\']$', fn.all.abs)]) 
            file.remove(fn.rm)
            
            # keep only the latest zip file
            all.zip <- dir(global.param$session.dir, pattern = '\\.zip$', full.names = T)
            if(length(all.zip) > 1){
              all.zip <- all.zip[ -which.max(file.info(all.zip)$ctime) ]
              file.remove(all.zip)
            }
            #

            # remove Rdata file, if session won;t be saved            
            if(!input$export.save.session)
                file.remove(fn.tmp)

             ###############################################################
             ## flag export
             global.results$export.results=T

             ## redirect to the same panel
             updateTabsetPanel(session, 'mainPage', selected='Export')
        })


        ###################################################################################
        ##
        ##                             Do the actual computation
        ##
        ###################################################################################

        ###############################################
        ## 4) ID column
        ##  - make unique ids
        ##  - determine id type
        ##  - map to gene names
        ##  - initialize group assignment
        observeEvent( input$id.col ,{

            if( is.null( global.input$table) | is.null(input$id.col.value) ) return()

            ## store name of id column
            global.param$id.col.value <- input$id.col.value

            ## update 'input$id.col'
            global.input$id.col <- input$id.col

            ## ###########################################
            ## check the id column
            tab <- global.input$table

            ## ###########################################
            ## make sure the ids are unique
            ##ids <- make.names( make.unique(as.character(tab[, global.param$id.col.value]), sep='_') )
            ids <- make.unique( as.character(tab[, global.param$id.col.value] ), sep='_')

            ## replace values in id column
            tab[, global.param$id.col.value] <- ids
            ##tab <- data.frame(id=ids, tab)
            ## use id as rownames
            rownames(tab) <- ids

            ## ############################################
            ## map to gene names
            ## ############################################
            map.res <- mapIDs(ids)
            global.results$keytype <- map.res$keytype
            #global.results$id.map <- data.frame(id=ids, map.res$id.map)
            global.results$id.map <- map.res$id.map
            
            #View(data.frame(id=ids, map.res$id.map))
            
            ########################################
            ## store
            global.input$table <- tab
            
            # #######################
            # flag
            global.param$id.done <- T

        })

        #################################################################################
        # IMPORT DATA
        ## 2) - upload file
        ##    - determine file type
        ##    - import file
        ##    - user folders are create here
        ##    - extract label from filename
        #################################################################################
        observeEvent( input$file, {

            global.input$id.col <- 0

            ########################################
            ## generate session ID and prepare data
            ## directory
            #########################################
            ## ## if there is an user ID...
            if(!is.null( global.param$user )){
                 ## create user directory, if not present already
                 if(!dir.exists(paste(DATADIR, global.param$user, sep='')))
                     dir.create(paste(DATADIR, global.param$user, sep=''))

                 global.param$session.dir <- paste(DATADIR, global.param$user,'/' ,global.param$session, '/', sep='')

            } else{
                 global.param$session.dir <- paste(DATADIR, global.param$session, '/' ,sep='')
            }
            ## create directory on server to store the results
            dir.create(global.param$session.dir)

            ## ############################################
            ## copy the input file to the session folder
            fn <- paste0( global.param$session.dir, input$file$name)
            file.copy(input$file$datapath, fn)

            ## ###############################
            ## generate label
            fn.split <- unlist(strsplit( sub('.*/', '', fn), '_'))
            if(length(fn.split) > 1){
                label <- fn.split[1]
                if(nchar(label) > 10)
                    label <- chopString(label, 10, add.dots=F)
                global.param$label <- label
            } else {
                global.param$label <- chopString( sub('.*/', '', fn) , 10, add.dots=F)
            }

            ########################################################
            # file import
            
            # ###############################################
            #                     GCT 1.2
            if( length( grep( '^\\#1\\.2', readLines(fn,n=1))) > 0){
                  tab <- read.delim( fn, stringsAsFactors=F, na.strings=NASTRINGS, skip=2)
                  
                  ## shorten column names and store together with the original names
                  colnames.tmp <- chopString(colnames(tab), STRLENGTH)
                  names(colnames.tmp) <- colnames(tab)
                  
                  ## store values
                  global.input$table <- global.input$table.org <- tab
                  global.input$file <- input$file
                  global.input$table.colnames <- colnames.tmp
                  
                  rm(tab, colnames.tmp)
            # ###################################################      
            #                     GCT 1.3
            } else if( length( grep( '^\\#1\\.3', readLines(fn,n=1))) > 0){
              
              # parse gct file
              gct <- parse.gctx(fn)
              
              # expression table
              tab <- data.frame(id=gct@rid, gct@rdesc, gct@mat, stringsAsFactors = F)
              
              colnames.tmp <- chopString(colnames(tab), STRLENGTH)
              names(colnames.tmp) <- colnames(tab)
              
              ## store values
              global.input$table <- global.input$table.org <- tab
              global.input$file <- input$file
              global.input$table.colnames <- colnames.tmp
              
              # meta data
              global.input$rdesc <- gct@rdesc
              
              global.input$cdesc <- gct@cdesc
              rownames(global.input$cdesc) <- make.names(rownames(global.input$cdesc)) # convert to proper names
              
              global.param$cdesc.all <- global.param$cdesc.selection <- colnames(global.input$cdesc)
              
              # id column 
              global.param$id.col.value='id'
              global.param$id.done=T
              
              # robustify ids
              ids <- make.unique( as.character(tab[, global.param$id.col.value] ), sep='_')
              ## replace values in id column
              tab[, global.param$id.col.value] <- ids
              ## use id as rownames
              rownames(tab) <- ids
              
              ## #################
              ## map to gene names
              map.res <- mapIDs(ids)
              global.results$keytype <- map.res$keytype
              global.results$id.map <- map.res$id.map
              
              # flag
              global.param$file.gct3 <- T
              
              rm(tab, colnames.tmp)
              
            # ##########################################################  
            #                    other text file
            } else {

                ## ################################
                ## determine the separator
                tab.sep=NULL
                ## try to figure out the separator, DON'T USE THE HEADER FOR THAT
                ## use the fourth row instead (should be data)
                for(s in SEPARATOR){

                    ##tab <- read.table(fn, sep=s, header=T, stringsAsFactors=F, nrows=1, skip=3)
                    tab <- read.table(fn, sep=s, header=F, stringsAsFactors=F, nrows=1, skip=4)

                    if(length(tab) > 1){
                        global.param$tabsep <- s
                        break;
                    }
                }
                ## #########################################################
                ## import the table: txt
                if( global.param$tabsep == '\t'){
                    tab <- read.delim( fn, stringsAsFactors=F, na.strings=NASTRINGS)
                } else {
                    tab <- read.table( fn, sep=global.param$tabsep, header=T, stringsAsFactors=F, na.strings=NASTRINGS, quote = "\"", dec = ".", fill = TRUE, comment.char = "")
                }
                ## shorten column names and store together with the original names
                colnames.tmp <- chopString(colnames(tab), STRLENGTH)
                names(colnames.tmp) <- colnames(tab)
                
                ## store values
                global.input$table <- global.input$table.org <- tab
                global.input$file <- input$file
                global.input$table.colnames <- colnames.tmp
                
                rm(tab, colnames.tmp)
            }## end if GCT

            # flag
            global.param$file.done <- T
            
        })

   
        ## ###############################################################
        ##
        ## import saved session from the server via drop down menu
        ##
        ## ###############################################################
        observeEvent(input$session.browse.import, {

            ## ###############################
            ## import workspace
            ## identify selected session file
            load(global.param$saved.sessions[[ input$session.browse ]])
           
            # for backwards compatibility
            if(!exists('global.plotparam.imp')){
              global.plotparam.imp <- plotparams.imp 
            }
            
            ## ###############################
            ## assign the imported values to the global reactive variables
            for(i in names(global.input.imp)){
                global.input[[i]] <- global.input.imp[[i]]
            }
            for(i in names(global.param.imp)){
                global.param[[i]] <- global.param.imp[[i]]
                ##cat(i,'\t',global.param[[i]], '\n')
            }
            for(i in names(global.results.imp)){
                global.results[[i]] <- global.results.imp[[i]]
            }
            for(i in names(global.plotparam.imp)){
                global.plotparam[[i]] <- global.plotparam.imp[[i]]
            }
            for(i in names(volc.imp)){
                volc[[i]] <- volc.imp[[i]]
            }
            
            ## backwards compatibility II
            if(is.null(names(global.param$grp.colors)))
              names(global.param$grp.colors) <- names(global.param$grp.colors)
          
          

            ## ##############################################################
            ##                       update filter
            ## ##############################################################
            if(global.param$filter.type == 'top.n'){
                updateNumericInput(session, inputId='filter.value.top.n', value=global.param$filter.value)
            }
            if(global.param$filter.type == 'nom.p'){
                updateNumericInput(session, inputId='filter.value.nom.p', value=global.param$filter.value)
            }
            if(global.param$filter.type == 'adj.p'){
                updateNumericInput(session, inputId='filter.value.adj.p', value=global.param$filter.value)
            }

            ## ##############################################################
            ##                 update plotting parameters
            ## ##############################################################

            ## heatmap
            updateNumericInput(session, inputId='cexCol', value=global.plotparam$hm.cexCol)
            updateNumericInput(session, inputId='cexRow', value=global.plotparam$hm.cexRow)
            updateSelectInput(session, inputId='hm.scale', selected=global.plotparam$hm.scale)
            updateSelectInput(session, inputId='hm.clust', selected=global.plotparam$hm.clust)
            updateNumericInput(session, inputId='hm.max.val', value=global.plotparam$hm.max.val)
            updateCheckboxInput(session, inputId='hm.max', value=global.plotparam$hm.max)
            ## hm clust missing

            ## PCA
            updateSelectInput(session, inputId='pca.x', selected=global.plotparam$pca.x)
            updateSelectInput(session, inputId='pca.y', selected=global.plotparam$pca.y)
            updateSelectInput(session, inputId='pca.z', selected=global.plotparam$pca.z)

            ## multiscatter
            updateCheckboxInput(session, inputId='ms.max', value=global.plotparam$ms.max)
            updateNumericInput(session, inputId='ms.max.val', value=global.plotparam$ms.max.val)
            updateNumericInput(session, inputId='ms.min.val', value=global.plotparam$ms.min.val)

            ## correlation matrix
            updateCheckboxInput(session, inputId='cm.numb', value=global.plotparam$cm.numb)
            updateSelectInput(session, inputId='cm.upper', selected=global.plotparam$cm.upper)
            updateSelectInput(session, inputId='cm.lower', selected=global.plotparam$cm.lower)



            ##################################
            ## set flags
            global.param$file.done=T
            global.param$id.done=T
            global.param$session.imported=T
            global.param$analysis.run=T

            global.param$session.import.init=T

            global.results$export.results=F

            ##################################
            ## clean up
            rm(global.input.imp, global.param.imp, global.results.imp, global.plotparam.imp, volc.imp)

            ###################################################################
            ##            insert the panels for the volcanos
            ###################################################################
            if(!(global.param$which.test %in% c('mod F', 'none'))){
                ins.volc()
            }
            ###################################################################
            ##            insert the panels for the scatterplots
            ###################################################################
            ins.scat()

            ## #####################################
            ## generate id.map for compatibility with < v0.7.0
            if(is.null(global.results$id.map)){

                ## ############################################
                ## map to gene names
                ## ############################################
                ##ids <- rownames( global.results$data$output)
                ##ids <- global.results$data$output[, global.param$id.col.val]
                ids <- global.results$data$output[, 'id']

                map.res <- mapIDs( ids )

                global.results$keytype <- map.res$keytype
                global.results$id.map <- data.frame(id=ids,  map.res$id.map )
                global.results$data$output <- left_join(global.results$data$output, global.results$id.map, 'id')
                ##rownames(global.results$data$output) <- global.results$id.map$id.concat
            }
        })

        #
        
        ## ###############################################################
        ## 4c)       upload experimental design file
        ##
        ## - divide input table into expression and annotation columns
        observeEvent( input$exp.file, {

            ## reset error message
            error$msg <- NULL

            ## ###########################
            ## copy file into sessions folder
            fn <- paste0( global.param$session.dir, input$exp.file$name)
            file.copy(input$exp.file$datapath, fn)


            ## ##########################
            ## read the experimental design file
            ##grp.file <- read.delim(input$exp.file$datapath, header=T, stringsAsFactors=F)
            grp.file <- read.delim(input$exp.file$datapath, header=T, stringsAsFactors=F)
            Column.Name <- grp.file$Column.Name
            Experiment <- as.character(grp.file$Experiment)
            
            ## #############################################################
            ## index on non-empty 'Experiment' rows
            exprs.idx <- which(nchar(Experiment) > 0 )


            ## ###############################
            ## update label
            fn.split <- unlist(strsplit( sub('.*/', '', fn), '_'))
            if(length(fn.split) > 1){
                label <- fn.split[1]
                if(nchar(label) > 10)
                    label <- chopString(label, 10, add.dots=F)
                global.param$label <- paste(global.param$label, label )
            } else {
                label <- chopString( sub('.*/', '', fn) , 10, add.dots=F)
                global.param$label <- paste(global.param$label, label, sep='-')
            }

           #cat('\n',nchar(Experiment),'\n')
            ###############################################################
            ##
            ## - do some sanity checks
            ## - separate expression data from annotation columns
            ##
            ###############################################################

            ## names in the exp design file do not match to the table
            if( sum( Column.Name != colnames(global.input$table)) > 0 ){
                error$msg <- 'Experimental design files does not match the table you have uploaded!'
                return()
            }

            ## not an experimental design file
            if( sum( colnames(grp.file) %in% c('Column.Name', 'Experiment'), na.rm=T) != 2 )  {
                error$msg <- 'This is not an experimental desgin file!\n\nThe file should contain two columns (Column.Name, Experiment)!'
                return()
            }
            ## 'empty' file
            if( sum( nchar(Experiment) > 0, na.rm=T ) == 0 | sum(!is.na( Experiment) == 0) ){
                error$msg <- 'No experiments defined!'
                return()
            }
            ##cat('L=', sum( Column.Name[ exprs.idx ] %in%  colnames(global.input$table)))
            ## column names specified in exp design file not found in table
            ##if( sum( Column.Name[ exprs.idx ] %in%  colnames(global.input$table)) != length(exprs.idx) ){
            ##    error$msg <- 'Column names in the experimental design file cannot be found in the data table!'
            ##    return()
           ## }
            ## check whether there are at least 2 replicates per group
            num.rep=table(Experiment[exprs.idx])
            ##if(min(num.rep) == 1){
            ##    error$msg <- paste('No replicate measurements defined!')
            ##    return()
            ##}

            ## ################################
            ## ANNOTATION: extract empty cells
            ## - corresponding columns will be carried over as
            ##   annotation columns in the result file
            grp.anno <- grp.file[which(nchar( Experiment) == 0 ), ]
            grp.anno <- setdiff( grp.anno$Column.Name, global.param$id.col.value )

            if(length(grp.anno)>0)
                global.input$table.anno <- data.frame(id=global.results$id.map[, 'id'], global.input$table[ , grp.anno])

            ## ################################
            ## EXPRESSION
            ## - extract all non-empty cells in the 'Experiment' column
            grp.exprs <- grp.file[exprs.idx, ]

            ## order alphabetically to make coloring consistent
            grp.exprs <- grp.exprs[order(grp.exprs$Experiment), ]

            ## class vector
            grp=grp.exprs$Experiment
            names(grp)=grp.exprs$Column.Name

            ## update input table, keep id and expression columns
            global.input$table <- global.input$table[ , c(global.param$id.col.value, names(grp))]

            ################################
            ## update number of groups
            global.param$N.grp <- length(unique( na.omit(grp)) )
            
            ## store group assignment
            global.param$grp <- global.param$grp.all <- grp
            global.param$grp.comp.all <- global.param$grp.comp <- unique(grp)
            
            ## group colors
            grp.col <- rep(GRPCOLORS[1], length(grp))
            names(grp.col) <- names(grp)
            
            for(i in 2:length(unique(grp))) grp.col[ which(grp == unique(grp)[i]) ] <- GRPCOLORS[i]
            global.param$grp.colors <- grp.col

            ## group colors for figure legend
            idx <- !duplicated(grp)
            grp.col.legend = grp.col[idx]
            names(grp.col.legend) <- grp[idx]
            global.param$grp.colors.legend <- grp.col.legend

            ## all done
            global.param$grp.done = T

            ## save column name used as 'id'
            ##global.param$id.col.value = input$id.col.value

        })

        # #####################################
        #     determine groups to test
        # 
        observeEvent(input$which.test, {
          
          test <- input$which.test
          
          ## #############################################
          ## specify which comparisons should be performed
          if(test %in% c('One-sample mod T', 'mod F', 'none')){
            ## each group separetely
            groups.comp <- global.param$grp.all
            
            #if( global.param$run.test == 0  ){
            global.param$grp.comp.all <- groups.comp
            #}
            global.param$grp.comp.selection <- groups.comp
            global.param$grp.selection <- groups.comp
          }
          ## #############################################
          ## pairwise combinations for 2-sample test
          if(test == 'Two-sample mod T'){
            
            ## all pairwise combinations
            groups.unique <- unique(global.param$grp.all)
            
            groups.comp <- c()
            count=1
            for(i in 1:(length(groups.unique)-1))
              for(j in (i+1):length(groups.unique)){
                ## order alphabetically
                groups.tmp <- sort( groups.unique[c(i,j)] )
                ##groups.comp[count] <- paste(groups.unique[i], groups.unique[j], sep='.vs.')
                groups.comp[count] <- paste(groups.tmp[1], groups.tmp[2], sep='.vs.')
                count <- count+1
              }
            #if( global.param$run.test == 0  ){
              global.param$grp.comp.all <- groups.comp
            #}
            global.param$grp.comp.selection <- groups.comp
            #global.param$grp.selection <- 
          }
          
          
        })
        
        
        
        ################################################################################
        ##
        ##                once the 'run test' button was pressed...
        ##
        ## - log-transform (optional)
        ## - normalization (optional)
        ## - filter data(optional)
        ## - test (optional)
        ##
        ################################################################################
        observeEvent(input$run.test, {

            ## reset any error messages
            error$msg <- NULL

            ##global.param$update.ppi.select <- FALSE

            global.input$run.test <- input$run.test

            ##Vie(global.input$table)

            ###########################################
            ## - store which test has been used
            ## - needed to supress volcanos for mod F test
            global.param$which.test <- input$which.test
            global.param$filt.data <- input$filt.data
            global.param$repro.filt.val <- input$repro.filt.val
            global.param$sd.filt.val <- input$sd.filt.val
            global.param$norm.data <- input$norm.data
            global.param$log.transform <- input$log.transform


            #####################################################################
            ## if the 'Run test' - button has been pressed for the first time,
            ## store a copy of the original input matrix (unnormalized, unfiltered)
            if(!(global.param$analysis.run)){
                global.input$table.org <- global.input$table
                tab <- data.frame(global.input$table)
            }
            if(global.param$analysis.run){
                tab <- data.frame(global.input$table.org)
            }
            ##View(tab)
            ## id column
            id.col = global.param$id.col.value
            
            ## all group labels
            groups=global.param$grp.selection
            global.param$grp <- groups
            ##View(tab[, id.col])

            ###############################################
            ## initialize values for normalized and filtered matrix
            global.results$table.log=NULL
            global.results$table.norm=NULL
            global.results$table.repro.filt=NULL
            global.results$table.sd.filt=NULL
            global.results$values.filtered=NULL
            global.results$pca=NULL

            ## ##############################################################################
            ##
            ##           set up the analysis pipeline
            ##
            ## ##############################################################################
            test = global.param$which.test
            norm.data = global.param$norm.data
            filt.data = global.param$filt.data
            repro.filt.val = global.param$repro.filt.val
            sd.filt.val = global.param$sd.filt.val
            log.trans = global.param$log.transform

            # update group selection
            groups.comp = global.param$grp.comp.selection
            global.param$grp.comp = groups.comp
            
            # ## #############################################
            # ## specify which comparisons should be performed
            # if(test %in% c('One-sample mod T', 'mod F', 'none')){
            #     ## each group separetely
            #     groups.comp <- global.param$grp
            #     global.param$grp.comp <- groups.comp
            # }
            # ## #############################################
            # ## pairwise combinations for 2-sample test
            # if(test == 'Two-sample mod T'){
            # 
            #     ## all pairwise combinations
            #     groups.unique <- unique(global.param$grp)
            # 
            #     groups.comp <- c()
            #     count=1
            #     for(i in 1:(length(groups.unique)-1))
            #         for(j in (i+1):length(groups.unique)){
            #             ## order alphabetically
            #             groups.tmp <- sort( groups.unique[c(i,j)] )
            #             ##groups.comp[count] <- paste(groups.unique[i], groups.unique[j], sep='.vs.')
            #             groups.comp[count] <- paste(groups.tmp[1], groups.tmp[2], sep='.vs.')
            #             count <- count+1
            #         }
            #     global.param$grp.comp <- groups.comp
            # }

            ## ###########################################
            ## log transformation
            ## ###########################################
            if(log.trans != 'none'){
                ids.tmp <- tab[, id.col]
                dat.tmp <- tab[, -which(colnames(tab) == id.col)]
                dat.tmp <- data.matrix(dat.tmp)

                dat.tmp[dat.tmp == 0] <- NA

                if(log.trans == 'log2'){
                    dat.tmp <- log(dat.tmp, 2)
                }
                if(log.trans == 'log10'){
                    dat.tmp <- log(dat.tmp, 10)
                }
                ## putting a dataframe around here turns out to be ESSENTIAL!! DON'T USE CBIND HERE!
                tab <- data.frame( ids.tmp, dat.tmp)
                colnames(tab)[1] <- id.col
                ## store
                global.results$table.log <- tab
            }

            ## ###########################################
            ##
            ## normalization
            ##
            ## ###########################################
            if(norm.data != 'none'){
                ##tab.org = tab
                withProgress(message='Applying normalization...', {
                    tab <- normalize.data(tab, id.col, norm.data)
                    if(is.null(dim(tab)) | unlist(tab)[1] == 'No_success'){

                        error$msg <- paste('Error in Two-component normalization:<br><br>Could not fit mixture model for data column:<br><br>', tab, '<br><br>Giving up...')

                        ## if not successful skip the rest
                        test='none'
                        repro.filt='no'

                    } else {
                        global.results$table.norm <- tab
                        ##View(tab)
                    }
                })
            }
            ## ############################################
            ##
            ## reproducibility filter
            ## - only for one-sample test
            ##
            ## ############################################
            if(filt.data == 'Reproducibility'){

                    withProgress(message='Applying reproducibility filter',  {
                        repro = my.reproducibility.filter(tab, id.col=id.col, groups, alpha=global.param$repro.filt.val)


                        tab = repro$table
                        ##View(repro$table)
                    })
                    ## store indices of filtered values in the original table
                    global.results$values.filtered <- repro$values.filtered
                    global.results$table.filt <- tab

            }
            ## #############################################
            ##
            ##   Standard deviation filter
            ##
            ## #############################################
            if(filt.data == 'StdDev'){

                 withProgress(message='Applying StdDev filter',  {
                        filt.tmp = sd.filter(tab, id.col=id.col, groups, sd.perc=global.param$sd.filt.val)
                        tab = filt.tmp$table
                        ##View(tab)
                    })
                    ## store indices of filtered values in the original table
                    global.results$values.filtered <- filt.tmp$values.filtered
                    global.results$table.filt <- tab



            }

            ## #############################################################################
            ##
            ##                                     TEST
            ##
            ###############################################################################

            ##################################
            ## two sample
            if(test == 'Two-sample mod T'){
                ##View(tab)
                withProgress(message='Two-sample test', value=0, {

                    count=0
                    ## loop over groups
                    for(g in unique(groups.comp)){

                        ## extract current groups
                        groups.tmp <- groups[ c(grep( paste('^',unlist( strsplit(g, '\\.vs\\.'))[1],'$', sep='') , groups), grep(paste('^',unlist( strsplit(g, '\\.vs\\.'))[2], '$', sep='') , groups)) ]

                        ## extract table of current group
                        tab.group <- cbind(tab[, id.col], tab[, names(groups.tmp)])
                        colnames(tab.group)[1] <- id.col


                        #############################
                        ## the actual test
                        #############################
                        res.tmp <-  modT.test.2class( tab.group, groups=groups.tmp, id.col=id.col, label=g )$output

                        if(count == 0){
                            res.comb <- res.tmp
                        } else {
                            ## make sure the order is correct
                            if(nrow(res.tmp ) != nrow(res.comb)) stop( "number of rows don't match!\n" )
                            res.tmp <- res.tmp[rownames(res.comb), ]
                            ##res.comb <- cbind(res.comb, res.tmp)
                            res.comb <- data.frame(res.comb, res.tmp, stringsAsFactors=F)
                        }
                        ##################################################
                        ## progress bar
                        incProgress(count/length(unique(groups.comp)), detail=g)
                        count=count + 1

                    }

                ##################################
                ## reorder table
                    #save(groups, res.comb, file='groups.RData')
                res.id <- res.comb$id ## id column
                res.exprs <- res.comb[, names(groups)] ## expression values
                res.test <- res.comb[, grep('^logFC|^AveExpr|^t\\.|^P\\.Value|^adj.P.Val|^Log\\.P\\.Value', colnames(res.comb))] ## test results
                res.test <- res.test[, order(colnames(res.test))]

                ## assemble new table
                res.comb <- data.frame(id=res.id, res.test, res.exprs)
                res.comb <- res.comb[rownames(tab),]

                ##global.results$data$output <- res.comb


                })
            }
            ##################################
            ## one sample
            if(test == 'One-sample mod T'){
               ## cat('test...')
                withProgress(message='One-sample test', value=0, {

                    count=0
                    ## loop over groups
                    for(g in unique(groups.comp)){

                       ## setProgress(detail=g)

                        ## extract table of current group
                        tab.group <- cbind(tab[, id.col], tab[, names(groups)[which(groups == g)]])
                        colnames(tab.group)[1] <- id.col

                        ##View(tab.group)

                        res.tmp <- modT.test( tab.group, id.col=id.col, plot=F, nastrings=NASTRINGS, label=g, na.rm=FALSE)$output
                        ##View(res.tmp)
                        if(count == 0){
                            res.comb <- res.tmp
                        } else {
                            if(nrow(res.tmp ) != nrow(res.comb)) stop( "number of rows don't match!\n" )
                            res.tmp <- res.tmp[rownames(res.comb), ]
                            res.comb <- cbind(res.comb, res.tmp)
                            ##res.comb <- merge(res.comb, res.tmp, sort=F)
                            ##res.comb <- merge(res.comb, res.tmp, by='id')
                        }

                        #############################################
                        ## progress bar
                        ##incProgress(1/length(unique(groups.comp)))
                        ##incProgress( count/length(unique(groups.comp) ), detail=g)
                        incProgress( 1/length(unique(groups.comp) ), detail=g)
                        count=count + 1
                    }
                })

                ##################################
                ## reorder columns of the table
                res.id <- res.comb$id ## id column
                res.exprs <- res.comb[, names(groups)] ## expression values
                ##View(res.exprs)
                res.test <- res.comb[, grep('^logFC|^AveExpr|^t\\.|^P\\.Value|^adj.P.Val|^Log\\.P\\.Value', colnames(res.comb))] ## test results
                ##View(res.test)
                res.test <- res.test[, order(colnames(res.test))]
                ## assemble new table
                res.comb <- data.frame(id=res.id, res.test, res.exprs, stringsAsFactors=F)
                ##res.comb <- res.comb[tab[, id.col], ]
                res.comb <- res.comb[rownames(tab),]
                ###########################################
                ## store the results
                ##global.results$data$output <- res.comb
            }

            ##################################
            ## moderated F test
            if(test == 'mod F'){
                ##cat('test F...')
                withProgress(message='moderated F-test', value=0, {
                    tab.group <- cbind(tab[, id.col], tab[, names(groups)])
                    colnames(tab.group)[1] <- id.col
                  
                    #save(groups, file='groups.RData')
                    res.comb <- modF.test( tab.group, id.col=id.col, class.vector=groups, nastrings=NASTRINGS, na.rm=FALSE)$output
                    colnames(res.comb) <- sub('^X', '', colnames(res.comb))
                })
                ##################################
                ## reorder table
                res.id <- res.comb$id ## id column
                ##res.exprs <- res.comb[, names(groups)] ## expression values
                res.exprs <- tab[, names(groups)]
                res.test <- res.comb[, grep('^logFC|^AveExpr|^F$|^P\\.Value$|^adj.P.Val$|^Log\\.P\\.Value', colnames(res.comb))] ## test results
                res.test <- res.test[, order(colnames(res.test))]

                ## assemble new table
                res.comb <- data.frame(id=res.id, res.test, res.exprs, stringsAsFactors=F)
                ##res.comb <- res.comb[tab[, id.col], ]
                res.comb <- res.comb[rownames(tab),]
                ###########################################
                ## store the results
##                global.results$data$output <- res.comb
            }

            ###################################################################
            ##
            ##                       no test
            ##
            ###################################################################
            if(test == 'none'){

                ###########################################
                ## store data matrix as test results
                ## - values are log

                tab <- data.frame(id=tab[ , id.col], tab)
                res.comb <- tab

            }



            ## #########################################
            ##   add id-mapping if not present already
            ##
            ## #########################################
            ##if(!is.null(global.results$id.map)){
            if( sum(c('id.concat', 'id.mapped', 'id.query') %in% colnames(res.comb)) < 3 ){
                ##cat('test\n')
                ##cat(res.comb$id[1:3], '\n')
                ##res.comb <- data.frame(res.comb, stringsAsFactors=F)
                ##cat(res.comb$id[1:3], '\n')
                res.comb <- left_join(res.comb, global.results$id.map, 'id')
                ##View(res.comb)
                ##rownames(res.comb) <- res.comb[, 'id']
                ##cat(res.comb$id[1:3], '\n')
                ##cat('test2\n')
            }

            ##save(res.comb, file='tmp.RData')
            ##rownames(res.comb) <- make.names(res.comb$id.concat)
            rownames(res.comb) <- make.names(res.comb$id)

            ## #########################################
            ## store the results
            global.results$data$output <- res.comb

            ##cat(dim(global.results$data$output))

            ## #####################################
            ## set some flags
            global.param$analysis.run <- T
            global.results$export.results <- F

            ## #################################################################
            ##            insert the panels for the volcanos
            ## #################################################################
            if(!(global.param$which.test %in% c('mod F', 'none'))){
                ins.volc()
            }

            ## #################################################################
            ##            insert panels for scatterplots
            ## #################################################################
            ins.scat()


            ## #################################################################
            ## increment counter: will invoke the filter
            global.param$run.test <- global.param$run.test + 1
        })


        ###################################################################################
        ##
        ##
        ##                   filter the test results accross multiple groups
        ##
        ##
        ###################################################################################
        ##filter.res  <-  reactive({
        observeEvent(c(input$filter.type, input$filter.value.top.n, input$filter.value.nom.p, input$filter.value.adj.p, global.param$run.test),{

            ## only run after one analysis has been completed
            if(global.param$run.test == 0 | is.null(input$filter.type)) return()

            cat('\n-- filter.res --\n')
            cat('filter.type: ', global.param$filter.type, '\nfilter.value:', global.param$filter.value, '\n')

            ## ######################################
            ## get the current tab
            tab.select <- input$mainPage

            groups.comp=unique(global.param$grp.comp)

            ## test results
            res <- data.frame(global.results$data$output, stringsAsFactors=F )
            ##View(data.frame(global.results$data$output))

            ## get the filter type
            global.param$filter.type=input$filter.type

            #################################
            ## top N
            if(global.param$filter.type=='top.n'){

                ######################################
                ## multiple comparisons
                if(length(groups.comp) > 1){

                    ###########################################
                    ## one/two sample tests
                    if(global.param$which.test != 'mod F'){

                        ## order according to minimal p-value
                        res <- res[ order( unlist(apply( res[, grep('^P.Value', colnames(res) )], 1, min))  ), ]
                        res <- res[ 1:input$filter.value.top.n, ]

                        ## order according to FC
                        res <- res[order( unlist(apply( res[, grep('^logFC', colnames(res) )], 1, min))  ), ]

                        ## now separate for each group comparison
                        res.groups <- lapply(groups.comp, function(g){
                            res.filt=res[ which(!is.na(paste('P.Value.', g, sep='') )), ]
                            res.filt=res.filt[ order( res.filt[, paste('P.Value.', g, sep='') ], decreasing=F) , ]
                            res.filt[1:input$filter.value.top.n, ]
                        })
                        names(res.groups) <- groups.comp

                    #########################################
                    ## F test
                    } else {

                        ## order according to
                        res <- res[order(res[, 'P.Value'], decreasing=F)[1:input$filter.value.top.n], ]

                        res.groups <- list(res)
                        names(res.groups) <- paste(groups.comp, collapse='|')
                    }

                    #################################################
                    ## one comparison only
                } else {
                     res <- res[order(res[, paste('P.Value.', groups.comp, sep='')]) , ]
                     res <- res[1:input$filter.value.top.n, ]
                     res <- res[order(res[, paste('logFC.', groups.comp, sep='')]) , ]

                     res.groups <- list(res)
                     names(res.groups) <- groups.comp
                }
                global.param$filter.value <- input$filter.value.top.n
            }
            #################################
            ## nominal p-value
            if(global.param$filter.type=='nom.p'){

                ###############################################
                ## multiple group comparisons
                if(length(groups.comp) > 1){

                    ###########################################
                    ## one/two sample tests
                    if(global.param$which.test != 'mod F'){

                        res <- res[ which( unlist(apply( res[, grep('^P.Value', colnames(res) )], 1, function(x) sum(x < input$filter.value.nom.p))) > 0), ]
                        res <- res[order( unlist(apply( res[, grep('^logFC', colnames(res) )], 1, min))  ), ]
                        ## now separate for each group comparison
                        res.groups <- lapply(groups.comp, function(g){
                            res[ which( res[, paste('P.Value.', g, sep='')] < input$filter.value.nom.p) , ]
                        })
                        names(res.groups) <- groups.comp
                    #########################################
                    ## F test
                    } else {
                        res <- res[which( res[, 'P.Value'] < input$filter.value.nom.p), ]

                        res.groups <- list(res)
                        names(res.groups) <- paste(groups.comp, collapse='|')
                    }
                ##############################################
                ## one group comparison
                } else {
                    res <- res[which(res[, paste('P.Value.', groups.comp, sep="")] < input$filter.value.nom.p), ]
                    res <-  res[order(res[, paste('logFC.', groups.comp, sep="")]), ]
                    res.groups <- list(res)
                    names(res.groups) <- groups.comp
                }
                global.param$filter.value <- input$filter.value.nom.p
            }

            #################################
            ## adjusted p-value
            if(global.param$filter.type=='adj.p'){

                ############################################################
                ## multiple group comparisons
                if(length(groups.comp) > 1){

                    ###########################################
                    ## one/two sample tests
                    if(global.param$which.test != 'mod F'){

                        res <- res[ which( unlist(apply( res[, grep('^adj.P.Val', colnames(res) )], 1, function(x) sum(as.numeric(x) < input$filter.value.adj.p ))) > 0), ]
                        res <- res[order( unlist(apply( res[, grep('^logFC', colnames(res) )], 1, min))  ), ]
                        ## now separate for each group comparison
                        res.groups <- lapply(groups.comp, function(g){
                            res[ which( res[, paste('adj.P.Val.', g, sep='')] < input$filter.value.adj.p) , ]
                        })
                        names(res.groups) <- groups.comp
                    ###########################################
                    ## F-test
                    } else {
                        res <- res[which( res[, 'adj.P.Val'] < input$filter.value.adj.p), ]

                        res.groups <- list(res)
                        names(res.groups) <- paste(groups.comp, collapse='|')
                    }
                ##############################################################
                ## one comparison
                } else {
                    res <- res[which(res[, paste('adj.P.Val.', groups.comp,sep='')] < input$filter.value.adj.p), ]
                    res <-  res[order(res[, paste('logFC.', groups.comp, sep="")]), ]
                    res.groups <- list(res)
                    names(res.groups) <- groups.comp
                }

                global.param$filter.value <- input$filter.value.adj.p
            }

            ###################################
            ## global FDR

            #################################
            ## no filter
            if(global.param$filter.type=='none'){
                global.param$filter.value <- 'none'
                res.groups <- lapply(groups.comp, function(g) res)
                names(res.groups) <- groups.comp
            }

            ###################################################
            ## global filter accross all experiments
            ##rownames(res) <- res[, 'id.concat']

            global.results$filtered <- res
            ##View(res)
            global.results$filtered.groups <- res.groups


            ## #####################################################
            ## suppress switching to the first tab
            updateNavbarPage(session, 'mainPage', selected=tab.select)

            ## trigger selectize update
            global.param$update.ppi.select <- TRUE
            global.param$update.ppi.select.scat <- TRUE

        })





        ##@################################################################################
        ##
        ##                                     output
        ##
        ##@################################################################################

        #############################################################
        ##
        ##                download zip file
        ## - append timestamp to zip-file
        ##
        #############################################################

        ## download handler for zip-file
        output$download.results <- downloadHandler(

	    filename = function(){paste('results', global.param$zip.name, sep='_')},
            content = function(file){
                file.copy( paste(global.param$session.dir, global.param$zip.name, sep='/'), file)
            }, contentType = "application/zip"
        )

        #########################################
        ## download handler for experimental design template
        output$exportTemplate <- downloadHandler(
            filename = function(){ 'experimentalDesign.txt' },
            content = function(file){
                tab <- global.input$table
                exp.design <- cbind(colnames(tab), rep('', ncol(tab)))
                colnames(exp.design) <- c('Column.Name', 'Experiment')
                write.table(  exp.design, file, sep='\t', quote=F, na='', row.names=F  )
            }
        )

        ###################################################################################
        ##
        ##                    summary dataset
        ##
        ###################################################################################
        output$summary.data <- renderTable({

            if(is.null(global.results$data)) return()
            if(!is.null(error$msg)) return()

            tab <- data.frame(global.input$table.org)
            ##View(tab)
            ##cat(dim(tab), '\n')
            grp <- global.param$grp
            N.grp <- global.param$N.grp

            ##cat(nrow(tab), ' - ', length(grp), ' - ', N.grp)

            sum.tab <- t(data.frame(
                ##id=c('Rows', 'Expression columns', 'Groups'),
                N.rows=nrow(tab),
                N.columns=length(grp),
                N.groups=N.grp
            ))
            ##View(sum.tab)
            sum.tab <- data.frame(id=c('Rows', 'Expression columns', 'Groups'), sum.tab)
            colnames(sum.tab) <- c('', 'Number')
            ##rownames(sum.tab) <- c('Rows', 'Expression columns', 'Groups')

            sum.tab
        })
        #####################################################################################
        ##
        ##                summary workflow
        ##
        #####################################################################################
        output$summary.workflow <- renderTable({

            if(is.null(global.results$data)) return()
            if(!is.null(error$msg)) return()

            wf.tab.ids <- c('Log scale', 'Normalization', 'Filter data', 'Test', 'Filter results')

            ## #########################
            ## Reproducibility filter
            if(global.param$filt.data == 'Reproducibility'){

                wf.tab <- t(data.frame( global.param$log.transform, global.param$norm.data, paste( global.param$filt.data, ' (alpha=',global.param$repro.filt.val, ')',sep=''), global.param$which.test,
                                       paste( global.param$filter.type, ' < ', global.param$filter.value)))

                wf.tab <- data.frame(id=wf.tab.ids, wf.tab, stringsAsFactors=F)

            }

            ## #######################
            ## SD filter
            if(global.param$filt.data == 'StdDev'){
                wf.tab <- t(data.frame( global.param$log.transform, global.param$norm.data, paste( global.param$filt.data, ' (SD=',global.param$sd.filt.val, '%)',sep=''), global.param$which.test,
                                       paste( global.param$filter.type, ' < ', global.param$filter.value) ))
                wf.tab <- data.frame(id=wf.tab.ids, wf.tab, stringsAsFactors=F)
            }

            ## #############################
            ## no data filter
            if(global.param$filt.data == 'none'){
                wf.tab <- t(data.frame( global.param$log.transform, global.param$norm.data, paste( global.param$filt.data, sep=''), global.param$which.test,
                                       paste( global.param$filter.type, ' < ', global.param$filter.value) ))
                wf.tab <- data.frame(id=wf.tab.ids, wf.tab, stringsAsFactors=F)
            }

            ## special case: no filter
            if(global.param$filter.type == 'none')
                wf.tab[which(wf.tab$id == 'Filter results'), 2] <- 'none'

            ## special case: top N
            if(global.param$filter.type == 'top.n')
                wf.tab[which(wf.tab$id == 'Filter results'), 2] <-  paste('top', global.param$filter.value)

            ## suppress column names
            colnames(wf.tab) <- c('', '')

            ## insert
            wf.tab

        })
        ###################################################################################
        ##
        ##                summary test results
        ##
        ###################################################################################
        output$summary.test <- renderTable({

            if(is.null(global.results$data)) return()
            if(!is.null(error$msg)) return()

            ## extract results ## kk 20161025
            ##tab.select <- input$mainPage
            ###filter.res()#--
            ##updateNavbarPage(session, 'mainPage', selected=tab.select)

            ##filter.res()

            ##res = global.results$filtered
            res <- global.results$data$output

            ## tested groups
            grp.comp=unique(global.param$grp.comp)

            ## extract filter type
            filter.type=global.param$filter.type
            filter.value=global.param$filter.value

            ######################################
            ## one/two sample mod T
            if(global.param$which.test != 'mod F'){

                if(filter.type == 'adj.p')
                    test.tab=unlist(lapply(paste('adj.P.Val', grp.comp, sep='.'), function(x) sum(res[, x] < filter.value, na.rm=T) ))
                if(filter.type == 'nom.p')
                    test.tab=unlist(lapply(paste('P.Value', grp.comp, sep='.'), function(x) sum(res[, x] < filter.value, na.rm=T) ))
                if(filter.type == 'top.n')
                    return(NULL)

                if(filter.type == 'none')
                    return(NULL)

                test.tab=data.frame(test.tab)

                ##rownames(test.tab) <- chopString(as.character(make.names(make.unique(grp.comp))), nChar=20)
                test.tab <- data.frame(id=chopString(as.character(make.names(make.unique(grp.comp))), nChar=20), test.tab)
                colnames(test.tab) <- c('','Number significant')
                return(test.tab)

            #######################################
            ## moderated F
            } else {
                if(filter.type == 'adj.p'){
                    test.tab=data.frame(  sum(res[, 'adj.P.Val'] < filter.value, na.rm=T) )
                    ##sum.na <- sum(is.na(res[, 'adj.P.Val']))
                    ##if( sum.na> 0)
                    ##    test.tab=paste(test.tab, '(', ,')')
                }
                if(filter.type == 'nom.p'){
                    test.tab=data.frame(  sum(res[, 'P.Value'] < filter.value, na.rm=T) )
                }
                if(filter.type == 'top.n')
                    return(NULL)
                if(filter.type == 'none')
                    return(NULL)

                test.tab <- data.frame(test.tab)

                test.tab <- data.frame(id=chopString( paste(unique(global.param$grp), collapse=' vs. '), 20), test.tab)
                colnames(test.tab) <- c('','Number significant')
                return(test.tab)
            }

        })
        ##@###############################################################
        ##
        ##          missing data distribution
        ##
        ##@###############################################################

        ##@################################
        ## per row
        output$summary.missing.data.row <- renderPlotly({

            if(is.null(global.results$data)) return()
            ##if(!is.null(error$msg)) return()

            if(global.param$log.transform == 'none')
                tab <- data.frame(global.input$table.org)
            else
                tab <- data.frame(global.results$table.log)

            #grp <- global.param$grp
            #N.grp <- global.param$N.grp
            #grp.colors.legend <- global.param$grp.colors.legend

            ## extract expression values
            #dat <- tab[, -which(colnames(tab) == global.param$id.col.value)]
            dat <- tab[, names(global.param$grp)]
            dat <- data.matrix(dat)

            dat[is.infinite(dat)] <- NA

            ## number of missing values per row
            na.row.idx <- table(apply(dat, 1, function(x) sum(is.na(x))))

            
            #n.miss <- seq(as.numeric(names(na.row.idx)[1]):as.numeric(names(na.row.idx)[length(na.row.idx)]) )
            n.miss <- rep(0, ncol(dat)+1)
            names(n.miss) <- 0:ncol(dat)
            n.miss[names(na.row.idx)] <- na.row.idx
            n.miss <- cumsum(n.miss)
            
            p <- plot_ly( x=as.numeric(names(n.miss)), y=n.miss, type='scatter', mode='lines+markers', marker=list(color = 'black'), line=list(color='black') )
            p <- layout(p, title=paste('Fully quantified features:', n.miss[1]), xaxis=list(title=paste('# missing values')), yaxis=list(title=paste('# quantified features')))
            
           # save(na.row.idx, dat, n.miss, file='tmp.RData')
            
            #p <- plot_ly( x=names(na.row.idx)[2:length(na.row.idx)], y=na.row.idx[2:length(na.row.idx)], type='bar' )
            #p <- layout(p, title=paste('Fully quantified features:', na.row.idx[1]), xaxis=list(title=paste('# missing values')), yaxis=list(title=paste('# data rows')))
            p
        })
        ## #############################################################################
        ##
        ## per column: Non-missing values
        ##
        ## #############################################################################
        output$summary.nonmissing.data <- renderPlotly({
            if(is.null(global.results$data)) return()


            if(global.param$log.transform == 'none')
                tab <- data.frame(global.input$table.org)
            else
                tab <- data.frame(global.results$table.log)

            grp <- global.param$grp
            N.grp <- global.param$N.grp
            grp.colors.legend <- global.param$grp.colors.legend
            grp.colors <- global.param$grp.colors

            ## extract expression values
            dat <- tab[, -which(colnames(tab) == global.param$id.col.value)]
            dat <- data.matrix(dat)

            ## ############################################################
            ## order columns
            ord.idx <- order(grp)
            dat <- dat[, ord.idx]

            grp.colors <- grp.colors[ord.idx]
            grp <- grp[ord.idx]

            grp.colors.legend <- grp.colors.legend[order(names(grp.colors.legend))]
            names(grp.colors) <- colnames(dat)


            ## #############################################################
            ## number of non-missing values per row
            na.col.idx <- apply(dat, 2, function(x) sum(!is.na(x)))

            dat.plot <- data.frame(x=names(na.col.idx), y=na.col.idx)
            dat.plot$x <- factor(dat.plot$x, levels=dat.plot[['x']])

            ## plot
            ##p <- plot_ly( x=names(na.col.idx), y=na.col.idx,  color=grp, colors=grp.colors.legend, type='bar')
            p <- plot_ly(dat.plot, x=~x, y=~y, color=grp, colors=grp.colors.legend, type='bar')
            p <- layout(p, title='Number of valid data points per data column', xaxis=list(title=paste('Data columns')), yaxis=list(title=paste('# data points')))

            ##################################################
            ##global.param$session.import.init <- F
            p
        })


        ##@################################################################################
        ##
        ##             display the corresponding part of the table
        ##
        ##@################################################################################
        output$tableprev <- DT::renderDataTable({

            if(is.null(global.results$data)) return()
            if(!is.null(error$msg)) return()

  
            tab <- global.results$filtered
            colnames(tab) <- sub('^X','',colnames(tab))
            ##rownames(tab) <- tab[, input$id.col.value]

            ## append annotation columns
            table.anno <- global.input$table.anno


            ## check whether there is annotation stored
            if(!is.null(table.anno)){
                if(is.null(dim(table.anno)))
                    table.anno <- data.frame(table.anno)
               ## View(tab)
               ## View(table.anno)
                ##tab <- cbind(tab, table.anno[ rownames(tab), ])
                tab <- left_join(tab, table.anno, 'id')
            }

            if(nrow(tab) > 0){
                ## add links to GeneCard
                up.id <- tab$id
                up.link <- link.db(up.id, global.results$keytype)
                tab[, 'id'] <- up.link
            }
            tab

        }, options = list( pageLength = 20, scrollX = T), escape=F, filter='top', rownames=F)


        ## ################################################################################
        ##
        ##                             Volcano plot
        ##
        ## ################################################################################


        ## ##################################################################
        ## obserever to trigger  'updateSelectizeInput' for volcano plots
        ## - server side  rendering of selectizeInput
        ## ##################################################################
        observe({
            if( !global.param$update.ppi.select ) return()
            ##cat('now:', global.param$update.ppi.select, '\n')

            grp.comp <- unique( global.param$grp.comp )

            for(i in 1:length(grp.comp)){

                ## commenting out the if-statemnet made select input work after multiple rounds of analysis
                ##if(!is.null(input[[ gsub('\\.','', paste0('ppi.bait.',  grp.comp[i])) ]] )){

                ## #################################
                ## all ids found in data
                choices=unlist(global.results$id.map$id.concat)


                ## server-side rendering of selectizeInput
                updateSelectizeInput( session=session, inputId = gsub('\\.','',paste0('ppi.bait.',  grp.comp[i])), label='Bait protein', choices=choices, selected=NULL, server=T)

                ##}
            }
            global.param$update.ppi.select <- FALSE
        })
        ## ##################################################################
        ## obserever to trigger  'updateSelectizeInput' for SCATTER plots
        ## - server side rendering of selectizeInput
        ## ##################################################################
        observe({
            if( !global.param$update.ppi.select.scat ) return()

            ##cat('now:', global.param$update.ppi.select.scat, '\n')

            grp.comp <- unique( global.param$grp )

            for(i in 1:length(grp.comp)){

                ## #################################
                ## all ids found in data
                choices=unlist(global.results$id.map$id.concat)

                ## server-side rendering of selectizeInput
                updateSelectizeInput( session=session, inputId = gsub('\\.','',paste0('ppi.bait.scat.',  grp.comp[i])), label='Bait protein', choices=choices, selected=NULL, server=T)
            }
            global.param$update.ppi.select.scat <- FALSE
        })



        ## ##############################################################
        ##
        ##            insert scatterplots
        ##
        ## ##############################################################
        ins.scat <- reactive({

            withProgress({
                setProgress( message='Preparing scatterplots...')

                ## volcanos for each group comparison
                ##if(global.param$which.test == 'Two-sample mod T')
                ##    grp.comp <- unique( global.param$grp.comp )
                ##else
                    grp.comp <- unique( global.param$grp )

                ## #########################################
                ## loop over group comparsions
                for(i in 1:length(grp.comp)){
                    local({

                        my_i <- i
                        ## ########################
                        ## the actual plots
                        output[[ paste("scatterplot", grp.comp[my_i], sep='.') ]] <- renderPlotly({
                            plotScatter( grp.comp[my_i] )
                        })
                    })
                }})
            global.param$update.ppi.select.scat <- TRUE
        }) ## end reactive


        ## #############################################################
        ##
        ##            insert volcano panels
        ##
        ## #############################################################
        ins.volc <- reactive({

            if(global.param$which.test %in% c('mod F', 'none')) return()

            withProgress({
                setProgress( message='Preparing volcanos...')

            ## volcanos for each group comparison
            grp.comp <- unique( global.param$grp.comp )

            ## #########################################
            ## loop over group comparsions
            for(i in 1:length(grp.comp)){


                local({

                    my_i <- i
                    ## ########################
                    ## the actual plots
                    output[[ paste("volcano", grp.comp[my_i], sep='.') ]] <- renderPlot({
                        plotVolcano( grp.comp[my_i] )
                    })

                  ## ###############################
                  ## observe brush
                  observeEvent(input[[paste('plot_brush', grp.comp[my_i], sep='.')]], {

                      tmp <- input[[paste('plot_brush', grp.comp[my_i], sep='.')]]

                      volc.brush[[paste('xmin', grp.comp[my_i], sep='.')]] <- tmp$xmin
                      volc.brush[[paste('xmax', grp.comp[my_i], sep='.')]] <- tmp$xmax
                      volc.brush[[paste('ymax', grp.comp[my_i], sep='.')]] <- tmp$ymax
                      volc.brush[[paste('ymin', grp.comp[my_i], sep='.')]] <- tmp$ymin

                  })

                  ## ###############################
                  ## observe double click
                  observeEvent(input[[paste('plot_dblclick', grp.comp[my_i], sep='.')]], {
                      volc.brush[[paste('xmin', grp.comp[my_i], sep='.')]] <- NULL
                      volc.brush[[paste('xmax', grp.comp[my_i], sep='.')]] <- NULL
                      volc.brush[[paste('ymax', grp.comp[my_i], sep='.')]] <- NULL
                      volc.brush[[paste('ymin', grp.comp[my_i], sep='.')]] <- NULL
                  })


                  ## ################################
                  ## observe clicks
                  observeEvent( input[[paste('plot_click', grp.comp[my_i], sep='.')]], {

                      ## the data set
                      res <- as.data.frame( global.results$data$output )

                      ## determine what to show in the plot, i.e. 'id' or mapped gene names
                      if( !is.null(global.results$id.map ))
                          txt.col <- 'id.concat'
                      else
                          txt.col <- 'id'

                      ## identify the points clicked
                      text.tmp <- nearPoints(res, input[[paste('plot_click', grp.comp[my_i], sep='.')]], threshold=10, maxpoints =  1, xvar=paste('logFC', grp.comp[my_i], sep='.'), yvar=paste('Log.P.Value', grp.comp[my_i], sep='.'))

                      ## if there are any
                      if(nrow(text.tmp) == 1){

                          ################################################
                          ## first click
                          if(is.null(volc[[ paste('x', grp.comp[my_i], sep='.')]] )){
                              volc[[paste('x', grp.comp[my_i], sep='.')]] = text.tmp[paste('logFC', grp.comp[my_i], sep='.')]
                              volc[[paste('y', grp.comp[my_i], sep='.')]] = text.tmp[paste('Log.P.Value', grp.comp[my_i], sep='.')]
                              ## volc[[paste('text', grp.comp[my_i], sep='.')]] = text.tmp['id']
                               volc[[paste('text', grp.comp[my_i], sep='.')]] = text.tmp[ txt.col ]
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
                                  ##volc[[paste('text', grp.comp[my_i], sep='.')]]=c(volc[[paste('text', grp.comp[my_i], sep='.')]],  text.tmp[ 'id'] )
                                  volc[[paste('text', grp.comp[my_i], sep='.')]]=c(volc[[paste('text', grp.comp[my_i], sep='.')]],  text.tmp[ txt.col ] )
                                  volc[[paste('xy', grp.comp[my_i], sep='.')]] = c(volc[[paste('xy', grp.comp[my_i], sep='.')]], paste(text.tmp[paste('logFC', grp.comp[my_i], sep='.')], text.tmp[paste('Log.P.Value', grp.comp[my_i], sep='.')]) )

                                  volc[[paste('P.Value', grp.comp[my_i], sep='.')]]=c(volc[[paste('P.Value', grp.comp[my_i], sep='.')]],  text.tmp[paste('P.Value', grp.comp[my_i], sep='.')] )
                                  volc[[paste('adj.P.Val', grp.comp[my_i], sep='.')]]=c(volc[[paste('adj.P.Val', grp.comp[my_i], sep='.')]],  text.tmp[paste('adj.P.Val', grp.comp[my_i], sep='.')] )
                              }
                          }
                      }
                  }) ## end observe clicks

                  ## ###################################################
                  ## reset volcano annotations
                  observeEvent(input[[paste('volc.tab.reset', grp.comp[my_i], sep='.')]],{
                      volc[[paste('x', grp.comp[my_i], sep='.')]] <- NULL
                      volc[[paste('y', grp.comp[my_i], sep='.')]] <- NULL
                      volc[[paste('xy', grp.comp[my_i], sep='.')]] <- NULL
                      volc[[paste('P.Value', grp.comp[my_i], sep='.')]] <- NULL
                      volc[[paste('adj.P.Val', grp.comp[my_i], sep='.')]] <- NULL
                  })

                  ## ####################################################
                  ## table of selected features
                  output[[paste('volc.tab.selected', grp.comp[my_i], sep='.')]] <- renderTable({

                      if(is.null(volc[[paste('x', grp.comp[my_i], sep='.')]])) return()
                      if(length(volc[[paste('x', grp.comp[my_i], sep='.')]]) == 0) return()

                      tags$h4('Selection:')

                      id.tmp <- volc[[paste('text', grp.comp[my_i], sep='.')]]

                      dat.select = data.frame(id=unlist(volc[[paste('text', grp.comp[my_i], sep='.')]]), logFC=unlist(volc[[paste('x', grp.comp[my_i], sep='.')]]), P.Value=unlist(volc[[paste('P.Value', grp.comp[my_i], sep='.')]]), adj.P.Value=unlist(volc[[paste('adj.P.Val', grp.comp[my_i], sep='.')]]) )
                      up.id <- dat.select[, 'id']
                      #up.link <- paste("<a href='http://www.uniprot.org/uniprot/", sub('(_|,|;|\\.).*', '', up.id),"' target='_blank'>", up.id, "</a>", sep='')
                      up.link <- link.db(up.id, global.results$keytype) 
                      dat.select[, 'id'] <- up.link

                      dat.select
                  }, sanitize.text.function = function(x) x)

              }) ## end local

                ##incProgress( 1/length(grp.comp))

            } ## end for loop

        }) ## end withProgress

            ## trigger selectize update
            global.param$update.ppi.select <- TRUE

        }) ## end 'ins.volc'


        ## #############################################################
        ##
        ##                   SCATTERPLOTS
        ##  - actual plot is generated here
        ##
        ## #############################################################
        plotScatter <- function(group){

            cat('\n-- plotScatter --\n')
            if(!is.null(error$msg)) return()

            ## #############################
            ## unfiltered data set
            res = as.data.frame( global.results$data$output )

            ## #############################
            ## p-values one-sample T
            if(global.param$which.test == 'One-sample mod T'){

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
            ## p-values two-sample
            if(global.param$which.test == 'Two-sample mod T') {
                pval <- rep(1, nrow(res))
            }
            ## ###############################
            ## p-values no test
            if(global.param$which.test == 'none') {
                pval <- rep(1, nrow(res))
            }

            ## ##############################
            ## extract data columns
            ##cat(paste0('scat.x.', group), '\n')
            ##cat(paste0('scat.y.', group), '\n')
            x.ax <- res[ , input[[ paste0('scat.x.', group) ]] ]
            y.ax <- res[ , input[[ paste0('scat.y.', group) ]] ]

            ## ##############################
            ## ids
            #IDs <- global.results$id.map$id.concat
            IDs <- res[ , 'id.concat'] 
            names(x.ax) <- names(y.ax) <- names(pval) <- IDs

            ## index of missing values
            rm.idx <- union( which(is.na(x.ax)), which(is.na(y.ax)) )
            if(length(rm.idx) > 0){
                res <- res[-rm.idx, ]
                x.ax <- x.ax[-rm.idx]
                y.ax <- y.ax[-rm.idx]
                IDs <- IDs[-rm.idx]
                pval <- pval[-rm.idx]
            }

            ## ###############################
            ## color significant values
            col <- rep('black', nrow(res))
            sig.idx <- c()
            if(global.param$filter.type == 'adj.p' | global.param$filter.type == 'nom.p'){
                sig.idx <- which(pval < as.numeric(global.param$filter.value ))
                col[ sig.idx ] <- 'red'
            }

            ## ################################
            ## data frame
            dat.plot <- data.frame(x.ax, y.ax, pval, IDs, stringsAsFactors=F)

            ## check whether to draw PPI stuff
            ppi.bait <- input[[ gsub('\\.', '', paste0('ppi.bait.scat', group)) ]]
            PPI <- FALSE

            if( nchar(ppi.bait) > 0 & toupper(ppi.bait) %in% toupper(IDs))
                PPI <- TRUE

            ## ##############################
            ##        set up the  plot
            ## ##############################
            if(PPI)
                p <- plot_ly(x=dat.plot$x.ax, y=dat.plot$y.ax, type='scatter', mode='markers', marker=list(size=10, color=col), text=IDs, opacity=.2, showlegend=F, name='non-interactors' )
                ##p <- plot_ly(x=dat.plot$x.ax, y=dat.plot$y.ax, type='scatter', mode='markers', marker=list(size=10, color=col, showlegend=FALSE), text=IDs, alpha=.1 )
            else
                p <- plot_ly(x=dat.plot$x.ax, y=dat.plot$y.ax, type='scatter', mode='markers', marker=list(size=10, color=col), text=IDs, showlegend=F )
            ## ########################################################
            ##                  add PPI stuff
            ## ##########################################################
            ppi.idx <- c()
            ##if(toupper(ppi.bait) %in% toupper(IDs)) {
            if(PPI){
                ## ##############################################
                ##
                ##            extract interactors
                ##
                ## ###############################################
                ppi.map <- get.interactors( ppi.bait=ppi.bait,
                                           IDs=IDs,
                                           sig.idx=sig.idx,
                                           db=input[[ paste0('ppi.db.scat.', group) ]],
                                           ppi.db=ppi,
                                           ppi.db.col=ppi.db.col
                                           )

                ## ################################
                ## extract results
                ppi.int.idx <- ppi.map$ppi.int.idx
                ppi.bait.idx <- ppi.map$ppi.bait.idx
                leg <- ppi.map$leg
                leg.col <- ppi.map$leg.col
                ppi.col <- ppi.map$ppi.col

                ## ###############################
                ## add interactors
                if(length(ppi.int.idx) > 0){
                    col[ppi.int.idx] <- ppi.col[ nchar(ppi.col) > 0 ]
                    ppi.idx <- ppi.int.idx

                    p <- p %>% add_trace(x=dat.plot$x.ax[ppi.int.idx], y=dat.plot$y.ax[ppi.int.idx], type='scatter', mode = 'markers', marker=list( size=10, color=col[ppi.int.idx]), text=IDs[ppi.int.idx], name='Interactors', opacity=1 )
                }
                ## ################################
                ## add bait
                if(length(ppi.bait.idx) > 0){
                    col[ppi.bait.idx] <- 'green'
                    ppi.idx <- c(ppi.idx, ppi.int.idx)

                    p <- p %>% add_trace(x=dat.plot$x.ax[ppi.bait.idx], y=dat.plot$y.ax[ppi.bait.idx], type='scatter', mode = 'markers', marker=list( size=15, color=col[ppi.bait.idx]), text=IDs[ppi.bait.idx], name=paste('Bait protein:', ppi.bait), opacity=1 )
                }
            }

            ## some astethics
            p <- layout(p,
                        title='',
                        xaxis=list(
                            title=paste( input[[ paste0('scat.x.', group) ]])
                        ),
                        yaxis=list(
                            title=paste( input[[ paste0('scat.y.', group) ]] )
                        )
                        )

            return(p)
        }

        ################################################################
        ##
        ##                 VOLCANO PLOTS
        ##
        ##  - actual plot is generetad here
        ##
        ################################################################
        plotVolcano <- function(group, interactors=NULL){

            cat('\n-- plotVolcano --\n')
            if(!is.null(error$msg)) return()

            ## hyperbolic curve filter?
            hyperbol <- input[[paste('ppi.hyper.filt', group, sep='.' )]]

            ## maximal log10 p-value
            max.logP <- input[[paste('max.logP', group, sep='.')]]

            ## pch for significant points
            sig.pch=23

            ## vectors for selected points/PP interactors
            volc.add.X <- volc.add.Y <- volc.add.text <- volc.add.col <- c()

            ## unfiltered data set
            res = as.data.frame( global.results$data$output )

            ## #############################
            ## p-values
            if(global.param$which.test != 'mod F'){
                logPVal <- res[, paste('Log.P.Value.', group, sep='')]
            } else {
                logPVal <- res[, paste('Log.P.Value',  sep='')]
            }

            # ###############################
            # adjusted p-values
            adjPVal <- res[, paste('adj.P.Val.', group, sep='')]
            PVal <- res[, paste('P.Value.', group, sep='')]
            
            ## ##############################
            ## log fold change
            logFC <- res[, paste('logFC.', group, sep='')]

            ## ##############################
            ## ids
            #IDs <- global.results$id.map$id.concat
            IDs <- res[ , 'id.concat']
            #View(head(res))
            #cat(IDs[1:4])
            
            ## ###################################################
            ##             use IDs as vector names
            names(logPVal) <- names(adjPVal) <- names(logFC) <- names(PVal) <- IDs

            ## index of missing values
            rm.idx <- union( which(is.na(logFC)), which(is.na(logPVal)) )

            if(length(rm.idx) > 0){
                res <- res[-rm.idx, ]
                logFC <- logFC[-rm.idx]
                logPVal <- logPVal[-rm.idx]
                PVal <- PVal[-rm.idx]
                IDs <- IDs[-rm.idx]
                adjPVal <- adjPVal[ -rm.idx ]
            }
            ## store a copy of IDs before zoom
            IDs.all <- IDs

            ## which filter has been used?
            filter.str <- paste('filter:', global.param$filter.type, '\ncutoff:', global.param$filter.value)

            ## ###################################################################
            ##                 set maximal log p value
            ## ###################################################################
            if(!is.null( max.logP))
                logPVal[which(logPVal > max.logP)] <- max.logP

            ## ###################################################################
            ##                         limits
            ## ###################################################################
            xlim = max(abs(logFC), na.rm=T)
            xlim = xlim + xlim*.1

            ylim = ifelse(is.null(max.logP), max(logPVal, na.rm=T), max.logP)
            ylim = ylim + .2*ylim

            ## ###################################################################
            ##                    hyperbolic curve
            ## ###################################################################
            if( hyperbol ){

                x0 <- as.numeric(input[[paste( "ppi.min.fc", group, sep='.')]])
                c <- as.numeric(input[[paste( "ppi.curve", group, sep='.')]])

                y.hc=function(x, x0, y0, c) return(c/(x-x0) + y0)
                x.hc <- seq(x0, xlim, 0.1)
            }

            ######################################################################
            ## extract significant proteins of current group/test
            ######################################################################
            ## one/two sample
            if( global.param$which.test != 'mod F'){
                if(global.param$filter.type == 'top.n'){
                    #PVal <- res[, paste('P.Value.', group, sep='')]
                    sig.idx = order(PVal, decreasing=F)[1:global.param$filter.value]

                }
                if(global.param$filter.type == 'nom.p'){
                    #PVal <- res[, paste('P.Value.', group, sep='')]
                    sig.idx = which(PVal < global.param$filter.value)
                }
                ## ##################################
                ## adjusted p
                if(global.param$filter.type == 'adj.p'){
                  
                    ## hyperbol
                    if(hyperbol){
                        #adjPVal <- res[, paste('adj.P.Val.', group, sep='')]

                        sig.idx = which(adjPVal < global.param$filter.value)
                        #names(sig.idx) <- IDs[sig.idx]

                        y0 <- min(logPVal[names(sig.idx)], na.rm=T)

                        sig.idx <- which( logPVal > y.hc( abs(logFC), x0=x0, y0=y0, c=c) & abs(logFC) > x0 )

                    ## adjusted p only
                    } else {
                        #adjPVal <- res[, paste('adj.P.Val.', group, sep='')]
                        sig.idx = which(adjPVal < global.param$filter.value)
                        #names(sig.idx) <- IDs[sig.idx]
                        
                    }
                }
            ######################################
            ## F-test
           ## } else {
           ##     if(global.param$filter.type == 'top.n'){
           ##         PVal <- res[, paste('P.Value', sep='')]
           ##         sig.idx = order(PVal, decreasing=F)[1:global.param$filter.value]
           ##     }
           ##     if(global.param$filter.type == 'nom.p'){
           ##         PVal <- res[, paste('P.Value', sep='')]
           ##         sig.idx = which(PVal < global.param$filter.value)
           ##     }
           ##     if(global.param$filter.type == 'adj.p'){
           ##         adjPVal <- res[, paste('adj.P.Val', sep='')]
           ##         sig.idx = which(adjPVal < global.param$filter.value)
           ##     }
            }

            if(global.param$filter.type == 'none')
                sig.idx <- c()
                ##sig.idx = 1:length(logFC)

            ## use IDs as names
            if(length(sig.idx) >0 )
                names(sig.idx) <- IDs[sig.idx]

            sig.idx.all <- sig.idx

            ## ####################################################
            ## PPI?
            ppi.bait <- input[[ gsub('\\.', '', paste0('ppi.bait.', group)) ]]
            PPI <- toupper(ppi.bait) %in% toupper(IDs.all)
            
            ## #################################################
            ##
            ##              pch and cex
            ##
            ## ##################################################
            pch.vec=rep(21, nrow(res))
            cex.vec=rep( input[[paste('cex.volcano', group, sep='.')]], nrow(res))

            if(length(sig.idx) > 0){
                pch.vec[ sig.idx ] <- sig.pch
                #cex.vec[ sig.idx ] <- cex.vec[1]+.5
            }

            #############################
            ## color gradient
            if(PPI){
              opac1 = 0.3
              opac2 = 0.2
            } else{
              opac1 = 1
              opac2 = 0.2
            }

            col=myColorRamp(c('grey30', 'grey50', 'darkred', 'red', 'deeppink'), na.omit(logPVal), range=c(0, max.logP), opac=opac1)
            col.opac=myColorRamp(c('grey30', 'grey50', 'darkred', 'red', 'deeppink'), na.omit(logPVal), range=c(0, max.logP), opac=opac2)



            ## ########################################################
            ##
            ##                Query PPI databases
            ##
            ## ##########################################################
           
            ## ##############################################
            ##
            ## if the bait was detected
            ##
            ## ###############################################
            
            if(PPI) {
                # extract interactors
                ppi.map <- get.interactors( ppi.bait=ppi.bait,
                                           IDs=IDs.all,
                                           sig.idx=sig.idx.all,
                                           db=input[[ paste0('ppi.db.', group) ]],
                                           ppi.db=ppi,
                                           ppi.db.col=ppi.db.col
                                           )

                ## ################################
                ## extract results
                ppi.int.idx <- ppi.map$ppi.int.idx

                ##ppi.bait.idx <- ppi.map$ppi.bait.idx
                leg <- ppi.map$leg
                leg.col <- ppi.map$leg.col
                ppi.col <- ppi.map$ppi.col

                #ppi.bait.idx <- which(toupper(IDs) == toupper(ppi.bait))
                ppi.bait.idx <- ppi.map$ppi.bait.idx
                ppi.col[ ppi.bait.idx ] <- 'green'
                
                cex.vec[ ppi.bait.idx ] <- cex.vec[ ppi.bait.idx ] + 1
                cex.vec[ ppi.int.idx ] <- cex.vec[ ppi.int.idx ] + 1
                
                ppi.int.vec <- rep(FALSE, length(IDs))
                ppi.int.vec[ppi.int.idx] <- TRUE

                ## udpate color vector
                col[ nchar(ppi.col) > 0] <- ppi.col[ nchar(ppi.col) > 0]
                col.opac[ nchar(ppi.col) > 0] <- ppi.col[ nchar(ppi.col) > 0 ]
                
                

                ## ###############################
                ## plot if there are interactors
                if(length(ppi.int.idx) > 0) {
                    ## labels
                    if(input[[paste('ppi.show.labels', group, sep='.')]]){

                        volc.add.X <- c(volc.add.X, logFC[ppi.int.idx])
                        volc.add.Y <- c(volc.add.Y, logPVal[ppi.int.idx])
                        volc.add.text <- c( volc.add.text, as.character(IDs[ppi.int.idx]))
                        volc.add.col <- c(volc.add.col, col[ppi.int.idx])
                    }
                }
            }

            ## ##################################################################
            ##
            ##                   zoomed vs not zoomed
            ##
            ## ##################################################################
            if( is.null( volc.brush[[ paste('xmin', group, sep='.') ]] ) ){
                xlim = c(-xlim, xlim)
                ## y-limits
                ylim = c(0, ylim)

            } else {

                ## ##################################
                ##        zoomed
                ## x-axis
                xlim = c(volc.brush[[ paste('xmin', group, sep='.') ]], volc.brush[[ paste('xmax', group, sep='.') ]])
                logFC[ logFC < xlim[1] | logFC > xlim[2] ] <-  NA

                ## y-axis
                ylim = c( max(0, volc.brush[[ paste('ymin', group, sep='.') ]]), volc.brush[[ paste('ymax', group, sep='.') ]])
                logPVal[ logPVal < ylim[1] | logPVal > ylim[2] ] <-  NA

                ## remove missing  values
                rm.idx <- union( which(is.na(logFC)), which(is.na(logPVal)) )

                if(length(rm.idx) > 0){
                    res <- res[-rm.idx, ]
                    logFC <- logFC[-rm.idx]
                    logPVal <- logPVal[-rm.idx]
                    IDs <- IDs[-rm.idx]

                    ## colors, pch, cex
                    col <- col[-rm.idx]
                    col.opac <- col.opac[-rm.idx]
                    pch.vec <- pch.vec[-rm.idx]
                    cex.vec <- cex.vec[-rm.idx]

                    ## update PPI stuff
                    if(PPI){
                        ppi.int.vec <- ppi.int.vec[-rm.idx]
                        ppi.int.idx <- which(ppi.int.vec)
                        ppi.bait.idx <- which(toupper(IDs) == toupper(ppi.bait))
                    }

                    ## update index of significant stuff
                    sig.idx <- sig.idx[ !(sig.idx %in% rm.idx) ]
                    ##sig.idx <- sig.idx[ -rm.idx ]
                }
            }


            ## ##################################################
            ##               set up the plot
            ## ##################################################
            par(mar=c(4,5,5,2))
            plot.new()
            plot.window( xlim=xlim, ylim=ylim, cex.axis=1.8, cex.lab=1.8)
            ## title
            mtext(group, side=3, cex=2, line=2)
            ## label axes
            if(global.param$which.test == 'Two-sample mod T')
                mtext( paste("log(", sub('.*\\.vs\\.', '', group), "/", sub('\\.vs.*', '', group),")"), side=1, cex=1.8, line=3)
            else
                mtext(expression(log(FC)), side=1, cex=1.8, line=3)
            mtext(expression(-10*log[10](P-value)), side=2, cex=1.8, line=3)
            ## draw axes
            axis(1, cex.axis=1.8)
            axis(2, las=2, cex.axis=1.8)
            ## grid
            if( input[[paste('grid.volcano', group, sep='.')]] )
                grid()

            ## actual plot
            points(logFC, logPVal, col=col, bg=col.opac, pch=pch.vec, cex=cex.vec, lwd=2)


            ## #######################################
            ## hyperbolic curve
            if( hyperbol & global.param$filter.type == 'adj.p'){
                lines( x.hc, y.hc(x.hc, x0, y0, c), col='grey30', lty='dashed')
                lines(-x.hc, y.hc(x.hc, x0, y0, c), col='grey30', lty='dashed')
                text( xlim[2]-(xlim[2]*.05), y0, paste(global.param$filter.type, global.param$filter.value, sep='='), pos=3, col='grey30')
            }

            ## ###################################
            ## add filter
            ## minimal log P-value for given filter
            if(length(sig.idx) > 0 & !(input[[paste('ppi.hyper.filt', group, sep='.' )]])){
                
              filt.minlogPVal <- min(logPVal[names(sig.idx)], na.rm=T)
            
                abline(h=filt.minlogPVal, col=my.col2rgb('grey30', 50), lwd=2, lty='dashed')
                text( xlim[2]-(xlim[2]*.05), filt.minlogPVal, paste(global.param$filter.type, global.param$filter.value, sep='='), pos=3, col='grey30')
            }

            ## number of significant
            legend(ifelse(PPI, 'topright', 'top'), bty='n', legend=paste(filter.str, '\nsig / tot: ', length(sig.idx),' / ', sum(!is.na(logFC) & !is.na(logPVal)), sep=''), cex=1.5)

            ## ############################
            ## indicate directionality for two-sample tests
            if(global.param$which.test == 'Two-sample mod T'){
                #legend('topleft', legend=sub('\\.vs.*', '', group), cex=2, text.col='darkblue', bty='n')
                #legend('topright', legend=sub('.*\\.vs\\.', '', group), cex=2, text.col='darkblue', bty='n')
                mtext(sub('\\.vs.*', '', group), side=3, line=1, at=(xlim[1]+abs(xlim[1])*0.05), cex=2, col='darkblue')
                mtext(sub('.*\\.vs\\.', '', group), side=3, line=1, at=(xlim[2]-abs(xlim[2])*0.05), cex=2, col='darkblue')
                
                }


            ## ###########################################################
            ##
            ##               selected points
            ## ###########################################################
            if(!is.null( volc[[paste('x', group, sep='.')]] ) & length(volc[[paste('x', group, sep='.')]]) ){

                for(i2 in 1:length( unlist( volc[[paste('x', group, sep='.')]]))){

                    volc.add.X[i2] <- as.numeric( unlist( volc[[paste('x', group, sep='.')]][i2] ) )
                    volc.add.Y[i2] <- as.numeric( unlist( volc[[paste('y', group, sep='.')]][i2] ) )
                    volc.add.text[i2] <- as.character( unlist( volc[[paste('text', group, sep='.')]][i2]) )
                    volc.add.col[i2] <- 'black'
                }
            }

            ## ###########################################################
            ##
            ##                 plot PPI stuff
            ##
            ## ###########################################################
            if(PPI){
                if(length(ppi.int.idx) > 0)
                    points(logFC[ppi.int.idx], logPVal[ppi.int.idx], col=col[ppi.int.idx], bg=col[ppi.int.idx], pch=pch.vec[ppi.int.idx], cex=cex.vec[ppi.int.idx])

                if(length(ppi.bait.idx) > 0)
                    points(logFC[ppi.bait.idx], logPVal[ppi.bait.idx], col='green', bg='green', pch=pch.vec[ppi.bait.idx], cex=cex.vec[ppi.bait.idx])

                legend('topleft', legend=leg, col=leg.col, pch=16, bty='n', cex=1.5, title=paste('Known interactors sig/det/tot'))

            }

            ## ########################################
            ##  draw ids of selected points
            if(length(volc.add.X) > 0)
                pointLabel(as.numeric(unlist(volc.add.X)), as.numeric(unlist(volc.add.Y)), labels=as.character(unlist(volc.add.text)), col=volc.add.col, offset=20, method='SANN', cex=input[[paste('cex.volcano.lab', group, sep='.')]])

        } ## end plotVolcano






        #######################################################################################
        ##
        ##                                profile plots
        ##
        #######################################################################################
        output$expr.profile <- renderPlot({
            if(is.null(global.results$data)) return()

            ## dataset
            if(is.null(global.results$table.log))
                tab <- data.frame(global.input$table)
            else
                tab <- data.frame(global.results$table.log)

            ## id column
            id.col.value <- global.param$id.col.value
            ## group vector
            grp <- global.param$grp
            ## group colors
            grp.col <- global.param$grp.colors
            grp.col.leg <- global.param$grp.colors.legend

            # update selection
            tab <- tab[, c(id.col.value, names(grp))]
            grp.col <- grp.col[names(grp)]
            grp.col.leg <- grp.col.leg[unique(grp)]
              
            withProgress({
                   setProgress(message = 'Processing...', detail= 'Generating profile plots')
                   makeProfileplot(tab, id.col.value, grp, grp.col, grp.col.leg, main='Before normalization')
            })

        })
        ###########################
        ## normalized
        output$expr.profile.norm <- renderPlot({
            if(is.null(global.results$data)) return()

            ## dataset
            tab <- data.frame(global.results$table.norm)

            ## id column
            id.col.value <- global.param$id.col.value
            ## group vector
            grp <- global.param$grp
            ## group colors
            grp.col <- global.param$grp.colors
            grp.col.leg <- global.param$grp.colors.legend

            
            # update selection
            tab <- tab[, c(id.col.value, names(grp))]
            grp.col <- grp.col[names(grp)]
            grp.col.leg <- grp.col.leg[unique(grp)]
            
            
            withProgress({
                   setProgress(message = 'Processing...', detail= 'Generating profile plots')
                   makeProfileplot(tab, id.col.value, grp, grp.col, grp.col.leg, main=paste(global.param$norm.data, 'normalized'))
            })

        })

        #######################################################################################
        ##
        ##                                    boxplots
        ##
        #######################################################################################

        ######################################
        ## without normalization
        output$expr.boxplot <- renderPlot({

            if(is.null(global.results$data)) return()

            ## dataset
            if(is.null(global.results$table.log))
                tab <- data.frame(global.input$table)
            else
                tab <- data.frame(global.results$table.log)

            ## id column
            id.col.value <- global.param$id.col.value
            ## group vector
            grp <- global.param$grp
            ## group colors
            grp.col <- global.param$grp.colors
            grp.col.leg <- global.param$grp.colors.legend
            
            
            # update selection
            tab <- tab[, c(id.col.value, names(grp))]
            grp.col <- grp.col[names(grp)]
            grp.col.leg <- grp.col.leg[unique(grp)]
            
               ##withProgress({
               ##    setProgress(message = 'Processing...', detail= 'Generating Boxplots')
            makeBoxplot(tab, id.col.value, grp, grp.col, grp.col.leg)
               ##})
        })

        #################################################
        ## with normalization
        output$expr.boxplot.norm <- renderPlot({

            if(is.null(global.results$data)) return()
            if(is.null(global.results$table.norm)) return()
            if(!is.null(error$msg)) return()

            ## dataset
            tab <- data.frame(global.results$table.norm)
            ## id column
            id.col.value <- global.param$id.col.value
            ## group vector
            grp <- global.param$grp
            ## group colors
            grp.col <- global.param$grp.colors
            grp.col.leg <- global.param$grp.colors.legend
            
            
            # update selection
            tab <- tab[, c(id.col.value, names(grp))]
            grp.col <- grp.col[names(grp)]
            grp.col.leg <- grp.col.leg[unique(grp)]
            

            makeBoxplot(tab, id.col.value, grp, grp.col, grp.col.leg, legend=T)
        })

        ######################################################################################
        ##
        ##                                   Correlation
        ##
        ######################################################################################
        output$multi.scatter <- renderPlot({
            ##cat('tesssttt')
            if(is.null(global.results$data)) return()
            global.plotparam$ms.max <- input$ms.max

            global.plotparam$ms.max.val <- isolate( input$ms.max.val )
            global.plotparam$ms.min.val <- isolate( input$ms.min.val )

            ##cat('maxval: ', global.plotparam$ms.max.val,'\n')
            plotMultiScatter( define.max=global.plotparam$ms.max, min.val=global.plotparam$ms.min.val, max.val=global.plotparam$ms.max.val )
            ##plotMultiScatter( define.max=input$ms.max, min.val=input$ms.min.val, max.val=input$ms.max.val )
        },
        width = function(){120*(ncol(data.frame(global.input$table))-1)},
       height= function(){120*(ncol(data.frame(global.input$table))-1)}
        )

        ###############################
        ## actual plot
        plotMultiScatter <- function(define.max, min.val, max.val){

            cat('\n-- plotMultiScatter  --\n')
            ## dataset
            ##tab <- data.frame(global.input$table)
            ## dataset
            if(is.null(global.results$table.log))
                tab <- data.frame(global.input$table)
            else
                tab <- data.frame(global.results$table.log)

            ##tab <- data.frame(global.input$table)
            ##View(tab)
            ## id column
            id.col <- global.param$id.col.value
            rownames(tab) <- tab[, id.col]
            ## table
            tab <- tab[, setdiff(colnames(tab), id.col)]
            ## get groups
            grp <-  global.param$grp
            grp <- sort(grp)
            tab <- tab[, names(grp)]

            ## mapping to colors
            ##grp.col <- rep('grey10', length(grp))
            ##grp.col[which(grp == input$label.g2)] <- 'darkblue'

            colnames(tab) <- chopString(colnames(tab), STRLENGTH)

           ###############################
            ## plot
            withProgress({
                setProgress(message = 'Processing...', detail= 'Calculating correlations')
                my.multiscatter(tab, repro.filt=global.results$values.filtered, grp=grp,  grp.col.legend=global.param$grp.colors.legend, define.max=define.max, max.val=max.val, min.val=min.val)
            })
        }
        ################################
        ## download image, Multiscatter
        ##output$downloadMS <- downloadHandler(
        ##    filename =  paste( 'multiscatter.pdf'),
        ##    content = function(file){
        ##        pdf(file, height=100*ncol(global.input$table)*(11/800), width=100*ncol(global.input$table)*(11/800))
        ##        withProgress({
        ##            ##setProgress(message = 'Processing...', detail= 'Calculation correlations')
        ##            plotMultiScatter(define.max=input$ms.max, max.val=input$ms.max.val, min.val=input$ms.min.val)
        ##        })
        ##        dev.off()
        ##    }
        ##)

        #####################################################
        ## correlation matrix
        output$correlation.matrix <- renderPlot({
            if(is.null(global.results$data)) return()
             withProgress({
                 setProgress(message = 'Processing...', detail= 'Generating Heatmap')
                 ##plotCorrMat(lower=input$cm.lower, upper=input$cm.upper, display_numbers=input$cm.numb, width=dynamicWidthHM( length(global.param$grp), unit='in'), height=dynamicWidthHM( length(global.param$grp), unit='in' ))
                 plotCorrMat(lower=input$cm.lower, upper=input$cm.upper, display_numbers=input$cm.numb, width=dynamicWidthHM( length(global.param$grp), unit='in'), height=dynamicWidthHM( length(global.param$grp), unit='in' ))
             })
        }, width=1200, height=1000)

        #####################################################
        ## correlation matrix transposed
       ## output$correlation.matrix.trans <- renderPlot({
       ##     if(is.null(global.results$data)) return()
       ##      withProgress({
       ##          setProgress(message = 'Processing...', detail= 'Generating Heatmap')
       ##          plotCorrMat(lower=input$cm.lower, upper=input$cm.upper, trans=T, display_numbers=input$cm.numb)
       ##      })
       ## }, width=1200, height=1000)




        ###################################################
        ## correlation matrix
        ###################################################
        plotCorrMat <- function(filename=NA, lower=c('pearson', 'spearman', 'kendall', 'pcor'), upper=c('pearson', 'spearman', 'kendall', 'pcor'), trans=F, display_numbers=T, width=12, height=12){


            cat('\n-- plotCorrMat --\n')

            ## dataset
            tab <- data.frame(global.input$table)
            ## id column
            id.col <- global.param$id.col.value
            ## class vector
            grp <- sort(global.param$grp)
            grp.col.legend <- global.param$grp.colors.legend

            ## table
            tab <- tab[, setdiff(colnames(tab), id.col)]
            tab <- tab[, names(grp)]

            ## transpose
            if(trans)
                tab=t(tab)

            ###########################
            ## calculate correlations
            ## withProgress({
            ##     setProgress(message = 'Processing...', detail= 'Calculation correlations')
            cm.upper <- cor(tab, use='pairwise', method=match.arg(upper))
            cm.lower <- cor(tab, use='pairwise', method=match.arg(lower))
            ##})
            ###########################
            ## initialize correlation matrix
            cm <- matrix(NA, ncol=ncol(cm.upper),nrow=nrow(cm.upper), dimnames=dimnames(cm.upper))
            cm[ lower.tri(cm, diag=T) ] <- cm.lower[lower.tri(cm.lower, diag=T)]
            cm[ upper.tri(cm, diag=F) ] <- cm.upper[upper.tri(cm.upper, diag=F)]

            ## colors
            color.breaks = seq(-1, 1, length.out=20)
            color.hm=colorRampPalette(c('blueviolet','blue','cyan', 'aliceblue', 'white' , 'blanchedalmond', 'orange', 'red', 'darkred'))(length(color.breaks))

            ## gaps between groups
            gaps.column=cumsum(table(grp))
            gaps.row=gaps.column

            ## annotation of rows/columns
            anno=data.frame(Group=grp)
            anno.color=list(Group=grp.col.legend)

            Rowv=F
            Colv=F

            ##setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=0.9, name="vp", just=c("right","top"))), action="prepend")
            setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=0.9, name="vp", just=c("right","top"))))
            if(trans)
                pheatmap(cm, fontsize_row=10, fontsize_col=10,
                     cluster_rows=Rowv, cluster_cols=Colv, border_col=NA, col=color.hm, filename=filename, labels_col=chopString(colnames(cm), STRLENGTH), labels_row=chopString(rownames(cm), STRLENGTH), main='', display_numbers=display_numbers, fontsize_number=100/ncol(cm)+10, breaks=color.breaks, width=width, height=height)
            else
                pheatmap(cm, fontsize_row=10, fontsize_col=10,
                     cluster_rows=Rowv, cluster_cols=Colv, border_col=NA, col=color.hm, filename=filename, labels_col=chopString(colnames(cm), STRLENGTH), labels_row=chopString(rownames(cm), STRLENGTH), main='', annotation_col=anno, annotation_colors=anno.color,  annotation_row=anno, display_numbers=display_numbers, fontsize_number=100/ncol(cm)+10, breaks=color.breaks, gaps_col=gaps.column, gaps_row=gaps.row, width=width, height=height)
            setHook("grid.newpage", NULL, "replace")

            ## add corr coeff
            grid.text(paste(match.arg(upper)), y=.995, x=.4, gp=gpar(fontsize=25))
            grid.text(paste(match.arg(lower)), x=-0.01, rot=90, gp=gpar(fontsize=25))
        }



        ## ################################################################################
        ## boxplots: correlations per group
        ## ################################################################################
        output$corr.box.group <- renderPlot({

            if(is.null(global.results$data)) return()
            if(!is.null(error$msg)) return()

            ## dataset
            tab <- data.frame(global.input$table)
            ## id column
            id.col <- global.param$id.col.value
            ## class vector
            grp <- sort(global.param$grp)
            grp.col.legend <- global.param$grp.colors.legend

            ## table
            tab <- tab[, setdiff(colnames(tab), id.col)]
            tab <- tab[, names(grp)]


            if( is.null(global.results$cm) ) ## will always be NULL atm
                cm <- cor(tab, use='pairwise', method='pearson')



        })

        ## #################################################################################
        ## 
        ##                               Morpheus
        ##
        ## #################################################################################
        output$HM.morpheus <- renderMorpheus({
          
          if(!is.null(error$msg)) return()
          
          ## extract results
          res = global.results$filtered
          
        
          ## require at least three significant hits
          validate(need(nrow(res) > 1, 'Need at least 2 features to draw a heatmap!'))
          
          
          ## extract expression values
          res <- res[, names(global.param$grp)]
          res <- data.matrix(res)
          
          ## morpheus
          withProgress({
            setProgress(message = 'Creating morpheus widget...', detail= 'hold on')
            morpheus(res, Rowv = F, Colv = F)
          })
          
        })

        ####################################################################################
        ##
        ##                                  Heatmap
        ##
        ####################################################################################
        output$HM <- renderPlot({

            ## if(is.null(global.results$data)) return()
            if(!is.null(error$msg)) return()

            ######################################
            ## extract results
            ##filter.res()
            ##tab.select <- input$mainPage
            ###filter.res()#--
            ##updateNavbarPage(session, 'mainPage', selected=tab.select)

            res = global.results$filtered



            ######################################
            ## require at least three significant hits
            ##if(nrow(res) < 3) return()
            validate(need(nrow(res) > 1, 'Need at least 2 features to draw a heatmap!'))

            #######################################
            ## heatmap title
            hm.title <- paste('filter:', global.param$filter.type, ' / cutoff:', global.param$filter.value, sep='')
            hm.title = paste(hm.title, '\nsig / total: ', nrow(res), ' / ', nrow( global.results$data$output ), sep='')

            #######################################
            ## extract expression values
            res = res[, names(global.param$grp)]

            ##@#####################################
            ##  dimensions depending on no. rows/columns
            cw <- cwHM(ncol(res))
         
            #if(!is.null(global.input$cdesc)){
            if(!is.null(global.param$anno.col)){
              #hm.cdesc <- data.frame( global.input$cdesc[names(global.param$grp),  global.param$cdesc.selection], stringsAsFactors = F )
              anno.col=global.param$anno.col
              anno.col.color=global.param$anno.col.color
            } else {
              anno.col=data.frame(Group=global.param$grp)
              anno.col.color=list(Group=global.param$grp.colors.legend)
            }
              #hm.cdesc <- global.param$grp
            #save(hm.cdesc, file='hm.cdesc.RData')
            ######################################
            ## plot
            if(input$hm.max){
                withProgress({
                     setProgress(message = 'Processing...', detail= 'Generating Heatmap')
                     #plotHM(res=res, grp=global.param$grp, grp.col=global.param$grp.colors, grp.col.legend=global.param$grp.colors.legend,  hm.clust=input$hm.clust, hm.title=hm.title, hm.scale=input$hm.scale, cellwidth=cw, fontsize_row=input$cexRow, fontsize_col=input$cexCol, max.val=input$hm.max.val, style=global.param$which.test, cdesc=hm.cdesc, cdesc.grp=global.param$grp.gct3)
                  plotHM(res=res, grp=global.param$grp, grp.col=global.param$grp.colors, grp.col.legend=global.param$grp.colors.legend,  hm.clust=input$hm.clust, hm.title=hm.title, hm.scale=input$hm.scale, cellwidth=cw, fontsize_row=input$cexRow, fontsize_col=input$cexCol, max.val=input$hm.max.val, style=global.param$which.test, anno.col=anno.col, anno.col.color=anno.col.color)
                  
                  })
            } else {
                 withProgress({
                    setProgress(message = 'Processing...', detail= 'Generating Heatmap')
                    #plotHM(res=res, grp=global.param$grp, grp.col=global.param$grp.colors, grp.col.legend=global.param$grp.colors.legend,  hm.clust=input$hm.clust, hm.title=hm.title, hm.scale=input$hm.scale, cellwidth=cw, fontsize_row=input$cexRow, fontsize_col=input$cexCol, style=global.param$which.test, cdesc=hm.cdesc, cdesc.grp=global.param$grp.gct3)
                   plotHM(res=res, grp=global.param$grp, grp.col=global.param$grp.colors, grp.col.legend=global.param$grp.colors.legend,  hm.clust=input$hm.clust, hm.title=hm.title, hm.scale=input$hm.scale, cellwidth=cw, fontsize_row=input$cexRow, fontsize_col=input$cexCol, style=global.param$which.test, anno.col=anno.col, anno.col.color=anno.col.color)
                   })
            }
        },
        width = function(){
            width=dynamicWidthHM(length(global.param$grp))
            return(width)
        },
        height= function(){
            ##filter.res()
            height=min( dynamicHeightHM( nrow(global.results$filtered)), 1200 )
            return(height)}
        )
        ####################################################################################
        ##
        ##                                interactive  Heatmap
        ##
        ####################################################################################
        output$HM.int <- renderPlotly({
          
          ## if(is.null(global.results$data)) return()
          if(!is.null(error$msg)) return()
          
          ######################################
          ## extract results
          res = global.results$filtered
  
          ######################################
          ## require at least three significant hits
          validate(need(nrow(res) > 1, 'Need at least 2 features to draw a heatmap!'))
          
          #######################################
          ## heatmap title
          hm.title <- paste('filter:', global.param$filter.type, ' / cutoff:', global.param$filter.value, sep='')
          hm.title <- paste(hm.title, '\nsig / total: ', nrow(res), ' / ', nrow( global.results$data$output ), sep='')
          
          #######################################
          ## extract expression values
          res = res[, names(global.param$grp)]
          
          ##@#####################################
          ##  dimensions depending on no. rows/columns
          cw <- cwHM(ncol(res))
          
          if(!is.null(global.param$anno.col)){
            #hm.cdesc <- data.frame( global.input$cdesc[names(global.param$grp),  global.param$cdesc.selection], stringsAsFactors = F )
            anno.col=global.param$anno.col
            anno.col.color=global.param$anno.col.color
          } else {
            anno.col=data.frame(Group=global.param$grp)
            anno.col.color=list(Group=global.param$grp.colors.legend)
          }
          
          #if(!is.null(global.input$cdesc))
          #  hm.cdesc <- data.frame(global.input$cdesc[names(global.param$grp), global.param$cdesc.selection], stringsAsFactors = F)
          #else
          #  hm.cdesc <- NULL
          #  hm.cdesc <- global.param$grp
          
          ######################################
          ## plot
          if(input$hm.int.max){
            withProgress({
              setProgress(message = 'Processing...', detail= 'Generating Heatmap')
              #plotHM(res=res, grp=global.param$grp, grp.col=global.param$grp.colors, grp.col.legend=global.param$grp.colors.legend,  hm.clust=input$hm.int.clust, hm.title=hm.title, hm.scale=input$hm.int.scale, cellwidth=cw, max.val=input$hm.int.max.val, style=global.param$which.test, cdesc=hm.cdesc, cdesc.grp=global.param$grp.gct3, plotly = T)
              plotHM(res=res, grp=global.param$grp, grp.col=global.param$grp.colors, grp.col.legend=global.param$grp.colors.legend,  hm.clust=input$hm.int.clust, hm.title=hm.title, hm.scale=input$hm.int.scale, cellwidth=cw, max.val=input$hm.int.max.val, style=global.param$which.test, anno.col=anno.col, anno.col.color=anno.col.color, plotly = T)
            })
          } else {
            withProgress({
              setProgress(message = 'Processing...', detail= 'Generating Heatmap')
              #plotHM(res=res, grp=global.param$grp, grp.col=global.param$grp.colors, grp.col.legend=global.param$grp.colors.legend,  hm.clust=input$hm.int.clust, hm.title=hm.title, hm.scale=input$hm.int.scale, cellwidth=cw, style=global.param$which.test, cdesc=hm.cdesc, cdesc.grp=global.param$grp.gct3, plotly = T)
              plotHM(res=res, grp=global.param$grp, grp.col=global.param$grp.colors, grp.col.legend=global.param$grp.colors.legend,  hm.clust=input$hm.int.clust, hm.title=hm.title, hm.scale=input$hm.int.scale, cellwidth=cw, style=global.param$which.test, anno.col=anno.col, anno.col.color=anno.col.color, plotly = T)
              
            })
          }
        }
        )
        ## ################################################################
        ##   
        ##                         fanplot 
        ##
        ## ################################################################
        output$HC.fan <- renderPlot({
          
          if(!is.null(error$msg)) return()
          
          ######################################
          ## extract results
          res = global.results$filtered
          
          ######################################
          ## require at least three significant hits
          validate(need(nrow(res) > 1, 'Need at least 2 features to draw a heatmap!'))
        
          #######################################
          ## extract expression values
          res = res[, names(global.param$grp)]
          
          if(!is.null(global.input$cdesc))
            hm.cdesc <- global.input$cdesc
          else
            hm.cdesc <- NULL
            #hm.cdesc <- global.param$grp
          
          ######################################
          ## plot
            withProgress({
              setProgress(message = 'Processing...', detail= 'Fanplot')
              plotFAN(res=res, grp=global.param$grp, grp.col=global.param$grp.colors, grp.col.legend=global.param$grp.colors.legend)
            })

          
        })
        
        
        
        
        ##@################################################################################
        ## histogram of p-values
        ##@################################################################################
        output$pval.hist <- renderPlot({

            if(is.null(global.results$data)) return()
            if(!is.null(error$msg)) return()

            groups.comp <- unique(global.param$grp.comp)

            res = global.results$data$output

            ############################################
            ## mod T
            if(global.param$which.test != 'mod F'){
                par(mfrow=c(length(groups.comp),1))
                for(g in groups.comp){

                    pval <- res[, paste('P.Value', g, sep='.')]
                    hist(pval, breaks=50, main=paste('Number of tests: ', sum(!is.na(pval)), '',sep=''), xlab='P-value', cex.main=2.5, cex.axis=1.8, cex.lab=1.8, col='darkblue', border=NA)
                    legend('top', legend=g, cex=2)
                }
            ############################################
            ## mod F
            } else {
                pval <- res[, paste('P.Value')]
                hist(pval, breaks=50, main=paste('Number of tests: ', sum(!is.na(pval)), '',sep=''), xlab='P-value', cex.main=2.5, cex.axis=1.8, cex.lab=1.8, col='darkblue', border=NA)
                ##legend('top', legend=paste(groups), cex=2)
            }

        },
        width = function(){ width=1000},
        height= function(){ height=500*ifelse( global.param$which.test != 'mod F', length(unique(global.param$grp.comp)), 1 )} )

        ######################################################################################
        ##
        ##                                 PCA
        ##
        ######################################################################################
        output$run.pca <- renderText({
            if(is.null(global.results$data) | is.na(global.param$filter.value)) return()

            res <- global.results$filtered

            validate(need(nrow(res) > 2, 'Need at least 3 features to perform PC.'))

            grp <- global.param$grp

            ## run PCA
            withProgress(message = 'PCA...',{
                    pca=my.prcomp2( res, grp )
                  })

            ## store results
            global.results$pca <- pca

            ## short summary, same as 'PCA' is generating
            txt = paste('<p><font size=\"5\">PCA model of a mean-centered and scaled matrix of ',length(grp), ' by ', nrow(pca$loadings), '.</font></p>')
            txt = paste(txt, '<p><font size=\"5\">Number of PCs to cover 90% of the variance:', min(which((cumsum(pca$var)/pca$totalvar) > .9)), '.</font></p>')

            HTML(txt)


        })


        ################################################
        ## PCA variance plots
        output$pca.var <- renderPlot({

          if(is.null(global.results$pca)) return()

          pca <- global.results$pca

          plotPCAvar(pca)

        })

        ## ###########################################################
        ## PCA plotly 2D scatter
        output$pcaxy.plotly <- renderPlotly({
          if(is.null(global.results$data) | is.na(global.param$filter.value) | is.null(global.results$pca)) return()
          ##if(!is.null(error$msg)) return()

          pca <- global.results$pca
          grp <- global.param$grp
          grp.unique <- unique(grp)
          grp.colors <- global.param$grp.colors[names(grp)]

          # selected PCs
          pca.x <- as.numeric(sub('PC ','', input$pca.x))
          pca.y <- as.numeric(sub('PC ','', input$pca.y))

          # build a data frame for plotly
          pca.mat = data.frame(
            PC1=pca$scores[, pca.x],
            PC2=pca$scores[, pca.y]
          )
          rownames(pca.mat) <- rownames(pca$scores)

          p <- plot_ly( pca.mat, type='scatter', mode='markers' )
          
          for(g in grp.unique){
            grp.tmp <- names(grp)[grp == g]
            p <-  add_trace(p, x=pca.mat[grp.tmp , 'PC1'], y=pca.mat[grp.tmp, 'PC2'], type='scatter', mode='markers', marker=list(size=15, color=grp.colors[grp.tmp]), text=grp.tmp, name=g  )
            #p <- plot_ly( x=pca.mat$PC1, y=pca.mat$PC2, type='scatter', mode='markers', marker=list(size=20, color=grp.colors), text=rownames(pca.mat) )
          }
          
          p <- layout(p, title=paste('PC', pca.x,' vs. PC', pca.y, sep=''), xaxis=list(title=paste('PC', pca.x)), yaxis=list(title=paste('PC', pca.y)) )

        })

        ################################################
        ## PCA plotly 3D scatterplot
        output$pcaxyz.plotly <- renderPlotly({
            if(is.null(global.results$data) | is.na(global.param$filter.value) | is.null(global.results$pca)) return()
            ##if(!is.null(error$msg)) return()

            validate(need(length(global.param$grp) >= 3, 'Cannot generate 3D plots from data with dimension < 3.'))

            ##if(length(global.param$grp) < 3) return()

            pca <- global.results$pca
            grp <- global.param$grp
            grp.unique <- unique(grp)
            grp.colors <- global.param$grp.colors[names(grp)]
            
            # selected PCs
            pca.x <- as.numeric(sub('PC ','', input$pca.x))
            pca.y <- as.numeric(sub('PC ','', input$pca.y))
            pca.z <- as.numeric(sub('PC ','', input$pca.z))

            pca.mat = data.frame(
                PC1=pca$scores[, pca.x],
                PC2=pca$scores[, pca.y],
                PC3=pca$scores[, pca.z]
            )
            rownames(pca.mat) <- rownames(pca$scores)

            ###########################
            ## plot
            p <- plot_ly( pca.mat, type='scatter3d', mode='markers' )
            for(g in grp.unique){
              grp.tmp <- names(grp)[grp == g]
              p <- add_trace(p, x=pca.mat[grp.tmp, 'PC1'], y=pca.mat[grp.tmp,'PC2'], z=pca.mat[grp.tmp,'PC3'], type='scatter3d', mode='markers', marker=list(size=15, color=grp.colors[grp.tmp]), text=grp.tmp, name=g  )
            }
            #p <- plot_ly(  x=pca.mat$PC1, y=pca.mat$PC2, z=pca.mat$PC3, type='scatter3d', mode='markers', marker=list(size=15, color=grp.colors), text=rownames(pca.mat) )
            p <- layout(p, title=paste('PC', pca.x,' vs. PC', pca.y, 'vs. PC', pca.z, sep=''), scene=list( xaxis=list(title=paste('PC', pca.x)), yaxis=list(title=paste('PC', pca.y)), zaxis=list(title=paste('PC', pca.z))) )

        })
        ####################################################
        ## PCA loadings
        output$pca.loadings <- renderPlot({

          if(is.null(global.results$data) | is.na(global.param$filter.value) | is.null(global.results$pca)) return()
          if(!is.null(error$msg)) return()
          if(length(global.param$grp) < 3) return()

          pca <- global.results$pca

          pca.x <- as.numeric(sub('PC ','', input$pca.x))
          pca.y <- as.numeric(sub('PC ','', input$pca.y))
          pca.z <- as.numeric(sub('PC ','', input$pca.z))


          topn=input$pca.load.topn


          plotPCAloadings( pca, topn, pca.x, pca.y, pca.z )

       })

        ####################################################
        ##  PCA loadings as a scatter
        ####################################################
        output$scatter.pca.loadings <- renderPlot({

                if(is.null(global.results$data) | is.na(global.param$filter.value) | is.null(global.results$pca)) return()
                if(!is.null(error$msg)) return()
                if(length(global.param$grp) < 3) return()

                pca <- global.results$pca
                pca.x <- as.numeric(sub('PC ','', input$pca.x))
                pca.y <- as.numeric(sub('PC ','', input$pca.y))
                pca.z <- as.numeric(sub('PC ','', input$pca.z))
                topn=input$pca.load.topn

                scatterPlotPCAloadings( pca, topn, pca.x, pca.y, pca.z )

        })

        ###############################################
        ## static PCA plot
        ###############################################
        plotPCA <- function(pca.x, pca.y, pca.z, plot=T){

            ##filter.res()
            ##tab.select <- input$mainPage
            ##filter.res()#--
            ##updateNavbarPage(session, 'mainPage', selected=tab.select)

            res <- global.results$filtered

            ##View(res)

            ## require at least three rows

            validate(need(nrow(res) > 2, 'Need at least 3 features to perform PC.'))

            ##if(nrow(res) < 3) return()

            ## get groups
            grp <- global.param$grp

            ## require at least 2 columns (2D PCA)
            if(length(names(grp)) < 2) return()

            ## mapping to colors
            grp.col <- global.param$grp.colors
            grp.col.leg <- global.param$grp.colors.legend

            # update selection
            grp.col <- grp.col[names(grp)]
            grp.col.leg <- grp.col.leg[unique(grp)]
            #res <- res[, names(grp)]
            
            ## remove missing values
            rm.idx <- apply(res, 1, function(x) sum(is.na(x)) + sum(is.infinite(x)))
            rm.idx <- which(rm.idx > 0)
            if(length(rm.idx)>0) res <- res[-rm.idx, ]
            if(nrow(res) < 3) return()

            ##View(res[, names(grp)])
            ##View(t(res[, names(grp)]))

            ## plot
            pca <- my.prcomp(t(res[, names(grp)]), col=grp.col, plot=plot, rgl=F, main='', cex.points=5, leg.vec=names(grp.col.leg), leg.col=grp.col.leg)

            return(pca)
        }


})


###########################################################
