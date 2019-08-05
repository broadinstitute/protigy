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

        #####################################
        ## reactive variables to store
        ## data accross a session
        #####################################

        ## error messages
        error <- reactiveValues()

        ## test results
        global.results <-  reactiveValues(

            ## data tables
            data=NULL,
            table.norm=NULL,
            table.log=NULL,
            table.repro.filt=NULL,
            filtered=NULL,

            #export.results=F,
            export.rmd=F,
            export.xls=F,
            export.gct=F,
            export.results=F, 
            
            pca=NULL,
            cm=NULL,
            
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
            session.saved=F,             ## flag whether session has been saved (v0.8.0.2)
            
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
            
            na.filt.val=100,            ## max. % missin  values

            
            filter.type='adj.p',         ## default filter
            filter.value=0.05,           ## default filter value

            run.test=0,                   ## number of times the 'Run analysis' button has been pressed during a user session

            update.ppi.select= FALSE,     ## trigger selectize, volcano
            update.ppi.select.scat=FALSE, ## trigger selectize, scatterplot
            collapse.ppi=TRUE,             ## should PPI query panel be collapsed?
            
            update.pca=FALSE,
            update.cm=FALSE

        )

        #####################################################
        ## plotting parameter
        global.plotparam <- reactiveValues(

            ## multiscatter
            ms.max=FALSE,
            ms.min.val=-4,
            ms.max.val=4,
            ms.robustify=T,

            ## volcano
            volc.ps=1.8,         ## point size
            volc.ls=1,         ## label size
            volc.grid=T,       ## grid
            volc.maxp=100,     ## max. -log10 p-value
            volc.hyper.fc=1,    ## min. FC for hyperbolic curve
            volc.hyper.curv=3,  ## curvation parameter for hyperbol. curve
            
            #volc.reset=F,    ## flag to trigger reset of annotated points
            volc.init=T,
            
            ## heatmap
            hm.cexCol=8,
            hm.cexRow=3,
            hm.scale="none",
            hm.max=FALSE,
            hm.max.val=4,
            hm.show.rownames=T,
            hm.show.colnames=T,
            
            ## PCA
            pca.x='PC 1',
            pca.y='PC 2',
            pca.z='PC 3',
            pca.load.topn=20,

            ## correlation matrix
            cm.upper='pearson',
            cm.lower='spearman',
            cm.numb=FALSE,
            
            ## fanplot
            HC.fan.show.tip.label=T,
            HC.fan.tip.cex=1
        )

        ## coordinates in volcano plot
        volc <- reactiveValues()
        volc.brush <- reactiveValues()

        ################################################################################
        ##
        ##                                instructions / help pages
        ##
        ################################################################################
        callModule(printHTML, id='getting.started', what='gs', global.input = global.input, error=error)
        callModule(printHTML, id='change.log', what='cl', global.input = global.input, error=error)
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
          
            if(is.null(error$title))
              etitle <- 'Error'
            else
              etitle <- error$title
            
            shinyalert(etitle, error$msg, type = "error")
            #HTML(paste('<p align=\"center\"><font size=\"5\" color=\"red\">', error$msg,'</font></p>'))
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
            updateCheckboxInput(session, 'export.gct.file', 'GCT files: 1) original data and 2) singed log P-values', value=!input$export.toggle.all)
            updateCheckboxInput(session, 'export.cm', 'Correlation matrix', value=!input$export.toggle.all)
            updateCheckboxInput(session, 'export.profile', 'Profile plot', value=!input$export.toggle.all)

        })


        ## ########################################################
        ## session name
        output$session.label <- renderMenu({
          
          #if(is.null(session$user)) return()
          if(!global.param$session.saved) return()
          label <- paste("Session:", global.param$label)
          
          notificationItem(label, status='success', shiny::icon("folder", "fa-1x", lib='glyphicon')) 
           
          
        })
        
        
        ## ########################################################
        ## logged user
        output$logged.user <-renderMenu({

            if(is.null(session$user)) return()
            user <- session$user
          
            notificationItem(user, shiny::icon("user", "fa-1x", lib='glyphicon'), status='success') ##, ic='users', status='info', text='test'),
        })

        ## #########################################################
        ## logout user
        output$logout <- renderMenu({
            if(is.null(session$user)) return()
            notificationItem('Logout', icon=shiny::icon("sign-out", "fa-1x"), status='success', href="__logout__")
        })


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
                                   
                                     fluidRow(
                                       column(width=6,
                                              box(title="Session name",
                                                  fluidPage(
                                                    fluidRow(
                                                      HTML('Specify a name for the current session.')
                                                      ),
                                                    fluidRow(
                                                      textInput( 'label', '', value=global.param$label, width=300)
                                                    ),
                                                    fluidRow(
                                                      HTML('Save the current state of your session (SSPro).')
                                                    ),
                                                    fluidRow(
                                                      actionButton('export.save.session', 'Save session')
                                                    )
                                                  ),
                                                  status = "primary",
                                                  solidHeader = T,
                                                  width=NULL),
                                              
                                              ## ################################################
                                              ## action and download button for rmd, xls and zip
                                              box(title="Export results", 
                                                  fluidPage(
                                                    fluidRow(column(3, HTML('<a href="https://rmarkdown.rstudio.com/" target="_blank_">Markdown</a> report')), 
                                                             column(3, HTML('Spreadsheet')),
                                                             column(3, HTML('GCT')),
                                                             column(3, HTML('Gimme all!'))),
                                                    fluidRow(column(3, 
                                                           if(!global.results$export.rmd)
                                                             actionButton('export.rmd', 'html', icon = icon("code", lib="font-awesome"))
                                                            #actionButton('export.rmd', 'html', icon = icon("code", lib="font-awesome"), style="color: #000000; background-color: #ffa09b; border-color: #c7c9c7")
                                                           else
                                                             downloadButton('download.rmd', 'html', style="color: #00a805") 
                                                            #downloadButton('download.rmd', 'html', style="color: #000000; background-color: #b7ff9b; border-color: #c7c9c7") 
                                                           ),
                                                    column(3,
                                                           if(!global.results$export.xls)
                                                            actionButton('export.xls', 'xlsx', icon = icon("table", lib="font-awesome"))
                                                           else
                                                            downloadButton('download.xls', 'xlsx', style="color: #00a805")
                                                           ),
                                                    column(3,
                                                           if(!global.results$export.gct)
                                                             actionButton('export.gct', 'GCT', icon = icon("table", lib="font-awesome"))
                                                           else
                                                             downloadButton('download.gct', 'GCT', style="color: #00a805")
                                                    ),
                                                    
                                                    column(3, if(!global.results$export.results)
                                                      actionButton('export.results', 'zip', icon = icon("archive", lib="font-awesome"))
                                                      else
                                                        downloadButton('download.results', 'zip', style="color: #00a805"))
                                                  )),
                                                  status = "primary",
                                                  solidHeader = T,
                                                  width=NULL)#,
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
                                                  checkboxInput('export.gct.file', 'GCT files: 1) original data and 2) singed log P-values',value=T),
                                                  checkboxInput('export.cm', 'Correlation matrix',value=T),
                                                  checkboxInput('export.profile', 'Profile plot',value=T),
                                                  
                                                  status = "primary",
                                                  solidHeader = T,
                                                  width=NULL
                                                  
                                              )
                                       )
                                     ) # end fluidRow

            ) # end tabPanel
            
         
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
                                        box(title="Missing values", solidHeader = T, status = "primary",width = 12,
                                            plotlyOutput('summary.missing.data.row'))
                                    ),
                                    fluidRow(
                                      column(6, 
                                             box(title="Significant features (down-regulated)", solidHeader = T, status = "primary",width = 12,
                                         
                                             plotOutput('summary.upset.dn'))
                                            
                                          ),
                                      column(6,
                                         box( title="Significant features (up-regulated)", solidHeader = T, status = "primary",width = 12,  
                                              plotOutput('summary.upset.up')) 
                                         )   
                                      
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
                                                                                maxOptions=10,
                                                                                placeholder='Protein/Gene'#,
                                                                                #onInitialize = I('function() { this.setValue(""); }')
                                                                            )

                                                                            )
                                                             ),
                                                      column(3, checkboxGroupInput(paste('ppi.db.scat', grp.unique[i], sep='.' ), 'Source data', choices=c('BioGRID (human)' = 'bg', 'InWeb' = 'iw', 'Reactome (human)' = 'react'), selected=c('bg') ))#,
                                                      #column(3, checkboxInput( paste('ppi.show.labels.scat', grp.unique[i], sep='.' ), 'Show labels', value=F) )

                                                      ##column(1)
                                                  )),


                                              ## ###############################
                                              ## plot
                                              box( title=grp.unique[i], status = 'primary', solidHeader = T, width=12,
                                                  fluidRow(
                                                       column(12, align='center', plotlyOutput(  paste0('scatterplot.', grp.unique[i]) , width=800, height=600) )
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
                                                                  selectizeInput( inputId=gsub('\\.','', gsub('\\.', '', paste0('ppi.bait.', groups.comp[i])) ), 
                                                                                  label=NULL,
                                                                                  choices=NULL,
                                                                                  selected=NULL,
                                                                                  multiple=F,
                                                                                  options = list(
                                                                                    maxOptions=10,
                                                                                    placeholder='Protein/Gene'#,
                                                                                    #oninit = I('function() { this.setValue(" "); }')
                                                                                    #onInitialize = I('function() { this.setValue("Testtststs"); }')
                                                                                )
    
                                                                                )
                                                                 ),
                                                          column(3, checkboxGroupInput(paste('ppi.db', groups.comp[i], sep='.' ), 'Source data', choices=c('BioGRID (human)' = 'bg', 'InWeb' = 'iw', 'Reactome (human)' = 'react'), selected=c('bg') )),
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
                                                             box(width=NULL,  title='Selection', status = 'primary', solidHeader = T,
                                                                 fluidRow(
                                                                 column(6, actionButton(inputId=paste('volc.tab.reset', groups.comp[i], sep='.'), label='Remove all')),
                                                                 column(6, actionButton(inputId=paste('volc.tab.reset.select', groups.comp[i], sep='.'), label='Remove selected'))
                                                                 ),
                                                                 tags$hr(),
                                                                 dataTableOutput(paste('volc.tab.selected', groups.comp[i], sep='.'))
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
                                          
                                          column(2, checkboxInput('hm.show.colnames', 'Show column labels', value=global.plotparam$hm.show.colnames)),
                                          column(2, checkboxInput('hm.show.rownames', 'Show row labels', value=global.plotparam$hm.show.rownames)),
                                          column(8)
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
                                          column(1, checkboxInput('HC.fan.show.tip.label', label = 'Show names', value = global.plotparam$HC.fan.show.tip.label)),
                                          column(2, numericInput('HC.fan.tip.cex', label = 'Label size', value = global.plotparam$HC.fan.tip.cex, min = 0.1, max=10, step = 0.2, width = '40%')),
                                          column(9)
                                        ),
                                        fluidRow(
                                          column(12, align='center', plotOutput("HC.fan", height = 600, width=1200 ))
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
            # morph.tab <- tabPanel('Morpheus',
            #                       box(title='Morpheus', status = 'primary', solidHeader = T, width="100%", height="100%",
            #                        fluidRow(
            #                          
            #                          column(12, morpheusOutput(outputId = "HM.morpheus"#,
            #                                                    #height=min( dynamicHeightHM( nrow(global.results$filtered)), 1200 ), 
            #                                                    #width=dynamicWidthHM(length(global.param$grp)) 
            #                                                     )
            #                                 )
            #                        )
            #                        )
            # )
            # 
            # 

            #############################################
            ##
            ##               PCA
            ##
            #############################################
            pca.tab2 <- list()
            
            ## ##################################################
            ## Run PCA
            ## ##################################################
            pca.tab2[[1]] <- tabPanel('Explained variance',
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
                                                 column(12, align='center', plotlyOutput("pcaxy.plotly", width=800, height=600))
                                            )
                                          ),

                                          fluidRow(
                                            box( title='3D', status = 'primary', solidHeader = T, width = 800, height = 900,
                                                 column(12, align='center', plotlyOutput("pcaxyz.plotly", width=800, height=800))
                                            )
                                          )

                                      )
                                      )
            ## ####################################################################
            ## Loadings plot
            ## ####################################################################
            #pca.tab2[[3]] <- tabPanel('Loadings',
            #                          fluidPage(
            #                            fluidRow(
            #                                column(width=12,
            #                                box(title='Loadings by Ozan Aygun', solidHeader=T, status='primary', width=1000,## height=min( nrow(global.results$filtered), global.plotparam$pca.load.topn )*20+50,
            #                                    sliderInput("pca.load.topn", "Choose number of loadings", 1, 100, 20),
            #          ##                          plotOutput("pca.loadings")##, width=1000, height=min(nrow(global.results$filtered), global.plotparam$pca.load.topn )*20 )
            #          ##                        ),
            #          ##                      box(title = "PCA loadings scatterplots",solidHeader = T, status = 'primary',width = 1000,
            #                                    background = "navy",
            #                                    plotOutput("scatter.pca.loadings")
            #                                ))
            #
            #                            )
            #                          )
            #)


            ## ###########################################
            ## TABLE: filtered result table
            ##
            ## ###########################################
            table.tab <- tabPanel('Table',
                     fluidPage(
                         #fluidRow(column(12, tags$h3(paste('Result table')))),
                         #fluidRow(column(12, tags$br())),
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
                                                           box(title="Parameters", solidHeader=T, status="primary", width='100%',
                                                               column(4, checkboxInput('ms.robustify', 'Robust scale (99.99 % of data)', value=global.plotparam$ms.robustify)),
                                                               column(1, checkboxInput('ms.max', 'Define limits', value=global.plotparam$ms.max)),
                                                               column(2, numericInput( "ms.min.val", "min.", value=global.plotparam$ms.min.val, step=1)),
                                                               column(2, numericInput( "ms.max.val", "max.", value=global.plotparam$ms.max.val, step=1)),
                                                               column(3, actionButton('ms.update',label = 'Update plot'))
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
                                                                     column(12, plotOutput("correlation.matrix", 800, 800))
                                                                 )
                                                             )
                                                         )
                                                         )
 
            ##############################
            ## correlation boxplots
            qc.tabs[[6]] <- tabPanel('Correlation boxplots',
                                     
                                     fluidPage(
                                       fluidRow(
                                         box(title='Correlation boxplots', solidHeader=T, status="primary", width=1000, height=800,
                                             column(12, plotOutput("corr.box.group", 800, 600))
                                        )
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
                        
                        #######################################
                        ##              insert volcanos
                        #######################################
                        #if( !(global.param$which.test %in% c('mod F', 'none'))){
                        do.call(navbarMenu, volc.tabs),
                      #  } else {
                      #       list()
                      #     },
                       
                       #morph.tab,
                      
                       ## ####################################
                       ##          insert scatterplots
                       ## ####################################
                       do.call(navbarMenu, scat.tabs),

                       #######################################
                       ##              insert PCA
                       #######################################
                       #navbarMenu('PCA', pca.tab2[[1]], pca.tab2[[2]], pca.tab2[[3]]),
                       #navbarMenu('PCA', pca.tab2[[1]], pca.tab2[[2]], pca.tab2[[3]]),
                       navbarMenu('PCA', pca.tab2[[1]], pca.tab2[[2]]),
                       #######################################
                       ##           insert table preview
                       #######################################
                       table.tab,

                       #######################################
                       ## QC
                       if(global.param$which.test != 'none') {
                           navbarMenu('QC', qc.tabs[[1]], qc.tabs[[2]], qc.tabs[[3]], qc.tabs[[4]], qc.tabs[[5]], qc.tabs[[6]])
                       } else {
                           navbarMenu('QC', qc.tabs[[1]], qc.tabs[[2]], qc.tabs[[4]], qc.tabs[[5]], qc.tabs[[6]])
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
                      navbarMenu('PCA', pca.tab2[[1]], pca.tab2[[2]]),
                      
                      #######################################
                      ##           insert table preview
                      #######################################
                      table.tab,
                      
                      #######################################
                      ## QC
                      if(global.param$which.test != 'none') {
                        navbarMenu('QC', qc.tabs[[1]], qc.tabs[[2]], qc.tabs[[3]], qc.tabs[[4]], qc.tabs[[5]], qc.tabs[[6]])
                      } else {
                        navbarMenu('QC', qc.tabs[[1]], qc.tabs[[2]], qc.tabs[[4]], qc.tabs[[5]], qc.tabs[[6]])
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
            cat( '---------------------------------\n')
          
            ## ######################################
            ## generate session ID and prepare data
            ## directory
            #########################################
  
            ## generate 'session id'
            if(is.null(global.param$session)){
              
              cat('No session id found:', global.param$session, '\n')
              
              cat('creating session id...\n')
              
              ## set a seed for 'sample'
              ## it happened that the same session ids were created...
              session.id.ok <- FALSE
              while(!session.id.ok){
                
                set.seed(as.numeric(Sys.time()))
                global.param$session <- paste(paste(letters[sample(26, 5)], collapse=''), paste(sample(100,5), collapse=''), sep='')
              
                ## check whether the session ids exists as directory on the server
                if(!dir.exists( paste(DATADIR, global.param$user,'/' ,global.param$session, '/', sep='') )) session.id.ok <- TRUE
              }
              cat('session:', global.param$session, '\n')
            }
            
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
              HTML('<font size=\"3\"><b>Upload file (txt, csv, gct, gctx):</b></font>'),
              fileInput("file", "", accept=c('text/csv',
                       'text/comma-separated-values,text/plain',
                       '.csv', '.txt', '.tsv', '.gct', '.gctx')),
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
                #actionButton('session.manage', 'Manage sessions')
                
                )
        })

        ## ##################################################
        ##
        ##       observer for manage sessions module
        ##
        ## ##################################################
        #observeEvent(input$session.manage, {
        #  callModule(manageSessions, id = 'manageSessions', data.dir=DATADIR,  session = session )
        #})
        
        
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
          
          #save(tab, file='tab.RData')
          
          tab <- data.frame(Level=names(tab), Freq=as.character(unlist(tab)), stringsAsFactors = F )
          
          list(
            tab
          )
          
        })
        
        # #####################################
        # define groups for GCTv3
        observeEvent(input$update.grp.gct3, {
          
          tab <- global.input$table
          cdesc <- data.frame(global.input$cdesc)
          #if(nrow)
          #View(head(cdesc))
          # store grp coloum
          global.param$grp.gct3 <- input$grp.gct3
          
          # initialize grp file
          Column.Name <- colnames(tab)
          Experiment <- rep('', length(Column.Name))
          names(Experiment) <- Column.Name
          Experiment[ rownames(cdesc) ] <- make.names(cdesc[, input$grp.gct3])
          

          global.param$cdesc.all <- global.param$cdesc.selection <- setdiff(colnames(cdesc),  input$grp.gct3)
        
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
          if(ncol(cdesc) > 1){
            col.tmp <- cdesc.colors(cdesc, global.param$grp.gct3, grp.col.legend)
            global.param$anno.col.all <- global.param$anno.col <- col.tmp$anno.col
            global.param$anno.col.color.all <- global.param$anno.col.color <- col.tmp$anno.col.color
          }

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
            tab.colnames.names <- names(tab.colnames)
            names(tab.colnames) <- NULL
           # cat(tab.colnames[1:4], '\n')
            ## radio button to pick id column
            list(
                ## experimental design
                HTML('<font size=\"3\"><b>Export experimental design file:</b></font>\n<br>'),
                HTML('<br>'),
                downloadButton("exportTemplate", 'Export',  style="color: #000000", width='100'),
                HTML('<br><hr size=\"5\">'),
                actionButton("id.col", 'Next', width='100'),
                HTML('<hr size=\"5\">'),
                radioButtons( inputId = "id.col.value", label = "Choose ID column", choiceValues = tab.colnames.names, choiceNames = tab.colnames )
                
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
                actionButton( 'update.grp', 'Next', width='100')
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
          
              showModal(modalDialog(
          
                size='m',
                title = "Modify selection",
                footer = fluidRow(
                  column(6),
                  column(3, actionButton(inputId = 'update.groups.button.modal' , label='Update')),
                  column(3, modalButton(label='Close'))
                  ),
                fluidPage(
                  
                  fluidRow(
                    column(6, HTML('<b>Select groups</b>')),
                    column(6, HTML('<b>Select annotation tracks</b>'))
                    
                  ),
                  fluidRow(
                    column(6, actionButton(inputId = 'toggle.select.groups', label = 'Toggle selection')),
                    column(6, actionButton(inputId = 'toggle.select.anno', label = 'Toggle selection'))
                  ),
                  
                  fluidRow(
                    #if(length(unique(global.param$grp.comp.all)) > 10){
                    column( 6, checkboxGroupInput('select.groups', label=' ', 
                                                  choices = unique(global.param$grp.comp.all), 
                                                  selected = unique(global.param$grp.comp.selection) 
                                                  #selected = NULL
                            )),
                    #column( 6, checkboxGroupInput('select.groups', label=' ', choices = unique(global.param$grp.all), selected = unique(global.param$grp.selection))),
                    #} else{ column( 6, checkboxGroupInput('select.groups', label=' ', choices = unique(global.param$grp.comp.all), selected = unique(global.param$grp.comp.selection)))},
                    column( 6, checkboxGroupInput('select.anno', label=' ', 
                                                  choices = unique(global.param$cdesc.all), 
                                                  selected = unique(global.param$cdesc.selection)
                                                  #selected=NULL
                                                  ))
                  )
                ),
                easyClose = FALSE
              ))
          jqui_draggable(ui = '.modal-content')
        })
        ####################################################
        ## toggle group selection
        observeEvent(input$toggle.select.groups, {
          #if(input$select.groups.button == 0) return()
          updateCheckboxGroupInput(inputId = 'select.groups', session = session,
                                   selected = setdiff(unique(global.param$grp.comp.all), input$select.groups)
                                   #selected = setdiff( global.param$grp.comp.all, unique(global.param$grp.comp.selection))
                                   )
        })
        ####################################################
        ## toggle selection of annotation tracks
        observeEvent(input$toggle.select.anno, {
          #if(input$select.groups.button == 0) return()
          updateCheckboxGroupInput(inputId = 'select.anno', session = session,
                                   selected = setdiff(unique(global.param$cdesc.selection), input$select.anno)
                                   #selected = setdiff( global.param$grp.comp.all, unique(global.param$grp.comp.selection))
          )
        })
        
        
        # ##################################################
        # reset group selection if the type of test changed
        observeEvent( input$which.test, {
          global.param$grp.comp.selection <- global.param$grp.comp.all
          global.param$grp.selection <- global.param$grp.all
        })
        
        # ##################################################
        # update based on selection
        #observeEvent(c(input$update.groups.button.modal, input$run.test),{  
        observeEvent(input$update.groups.button.modal,{  
          
          ## ####################
          ## update class vector
          
          # super important! to memorize previous selections in the modal window
          global.param$grp.comp.selection <- input$select.groups
          #global.param$grp.selection <- input$select.groups
          grp.selection <- global.param$grp.all
          
          ## extract groups selected in the modal window
          grp.unique <- unique( unlist( strsplit( sub('\\.vs\\.', ' ', input$select.groups), ' ')))
          
          ####################################
          ## make sure at least one group 
          ## has been selected
          if(length(grp.unique) == 0){
            shinyalert("No data selected!", "Please slect at least one group.", type = "error")
          }
          
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
          
        })
        
       
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

                     #radioButtons('filt.data', 'Filter data', choices=c('Reproducibility', 'StdDev', 'MissingValues', 'none'), selected=global.param$filt.data ),
                     radioButtons('filt.data', 'Filter data', choices=c('Reproducibility', 'StdDev', 'none'), selected=global.param$filt.data ),
                     #sliderInput('na.filt.val', 'Max. % missing values', min=0, max=100, value=global.param$na.filt.val),
                     radioButtons('which.test', 'Select test', choices=c('One-sample mod T', 'Two-sample mod T', 'mod F', 'none'), selected=global.param$which.test),

                     actionButton('run.test', 'Run analysis!'),
                     br(),
                     hr(),
                     actionButton('select.groups.button', 'Selected groups')
      
                )
            }
            ## ###################################################
            ## no filter
            else if(input$filt.data == 'none'){

                list(
                  
                     radioButtons('log.transform', 'Log-transformation', choices=c('none', 'log10', 'log2'), selected=global.param$log.transform),
                     radioButtons('norm.data', 'Data normalization', choices=c('Median', 'Median-MAD', '2-component', 'Quantile', 'none'), selected=global.param$norm.data),

                     radioButtons('filt.data', 'Filter data', choices=c('Reproducibility', 'StdDev', 'none'), selected='none'),
                     
                   #  sliderInput('na.filt.val', 'Max. % missing values', min=0, max=100, value=global.param$na.filt.val),
                     
                     radioButtons('which.test', 'Select test', choices=c('One-sample mod T', 'Two-sample mod T', 'mod F', 'none'), selected=global.param$which.test),

                     actionButton('run.test', 'Run analysis!'),
                     br(),
                     hr(),
                    # fluidRow(
                    #   column(6, actionButton('select.groups.button', 'Select Groups')),
                    #   column(6, actionButton('select.anno.button', 'Select Annotation Tracks'))
                     #)
                     actionButton('select.groups.button', 'Select groups')
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
                    # sliderInput('na.filt.val', 'Max. % missing values', min=0, max=100, value=global.param$na.filt.val),
                     radioButtons('which.test', 'Select test', choices=c('One-sample mod T', 'none'), selected='One-sample mod T'),

                     actionButton('run.test', 'Run analysis!'),
                     br(),
                     hr(),
                     actionButton('select.groups.button', 'Select groups')
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

                    #sliderInput('na.filt.val', 'Max. % missing values', min=0, max=100, value=global.param$na.filt.val),
                    
                    radioButtons('which.test', 'Select test', choices=c('One-sample mod T', 'Two-sample mod T', 'mod F', 'none'), selected=global.param$which.test),
                    actionButton('run.test', 'Run analysis!'),
                    br(),
                    hr(),
                    actionButton('select.groups.button', 'Select groups')
                )
            }

        })
        
        ## #####################################################
        ##              savePlotparams
        ## - save input parameters to reactive values
        ## 
        savePlotparams <- function(){
          ## heatmap
          global.plotparam$hm.clust <- input$hm.clust
          global.plotparam$hm.scale <- input$hm.scale
          global.plotparam$cexRow <- input$cexRow
          global.plotparam$cexCol <- input$cexCol
          global.plotparam$hm.max.val <- input$hm.max.val
          global.plotparam$hm.max <- input$hm.max
          global.plotparam$hm.show.rownames <- input$hm.show.rownames
          global.plotparam$hm.show.colnames <- input$hm.show.colnames
          ## volcano
          #global.plotparam$volc.ps <- input$volc.ps
          #global.plotparam$volc.ls <- input$volc.ls
          #global.plotparam$volc.grid <- input$volc.grid
          #global.plotparam$volc. <- input$volc.ps
          
          ## pca
          global.plotparam$pca.x <- input$pca.x
          global.plotparam$pca.y <- input$pca.y
          global.plotparam$pca.z <- input$pca.z
          ## multiscatter
          global.plotparam$ms.max <- input$ms.max
          global.plotparam$ms.min <- input$ms.min
          global.plotparam$ms.max <- input$ms.max
          global.plotparam$ms.robustify <- input$ms.robustify
          ## correlation matrix
          global.plotparam$cm.upper <- input$cm.upper
          global.plotparam$cm.lower <- input$cm.lower
          global.plotparam$cm.numb <- input$cm.numb
          ## fanplot
          global.plotparam$HC.fan.show.tip.label <- input$HC.fan.show.tip.label
          global.plotparam$HC.fan.tip.cex <- input$HC.fan.tip.cex
        }
        ## ##############################################################
        ##                 update plotting parameters
        ## ##############################################################
        updatePlotparams <- function(){
          ## heatmap
          updateSelectInput(session, inputId='hm.clust', selected=global.plotparam$hm.clust)
          updateSelectInput(session, inputId='hm.scale', selected=global.plotparam$hm.scale)
          updateNumericInput(session, inputId='cexCol', value=global.plotparam$hm.cexCol)
          updateNumericInput(session, inputId='cexRow', value=global.plotparam$hm.cexRow)
          updateNumericInput(session, inputId='hm.max.val', value=global.plotparam$hm.max.val)
          updateCheckboxInput(session, inputId='hm.max', value=global.plotparam$hm.max)
          updateCheckboxInput(session, inputId='hm.show.rownames', value=global.plotparam$hm.show.rownames)
          updateCheckboxInput(session, inputId='hm.show.colnames', value=global.plotparam$hm.show.colnames)
          
          ## PCA
          updateSelectInput(session, inputId='pca.x', selected=global.plotparam$pca.x)
          updateSelectInput(session, inputId='pca.y', selected=global.plotparam$pca.y)
          updateSelectInput(session, inputId='pca.z', selected=global.plotparam$pca.z)
          ## multiscatter
          updateCheckboxInput(session, inputId='ms.max', value=global.plotparam$ms.max)
          updateNumericInput(session, inputId='ms.max.val', value=global.plotparam$ms.max.val)
          updateNumericInput(session, inputId='ms.min.val', value=global.plotparam$ms.min.val)
          updateCheckboxInput(session, inputId='ms.robustify', value=global.plotparam$ms.robustify)
          ## correlation matrix
          updateCheckboxInput(session, inputId='cm.numb', value=global.plotparam$cm.numb)
          updateSelectInput(session, inputId='cm.upper', selected=global.plotparam$cm.upper)
          updateSelectInput(session, inputId='cm.lower', selected=global.plotparam$cm.lower)
          
          ## fanplot
          updateCheckboxInput(session, inputId='HC.fan.show.tip.label', value= global.plotparam$HC.fan.show.tip.label)
          updateNumericInput(session, inputId='HC.fan.tip.cex', value=global.plotparam$HC.fan.tip.cex)
        }
        
        ## ######################################################################
        ##
        ##               RMarkdown analysis report
        ##
        observeEvent(input$export.rmd,{
          
             
          ## label
          global.param$label <- gsub('_| |,|;|\\:|\\+|\\*', '-', input$label)
          #global.param$label <- input$label
          
          ## filename
          fn.rmd <- paste('report', global.param$label, sep='-')
          
          ## extract results  
          res = global.results$filtered
          grp <- global.param$grp
          
          ## ids to show in heatmap
          hm.rownames <- res[, 'id.concat']
          
          ## groups to compare
          grp.comp <- unique( global.param$grp.comp )
          
          ## ####################################
          ##         what to export
          export.volc <- input$export.volc & !(global.param$which.test %in% c('mod F', 'none'))
          export.cm <- input$export.cm
          export.hm <- input$export.hm
          export.pca <- input$export.pca
          export.box <- input$export.box
          
          #n.chunks <- sum( )
          
          
          #######################################
          ##    extract expression values
          res = res[, names(grp)]
          
          ## ####################################
          ##        perform pca
          if(export.pca & nrow(res) > 2){
            pca=my.prcomp2( res, grp )
            global.results$pca <- pca
          }
          
          ## ################################################
          ##        update plotting parameters 
          savePlotparams()
          

          ## ##########################################################
          ##                       header
          ## ##########################################################
          withProgress(message='Rmarkdown report', value=0, min = 0, max = 1, detail='hold on...',{
            
          rmd <- paste('<a name="top"></a>
            \n#`r global.param$label` - analysis report
            \n***
            \n')
          ## ##########################################################
          ##                   init
          ## ##########################################################
          rmd <- paste(rmd, '
                      \n### Table of contents <a name="toc"></a>
                      \n* [Summary](#summary)
                      \n     + [Data set](#dataset)
                      \n     + [Quantified features](#quantifiedfeatures)
                      \n     + [Workflow](#workflow)
                      \n     + [Test results](#testresults)
                      \n',ifelse(export.hm, '\n* [Heatmap](#heatmap)','' ),'
                      \n',ifelse(export.volc, '\n* [Volcano plots](#volcanos)','' ),'
                      \n',ifelse(export.pca, '\n* [Principle Components Analysis](#pca)','' ),'
                      \n* [QC-metrics](#qc)
                      ',ifelse(export.cm, '\n     - [Correlation matrix](#corrmat)','' ),'
                      ',ifelse(export.box, '\n    - [Box-and-whisker plots](#boxplot)','' ),'
                       ', sep='')
          #setProgress(value=0.1)
          ## ##########################################################
          ##                   data set
          ## ##########################################################
          rmd <- paste(rmd, '
                  \n***
                   \n### <font color="black">Summary</font><a name="summary"></a>
                  \n#### Data set <a name="dataset"></a>
                  \nData file: ```r global.input$file[1]```
                  \n```{r sumtab, echo=F, width="30%"}
                  \ntab <- data.frame(global.input$table.org)
                  \ngrp <- global.param$grp
                  \nN.grp <- global.param$N.grp
                  \nsum.tab <- t(data.frame(global.results$N.feat, length(grp), N.grp, ncol(global.input$table.anno)))
                  \nsum.tab <- data.frame(id=c("No. rows", "No. expression columns", "No. expression groups", "No. annotation columns"), sum.tab)
                  \ncolnames(sum.tab) <- c("Data", "Number")
                  \nkable(sum.tab, row.names=F)
                  \n```
                  \n[Back to top](#top)
                   \n
                   ', sep='\n')
          ## ##########################################################
          ##                  workflow
          ## ##########################################################
          rmd <- paste(rmd, "\n***
                       \n####  Workflow <a name='workflow'></a>
                       \n```{r workflow, echo=F}
                       \nwf.tab.ids <- c('Log scale', 'Normalization', 'Filter data', 'Test', 'Filter results')
                       \n## Reproducibility filter
                       \nif(global.param$filt.data == 'Reproducibility'){
                       \nwf.tab <- t(data.frame( global.param$log.transform, global.param$norm.data, paste( global.param$filt.data, ' (alpha=',global.param$repro.filt.val, ')',sep=''), global.param$which.test,
                       \npaste( global.param$filter.type, ' < ', global.param$filter.value)))
                       \nwf.tab <- data.frame(id=wf.tab.ids, value=wf.tab, stringsAsFactors=F)
                       \n}
                       \n## SD filter
                       \nif(global.param$filt.data == 'StdDev'){
                       \nwf.tab <- t(data.frame( global.param$log.transform, global.param$norm.data, paste( global.param$filt.data, ' (SD=',global.param$sd.filt.val, '%)',sep=''), global.param$which.test,
                       \npaste( global.param$filter.type, ' < ', global.param$filter.value) ))
                       \nwf.tab <- data.frame(id=wf.tab.ids, wf.tab, stringsAsFactors=F)
                       \n}
                       \n## no data filter
                       \nif(global.param$filt.data == 'none'){
                       \nwf.tab <- t(data.frame( global.param$log.transform, global.param$norm.data, paste( global.param$filt.data, sep=''), global.param$which.test,
                       \npaste( global.param$filter.type, ' < ', global.param$filter.value) ))
                       \nwf.tab <- data.frame(id=wf.tab.ids, wf.tab, stringsAsFactors=F)
                       \n}
                       \n## special case: no filter
                       \nif(global.param$filter.type == 'none')
                       \nwf.tab[which(wf.tab$id == 'Filter results'), 2] <- 'none'
                       \n## special case: top N
                       \nif(global.param$filter.type == 'top.n')
                       \nwf.tab[which(wf.tab$id == 'Filter results'), 2] <-  paste('top', global.param$filter.value)
                       \n
                       \n## suppress column names
                       #colnames(wf.tab) <- c('Module', 'Value')
                       \nkable(wf.tab, row.names=F)
                       \n```
                       \n[Back to top](#top)   
                       ", sep='\n')
          
          ## #######################################################
          ##                  results
          ## #######################################################
          #setProgress(0.2)
          ## ##########################################################
          ## Quantified features
          rmd <- paste(rmd, "
                      \n***
                      \n#### Quantified features <a name='quantifiedfeatures'></a>
                      \n```{r quantfeatures, echo=F, fig.width=10, fig.height=3}
                      \nif(global.param$log.transform == 'none'){
                      \n   tab <- data.frame(global.input$table.org)
                      \n} else {
                      \n    tab <- data.frame(global.results$table.log)
                      \n} 
                      \ngrp <- global.param$grp
                      \nN.grp <- global.param$N.grp
                      \ngrp.colors.legend <- global.param$grp.colors.legend
                      \ngrp.colors <- global.param$grp.colors
                       
                      \n## extract expression values
                      \ndat <- tab[, -which(colnames(tab) == global.param$id.col.value)]
                      \ndat <- data.matrix(dat)
                       
                       ## ############################################################
                       ## order columns
                      \nord.idx <- order(grp)
                      \ndat <- dat[, ord.idx]
                       
                      \ngrp.colors <- grp.colors[ord.idx]
                      \ngrp <- grp[ord.idx]
                       
                      \ngrp.colors.legend <- grp.colors.legend[order(names(grp.colors.legend))]
                      \nnames(grp.colors) <- colnames(dat)
                       
                       
                       ## #############################################################
                       ## number of non-missing values per row
                      \nna.col.idx <- apply(dat, 2, function(x) sum(!is.na(x)))
                       
                      \ndat.plot <- data.frame(x=names(na.col.idx), y=na.col.idx)
                      \ndat.plot$x <- factor(dat.plot$x, levels=dat.plot[['x']])
                       
                       ## plot
                       ##p <- plot_ly( x=names(na.col.idx), y=na.col.idx,  color=grp, colors=grp.colors.legend, type='bar')
                      \np <- plot_ly(dat.plot, x=~x, y=~y, color=grp, colors=grp.colors.legend, type='bar')
                      \np <- layout(p, title='Number of quantified features per data column', xaxis=list(title=paste('Data columns')), yaxis=list(title=paste('# quant features')))
                       
                       ##################################################
                       ##global.param$session.import.init <- F
                      \np
                       \n```
                      \n[Back to top](#top)
                       \n", sep='\n')
          
          #setProgress(0.3)
          ## ##########################################################
          ##                 Heatmap
          ## ##########################################################
          if(export.hm){
            
            
           rmd <- paste(rmd, '
                        \n***
                        \n### Heatmap <a name="heatmap"></a>
                         ', sep='\n')
            if(nrow(res) < 3){
              rmd <- paste(rmd, '\nCould not draw a heatmap: too few features passed the significance threshold.', sep='\n')
            } else {
              
             
                rmd <- paste(rmd, "
                           \n```{r heatmap, echo=F, fig.width=8, fig.height=8}
                           \nwithProgress(message='Exporting', detail='heatmap',{
                           \n# heatmap title
                           \nhm.title <- paste('filter:', global.param$filter.type, ' / cutoff:', global.param$filter.value, sep='')
                           \nhm.title <- paste(hm.title, '\nsig / total: ', nrow(res), ' / ', nrow( global.results$data$output ), sep='')
                           \n# column annotation
                           \nif(!is.null(global.input$cdesc)){
                           \n  hm.cdesc <- global.input$cdesc[global.param$cdesc.selection,  ]
                           \n} else {
                           \n   hm.cdesc <- NULL
                           \n}
                           \nif(!is.null(global.param$anno.col)){
                           \n  anno.col=global.param$anno.col
                           \n     anno.col.color=global.param$anno.col.color
                           \n} else {
                           \n  anno.col=data.frame(Group=global.param$grp)
                           \n  anno.col.color=list(Group=global.param$grp.colors.legend)
                           \n}
                           \nif(global.plotparam$hm.max){
                           \n  plotHM(res=res, hm.rownames=hm.rownames, grp=global.param$grp, grp.col=global.param$grp.colors, grp.col.legend=global.param$grp.colors.legend,  hm.clust=global.plotparam$hm.clust, hm.title=hm.title, hm.scale=global.plotparam$hm.scale , fontsize_row= global.plotparam$cexRow, fontsize_col= global.plotparam$cexCol, max.val=global.plotparam$hm.max.val, style=global.param$which.test, anno.col=anno.col, anno.col.color=anno.col.color, show.rownames=global.plotparam$hm.show.rownames, show.colnames=global.plotparam$hm.show.colnames)
                           \n} else {
                           \n  plotHM(res=res, hm.rownames=hm.rownames, grp=global.param$grp, grp.col=global.param$grp.colors, grp.col.legend=global.param$grp.colors.legend,  hm.clust=global.plotparam$hm.clust, hm.title=hm.title, hm.scale=global.plotparam$hm.scale , fontsize_row= global.plotparam$cexRow, fontsize_col= global.plotparam$cexCol, style=global.param$which.test, anno.col=anno.col, anno.col.color=anno.col.color, show.rownames=global.plotparam$hm.show.rownames, show.colnames=global.plotparam$hm.show.colnames)
                           \n}
                           \n}) # end withProgress
                           \n```   
                           \n[Back to top](#top)      
                           ", sep='\n')
                
              }
        
          }
          setProgress(0.5)
          ## ##########################################################
          ##                       volcano
          ## ##########################################################
          if(export.volc){
           
              
              rmd <- paste(rmd, "
              \n***
              \n### Volcano Plot <a name=\"volcanos\"></a>
              \n```{r volcano, echo=F, fig.width=7, fig.height=7}
              \nwithProgress(message='Exporting', detail='volcano plot',{
              \nfor(j in 1:length(grp.comp)){
              \n  local({
              \n    my_j=j
              \n    plotVolcano(grp.comp[my_j], verbose=F, cex.main=1.5, cex.axis=1, cex.leg=1)
              \n  })
              \n}
              \n})
              \n````
              \n[Back to top](#top)
              \n", sep='\n')
            
          }
          #setProgress(0.6)
          ## ############################################################
          ##                       PCA
          ## ############################################################
          if(export.pca){
            rmd <- paste(rmd, '
                        \n***
                        \n<a name="pca"></a>
                        \n### Principle Component Analysis
                       ', sep='\n')
          
            if(nrow(res) < 3){
              rmd <- paste(rmd, '\nCould not calculate principle components: Too few features passed the significance threshold.', sep='\n')
            
            } else {
            
              ## variance plot
              rmd <- paste(rmd, '\n
                      \nPCA model of a mean-centered and scaled matrix of ',length(grp), ' by ', nrow(pca$loadings), '. Number of PCs to cover 90% of the variance:', min(which((cumsum(pca$var)/pca$totalvar) > .9)), '.
                      \n```{r pca-varplot, echo=F}
                      \nplotPCAvar(pca, main="Variance explained by PCs", cex=1.5)
                      \n```
                      \n[Back to top](#top)     
                      \n', sep='')
            
              ## pc plot
              ## variance plot
              rmd <- paste(rmd, "\n
                         \n Scatterplot of ```r global.plotparam$pca.x``` and ```r global.plotparam$pca.y```.
                         \n```{r pca-scatter, echo=F}
                         \ngrp.unique <- unique(grp)
                         \ngrp.colors <- global.param$grp.colors[names(grp)]
                         \n# selected PCs
                         \npca.x <- as.numeric(sub('PC ','', input$pca.x))
                         \npca.y <- as.numeric(sub('PC ','', input$pca.y))
                         \n# build a data frame for plotly
                         \npca.mat = data.frame(
                         \n   PC1=pca$scores[, pca.x],
                         \n   PC2=pca$scores[, pca.y]
                         \n)
                         \nrownames(pca.mat) <- rownames(pca$scores)
                         \np <- plot_ly( pca.mat, type='scatter', mode='markers' )
                         \nfor(g in grp.unique){
                         \n   grp.tmp <- names(grp)[grp == g]
                         \n   p <-  add_trace(p, x=pca.mat[grp.tmp , 'PC1'], y=pca.mat[grp.tmp, 'PC2'], type='scatter', mode='markers', marker=list(size=15, color=grp.colors[grp.tmp]), text=grp.tmp, name=g  )
                         \n}
                         \np <- layout(p, title=paste('PC', pca.x,' vs. PC', pca.y, sep=''), xaxis=list(title=paste('PC', pca.x)), yaxis=list(title=paste('PC', pca.y)) )
                         \np
                          \n```
                         \n[Back to top](#top)     
                         \n", sep='')
            
            }
          } # end if export.pca
          #setProgress(0.7)
          ###############################################################
          ##            QC
          rmd <- paste(rmd, '
                      \n### QC-metrics <a name="qc"></a>
                      ', sep='')
          
          ## ##########################################################
          ##                correlation matrix
          ## ##########################################################
          if(export.cm){
            

            rmd <- paste(rmd, 
                 '\n***
              \n<a name="corrmat">
              \n#### Correlation matrix</a>
              \n<br>
              \n```{r corrmat, echo=F, warning=F, message=F, fig.width=8, fig.height=8}
              \nwithProgress(message="Exporting", detail="correlation matrix",{
              \nif(is.null(global.results$table.log)){
              \n  tab <- data.frame(global.input$table)
              \n} else{
              \n  tab <- data.frame(global.results$table.log)
              \n}
              \nid.col.value <- global.param$id.col.value
              \ngrp <- global.param$grp
              \n## group colors
              \ngrp.col <- global.param$grp.colors
              \ngrp.col.leg <- global.param$grp.colors.legend
              \n## update selection
              \ntab <- tab[, c(id.col.value, names(grp))]
              \ngrp.col <- grp.col[names(grp)]
              \ngrp.col.leg <- grp.col.leg[unique(grp)]

              \n## get correlation matrix
              \nif(is.null(global.results$cm) | global.param$update.cm == TRUE){
                   \n## calculate correlatio matrix
                   \nwithProgress(message = "Correlation matrix...",{
                   \ncm=calculateCorrMat( tab=tab,
                                         \ngrp=grp,
                                         \nlower= global.plotparam$cm.lower, upper= global.plotparam$cm.upper
                  \n)
                  \n})
                  
                  \n## store results
                  \nglobal.results$cm <- cm
                  \nglobal.param$update.cm <- FALSE
                  
                \n} else {
                \n  cm <- global.results$cm
                \n}
                \nplotCorrMat(#tab=tab,
                \n            #id.col=id.col.value,
                \n            cm=cm,
                \n            grp=grp,
                \n            grp.col.legend=grp.col.leg,
                \n  lower=input$cm.lower, upper=input$cm.upper, 
                \n  display_numbers=input$cm.numb, 
                \n  verbose=F)
              \n}) # end with progress
              \n```\n
              \n 
              \n[Back to top](#top)'          
                           , sep='\n')
           
          } # end if export.cm
          #setProgress(0.8)
          ## ##########################################################
          ##                Boxplots
          ## ##########################################################
          if(export.box){
            
            rmd <- paste(rmd, 
                         '\n***
              \n<a name="boxplot"></a>
              \n#### Box-and-whisker plots <a name="corrmat"></a>
              \nThe box plots depict the distribution of data points per sample column. 
              \n<br>
              \n```{r boxplot, echo=F, warning=F, message=F, fig.width=10}
              \nwithProgress(message="Exporting", detail="boxplots",{
              \nif(is.null(global.results$table.log)){
              \n  tab <- data.frame(global.input$table)
              \n} else{
              \n  tab <- data.frame(global.results$table.log)
              \n}
              \nid.col.value <- global.param$id.col.value
              \ngrp <- global.param$grp
              \n## group colors
              \ngrp.col <- global.param$grp.colors
              \ngrp.col.leg <- global.param$grp.colors.legend
              \n## update selection
              \ntab <- tab[, c(id.col.value, names(grp))]
              \ngrp.col <- grp.col[names(grp)]
              \ngrp.col.leg <- grp.col.leg[unique(grp)]
              \np1 <- makeBoxplotly(tab, id.col.value, grp, grp.col, verbose=F, title="No normalization")
              \nif(!is.null(global.results$table.norm) ){
              \n    # normalized ratios
              \n      tab <- data.frame(global.results$table.norm)
              \n      p2 <- makeBoxplotly(tab, id.col.value, grp, grp.col, verbose=F,title="Before and after normalization")
              \n      plotly::subplot( p1, p2, shareY=T, titleX = TRUE, titleY = TRUE)              
              \n} else{
              \n      p1
              \n}
              \n}) # end with progress
              \n```\n
              \n 
              \n[Back to top](#top)'          
                         , sep='\n')
            
          } # end if export.cm
          
          
          #setProgress(0.9)
          ## ##########################################################
          ##                    session info
          ## ##########################################################
          rmd <- paste(rmd, '
                          \n***
                          \n## Session info <a name="info"></a>
                          \n```{r sessioninfo, echo=F, width="50%"}
                          \npar.tab <- data.frame(param=c("version", "hostname", "user", "session id", "session name", "time stamp"), 
                          \n       value=c( paste(APPNAME," v",VER, sep=""), session$clientData$url_hostname, ifelse(!is.null(session$user), session$user, ""), global.param$session, global.param$label, format(Sys.time())) )   
                          \nkable(par.tab)
                          \n```
                          \n[Back to top](#top)', sep='\n')

          ## ##########################################################
          ##                  render document
          cat('## Rendering Rmarkdown file\n')
          writeLines(rmd, con=paste(global.param$session.dir, paste(fn.rmd, 'rmd', sep='.'), sep='/'))
          
          render.rmd <- try(
            rmarkdown::render(paste(global.param$session.dir, paste(fn.rmd, 'rmd', sep='.'), sep='/'))
          )
            
          }) # end withProgress
          
          if(class(render.rmd) == 'try-error'){
            #save(render.rmd, file='rmd.RData')
            
            showModal(modalDialog(
              size='m',
              title = "Problem generating R Markdown report",
              render.rmd[[1]],
              HTML('<br><br>Please install pandoc for your OS from <a href="https://github.com/jgm/pandoc/releases/tag/2.1.1
" target="_blank_">https://github.com/jgm/pandoc/releases/tag/2.1.1/</a><br>Please restart R after installing pandoc.')
            ))
            
          } else {
            global.param$rmd.name=paste(fn.rmd, 'html', sep='.')
            global.results$export.rmd <- TRUE
          }
          
          ## trigger selectize update
          global.param$update.ppi.select <- TRUE
          global.param$update.ppi.select.scat <- TRUE
          
          updateTabsetPanel(session, 'mainPage', selected='Export')
        })
        
        ####################################################
        ##         function to save current state
        ## #################################################
        saveSession <- function(){
          
          #########################################################
          ## save current plotting parameters
          savePlotparams()
          global.param$session.saved <- TRUE
          
          #########################################################
          ##          save session parameters
          global.input.imp <- reactiveValuesToList(global.input)
          global.param.imp <- reactiveValuesToList(global.param)
          global.results.imp <- reactiveValuesToList(global.results)
          global.plotparam.imp <- reactiveValuesToList(global.plotparam)
          volc.imp <-  reactiveValuesToList(volc)
          
          #################################
          ## save as R-object
          ## no label present
          if(is.null(global.param$label)){
            fn.tmp <- paste(global.param$session.dir, paste('session_', gsub('( |\\:)', '-', Sys.time()), '.RData', sep=''), sep='/')
            save(global.input.imp, global.param.imp, global.results.imp, global.plotparam.imp, volc.imp, file=fn.tmp)
          }
          ## label present
          if(!is.null(global.param$label) | nchar(global.param$label) == 0){
            fn.tmp <- paste(global.param$session.dir, paste(global.param$label, paste('_session_', gsub('( |\\:)', '-', Sys.time()), '.RData', sep=''), sep=''), sep='/')
            save(global.input.imp, global.param.imp, global.results.imp, global.plotparam.imp, volc.imp, file=fn.tmp)
          }
        }
        ## #############################################################
        ## observer
        ##            save current state as session on the server 
        ##
        ## 
        observeEvent(input$export.save.session, {
          
          ## #####################################
          ## update label
          global.param$label <- gsub('_| |,|;|\\:|\\+|\\*', '-', input$label)
          
          ## #####################################
          ## check wether a session with the same name exists
          double.check <- T
          check.session <- dir( global.param$session.dir, pattern=paste('^', global.param$label, '_session_[0-9]{4}-[0-9]{2}-[0-9]{2}-.*\\.RData$', sep=''), full.names = T)
          if(length(check.session) > 0){
            double.check <- F
            showModal(modalDialog(
              size='s',
              title = "Overwrite existing session?",
              ##fluidPage(
            ##    fluidRow(column(12, 'Overwrite session?'))
            ##  ),
            footer = NULL,
            fluidRow(    
              column(6, actionButton(inputId = 'session.overwrite' , label='Yes')),
                column(6, modalButton(label='Cancel'))
              )
            ))
            ## observe "yes"-button
            observeEvent(input$session.overwrite,{
              double.check <- T
              if(file.exists(check.session))
                file.remove(check.session)
              saveSession()
              removeModal()
            })
          }
          if(double.check)
            saveSession()
          
          ## trigger selectize update
          global.param$update.ppi.select <- TRUE
          global.param$update.ppi.select.scat <- TRUE
          
          
          updateTabsetPanel(session, 'mainPage', selected='Export')
          
        })
        
        ## ##########################################################################
        ## observer
        ##                      export GCT file
        ##
        ## ##########################################################################
        observeEvent(input$export.gct, {
          
          if(!is.null(error$msg)) return()
          
          ## update label
          global.param$label <- gsub('_| |,|;|\\:|\\+|\\*', '-', input$label)
          
          res.comb <- global.results$data$output
          grp.srt <- sort(global.param$grp)
          
          ## append annotation columns
          if(!is.null(global.input$table.anno))
            res.comb <- left_join(res.comb,  global.input$table.anno, 'id')
          
          
          ## #########################################################
          ## Two sample moderated T- test:
          ## adjust labels in the header of the Excel sheet
          ## WT.vs.KO -> KO.over.WT
          if(global.param$which.test == 'Two-sample mod T'){
            
            colnames.tmp <- colnames(res.comb)
            grp.comp <- unique(global.param$grp.comp)
            
            for(g in grp.comp){
              g.new <- strsplit(g, '\\.vs\\.') %>% unlist
              g.new <- paste(g.new[2], '.over.', g.new[1], sep='')
              colnames.tmp <- gsub(g, g.new, colnames.tmp)
            }
            colnames(res.comb) <- colnames.tmp
          }
          
          ##expDesign <- data.frame(Column=names(grp.srt), Experiment=grp.srt)
          
          fn.tmp <- sub(' ','_',
                        paste(
                          global.param$label, '_',
                          sub(' ', '_',global.param$which.test),
                          ifelse(global.param$log.transform != 'none', paste( '_', global.param$log.transform, '_', sep=''), '_'),
                          #ifelse(global.param$norm.data != 'none', paste( global.param$norm.data, '_', sep=''), '_'),
                          ifelse(input$repro.filt=='yes', paste(global.param$filt.data, sep=''), '_'),
                          sub(' .*', '', Sys.time()), sep='') 
          )
          
          ## assemble gct file
          withProgress(message='Exporting', detail='GCT file',{
            ##rdesc <- global.results$data$output
            rdesc <- res.comb
            mat <- rdesc[, names(grp.srt)] %>% data.matrix
            rdesc <- rdesc[ ,-which(colnames(rdesc) %in% names(grp.srt))]
            if(global.param$file.gct3){
              #View( global.input$cdesc)
              cdesc <- global.input$cdesc[ names(grp.srt),]
            } else {
              cdesc <- data.frame(experiment=grp.srt)
            }
            #cat(global.param$id.col.value, '\n')
            cid <- names(grp.srt)
            
            #View(rdesc)
            
            rid <- rdesc[, global.param$id.col.value] #rownames(mat)
            
            #cat('test\n')
            
            res.gct <- new('GCT')
            res.gct@mat <- mat
            res.gct@rid <- rid
            res.gct@cid <- cid
            res.gct@rdesc <- rdesc
            res.gct@cdesc <- cdesc
            
            #render.gct <- try(
              write.gct(res.gct, ofile =  sub('\\.xlsx','', paste(global.param$session.dir, fn.tmp, sep='/')) )
            #)
          })
          global.param$gct.name <- paste(fn.tmp, '_n', ncol(mat),'x',nrow(mat) ,'.gct', sep='')
          
          #if(class(render.gct) == 'try-error' | !render.gct){
          #  showModal(modalDialog(
          #    size='m',
          #    title = "Problem generating GCT file",
          #    #render.xlsx[[1]],
          #    HTML('cmapR package required.')
          #  ))
          #} else {
            global.results$export.gct <- TRUE
          #}
          
          ## trigger selectize update
          global.param$update.ppi.select <- TRUE
          global.param$update.ppi.select.scat <- TRUE
          
          updateTabsetPanel(session, 'mainPage', selected='Export')
        })
        
        
        ## ##########################################################################
        ## observer
        ##                      export Excel sheet
        ##
        ## ##########################################################################
        observeEvent(input$export.xls, {
          ##cat('User:', session$user, '\n')
          if(!is.null(error$msg)) return()
          
          ## update label
          global.param$label <- gsub('_| |,|;|\\:|\\+|\\*', '-', input$label)
          
          
          
          #########################################################
          ##               Excel sheet
          #########################################################
          withProgress(message='Exporting', detail='Excel sheet',{
              
            res.comb <- global.results$data$output
            tmp <- sort(global.param$grp)
              
            ## append annotation columns
            if(!is.null(global.input$table.anno))
              res.comb <- left_join(res.comb,  global.input$table.anno, 'id')
            
            ## #########################################################
            ## Two sample moderated T- test:
            ## adjust labels in the header of the Excel sheet
            ## WT.vs.KO -> KO.over.WT
            if(global.param$which.test == 'Two-sample mod T'){
              
              colnames.tmp <- colnames(res.comb)
              grp.comp <- unique(global.param$grp.comp)
              
              for(g in grp.comp){
                g.new <- strsplit(g, '\\.vs\\.') %>% unlist
                g.new <- paste(g.new[2], '.over.', g.new[1], sep='')
                colnames.tmp <- gsub(g, g.new, colnames.tmp)
              }
              colnames(res.comb) <- colnames.tmp
            }
            
            expDesign <- data.frame(Column=names(tmp), Experiment=tmp)
              
            ## generate_filename
            fn.tmp <- sub(' ','_',
                            paste(
                                global.param$label, '_',
                                sub(' ', '_',global.param$which.test),
                                ifelse(global.param$log.transform != 'none', paste( '_', global.param$log.transform, '_', sep=''), '_'),
                                #ifelse(global.param$norm.data != 'none', paste( global.param$norm.data, '_', sep=''), '_'),
                                ifelse(input$repro.filt=='yes', paste(global.param$filt.data, sep=''), '_'),
                                sub(' .*', '', Sys.time()),".xlsx", sep='') 
                            )
            global.param$xls.name <- fn.tmp
              
            ## Excel
            render.xlsx <- try(
              WriteXLS(c('res.comb', 'expDesign'), ExcelFileName=paste(global.param$session.dir, fn.tmp, sep='/'), FreezeRow=1, FreezeCol=1, SheetNames=c(global.param$which.test, 'class vector'), row.names=F, BoldHeaderRow=T, AutoFilter=T)
            )
          })
          
          if(class(render.xlsx) == 'try-error' | !render.xlsx){
            showModal(modalDialog(
              size='m',
              title = "Problem generating Excel sheet",
              #render.xlsx[[1]],
              HTML('Perl required to generate xlsx files. For Windows OS please install Strawberry Perl (<a href="http://strawberryperl.com/" target="_blank_">http://strawberryperl.com/</a><br>R needs to be restarted after installing Perl.)')
            ))
          } else {
            global.results$export.xls <- TRUE
          }
          
          ## trigger selectize update
          global.param$update.ppi.select <- TRUE
          global.param$update.ppi.select.scat <- TRUE
          
          updateTabsetPanel(session, 'mainPage', selected='Export')
          
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

            ## #####################################
            ## update label
            global.param$label <- gsub('_| |,|;|\\:|\\+|\\*', '-', input$label)

            ## #####################################
            ## update plot params
            savePlotparams()

            ###############################
            ## apply filter
            res = global.results$filtered

            ################################
            ## ids to show in heatmap
            hm.rownames <- res[, 'id.concat']
            
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
               
                if(is.null(global.results$table.log))
                  tab <- data.frame(global.input$table)
                else
                  tab <- data.frame(global.results$table.log)
                
                ## dataset
                #tab <- data.frame(global.results$table.log)
                
                ## id column
                id.col.value <- global.param$id.col.value
                
                #tab[, c(id.col.value, names(grp))]
                ## group vector
                grp <- global.param$grp
                ## group colors
                grp.col <- global.param$grp.colors
                grp.col.leg <- global.param$grp.colors.legend
                
                
                ## update selection
                tab <- tab[, c(id.col.value, names(grp))]
                grp.col <- grp.col[names(grp)]
                grp.col.leg <- grp.col.leg[unique(grp)]
                
                ## get correlation matrix
                if(is.null(global.results$cm) | global.param$update.cm == TRUE){
                  
                  
                  ## calculate correlatio matrix
                  withProgress(message = 'Correlation matrix...',{
                    cm=calculateCorrMat( tab=tab,
                                         grp=grp,
                                         #id.col=id.col.value,
                                         lower= global.plotparam$cm.lower, upper= global.plotparam$cm.upper
                    )
                  })
                  
                  ## store results
                  global.results$cm <- cm
                  global.param$update.cm <- FALSE
                  
                } else {
                  cm <- global.results$cm
                }
                
                ## plot
                plotCorrMat(#tab=tab,
                            #id.col=id.col.value,
                            cm=cm,
                            grp=grp,
                            grp.col.legend=grp.col.leg,
                  lower=input$cm.lower, upper=input$cm.upper, 
                  display_numbers=input$cm.numb, 
                  filename=fn.cm )

                })
            }
            ############################################################
            ##                   heatmap
            ## require at least three significant hits
            ############################################################
            if(input$export.hm){

                if(nrow(res) >= 3){
                  
                    withProgress(message='Exporting', detail='heatmap',{
                      
                    fn.hm <- paste(global.param$session.dir, 'heatmap.pdf', sep='/')
                    
                    ## heatmap title
                    hm.title <- paste('filter:', global.param$filter.type, ' / cutoff:', global.param$filter.value, sep='')
                    hm.title <- paste(hm.title, '\nsig / total: ', nrow(res), ' / ', nrow( global.results$data$output ), sep='')

                    # column annotation
                    if(!is.null(global.param$anno.col)){
                      anno.col=global.param$anno.col
                      anno.col.color=global.param$anno.col.color
                    } else {
                      anno.col=data.frame(Group=global.param$grp)
                      anno.col.color=list(Group=global.param$grp.colors.legend)
                    }
                    
                    if(input$hm.max){
                     pdf(fn.hm)
                      plotHM(res=res, 
                             hm.rownames=hm.rownames,
                             grp=global.param$grp, 
                             grp.col=global.param$grp.colors, 
                             grp.col.legend=global.param$grp.colors.legend,  
                             hm.clust=global.plotparam$hm.clust, 
                             hm.title=hm.title, 
                             hm.scale=global.plotparam$hm.scale, 
                             fontsize_row= global.plotparam$cexRow, 
                             fontsize_col= global.plotparam$cexCol, 
                             max.val=global.plotparam$hm.max.val, 
                             style=global.param$which.test, 
                             anno.col=anno.col, 
                             anno.col.color=anno.col.color, 
                             show.rownames=global.plotparam$hm.show.rownames, 
                             show.colnames=global.plotparam$hm.show.colnames)#,
                            # fn = fn.hm)
                      dev.off()
                   
                    } else {
                     pdf(fn.hm)
                      plotHM(res=res,
                             hm.rownames=hm.rownames,
                             grp=global.param$grp,
                             grp.col=global.param$grp.colors,
                             grp.col.legend=global.param$grp.colors.legend,
                             hm.clust=global.plotparam$hm.clust,
                             hm.title=hm.title,
                             hm.scale=global.plotparam$hm.scale,
                             fontsize_row= global.plotparam$cexRow,
                             fontsize_col= global.plotparam$cexCol,
                             style=global.param$which.test,
                             anno.col=anno.col,
                             anno.col.color=anno.col.color,
                             show.rownames=global.plotparam$hm.show.rownames,
                             show.colnames=global.plotparam$hm.show.colnames)#,
                             #fn = fn.hm)
                      dev.off()
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
                        
                        
                        # update selection
                        tab <- tab[, c(id.col.value, names(grp))]
                        grp.col <- grp.col[names(grp)]
                        grp.col.leg <- grp.col.leg[unique(grp)]
                        
                        
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
                    plotMultiScatter( define.max=input$ms.max, min.val=input$ms.min.val, max.val=input$ms.max.val, update='', robustify=isolate(input$ms.robustify) )
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
                    grp.srt <- sort(global.param$grp)

                    ## append annotation columns
                    if(!is.null(global.input$table.anno))
                        res.comb <- left_join(res.comb,  global.input$table.anno, 'id')
                       

                    ## #########################################################
                    ## Two sample moderated T- test:
                    ## adjust labels in the header of the Excel sheet
                    ## WT.vs.KO -> KO.over.WT
                    if(global.param$which.test == 'Two-sample mod T'){
                      
                      colnames.tmp <- colnames(res.comb)
                      grp.comp <- unique(global.param$grp.comp)
                      
                      for(g in grp.comp){
                        g.new <- strsplit(g, '\\.vs\\.') %>% unlist
                        g.new <- paste(g.new[2], '.over.', g.new[1], sep='')
                        colnames.tmp <- gsub(g, g.new, colnames.tmp)
                      }
                      colnames(res.comb) <- colnames.tmp
                    }
                    
                    expDesign <- data.frame(Column=names(grp.srt), Experiment=grp.srt)

                    fn.tmp <- sub(' ','_',
                                  paste(
                                    global.param$label, '_',
                                    sub(' ', '_',global.param$which.test),
                                    ifelse(global.param$log.transform != 'none', paste( '_', global.param$log.transform, '_', sep=''), '_'),
                                    #ifelse(global.param$norm.data != 'none', paste( global.param$norm.data, '_', sep=''), '_'),
                                    ifelse(input$repro.filt=='yes', paste(global.param$filt.data, sep=''), '_'),
                                    sub(' .*', '', Sys.time()),".xlsx", sep='') 
                    )
                    global.param$xls.name <- fn.tmp
                    
                    ## Excel
                    render.xlsx <- try(
                      WriteXLS(c('res.comb', 'expDesign'), ExcelFileName=paste(global.param$session.dir, fn.tmp, sep='/'), FreezeRow=1, FreezeCol=1, SheetNames=c(global.param$which.test, 'class vector'), row.names=F, BoldHeaderRow=T, AutoFilter=T)
                    )
                })
            }
            
            
            
            ## ###############################################################                    
            ##                       GCT v1.3
            ## export two GCT files   
            ## - input data plus test results
            ## - signed, log-transformed p-values plus row-meta data
            ## ###############################################################
            if(input$export.gct.file){
              
              res.comb <- global.results$data$output
              grp.srt <- sort(global.param$grp)
              
              ## append annotation columns
              if(!is.null(global.input$table.anno))
                res.comb <- left_join(res.comb,  global.input$table.anno, 'id')
              
              
              ## group comparisons: required for p-value GCT
              grp.comp <- unique(global.param$grp.comp)
              
              ## column names for log p valueas and fc
              ## depend on the type of test perfromed
              logp.colnames <- logfc.colnames <- grp.comp
              names(logp.colnames) <- names(logfc.colnames) <- grp.comp
              
              #############################################################
              ## moderated F: single P-value column
              ## 
              # if(global.param$which.test == 'mod F'){
              # 
              #   logp.colnames <- c('Log.P.Value')
              #   logfc.colnames <- c('')
              #   
              # }
               
              ## #########################################################
              ## Two sample moderated T- test:
              ## adjust labels in the header of the Excel sheet
              ## WT.vs.KO -> KO.over.WT
              if(global.param$which.test == 'Two-sample mod T'){
                
                colnames.tmp <- colnames(res.comb)
                
                for(g in grp.comp){
                  g.new <- strsplit(g, '\\.vs\\.') %>% unlist
                  g.new <- paste(g.new[2], '.over.', g.new[1], sep='')
                  colnames.tmp <- gsub(g, g.new, colnames.tmp)
                  
                 # grp.comp.pval.gct[g] <- g.new
                  logp.colnames[g] <- paste0('Log.P.Value.', g.new)
                  logfc.colnames[g] <- paste0('logFC.', g.new)
                }
                colnames(res.comb) <- colnames.tmp
              }
            
              ###########################################################
              ## One-sample mod T
              if(global.param$which.test == 'One-sample mod T'){
                logp.colnames <- paste0('Log.P.Value.', logp.colnames)
                logfc.colnames <- paste0('logFC.', logfc.colnames)
              }
              
              ####################
              ## filename 
              fn.tmp <- sub(' ','_',
                            paste(
                              global.param$label, '_',
                              sub(' ', '_',global.param$which.test),
                              ifelse(global.param$log.transform != 'none', paste( '_', global.param$log.transform, '_', sep=''), '_'),
                              ifelse(input$repro.filt=='yes', paste(global.param$filt.data, sep=''), '_'),
                              sub(' .*', '', Sys.time()),".xlsx", sep='') 
              )
              
              
              ################################################
              ## GCT file with transformed p-values as data
              if( !global.param$which.test == 'mod F'){
                    withProgress(message='Exporting', detail='GCT file (transformed P-values)',{
                      ##rdesc <- global.results$data$output
                      rdesc <- res.comb
                      
                      #if(global.param$which.test == 'Two-sample mod T'){
                        logp <- rdesc[, logp.colnames] %>% data.matrix
                        fc <- rdesc[, logfc.colnames ] %>% data.matrix
                      #} else {
                      #  logp <- rdesc[, paste0('Log.P.Value.', ) ] %>% data.matrix
                      #  fc <- rdesc[, paste0('logFC.', grp.comp) ] %>% data.matrix
                      #} 
                      
                      ## transformed and signed p-values
                      mat <- logp*sign(fc)
                      colnames(mat) <- paste0('signed.',colnames(logp))
                      
                      #mat <- rdesc[, names(grp.srt)] %>% data.matrix
                      
                      rdesc <- rdesc[ ,-which(colnames(rdesc) %in% names(grp.srt))]
                      #if(global.param$file.gct3){
                      #  cdesc <- global.input$cdesc[ names(grp.srt),]
                      #} else {
                      #  cdesc <- data.frame(experiment=grp.srt)
                      #}
                      cid <- colnames(mat)
                      rid <- rdesc[, global.param$id.col.value] #rownames(mat)
                      res.gct <- new('GCT')
                      res.gct@mat <- mat
                      res.gct@rid <- rid
                      res.gct@cid <- cid
                      res.gct@rdesc <- rdesc
                      #res.gct@cdesc <- cdesc
                      write.gct(res.gct, ofile =  sub('\\.xlsx','-transformed-p-val', paste(global.param$session.dir, fn.tmp, sep='/')) )
                    })
              } ## end if not mod F
              
              
              #####################################################
              ## GCT file
              
              ## assemble gct file
              withProgress(message='Exporting', detail='GCT file',{
                    ##rdesc <- global.results$data$output
                    rdesc <- res.comb
                    
                    mat <- rdesc[, paste(names(grp.srt) ) ] %>% data.matrix
                  
                    rdesc <- rdesc[ ,-which(colnames(rdesc) %in% names(grp.srt))]
                    if(global.param$file.gct3){
                      #View( global.input$cdesc)
                      cdesc <- global.input$cdesc[ names(grp.srt),]
                    } else {
                      cdesc <- data.frame(experiment=grp.srt)
                    }
                    cid <- names(grp.srt)
                    rid <- rdesc[, global.param$id.col.value] #rownames(mat)
                    res.gct <- new('GCT')
                    res.gct@mat <- mat
                    res.gct@rid <- rid
                    res.gct@cid <- cid
                    res.gct@rdesc <- rdesc
                    res.gct@cdesc <- cdesc
                    write.gct(res.gct, ofile =  sub('\\.xlsx','', paste(global.param$session.dir, fn.tmp, sep='/')) )
              })
      
                   
                
            }

            #################################
            ##       name of zip archive
            ## no label present
            if(is.null(global.param$label)){
                #fn.tmp <- paste(global.param$session.dir, paste('session_', gsub('( |\\:)', '-', Sys.time()), '.RData', sep=''), sep='/')
                #save(global.input.imp, global.param.imp, global.results.imp, global.plotparam.imp, volc.imp, file=fn.tmp)
                fn.zip <- paste( gsub('\\:', '', gsub(' ','-', gsub('-','',Sys.time()))),'.zip', sep='')
            }
            ## label present
            if(!is.null(global.param$label) | nchar(global.param$label) == 0){
                #fn.tmp <- paste(global.param$session.dir, paste(global.param$label, paste('_session_', gsub('( |\\:)', '-', Sys.time()), '.RData', sep=''), sep=''), sep='/')
                #save(global.input.imp, global.param.imp, global.results.imp, global.plotparam.imp, volc.imp, file=fn.tmp)
                fn.zip <- paste( global.param$label, '_', gsub('\\:', '', gsub(' ','-', gsub('-','',Sys.time()))),'.zip', sep='')
            }

            #########################################################
            ##   export parameters as small text file
            #########################################################
            global.param.imp <- reactiveValuesToList(global.param)
            params <- global.param.imp[c('log.transform', 'norm.data', 'filt.data',  'repro.filt.val', 'sd.filt.val', 'which.test', 'filter.type', 'filter.value')]
            #params <- global.param[[c('log.transform', 'norm.data', 'filt.data',  'repro.filt.val', 'sd.filt.val', 'which.test', 'filter.type', 'filter.value')]]
            
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
            fn.all <- grep('pdf$|xlsx$|txt$|gct$|RData$|html$',  dir(global.param$session.dir) , value=T, ignore.case=T)
	          fn.all.abs <- grep('pdf$|xlsx$|txt$|gct$|RData$html$', dir(global.param$session.dir, full.names=T, ignore.case=T), value=T)

            ## #################################################
            ## handle special characters in file names
            ## #################################################

            ## for windows: "
            fn.all.abs <- paste('"',fn.all.abs,'"', sep='')
            ## unix/linux: '
            fn.all <- paste('\'',fn.all,'\'', sep='')

            ###################################################
            ## run command
	          if(OS == 'Windows') {
                #system( paste('zip -0 ', paste(global.param$session.dir, fn.zip, sep='/'), ' ', paste(fn.all.abs, collapse=' '),sep='') )
	              system( paste('bin/zip.exe -0 ', paste(global.param$session.dir, fn.zip, sep='/'), ' ', paste(fn.all.abs, collapse=' '),sep='') )
            } else {
                system( paste('cd ', global.param$session.dir,' && zip -0 ', fn.zip, ' ', paste(fn.all, collapse=' '), sep='') )
            }
            ## store file.name
	          global.param$zip.name=fn.zip
            ## cat('test16\n')
	          
            #----------------------------------------------------------------------
            #             clean up
	          # remove archived files: all files except RData
	          rdata.idx <- grep('\\.RData[\\"|\']$', fn.all.abs)
	          fn.rm <- fn.all.abs
	          if(length(rdata.idx) > 0)
	            fn.rm <- fn.rm[-rdata.idx]
	          
	          #fn.rm <- gsub('"|\'', '', fn.all.abs[]) 
	          fn.rm <- gsub('"|\'', '', fn.rm)
	          
	        #  cat('ALL FILES:', fn.all.abs, '\n')
	        #  cat('CLEAN UP:', fn.rm, '\n')
	          
            file.remove(fn.rm)
            
            # keep only the latest zip file
            all.zip <- dir(global.param$session.dir, pattern = '\\.zip$', full.names = T)
            if(length(all.zip) > 1){
              all.zip <- all.zip[ -which.max(file.info(all.zip)$ctime) ]
              file.remove(all.zip)
            }
            # 
  
            ## trigger selectize update
            global.param$update.ppi.select <- TRUE
            global.param$update.ppi.select.scat <- TRUE
            
            
             ###############################################################
             ## flag export
             global.results$export.results=TRUE

             ## redirect to the same panel
             updateTabsetPanel(session, 'mainPage', selected='Export')
        })


        ###################################################################################
        ##
        ##                             Do the actual computation
        ##
        ###################################################################################
        
        observeEvent(input$ms.robustify,{
          updateCheckboxInput('ms.max', session = session, value=!input$ms.robustify)
        })
        observeEvent(input$ms.max,{
          updateCheckboxInput('ms.robustify', session=session, value=!input$ms.max)
        })
        
        
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
cat('id: ', global.param$id.col.value, '\n')
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
            
            ## ##############################################
            ##                      GCTX
            if(grepl('\\.gctx$', fn)){
              
                gct <- try( parse.gctx(fn) )

            # ###############################################
            #                     GCT 1.2
            } else if( length( grep( '^\\#1\\.2', readLines(fn,n=1))) > 0){
            

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
              #gct <- parse.gctx(fn)
              gct <- try( parse.gctx2(fn) )
              #gct <- try( parse.gctx(fn) )
              #cat('test2\n')

              if(class(gct) == 'try-error'){
                #error$msg <- paste('<p>Error importing GCT 1.3 file:<br>', gct[1],'<p>')
                error$title <- "Error importing GCT 1.3 file"
                error$msg <- gct[1]
                
                validate(need(class(gct) != 'try-error', 'Error importing GCT 1.3 file.'))
              }
              
              ## #################################################
              ## robustify ids
              #gct@rid <- make.unique(make.names(gct@rid))
              gct@rid <- make.unique( gct@rid )
              rownames(gct@mat) <- gct@rid
              if(nrow(gct@rdesc) > 0){
                rownames(gct@rdesc) <- gct@rid   
              }
              
              gct@cid <- make.unique(make.names(gct@cid))
              if(nrow(gct@cdesc) > 0)
                rownames(gct@cdesc) <- gct@cid
              colnames(gct@mat) <- gct@cid
              
              ## error checking for column meta data
              if(nrow(gct@cdesc) == 0){
                error$title <- "Error parsing GCT 1.3 column meta data"
                error$msg <- paste('No column meta data tracks defined! Need at least one column meta data track to use as class vector.') 
                validate( need(nrow(gct@cdesc) > 0, 'Error parsing GCT 1.3 column meta data.'))
              }
              
              ## remove 'id'
              if('id' %in% colnames(gct@cdesc)){
                cn.tmp <- colnames(gct@cdesc)
                rm.idx <- which(colnames(gct@cdesc) == 'id')
                #gct@cdesc <- data.frame(gct@cdesc[ ,-which(colnames(gct@cdesc) == 'id') ] )
                gct@cdesc <- data.frame(gct@cdesc[ ,-rm.idx] )
                cn.tmp <- cn.tmp[-rm.idx]
                colnames(gct@cdesc) <- cn.tmp
              }
              if(ncol(gct@cdesc) == 1)
                rownames(gct@cdesc) <- gct@cid
              
              
              
              
              # expression table
              if(nrow(gct@rdesc) > 0 ){
                tab <- data.frame(id=gct@rid, gct@rdesc, gct@mat, stringsAsFactors = F)
              } else {
                tab <- data.frame(id=gct@rid, gct@mat, stringsAsFactors = F)
              }
              rownames(tab) <-tab$id
              
              ## sample names
              colnames.tmp <- chopString(colnames(tab), STRLENGTH)
              names(colnames.tmp) <- colnames(tab)
              
               
              # id column 
              global.param$id.col.value='id'
              global.param$id.done=T
              
              # robustify ids
              #ids <- make.unique( as.character(tab[, global.param$id.col.value] ), sep='_')
              
              ## replace values in id column
              #tab[, global.param$id.col.value] <- ids
              ## use id as rownames
              #rownames(tab) <- ids
              
              ## #################
              ## map to gene names
              map.res <- mapIDs(tab$id)
              global.results$keytype <- map.res$keytype
              global.results$id.map <- map.res$id.map
             
              ## store values
              global.input$table <- global.input$table.org <- tab
              global.input$file <- input$file
              global.input$table.colnames <- colnames.tmp
              
              # meta data
              #rdesc <- data.frame(gct@rdesc)
              #if(ncol(rdesc) == 1){
              # if(length(rdesc) > 0)
              #   rownames(rdesc) <- gct@rid
              #}
              #cdesc <- data.frame(gct@cdesc)
              
              global.input$rdesc <- gct@rdesc
              global.input$cdesc <- gct@cdesc

              
              #rownames(global.input$cdesc) <- make.names(make.unique( rownames(global.input$cdesc))) # convert to proper names
              
              global.param$cdesc.all <- global.param$cdesc.selection <- colnames(global.input$cdesc)
              
               
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
           
            ## ###################################            
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
            
            ## #####################################################
            ## backwards compatibility II
            if(is.null(names(global.param$grp.colors)))
              names(global.param$grp.colors) <- names(global.param$grp.colors)
            if(is.null(global.param$session.saved))
              global.param$session.saved <- T
            if(is.null( global.param$update.pca) )
              global.param$update.pca <- TRUE
            if(is.null( global.param$update.cm) )
              global.param$update.cm <- TRUE
          

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

            ## ###############################################################
            ## update plotting paramaters
            updatePlotparams()
            
            ##################################
            ## set flags
            global.param$file.done=T
            global.param$id.done=T
            global.param$session.imported=T
            global.param$analysis.run=T

            global.param$session.import.init=T

            global.results$export.rmd=F
            global.results$export.results=F
            global.results$export.xls=F
            
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
                error$title <- "Problem parsing experimental design file."
                error$msg <- 'Experimental design file does not match the table you have uploaded!'
                return()
            }

            ## not an experimental design file
            if( sum( colnames(grp.file) %in% c('Column.Name', 'Experiment'), na.rm=T) != 2 )  {
                error$title <- "Problem parsing experimental design file."
                error$msg <- 'This is not an experimental desgin file! The file should contain two columns (Column.Name, Experiment)!'
                return()
            }
            ## 'empty' file
            if( sum( nchar(Experiment) > 0, na.rm=T ) == 0 | sum(!is.na( Experiment) == 0) ){
                error$title <- "Problem parsing experimental design file."
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
            tab <- global.input$table[ , c(global.param$id.col.value, names(grp))]
            
            ## #################################
            ## remove NA rows
            na.row.idx <- apply(tab[, names(grp)], 1, function(x) sum(is.na(x))/length(x) )
            na.row.idx <- which(na.row.idx == 1)
            if(length(na.row.idx) > 0){
              tab <- tab[-na.row.idx, ]
              global.input$NA.rows <- length(na.row.idx)
            }
            
            #global.input$table <- global.input$table[ , c(global.param$id.col.value, names(grp))]
            global.input$table <- tab
            
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
        #observeEvent(global.param$grp.done, {
          
          test <- input$which.test
          
          ## #############################################
          ## specify which comparisons should be performed
          if(test %in% c('One-sample mod T', 'mod F', 'none')){
            
            ## each group separetely
            groups.comp <- global.param$grp.all
            
            
            global.param$grp.comp.all <- groups.comp
            
            #if(global.param$run.test == 0) {
            if( is.null( global.param$grp.comp.selection ) | sum( global.param$grp.comp.selection %in% global.param$grp.comp.all) == 0 )
               global.param$grp.comp.selection <- groups.comp
            if( is.null( global.param$grp.selection)  | sum( global.param$grp.selection %in% global.param$grp.comp.all) == 0 )   
              global.param$grp.selection <- groups.comp
            
          }
          ## #############################################
          ## pairwise combinations for 2-sample test
          if(test == 'Two-sample mod T'){
            
            ## all pairwise combinations
            groups.unique <- unique(global.param$grp.all)
            #groups.unique <- unique(global.param$grp.selection)
            
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
          }
            #if( global.param$run.test == 0  ){
          global.param$grp.comp.all <- groups.comp
            #}
          #if( global.param$run.test == 0  )  
          if(is.null(global.param$grp.comp.selection) | sum( global.param$grp.comp.selection %in% global.param$grp.comp.all) == 0 )
            global.param$grp.comp.selection <- groups.comp
            #global.param$grp.selection <- 
          #}
          
          
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

                      error$title <- "'Error applying two-component normalization."
                      error$msg <- paste('Could not fit mixture model for data column: ', tab, ' Giving up...')

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
                res.comb <- left_join(res.comb, global.results$id.map, 'id')
              
            }

            rownames(res.comb) <- make.names(res.comb$id)

            ## #########################################
            ## store the results
            global.results$data$output <- res.comb
            
            ## number of features with at least one non-missing value
            global.results$N.feat <- sum(apply( tab[, names(groups)], 1, function(x) sum(is.na(x)/length(x))) < 1)

            ##cat(dim(global.results$data$output))

            ## #####################################
            ## set some flags
            global.param$analysis.run <- T
            
            #global.results$export.results <- F
            global.results$export.results <- F
            global.results$export.xls <- F
            global.results$export.rmd <- F
            
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
            
            ## #################################
            ## flag to update PCA
            global.param$update.pca <- TRUE
            ## flag to update correlation matrix
            global.param$update.cm <- TRUE
        })


        ###################################################################################
        ## observer
        ##
        ##                   filter the test results accross multiple groups
        ##
        ##
        ###################################################################################
        observeEvent(c(input$filter.type, input$filter.value.top.n, input$filter.value.nom.p, input$filter.value.adj.p, global.param$run.test),{

            ## only run after one analysis has been completed
            if(global.param$run.test == 0 | is.null(input$filter.type)) return()

            ## ######################################
            ## get the current tab
            tab.select <- input$mainPage

            groups.comp=unique(global.param$grp.comp)

            ## test results
            res <- data.frame(global.results$data$output, stringsAsFactors=F )

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
                        res <- res[ order( unlist(apply( res[, grep('^P.Value', colnames(res) )], 1, min, na.rm=T))  ), ]
                        res <- res[ 1:input$filter.value.top.n, ]

                        ## order according to FC
                        #res <- res[order( unlist(apply( res[, grep('^logFC', colnames(res) )], 1, min, na.rm=T))  ), ]

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

                        #res <- res[ which( unlist(apply( res[, grep('^P.Value', colnames(res) )], 1, function(x) sum(x < input$filter.value.nom.p))) > 0), ]
                      # keep all features that are significant in at least one group comparison
                      res <- res[ which( unlist(apply( res[, grep('^P.Value', colnames(res) )], 1, 
                                                       function(x) sum(x < input$filter.value.nom.p, na.rm=T))
                                                ) > 0), ]
                      #res <- res[order( unlist(apply( res[, grep('^logFC', colnames(res) )], 1, min, na.rm=T))  ), ]
                      
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
                    res <-  res[order( res[, paste('logFC.', groups.comp, sep="") ]), ]
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

                        #res <- res[ which( unlist(apply( res[, grep('^adj.P.Val', colnames(res) )], 1, function(x) sum(as.numeric(x) < input$filter.value.adj.p ))) > 0), ]
                      res <- res[ which( unlist(apply( res[, grep('^adj.P.Val', colnames(res) )], 1, 
                                                       #function(x) sum(as.numeric(x) < input$filter.value.adj.p, na.rm=T ))
                                                       #function(x) sum(as.numeric(x) < input$filter.value.adj.p, na.rm=T ))
                                                       function(x) sum( x < input$filter.value.adj.p, na.rm=T ))
                                                ) > 0), ]
                      
                      #res <- res[order( unlist(apply( res[, grep('^logFC', colnames(res) )], 1, min,  na.rm=T))  ), ]
                        
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

            cat('\n-- filter.res --\n')
            cat('filter.type: ', global.param$filter.type, '\nfilter.value:', global.param$filter.value, '\n')
            
            
            ###################################################
            ## global filter accross all experiments
            ##rownames(res) <- res[, 'id.concat']

            global.results$filtered <- res
            ##View(res)
            global.results$filtered.groups <- res.groups
            
            N.feat.filt <- sum( apply(res[, names(groups)], 1, function(x) sum(is.na(x))/length(x) ) < 1, na.rm=T)
            
            global.results$N.feat.filtered <- N.feat.filt


            ## #####################################################
            ## suppress switching to the first tab
            updateNavbarPage(session, 'mainPage', selected=tab.select)

            ## trigger selectize update
            global.param$update.ppi.select <- TRUE
            global.param$update.ppi.select.scat <- TRUE
            
            ## flag to update PCA
            global.param$update.pca <- TRUE
            ## flag to update correlation matrix
            global.param$update.cm <- TRUE
            
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

        ## download handler for Rmarkdown html file
        output$download.rmd <- downloadHandler(
          
          filename = function(){paste('results', global.param$rmd.name, sep='_')},
          content = function(file){
            file.copy( paste(global.param$session.dir, global.param$rmd.name, sep='/'), file)
          }, contentType = "application/html"
        )
        ## download handler for Excel file
        output$download.xls <- downloadHandler(
          
          filename = function(){paste('results', global.param$xls.name, sep='_')},
          content = function(file){
            file.copy( paste(global.param$session.dir, global.param$xls.name, sep='/'), file)
          }, contentType = "application/xlsx"
        )
        ## download handler for GCT file
        output$download.gct <- downloadHandler(
          filename = function(){paste('results', global.param$gct.name, sep='_')},
          content = function(file){
            file.copy( paste(global.param$session.dir, global.param$gct.name, sep='/'), file)
          }, contentType = "application/gct"
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
            
            grp <- global.param$grp
            N.grp <- global.param$N.grp

           # N.feat <- sum(apply(tab, 1, function(x) sum(is.na(x)/length(x))) < 1)
          #  global.input$N.feat <- N.feat
            
            ##cat(nrow(tab), ' - ', length(grp), ' - ', N.grp)

            sum.tab <- t(data.frame(
                ##id=c('Rows', 'Expression columns', 'Groups'),
                #N.rows=nrow(tab),
                global.results$N.feat,
                N.columns=length(grp),
                N.groups=N.grp
            ))
            
            if(!is.null(global.input$NA.row)){
              sum.tab[1,1] <- paste(sum.tab[1,1], '(no quant: ',global.input$NA.row,')')
            }
            
            ##View(sum.tab)
            sum.tab <- data.frame(id=c('No. features', 'No. expression columns', 'No. groups'), sum.tab)
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
        ###################################################################################
        ##
        ##                upset plot comparing significant hits
        ##
        ###################################################################################
        output$summary.upset.up <- renderPlot({
          
          if(is.null(global.results$data)) return()
          if(!is.null(error$msg)) return()
          
          ## tested groups
          grp.comp=unique(global.param$grp.comp)
          
          # results 
          res <- global.results$data$output
          
          
          ## extract filter type
          filter.type=global.param$filter.type
          filter.value=global.param$filter.value
          
          
          
          validate( need( global.param$which.test %in% c('One-sample mod T', 'Two-sample mod T'), paste('Only available for multi-group one-or two-sample moderated T-tests!')) )
          validate( need( length(grp.comp) > 1, paste('Need at least 2 groups!')) )
          validate( need( nrow(global.results$filtered) > 0 , 'No significant features to compare!') )
          validate( need( filter.type %in% c('adj.p', 'nom.p'), paste('No significance filter defined!')) )
          
         
          
          if(filter.type == 'adj.p')
            test.tab <- res[, grep( paste(paste('^adj.P.Val', grp.comp, sep='.'), collapse='|' ), colnames(res) ) ] %>% data.matrix()

          if(filter.type == 'nom.p')
            test.tab <- res[, grep( paste(paste('^P.Value', grp.comp, sep='.'), collapse='|' ), colnames(res)) ]  %>% data.matrix()
          
          # fold change
          fc.tab  <- res[, grep( paste(paste('^logFC', grp.comp, sep='.'), collapse='|' ), colnames(res)) ]  %>% data.matrix()
           
          
          ## binary 'upset' matrix 
          na.idx <- is.na(test.tab)
          ## up-regulated
          sig.up.idx <- test.tab < filter.value & !is.na(test.tab) & fc.tab > 0 & !is.na(fc.tab)
          #notsig.up.idx <- test.tab >= filter.value & !is.na(test.tab) # & fc.tab > 0 & !is.na(fc.tab)
          
          #test.tab.up <- test.tab
          test.tab.up <- matrix(0, nrow=nrow(test.tab), ncol=ncol(test.tab), dimnames=dimnames(test.tab))
          colnames(test.tab.up) <- sub('^adj.P.Val\\.|^P.Value\\.', '', colnames(test.tab.up))
          
          #test.tab.up[fc.tab < 0 & !is.na(fc.tab)] <- 0
          test.tab.up[sig.up.idx] <- 1
          #test.tab.up[notsig.up.idx] <- 0
          #test.tab.up[na.idx] <- 0
          
         
          #save(test.tab, test.tab.dn, fc.tab, sig.dn.idx, filter.value, file='upset.RData')
          
          ## need at least two groups with significant features
          tmp <- apply(test.tab.up, 2, sum)
          tmp.idx <- which(tmp > 0)
          
          validate( need( length(tmp.idx) > 1, paste('Nothing to show here!\n\nOnly one group with significant features (', colnames(test.tab.up)[tmp.idx], ').') ) )
          
          test.tab.up <- test.tab.up[, tmp.idx]
          
          ## colors one-sample -> use color vector of different groups
          ## colors two-sample -> no colors
          #sets.bar.color <- rep('grey', ncol(test.tab.up))
          sets.bar.color <- 'grey'
          
          #save(test.tab.up, sets.bar.color, fc.tab, filter.value, file='upset.RData')
          
          
          ## plot
          #par(mfrow=c(1,2))
          upset(data.frame( test.tab.up ), 
                order.by='degree', 
                nintersects=NA,
                text.scale=c(2,2,1.5, 1.2, 2, 2),
                point.size=6, main.bar.color='darkblue',
                sets.bar.color=sets.bar.color
                
                )
          
     
        })
        ## upset plot: down regulated
        output$summary.upset.dn <- renderPlot({
          
          if(is.null(global.results$data)) return()
          if(!is.null(error$msg)) return()
          
          ## tested groups
          grp.comp=unique(global.param$grp.comp)
          
          # results 
          res <- global.results$data$output
          
          
          ## extract filter type
          filter.type=global.param$filter.type
          filter.value=global.param$filter.value
          
          
          validate( need( global.param$which.test %in% c('One-sample mod T', 'Two-sample mod T'), paste('Only available for multi-group one-or two-sample moderated T-tests.')) )
          validate( need( length(grp.comp) > 1, paste('Need at least 2 groups.')) )
          validate( need( nrow(global.results$filtered) > 0 , 'No significant features to compare!') )
          validate( need( filter.type %in% c('adj.p', 'nom.p'), paste('No significance filter defined!')) )
          
        
          if(filter.type == 'adj.p')
            test.tab <- res[, grep( paste(paste('^adj.P.Val', grp.comp, sep='.'), collapse='|' ), colnames(res) ) ] %>% data.matrix()
          
          if(filter.type == 'nom.p')
            test.tab <- res[, grep( paste(paste('^P.Value', grp.comp, sep='.'), collapse='|' ), colnames(res)) ]  %>% data.matrix()
          
          # fold change
          fc.tab  <- res[, grep( paste(paste('^logFC', grp.comp, sep='.'), collapse='|' ), colnames(res)) ]  %>% data.matrix()
          
          
          ## binary 'upset' matrix 
          na.idx <- is.na(test.tab)
        
          ## down regulated
          sig.dn.idx <- test.tab < filter.value & !is.na(test.tab) & fc.tab < 0 & !is.na(fc.tab)
          #notsig.dn.idx <- test.tab >= filter.value & !is.na(test.tab)# & fc.tab < 0 & !is.na(fc.tab)
          
          
          
          test.tab.dn <- matrix(0, nrow=nrow(test.tab), ncol=ncol(test.tab), dimnames=dimnames(test.tab))
          colnames(test.tab.dn) <- sub('^adj.P.Val\\.|^P.Value\\.', '', colnames(test.tab.dn))
          test.tab.dn[sig.dn.idx] <- 1
         
          ## need at least two groups with significnat features
          tmp <- apply(test.tab.dn, 2, sum)
          tmp.idx <- which(tmp > 0)
          
          validate( need( length(tmp.idx) > 1, paste( 'Nothing to show here!\n\nOnly one group with significant features (', colnames(test.tab.dn)[tmp.idx], ').') ) )
          
          test.tab.dn <- test.tab.dn[, tmp.idx]
          #save(test.tab, test.tab.dn, fc.tab, sig.dn.idx, filter.value, file='upset.RData')
          
          ## colors one-sample -> use color vector of different groups
          ## colors two-sample -> no colors
          #sets.bar.color <- rep('grey', ncol(test.tab.dn))
          sets.bar.color <- 'grey'
          
          ## plot
          upset(data.frame( test.tab.dn), 
                          order.by='degree', 
                          nintersects=NA,
                          text.scale=c(2,2,1.5, 1.2, 2, 2),
                          point.size=6, main.bar.color='darkblue',
                          sets.bar.color=sets.bar.color
          )
          
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
            
            p <- plot_ly( x=as.numeric(names(n.miss))/ncol(dat)*100, y=n.miss, type='scatter', mode='lines+markers', marker=list(color = 'black'), line=list(color='black') )
            p <- layout(p, title=paste('Fully quantified features:', n.miss[1]), xaxis=list(title=paste('Max. percent missing'), showspikes=T, spikecolor='grey60'), yaxis=list(title=paste('# quant. features'), showspikes=T, spikecolor='grey60'))
            
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
            #dat <- tab[, -which(colnames(tab) == global.param$id.col.value)]
            dat <- tab[, names(grp)]
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
            p <- layout(p, title='Number of quantified features per sample column', xaxis=list(title=paste('Sample columns')), yaxis=list(title=paste('# quant. features')))

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
            
            ## remove all NA features
            keep.idx <- which(apply(tab[, names(global.param$grp)], 1, function(x) sum(is.na(x))/length(x) ) < 1)
            tab <- tab[keep.idx, ]
            
            ## column names
            colnames(tab) <- sub('^X','',colnames(tab))
            
            ## append annotation columns
            table.anno <- global.input$table.anno


            ## check whether there is annotation stored
            if(!is.null(table.anno)){
                if(is.null(dim(table.anno)))
                    table.anno <- data.frame(table.anno)

                tab <- left_join(tab, table.anno, 'id')
            }

            
            
            if(nrow(tab) > 0){
              
                ## add links to GeneCard/UniProt
                up.id <- tab$id
                up.link <- link.db(up.id, global.results$keytype)
                tab[, 'id'] <- up.link
            }
            
            #View(tab)
            
            datatable(tab, style = 'bootstrap', 
                      width = 1000, escape = F,  filter='top', 
                      rownames=F, options = list( pageLength = 10, scrollX = T, selection='none',
                                                  autoWidth = TRUE#,
                                                  #fnRowCallback = JS("function( nRow, aData, iDisplayIndex, iDisplayIndexFull ) {ind = 2; $('td:eq('+ind+')', nRow).html( (aData[ind]).toFixed(2) );}")
                                                  )
                      )
            #tab
        #}, options = list( pageLength = 10, scrollX = T), escape=F, filter='top', rownames=F, server = T)
        }, server=T)

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
            #cat('now:', global.param$update.ppi.select, '\n')
            if( !global.param$update.ppi.select ) return()
            #cat('now:', global.param$update.ppi.select, '\n')

            grp.comp <- unique( global.param$grp.comp )

            for(i in 1:length(grp.comp)){

                ## commenting out the if-statemnet made select input work after multiple rounds of analysis
                ##if(!is.null(input[[ gsub('\\.','', paste0('ppi.bait.',  grp.comp[i])) ]] )){

                ## #################################
                ## all ids found in data
                choices=c( unlist(global.results$id.map$id.concat) )


                ## server-side rendering of selectizeInput
                 #updateSelectizeInput( session=session, inputId = gsub('\\.','',paste0('ppi.bait.',  grp.comp[i])), label='Bait protein', choices=choices, selected=NULL, server=T)#,
                updateSelectizeInput( session=session, inputId = gsub('\\.','',paste0('ppi.bait.',  grp.comp[i])), choices=choices, selected="", server=T)#,
                                       # options=list( 
                                      #    placeholder = 'test'#,
                                        #  oninit = I('function() { this.setValue(""); }')
                                      #  ))
                 #updateSelectizeInput( session=session, inputId = gsub('\\.','',paste0('ppi.bait.',  grp.comp[i])), label='Bait protein', choices=choices, selected=NULL, server=T, 
                #                       options= list(
                #                          onInitialize = I('function() { this.setValue(""); }')
                #                       )
                 #)
                # 
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
                updateSelectizeInput( session=session, inputId = gsub('\\.','',paste0('ppi.bait.scat.',  grp.comp[i])), choices=choices, selected="", server=T)
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
            #grp.comp <- unique( global.param$grp.comp )
            grp.comp <- unique( global.param$grp.comp.selection )
              
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
                      #save(text.tmp, file='volc.RData')
                      
                      #cat('yes --- ', text.tmp[[1]], '\n')
                      
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
                      volc[[paste('text', grp.comp[my_i], sep='.')]] <- NULL
                      volc[[paste('P.Value', grp.comp[my_i], sep='.')]] <- NULL
                      volc[[paste('adj.P.Val', grp.comp[my_i], sep='.')]] <- NULL
                  })
                    
                  
                  ## ######################################################
                  ##       remove selected rows from table and volcano labels  
                  observeEvent(input[[paste('volc.tab.reset.select', grp.comp[my_i], sep='.')]],{
                    
                    if(is.null(input[[ paste( paste('volc.tab.selected', grp.comp[my_i], sep='.'), 'rows_selected', sep='_' )]])) return()
                    
                      ## if so remove from the list
                      #ids <- volc[[paste('text', grp.comp[my_i], sep='.')]]
                      
                      #idx <- c()
                      #for(l in 1:length(ids))
                      #  idx[l] = input[[ paste('rm', ids[l], grp.comp[my_i], sep='.') ]] 
                      
                      idx=input[[ paste( paste('volc.tab.selected', grp.comp[my_i], sep='.'), 'rows_selected', sep='_' )]]
                      #cat('test', idx, '---\n')
                      volc[[paste('x', grp.comp[my_i], sep='.')]] <- volc[[paste('x', grp.comp[my_i], sep='.')]][-idx]
                      volc[[paste('y', grp.comp[my_i], sep='.')]] <- volc[[paste('y', grp.comp[my_i], sep='.')]][-idx]
                      volc[[paste('text', grp.comp[my_i], sep='.')]] <- volc[[paste('text', grp.comp[my_i], sep='.')]][-idx]
                      volc[[paste('xy', grp.comp[my_i], sep='.')]] <- volc[[paste('xy', grp.comp[my_i], sep='.')]][-idx]
                      volc[[paste('P.Value', grp.comp[my_i], sep='.')]] <- volc[[paste('P.Value', grp.comp[my_i], sep='.')]][-idx]
                      volc[[paste('adj.P.Val', grp.comp[my_i], sep='.')]] <- volc[[paste('adj.P.Val', grp.comp[my_i], sep='.')]][-idx]
                      
                      })
                    
                    
                  ## ####################################################
                  ## table of selected features
                  output[[paste('volc.tab.selected', grp.comp[my_i], sep='.')]] <- DT::renderDataTable({

                      if(is.null(volc[[paste('x', grp.comp[my_i], sep='.')]])) return()
                      if(length(volc[[paste('x', grp.comp[my_i], sep='.')]]) == 0) return()

                      tags$h4('Selection:')

                      id.tmp <- volc[[paste('text', grp.comp[my_i], sep='.')]]

                      dat.select = data.frame(id=unlist(volc[[paste('text', grp.comp[my_i], sep='.')]]), logFC=round( unlist(volc[[paste('x', grp.comp[my_i], sep='.')]]), 2), P.Value=round( unlist(volc[[paste('P.Value', grp.comp[my_i], sep='.')]]),3), adj.P.Value=round( unlist(volc[[paste('adj.P.Val', grp.comp[my_i], sep='.')]]), 3) )
                      up.id <- dat.select[, 'id']
                      up.link <- link.db(up.id, global.results$keytype)
                      dat.select[, 'id'] <- up.link

                      #rm.button <- c()
                      #for(l in 1:length(up.id)){
                      #  rm.button[l] <- sprintf(
                      #  "<input type=\"checkbox\" name=\"%s\" value=''/>",
                      #  paste('rm', up.id[l], grp.comp[my_i], sep='.')
                      #)
                      #}
                      #dat.select <- data.frame(Remove=rm.button, dat.select)
                      #colnames(dat.select)[1] <- ""
                      datatable(dat.select, escape=F, rownames = F, filter = 'none',autoHideNavigation = T,
                                options = list( 
                                  searching=F, paging=F, selection='multiple'#, dom='t', server=FALSE
                                  #,callback = JS("table.rows().every(function(i, tab, row) {
                                  #                var $this = $(this.node());
                                  #              $this.attr('id', this.data()[0]);
                                  #              $this.addClass('shiny-input-radiogroup');
                                  #            });
                                  #            Shiny.unbindAll(table.table().node());
                                  #            Shiny.bindAll(table.table().node());")
                                
                                #                fnRowCallback = I("function( nRow, aData, iDisplayIndex, iDisplayIndexFull ) {ind = 2; $('td:eq('+ind+')', nRow).html( (aData[ind]).toFixed(2) );}")
                                  )
                                )
                      #dat.select
                  } )#, sanitize.text.function = function(x) x)

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

            ## #################################
            ## labels
            non.sig.txt <- 'not signif.'
            sig.txt <- 'signif.'
            
            ## ########################################################
            ##                  add PPI stuff
            ## ##########################################################
            ppi.idx <- c()

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
                ## if NO interactors have been found
                if(length(ppi.int.idx) == 0){
                  
                  non.int.idx <- 1:nrow(dat.plot)
                  if(length(sig.idx) > 0)
                    non.int.idx <- non.int.idx[ -c(sig.idx) ]  
                  
                  ## set up plot  
                  p <- plot_ly(x=dat.plot$x.ax[ non.int.idx ], y=dat.plot$y.ax[ non.int.idx], type='scatter', mode='markers', marker=list(size=10, color=col[non.int.idx]), text=IDs[non.int.idx], opacity=.2, showlegend=T, name=non.sig.txt )
                  
                  ## add significant proteins
                  if(length(sig.idx) > 0){
                    p <- p %>% add_trace(x=dat.plot$x.ax[ sig.idx ], y=dat.plot$y.ax[ sig.idx], type='scatter', mode='markers', marker=list(size=10, color=col[sig.idx]), text=IDs[sig.idx], opacity=.2, showlegend=T, name=sig.txt  )
                  }
                  
                  ## #######################################
                  ## if interactors have been found
                } else {
                    col[ppi.int.idx] <- ppi.col[ nchar(ppi.col) > 0 ]
                    #ppi.idx <- ppi.int.idx
  
                    
                    #cat(ppi.int.idx, '\n')
                    ## ###################################################
                    ## plot non-interactors
                    #non.int.idx <- 1:nrow(dat.plot)[ -c(ppi.int.idx) ]
                    non.int.idx <- setdiff( 1:nrow(dat.plot), ppi.int.idx)
                    
                   # cat(non.int.idx, '\n')
                    
                    if(length(sig.idx) > 0)
                      non.int.idx <- non.int.idx[ -c(sig.idx) ]  
                    
                    p <- plot_ly(x=dat.plot$x.ax[ non.int.idx ], y=dat.plot$y.ax[ non.int.idx], type='scatter', mode='markers', marker=list(size=10, color=col[non.int.idx]), text=IDs[non.int.idx], opacity=.2, showlegend=T, name=non.sig.txt )
                
                    ## add significant proteins
                    if(length(sig.idx) > 0){
                      p <- p %>% add_trace(x=dat.plot$x.ax[ sig.idx ], y=dat.plot$y.ax[ sig.idx], type='scatter', mode='markers', marker=list(size=10, color=col[sig.idx]), text=IDs[sig.idx], opacity=.2, showlegend=T, name=sig.txt  )
                    }
                    
                    ## plot interactors    
                    p <- p %>% add_trace(x=dat.plot$x.ax[ppi.int.idx], y=dat.plot$y.ax[ppi.int.idx], type='scatter', mode = 'markers', marker=list( size=10, color=col[ppi.int.idx]), text=IDs[ppi.int.idx], name='Interactors', opacity=1 )
                
                }
                ## ################################
                ## add bait
                if(length(ppi.bait.idx) > 0){
                  
                  col[ppi.bait.idx] <- 'green'
                  
                  p <- p %>% add_trace(x=dat.plot$x.ax[ppi.bait.idx], y=dat.plot$y.ax[ppi.bait.idx], type='scatter', mode = 'markers', marker=list( size=10, color=col[ppi.bait.idx]),  text=IDs[ppi.bait.idx], name=paste('Bait protein:', toupper(sub('.*_(.*)$', '\\1', ppi.bait) )), opacity=1 )
                }
                
            }  else { 
              ## ################################################################
              ## no PPI
              ## different traces for significnat and not-significant features
             non.sig.idx <- 1:nrow(dat.plot)
              if(length(sig.idx) > 0)
                non.sig.idx <- non.sig.idx[ -c(sig.idx) ]  
              
              p <- plot_ly(x=dat.plot$x.ax[non.sig.idx], y=dat.plot$y.ax[non.sig.idx], type='scatter', mode='markers', marker=list(size=10, color=col[non.sig.idx]), text=IDs[non.sig.idx], showlegend=T, name=non.sig.txt )
            
              if(length(sig.idx) > 0){
                p <- p %>% add_trace(x=dat.plot$x.ax[ sig.idx ], y=dat.plot$y.ax[ sig.idx], type='scatter', mode='markers', marker=list(size=10, color=col[sig.idx]), text=IDs[sig.idx], opacity=1, showlegend=T, name='signif.'  )
              }
              
            }
            
             
            ## ##############################
            ##  reproducibility filter?
            values.filtered <- global.results$values.filtered[[group]]  
            #View(values.filtered)
            
            if(length(values.filtered) > 0){

              ## get the entire dataset
              if(is.null(global.results$table.log))
                tab <- data.frame(global.input$table)
              else
                tab <- data.frame(global.results$table.log)
              
              ## extract filtered values
              filt.x <- tab[values.filtered, input[[ paste0('scat.x.', group) ]]]
              filt.y <- tab[values.filtered, input[[ paste0('scat.y.', group) ]]]
              
              filt.dat <- data.frame(x=filt.x, y=filt.y)
              
              ## add
              filt.txt <-  paste( global.param$filt.data, '\n(alpha=',global.param$repro.filt.val, ')',sep='')
              p <- p %>% add_trace(
                x=filt.dat$x,
                y=filt.dat$y,
                type='scatter', mode='markers', marker=list(size=10, color=rep('blue', nrow(filt.dat))), opacity=.2, text=values.filtered, name= filt.txt
              )
              
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
        plotVolcano <- function(group, interactors=NULL, 
                                verbose=T, cex.axis=1.5, cex.lab=1.5, cex.leg=1.2, cex.main=1.8,
                                sig.pch=21,
                                bg.pch=21,
                                sig.col='darkred',
                                bg.col='black'){

            if(verbose)
              cat('\n-- plotVolcano --\n')
            if(!is.null(error$msg)) return()

            ## hyperbolic curve filter?
            hyperbol <- input[[paste('ppi.hyper.filt', group, sep='.' )]]

            ## maximal log10 p-value
            max.logP <- input[[paste('max.logP', group, sep='.')]]

            ## pch for significant points
            #sig.pch=23

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
            pch.vec=rep(bg.pch, nrow(res))
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
              opac1 = 0.3
              opac2 = 0.2
            }

            #col=myColorRamp(c('grey30', 'grey50', 'darkred', 'red', 'deeppink'), na.omit(logPVal), range=c(0, max.logP), opac=opac1)
            #col.opac=myColorRamp(c('grey30', 'grey50', 'darkred', 'red', 'deeppink'), na.omit(logPVal), range=c(0, max.logP), opac=opac2)
            #col=myColorRamp(c('grey30', 'grey50', 'darkred', 'red', 'deeppink'), na.omit(logPVal), range=c(0, max.logP), opac=opac1)
            #col.opac=myColorRamp(c('grey30', 'grey50', 'darkred', 'red', 'deeppink'), na.omit(logPVal), range=c(0, max.logP), opac=opac2)
                   
            col <- rep(my.col2rgb(bg.col, alpha = 255*opac1), length(IDs))
            col[sig.idx] <- my.col2rgb(sig.col, alpha = 255*opac1)
            
            col.opac <- rep(my.col2rgb(bg.col, alpha = 255*opac2), length(IDs))
            col.opac[sig.idx] <- my.col2rgb(sig.col, alpha = 255*opac2)
            
            
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
                ppi.bait <-  toupper(sub('.*_(.*)$', '\\1', ppi.bait))
                #debug(get.interactors)
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

                leg <- ppi.map$leg
                leg.col <- ppi.map$leg.col
                ppi.col <- ppi.map$ppi.col

                ppi.bait.idx <- ppi.map$ppi.bait.idx
                if(length(ppi.bait.idx) > 0)
                  ppi.col[ ppi.bait.idx ] <- 'green'
                
                ppi.int.vec <- rep(FALSE, length(IDs))
                ppi.int.vec[ppi.int.idx] <- TRUE
  
                ## udpate color vector
                if( sum(nchar(ppi.col) > 0) > 0){
                  col[ nchar(ppi.col) > 0] <- ppi.col[ nchar(ppi.col) > 0]
                  col.opac[ nchar(ppi.col) > 0] <- ppi.col[ nchar(ppi.col) > 0 ]
                }
                #save(ppi.map, IDs,file='test.tmp.RData')
                

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
                        ppi.bait.idx <- which(toupper(sub('.*_(.*)$', '\\1', IDs) ) == toupper(ppi.bait))
                        
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
            plot.window( xlim=xlim, ylim=ylim, cex.axis=cex.axis, cex.lab=cex.lab)
            ## title
            mtext(group, side=3, cex=cex.main, line=2)
            ## label axes
            if(global.param$which.test == 'Two-sample mod T')
                mtext( paste("log(", sub('.*\\.vs\\.', '', group), "/", sub('\\.vs.*', '', group),")"), side=1, cex=cex.axis, line=3)
            else
                mtext(expression(log(FC)), side=1, cex=cex.axis, line=3)
            mtext(expression(-10*log[10](P-value)), side=2, cex=cex.axis, line=3)
            ## draw axes
            axis(1, cex.axis=cex.axis)
            axis(2, las=2, cex.axis=cex.axis)
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
                text( xlim[2]-(xlim[2]*.1), filt.minlogPVal, paste(global.param$filter.type, global.param$filter.value, sep='='), pos=3, col='grey30')
            }

            ## number of significant
            legend(ifelse(PPI, 'topright', 'top'), bty='n', legend=paste(filter.str, '\nsig / tot: ', length(sig.idx),' / ', sum(!is.na(logFC) & !is.na(logPVal)), sep=''), cex=cex.leg)

            ## ############################
            ## indicate directionality for two-sample tests
            if(global.param$which.test == 'Two-sample mod T'){
                #legend('topleft', legend=sub('\\.vs.*', '', group), cex=2, text.col='darkblue', bty='n')
                #legend('topright', legend=sub('.*\\.vs\\.', '', group), cex=2, text.col='darkblue', bty='n')
                mtext(sub('\\.vs.*', '', group), side=3, line=1, at=(xlim[1]+abs(xlim[1])*0.05), cex=cex.main, col='darkblue')
                mtext(sub('.*\\.vs\\.', '', group), side=3, line=1, at=(xlim[2]-abs(xlim[2])*0.05), cex=cex.main, col='darkblue')
                
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

                leg <- c(paste(ppi.bait), leg)
                leg.col <- c('green', leg.col)
                legend('topleft', legend=leg, col=leg.col, pch=20, bty='n', cex=cex.leg, title=ifelse( length(ppi.int.idx) > 0 ,paste('interactors sig/det/tot') , ''), pt.cex=cex.vec[1])

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
            
            if(is.null(global.results$data)) return()
  
            plotMultiScatter( 
                define.max=isolate(input$ms.max),
                min.val=isolate(input$ms.min.val), 
                max.val=isolate(input$ms.max.val), 
                robustify=isolate(input$ms.robustify),
                update=input$ms.update)
        },
        width = function(){120*(ncol(data.frame(global.input$table))-1)},
        height= function(){120*(ncol(data.frame(global.input$table))-1)}
        )

        ###############################
        ## actual plot
        plotMultiScatter <- function(define.max, min.val, max.val, robustify, update){

            cat('\n-- plotMultiScatter  --\n')
            
            ## dataset
            if(is.null(global.results$table.log))
                tab <- data.frame(global.input$table)
            else
                tab <- data.frame(global.results$table.log)

            id.col <- global.param$id.col.value
            rownames(tab) <- tab[, id.col]
            
            ## get groups
            grp <-  global.param$grp
            grp <- sort(grp)
            
            
            #############################################
            ## get correlation matrix
            if(is.null(global.results$cm) | global.param$update.cm == TRUE){
              
              ## calculate correlatio matrix
              withProgress(message = 'Correlation matrix...',{
                cm=calculateCorrMat( tab=tab,
                                     grp=grp,
                                     #id.col=id.col,
                                     lower= global.plotparam$cm.lower, upper= global.plotparam$cm.upper
                )
              })
              
              ## store results
              global.results$cm <- cm
              global.param$update.cm <- FALSE
              
            } else {
              cm <- global.results$cm
            }
            tab <- tab[, names(grp)]
            #colnames(tab) <- chopString(colnames(tab), STRLENGTH)
            
            ## table
            #tab <- tab[, setdiff(colnames(tab), id.col)]
            
            ###############################
            ## plot
            repro.filt=global.results$values.filtered
            
            withProgress({
                setProgress(message = 'Processing...', detail= 'Calculating correlations')
                my.multiscatter(tab, 
                                cm=cm,
                                repro.filt=global.results$values.filtered, 
                                grp=grp,  
                                grp.col.legend=global.param$grp.colors.legend[unique(grp)], 
                                define.max=define.max, 
                                max.val=max.val, 
                                min.val=min.val, 
                                robustify=robustify,
                                update=update)
            })
        }
 
        
        ## ##############################################
        ## observe input to trigger update 
        observeEvent(c(input$cm.lower, input$cm.upper),{
          
          global.plotparam$cm.lower <- input$cm.lower
          global.plotparam$cm.upper <- input$cm.upper
          
          global.param$update.cm <- TRUE
        })
        
        #####################################################
        ## correlation matrix
        output$correlation.matrix <- renderPlot({
            if(is.null(global.results$data)) return()
             withProgress({
                 setProgress(message = 'Processing...', detail= 'Generating Heatmap')
                  
                  if(is.null(global.results$table.log))
                    tab <- data.frame(global.input$table)
                  else
                    tab <- data.frame(global.results$table.log)
               
                  ## dataset
                  #tab <- data.frame(global.results$table.log)
                  
                  ## id column
                  id.col.value <- global.param$id.col.value
                  
                  #tab[, c(id.col.value, names(grp))]
                  ## group vector
                  grp <- global.param$grp
                  ## group colors
                  grp.col <- global.param$grp.colors
                  grp.col.leg <- global.param$grp.colors.legend
               
              
                  ## update selection
                  tab <- tab[, c(id.col.value, names(grp))]
                  grp.col <- grp.col[names(grp)]
                  grp.col.leg <- grp.col.leg[unique(grp)]
                
                  
                  ## get correlation matrix
                  if(is.null(global.results$cm) | global.param$update.cm == TRUE){
                    
                    
                    ## calculate correlatio matrix
                    withProgress(message = 'Correlation matrix...',{
                      cm=calculateCorrMat( tab=tab,
                                           grp=grp,
                                           #id.col=id.col.value,
                                           lower= global.plotparam$cm.lower, upper= global.plotparam$cm.upper
                                           )
                    })
                    
                    ## store results
                    global.results$cm <- cm
                    global.param$update.cm <- FALSE
                    
                  } else {
                    cm <- global.results$cm
                  }
                  
                 plotCorrMat(#tab=tab,
                             #id.col=id.col.value,
                            cm=cm,
                             grp=grp,
                             grp.col.legend=grp.col.leg,
                             #lower=input$cm.lower, upper=input$cm.upper, display_numbers=input$cm.numb)
                             lower= global.plotparam$cm.lower, upper= global.plotparam$cm.upper, display_numbers=input$cm.numb, width=dynamicWidthHM( length(global.param$grp), unit='in'), height=dynamicWidthHM( length(global.param$grp), unit='in' ))
             })
        })#, width=1200, height=1000)

        #####################################################
        ## correlation matrix transposed
       ## output$correlation.matrix.trans <- renderPlot({
       ##     if(is.null(global.results$data)) return()
       ##      withProgress({
       ##          setProgress(message = 'Processing...', detail= 'Generating Heatmap')
       ##          plotCorrMat(lower=input$cm.lower, upper=input$cm.upper, trans=T, display_numbers=input$cm.numb)
       ##      })
       ## }, width=1200, height=1000)




        


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

            
            ## update selection
            tab <- tab[, c(id.col, names(grp))]
            #grp.col <- grp.col[names(grp)]
            grp.col.legend <- grp.col.legend[unique(grp)]
            
            ## table
            #tab <- tab[, setdiff(colnames(tab), id.col)]
            tab <- tab[, names(grp)]

            ## ############################################
            ## get correlation matrix
            if(is.null(global.results$cm) | global.param$update.cm == TRUE){
              
              
              ## calculate correlatio matrix
              withProgress(message = 'Correlation matrix...',{
                cm=calculateCorrMat( tab=tab,
                                     grp=grp,
                                    #id.col=id.col,
                                     lower= global.plotparam$cm.lower, upper= global.plotparam$cm.upper
                )
              })
              
              ## store results
              global.results$cm <- cm
              global.param$update.cm <- FALSE
              
            } else {
              cm <- global.results$cm
            }
           
            ## #############################################
            ## extract correlations for each group
            cor.group <- lapply(names(grp.col.legend), function(x){
              cm.grp=cm[ names(grp) [grp == x], names(grp) [grp == x]]
              cm.grp[upper.tri(cm.grp, diag = FALSE)]
              }  )
            #save(cor.group, file='cor.group.RData')
            
            #ylim <- c(min(unlist(cor.group))-0.1*min(unlist(cor.group)), max( unlist(cor.group))+0.1*max(unlist(cor.group)) )
            ylim <- c(min(unlist(cor.group))-0.1*min(unlist(cor.group)), 1)
            
            # plot
            #debug(fancyBoxplot)
            par(mar=c(8, 5, 2, 1))
            fancyBoxplot(cor.group, col='white', box.border=unlist(grp.col.legend), vio.alpha = 0, lwd=2.5,
                         names=names(grp.col.legend), grid=F,
                         show.numb = 'median', numb.col='black', 
                         main=paste('Pairwise intra-group correlations (', global.plotparam$cm.upper,')'),
                         ylab='Correlation coeffient',
                         ylim=ylim,
                         numb.cex=1,
                         las=2)
            
            legend('topright', legend=names(grp.col.legend), fill=grp.col.legend, bty='n', ncol = ifelse(length(grp.col.legend) > 4, 2, 1))
            
        })

        ## #################################################################################
        ## 
        ##                               Morpheus
        ##
        ## #################################################################################
        # output$HM.morpheus <- renderMorpheus({
        # 
        #   if(!is.null(error$msg)) return()
        # 
        #   ## extract results
        #   res = global.results$filtered
        # 
        # 
        #   ## require at least three significant hits
        #   validate(need(nrow(res) > 1, 'Need at least 2 features to draw a heatmap!'))
        # 
        # 
        #   ## extract expression values
        #   res <- res[, names(global.param$grp)]
        #   res <- data.matrix(res)
        # 
        #   ## morpheus
        #   withProgress({
        #     setProgress(message = 'Creating morpheus widget...', detail= 'hold on')
        #     morpheus(res, Rowv = F, Colv = F, labCol = chopString(colnames(res), nChar = 20))
        #   })
        # 
        # } )

        ####################################################################################
        ##
        ##                                  Heatmap
        ##
        ####################################################################################
        output$HM <- renderPlot({

            if(!is.null(error$msg)) return()

            ######################################
            ## extract results
            res = global.results$filtered

            ######################################
            ## require at least three significant hits
            validate(need(nrow(res) > 1, 'Need at least 2 features to draw a heatmap!'))

            ## number of features with at least one non-missing value
            #n.feat <- sum( apply() )
            
            #######################################
            ## heatmap title
            hm.title <- paste('filter:', global.param$filter.type, ' / cutoff:', global.param$filter.value, sep='')
            #hm.title = paste(hm.title, '\nsig / total: ', nrow(res), ' / ', nrow( global.results$data$output ), sep='')
            hm.title = paste(hm.title, '\nsig / total: ', global.results$N.feat.filtered, ' / ',  global.results$N.feat, sep='')
            
            ###################################
            ## ids to show in heatmap
            hm.rownames <- res[, 'id.concat']
            
            #######################################
            ## extract expression values
            res = res[, names(global.param$grp)]

            ##@#####################################
            ##  dimensions depending on no. rows/columns
            cw <- cwHM(ncol(res))
         
            #if(!is.null(global.input$cdesc)){
            if(!is.null(global.param$anno.col)){
              anno.col=global.param$anno.col
              anno.col.color=global.param$anno.col.color
            } else {
              anno.col=data.frame(Group=global.param$grp)
              anno.col.color=list(Group=global.param$grp.colors.legend)
            }

            ######################################
            ## plot
            if(input$hm.max){
                withProgress({
                    setProgress(message = 'Processing...', detail= 'Generating Heatmap')
                    #plotHM(res=res, grp=global.param$grp, grp.col=global.param$grp.colors, grp.col.legend=global.param$grp.colors.legend,  hm.clust=input$hm.clust, hm.title=hm.title, hm.scale=input$hm.scale, cellwidth=cw, fontsize_row=input$cexRow, fontsize_col=input$cexCol, max.val=input$hm.max.val, style=global.param$which.test, cdesc=hm.cdesc, cdesc.grp=global.param$grp.gct3)
                    plotHM(res=res, hm.rownames=hm.rownames, grp=global.param$grp, grp.col=global.param$grp.colors, grp.col.legend=global.param$grp.colors.legend,  hm.clust=input$hm.clust, hm.title=hm.title, hm.scale=input$hm.scale, cellwidth=cw, fontsize_row=input$cexRow, fontsize_col=input$cexCol, max.val=input$hm.max.val, style=global.param$which.test, anno.col=anno.col, anno.col.color=anno.col.color, show.rownames=input$hm.show.rownames, show.colnames=input$hm.show.colnames)
                  
                  })
            } else {
                 withProgress({
                   setProgress(message = 'Processing...', detail= 'Generating Heatmap')
                   
                   plotHM(res=res, hm.rownames=hm.rownames, grp=global.param$grp, grp.col=global.param$grp.colors, grp.col.legend=global.param$grp.colors.legend,  hm.clust=input$hm.clust, hm.title=hm.title, hm.scale=input$hm.scale, cellwidth=cw, fontsize_row=input$cexRow, fontsize_col=input$cexCol, style=global.param$which.test, anno.col=anno.col, anno.col.color=anno.col.color, show.rownames=input$hm.show.rownames, show.colnames=input$hm.show.colnames)
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
          
          ###################################
          ## ids to show in heatmap
          hm.rownames <- res[, 'id.concat']
          
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
          #debug(plotHM)
          ######################################
          ## plot
          if(input$hm.int.max){
            withProgress({
              setProgress(message = 'Processing...', detail= 'Generating Heatmap')
              #plotHM(res=res, grp=global.param$grp, grp.col=global.param$grp.colors, grp.col.legend=global.param$grp.colors.legend,  hm.clust=input$hm.int.clust, hm.title=hm.title, hm.scale=input$hm.int.scale, cellwidth=cw, max.val=input$hm.int.max.val, style=global.param$which.test, cdesc=hm.cdesc, cdesc.grp=global.param$grp.gct3, plotly = T)
              plotHM(res=res, hm.rownames=hm.rownames, grp=global.param$grp, grp.col=global.param$grp.colors, grp.col.legend=global.param$grp.colors.legend,  hm.clust=input$hm.int.clust, hm.title=hm.title, hm.scale=input$hm.int.scale, cellwidth=cw, max.val=input$hm.int.max.val, style=global.param$which.test, anno.col=anno.col, anno.col.color=anno.col.color, plotly = T)
            })
          } else {
            withProgress({
              setProgress(message = 'Processing...', detail= 'Generating Heatmap')
              #plotHM(res=res, grp=global.param$grp, grp.col=global.param$grp.colors, grp.col.legend=global.param$grp.colors.legend,  hm.clust=input$hm.int.clust, hm.title=hm.title, hm.scale=input$hm.int.scale, cellwidth=cw, style=global.param$which.test, cdesc=hm.cdesc, cdesc.grp=global.param$grp.gct3, plotly = T)
              plotHM(res=res, hm.rownames=hm.rownames, grp=global.param$grp, grp.col=global.param$grp.colors, grp.col.legend=global.param$grp.colors.legend,  hm.clust=input$hm.int.clust, hm.title=hm.title, hm.scale=input$hm.int.scale, cellwidth=cw, style=global.param$which.test, anno.col=anno.col, anno.col.color=anno.col.color, plotly = T)
              
            })
          }
        }
        )
        ## ################################################################
        ##   
        ##                         fanplot 
        ##
        ## ################################################################
        observeEvent( input$HC.fan.show.tip.label, {
          if(!input$HC.fan.show.tip.label)
            updateNumericInput(session = session, inputId = 'HC.fan.tip.cex', value = 4)
          else
            updateNumericInput(session = session, inputId = 'HC.fan.tip.cex', value = 1)
          })
        ## plot
        output$HC.fan <- renderPlot({
          
          if(!is.null(error$msg)) return()
          
          ######################################
          ## extract results
          res = global.results$filtered
          
          ######################################
          ## require at least three significant hits
          validate(need(nrow(res) > 1, 'Need at least 2 features to draw a fanplot!'))
        
          #######################################
          ## extract expression values
          res = res[, names(global.param$grp)]
          
          if(!is.null(global.input$cdesc))
            hm.cdesc <- global.input$cdesc
          else
            hm.cdesc <- NULL
            #hm.cdesc <- global.param$grp
          
          # groups and colors
          grp <- global.param$grp
          grp.col <- global.param$grp.colors
          grp.col.leg <- global.param$grp.colors.legend
          
          # update selection
          grp.col <- grp.col[names(grp)]
          grp.col.leg <- grp.col.leg[unique(grp)]
          
          # check
          HC.fan.show.tip.label <- input$HC.fan.show.tip.label
          HC.fan.tip.cex <- input$HC.fan.tip.cex
          validate(need(!is.na(HC.fan.tip.cex), 'Label size must be specified.'))
          #if(!HC.fan.show.tip.label){
            #global.plotparam$HC.fan.tip.cex <- 4
          #  updateNumericInput(session = session, inputId = 'HC.fan.tip.cex', value = 4)
          #} else {
            #global.plotparam$HC.fan.tip.cex <- 1
           # updateNumericInput(session = session, inputId = 'HC.fan.tip.cex', value = 1)
          #}
          
          ######################################
          ## plot
            withProgress({
              setProgress(message = 'Processing...', detail= 'Fanplot')
              plotFAN(res=res, grp=grp, grp.col=grp.col, grp.col.legend=grp.col.leg, 
                      show.tip.label=HC.fan.show.tip.label, 
                      tip.cex=HC.fan.tip.cex)
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

            #res <- global.results$filtered

            #validate(need(nrow(res) > 2, 'Need at least 3 features to perform PC.'))

            grp <- global.param$grp

            
            if(is.null(global.results$pca) | global.param$update.pca == TRUE){
              
              res <- global.results$filtered
              
              validate(need(nrow(res) > 2, 'Need at least 3 features to perform PC.'))
              
              ## run PCA
              withProgress(message = 'PCA...',{
                pca=my.prcomp2( res, grp )
              })
              ## store results
              global.results$pca <- pca
              global.param$update.pca <- FALSE
              
            } else {
              pca <- global.results$pca
            }
            
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
          #if(is.null(global.results$data) | is.na(global.param$filter.value) | is.null(global.results$pca)) return()
          
          if(is.null(global.results$data) | is.na(global.param$filter.value)) return()
          ##if(!is.null(error$msg)) return()

          grp <- global.param$grp
          grp.unique <- unique(grp)
          grp.colors <- global.param$grp.colors[names(grp)]
          
          if(is.null(global.results$pca) | global.param$update.pca == TRUE){
 
            res <- global.results$filtered
          
            validate(need(nrow(res) > 2, 'Need at least 3 features to perform PC.'))
            
            ## run PCA
            withProgress(message = 'PCA...',{
              pca=my.prcomp2( res, grp )
            })
            ## store results
            global.results$pca <- pca
            global.param$update.pca <- FALSE
            
          } else {
            pca <- global.results$pca
          }
          
          
          
         # pca <- global.results$pca
          
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
            
          if(is.null(global.results$data) | is.na(global.param$filter.value)) return()
            
            grp <- global.param$grp
            grp.unique <- unique(grp)
            grp.colors <- global.param$grp.colors[names(grp)]
            
            
            
            if(is.null(global.results$pca) | global.param$update.pca == TRUE){
              
              res <- global.results$filtered
              
              validate(need(nrow(res) > 2, 'Need at least 3 features to perform PC.'))
              
              ## run PCA
              withProgress(message = 'PCA...',{
                pca=my.prcomp2( res, grp )
              })
              ## store results
              global.results$pca <- pca
              global.param$update.pca <- FALSE
              
            } else {
              pca <- global.results$pca
            }
            
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
