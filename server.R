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
## Last updated January 18, 2023 by Natalie Clark (nclark@broadinstitute.org) - v 1.1.2
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
            gprof=NULL,
            
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
            
            file.gct3=F,                 ## GCT v1.3?
            
            session=NULL,                ## session id
            user=NULL,                   ## user
            grp=NULL,                    ## group to use for statistical testing
            N.grp=NULL,                  ## number of defined groups for statistical testing
            grp.norm = NULL,             ## group to use for group-wise normalization
            grp.colors=NULL,             ## group color assignment
            grp.colors.legend=NULL,      ## group colors, names are group names
            grp.colors.norm=NULL,        ## group color assignment for normalization
            grp.colors.legend.norm=NULL, ## group colors for normalization, names are group names
            
            tabsep = '\t',               ## default separator for non-gct files
            tabsep.anno = '\t',          ## default separator for non-gct files, used for annotation file
            
            which.test='One-sample mod T', ## specify test
            
            log.transform='none',        ## log transformation
            norm.data='none',            ## data normalization
            norm.per.group=FALSE,        ## normalize per group? 
            filt.data='none',            ## data filtering

            ##repro.filt='no',           ## reproducibility filter
            repro.filt.val=0.001,
            ##sd.filt='no',              ## sd filter
            sd.filt.val=10,              ## remove lower 10 percent of features with lowest sd
            
            na.filt.val=100,             ## max. % missing  values

            
            filter.type='adj.p',         ## default filter
            filter.value=0.05,           ## default filter value

            run.test=0,                   ## number of times the 'Run analysis' button has been pressed during a user session

            update.ppi.select= FALSE,     ## trigger selectize, volcano
            update.ppi.select.scat=FALSE, ## trigger selectize, scatterplot
            collapse.ppi=TRUE,             ## should PPI query panel be collapsed?
            
            update.pca=FALSE,
            update.cm=FALSE,
            update.gprof.up=TRUE,
            update.gprof.dn=TRUE

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
            volc.label="ID_Symbol", #default text label
            
            volc.init=T,
            
            ## heatmap
            hm.cexCol=8,
            hm.cexRow=3,
            hm.scale="none",
            hm.max=FALSE,
            hm.max.val=4,
            hm.show.rownames=T,
            hm.show.colnames=T,
            hm.clust="none",
            
            ## PCA
            pca.x='PC 1',
            pca.y='PC 2',
            pca.z='PC 3',
         
            ## correlation matrix
            cm.upper='pearson',
            cm.lower='spearman',
            cm.numb=FALSE,
            
            ## gprofiler
            gprof.source.all=c('GO:BP', 'GO:MF', 'GO:CC', 'KEGG', 'WP'),
            gprof.source.selected=c('GP:BP', 'KEGG'),
            gprof.fdr=0.01,
          
            ## fanplot
            HC.fan.show.tip.label=T,
            HC.fan.tip.cex=1,
            
            ## profile plot
            profile.plot.xlim='as-is'  ## 'symmetric' or 'as-is'
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
            HTML('<p align=\"center\"><font size=\"5\" color=\"red\">To analyze another data set or to start over, refresh this page.</font></p>' )
        })

        #####################################
        ## Error messages
        output$error <- renderText({
         
            if( is.null(error$msg) ) return()
          
            if(is.null(error$title))
              etitle <- 'Error'
            else
              etitle <- error$title
            
            shinyalert(etitle, error$msg, type = "error", html=TRUE)
        })

        ## #########################################################
        ## toggle all options for export
        ## #########################################################
        observeEvent(input$export.toggle.all, {

            updateCheckboxInput(session, "export.hm", "Heatmap", value=!input$export.toggle.all)
            updateCheckboxInput(session, 'export.box', 'Boxplots', value=!input$export.toggle.all)
            updateCheckboxInput(session, 'export.volc', 'Volcano plot', value=!input$export.toggle.all)
            updateCheckboxInput(session, 'export.phist', 'p-value histogram', value=!input$export.toggle.all)
            updateCheckboxInput(session, 'export.pca', 'PCA', value=!input$export.toggle.all)
            #updateCheckboxInput(session, 'export.pca.loadings', 'PCA loadings (xls)', value=!input$export.toggle.all)
            updateCheckboxInput(session, 'export.ms', 'Multiscatter', value=!input$export.toggle.all)
            updateCheckboxInput(session, 'export.excel', 'Excel sheet', value=!input$export.toggle.all)
            updateCheckboxInput(session, 'export.gct.file', 'GCT files: 1) normalized data and 2) signed log p-values', value=!input$export.toggle.all)
            updateCheckboxInput(session, 'export.cm', 'Correlation matrix', value=!input$export.toggle.all)
            updateCheckboxInput(session, 'export.cb', 'Correlation boxplot', value=!input$export.toggle.all)
            updateCheckboxInput(session, 'export.profile', 'Profile plot', value=!input$export.toggle.all)

        })


        ## ########################################################
        ## session name
        output$session.label <- renderMenu({
          
          if(!global.param$session.saved) return()
          label <- paste("Session:", global.param$label)
          
          notificationItem(label, status='success', shiny::icon("folder", "fa-1x", lib='glyphicon')) 
        })
        
        
        ## ########################################################
        ## logged user
        output$logged.user <- renderMenu({

            if(is.null(session$user)) return()
            user <- session$user
          
            notificationItem(user, shiny::icon("user", "fa-1x", lib='glyphicon'), status='success') 
        })

        ## #########################################################
        ## logout user
        output$logout <- renderMenu({
            if(is.null(session$user)) return()
            notificationItem('Logout', icon=shiny::icon("sign-out", "fa-1x"), status='success', href="__logout__")
        })


        ################################################################################
        ##
        ##                      navbar - render UI
        ##
        ################################################################################
        output$navbar <- renderUI({

            if(!global.param$analysis.run) return()

            ##############################################
            ## determine the number of group comparisons,
            ## e.g. for the number of volcano plots to draw
            ##############################################
            groups.comp <- unique(global.param$grp.comp)

            ## ##########################################
            ## class vector
            grp <- global.param$grp
            grp.norm <- global.param$grp.norm
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
                                                             column(3, HTML('Export all checked files'))),
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
                                                  checkboxInput('export.phist', 'p-value histogram',value=T),
                                                  checkboxInput('export.pca', 'PCA',value=T),
                                                  #checkboxInput('export.pca.loadings', "PCA loadings (xls)", value = T),
                                                  checkboxInput('export.ms', 'Multiscatter',value=T),
                                                  checkboxInput('export.excel', 'Excel sheet',value=T),
                                                  checkboxInput('export.gct.file', 'GCT files: 1) normalized data and 2) signed log p-values',value=T),
                                                  checkboxInput('export.cm', 'Correlation matrix',value=T),
                                                  checkboxInput('export.cb', 'Correlation boxplot',value=T),
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
                                        box(title="Dataset:", solidHeader = T, status = "primary", width = 4,# color='purple',
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

            ############################################
            ## gProfiler
            gprof.tabs <- list()
            gprof.tabs[[1]] <- 'Functional enrichment analysis'
            for(i in 1:length( grp.unique )){
              
              gprof.tabs[[i+1]] <- tabPanel(paste0( grp.unique[ i ] ),
                                            fluidPage(
                                              ## ###############################
                                              ## plot
                                              box( title=grp.unique[i], status = 'primary', solidHeader = T, width=12,
                                                   fluidRow(column(2, selectInput(paste0('gprof.src.', grp.unique[i]), label='Source' , 
                                                                                  choices=unlist(global.plotparam$gprof.source.all), 
                                                                                  multiple = T, 
                                                                                  selected = unlist(global.plotparam$gprof.source.selected )) ),
                                                            column(2, numericInput(paste0('gprof.fdr.', grp.unique[i]), label='FDR', value = global.plotparam$gprof.fdr, min = 0, max=1, step = 0.01)),
                                                            column(8)
                                                            ),
                                                   fluidRow(
                                                     column(12, align='center', plotlyOutput(  paste0('gprof.up.', grp.unique[i]) , width=800, height=600) )
                                                   ), ## end fluiRow
                                                   fluidRow(
                                                     column(12, align='center', plotlyOutput(  paste0('gprof.dn.', grp.unique[i]) , width=800, height=600) )
                                                   )
                                                   
                                              ) ## end box
                                            )
                                            
                                            )
            }
            
            ## ##########################################
            ## SCATTERPLOTS plotly
            ## - for each group comparison
            ##
            ## ##########################################
            scat.tabs <- list()
            scat.tabs[[1]] <- 'Scatterplots'

            for(i in 1:length( grp.unique )){

                ## extract data columns of current experiment
                grp.idx <- which(grp == grp.unique[i])
                
                ## index of other experiments
                other.idx <- setdiff(1:length(grp), grp.idx)
                grp.idx <- c(grp.idx, other.idx)
                scat.x <- names(grp)[grp.idx]
                scat.y <- names(grp)[grp.idx]

                scat.tabs[[i+1]]=tabPanel(paste0( grp.unique[ i ] ),
                                          fluidPage(

                                              ## ###############################
                                              ## select data columns to plot
                                              box( title='Choose expression columns', status = 'primary', solidHeader = T, width=12,

                                                  fluidRow(
                                                      column(1),
                                                      #column(1, checkboxInput(paste0('scat.showall'), label='Show all columns', value=FALSE)),
                                                      column(5, selectInput( paste0('scat.x.', grp.unique[i]), 'x-axis', choices=scat.x, selected=scat.x[1], multiple=FALSE, selectize=TRUE)),
                                                      column(5, selectInput( paste0('scat.y.', grp.unique[i]), 'y-axis', choices=scat.y, selected=scat.y[2], multiple=FALSE, selectize=TRUE)),
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
    
    
                                                              column(2, numericInput( paste("cex.volcano",groups.comp[i], sep='.'), "Point size", value=global.plotparam$volc.ps, min=1, step=1, width='100px')),
                                                              ##column(1, numericInput( paste("opac.volcano",groups.comp[i],sep='.'), "Opacity %", value=50, min=0, max=100, step=10)),
                                                              column(2, numericInput( paste("cex.volcano.lab",groups.comp[i],sep='.'), "Label size", value=global.plotparam$volc.ls, min=.1, step=.1, width='100px')),
                                                              column(2, selectInput( paste("grid.volcano",groups.comp[i],sep='.'), "Grid", c(T, F), selected=global.plotparam$volc.grid, width='100px')),
                                                              column(2, numericInput( paste( "max.logP", groups.comp[i], sep='.'), "Max. Log10(p-value)", value=global.plotparam$volc.maxp, min=20, max=300, step=10, width='100px') ),
                                                              column(2, selectInput( paste("volc.label",groups.comp[i],sep='.'), "Labels", c("ID_Symbol", "ID","Symbol"), selected=global.plotparam$volc.label, width='100px'))
    
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
                                                      column(width=7,
                                                             box( width=NULL,  title='Volcano plot', status = 'primary', solidHeader = T,
                                                                 plotOutput( paste("volcano", groups.comp[i], sep='.'), height=600, click=paste('plot_click', groups.comp[i], sep='.'), hover=hoverOpts(id=paste('plot_hover', groups.comp[i], sep='.'), delay=10), brush=brushOpts(id=paste('plot_brush', groups.comp[i], sep='.'), resetOnNew=T, delayType='debounce', delay='1000' ), dblclick=paste('plot_dblclick', groups.comp[i], sep='.'))
                                                                 )),
                                                      ## table
                                                      column(width=5,
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
            # hm.int.tab <-  tabPanel('Interactive heatmap',
            #                     box(title='Interactive Heatmap', status = 'primary', solidHeader = T, width="100%", height="100%",
            #                         fluidRow(
            #                           column(2, selectInput( "hm.int.scale", "Scale", c("row","column","none"), selected=global.plotparam$hm.scale)),
            #                           column(2, selectInput( "hm.int.clust", "Cluster", c("column","row","both","none"), selected=ifelse(global.param$which.test != "mod F", "none" ,"both"))),
            #                           column(2, checkboxInput('hm.int.max', 'Cap values', value=global.plotparam$hm.max)),
            #                           column(2, numericInput( "hm.int.max.val", "Max. value", value=global.plotparam$hm.max.val, step=1, min=2))
            #                         ),
            #                         fluidRow(
            #                           column(12, align='center', plotlyOutput("HM.int", height=1000, width=1000 ))
            #                         )
            #                     )
            # )
            # 
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
            
            clust.tab <- vector('list', 2)
            clust.tab[[1]] <- hm.tab
            #clust.tab[[2]] <- hm.int.tab
            #clust.tab[[3]] <- hc.fanplot
            clust.tab[[2]] <- hc.fanplot
            
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
            if(global.param$file.gct3){    ## CGT v1.3 - the user can user different annotations to color the plots
                
                    pca.tab2[[2]] <- tabPanel('PC plots',
                                              fluidPage(
        
                                                  fluidRow(
                                        
                                                        box( title = 'Select principle components', status='primary', solidHeader = T, align='center',
                                                        column(width=4, selectInput('pca.x', 'x-axis', paste('PC', 1:10)), selected=global.plotparam$pca.x),
                                                        column(width=4, selectInput('pca.y', 'y-axis', paste('PC', 1:10), selected=global.plotparam$pca.y)),
                                                        column(width=4, selectInput('pca.z', 'z-axis', paste('PC', 1:10), selected=global.plotparam$pca.z))
                                                    
                                                    ),
                                                         box(title='Select annotation', status='primary', solidHeader = T, align='center',
                                                             column(6, selectInput('pca.grp.col', label = 'Color by', choices = colnames(global.input$cdesc), selected = global.plotparam$pca.grp.col),
                                                             column(6),         
                                                               )
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
            } else { ## no GCT v1.3
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
            }
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
                         fluidRow(column(12, dataTableOutput("tableprev")))
                     )
             )
            #############################################
            ## QC tabs
            ##
            #############################################
            qc.tabs <- vector('list', 5)
            ##names(qc.tabs) <- c('Boxplots', 'p-values', 'Multi scatter', 'Correlation matrix', 'Correlation matrix transposed')

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
                                                                   box(title='Data range', solidHeader=T, status='primary',
                                                                       column(width=6, selectInput('profile.plot.xlim', label='x-axis', choices=c('symmetric', 'as-is'), selected = global.plotparam$profile.plot.xlim))
                                                                       )
                                                               ),
                                                               fluidRow(
                                                                   box(title="Before normalization", solidHeader=T, status="primary",
                                                                       column(width=12, plotlyOutput("expr.profile"))
                                                               ))
                                                            )
                                                       } else {
                                                           fluidPage(
                                                               fluidRow(
                                                                   box(title='Data range', solidHeader=T, status='primary',
                                                                       column(width=6, selectInput('profile.plot.xlim', label='x-axis', choices=c('symmetric', 'as-is'), selected = global.plotparam$profile.plot.xlim))
                                                                   )
                                                               ),
                                                               fluidRow(
                                                                   box(title="Before normalization", solidHeader=T, status="primary",
                                                                      column(width=12, plotlyOutput("expr.profile"))
                                                                      ),
                                                                   box(title="After normalization", solidHeader=T, status="primary",
                                                                       column(width=12, plotlyOutput("expr.profile.norm"))
                                                                       )
                                                               )
                                                           )
                                                       }
                                                   )

            ###########################
            ## P-value distribution
            ##qc.tabs[['P-values']] <- tabPanel('P-values',
            qc.tabs[[3]] <- tabPanel('p-values',
                                              fluidPage(
                                                  fluidRow(
                                                      box(title="Distribution of p-values", solidHeader=T, status="primary", width=1000, height=600*ifelse( global.param$which.test != 'mod F', length(unique(global.param$grp.comp)), 1 ),
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
                       # navbarMenu('Clustering', clust.tab[[1]], clust.tab[[2]], clust.tab[[3]]),
                       navbarMenu('Clustering', clust.tab[[1]], clust.tab[[2]]),
                       
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

                       ## ####################################
                       ##          insert gprofiler plots
                       ## ####################################
                       #do.call(navbarMenu, gprof.tabs),
                      
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
                      #navbarMenu('Clustering', clust.tab[[1]], clust.tab[[2]], clust.tab[[3]]),
                      navbarMenu('Clustering', clust.tab[[1]], clust.tab[[2]]),
                      
                      #do.call(navbarMenu, hm.tab),
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
            ## only works in SSP 
            if(!is.null(session$user) & file.exists(USERDB)){
                
                cat('\n---------------\n', 
                    paste("Found session db at \"", USERDB),
                    paste("Parsing session db for user: ", session$user, collapse='\n'), 
                    '\n-----------------')
                
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
                user.roles <- read.delim( USERDB, stringsAsFactors=F)
                
                ## check whether the user appears as collaborator in a project
                ## 20210830: - session$user was expected to be an email address (authenification on SSP via Google Auth) 
                ##           - session$user on RSC is just the user name (up to the @ sign)
                idx <- grep( paste('(^|;)', session$user, '(@|$|;)', sep=''), user.roles$collaborator)

                ## if so add the project path to the search path
                if(length(idx) > 0){
                    for(i in 1:length(idx)){

                        ## folder of the project OWNER to be parsed
                        dir.owner <- paste(DATADIR, sub('@.*', '', user.roles$owner[idx[i]]), sep='')

                        ## check if the folder exists (if not, 'user-roles.txt' has not been updated)
                        if(dir.exists(dir.owner)){
                            tmp <- grep( paste(user.roles$project[ idx[i] ], '_session.*RData$', sep='' ),
                                        dir( dir.owner, full.names=T, recursive=T), value=T)
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

            ## don't show the panel if there is no saved session
            if(length(saved.sessions) == 0) return()

            ## get the time stamp of the files
            #time.tmp <- file.info(saved.sessions)$ctime
            time.tmp <- file.info(saved.sessions)$mtime
            
            names(saved.sessions) <-  paste( sub('_.*','', sub('.*/','',saved.sessions)), time.tmp, sep='_' )

            ## order by time
            saved.sessions <- saved.sessions[ order(time.tmp) ]

            ## store saved sessions
            global.param$saved.sessions <- saved.sessions

            list(
                selectizeInput(inputId = 'session.browse', 
                               label = 'Saved sessions:',
                               choices=names(saved.sessions), 
                               
                               options=list( 
                                 placeholder = 'Search sessions',
                                 onInitialize = I('function() { this.setValue(""); }')
                               )
                               ),
              
                actionButton('session.browse.import', 'Import'),
                actionButton('session.manage', 'Manage sessions', onclick =paste("window.open('", CONFAPP,"', 'newwindow', 'width=500 height=600'); return false;", sep=''))
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
        #  UI:   group assignment for gct 1.3 files
        # - id column defined as first column in gct
        # - group assignment from column annotations
        output$define.groups.gct3 <- renderUI({
          
          if(!global.param$file.gct3) return()
          if(global.param$grp.done) return()
          
          ## list all column description columns (cdesc)    
          if( is.null(input$grp.norm.check)){
          list(
            checkboxInput("QC.filter","Filter out QC.fail samples (requires QC.status column)",FALSE),
            radioButtons("grp.gct3", "Choose column for statistical testing", colnames(global.input$cdesc)),
            checkboxInput("grp.norm.check","Perform group-wise normalization",FALSE),
            actionButton("update.grp.gct3", 'OK')
          )
          }else if(input$grp.norm.check){
            list(
              checkboxInput("QC.filter","Filter out QC.fail samples (requires QC.status column)",FALSE),
              radioButtons("grp.gct3", "Choose column for statistical testing", colnames(global.input$cdesc)),
              checkboxInput("grp.norm.check","Perform group-wise normalization",TRUE),
              radioButtons("grp.norm", 'Choose column for group-wise normalization.', colnames(global.input$cdesc)),
              actionButton("update.grp.gct3", 'OK')
            )
          }else if(!input$grp.norm.check){
            list(
              checkboxInput("QC.filter","Filter out QC.fail samples (requires QC.status column)",FALSE),
              radioButtons("grp.gct3", "Choose column for statistical testing", colnames(global.input$cdesc)),
              checkboxInput("grp.norm.check","Perform group-wise normalization",FALSE),
              actionButton("update.grp.gct3", 'OK')
            )
          }
          
        })
        # preview current selection
        output$grp.gct3.prev <- renderText({
          if(is.null(input$grp.gct3)) return()
          if(global.param$grp.done) return()
          
          if (is.null(input$grp.norm.check)){
            HTML(paste('<br><p><font size=\"4\"><b>Current selection for statistical testing:</b><b>', input$grp.gct3,'</b><br>','<br><font size=\"4\"><b>Perform group-wise normalization:</b><b>', input$grp.norm.check,'</b><br>','<br><font size=\"4\"><b>Filter QC.fail samples:</b><b>', input$QC.filter,'</b></p><br>'))
            }else if(input$grp.norm.check){
            HTML(paste('<br><p><font size=\"4\"><b>Current selection for statistical testing:</b><b>', input$grp.gct3,'</b><br>','<br><font size=\"4\"><b>Perform group-wise normalization:</b><b>', input$grp.norm.check,'</b><br>','<br><font size=\"4\"><b>Current selection for group-wise normalization:</b><b>', input$grp.norm,'</b><br>','<br><font size=\"4\"><b>Filter QC.fail samples:</b><b>', input$QC.filter,'</b></p><br>'))
          }else if(!input$grp.norm.check){
            HTML(paste('<br><p><font size=\"4\"><b>Current selection for statistical testing:</b><b>', input$grp.gct3,'</b><br>','<br><font size=\"4\"><b>Perform group-wise normalization:</b><b>', input$grp.norm.check,'</b><br>','<br><font size=\"4\"><b>Filter QC.fail samples:</b><b>', input$QC.filter,'</b></p><br>'))
          }
        })
        # preview levels current selection for statistical testing
        output$grp.gct3.prev.tab <- renderTable({
          if(is.null(input$grp.gct3)) return()
          if(global.param$grp.done) return()
          
          #if QC.filter is checked but QC.status column does not exist inform the user
          if(input$QC.filter & "QC.status"%in%colnames(global.input$cdesc)){
            tab <- table(global.input$cdesc[global.input$cdesc$QC.status!="QC.fail", input$grp.gct3])
            global.param$QC.filter=input$QC.filter
          }else if (input$QC.filter & !"QC.status"%in%colnames(global.input$cdesc)){
            shinyalert("No QC.status column detected!","No QC.status column detected for filtering. Analysis will proceed using all samples.", type="warning")
            tab <- table(global.input$cdesc[, input$grp.gct3])
            global.param$QC.filter=FALSE
          }else{
            tab <- table(global.input$cdesc[, input$grp.gct3])
            global.param$QC.filter=input$QC.filter
          }
          
          tab <- data.frame(Stat.Test.Level=names(tab), Stat.Test.Freq=as.character(unlist(tab)), stringsAsFactors = F )
          
          list(tab)
        })
        
        # preview levels current selection for group normalization
        output$grp.norm.prev.tab <- renderTable({
          if(is.null(input$grp.norm)) return()
          if(global.param$grp.done) return()
          if(!input$grp.norm.check) return()
          
          if(input$QC.filter & "QC.status"%in%colnames(global.input$cdesc)){
            tab <- table(global.input$cdesc[global.input$cdesc$QC.status!="QC.fail", input$grp.norm])
          }else if (input$QC.filter & !"QC.status"%in%colnames(global.input$cdesc)){
            shinyalert("No QC.status column detected!","No QC.status column detected for filtering. Analysis will proceed using all samples.", type="warning")
            tab <- table(global.input$cdesc[, input$grp.norm])
          }else{
            tab <- table(global.input$cdesc[, input$grp.norm])
          }
          tab <- data.frame(Group.Norm.Level=names(tab), Group.Norm.Freq=as.character(unlist(tab)), stringsAsFactors = F )
          
          list(tab)
        })
        
        
        
        #######################################
        ## OBSERVER: define groups for GCT v1.3
        observeEvent(input$update.grp.gct3, {
          
          tab <- global.input$table
          cdesc <- data.frame(global.input$cdesc)
         
          ## store grp column
          global.param$grp.gct3 <- input$grp.gct3
          
          
          ## store grp normalization column
          if(!input$grp.norm.check | is.null(input$grp.norm.check)){
            global.param$grp.norm <- input$grp.gct3
            global.param$norm.per.group <- FALSE
          }else{
            global.param$grp.norm <- input$grp.norm
            global.param$norm.per.group <- TRUE
          }
          
          ## store grp column for PCA colors
          global.plotparam$pca.grp.col <- input$grp.gct3
          
          if(!global.param$grp.gct3 %in% colnames(cdesc)){
              error$title <- paste("Parsing error")
              error$msg <- paste("Can't find column'",global.param$grp.gct3, "'in the sample meta data. Does it contain special characters (e.g. blanks)?")
              return()
          }
          
          if(!global.param$grp.norm %in% colnames(cdesc)){
            error$title <- paste("Parsing error")
            error$msg <- paste("Can't find column'",global.param$grp.norm, "'in the sample meta data. Does it contain special characters (e.g. blanks)?")
            return()
          }
          
          ## robustify levels of group variable
          ## prevent NA from being converted (TEXT FILE INPUT ONLY)
          cdesc[!is.na(cdesc[,input$grp.gct3]), input$grp.gct3] <- make.names(cdesc[!is.na(cdesc[,input$grp.gct3]), input$grp.gct3])
          cdesc[!is.na(cdesc[,input$grp.norm]), input$grp.norm] <- make.names(cdesc[!is.na(cdesc[,input$grp.norm]), input$grp.norm])
          global.input$cdesc <- cdesc
          
          # initialize grp file
          Column.Name <- colnames(tab)
          Experiment <- rep('', length(Column.Name))
          names(Experiment) <- Column.Name
          Group <- rep('', length(Column.Name))
          names(Group) <- Column.Name
          QC <- rep('QC.pass',length(Column.Name))
          names(QC) <- Column.Name

          Experiment[ rownames(cdesc) ] <- cdesc[, global.param$grp.gct3]
          Group[ rownames(cdesc) ] <- cdesc[, global.param$grp.norm]
          #add QC.status column if it exists
          if(input$QC.filter & "QC.status"%in%colnames(cdesc)){
            QC[ rownames(cdesc) ] <- cdesc[, "QC.status"]
          }
          
          #don't need to keep QC column for analysis, just for filtering, so this should work fine
          global.param$cdesc.all <- global.param$cdesc.selection <- setdiff(colnames(cdesc),  c(global.param$grp.gct3,global.param$grp.norm))
        
          grp.file=data.frame(
            Column.Name,
            Experiment,
            Group,
            QC,
            stringsAsFactors = F
              )
          
          ## ################################
          ## ANNOTATION: extract empty cells
          ## - corresponding columns will be carried over as
          ##   annotation columns in the result file
          grp.anno <- grp.file[which(nchar( Experiment) == 0 | is.na(nchar( Experiment))), ]
          grp.anno <- setdiff( grp.anno$Column.Name, global.param$id.col.value )
          
          if(length(grp.anno)>0)
            global.input$table.anno <- data.frame(id=global.results$id.map[, 'id'], global.input$table[ , grp.anno])
          
          #replace NA strings with actual NA values
          grp.file$Experiment[grp.file$Experiment%in%NASTRINGS]=NA
          grp.file$Group[grp.file$Group%in%NASTRINGS]=NA
          
          #filter out QC.fail samples if requested
          if(input$QC.filter){
            grp.file$Experiment[grp.file$QC=="QC.fail"]=NA
            grp.file$Group[grp.file$QC=="QC.fail"]=NA
          }
          
          #remove samples with missing annotations from the table, group file, and cdesc
          grp.file = grp.file[!is.na(grp.file$Experiment) & !is.na(grp.file$Group),]
          #robustify levels of group variables
          grp.file$Experiment <- make.names(grp.file$Experiment)
          grp.file$Group <- make.names(grp.file$Group)
          global.input$table <- tab <- tab[,colnames(tab)%in%grp.file$Column.Name]
          #replace the Group annotation column if needed
          if(!input$grp.norm.check & "Group"%in%colnames(cdesc)){
            cdesc$Group <- rep("None",dim(cdesc)[1])
          }
          global.param$cdesc.all <- global.param$cdesc.selection <- cdesc[rownames(cdesc)%in%grp.file$Column.Name,colSums(is.na(cdesc))<nrow(cdesc),drop=F]
          cdesc <- global.param$cdesc.all
          
          #if an annotation only appears once, throw an error
          tab2 <- table(cdesc[, global.param$grp.gct3]) 
          tab2 <- data.frame(Stat.Test.Level=names(tab2), Stat.Test.Freq=as.character(unlist(tab2)), stringsAsFactors = F )
          if("1"%in%tab2$Stat.Test.Freq){
            error$title <- paste("All column levels must have more than one sample!")
            error$msg <- paste("At least one column level has only one sample assigned to it. Analysis cannot proceed. Please remove this sample or choose another column to use for statistical testing.")
            return()
          }
          
          if(global.param$norm.per.group){
            tab3 <- table(cdesc[, global.param$grp.norm]) 
            tab3 <- data.frame(Stat.Test.Level=names(tab3), Stat.Test.Freq=as.character(unlist(tab3)), stringsAsFactors = F )
            if("1"%in%tab3$Stat.Test.Freq){
              error$title <- paste("All column levels must have more than one sample!")
              error$msg <- paste("At least one column level has only one sample assigned to it. Analysis cannot proceed. Please remove this sample or choose another column to use for group-wise normalization.")
              return()
            }
          }
          
          ## ################################
          ## EXPRESSION
          ## - extract all non-empty cells in the 'Experiment' column
          exprs.idx <- rownames(cdesc)
          grp.exprs <- grp.file[exprs.idx, ]

          ## - extract all non-empty cells in the 'Group' column
          norm.idx <- rownames(cdesc)
          grp.norms <- grp.file[norm.idx, ]
          
          ## order alphabetically to make coloring consistent
          #grp.exprs <- grp.exprs[order(grp.exprs$Experiment), ]
          #grp.norms <- grp.norms[order(grp.norms$Group), ]
          
          ## class vector
          grp=grp.exprs$Experiment
          names(grp)=grp.exprs$Column.Name
          grp.norm=grp.norms$Group
          names(grp.norm)=grp.norms$Column.Name
          
          ## update input table, keep id and expression columns
          global.input$table <- global.input$table[ , c(global.param$id.col.value, names(grp), names(grp.norm))]
          
          ################################
          ## update number of groups
          global.param$N.grp <- length(unique( na.omit(grp)) )
          
          ## store group assignment
          global.param$grp <- global.param$grp.all <- grp
          global.param$grp.norm <- grp.norm
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
          
          ## group colors for normalization
          grp.col.norm <- rep(GRPCOLORS[1], length(grp.norm))
          names(grp.col.norm) <- names(grp.norm)
          
          for(i in 2:length(unique(grp.norm))) grp.col.norm[ which(grp.norm == unique(grp.norm)[i]) ] <- GRPCOLORS[i]
          global.param$grp.colors.norm <- grp.col.norm
          
          ## group colors for figure legend
          idx <- !duplicated(grp.norm)
          grp.col.legend.norm = grp.col.norm[idx]
          names(grp.col.legend.norm) <- grp.norm[idx]
          
          global.param$grp.colors.legend.norm <- grp.col.legend.norm
          
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
            id.idx <- grep('id', tab.colnames, ignore.case = T)
            if(length(id.idx) > 0){
                tab.colnames <- c(tab.colnames[id.idx], tab.colnames[-id.idx])
            }
            tab.colnames.names <- names(tab.colnames)
            names(tab.colnames) <- NULL
           
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

            if(is.null(input$id.col)) return() ## no id column selected
            if( !is.null(input$id.col))
                if( input$id.col == 0 ) return() ## not pressed yet
            if(!is.null(global.results$data)) return() ## test has been run
            if(!is.null(global.param$grp)){            ## group assignment has been RUN
                if(sum(is.na(global.param$grp)) == 0) return() ## group assignemnt has been DONE
            }

            list(
                ## upload template
                fileInput("exp.file", "Upload experimental design file", accept=c('text/plain','.txt'))
            )
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
            ##                file import
            ########################################################
            
            
            ################################################
            ##                      GCTX
            if(grepl('\\.gctx$', fn)){
                
                gct <- try( parse.gctx(fn) )
                
                #################################################
                ##                     GCT 1.2
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
                
                #####################################################      
                ##                     GCT 1.3
            } else if( length( grep( '^\\#1\\.3', readLines(fn,n=1))) > 0){
                
                # parse gct file
                gct <- try( parse.gctx2(fn, show_modal = TRUE) )
                if(class(gct) == 'try-error'){
                    error$title <- "Error importing GCT 1.3 file"
                    error$msg <- gct[1]
                    
                    validate(need(class(gct) != 'try-error', 'Error importing GCT 1.3 file.'))
                }
                
                ## #################################################
                ## robustify rids
                if(sum(duplicated(gct@rid)) > 0)
                    gct@rid <- de_duplicate_ids(gct@rid, reactiveValuesToList(global.param))
                rownames(gct@mat) <- gct@rid
                
                ## rdesc
                if(nrow(gct@rdesc) > 0){
                    rownames(gct@rdesc) <- gct@rid   
                }
                ## cid
                gct@cid <- make.unique(make.names(gct@cid))
                #cdesc
                if(nrow(gct@cdesc) > 0)
                    rownames(gct@cdesc) <- gct@cid
                colnames(gct@mat) <- gct@cid
                
                ## error checking for column meta data
                if(nrow(gct@cdesc) == 0){
                    error$title <- "Error parsing GCT 1.3 column meta data"
                    error$msg <- paste('No column meta data tracks defined! Need at least one column meta data track to use as class vector.') 
                    validate( need(nrow(gct@cdesc) > 0, 'Error parsing GCT 1.3 column meta data.'))
                }
                ## robustify cdesc column names
                colnames(gct@cdesc) <- make.names(colnames(gct@cdesc))
                
                ## remove 'id' from cdesc
                if('id' %in% colnames(gct@cdesc)){
                    cn.tmp <- colnames(gct@cdesc)
                    rm.idx <- which(colnames(gct@cdesc) == 'id')
                    gct@cdesc <- data.frame(gct@cdesc[ ,-rm.idx] )
                    cn.tmp <- cn.tmp[-rm.idx]
                    colnames(gct@cdesc) <- cn.tmp
                }
                if(ncol(gct@cdesc) == 1)
                    rownames(gct@cdesc) <- gct@cid
                
                #####################
                # expression table
                if(nrow(gct@rdesc) > 0 ){
                    tab <- data.frame(id=gct@rid, gct@rdesc, gct@mat, stringsAsFactors = F)
                } else {
                    tab <- data.frame(id=gct@rid, gct@mat, stringsAsFactors = F)
                }
                rownames(tab) <- tab$id
                
                ## sample names
                colnames.tmp <- chopString(colnames(tab), STRLENGTH)
                names(colnames.tmp) <- colnames(tab)
                
                # id column 
                global.param$id.col.value='id'
                global.param$id.done=T
                
                ## #################
                ## map to gene names
                map.res <- mapIDs(tab$id,gct@rdesc)
                global.results$keytype <- map.res$keytype
                global.results$id.map <- map.res$id.map
                
                ## store values
                global.input$table <- global.input$table.org <- tab
                global.input$file <- input$file
                global.input$table.colnames <- colnames.tmp
                
                global.input$rdesc <- gct@rdesc
                global.input$cdesc <- gct@cdesc
                
                global.param$cdesc.all <- global.param$cdesc.selection <- colnames(global.input$cdesc)
                
                # flag
                global.param$file.gct3 <- T
                
                rm(tab, colnames.tmp)
                
                # ##########################################################  
                #                    other text file
            } else { ## end if GCT 1.3
                
                ## ################################
                ## determine the separator
                ## try to figure out the separator, DON'T USE THE HEADER FOR THAT
                ## use the fourth row instead (should be data)
                for(s in SEPARATOR){
                    
                    tab <- read.table(fn, sep=s, header=F, stringsAsFactors=F, nrows=1, skip=4)
                    
                    if(length(tab) > 1){
                        global.param$tabsep <- s
                        break;
                    }
                }
                
                ## #########################################################
                ## import the table
                if( global.param$tabsep == '\t'){
                    tab <- read.delim( fn, stringsAsFactors=F, na.strings=NASTRINGS)
                } else {
                    tab <- read.table( fn, sep=global.param$tabsep, header=T, 
                                       stringsAsFactors=F, na.strings=NASTRINGS, quote = "\"", 
                                       dec = ".", fill = TRUE, comment.char = "")
                }
                
                ## shorten column names and store together with the original names
              colnames.tmp <- chopString(colnames(tab), STRLENGTH)
              names(colnames.tmp) <- colnames(tab)
                
                ## store values
                global.input$table <- global.input$table.org <- tab ## need to make sure that 'global.input$table.org' never gets overwritten 
                global.input$file <- input$file
                global.input$table.colnames <- colnames.tmp
                
                rm(tab, colnames.tmp)
            } 
            
            # flag
            global.param$file.done <- T
        })
        
        ###############################################
        ## observer
        ## 4) ID column
        ##  - make unique ids
        ##  - determine id type
        ##  - map to gene names
        ##  - initialize group assignment
        ##
        ## - for txt and gct 1.2 file formats
        observeEvent( input$id.col ,{
            
            if( is.null( global.input$table) | is.null(input$id.col.value) ) return()
            
            ## store name of id column
            global.param$id.org.col <- global.param$id.col.value <- input$id.col.value
            cat('selected id column: ', global.param$id.col.value, '\n')
            
            ## update 'global.input$id.col'
            global.input$id.col <- input$id.col
            
            ## ###########################################
            ## check the id column
            tab <- global.input$table
            
            ## ###########################################
            ## make sure the ids are unique
            ids <- as.character(tab[, global.param$id.col.value])
            if(sum(duplicated(ids)) > 0)
                ids <- de_duplicate_ids(ids, reactiveValuesToList(global.param))
            
            #####################################    
            ## - update 'id.col.value' with 'id'
            ## - replace values in 'id' column
            global.param$id.col.value.updated <- FALSE
            global.param$id.col.renamed <- FALSE
            
            if(global.param$id.col.value != 'id'){
                
                cat('updated id column: ', global.param$id.col.value, ' --> ')
                cat('id\n\n')
                
                if('id' %in% colnames(tab)){
                    colnames(tab)[ colnames(tab) == 'id' ] <- 'id_org'
                    cat('renamed column "id" --> "id_org"\n\n')
                    
                    global.param$id.col.renamed <- TRUE
                }
                    
                tab <- data.frame(id=ids, tab)
                
                global.param$id.col.value <- "id"
                global.param$id.col.value.updated <- TRUE
                
               # global.input$table.colnames <- colnames(tab)

            }
            tab[, 'id'] <- ids
            
            ## use id as rownames
            rownames(tab) <- ids
            
            ## ############################################
            ## map to gene names
            ## ############################################
            map.res <- mapIDs(ids)
            global.results$keytype <- map.res$keytype
            global.results$id.map <- map.res$id.map
            
            ########################################
            ## store
            global.input$table <- tab
            
            
            # #######################
            # flag
            global.param$id.done <- T
            
        })
        
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
            
            ## ################################
            ## determine the separator
            ## try to figure out the separator, DON'T USE THE HEADER FOR THAT
            ## use the fourth row instead (should be data)
            for(s in SEPARATOR){
              
              grp.file<- read.table(input$exp.file$datapath, sep=s, header=F, stringsAsFactors=F, nrows=1, skip=4)
              
              if(length(grp.file) > 1){
                global.param$tabsep.anno <- s
                break;
              }
            }
            
            ## #########################################################
            ## import the table: txt
            if( global.param$tabsep.anno == '\t'){
              grp.file<- read.delim( input$exp.file$datapath, stringsAsFactors=F, na.strings=NASTRINGS)
            } else {
              grp.file <- read.table( input$exp.file$datapath, sep=global.param$tabsep.anno, header=T, 
                                 stringsAsFactors=F, na.strings=NASTRINGS, quote = "\"", 
                                 dec = ".", fill = TRUE, comment.char = "")
            }
            
            #remove the ID column (if present)
            colnames.noid <- colnames(global.input$table)[!colnames(global.input$table)%in%global.param$id.col.value & !colnames(global.input$table)%in%global.param$id.org.col]
            grp.file <- grp.file[grp.file[,1]%in%colnames.noid,]
            
            #check that the first column contains the column names of the table
            #if not, throw an error
            if(sum(colnames.noid %in% grp.file[,1])!=length(colnames.noid)){
              error$title <- paste("Parsing error")
              error$msg <- paste("The first column of the experimental design file must be the sample (column) names, and they must match the column names in the table exactly! All columns (except the ID column) must be present in the experimental design file.")
              return()
            }
            
            
            #after this, the grp.file will be used as the cdesc file, and non-gct table will be processed as a gct file.
            #assume that the first column of the annotation file is the id
            my.cdesc <- data.frame(grp.file[,-1])
            row.names(my.cdesc) <- grp.file[,1]
            global.input$cdesc <- my.cdesc
            global.param$file.gct3 <- T
            
            ###############################################################
            
            #grp.file <- read.delim(input$exp.file$datapath, header=T, stringsAsFactors=F)
            
            ##OLD ANNOTATION FILE IMPORT CODE
            # ############################# 
            # ## if existing 'id' column has been renamed to 'id_org'
            # if(global.param$id.col.renamed){
            #     
            #     id_idx <- which(grp.file$Column.Name == 'id')
            #     grp.file$Column.Name[id_idx] <- 'id_org'
            # }
            
            # #############################
            # ## if id column has been updated
            # ## to 'id':
            # if(global.param$id.col.value.updated){
            #     id_row <- data.frame(Column.Name='id', Experiment='', Group="")
            #     grp.file <- rbind(id_row, grp.file) 
            # }
            
            
            # Column.Name <- grp.file$Column.Name
            # Experiment <- as.character(grp.file$Experiment)
            # Group <- as.character(grp.file$Group)
            # 
            # ## #############################################################
            # ## index on non-empty 'Experiment' rows
            # exprs.idx <- which(nchar(Experiment) > 0 )
            # norm.idx <- which(nchar(Group) > 0)
            # 
            # ## ###############################
            # ## update label
            # fn.split <- unlist(strsplit( sub('.*/', '', fn), '_'))
            # if(length(fn.split) > 1){
            #     label <- fn.split[1]
            #     if(nchar(label) > 10)
            #         label <- chopString(label, 10, add.dots=F)
            #     global.param$label <- paste(global.param$label, label )
            # } else {
            #     label <- chopString( sub('.*/', '', fn) , 10, add.dots=F)
            #     global.param$label <- paste(global.param$label, label, sep='-')
            # }
            # 
            # ###############################################################
            # ##
            # ## - do some sanity checks
            # ## - separate expression data from annotation columns
            # ##
            # ###############################################################
            # 
            # ## Number of rows in exp design file does not match number of columns in data file
            # if( length(Column.Name) != ncol(global.input$table)  ){
            #     error$title <- "Problem parsing experimental design file."
            #     error$msg <- 'Experimental design file does not match the table you have uploaded (different number of rows/columns)!'
            #     return()
            # }
            # #table <- global.input$table
            # #save(Column.Name, table, file='debug.RData')
            # ## names in the exp design file do not match to the table
            # if( sum( Column.Name != colnames(global.input$table)) ){
            #   
            #     error$title <- "Problem parsing experimental design file."
            #     error$msg <- 'Experimental design file does not match the table you have uploaded!'
            #     return()
            # }
            # 
            # ## not an experimental design file
            # if( sum( colnames(grp.file) %in% c('Column.Name', 'Experiment'), na.rm=T) < 2 )  {
            #     error$title <- "Problem parsing experimental design file."
            #     error$msg <- 'This is not an experimental design file! The file should contain at least two columns named (Column.Name, Experiment)!'
            #     return()
            # }
            # ## 'empty' file
            # if( sum( nchar(Experiment) > 0, na.rm=T ) == 0 | sum(!is.na( Experiment) == 0) ){
            #     error$title <- "Problem parsing experimental design file."
            #     error$msg <- 'No experiments defined!'
            #     return()
            # }
            # ##cat('L=', sum( Column.Name[ exprs.idx ] %in%  colnames(global.input$table)))
            # ## column names specified in exp design file not found in table
            # ##if( sum( Column.Name[ exprs.idx ] %in%  colnames(global.input$table)) != length(exprs.idx) ){
            # ##    error$msg <- 'Column names in the experimental design file cannot be found in the data table!'
            # ##    return()
            # ## }
            # ## check whether there are at least 2 replicates per group
            # num.rep=table(Experiment[exprs.idx])
            # if(min(num.rep) == 1){
            #   error$title <- "Not enough replicates."
            #         error$msg <- paste('Warning! All groups must have at least 2 replicates!')
            #         return()
            # }
            # 
            # ## ################################
            # ## ANNOTATION: extract empty cells
            # ## - corresponding columns will be carried over as
            # ##   annotation columns in the result file
            # grp.anno <- grp.file[which(nchar( Experiment) == 0 ), ]
            # grp.anno <- setdiff( grp.anno$Column.Name, global.param$id.col.value )
            # 
            # if(length(grp.anno)>0)
            #     global.input$table.anno <- data.frame(id=global.results$id.map[, 'id'], global.input$table[ , grp.anno])
            # 
            # ## ################################
            # ## EXPRESSION
            # ## - extract all non-empty cells in the 'Experiment' column
            # grp.exprs <- grp.file[exprs.idx, ]
            # grp.norms <- grp.file[norm.idx, ]
            # 
            # ## order alphabetically to make coloring consistent
            # grp.exprs <- grp.exprs[order(grp.exprs$Experiment), ]
            # grp.norms <- grp.norms[order(grp.norms$Group), ]
            # 
            # ## class vector
            # grp=grp.exprs$Experiment
            # grp.norm = grp.exprs$Group
            # ## robustify experiment names
            # grp <- gsub('^ {1,10}', '', grp)
            # grp <- gsub(' {1,10}$', '', grp)
            # grp <- make.names(grp)
            # names(grp)=grp.exprs$Column.Name
            # 
            # grp.norm <- gsub('^ {1,10}', '', grp.norm)
            # grp.norm <- gsub(' {1,10}$', '', grp.norm)
            # grp.norm <- make.names(grp.norm)
            # names(grp.norm)=grp.exprs$Column.Name
            # 
            # ## update input table, keep id and expression columns
            # #tab <- global.input$table[ , c(global.param$id.col.value, names(grp))]
            # tab <- global.input$table[ , c('id', names(grp))]
            # 
            # ## #################################
            # ## remove NA rows
            # na.row.idx <- apply(tab[, names(grp)], 1, function(x) sum(is.na(x))/length(x) )
            # na.row.idx <- which(na.row.idx == 1)
            # if(length(na.row.idx) > 0){
            #     tab <- tab[-na.row.idx, ]
            #     global.input$NA.rows <- length(na.row.idx)
            # }
            # 
            # global.input$table <- tab
            # 
            # ################################
            # ## update number of groups
            # global.param$N.grp <- length(unique( na.omit(grp)) )
            # 
            # ## store group assignment
            # global.param$grp <- global.param$grp.all <- grp
            # global.param$grp.norm <- grp.norm
            # global.param$grp.comp.all <- global.param$grp.comp <- unique(grp)
            # 
            # ## group colors
            # grp.col <- rep(GRPCOLORS[1], length(grp))
            # names(grp.col) <- names(grp)
            # 
            # for(i in 2:length(unique(grp))) grp.col[ which(grp == unique(grp)[i]) ] <- GRPCOLORS[i]
            # global.param$grp.colors <- grp.col
            # 
            # ## group colors for figure legend
            # idx <- !duplicated(grp)
            # grp.col.legend = grp.col[idx]
            # names(grp.col.legend) <- grp[idx]
            # global.param$grp.colors.legend <- grp.col.legend
            # 
            # ## group colors for normalization
            # grp.col.norm <- rep(GRPCOLORS[1], length(grp.norm))
            # names(grp.col.norm) <- names(grp.norm)
            # 
            # for(i in 2:length(unique(grp.norm))) grp.col.norm[ which(grp.norm == unique(grp.norm)[i]) ] <- GRPCOLORS[i]
            # global.param$grp.colors.norm <- grp.col.norm
            # 
            # ## group colors for figure legend
            # idx <- !duplicated(grp.norm)
            # grp.col.legend.norm = grp.col.norm[idx]
            # names(grp.col.legend.norm) <- grp.norm[idx]
            # 
            # global.param$grp.colors.legend.norm <- grp.col.legend.norm
            # 
            # ## all done
            # global.param$grp.done = T
            # 
            # ## save column name used as 'id'
            # ##global.param$id.col.value = input$id.col.value
            # 
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
                  conditionalPanel(condition = "input['filter.type'] == 'nom.p'", numericInput( "filter.value.nom.p", "p-value filter", value=0.01, min=0, max=1, step=1e-2)),
                  conditionalPanel(condition = "input['filter.type'] == 'adj.p'", numericInput( "filter.value.adj.p", "Adj.p-value", value=0.05, min=0, max=1, step=1e-2))
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
                                                  choices = unique(colnames(global.param$cdesc.all)), 
                                                  selected = unique(colnames(global.param$cdesc.all))
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
                                   selected = setdiff(unique(colnames(global.param$cdesc.all)), input$select.anno)
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
            shinyalert("No data selected!", "Please select at least one group.", type = "error")
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
            global.param$anno.col <- anno.col[names(global.param$grp.selection), colnames(anno.col)%in%unique( c(cdesc.selection, global.param$grp.gct3)) ]
            global.param$anno.col.color <- anno.col.color[ unique(c(cdesc.selection, global.param$grp.gct3 )) ]
            if(is.character(global.param$anno.col)){
              global.param$anno.col <- data.frame(global.param$anno.col)
              colnames(global.param$anno.col) <- colnames(anno.col)[colnames(anno.col)%in%unique( cdesc.selection, global.param$grp.gct3)]
            }
          }
          shinyalert("Group Selection Updated!", "Press OK to close this window and proceed with analysis.", type = "success")
          
        })
        
       
        ## #####################################################################
        ## UI: set up analysis
        ##
        ## - the if-conditions regarding the type of filtering below are
        ##   probably obsolete and can be removed
        ## - initial reason to have the if-statements was to avoid the  
        ##   application of the reproducibility filter to 2-sample t and F-test
        ##
        ## #####################################################################
        output$list.groups <- renderUI({

            if( !global.param$grp.done ) return()

            ####################################
            ## initialize with default values from
            ## 'global.param'
            ## - show everything
            ## - only visible for a fraction of a second (until input$filt.data is set)
            ##  
            if( is.null(input$filt.data)){

                list(
                  checkboxInput('intensity','Intensity data',FALSE),
                     radioButtons('log.transform', 'Log-transformation', choices=c('none', 'log10', 'log2'), selected=global.param$log.transform),              #checkboxInput('norm.per.group', 'Normalize per group', value = global.param$norm.per.group ),
                     radioButtons('norm.data', 'Data normalization', choices=c('Median', 'Median (non-zero)', 'Median-MAD', 'Median-MAD (non-zero)', 'Upper-quartile', '2-component', 'Quantile', 'VSN', 'none'), selected=global.param$norm.data),                  #checkboxInput('norm.per.group', 'Normalize per group', value = global.param$norm.per.group),
                     #sliderInput('na.filt.val', 'Max. % missing values', min=0, max=100, value=global.param$na.filt.val),
                     numericInput('na.filt.val', 'Max. % missing values', min=0, max=100, step=5, value=global.param$na.filt.val),
                  actionButton("exp","Click for missing value info"), 
                  tippy_this("exp",tooltip="<span style='font-size:14px;'>Example: A value of 70 means a feature may have up to 70% missing values (quantified in at least 30% of samples). If you do not want to filter missing values, leave the value at 100 (99 for intensity-based data).<span>",allowHTML=TRUE,placement="auto",trigger="click",theme="light"),
                    
                     radioButtons('filt.data', 'Filter data', choices=c('Reproducibility', 'StdDev', 'none'), selected=global.param$filt.data ),
                     #radioButtons('which.test', 'Select test', choices=c('One-sample mod T', 'Two-sample mod T', 'Two-sample LM', 'mod F', 'none'), selected=global.param$which.test),
                     radioButtons('which.test', 'Select test', choices=c('One-sample mod T', 'Two-sample mod T', 'mod F', 'none'), selected=global.param$which.test),
                     
                     actionButton('run.test', 'Run analysis!'),
                     br(),
                     hr(),
                     actionButton('select.groups.button', 'Selected groups')
      
                )
            }
            ## ###################################################
            ## no filter
            ## - show everything
            else if(input$filt.data == 'none' & !input$intensity){

                list(
                  checkboxInput('intensity','Intensity data',FALSE),
                     #radioButtons('log.transform', 'Log-transformation', choices=c('none', 'log10', 'log2'), selected=global.param$log.transform),
                     radioButtons('log.transform', 'Log-transformation', choices=c('none', 'log10', 'log2'), selected=input$log.transform),             #checkboxInput('norm.per.group', 'Normalize per group', value = input$norm.per.group),
                     
                     #radioButtons('norm.data', 'Data normalization', choices=c('Median', 'Median (non-zero)', 'Median-MAD', 'Median-MAD (non-zero)', 'Upper-quartile', '2-component', 'Quantile', 'VSN', 'none'), selected=global.param$norm.data),
                     radioButtons('norm.data', 'Data normalization', choices=c('Median', 'Median (non-zero)', 'Median-MAD', 'Median-MAD (non-zero)', 'Upper-quartile', '2-component', 'Quantile', 'VSN', 'none'), selected=input$norm.data),
                     #checkboxInput('norm.per.group', 'Normalize per group', value = input$norm.per.group),
                     #sliderInput('na.filt.val', 'Max. % missing values', min=0, max=100, value=global.param$na.filt.val),
                     #sliderInput('na.filt.val', 'Max. % missing values', min=0, max=100, value=input$na.filt.val),
                  numericInput('na.filt.val', 'Max. % missing values', min=0, max=100, step=5, value=input$na.filt.val),
                  actionButton("exp","Click for missing value info"), 
                  tippy_this("exp",tooltip="<span style='font-size:14px;'>Example: A value of 70 means a feature may have up to 70% missing values (quantified in at least 30% of samples). If you do not want to filter missing values, leave the value at 100 (99 for intensity-based data).<span>",allowHTML=TRUE,placement="auto",trigger="click",theme="light"),
                     
                     radioButtons('filt.data', 'Filter data', choices=c('Reproducibility', 'StdDev', 'none'), selected='none'),
                     
                     #radioButtons('which.test', 'Select test', choices=c('One-sample mod T', 'Two-sample mod T', 'Two-sample LM', 'mod F', 'none'), selected=global.param$which.test),
                      
                     #radioButtons('which.test', 'Select test', choices=c('One-sample mod T', 'Two-sample mod T', 'mod F', 'none'), selected=global.param$which.test),
                     radioButtons('which.test', 'Select test', choices=c('One-sample mod T', 'Two-sample mod T', 'mod F', 'none'), selected=input$which.test),
                     
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
          ## no filter + intensity data
          ## - only show relevant normalization+filtering methods and statistical tests
          ## max missing values is 99%, otherwise statistics fails
          else if(input$filt.data == 'none' & input$intensity){
            
            list(
              checkboxInput('intensity','Intensity data',TRUE),
              #radioButtons('log.transform', 'Log-transformation', choices=c('none', 'log10', 'log2'), selected=global.param$log.transform),
              radioButtons('log.transform', 'Log-transformation', choices=c('none', 'log10', 'log2'), selected=input$log.transform),                    #checkboxInput('norm.per.group', 'Normalize per group', value = input$norm.per.group),
              
              #radioButtons('norm.data', 'Data normalization', choices=c('Median', 'Median (non-zero)', 'Median-MAD', 'Median-MAD (non-zero)', 'Upper-quartile', '2-component', 'Quantile', 'VSN', 'none'), selected=global.param$norm.data),
              radioButtons('norm.data', 'Data normalization', choices=c('Median (non-zero)', 'Median-MAD (non-zero)','Upper-quartile', 'Quantile','VSN', 'none'), selected=input$norm.data),
              #checkboxInput('norm.per.group', 'Normalize per group', value = input$norm.per.group),
              #sliderInput('na.filt.val', 'Max. % missing values', min=0, max=100, value=global.param$na.filt.val),
              #sliderInput('na.filt.val', 'Max. % missing values', min=0, max=100, value=input$na.filt.val),
              numericInput('na.filt.val', 'Max. % missing values', min=0, max=99, step=3, value=min(input$na.filt.val,99)),
              actionButton("exp","Click for missing value info"), 
              tippy_this("exp",tooltip="<span style='font-size:14px;'>Example: A value of 70 means a feature may have up to 70% missing values (quantified in at least 30% of samples). If you do not want to filter missing values, leave the value at 100 (99 for intensity-based data).<span>",allowHTML=TRUE,placement="auto",trigger="click",theme="light"),
              
              radioButtons('filt.data', 'Filter data', choices=c( 'StdDev', 'none'), selected='none'),
              
              #radioButtons('which.test', 'Select test', choices=c('One-sample mod T', 'Two-sample mod T', 'Two-sample LM', 'mod F', 'none'), selected=global.param$which.test),
              
              #radioButtons('which.test', 'Select test', choices=c('One-sample mod T', 'Two-sample mod T', 'mod F', 'none'), selected=global.param$which.test),
              radioButtons('which.test', 'Select test', choices=c('Two-sample mod T', 'mod F', 'none'), selected=input$which.test),
              
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
            ## - show only one sample T test
            ## no option for intensity data with reproducibility filter
            else if(input$filt.data == 'Reproducibility'){

                list(
                  #checkboxInput('intensity','Intensity data',FALSE),
                    #radioButtons('log.transform', 'Log-transformation', choices=c('none', 'log10', 'log2'), selected=global.param$log.transform),
                    radioButtons('log.transform', 'Log-transformation', choices=c('none', 'log10', 'log2'), selected=input$log.transform),                    #checkboxInput('norm.per.group', 'Normalize per group', value = input$norm.per.group),
                    
                    #radioButtons('norm.data', 'Data normalization', choices=c('Median', 'Median (non-zero)', 'Median-MAD', 'Median-MAD (non-zero)', 'Upper-quartile', '2-component', 'Quantile', 'VSN', 'none'), selected=global.param$norm.data),
                    radioButtons('norm.data', 'Data normalization', choices=c('Median', 'Median (non-zero)', 'Median-MAD', 'Median-MAD (non-zero)', 'Upper-quartile', '2-component', 'Quantile', 'VSN', 'none'), selected=input$norm.data),
                    
                    #sliderInput('na.filt.val', 'Max. % missing values', min=0, max=100, value=global.param$na.filt.val),
                    #sliderInput('na.filt.val', 'Max. % missing values', min=0, max=100, value=input$na.filt.val),
                    numericInput('na.filt.val', 'Max. % missing values', min=0, max=100, step=5, value=input$na.filt.val),
                    actionButton("exp","Click for missing value info"), 
                    tippy_this("exp",tooltip="<span style='font-size:14px;'>Example: A value of 70 means a feature may have up to 70% missing values (quantified in at least 30% of samples). If you do not want to filter missing values, leave the value at 100 (99 for intensity-based data).<span>",allowHTML=TRUE,placement="auto",trigger="click",theme="light"),
                    
                    radioButtons('filt.data', 'Filter data', choices=c('Reproducibility', 'StdDev', 'none'), selected='Reproducibility'),
                    
                    selectInput('repro.filt.val', 'alpha', choices=c(.1, .05, 0.01, 0.001 ), selected=global.param$repro.filt.val),
                     
                    # sliderInput('na.filt.val', 'Max. % missing values', min=0, max=100, value=global.param$na.filt.val),
                     # radioButtons('which.test', 'Select test', choices=c('One-sample mod T', 'none'), selected='One-sample mod T'),
                    #radioButtons('which.test', 'Select test', choices=c('One-sample mod T', 'Two-sample mod T',  'mod F', 'none'), selected=global.param$which.test),
                    radioButtons('which.test', 'Select test', choices=c('One-sample mod T', 'none'), selected=input$which.test),
                    

                     actionButton('run.test', 'Run analysis!'),
                     br(),
                     hr(),
                     actionButton('select.groups.button', 'Select groups')
                )
            }

            ## ###################################################
            ## StdDev filter
            ## - show everything
            else if(input$filt.data == 'StdDev' & !input$intensity){

                list(
                  checkboxInput('intensity','Intensity data',FALSE),
                    #radioButtons('log.transform', 'Log-transformation', choices=c('none', 'log10', 'log2'), selected=global.param$log.transform),
                    radioButtons('log.transform', 'Log-transformation', choices=c('none', 'log10', 'log2'), selected=input$log.transform),                    #checkboxInput('norm.per.group', 'Normalize per group', value = input$norm.per.group),
                    
                    #radioButtons('norm.data', 'Data normalization', choices=c('Median', 'Median (non-zero)', 'Median-MAD', 'Median-MAD (non-zero)', 'Upper-quartile', '2-component', 'Quantile', 'VSN', 'none'), selected=global.param$norm.data),
                    radioButtons('norm.data', 'Data normalization', choices=c('Median', 'Median (non-zero)', 'Median-MAD', 'Median-MAD (non-zero)', 'Upper-quartile', '2-component', 'Quantile', 'VSN', 'none'), selected=input$norm.data),
                    #sliderInput('na.filt.val', 'Max. % missing values', min=0, max=100, value=global.param$na.filt.val),
                    #sliderInput('na.filt.val', 'Max. % missing values', min=0, max=100, value=input$na.filt.val),
                    numericInput('na.filt.val', 'Max. % missing values', min=0, max=100, step=5, value=input$na.filt.val),
                  actionButton("exp","Click for missing value info"), 
                  tippy_this("exp",tooltip="<span style='font-size:14px;'>Example: A value of 70 means a feature may have up to 70% missing values (quantified in at least 30% of samples). If you do not want to filter missing values, leave the value at 100 (99 for intensity-based data).<span>",allowHTML=TRUE,placement="auto",trigger="click",theme="light"),
                    
                    radioButtons('filt.data', 'Filter data', choices=c('Reproducibility', 'StdDev', 'none'), selected='StdDev'),
                    
                    sliderInput('sd.filt.val', 'Percentile StdDev', min=10, max=90, value=global.param$sd.filt.val),

                    #sliderInput('na.filt.val', 'Max. % missing values', min=0, max=100, value=global.param$na.filt.val),
                    
                    #radioButtons('which.test', 'Select test', choices=c('One-sample mod T', 'Two-sample mod T', 'Two-sample LM', 'mod F', 'none'), selected=global.param$which.test),
                    #radioButtons('which.test', 'Select test', choices=c('One-sample mod T', 'Two-sample mod T', 'mod F', 'none'), selected=global.param$which.test),
                    radioButtons('which.test', 'Select test', choices=c('One-sample mod T', 'Two-sample mod T', 'mod F', 'none'), selected=input$which.test),
                    actionButton('run.test', 'Run analysis!'),
                    br(),
                    hr(),
                    actionButton('select.groups.button', 'Select groups')
                )
            }
          
          ## ###################################################
          ## StdDev filter plus intensity data
          ## - show only appropriate normalization+filtering methods and statistical tests
          ## max missing values is 99%, otherwise statistics fails
          else if(input$filt.data == 'StdDev' & input$intensity){
            
            list(
              checkboxInput('intensity','Intensity data',TRUE),
              #radioButtons('log.transform', 'Log-transformation', choices=c('none', 'log10', 'log2'), selected=global.param$log.transform),
              radioButtons('log.transform', 'Log-transformation', choices=c('none', 'log10', 'log2'), selected=input$log.transform),                    #checkboxInput('norm.per.group', 'Normalize per group', value = input$norm.per.group),
              
              #radioButtons('norm.data', 'Data normalization', choices=c('Median', 'Median (non-zero)', 'Median-MAD', 'Median-MAD (non-zero)', 'Upper-quartile', '2-component', 'Quantile', 'VSN', 'none'), selected=global.param$norm.data),
              radioButtons('norm.data', 'Data normalization', choices=c('Median (non-zero)', 'Median-MAD (non-zero)','Upper-quartile', 'Quantile','VSN', 'none'), selected=input$norm.data),
              #sliderInput('na.filt.val', 'Max. % missing values', min=0, max=100, value=global.param$na.filt.val),
              #sliderInput('na.filt.val', 'Max. % missing values', min=0, max=100, value=input$na.filt.val),
              numericInput('na.filt.val', 'Max. % missing values', min=0, max=99, step=3, value=min(input$na.filt.val,99)),
              actionButton("exp","Click for missing value info"), 
              tippy_this("exp",tooltip="<span style='font-size:14px;'>Example: A value of 70 means a feature may have up to 70% missing values (quantified in at least 30% of samples). If you do not want to filter missing values, leave the value at 100 (99 for intensity-based data).<span>",allowHTML=TRUE,placement="auto",trigger="click",theme="light"),
              
              radioButtons('filt.data', 'Filter data', choices=c( 'StdDev', 'none'), selected='StdDev'),
              
              sliderInput('sd.filt.val', 'Percentile StdDev', min=10, max=90, value=global.param$sd.filt.val),
              
              #sliderInput('na.filt.val', 'Max. % missing values', min=0, max=100, value=global.param$na.filt.val),
              
              #radioButtons('which.test', 'Select test', choices=c('One-sample mod T', 'Two-sample mod T', 'Two-sample LM', 'mod F', 'none'), selected=global.param$which.test),
              #radioButtons('which.test', 'Select test', choices=c('One-sample mod T', 'Two-sample mod T', 'mod F', 'none'), selected=global.param$which.test),
              radioButtons('which.test', 'Select test', choices=c('Two-sample mod T', 'mod F', 'none'), selected=input$which.test),
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
          global.plotparam$pca.grp.col <- input$pca.grp.col
          
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
          
          ## profile plot
          global.plotparam$profile.xlim.mode <- input$profile.xlim.mode
          
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
          updateSelectInput(session, inputId='pca.grp.col', selected=global.plotparam$pca.grp.col)
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
          ## profile plot
          updateSelectInput(session, inputId='profile.xlim.mode', selected=global.plotparam$profile.xlim.mode)
          
          
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
          ## data for heatmap
          hm.res <- res
          
          ## groups to compare
          grp.comp <- unique( global.param$grp.comp )
          
          ## ####################################
          ##         what to export
          export.volc <- input$export.volc & !(global.param$which.test %in% c('mod F', 'none'))
          export.cm <- input$export.cm
          export.cb <- input$export.cb
          export.hm <- input$export.hm
          export.pca <- input$export.pca
          export.box <- input$export.box
          

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
            \n# `r global.param$label` - analysis report
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
                      ',ifelse(export.cb, '\n     - [Correlation boxplot](#corrbox)','' ),'
                      ',ifelse(export.box, '\n    - [Box-and-whisker plots](#boxplot)','' ),'
                       ', sep='')
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
                       \n Summary about the processing steps specifified in Protigy.
                       \n```{r workflow, echo=F}
                       \nwf.tab.ids <- c('Intensity data', 'Group-wise normalization', 'Filter QC.fail','Log scale', 'Normalization', 'Filter data', 'Test', 'Filter results')
                       \n## Reproducibility filter
                       \nif(global.param$filt.data == 'Reproducibility'){
                       \nwf.tab <- t(data.frame( global.param$intensity, global.param$norm.per.group, global.param$QC.filter,global.param$log.transform, global.param$norm.data, paste( global.param$filt.data, ' (alpha=',global.param$repro.filt.val, ')',sep=''), global.param$which.test,
                       \npaste( global.param$filter.type, ' < ', global.param$filter.value)))
                       \nwf.tab <- data.frame(id=wf.tab.ids, value=wf.tab, stringsAsFactors=F)
                       \n}
                       \n## SD filter
                       \nif(global.param$filt.data == 'StdDev'){
                       \nwf.tab <- t(data.frame( global.param$intensity, global.param$norm.per.group, global.param$QC.filter,global.param$log.transform, global.param$norm.data, paste( global.param$filt.data, ' (SD=',global.param$sd.filt.val, '%)',sep=''), global.param$which.test,
                       \npaste( global.param$filter.type, ' < ', global.param$filter.value) ))
                       \nwf.tab <- data.frame(id=wf.tab.ids, wf.tab, stringsAsFactors=F)
                       \n}
                       \n## no data filter
                       \nif(global.param$filt.data == 'none'){
                       \nwf.tab <- t(data.frame( global.param$intensity, global.param$norm.per.group, global.param$QC.filter,global.param$log.transform, global.param$norm.data, paste( global.param$filt.data, sep=''), global.param$which.test,
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
                      \nThe barchart shows the number of non-missing data points in each data column, colored bu sample groups.
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
                              ######################################
                           \n# heatmap title
                           \nhm.title <- paste('filter:', global.param$filter.type, ' / cutoff:', global.param$filter.value, sep='')
                           \nhm.title <- paste(hm.title, '\nsig / total: ', nrow(res), ' / ', nrow( global.results$data$output ), sep='')

                           \nif(!is.null(global.param$anno.col)){
                           \n  anno.col=global.param$anno.col
                           \n     anno.col.color=global.param$anno.col.color
                           \n} else {
                           \n  anno.col=data.frame(Group=global.param$grp)
                           \n  anno.col.color=list(Group=global.param$grp.colors.legend)
                           \n}
                           \nif(global.plotparam$hm.max){
                           \n  plotHM(res=hm.res, hm.rownames=hm.rownames, grp=global.param$grp, grp.col=global.param$grp.colors, grp.col.legend=global.param$grp.colors.legend,  hm.clust=global.plotparam$hm.clust, hm.title=hm.title, hm.scale=global.plotparam$hm.scale , fontsize_row= global.plotparam$hm.cexRow, fontsize_col= global.plotparam$hm.cexCol, max.val=global.plotparam$hm.max.val, style=global.param$which.test, anno.col=anno.col, anno.col.color=anno.col.color, show.rownames=global.plotparam$hm.show.rownames, show.colnames=global.plotparam$hm.show.colnames)
                           \n} else {
                           \n   plotHM(res=hm.res, hm.rownames=hm.rownames, grp=global.param$grp, grp.col=global.param$grp.colors, grp.col.legend=global.param$grp.colors.legend,  hm.clust=global.plotparam$hm.clust, hm.title=hm.title, hm.scale=global.plotparam$hm.scale , fontsize_row= global.plotparam$hm.cexRow, fontsize_col= global.plotparam$hm.cexCol, style=global.param$which.test, anno.col=anno.col, anno.col.color=anno.col.color, show.rownames=global.plotparam$hm.show.rownames, show.colnames=global.plotparam$hm.show.colnames)
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
              rmd <- paste(rmd, "\n
                         \n Scatterplot of ```r global.plotparam$pca.x``` and ```r global.plotparam$pca.y```.
                         \n```{r pca-scatter, echo=F}
                         \np <- plotlyPCA( reactiveValuesToList(global.param), reactiveValuesToList(global.results), reactiveValuesToList(global.input), input$pca.x, input$pca.y, grp.other=input$pca.grp.col)
                         \np$p 
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
              \n<a name="corrmat"></a>
              \n#### Correlation matrix
              \n<br>
              \n```{r corrmat, echo=F, warning=F, message=F, fig.width=8, fig.height=8}
              \nwithProgress(message="Exporting", detail="correlation matrix",{
              \nif(is.null(global.results$table.log)){
              \n  tab <- data.frame(global.results$table.na.filt)
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
                   \n## calculate correlation matrix
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
          
          ## ##########################################################
          ##                correlation boxplots
          ## ##########################################################
          if(export.cm){
            
            rmd <- paste(rmd, 
                         '\n***
              \n<a name="corrbox"></a>
              \n#### Correlation boxplots
              \nThe boxplots depict the distribution of correlation coefficients within a sample group.
              \n<br>
              \n```{r corrbox, echo=F, warning=F, message=F, fig.width=8, fig.height=8}
              \nwithProgress(message="Exporting", detail="correlation matrix",{
              \nif(is.null(global.results$table.log)){
              \n  tab <- data.frame(global.results$table.na.filt)
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
                   \n## calculate correlation matrix
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
                \nplotCorrBox(cm, grp, grp.col.leg, global.plotparam$cm.upper)
              \n}) # end with progress
              \n```\n
              \n 
              \n[Back to top](#top)'          
                         , sep='\n')
            
          } # end if export.cm
          
          
          
          
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
            stat_header <-'^logFC\\.|^AveExpr\\.|^t\\.|^P\\.Value\\.|^adj\\.P\\.Val\\.|^Log\\.P\\.Value\\.|^RawlogFC\\.|^RawAveExpr\\.'
            
            #do this for each column
            for(i in 1:length(colnames.tmp)){
              #if this is not a stat-related column, skip
              if(!grepl(stat_header,colnames.tmp[i])){
                next
              }else{
                #remove the stat header to get the group comparison
                colnames.split <- unlist(strsplit(colnames.tmp[i],stat_header))
                groups <- colnames.split[2]
                #flip the order
                g.new <- unlist(strsplit(groups, '\\.vs\\.'))
                g.new <- paste0(g.new[2], '.over.', g.new[1])
                #put the newly ordered groups back into the column
                header.loc <- str_locate(colnames.tmp[i],stat_header)
                my.header <- substring(colnames.tmp[i],header.loc[1],header.loc[2])
                colnames.tmp[i] <- paste0(my.header,g.new)
              }
            }
            
            
            # old code which just matched substrings was inaccurate when a comparison was an exact substring match
            # for(g in grp.comp){
            #   g.new <- strsplit(g, '\\.vs\\.') %>% unlist
            #   g.new <- paste(g.new[2], '.over.', g.new[1], sep='')
            #   colnames.tmp <- gsub(paste0(g,'$'), g.new, colnames.tmp)
            # }
            colnames(res.comb) <- colnames.tmp
          }
          
          ## filename 
          fn.tmp <- sub(' ','_',
                        paste(
                          global.param$label, '_',
                          sub(' ', '_',global.param$which.test),
                          ifelse(global.param$log.transform != 'none', paste( '_', global.param$log.transform, '_', sep=''), '_'),
                          ifelse(input$repro.filt=='yes', paste(global.param$filt.data, sep=''), '_'),
                          sub(' .*', '', Sys.time()), sep='') 
          )
          
          ## assemble gct file
          withProgress(message='Exporting', detail='GCT file',{
            rdesc <- res.comb
            mat <- rdesc[, names(grp.srt)] %>% data.matrix
            rdesc <- rdesc[ ,-which(colnames(rdesc) %in% names(grp.srt))]
            if(global.param$file.gct3){

              cdesc <- global.input$cdesc[ names(grp.srt),]
            } else {
              cdesc <- data.frame(experiment=grp.srt)
            }
            cat( paste0(rep('-', 20), '\nid-column: ', global.param$id.col.value, '\n', rep('-', 20), '\n') )
            cid <- names(grp.srt)
            
            
            #rid <- rdesc[, global.param$id.col.value] #rownames(mat)
            rid <- rdesc[, 'id'] #rownames(mat)
            
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
          
          if(!is.null(error$msg)) return()
          
          ## update label
          global.param$label <- gsub('_| |,|;|\\:|\\+|\\*', '-', input$label)
          
          #########################################################
          ## export
          withProgress(message='Exporting', detail='Excel sheet',{
            
            ## generate_filename
            global.param$xls.name <- create.fn(global.param$label, global.param$which.test, global.param$log.transform, 
                                               input$repro.filt, global.param$filt.data, suffix='xlsx')
              
            ## Excel
            render.xlsx <- export2xlsx(res.comb=global.results$data$output, grp=global.param$grp, 
                              grp.comp=unique(global.param$grp.comp), rdesc=global.input$table.anno, 
                          which.test=global.param$which.test,  fn=global.param$xls.name, 
                          session.dir=global.param$session.dir, headerDesc=headerDesc)
              
           }) ## end with.progress
     
          if(!(class(render.xlsx) == 'try-error' | !render.xlsx))
            global.results$export.xls <- TRUE
          
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
                  #tab <- data.frame(global.input$table)
                  tab <- data.frame(global.results$table.na.filt)
                else
                  tab <- data.frame(global.results$table.log)
                
                ## id column
                id.col.value <- global.param$id.col.value
                
                ## group vector
                grp <- global.param$grp
                grp <- sort(grp)
                ## group colors
                grp.col <- global.param$grp.colors
                grp.col.leg <- global.param$grp.colors.legend
                
                
                ## update selection
                tab <- tab[, c(id.col.value, names(grp))]
                grp.col <- grp.col[names(grp)]
                grp.col.leg <- grp.col.leg[unique(grp)]
                
                ## get correlation matrix
                if(is.null(global.results$cm) | global.param$update.cm == TRUE){
                  
                  
                  ## calculate correlation matrix
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
            #############################################################
            ##                correlation boxplot
            #############################################################
            if(input$export.cb){
              withProgress(message='Exporting', detail='correlation boxplot',{
                
                fn.cm <- paste0(global.param$session.dir, '/correlation_boxplot.pdf')
                
                if(is.null(global.results$table.log))
                  tab <- data.frame(global.results$table.na.filt)
                  #tab <- data.frame(global.input$table)
                else
                  tab <- data.frame(global.results$table.log)
                
                ## id column
                id.col.value <- global.param$id.col.value
                
                ## group vector
                grp <- global.param$grp
                ## group colors
                grp.col <- global.param$grp.colors
                grp.col.legend <- global.param$grp.colors.legend
                
                
                ## update selection
                tab <- tab[, c(id.col.value, names(grp))]
                grp.col <- grp.col[names(grp)]
                grp.col.legend <- grp.col.legend[unique(grp)]
                
                ## get correlation matrix
                if(is.null(global.results$cm) | global.param$update.cm == TRUE){
                  
                  ## calculate correlation matrix
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
                pdf(fn.cm,  max(2*(length(unique(grp))), 4), 6)
                plotCorrBox(cm, grp, grp.col.legend, global.plotparam$cm.upper)
                dev.off()
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
                    
                    ######################################
                    ## extract results
                    res = global.results$filtered
                    
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
                    if (!is.null(global.param$grp.selection)){
                      res = res[, names(global.param$grp.selection)]
                      grp.hm <- global.param$grp.selection
                    }else{
                      res = res[, names(global.param$grp)]
                      grp.hm <- global.param$grp
                    }
                    
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

                    pdf(fn.hm, height=12, width=24)
                    if(input$hm.max){
                        print(plotHM(res=res, hm.rownames=hm.rownames, grp=grp.hm, grp.col=global.param$grp.colors, grp.col.legend=global.param$grp.colors.legend,  hm.clust=input$hm.clust, hm.title=hm.title, hm.scale=input$hm.scale, cellwidth=cw, fontsize_row=input$cexRow, fontsize_col=input$cexCol, max.val=input$hm.max.val, style=global.param$which.test, anno.col=anno.col, anno.col.color=anno.col.color, show.rownames=input$hm.show.rownames, show.colnames=input$hm.show.colnames,
                               height=min( dynamicHeightHM( nrow(global.results$filtered)), 1200 ),
                               width=dynamicWidthHM(length(global.param$grp))))
                        
                    } else {
                        print(plotHM(res=res, hm.rownames=hm.rownames, grp=grp.hm, grp.col=global.param$grp.colors, grp.col.legend=global.param$grp.colors.legend,  hm.clust=input$hm.clust, hm.title=hm.title, hm.scale=input$hm.scale, cellwidth=cw, fontsize_row=input$cexRow, fontsize_col=input$cexCol, style=global.param$which.test, anno.col=anno.col, anno.col.color=anno.col.color, show.rownames=input$hm.show.rownames, show.colnames=input$hm.show.colnames,
                               height=min( dynamicHeightHM( nrow(global.results$filtered)), 1200 ),
                               width=dynamicWidthHM(length(global.param$grp))))
                    }
                    dev.off()
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

            # if(input$export.pca.loadings){
            #         withProgress(message='Exporting', detail='PCA loadings as Excel sheet',{
            #                 if(is.null(global.results$data) | is.na(global.param$filter.value) | is.null(global.results$pca)) return()
            #                 if(!is.null(error$msg)) return()
            #                 if(length(global.param$grp) < 3) return()
            # 
            #                 pca <- global.results$pca
            #                 pca.loadings <- as.data.frame(pca$loadings)
            # 
            #                 ## generate_filename
            #                 fn.tmp <- sub(' ','_', paste(global.param$session.dir, '/', 'pcaLoadings_', sub(' ', '_',global.param$which.test),  ifelse(global.param$log.transform != 'none', paste('_', global.param$log.transform, sep=''), '_'), ifelse(global.param$norm.data != 'none', paste('_', global.param$norm.data, sep=''), '_'), ifelse(input$repro.filt=='yes', paste('_reprofilt', sep=''), '_'), sub(' .*', '', Sys.time()),".xlsx", sep=''))
            #                 global.param$pcaloadings.ExcelFileName <- fn.tmp
            #                 ## Excel
            #                 WriteXLS("pca.loadings", ExcelFileName= fn.tmp, SheetNames= "pcaLoadings", row.names=T, BoldHeaderRow=T, AutoFilter=T)
            # 
            #                 cat("-- writing pca_loadings --")
            # 
            #         })
            # 
            # }
            # ############################################################
            ##                    boxplots
            ############################################################
            if(input$export.box){
              
                withProgress(message='Exporting', detail='box plot',{
                    fn.box <- paste(global.param$session.dir, 'boxplots_unnormalized.pdf', sep='/')

                    ###############################################
                    ## unnormalized ratios
                    if(is.null(global.results$table.log))
                        tab <- data.frame(global.input$table)
                        #tab <- data.frame(global.results$table.na.filt)
                    else
                        tab <- data.frame(global.results$table.log)

                    ## id column
                    id.col.value <- global.param$id.col.value
                    ## group vector
                    if(global.param$norm.per.group){
                      grp <- global.param$grp.norm
                    }else{
                      grp <- global.param$grp
                    }
                    ## group colors
                    if(global.param$norm.per.group){
                      grp.col <- global.param$grp.colors.norm
                      grp.col.leg <- global.param$grp.colors.legend.norm
                    }else{
                      grp.col <- global.param$grp.colors
                      grp.col.leg <- global.param$grp.colors.legend
                    }

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
                    if(global.param$norm.per.group){
                      grp <- global.param$grp.norm
                    }else{
                      grp <- global.param$grp
                    }
                    ## group colors
                    if(global.param$norm.per.group){
                      grp.col <- global.param$grp.colors.norm
                      grp.col.leg <- global.param$grp.colors.legend.norm
                    }else{
                      grp.col <- global.param$grp.colors
                      grp.col.leg <- global.param$grp.colors.legend
                    }

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
                        makeProfileplot(tab, id.col.value, grp, grp.col, grp.col.leg, 
                                        xlim.mode=input$profile.plot.xlim, 
                                        main='Before normalization')
                        makeProfileplot(tab.norm, id.col.value, grp, grp.col, grp.col.leg, 
                                        xlim.mode=input$profile.plot.xlim,
                                        main=paste(global.param$norm.data, 'normalized'))
                        dev.off()

                    } else{
                        pdf(fn.profile, 7, 7)
                        par(mfrow=c(1,1))
                        makeProfileplot(tab, id.col.value, grp, grp.col, grp.col.leg, 
                                        xlim.mode=input$profile.plot.xlim,
                                        main='unnormalized')
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
                fn.pval <- paste(global.param$session.dir, 'histogram_p-values.pdf', sep='/')

                ## unfiltered results
                res.all = global.results$data$output

                pdf(fn.pval, 10, 5*ifelse( global.param$which.test != 'mod F', length(grp.comp), 1 ))

                ############################################
                ## mod T
                if(!(global.param$which.test %in% c('mod F'))){
                    par(mfrow=c(length(grp.comp),1))
                    for(g in grp.comp){
                          pval <- res.all[, paste('P.Value', g, sep='.')]
                          hist(pval, breaks=50, main=paste('Histogram of p-values (N=', sum(!is.na(pval)), ')',sep=''), xlab='p-value', cex.main=2.2, cex.axis=2, cex.lab=2, col='darkblue', border=NA)
                          legend('top', legend=g, cex=2)
                    }
                 ############################################
                 ## mod F
                } else {
                    pval <- res.all[, paste('P.Value')]
                    hist(pval, breaks=50, main=paste('Histogram of p-values (N=', sum(!is.na(pval)), ')',sep=''), xlab='p-value', cex.main=2.2, cex.axis=2, cex.lab=2, col='darkblue', border=NA)
                }
                dev.off()
                })

            }

            #########################################################
            ##               Excel sheet
            #########################################################
            if(input$export.excel){
              
                withProgress(message='Exporting', detail='Excel sheet',{

                  
                  ## generate_filename
                  global.param$xls.name <- create.fn(global.param$label, global.param$which.test, global.param$log.transform, 
                                                     input$repro.filt, global.param$filt.data, suffix='xlsx')
                  
                  ## Excel
                  render.xlsx <- export2xlsx(res.comb=global.results$data$output, grp=global.param$grp, 
                                             grp.comp=unique(global.param$grp.comp), rdesc=global.input$table.anno, 
                                             which.test=global.param$which.test, fn=global.param$xls.name, 
                                             session.dir=global.param$session.dir, headerDesc=headerDesc)
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
                stat_header <-'^logFC\\.|^AveExpr\\.|^t\\.|^P\\.Value\\.|^adj\\.P\\.Val\\.|^Log\\.P\\.Value\\.|^RawlogFC\\.|^RawAveExpr\\.'
                
                #do this for each column
                for(i in 1:length(colnames.tmp)){
                  #if this is not a stat-related column, skip
                  if(!grepl(stat_header,colnames.tmp[i])){
                    next
                  }else{
                    #remove the stat header to get the group comparison
                    colnames.split <- unlist(strsplit(colnames.tmp[i],stat_header))
                    groups <- colnames.split[2]
                    #flip the order
                    g.new <- unlist(strsplit(groups, '\\.vs\\.'))
                    g.new <- paste0(g.new[2], '.over.', g.new[1])
                    #put the newly ordered groups back into the column
                    header.loc <- str_locate(colnames.tmp[i],stat_header)
                    my.header <- substring(colnames.tmp[i],header.loc[1],header.loc[2])
                    colnames.tmp[i] <- paste0(my.header,g.new)
                  }
                }
                
                # old code which just matched substrings was inaccurate when a comparison was an exact substring match
                # for(g in grp.comp){
                #   g.new <- strsplit(g, '\\.vs\\.') %>% unlist
                #   g.new <- paste(g.new[2], '.over.', g.new[1], sep='')
                #   colnames.tmp <- gsub(paste0(g,'$'), g.new, colnames.tmp)
                #  # colnames.tmp <- gsub(g, g.new, colnames.tmp)
                #   
                #  # grp.comp.pval.gct[g] <- g.new
                #   logp.colnames[g] <- paste0('Log.P.Value.', g.new)
                #   logfc.colnames[g] <- paste0('logFC.', g.new)
                # }
                logp.colnames <- colnames.tmp[grepl('^Log\\.P\\.Value',colnames.tmp)]
                logfc.colnames <- colnames.tmp[grepl('^logFC',colnames.tmp)]
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
              if( !global.param$which.test %in% c('mod F','none')){
                    withProgress(message='Exporting', detail='GCT file (transformed P-values)',{
                      rdesc <- res.comb
                      global.param.list <- reactiveValuesToList(global.param)
                     
                      #save(rdesc, logp.colnames, logfc.colnames, file='debug.RData')
                      logp <- rdesc[, logp.colnames] %>% data.matrix
                      fc <- rdesc[, logfc.colnames ] %>% data.matrix
                      
                      ## transformed and signed p-values
                      mat <- logp*sign(fc)
                      #colnames(mat) <- paste0('signed.', colnames(logp))
                      colnames(mat) <- paste0('signed.', logp.colnames) 
                      
                      rdesc <- rdesc[ ,-which(colnames(rdesc) %in% names(grp.srt))]
                      
                      cid <- colnames(mat)
                      rid <- rdesc[, 'id']
                      
                      res.gct <- new('GCT')
                      res.gct@mat <- mat
                      res.gct@rid <- rid
                      res.gct@cid <- cid
                      res.gct@rdesc <- rdesc
                      
                      write.gct(res.gct, ofile =  sub('\\.xlsx','-transformed-p-val', paste(global.param$session.dir, fn.tmp, sep='/')) )
                    })
              } ## end if not mod F
              
              
              #####################################################
              ## GCT file
              
              ## assemble gct file
              withProgress(message='Exporting', detail='GCT file',{
                    
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
                    #rid <- rdesc[, global.param$id.col.value] ## 20200928
                    rid <- rdesc[, 'id']
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
                fn.zip <- paste( gsub('\\:', '', gsub(' ','-', gsub('-','',Sys.time()))),'.zip', sep='')
            }
            ## label present
            if(!is.null(global.param$label) | nchar(global.param$label) == 0){
                fn.zip <- paste( global.param$label, '_', gsub('\\:', '', gsub(' ','-', gsub('-','',Sys.time()))),'.zip', sep='')
            }

            #########################################################
            ##   export parameters as small text file
            #########################################################
            global.param.imp <- reactiveValuesToList(global.param)
            params <- global.param.imp[c('log.transform', 'norm.data', 'na.filt.val', 'filt.data',  'repro.filt.val', 'sd.filt.val', 'which.test', 'filter.type', 'filter.value')]
            
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

            ####################################################
            ##             clean up
            ## remove archived files: all files except RData
	        rdata.idx <- grep('\\.RData[\\"|\']$', fn.all.abs)
	        fn.rm <- fn.all.abs
	        if(length(rdata.idx) > 0)
	           fn.rm <- fn.rm[-rdata.idx]
	          
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
        ##################################################################################
        observeEvent(input$ms.robustify,{
          updateCheckboxInput('ms.max', session = session, value=!input$ms.robustify)
        })
        observeEvent(input$ms.max,{
          updateCheckboxInput('ms.robustify', session=session, value=!input$ms.max)
        })
        
        
        
   
        ## ###############################################################
        ##
        ## import saved session from the server via drop down menu
        ##
        ## ###############################################################
        observeEvent(input$session.browse.import, {

            ## if nothing has been selected...
            if( input$session.browse == '' ) return()
          
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
            if(is.null(global.param$QC.filter))
              global.param$QC.filter <- F
          
            ## backwards compatibility III
            if( is.null(global.results$table.na.filt) )
                global.results$table.na.filt <- global.input$table

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

            ###################################################################
            ##            insert the panels for the gprofiler plots
            ###################################################################
            #ins.gprof()
            
            ## #####################################
            ## generate id.map for compatibility with < v0.7.0
            if(is.null(global.results$id.map)){

                ## ############################################
                ## map to gene names
                ## ############################################
                ##ids <- global.results$data$output[, global.param$id.col.val]
                ids <- global.results$data$output[, 'id']

                map.res <- mapIDs( ids )

                global.results$keytype <- map.res$keytype
                global.results$id.map <- data.frame(id=ids,  map.res$id.map )
                global.results$data$output <- left_join(global.results$data$output, global.results$id.map, 'id')
            }
        }) ## end: observeEvent(input$session.browse.import

        
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
            
            global.param$grp.comp.all <- groups.comp
            
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
          
            ## get all distinct pairwise comparisons  
            groups.comp <- get_pairwise_group_comparisons(groups.unique, mode='sort')
            
            # groups.comp <- c()
            # count=1
            # for(i in 1:(length(groups.unique)-1))
            #   for(j in (i+1):length(groups.unique)){
            #       
            #     ## order alphabetically
            #     groups.tmp <- sort( groups.unique[c(i,j)] )
            #     groups.comp[count] <- paste(groups.tmp[1], groups.tmp[2], sep='.vs.')
            #     count <- count+1
            #   }
          }
          global.param$grp.comp.all <- groups.comp
          
          if(is.null(global.param$grp.comp.selection) | sum( global.param$grp.comp.selection %in% global.param$grp.comp.all) == 0 )
            global.param$grp.comp.selection <- groups.comp
        })
        
        ################################################################################
        ##
        ##                once the 'run test' button was pressed...
        ##
        ## - log-transform (optional)
        ## - normalization (optional)
        ## - filter missing values (optional)
        ## - filter data(optional)
        ## - test (optional)
        ## - flag for intensity data
        ##
        ################################################################################
        observeEvent(input$run.test, {

            ## reset any error messages
            error$msg <- NULL

            global.input$run.test <- input$run.test


            ###########################################
            ## - store use selections
            global.param$which.test <- input$which.test
            global.param$na.filt.val <- input$na.filt.val
            global.param$filt.data <- input$filt.data
            global.param$repro.filt.val <- input$repro.filt.val
            if(!is.null(input$sd.filt.val)) global.param$sd.filt.val <- input$sd.filt.val
            global.param$norm.data <- input$norm.data
            #global.param$norm.per.group <- input$norm.per.group
            global.param$log.transform <- input$log.transform
            global.param$intensity <-input$intensity

            #####################################################################
            ## if the 'Run test' - button has been pressed for the first time,
            ## store a copy of the original input matrix (un-normalized, unfiltered)
            if(!(global.param$analysis.run)){
                global.input$table.org <- global.input$table
                tab <- data.frame(global.input$table)
            }
            if(global.param$analysis.run){
                tab <- data.frame(global.input$table.org)
            }
            
            ## id column
            #id.col = global.param$id.col.value
            id.col <- 'id'
            
            ## all group labels
            groups <- global.param$grp.selection
            global.param$grp <- groups
            
            ###############################################
            ## initialize values for normalized and filtered matrices
            global.results$table.log=NULL
            global.results$table.norm=NULL
            global.results$table.repro.filt=NULL
            global.results$table.sd.filt=NULL
            global.results$table.na.filt=NULL
            global.results$values.filtered=NULL
            global.results$pca=NULL

            ## ##############################################################################
            ##
            ##           set up the analysis pipeline
            ##
            ## ##############################################################################
            test <- global.param$which.test
            norm.data <- global.param$norm.data
            intensity <- global.param$intensity
            
            ## NA filter: backwards compatibility
            if(!is.null(global.param$na.filt.val)){
                na.filt.val <- as.numeric(global.param$na.filt.val)
            } else {
                na.filt.val <- 100
            }
            
            filt.data <- global.param$filt.data
            repro.filt.val <- global.param$repro.filt.val
            sd.filt.val <- global.param$sd.filt.val
            log.trans <- global.param$log.transform

            # update group selection
            groups.comp <- global.param$grp.comp.selection
            global.param$grp.comp <- groups.comp
            
            ## ###########################################
            ## log transformation
            ## ###########################################
            if(log.trans != 'none'){
                
                ids.tmp <- tab[, id.col]
                dat.tmp <- tab[, -which(colnames(tab) == id.col)]
                dat.tmp <- data.matrix(dat.tmp)
                dat.tmp[dat.tmp == 0] <- NA
                
                #if there are negative values in the matrix, do not log transform!
                if(sum(na.omit(dat.tmp<0))>0){
                  shinyalert("Dataset contains negative values!", "Analysis will proceed WITHOUT log-transformation. If you wish to log-transform, please re-upload a dataset without negative values.", type = "warning")
                  log.trans="none"
                  global.param$log.transform = log.trans
                }else if(log.trans == 'log2'){
                    dat.tmp <- log(dat.tmp, 2)
                }else if(log.trans == 'log10'){
                    dat.tmp <- log(dat.tmp, 10)
                }
                ## putting a data frame around here turns out to be ESSENTIAL!! DON'T USE CBIND HERE!
                tab <- data.frame( ids.tmp, dat.tmp)
                colnames(tab)[1] <- id.col
                
                ## store
                global.results$table.log <- tab
            }

            ## ###########################################
            ## calculate FC between groups before 
            ## normalization
            ## ###########################################
            fc.before.norm <- calculate_fc(tab,global.param$grp, groups.comp, test)
            #add rownames
            rownames(fc.before.norm) <- rownames(tab)
            
            
            ## ###########################################
            ##
            ## normalization
            ##
            ## ###########################################
            if(norm.data != 'none'){
               
                withProgress(message='Applying normalization...', {
                    ## run normalization
                  if(!global.param$norm.per.group){
                    tab.norm <- normalize.data(tab, id.col, norm.data )
                  }else{
                    if(is.null(global.param$grp.norm) | input$grp.norm=="None"){
                      ## throw error
                      shinyalert("Error applying group-wise normalization!", "No Group column detected. Proceeding with global normalization. Press OK to close this window.", type = "warning")
                      global.param$norm.per.group=FALSE
                      tab.norm <- normalize.data(tab, id.col, norm.data )
                    }else{
                      tab.norm <- normalize.data(tab, id.col, norm.data, grp.vec=global.param$grp.norm )
                    }
                  }
                    
                    ## if two-component norm fails....
                    if(class(tab.norm) == 'try-error'){
                        
                      ## reset test
                      test='none'
                      global.param$which.test <- 'none'
                      updateRadioButtons('which.test', session=session, selected = 'none')
                      
                      ## reset normalization
                      norm.data='none'
                      global.param$norm.data <- 'none'
                      updateRadioButtons('norm.data', session=session, selected = 'none')
                      global.results$table.norm <- NULL
                      
                      ## reset significance filter
                      global.param$filter.type <- 'none'
                      updateSelectInput('filter.type', session=session, selected = 'none')
                      
                      ## throw error
                      error$title <- "'Error applying two-component normalization."
                      error$msg <- HTML('The two-component normalization failed to converge on at least one data column. Please note that this type of normalization expects <b>log-ratio</b> data that is approximately <b>centered around zero</b>. Please make sure this is the case by <b>inspecting the profile plots</b> under the QC tab.')
                      
                    } else {
                        ## store normalized table    
                        global.results$table.norm <- tab.norm
                    }
                    
                    ## Make sure the normalized table is used
                    ## in downstream analysis.
                    ## This was broken in v0.8.8 and v0.8.8.1 
                    tab <- tab.norm
                })
            }
            ## ############################################
            ##
            ## missing value filter
            ##
            ## ############################################
            if(na.filt.val < 100) {
                
                withProgress(message='Applying missing data filter',  {
                    na.filt = na.filter(tab, id.col=id.col, groups, na.filt.val=na.filt.val)
                    tab.na.filt <- na.filt$table
                })
                
                ## store
                global.results$table.na.filt <- tab.na.filt
                global.results$ids.na.filt <- na.filt$ids
                
                ## update other tables 
                global.results$table.norm <- global.results$table.norm[na.filt$ids, ]
                global.results$table.log <- global.results$table.log[na.filt$ids, ]
                
                ## fc.before.norm doesn't have rownames..
                fc.before.norm  <- fc.before.norm[ which(tab[, id.col] %in% na.filt$ids), ,drop=F]
                
                ## use na-filtered table downstream
                tab <- tab.na.filt
                
            } else {
                global.results$table.na.filt <- tab
                global.results$ids.na.filt <- tab$id
            }
            
            ## ############################################
            ##
            ## reproducibility filter
            ##
            ## ############################################
            if(filt.data == 'Reproducibility'){

                    withProgress(message='Applying reproducibility filter',  {
                        repro = my.reproducibility.filter(tab, id.col=id.col, groups, alpha=global.param$repro.filt.val)
                        tab = repro$table
                        
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
                    })
                
                    ## store indices of filtered values in the original table
                    global.results$values.filtered <- filt.tmp$values.filtered
                    global.results$table.filt <- tab
                    global.param$sd.perc.val <- filt.tmp$sd.perc.val
                    global.param$sd.filt.str <- paste( global.param$filt.data, ' (SD > ', round(global.param$sd.perc.val, 2),' / ', global.param$sd.filt.val, 'th percentile)',sep='')
            }

            ## #############################################################################
            ##
            ##                                     TEST
            ##
            ###############################################################################

            ##################################
            ## two sample
            if(test == 'Two-sample mod T'){
                
                withProgress(message='Two-sample test', value=0, {

                    count=0
                    res.comb <- tab
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
                        res.tmp <-  modT.test.2class( tab.group, groups=groups.tmp, id.col=id.col, label=g , intensity=intensity)$output
                        #previous code would incorrectly throw away a test result if a feature was missing
                        # if(count == 0){
                        #     res.comb <- res.tmp
                        # } else {
                        #     ## make sure the order is correct
                        #     if(nrow(res.tmp ) != nrow(res.comb)) stop( "number of rows don't match!\n" )
                        #     #res.tmp <- res.tmp[rownames(res.comb), ]
                        #     ##res.comb <- cbind(res.comb, res.tmp)
                        #     res.comb <- data.frame(res.comb, res.tmp, stringsAsFactors=F)
                        # }
                        
                        
                        #create data frame of expression values and test results
                        res.test <- res.tmp[, !colnames(res.tmp)%in%colnames(res.comb)]
                        res.comb <- base::merge(res.comb,res.test,by="row.names",all=T)
                        rownames(res.comb) <- res.comb[,1]
                        res.comb <- res.comb[,-1]
                        
                        ##################################################
                        ## progress bar
                        incProgress(count/length(unique(groups.comp)), detail=g)
                        count=count + 1

                    }
                
                ##################################    
                ## add FC before normalization
                res.comb <- merge(res.comb, fc.before.norm,by="row.names")
                rownames(res.comb) <- res.comb[,1]
                res.comb <- res.comb[,-1]
                    

                ##################################
                ## reorder columns in table
                res.id <- res.comb$id ## id column
                res.exprs <- res.comb[, names(groups)] ## expression values
                res.test <- res.comb[, grep('^logFC\\.|^AveExpr\\.|^t\\.|^P\\.Value\\.|^adj\\.P\\.Val\\.|^Log\\.P\\.Value\\.|^RawlogFC\\.|^RawAveExpr\\.', colnames(res.comb))] ## test results
                res.test <- res.test[, order(colnames(res.test))]

                ## assemble new table
                res.comb <- data.frame(id=res.id, res.test, res.exprs)
                res.comb <- res.comb[rownames(tab),]

                })
            }
            ##################################
            ## one sample
            if(test == 'One-sample mod T'){
               
                withProgress(message='One-sample T test', value=0, {

                    count=0
                    res.comb <- tab
                    ## loop over groups
                    for(g in unique(groups.comp)){

                        ## extract table of current group
                        tab.group <- cbind(tab[, id.col], tab[, names(groups)[which(groups == g)]])
                        colnames(tab.group)[1] <- id.col

                        res.tmp <- modT.test( tab.group, id.col=id.col, plot=F, nastrings=NASTRINGS, label=g, na.rm=FALSE)$output
                        #previous code would incorrectly throw away a test result if a feature was missing
                        # if(count == 0){
                        #     res.comb <- res.tmp
                        # } else {
                        #     if(nrow(res.tmp ) != nrow(res.comb)) stop( "number of rows don't match!\n" )
                        #     res.tmp <- res.tmp[rownames(res.comb), ]
                        #     res.comb <- cbind(res.comb, res.tmp)
                        # }
                        
                        #create data frame of expression values and test results
                        res.test <- res.tmp[, !colnames(res.tmp)%in%colnames(res.comb)]
                        res.comb <- base::merge(res.comb,res.test,by="row.names",all=T)
                        rownames(res.comb) <- res.comb[,1]
                        res.comb <- res.comb[,-1]

                        #############################################
                        ## update progress bar
                        incProgress( 1/length(unique(groups.comp) ), detail=g)
                        count=count + 1
                    }
                })

                ##################################    
                ## add FC before normalization
              res.comb <- merge(res.comb, fc.before.norm,by="row.names")
              rownames(res.comb) <- res.comb[,1]
              res.comb <- res.comb[,-1]
                
                ##################################
                ## reorder columns of the table
                res.id <- res.comb$id ## id column
                res.exprs <- res.comb[, names(groups)] ## expression values
                res.test <- res.comb[, grep('^logFC\\.|^AveExpr\\.|^t\\.|^P\\.Value\\.|^adj\\.P\\.Val\\.|^Log\\.P\\.Value\\.|^RawlogFC\\.|^RawAveExpr\\.', colnames(res.comb))] ## test results
                res.test <- res.test[, order(colnames(res.test))]
                ## assemble new table
                res.comb <- data.frame(id=res.id, res.test, res.exprs, stringsAsFactors=F)
                ##res.comb <- res.comb[tab[, id.col], ]
                res.comb <- res.comb[rownames(tab),]
                
            }

            ##################################
            ## moderated F test
            if(test == 'mod F'){
                
                withProgress(message='moderated F-test', value=0, {
                
                    tab.group <- cbind(tab[, id.col], tab[, names(groups)])
                    colnames(tab.group)[1] <- id.col
                  
                    res.comb <- modF.test( tab.group, id.col=id.col, class.vector=groups, nastrings=NASTRINGS, na.rm=FALSE, intensity=intensity)$output
                    colnames(res.comb) <- sub('^X', '', colnames(res.comb))
                })
                
                ##################################    
                ## add FC before normalization
              res.comb <- merge(res.comb, fc.before.norm,by="row.names")
              rownames(res.comb) <- res.comb[,1]
              res.comb <- res.comb[,-1]
                
                ##################################
                ## reorder columns of the data table
                res.id <- res.comb$id ## id column
                res.exprs <- tab[, names(groups)]
                res.test <- res.comb[, grep('^logFC\\.|^AveExpr|^t|^P\\.Value|^adj\\.P\\.Val|^Log\\.P\\.Value|^RawlogFC\\.|^RawAveExpr\\.', colnames(res.comb))] ## test results
                res.test <- res.test[, order(colnames(res.test))]

                ## assemble new table
                res.comb <- data.frame(id=res.id, res.test, res.exprs, stringsAsFactors=F)
                ##res.comb <- res.comb[tab[, id.col], ]
                res.comb <- res.comb[rownames(tab),]
           }
            
            ##########################################################
            ## two group linear model with repeated measurements
            ## -- NOT FUNCTIONAL --
            if(test == 'Two-sample LM'){
                
                withProgress(message='Two-sample LM', value=0, {
                    
                    count=0
                    ## loop over groups
                    for(g in unique(groups.comp)){
                        
                        ## extract current groups
                        groups.tmp <- groups[ c(grep( paste('^',unlist( strsplit(g, '\\.vs\\.'))[1],'$', sep='') , groups), grep(paste('^',unlist( strsplit(g, '\\.vs\\.'))[2], '$', sep='') , groups)) ]
                        
                        ## extract table of current group
                        tab.group <- cbind(tab[, id.col], tab[, names(groups.tmp)])
                        colnames(tab.group)[1] <- id.col
                        
            
                        groups.tmp <- groups.tmp %>% as.factor
                        reps.tmp <- cdesc[, repeats] %>% as.factor %>% as.numeric
                        
                        #############################
                        ## the actual test
                        #############################
                        
                        ## ###############################################
                        ## design matrix
                        design <- model.matrix(~ 0 + groups.tmp + reps.tmp)
                        colnames(design) <- make.names(colnames(design))
                        
                        
                        ## ###############################################
                        ## contrast matrix
                        contrasts <- makeContrasts( contrasts=glue('group{group.levels[1]}-group{group.levels[2]}') , levels=design)
                        
                        ## ################################################
                        ## run the model
                        res.tmp <- lm.repeats(tab.group, design, contrasts, repeats='rep')
                        
                        #res.tmp <-  modT.test.2class( tab.group, groups=groups.tmp, id.col=id.col, label=g )$output
                        
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
                    
                })
                
            }
            

            ###################################################################
            ##
            ##                       no test
            ##
            ###################################################################
            if(test == 'none'){

                ###########################################
                ## store data matrix as test results
                ## - values are log, if performed
                #tab <- data.frame(id=tab[ , id.col], tab)
                res.comb <- tab
                
                ##################################    
                ## add FC before normalization
                res.comb <- merge(res.comb, fc.before.norm,by="row.names")
                rownames(res.comb) <- res.comb[,1]
                res.comb <- res.comb[,-1]

            }

            ############################################
            ## add logFC.raw columns
            #res.comb <- data.frame(res.comb, fc.before.norm)
            
            ## #########################################
            ##   add id-mapping if not present already
            ##
            ## #########################################
            if( sum(c('id.concat', 'id.mapped', 'id.query') %in% colnames(res.comb)) < 3 ){
                res.comb <- left_join(res.comb, global.results$id.map, 'id')
            }
            rownames(res.comb) <- make.names(res.comb$id, unique = T)

            ## #########################################
            ## store the results
            global.results$data$output <- res.comb
            
            ## number of features with at least one non-missing value
            global.results$N.feat <- sum(apply( tab[, names(groups)], 1, function(x) sum(is.na(x)/length(x))) < 1)

            ## #####################################
            ## set some flags
            global.param$analysis.run <- T
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

            ###################################################################
            ##            insert panels fo grprofiler
            #ins.gprof()

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
        ##                   filter the test results across multiple groups
        ##
        ##
        ###################################################################################
        observeEvent(c(input$filter.type, input$filter.value.top.n, input$filter.value.nom.p, input$filter.value.adj.p, global.param$run.test),{

            ## only run after one analysis has been completed
            if(global.param$run.test == 0 | is.null(input$filter.type)) return()

            ## ######################################
            ## get the current tab
            tab.select <- input$mainPage

            groups.comp <- unique(global.param$grp.comp)
            groups <- global.param$grp

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

                      res <- res[ which( unlist(apply( res[, grep('^adj.P.Val', colnames(res) )], 1, 
                                                       function(x) sum( x < input$filter.value.adj.p, na.rm=T ))
                                                ) > 0), ]
                      
                        ## now separate for each group comparison
                        res.groups <- lapply(groups.comp, function(g){
                            res[ which( res[, paste('adj.P.Val.', g, sep='')] < input$filter.value.adj.p) , ]
                        })
                        names(res.groups) <- groups.comp
                    ###########################################
                    ## F-test
                    } else {
                        #View(res)
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
            global.results$filtered <- res
            global.results$filtered.groups <- res.groups
            
            #save(res, groups, res.groups, file='debug.RData')
            
            ## number of features after filtering
            #N.feat.filt <- nrow(res)
            N.feat.filt <- sum( apply(res[, names(groups)], 1, function(x) sum(is.na(x))/length(x) ) < 1, na.rm=T)
            global.results$N.feat.filtered <- N.feat.filt
          
            #save(res, groups, res.groups, groups.comp, file='debug.RData')

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
            
            ## flag to update gprofiler
            global.param$update.gprof.up <- TRUE
            global.param$update.gprof.dn <- TRUE
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
                exp.design <- cbind(colnames(tab), rep('NA', ncol(tab)), rep('NA', ncol(tab)))
                colnames(exp.design) <- c('Column.Name', 'Experiment', 'Group')
                write.table(  exp.design, file, sep='\t', quote=F, na='NA', row.names=F  )
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

            ## original input table
            tab <- data.frame(global.input$table.org)
            
            ## table after NA filter
            if(!is.null(global.results$table.na.filt)){
                tab.na.filt <- data.frame(global.results$table.na.filt)
                ## No. features after NA filter
                N.feat.na.filt <- nrow(tab.na.filt)
                global.results$N.feat.na.filt <- N.feat.na.filt
            }
            
            grp <- global.param$grp
            N.grp <- global.param$N.grp

            ## No. features with at least 1 non-missing value 
            N.feat <- sum(apply(tab, 1, function(x) sum(is.na(x)/length(x))) < 1)
            global.results$N.feat <- N.feat
            
            sum.tab <- t(data.frame(
                ifelse(is.null(global.results$N.feat), nrow(tab), N.feat),
                ifelse(is.null(global.results$N.feat.na.filt), NA, N.feat.na.filt),
                N.columns=length(grp),
                N.groups=N.grp
            ))
            
            sum.tab <- data.frame(id=c('No. features', 'No. features (NA-filtered)', 'No. expression columns', 'No. groups'), sum.tab)
            colnames(sum.tab) <- c('', 'Number')
            
            ## Number of features without quant in any sample
            if(!is.null(global.input$NA.rows)){
              row.tmp <-  c('No. features w/o quant', global.input$NA.rows)
              names(row.tmp) <- colnames(sum.tab)
              sum.tab <- rbind(sum.tab[1:2,], row.tmp, sum.tab[3:nrow(sum.tab),])
            }
            
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

            wf.tab.ids <- c('Intensity data', 'Group-wise normalization', 'Filter QC.fail', 'Log scale', 'Normalization', 'Max % missing (Min % quantified)', 'Filter data', 'Test', 'Filter results')

            ## #########################
            ## Reproducibility filter
            if(global.param$filt.data == 'Reproducibility'){

                wf.tab <- t(data.frame( global.param$intensity, global.param$norm.per.group, global.param$QC.filter, global.param$log.transform, global.param$norm.data, paste(global.param$na.filt.val," (",100-global.param$na.filt.val,")",sep=""), paste( global.param$filt.data, ' (alpha=',global.param$repro.filt.val, ')',sep=''), global.param$which.test,
                                       paste( global.param$filter.type, ' < ', global.param$filter.value)))

                wf.tab <- data.frame(id=wf.tab.ids, wf.tab, stringsAsFactors=F)

            }

            ## #######################
            ## SD filter
            if(global.param$filt.data == 'StdDev'){
                wf.tab <- t(data.frame( global.param$intensity, global.param$norm.per.group, global.param$QC.filter,global.param$log.transform, global.param$norm.data, paste(global.param$na.filt.val," (",100-global.param$na.filt.val,")",sep=""), global.param$sd.filt.str, global.param$which.test,
                                       paste( global.param$filter.type, ' < ', global.param$filter.value) ))
                wf.tab <- data.frame(id=wf.tab.ids, wf.tab, stringsAsFactors=F)
            }

            ## #############################
            ## no data filter
            if(global.param$filt.data == 'none'){
                wf.tab <- t(data.frame( global.param$intensity, global.param$norm.per.group, global.param$QC.filter,global.param$log.transform, global.param$norm.data, paste(global.param$na.filt.val," (",100-global.param$na.filt.val,")",sep=""), paste( global.param$filt.data, sep=''), global.param$which.test,
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

            res <- global.results$data$output

            ## tested groups
            grp.comp=unique(global.param$grp.comp)

            ## extract filter type
            filter.type=global.param$filter.type
            filter.value=global.param$filter.value

            ######################################
            ## one/two sample mod T
            #if(global.param$which.test != 'mod F'){
            if( !global.param$which.test %in% c('mod F', 'none')){
                
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
                nsets = length(grp.comp),
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
                nsets = length(grp.comp),
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
            
            if(global.param$log.transform == 'none')
                tab <- data.frame(global.results$table.na.filt)
                #tab <- data.frame(global.input$table.org)
            else
                tab <- data.frame(global.results$table.log)
                #tab <- data.frame(global.results$table.log)

            ## extract expression values
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
                tab <- data.frame(global.results$table.na.filt)
                #tab <- data.frame(global.input$table.org)
            else
                tab <- data.frame(global.results$table.log)
                #tab <- data.frame(global.results$table.log)

            grp <- global.param$grp
            N.grp <- global.param$N.grp
            grp.colors.legend <- global.param$grp.colors.legend
            grp.colors <- global.param$grp.colors

            ## extract expression values
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
            
            ## get column data types
            ## used for formatSignif() below
            column_type <- sapply(tab, class)
            numeric_column_idx <- which(column_type == 'numeric')
            
            dt <- datatable(tab, style = 'bootstrap', 
                           width = 1000, escape = F,  filter='top', 
                           rownames=F, 
                           options = list( pageLength = 10, scrollX = T, selection='none',
                                           autoWidth = TRUE#,
                                           #fnRowCallback = JS("function( nRow, aData, iDisplayIndex, iDisplayIndexFull ) {ind = 2; $('td:eq('+ind+')', nRow).html( (aData[ind]).toFixed(2) );}")
                           )
            )
            formatSignif( dt, digits=3, columns = numeric_column_idx )
 
        }, server=T)

        ## ################################################################################
        ##
        ##                             Volcano plot
        ##
        ## ################################################################################


        ## ##################################################################
        ## observer to trigger  'updateSelectizeInput' for volcano plots
        ## - server side  rendering of selectizeInput
        ## ##################################################################
        observe({

            if( !global.param$update.ppi.select ) return()

            grp.comp <- unique( global.param$grp.comp )

            for(i in 1:length(grp.comp)){

                ## commenting out the if-statemnet made select input work after multiple rounds of analysis
                ##if(!is.null(input[[ gsub('\\.','', paste0('ppi.bait.',  grp.comp[i])) ]] )){

                ## #################################
                ## all ids found in data
                choices=c( unlist(global.results$id.map$id.concat) )


                ## server-side rendering of selectizeInput
                updateSelectizeInput( session=session, inputId = gsub('\\.','',paste0('ppi.bait.',  grp.comp[i])), choices=choices, selected="", server=T)#,
            
            }
            global.param$update.ppi.select <- FALSE
        })
        
        ## ##################################################################
        ## obserever to trigger  'updateSelectizeInput' for SCATTER plots
        ## - server side rendering of selectizeInput
        ## ##################################################################
        observe({
            if( !global.param$update.ppi.select.scat ) return()

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
        ##            insert gprofile plots
        ##
        ## ##############################################################
        ins.gprof <- reactive({
          cat('gprof test...\n')
          withProgress({
            
            setProgress( message='Preparing gprofile plots...')
            
            grp.comp <- unique( global.param$grp )
            
            ## #########################################
            ## loop over group comparsions
            for(i in 1:length(grp.comp)){
              local({
                
                my_i <- i
                
                ## ########################
                ## the actual plots
                
                ## up regulated
                output[[ paste("gprof.up", grp.comp[my_i], sep='.') ]] <- renderPlotly({
                 ## cat('\ntes1...\n')
                   ## run gprofiler
              #if(is.null( global.results$gprof[[ paste("gprof.up", grp.comp[my_i], sep='.') ]]) | global.param$update.gprof.up ){
               ##cat('tesssttt:',global.param$update.gprof.up, '\n')
            #  if(global.param$update.gprof.up ){
                     cat('\ntes2...\n')         
                      gp <- run_gProfileR(reactiveValuesToList(global.results),
                                 reactiveValuesToList(global.param),
                                 direction='up',
                                 group=grp.comp[my_i],
                                 src=input[[paste0('gprof.src.', grp.comp[my_i])]],
                                 fdr=input[[paste0('gprof.fdr.', grp.comp[my_i])]])
                   
                      #if(class(gp) != 'try-error'){
                        global.results$gprof[[paste("gprof.up", grp.comp[my_i], sep='.')]] <- gp
                      #}
                     global.param$update.gprof.up <- FALSE
                      gostplot(gp$gp, interactive = T ) %>% layout(title=paste0('up-regulated (n=', length(gp$genes),')'))
                      
                   #} # else {
                   #  gp <- global.results$gprof[[ paste("gprof.up", grp.comp[my_i], sep='.') ]]  
                    ##cat('\ntest3...\n')
                    #save(gp, file='debug.RData')
                   #}
                   #validate(need(!is.null(gp), message = 'gp test'))
                   ## plot
                   #gostplot(gp$gp, interactive = T ) %>% layout(title=paste0('up-regulated (n=', length(gp$genes),')'))
                
                 })
                ## down regulated
                output[[ paste("gprof.dn", grp.comp[my_i], sep='.') ]] <- renderPlotly({
                  # gp <- run_gProfileR(reactiveValuesToList(global.results),
                  #                     reactiveValuesToList(global.param),
                  #                     direction='down',
                  #                     group=grp.comp[my_i],
                  #                     src=input[[paste0('gprof.src.', grp.comp[my_i])]],
                  #                     fdr=input[[paste0('gprof.fdr.', grp.comp[my_i])]])
                  # 
                  # validate(need(!is.null(gp), message = 'gp test'))
                  # 
                  # #gostplot(gp$gp, interactive = T ) %>% layout(title=paste0('down-regulated (n=', length(gp$genes),')'))
                  
                     })
              })
            }})
          
        }) ## end reactive
        
        


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
                      volc.label <- input[[paste('volc.label', grp.comp[my_i], sep='.')]]
                      if( volc.label=="ID_Symbol"&!is.null(global.results$id.map )){
                        txt.col <- 'id.concat' 
                      }else if(volc.label=="Symbol"&!is.null(global.results$id.map )){
                        txt.col <- 'id.mapped'
                      }else{
                        txt.col <- 'id'
                      }

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
                               volc[[paste('text', grp.comp[my_i], sep='.')]] = text.tmp[txt.col]
                               volc[[paste('id', grp.comp[my_i], sep='.')]] = text.tmp[ 'id' ]
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
                                      volc[[paste('id', grp.comp[my_i], sep='.')]] <- volc[[paste('id', grp.comp[my_i], sep='.')]][-idx]
                                      volc[[paste('xy', grp.comp[my_i], sep='.')]] <- volc[[paste('xy', grp.comp[my_i], sep='.')]][-idx]
                                      volc[[paste('P.Value', grp.comp[my_i], sep='.')]] <- volc[[paste('P.Value', grp.comp[my_i], sep='.')]][-idx]
                                      volc[[paste('adj.P.Val', grp.comp[my_i], sep='.')]] <- volc[[paste('adj.P.Val', grp.comp[my_i], sep='.')]][-idx]
                                  } else {
                                      volc[[paste('text', grp.comp[my_i], sep='.')]] <-volc[[paste('id', grp.comp[my_i], sep='.')]] <- volc[[ paste('y', grp.comp[my_i], sep='.') ]] <- volc[[paste('x', grp.comp[my_i], sep='.')]] <- volc[[paste('xy', grp.comp[my_i], sep='.')]] <- volc[[paste('adj.P.Val', grp.comp[my_i], sep='.')]] <- volc[[paste('P.Value', grp.comp[my_i], sep='.')]]<- NULL
                                                    }
                                  ################################################
                                  ## ADD: if selected point is not present add it to the list
                              } else{
                                  volc[[paste('x', grp.comp[my_i], sep='.')]]=c( volc[[paste('x', grp.comp[my_i], sep='.')]],
                                                                      text.tmp[paste('logFC', grp.comp[my_i], sep='.')])
                                  volc[[paste('y', grp.comp[my_i], sep='.')]]=c(volc[[paste('y', grp.comp[my_i], sep='.')]], text.tmp[paste('Log.P.Value', grp.comp[my_i], sep='.')])
                                  volc[[paste('text', grp.comp[my_i], sep='.')]]=c(volc[[paste('text', grp.comp[my_i], sep='.')]],  text.tmp[ txt.col] )
                                  volc[[paste('id', grp.comp[my_i], sep='.')]]=c(volc[[paste('id', grp.comp[my_i], sep='.')]],  text.tmp[ 'id' ] )
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
                      volc[[paste('id', grp.comp[my_i], sep='.')]] <- NULL
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
                      volc[[paste('id', grp.comp[my_i], sep='.')]] <- volc[[paste('id', grp.comp[my_i], sep='.')]][-idx]
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

                      dat.select = data.frame(id=unlist(volc[[paste('id', grp.comp[my_i], sep='.')]]), label=unlist(volc[[paste('text', grp.comp[my_i], sep='.')]]), logFC=round( unlist(volc[[paste('x', grp.comp[my_i], sep='.')]]), 2), P.Value=round( unlist(volc[[paste('P.Value', grp.comp[my_i], sep='.')]]),3), adj.P.Value=round( unlist(volc[[paste('adj.P.Val', grp.comp[my_i], sep='.')]]), 3) )
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
            x.ax <- res[ , input[[ paste0('scat.x.', group) ]] ]
            y.ax <- res[ , input[[ paste0('scat.y.', group) ]] ]
     
            cor.xy.pears <- cor(x.ax, y.ax, use='pairwise.complete', method = 'pearson')
            cor.xy.spear <- cor(x.ax, y.ax, use='pairwise.complete', method = 'spearman')
            
            
            ## ##############################
            ## ids
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
                  p <- plot_ly(x=dat.plot$x.ax[ non.int.idx ], y=dat.plot$y.ax[ non.int.idx], 
                               type='scatter', mode='markers', 
                               marker=list(size=10, color=col[non.int.idx]), 
                               text=IDs[non.int.idx], opacity=.2, showlegend=T, name=non.sig.txt )
                  
                  ## add significant proteins
                  if(length(sig.idx) > 0){
                    p <- p %>% 
                        add_trace(x=dat.plot$x.ax[ sig.idx ], y=dat.plot$y.ax[ sig.idx], 
                                  type='scatter', mode='markers', 
                                  marker=list(size=10, color=col[sig.idx]), 
                                  text=IDs[sig.idx], opacity=.2, showlegend=T, name=sig.txt  )
                  }
                  
                  ## #######################################
                  ## if interactors have been found
                } else {
                    col[ppi.int.idx] <- ppi.col[ nchar(ppi.col) > 0 ]
  
                    ## ###################################################
                    ## plot non-interactors
                    non.int.idx <- setdiff( 1:nrow(dat.plot), ppi.int.idx)
                 
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
              ## different traces for significant and not-significant features
             non.sig.idx <- 1:nrow(dat.plot)
              if(length(sig.idx) > 0)
                non.sig.idx <- non.sig.idx[ -c(sig.idx) ]  
              
              p <- plot_ly(x=dat.plot$x.ax[non.sig.idx], 
                           y=dat.plot$y.ax[non.sig.idx], 
                           type='scatter', mode='markers', 
                           marker=list(size=10, color=col[non.sig.idx]), 
                           text=IDs[non.sig.idx], showlegend=T, name=non.sig.txt )
            
              if(length(sig.idx) > 0){
                p <- p %>% add_trace(x=dat.plot$x.ax[ sig.idx ], 
                                     y=dat.plot$y.ax[ sig.idx], 
                                     type='scatter', mode='markers', 
                                     marker=list(size=10, color=col[sig.idx]), 
                                     text=IDs[sig.idx], opacity=1, showlegend=T, name='signif.'  )
              }
              
            }
            
             
            ## ##############################
            ##  reproducibility filter?
            values.filtered <- global.results$values.filtered[[group]]  
            
            if(length(values.filtered) > 0){

              ## get the entire dataset
              if(is.null(global.results$table.log))
                tab <- data.frame(global.results$table.na.filt)
                #tab <- data.frame(global.input$table)
              else
                tab <- data.frame(global.results$table.log)
              
              ## extract filtered values
              filt.x <- tab[values.filtered, input[[ paste0('scat.x.', group) ]]]
              filt.y <- tab[values.filtered, input[[ paste0('scat.y.', group) ]]]
              
              filt.dat <- data.frame(x=filt.x, y=filt.y)
              
              ## add
              if(global.param$filt.data == 'Reproducibility')
                filt.txt <-  paste( global.param$filt.data, '\n(alpha=',global.param$repro.filt.val, ')',sep='')
              if(global.param$filt.data == 'StdDev')
                filt.txt <-  global.param$sd.filt.str
              
              p <- p %>% add_trace(
                x=filt.dat$x,
                y=filt.dat$y,
                type='scatter', mode='markers', 
                marker=list(size=10, color=rep('blue', nrow(filt.dat))), 
                opacity=.2, text=values.filtered, name= filt.txt
              )
              
            }
            
            ## add Pearson correlation, based on all data points
            

            ## some astethics
            p <- layout(p,
                        title=paste0( 'Pearson r: ', round(cor.xy.pears, 3), ' / Spearman r: ', round(cor.xy.spear, 3)),
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
            IDs <- res[ , 'id.concat']
            
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
            
            ##change IDs depending on label
            ##needed to plot PPI interactors in the same style
            if(input[[paste('volc.label', group, sep='.')]]=="ID_Symbol" & !is.null(global.results$id.map ) ){
              IDs <- res[ , 'id.concat']
            }else if(input[[paste('volc.label', group, sep='.')]]=="Symbol" & !is.null(global.results$id.map ) ){
              IDs <- res[ , 'id.mapped']
            }else{
              IDs <- res[ , 'id']
            }

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
            if(global.param$which.test == 'Two-sample mod T'){
              mtext( paste("log(", sub('.*\\.vs\\.', '', group), "/", sub('\\.vs.*', '', group),")"), side=1, cex=cex.axis, line=3)
            }else{
              mtext(expression(log(FC)), side=1, cex=cex.axis, line=3)
            }
            # if(global.param$filter.type=="adj.p"){
            #   mtext(expression(-10*log[10](adj.p-value)), side=2, cex=cex.axis, line=3)
            # }else{
            #   mtext(expression(-10*log[10](p-value)), side=2, cex=cex.axis, line=3)
            # }
            mtext(expression(-10*log[10](p-value)), side=2, cex=cex.axis, line=3)
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
                text( xlim[2]-(xlim[2]*.2), y0, paste(global.param$filter.type, "=", global.param$filter.value, ",", "logFC", "=", input[[paste( "ppi.min.fc", group, sep='.')]], sep=''), pos=1, col='grey30')
            }

            ## ###################################
            ## add filter
            ## minimal log P-value for given filter
            if(length(sig.idx) > 0 & !(input[[paste('ppi.hyper.filt', group, sep='.' )]])){
                
              filt.minlogPVal <- min(logPVal[names(sig.idx)], na.rm=T)
            
                abline(h=filt.minlogPVal, col=my.col2rgb('grey30', 50), lwd=2, lty='dashed')
                text( xlim[2]-(xlim[2]*.1), filt.minlogPVal, paste(global.param$filter.type, global.param$filter.value, sep='='), pos=1, col='grey30')
            }

            ## number of significant
            legend(ifelse(PPI, 'topright', 'top'), bty='n', legend=paste(filter.str, '\nsig / tot: ', length(sig.idx),' / ', sum(!is.na(logFC) & !is.na(logPVal)), sep=''), cex=cex.leg)

            ## ############################
            ## indicate directionality for two-sample tests
            if(global.param$which.test == 'Two-sample mod T'){
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
                legend('topleft', legend=leg, col=leg.col, pch=20, cex=cex.leg, title=ifelse( length(ppi.int.idx) > 0 ,paste('interactors sig/det/tot') , ''), pt.cex=cex.vec[1])

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
        output$expr.profile <- renderPlotly({
            if(is.null(global.results$data)) return()

            ## dataset
            if(is.null(global.results$table.log)){
                tab <- data.frame(global.input$table)
                #tab <- data.frame(global.results$table.na.filt)
            } else {
                tab <- data.frame(global.results$table.log)
            }
            ## id column
            id.col.value <- global.param$id.col.value
            ## group vector
            if(global.param$norm.per.group){
              grp <- global.param$grp.norm
            }else{
              grp <- global.param$grp
            }
            ## group colors
            if(global.param$norm.per.group){
              grp.col <- global.param$grp.colors.norm
              grp.col.leg <- global.param$grp.colors.legend.norm
            }else{
              grp.col <- global.param$grp.colors
              grp.col.leg <- global.param$grp.colors.legend
            }

            # update selection
            tab <- tab[, c(id.col.value, names(grp))]
            grp.col <- grp.col[names(grp)]
            grp.col.leg <- grp.col.leg[unique(grp)]
              
            withProgress({
                   setProgress(message = 'Processing...', detail= 'Generating profile plots')
                   p <- makeProfileplot(tab, id.col.value, grp, grp.col, grp.col.leg, 
                                   xlim.mode = input$profile.plot.xlim,
                                   main='Before normalization',
                                   plotly = T)
            })
            ## plot
            p

        })
        ###########################
        ## normalized
        output$expr.profile.norm <- renderPlotly({
            if(is.null(global.results$data)) return()

            ## dataset
            tab <- data.frame(global.results$table.norm)

            ## id column
            id.col.value <- global.param$id.col.value
            ## group vector
            if(global.param$norm.per.group){
              grp <- global.param$grp.norm
            }else{
              grp <- global.param$grp
            }
            ## group colors
            if(global.param$norm.per.group){
              grp.col <- global.param$grp.colors.norm
              grp.col.leg <- global.param$grp.colors.legend.norm
            }else{
              grp.col <- global.param$grp.colors
              grp.col.leg <- global.param$grp.colors.legend
            }

            # update selection
            tab <- tab[, c(id.col.value, names(grp))]
            grp.col <- grp.col[names(grp)]
            grp.col.leg <- grp.col.leg[unique(grp)]
            
            
            withProgress({
                   setProgress(message = 'Processing...', detail= 'Generating profile plots')
                   p <- makeProfileplot(tab, id.col.value, grp, grp.col, grp.col.leg, 
                                   xlim.mode = input$profile.plot.xlim,   
                                   main=paste(global.param$norm.data, 'normalized'),
                                   plotly = T)
            })
            # plot
            p
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
                #tab <- data.frame(global.results$table.na.filt)
                tab <- data.frame(global.input$table)
            else
                tab <- data.frame(global.results$table.log)

            ## id column
            id.col.value <- global.param$id.col.value
            ## group vector
            if(global.param$norm.per.group){
              grp <- global.param$grp.norm
            }else{
              grp <- global.param$grp
            }
            ## group colors
            if(global.param$norm.per.group){
              grp.col <- global.param$grp.colors.norm
              grp.col.leg <- global.param$grp.colors.legend.norm
            }else{
              grp.col <- global.param$grp.colors
              grp.col.leg <- global.param$grp.colors.legend
            }
            
            
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
            if(global.param$norm.per.group){
              grp <- global.param$grp.norm
            }else{
              grp <- global.param$grp
            }
            ## group colors
            if(global.param$norm.per.group){
              grp.col <- global.param$grp.colors.norm
              grp.col.leg <- global.param$grp.colors.legend.norm
            }else{
              grp.col <- global.param$grp.colors
              grp.col.leg <- global.param$grp.colors.legend
            }
            
            
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
        width = function(){120*(ncol(data.frame(global.results$table.na.filt))-1)},
        height= function(){120*(ncol(data.frame(global.results$table.na.filt))-1)}
        #width = function(){120*(ncol(data.frame(global.input$table))-1)},
        #height= function(){120*(ncol(data.frame(global.input$table))-1)}
        )

        ###############################
        ## actual plot
        plotMultiScatter <- function(define.max, min.val, max.val, robustify, update){

            cat('\n-- plotMultiScatter  --\n')
            
            ## dataset
            if(is.null(global.results$table.log))
                tab <- data.frame(global.results$table.na.filt)
                #tab <- data.frame(global.input$table)
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
              
              ## calculate correlation matrix
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
            
            ###############################
            ## plot
            repro.filt=global.results$values.filtered
            
            withProgress({
                setProgress(message = 'Processing...', detail= 'Calculating correlations')
                suppressWarnings(my.multiscatter(tab, 
                                cm=cm,
                                repro.filt=global.results$values.filtered, 
                                grp=grp,  
                                grp.col.legend=global.param$grp.colors.legend[unique(grp)], 
                                define.max=define.max, 
                                max.val=max.val, 
                                min.val=min.val, 
                                robustify=robustify,
                                update=update))
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
                    tab <- data.frame(global.results$table.na.filt)
                    #tab <- data.frame(global.input$table)
                  else
                    tab <- data.frame(global.results$table.log)
               
                  ## id column
                  id.col.value <- global.param$id.col.value
               
                  ## group vector
                  grp <- global.param$grp
                  grp <- sort(grp)
                  ## group colors
                  grp.col <- global.param$grp.colors
                  grp.col.leg <- global.param$grp.colors.legend
               
              
                  ## update selection
                  tab <- tab[, c(id.col.value, names(grp))]
                  grp.col <- grp.col[names(grp)]
                  grp.col.leg <- grp.col.leg[unique(grp)]
                
                  
                  ## get correlation matrix
                  if(is.null(global.results$cm) | global.param$update.cm == TRUE){
                    
                    
                    ## calculate correlation matrix
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
            #tab <- data.frame(global.input$table)
            tab <- data.frame(global.results$table.na.filt)
            
            ## id column
            id.col <- global.param$id.col.value
            ## class vector
            grp <- sort(global.param$grp)
            grp.col.legend <- global.param$grp.colors.legend

            ## update selection
            tab <- tab[, c(id.col, names(grp))]
            grp.col.legend <- grp.col.legend[unique(grp)]
            
            ## table
            tab <- tab[, names(grp)]

            ## ############################################
            ## get correlation matrix
            if(is.null(global.results$cm) | global.param$update.cm == TRUE){
              
              
              ## calculate correlation matrix
              withProgress(message = 'Correlation matrix...',{
                cm=calculateCorrMat( tab=tab,
                                     grp=grp,
                                     lower= global.plotparam$cm.lower, upper= global.plotparam$cm.upper
                )
              })
              
              ## store results
              global.results$cm <- cm
              global.param$update.cm <- FALSE
              
            } else {
              cm <- global.results$cm
            }
            #########################################
            ## plot
            plotCorrBox(cm, grp, grp.col.legend, global.plotparam$cm.upper)
   
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
            if (!is.null(global.param$grp.selection)){
              res = res[, names(global.param$grp.selection)]
              grp.hm <- global.param$grp.selection
            }else{
              res = res[, names(global.param$grp)]
              grp.hm <- global.param$grp
            }

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
                    plotHM(res=res, hm.rownames=hm.rownames, grp=grp.hm, grp.col=global.param$grp.colors, grp.col.legend=global.param$grp.colors.legend,  hm.clust=input$hm.clust, hm.title=hm.title, hm.scale=input$hm.scale, cellwidth=cw, fontsize_row=input$cexRow, fontsize_col=input$cexCol, max.val=input$hm.max.val, style=global.param$which.test, anno.col=anno.col, anno.col.color=anno.col.color, show.rownames=input$hm.show.rownames, show.colnames=input$hm.show.colnames,
                           height=min( dynamicHeightHM( nrow(global.results$filtered)), 1200 ),
                           width=dynamicWidthHM(length(global.param$grp)))
                    
                  })
            } else {
                 withProgress({
                   setProgress(message = 'Processing...', detail= 'Generating Heatmap')
                   plotHM(res=res, hm.rownames=hm.rownames, grp=grp.hm, grp.col=global.param$grp.colors, grp.col.legend=global.param$grp.colors.legend,  hm.clust=input$hm.clust, hm.title=hm.title, hm.scale=input$hm.scale, cellwidth=cw, fontsize_row=input$cexRow, fontsize_col=input$cexCol, style=global.param$which.test, anno.col=anno.col, anno.col.color=anno.col.color, show.rownames=input$hm.show.rownames, show.colnames=input$hm.show.colnames,
                          height=min( dynamicHeightHM( nrow(global.results$filtered)), 1200 ),
                          width=dynamicWidthHM(length(global.param$grp)))
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
        # output$HM.int <- renderPlotly({
        #   
        #   ## if(is.null(global.results$data)) return()
        #   if(!is.null(error$msg)) return()
        #   
        #   ######################################
        #   ## extract results
        #   res = global.results$filtered
        # 
        #   ######################################
        #   ## require at least three significant hits
        #   validate(need(nrow(res) > 1, 'Need at least 2 features to draw a heatmap!'))
        #   
        #   #######################################
        #   ## heatmap title
        #   hm.title <- paste('filter:', global.param$filter.type, ' / cutoff:', global.param$filter.value, sep='')
        #   hm.title <- paste(hm.title, '\nsig / total: ', nrow(res), ' / ', nrow( global.results$data$output ), sep='')
        #   
        #   ###################################
        #   ## ids to show in heatmap
        #   hm.rownames <- res[, 'id.concat']
        #   
        #   #######################################
        #   ## extract expression values
        #   res = res[, names(global.param$grp)]
        #   
        #   ##@#####################################
        #   ##  dimensions depending on no. rows/columns
        #   cw <- cwHM(ncol(res))
        #   
        #   if(!is.null(global.param$anno.col)){
        #     anno.col=global.param$anno.col
        #     anno.col.color=global.param$anno.col.color
        #   } else {
        #     anno.col=data.frame(Group=global.param$grp)
        #     anno.col.color=list(Group=global.param$grp.colors.legend)
        #   }
        #   
        #   ######################################
        #   ## plot
        #   if(input$hm.int.max){
        #     withProgress({
        #       setProgress(message = 'Processing...', detail= 'Generating Heatmap')
        #       plotHM(res=res, hm.rownames=hm.rownames, grp=global.param$grp, grp.col=global.param$grp.colors, grp.col.legend=global.param$grp.colors.legend,  hm.clust=input$hm.int.clust, hm.title=hm.title, hm.scale=input$hm.int.scale, cellwidth=cw, max.val=input$hm.int.max.val, style=global.param$which.test, anno.col=anno.col, anno.col.color=anno.col.color, plotly = T)
        #     })
        #   } else {
        #     withProgress({
        #       setProgress(message = 'Processing...', detail= 'Generating Heatmap')
        #       plotHM(res=res, hm.rownames=hm.rownames, grp=global.param$grp, grp.col=global.param$grp.colors, grp.col.legend=global.param$grp.colors.legend,  hm.clust=input$hm.int.clust, hm.title=hm.title, hm.scale=input$hm.int.scale, cellwidth=cw, style=global.param$which.test, anno.col=anno.col, anno.col.color=anno.col.color, plotly = T)
        #       
        #     })
        #   }
        # }
        # )
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
            if(!global.param$which.test %in% c('mod F','none')){
                par(mfrow=c(length(groups.comp),1))
                for(g in groups.comp){
                      pval <- res[, paste('P.Value', g, sep='.')]
                      hist(pval, breaks=50, main=paste('Number of tests: ', sum(!is.na(pval)), '',sep=''), xlab='p-value', cex.main=2.5, cex.axis=1.8, cex.lab=1.8, col='darkblue', border=NA)
                      legend('top', legend=g, cex=2)
                }
            ############################################
            ## mod F
            } else {
                pval <- res[, paste('P.Value')]
                hist(pval, breaks=50, main=paste('Number of tests: ', sum(!is.na(pval)), '',sep=''), xlab='p-value', cex.main=2.5, cex.axis=1.8, cex.lab=1.8, col='darkblue', border=NA)
            }

        },
        width = function(){ width=1000},
        height= function(){ height=500*ifelse( global.param$which.test != 'mod F', length(unique(global.param$grp.comp)), 1 )} )

        ######################################################################################
        ##
        ##                                 PCA
        ##
        ######################################################################################
        
        ## shown under 'Explained variance'-tab
        output$run.pca <- renderText({

            if(is.null(global.results$data) | is.na(global.param$filter.value)) return()
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
        ##             PCA plotly 2D scatter
        ##
        ##
        output$pcaxy.plotly <- renderPlotly({
          
          if(is.null(global.results$data) | is.na(global.param$filter.value)) return()
          
          pca <- plotlyPCA( reactiveValuesToList(global.param), reactiveValuesToList(global.results), reactiveValuesToList(global.input),
                            input$pca.x, input$pca.y, grp.other=input$pca.grp.col)
          
          ## update if PCA has been recalculated
          if(!is.null(pca$pca)){
            global.results$pca <- pca$pca
            global.param$update.pca <- FALSE
          }
          
          ##plot
          pca$p
          
        })
  
        ################################################
        ## PCA plotly 3D scatterplot
        output$pcaxyz.plotly <- renderPlotly({
            
          if(is.null(global.results$data) | is.na(global.param$filter.value)) return()
          
          pca <- plotlyPCA(reactiveValuesToList(global.param), reactiveValuesToList(global.results), reactiveValuesToList(global.input), 
                           input$pca.x, input$pca.y, input$pca.z, grp.other=input$pca.grp.col)  
         
          ## update if PCA has been recalculated
          if(!is.null(pca$pca)){
            global.results$pca <- pca$pca
            global.param$update.pca <- FALSE
          }
          ## plot
          pca$p
          
           
        })
        
        ####################################################
        ## PCA loadings
       #  output$pca.loadings <- renderPlot({
       # 
       #    if(is.null(global.results$data) | is.na(global.param$filter.value) | is.null(global.results$pca)) return()
       #    if(!is.null(error$msg)) return()
       #    if(length(global.param$grp) < 3) return()
       # 
       #    pca <- global.results$pca
       # 
       #    pca.x <- as.numeric(sub('PC ','', input$pca.x))
       #    pca.y <- as.numeric(sub('PC ','', input$pca.y))
       #    pca.z <- as.numeric(sub('PC ','', input$pca.z))
       # 
       # 
       #    topn=input$pca.load.topn
       # 
       # 
       #    plotPCAloadings( pca, topn, pca.x, pca.y, pca.z )
       # 
       # })

        ####################################################
        ##  PCA loadings as a scatter
        ####################################################
        # output$scatter.pca.loadings <- renderPlot({
        # 
        #         if(is.null(global.results$data) | is.na(global.param$filter.value) | is.null(global.results$pca)) return()
        #         if(!is.null(error$msg)) return()
        #         if(length(global.param$grp) < 3) return()
        # 
        #         pca <- global.results$pca
        #         pca.x <- as.numeric(sub('PC ','', input$pca.x))
        #         pca.y <- as.numeric(sub('PC ','', input$pca.y))
        #         pca.z <- as.numeric(sub('PC ','', input$pca.z))
        #         topn=input$pca.load.topn
        # 
        #         scatterPlotPCAloadings( pca, topn, pca.x, pca.y, pca.z )
        # 
        # })

        ###############################################
        ## static PCA plot
        ###############################################
        plotPCA <- function(pca.x, pca.y, pca.z, plot=T){

            res <- global.results$filtered
          
            ## require at least three rows
            validate(need(nrow(res) > 2, 'Need at least 3 features to perform PC.'))

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

            ## plot
            pca <- my.prcomp.static(t(res[, names(grp)]), col=grp.col, plot=plot, rgl=F, main='', cex.points=5, leg.vec=names(grp.col.leg), leg.col=grp.col.leg)

            return(pca)
        }


})


###########################################################
