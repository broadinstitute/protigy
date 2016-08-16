################################################################################################################
## Filename: global.r
## Created: October 09, 2015
## Author(s): Karsten Krug
##
## Purpose: Shiny-app to perform differential expression analysis, primarily on proteomics data, to perform
##          simple data QC, to interactively browse through the results and to download high-quality result
##          figures.
##
## This file defines the two panels of the user interface.
##
################################################################################################################
library(shiny)

##################################################################
## user interface
##################################################################
##shinyUI(fluidPage(theme = "bootstrap.min.css",
shinyUI(
    fluidPage(

        includeCSS("www/style.css"),

##        tags$src("<script src=\"www/MailtoComposeInGMail.user.js\"></script>"),

        headerPanel('', windowTitle=paste("modT (v",VER,")",sep='')),

        ## Application title
        titlePanel(paste(APPNAME, " (v",VER,") - Differential Expression Analysis using Moderated T-Statistics", sep='')),
        tags$hr(), tags$br(),



        fluidRow(

            ############################################################################
            ##
            ##                            LEFT panel
            ##
            ############################################################################
            column(3, wellPanel(

                          ###############################
                          ## file upload
                          uiOutput("file.upload"),

                          ###############################
                          ## import session
                          uiOutput("import.session"),

                          ###############################
                          ## test data set
                          ##uiOutput("file.testdata"),

                          ###############################
                          ## id column
                          uiOutput("choose.id.column"),

                          ###############################
                          ## define groups
                          uiOutput("define.groups"),

                          ###############################
                          ## show experimental design
                          uiOutput("list.groups"),

                          ###############################
                          ## filter type
                          tags$br(),
                          uiOutput("filter.type"),

                          ## top N
                          ##conditionalPanel(condition = "input['filter.type'] == 'top.n'",  numericInput( "top.n", "Top N features", value=50, min=2, step=1)),
                          conditionalPanel(condition = "input['filter.type'] == 'top.n'",  numericInput( "filter.value.top.n", "Top N features", value=50, min=2, step=1)),

                          ## nom p
                          ##conditionalPanel(condition = "input['filter.type'] == 'nom.p' & input['run.test'] > '0'", numericInput( "p.val", "P-value filter", value=0.01, min=0, max=1, step=1e-3)),
                          conditionalPanel(condition = "input['filter.type'] == 'nom.p'", numericInput( "filter.value.nom.p", "P-value filter", value=0.01, min=0, max=1, step=1e-2)),
                          ## adj p
                          ##conditionalPanel(condition = "input['filter.type'] == 'adj.p'", numericInput( "adj.p", "Corrected P-Value (FDR)", value=0.05, min=0, max=1, step=1e-3)),
                          conditionalPanel(condition = "input['filter.type'] == 'adj.p'", numericInput( "filter.value.adj.p", "Corrected P-Value (FDR)", value=0.05, min=0, max=1, step=1e-2)),

                          ###############################
                          ## export session
                          uiOutput("export.session"),


                          ###############################
                          ## F5 hint
                          tags$br(),
                          htmlOutput('F5hint'),
                          tags$br(),



                          ##########################################################
                          ## Email footer: works for Gmail only
                          ##########################################################
                          HTML(paste('<footer>Ran into problems?  <a href=\"https://mail.google.com/mail/?view=cm&fs=1&to=karsten@broadinstitute.org&su=Shiny%20modTv',VER,'%20problem&body=%0D%0A%0D%0A%0D%0A',paste(rep('-', 30), collapse=''),'%0D%0AMachine:', Sys.info()['nodename'],'%0D%0AApp:', APPNAME, '%0D%0AVersion:', VER, '%0D%0A','\" target=\"_blank\">Send me an email</a></footer>', sep=''))

                      ) ## end wellPanel
                   ),

            ############################################################################
            ##
            ##                            RIGHT panel
            ##
            ############################################################################
            column(9,

                   htmlOutput('error'),

                   htmlOutput('help.start' ),
                   htmlOutput('help.id.column' ),
                   htmlOutput('help.exp.design' ),
                   htmlOutput('help.test' ),

                   htmlOutput('help.results'),

                   uiOutput('navbar')
                   ) ## end column

        ) ## end end fluidRow


    ) ## end fluidPage
) ## end shinyUI
