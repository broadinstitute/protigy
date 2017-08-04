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

shinyUI(

    ##fluidPage(theme="custom.css",
    dashboardPage( skin='blue',


                  ## ##########################################
                  ## header
                  dashboardHeader( title=paste("modT (v",VER,")", sep=""),

                        ## logged user
                        ##dropdownMenuOutput('memfree'),
                        dropdownMenuOutput('session.name'),
                        dropdownMenuOutput('logged.user'),
                        dropdownMenuOutput('logout'),

                        ## E-mail
                        tags$li(class = "dropdown",
                                a(href=paste('https://mail.google.com/mail/?view=cm&fs=1&to=', MAIL,'&su=Shiny%20modTv',VER,'%20help&body=%0D%0A%0D%0A%0D%0A',paste(rep('-', 30), collapse=''),'%0D%0AMachine:', Sys.info()['nodename'],'%0D%0AApp:', APPNAME, '%0D%0AVersion:', VER, '%0D%0A', sep=''), target="_blank",
                                  "Help me!",
                                  img(src="help.jpg", height="20px", alt="Help"),  style = "padding-top:15px; padding-bottom:0px"
                                  )
                                ),
                        ## logo
                       tags$li(class = "dropdown",
                                a(href="http://www.broadinstitute.org/proteomics", target="_blank",
                                  img(height = "50px", alt="Proteomics Logo", src="BroadProteomicsLogo.png"),
                                  style = "padding-top:0px; padding-bottom:0px"
                                  )
                                )),


                       ##),

            ############################################################################
            ##
            ##                            LEFT panel
            ##
            ############################################################################
        ##    column(3, fluidRow(wellPanel(
                  dashboardSidebar(
                      tags$head(tags$style(HTML('.shiny-server-account { display: none; }'))),
                      tags$head(tags$style(".wrapper {overflow: visible !important;}")),

                          ###############################
                          ## file upload
                          uiOutput("file.upload"),

                          ###############################
                          ## import session
                          uiOutput("import.session"),

                          ## #############################
                          ## browse saved sessions
                          uiOutput("browse.sessions"),

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
                          uiOutput("filter.value"),

                          ###############################
                          ## F5 hint
                          tags$br(),
                          htmlOutput('F5hint'),
                          tags$br()



                          ##########################################################
                          ## Email footer: works for Gmail only
                          ##########################################################
                          ##HTML(paste('<footer>Ran into problems?  <a href=\"https://mail.google.com/mail/?view=cm&fs=1&to=karsten@broadinstitute.org&su=Shiny%20modTv',VER,'%20problem&body=%0D%0A%0D%0A%0D%0A',paste(rep('-', 30), collapse=''),'%0D%0AMachine:', Sys.info()['nodename'],'%0D%0AApp:', APPNAME, '%0D%0AVersion:', VER, '%0D%0A','\" target=\"_blank\">Send me an email</a></footer>', sep=''))

               ),

            ############################################################################
            ##
            ##                            RIGHT panel
            ##
            ############################################################################
            dashboardBody(
            ##column(9,

                ##tags$head(includeCSS('www/custom.css')),

                useShinyjs(),
                inlineCSS(appCSS),

                ## #######################
                ## loading animantion
                div(
                    id = "loading-content",
                    h2("Loading...")
                ),
                ## #######################
                ## app content
               ## hidden(
                    div( id= "app-content",

                         htmlOutput('error'),

                         printHTMLUI('getting.started'),
                         printHTMLUI('change.log'),
                         printHTMLUI('id.column'),
                         printHTMLUI('exp.design'),
                         printHTMLUI('analysis'),
                         printHTMLUI('results'),

                         ## the actual content
                         uiOutput('navbar'),


                         tags$head(tags$style(HTML('

                              /* main sidebar */
                              .skin-blue .main-sidebar {
                                      background-color: #000000;
                              }
                              /* logo */
                              .skin-blue .main-header .logo {
                                      background-color: #1477C5; ## Broad blue

                              }
                              /* rest of the header*/
                              .skin-blue .main-header .navbar {
                                      background-color: #1477C5;
                              }

                              # /* scroll bar  */
                              # .sidebar {
                              #          height: 90vh; overflow-x: auto
                              # }



                              # /* active selected tab in the sidebarmenu */
                              # .skin-blue .main-sidebar .sidebar .sidebar-menu .active a{
                              #          background-color: #ECF0F5;
                              # }
                              # /*  */
                              # .skin-blue .content-wrapper .right-side {
                              #          background-color: #FFFFFF;
                              # }
                   ')))
                   ) ## end div

               ## ) ## end hidden

            )## dashboard body

        ) ## end dashboardPage

) ## end shinyUI
