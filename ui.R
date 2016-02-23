library(shiny)
##library(shinyRGL)


shinyUI(fluidPage(

    ## Application title
    titlePanel("Moderated T-test (developmental version)"),
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
        ## test data set
        uiOutput("file.testdata"),

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

        ## nom p
        conditionalPanel(condition = "input['filter.type'] == 'nom.p' & input['run.test'] > '0'", numericInput( "p.val", "P-value filter", value=0.01, min=0, max=1, step=1e-3)),
        ## adj p
        conditionalPanel(condition = "input['filter.type'] == 'adj.p'", numericInput( "adj.p", "Corrected P-Value (FDR)", value=0.05, min=0, max=1, step=1e-3)),
        ## top N
        conditionalPanel(condition = "input['filter.type'] == 'top.n'",  numericInput( "top.n", "Top N features", value=50, min=2, step=1)),


        ###############################
        ## F5 hint
        tags$br(),
        htmlOutput('help')

        )),

      ############################################################################
      ##
      ##                            RIGHT panel
      ##
      ############################################################################
      column(9,
             htmlOutput('help.start' ),
             htmlOutput('help.id.column' ),
             htmlOutput('help.exp.design' ),
             htmlOutput('help.test' ),

             htmlOutput('help.results'),

             uiOutput('navbar')
      ) ## end column
  ) ## end row
))
