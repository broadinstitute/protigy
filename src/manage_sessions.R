########################################################################
## 	                    
## Filename: manage_session.R
## Created: January 10, 2018
## Author(s): Karsten Krug
##
## Purpose: shiny module to manage 'protigy' sessions saved on the server
########################################################################

## module ui
manageSessionsUI <- function(id){
  
  ns <- NS(id)
  
  uiOutput(ns("list.sessions"))

}

## module server funciton
manageSessions <- function(input, output, session,  data.dir){

  ## some variables
  data <- reactiveValues()
  sessions.select <- reactiveValues()
  param <- reactiveValues()
  param$init <- FALSE
  
  SHINYDATA <- data.dir
  BACKUP <- paste(SHINYDATA, 'TMP/', sep='')
  
  
  
  ## ####################################
  ## read database/config file
  getSharedSessions <- reactive({
    
    if(is.null(param$n)) return()
    sessions.select$shared <- rep(FALSE, param$n)  
    
  })
  
  
  ## ##############################
  ## get saved sessions
  getSessions <- reactive({
    
    user <- session$user
    search.path <- paste( SHINYDATA, sub('\\@.*', '',user), sep='')
    
    saved.sessions <- grep( '_session.*RData$', dir( search.path, full.names=T, recursive=T ), value=T )
    names(saved.sessions) <-  paste( sub('_.*','', sub('.*/','',saved.sessions)), file.info(saved.sessions)$ctime, sep='_' )
    
    ## order according to date
    time.tmp <- file.info(saved.sessions)$ctime
    saved.sessions <- saved.sessions[ order(time.tmp) ]
    
    data$saved.sessions <- saved.sessions
    param$n <- length(saved.sessions)
    
  })
  
  
  ## #####################################
  # ## observe checkboxes for sharing sessions
  observe({

    if(!param$init) return()
    
    if(is.null(input$shareID1)) return()  
    
    
    if(is.null(param$n)) return()

    n <- param$n
    #vec <- rep(FALSE, n)

    ns <- session$ns
    shared <- sapply(1:n, function(i) ifelse( input[[paste0("shareID",i)]], TRUE, FALSE))
    
    # shared <- sapply(1:n, function(i){
    #     local({
    #       my_i <- i
    #       res <- ifelse( input[[ns(paste0("shareID",my_i))]], 1, 0)
    #       return(res)
    #     })
    # })
    # #cat('test ')
    
    #shared <- rep(F, n)
    #if(!is.null(input$shareID1))
    #  shared[1] <- input$shareID1
    #shared[1] <-T
    
    sessions.select$shared <- unlist(shared)
    
    
    save(sessions.select, shared, n, file='tmp.RData')
    
   # sessionsShowModal()
    
  })

  
  
  sessionsShowModal <- reactive({
    
    
   # if(is.null(sessions.select$shared)) getSharedSessions()
    if(!is.null(param$n)) n <- param$n
  
    ns <- session$ns
      
    shared.sessions.select <- sessions.select$shared
    saved.sessions <- data$saved.sessions
    n <- param$n
    
    saved.sessions.select <- rep(FALSE, n)
    
    
    
    ## modal window
    showModal(modalDialog(
      
      size='m',
      title = paste0("Manage sessions for ", session$user),
      footer = fluidRow(
        column(6),
        column(3, actionButton(ns('delete.sessions'), 'Delete sessions')),
        column(3, modalButton(label='Close'))
      ),
      fluidPage(
        fluidRow(
          column(1, HTML('<b>No.</b>')),
          column(4, HTML('<b>Name</b>')),
          column(1, HTML('<b>Delete</b>')),
          column(1, HTML('<b>Share</b>')),
          if(sum(shared.sessions.select) > 0){
            column(5, HTML('<b>Share with (email)</b>'))
          } else {
            column(5, HTML(''))
          }
        ),
        
        #updateSelectedSessions()
        
        lapply( 1:length(saved.sessions), function(i){
          fluidRow(
            column(1, HTML(paste0(i))),
            column(4, paste0(names( data$saved.sessions )[i])),
            column(1, checkboxInput(inputId = ns(paste0('deleteID', i)), label = NULL, value = saved.sessions.select[i] )),
            column(1, checkboxInput(inputId = ns(paste0('shareID', i)), label =  NULL, value = shared.sessions.select[i])),
            if(shared.sessions.select[i]){
              column(5, textInput(inputId = ns(paste0('sharedWith', i)), label = NULL ))
            } else {
              column(5, HTML(''))
            }
          )
        })
      ),
      easyClose = TRUE,
      fade = FALSE
      #ifelse( param$init, fade = FALSE, fade = TRUE)
    ))
    
    
  })
  
  ## ############################################
  ## list all saved session of the current user
  ## ############################################
  output$list.sessions <- renderUI({
    
    if(is.null(session$user)) return()
    
    #ns <- session$ns
    
    getSessions()
    getSharedSessions()
    
    #saved.sessions <- data$saved.sessions
    #n <- length(saved.sessions) 
    # n <- param$n
    # 
    # saved.sessions.select <- rep(FALSE, n)
    # 
    # if(is.null(sessions.select$shared)) getSharedSessions()
    # 
    #shared.sessions.select <- sessions.select$shared
    #save(shared.sessions.select, file = 'tmp.RData')
    #shared.sessions.select <- rep(FALSE, n)
    #shared.sessions.select[c(2, 10)] <- TRUE
    
    sessionsShowModal()
    
    param$init <- TRUE
    
  })
  
  ## #############################################
  ## observer to delete sessions from the server
  ## #############################################
  observeEvent( input$delete.sessions, {
    
    user <- session$user
    
    
    saved.sessions <- data$saved.sessions
    
    sessions.to.delete <- input$saved.sessions
    
    data$saved.sessions <- saved.sessions[ setdiff(names(saved.sessions), sessions.to.delete) ]
    
    if(length(sessions.to.delete) == 0) return()
    
    sessions.to.delete <- saved.sessions[sessions.to.delete]
    
    sessions.to.delete.dir <- sub('^(.*/).*' , '\\1', sessions.to.delete) ## the entire directory
    
    ## create a user directory if not present already
    if(!file.exists( paste( BACKUP, sub('\\@.*', '', user), sep='') ))
      dir.create(paste( BACKUP, sub('\\@.*', '', user), sep=''))
    
    ## loop over selected sessions and move them
    for(i in 1:length(sessions.to.delete)){
      file.copy(sessions.to.delete.dir[i], paste( BACKUP, sub('\\@.*', '', user), sep=''), recursive=T)
      file.remove( sessions.to.delete[i] )
    }
  })
}

