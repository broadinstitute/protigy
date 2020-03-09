library(shiny)

shinyExportExcelUI <- function(id, signif.col, data) {
  
  ns <- NS(id)
  # fluidPage(
  #   fluidRow(
  #    column(6, 
  #           numericInput(inputId = ns('profile.pval'), label = 'max. p-value', min = 0, max = 1, value = 0.01, step = 0.01, width = 100)
  #    ),
  #    column(6,
  #           selectInput(inputId = ns('profile.signifColumn'), label = 'based on', choices = signif.col, width = 200)
  #    )
  #   ),
  #   fluidRow(column(12, selectizeInput(ns('profile.selected'), 
  #                                      choices=NULL, 
  #                                      label='Select protein/p-site/pathway', 
  #                                      selected=NULL,width = 1000
  #                                      ))),
  # 
  #   fluidRow(column(6, htmlOutput(ns('ext.link'))), column(6)),
  #   fluidRow(column(12, br())),
  #   fluidRow( plotlyOutput( ns("profile")) ) 
  # )
}

shinyExportExcel <- function(input, output, session, data, anno, label, signif.col){
  
  
  observeEvent(c(input$profile.pval, input$profile.signifColumn), {
    dat <- data[[label]]
    pval=as.numeric(input$profile.pval)
    if(is.na(pval)) pval=0
    sig.col = input$profile.signifColumn
    sig.idx <- dat[which( dat[, sig.col] < pval ),'id']
    updateSelectizeInput( session=session, inputId = 'profile.selected', label='Select protein/p-site/pathway', choices=sig.idx, selected=sig.idx[1], server=T)
  })
 
  
  
  # -----------------------------------------
  # define heatmap function
  plotProfile <- function(data, anno, label, selected, pval, sig.col, signif.col){
    dat <- data[[label]]
    rownames(dat) <- dat$id
    ann <- anno[[label]]
    rownames(ann) <- ann$expt
    ann <- ann[order(ann$repeated.measure), ]
    col.idx <- sig.col
    if(is.na(pval)) pval=0
    
    # extract p-values
    dat.p <- dat[, signif.col]
    
    # extract significant features
    #sig.idx <- which( dat[, sig.col] < pval )
    #dat <- dat[sig.idx, as.character( ann$expt)]
    dat <- dat[, as.character( ann$expt)]
    
    # find selected features
    selected.idx <- rownames(dat)[ rownames(dat) %in% selected ]
    
    validate(
      need(length(selected.idx) > 0, 'No significant features.')
    )
    
    dat <- dat[ selected.idx, ]
    dat.p <- dat.p[selected.idx, ]
    dat.p <- format(dat.p, scientific = T, digits = 3)
    dat.p <- sapply(names(dat.p), function(x) paste(x, dat.p[x], sep='='))
    
    dat.plot <- data.frame(time=ann$repeated.measure, feature=unlist(dat) )
    dat.plot.1 <- dat.plot[grep('Rep1', ann$expt),]
    dat.plot.2 <- dat.plot[grep('Rep2', ann$expt),]
    
    #View(dat.plot.1)
    # column annotation
    #ann.plot <- ann[, -c(1)]
    # set up the plot
    p <- plot_ly( x=dat.plot.1$time, y=dat.plot.1$feature, type = 'scatter', mode='line+markers', name='Rep1') 
    p <- p %>% add_trace(  x=dat.plot.2$time, y=dat.plot.2$feature, type = 'scatter', mode='line+markers', name='Rep2')
    p <- p %>% layout(xaxis=list( title='Time'), yaxis=list(title='Expression'), title=paste(dat.p, collapse='\n'))#, annotations=list(text=dat.p, yref='paper', xref='paper', y=1, x=1))
   # p <- p %>% add_annotations(p, text=dat.p, )
    p <- p %>% config(displayModeBar = FALSE) %>% config(showLink = FALSE)
    p
    }
  
  output$ext.link <- renderText({
    
    id <- input$profile.selected

    src <- 'UNIPROT'    
    # check whether id is UniProt or MSigDB
    if(length(grep('KEGG|PID|BIOCARTA|REACTOME', id)) > 0) 
      src <- 'MSIGDB'
  
    id.l <- unlist(strsplit(id, '_'))
    
    # UniProt
    if(src == 'UNIPROT')
      link <- paste("Link to  UniProt: <a href='https://www.uniprot.org/uniprot/", sub('(_|,|;|\\.).*', '', id.l[2]),"' target='_blank'>", id, "</a>", sep='')
    else
      link <- paste("Link to MSigDB: <a href='http://software.broadinstitute.org/gsea/msigdb/geneset_page.jsp?geneSetName=", id,"' target='_blank'>", id, "</a>", sep='')
    link
  })
  
  output$profile <- renderPlotly({
    plotProfile(data, anno, label=label, pval=as.numeric(input$profile.pval), sig.col = input$profile.signifColumn, selected=input$profile.selected,  signif.col= signif.col)
  
    })
  
}

