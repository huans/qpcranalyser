##---- A shiny webapp to input qPCR raw data and get Transcript Abundance and/or Fold Change

##---- load libraries ----
library(shiny)
library(DT)
library(data.table)
library(tidyverse)
library(magrittr)
library(stringr)

##---- UI ----
ui <- navbarPage('qPCR Analysis Tool',
## Input Tab ----                 
                 tabPanel('Input',
                          titlePanel("qPCR analyser"),
   sidebarLayout(
      sidebarPanel(
        helpText('Choose your raw data file and select your reference gene(s) and data columns to calculate relative Transcript Abundance and Fold Change',
                 br(),
                 'N.B.: Input file must have a gene, treatment, replicate and Ct column.'),
        actionButton("show", "Instructions for use"),
        br(),
        br(),
        fileInput("file",
                  label = ("Select a csv or txt file with your raw data"),
                  multiple=FALSE,
                  accept = c("text/csv",
                             "text/comma-separated-values,text/plain",
                             ".csv")),
        selectInput('genecol',
                    'Select Gene Column',
                    choices='',
                    selected='',
                    multiple=FALSE),
        radioButtons("refcount", label = ("Choose number of reference genes"),
                     inline=TRUE,
                     choices = list("1" = 1, "2" = 2, "3" = 3), 
                     selected = 1),
        helpText('Reference gene may not work if there are any spaces in name'),
## Selecting number of ref genes ----
        selectInput("refgene",
                      "Reference Gene",
                      choices="",
                      selected="",
                      multiple=FALSE),
        conditionalPanel(
          condition = "input.refcount == 2 || input.refcount == 3",
          selectInput("refgene2",
                      "Second Reference Gene",
                      choices="",
                      selected="",
                      multiple=FALSE)),
        conditionalPanel(
          condition = "input.refcount == 3",
          selectInput("refgene3",
                      "Third Reference Gene",
                      choices="",
                      selected="",
                      multiple=FALSE)),
## Select other features ----
        selectInput('treatcol',
                    'Select Treatment Column',
                    choices='',
                    selected='',
                    multiple=FALSE),
        selectInput("controlsample",
                    "Select Control Sample",
                    choices="",
                    selected="",
                    multiple=FALSE),
        selectInput('repcol',
                    'Select Replicate Column',
                    choices='',
                    selected='',
                    multiple=FALSE),
        selectInput('ctcol',
                    'Select Ct Column',
                    choices='',
                    selected='',
                    multiple=FALSE),
        radioButtons("disp", "Display",
                     choices = c(Head = "head",
                                 All = "all"),
                     selected = "head"),
        h5('Download the Delta Ct Summary Table'),
        downloadButton("dct_table_DL", "Download"),
        br(),
        h5('Download the Delta-Delta Ct Summary Table'),
        downloadButton("ddct_table_DL", "Download"),
        br(),
        br()
      ),
      mainPanel(
        tableOutput('tablehead'),
        br(),
        tableOutput('dct_table'),
        br(),
        tableOutput('ddct_table'),
        br()
      ),
   )
),
## Raw Data Tab ----
tabPanel('Raw Data Output',
         titlePanel("Raw Data Table"),
         sidebarPanel(
           helpText('Your raw data will display here once you have loaded a suitable file into the analyser.'),
           radioButtons("dispout", "Display",
                        choices = c(Head = "head",
                                    All = "all"),
                        selected = "head")),
         mainPanel(
           tableOutput('table'))),
## Graph Tab ----
tabPanel('Graph Outputs',
         titlePanel('Delta Ct and Delta-Delta Ct Graphs'),
         sidebarPanel(
           helpText('Once you load a data file and select appropriate columns a Transcript Abundance and a Fold Change graph should appear.'),
           br(),
           helpText('Download Transcript Abundance Graph'),
           downloadButton("downloadPlot1", "Download"),
           helpText('Download Fold Change Graph'),
         downloadButton("downloadPlot2", "Download"),
         br(),
         checkboxInput("geneshow", label = "Show all genes in graphs", value = TRUE),
         
         conditionalPanel(
           condition = "input.geneshow == false",
           selectInput("genestoshow", "Genes to Show in Graph",
                       list(""), multiple=TRUE)
         )),
         mainPanel(plotOutput('dctplot'),
                   plotOutput('ddctplot'),
                   br(),
                   br(),
                   br(),
                   br(),
                   br()
         )),
tags$footer("Created by Richard Browne and Charlotte Francois, La Trobe University, 2019. Last updated 26/06/2019.", align = "right", style = "
position:fixed;
min-height: 5vh;
bottom:0;
width:100%;
height:40px;
color: grey;
padding: 10px;
background-color: black;
z-index: 1000;"
)
)

##---- SERVER ----
server <- function(input, output, session) {
## Instructions and citation modals ----
  observeEvent(input$show, {
    showModal(modalDialog(
      title = "About this tool",
      tags$li("This tool determines transcript abundance and fold change using the 2^delta Ct method. Geometric Mean is used when using multiple reference genes"), br(),
      tags$li("Load a .csv or a .txt file with qPCR information containing a Gene, Treatment, Replicate and Cycle Threshold (Ct) column."), br(),
      tags$li("Select the number of reference genes and choose the correct columns using the dropdown menus."), br(),
      tags$li("The delta Ct and delta-delta Ct datatables should auto-fill and can be downloaded."), br(),
      tags$li("Click the Graph Outputs tab to view the relative transcript abundance and the relative fold changes for your genes and treatments."), br(),
      tags$li("If you found this tool useful, please make sure to reference its use and creators."),
      easyClose = TRUE
    ))
  })

## load data file  ----
  data_file <- reactive({
    req(input$file)
    data <- fread(input$file$datapath)
    return(data)
  })

## make data table and set 'head' function ---- 
  output$table <- renderTable({
    req(data_file)
    data <- data_file()
    if(input$dispout == "head") {
      return(head(data, 5))
    }
    else {
      return(data)
    }},caption = 'Raw Data Output Table', caption.placement = getOption("xtable.caption.placement", "top"))
  
## Head of raw data to visualise headers ----
  output$tablehead <- renderTable({
    req(data_file)
    data <- data_file()
      return(head(data, 3))
  },caption = 'Raw Data Output Table', caption.placement = getOption("xtable.caption.placement", "top"))
 
## Delta Ct Table ---- 
  output$dct_table <- renderTable({
    req(deltasummary)
    data <- deltasummary()
    colnames(data) <- c('Treatment','Gene','Delta Ct','Sample Size','Std Err','Transcript Abundance','Percent TrA','Lower TrA Error Pos','Upper TrA Error Pos')
    if(input$disp == "head") {
      return(head(data, 5))
    }
    else {
      return(data)
    }
  },caption = 'Delta Ct Summary', caption.placement = getOption("xtable.caption.placement", "top"))
## DeltaDelta Ct Table ----
  output$ddct_table <- renderTable({
    req(deltasummary)
    data <- ddCT_file()
    colnames(data) <- c('Treatment','Gene','Delta Delta Ct','Sample Size','Std Err','Fold Change','Fold Change Std Err','Lower FC Error Pos','Upper FC Error Pos')
    if(input$disp == "head") {
      return(head(data, 5))
    }
    else {
      return(data)
    }
  },caption = 'Delta Delta Ct Summary', caption.placement = getOption("xtable.caption.placement", "top"))
  
## Fill the 'column' dropdowns ----
  observeEvent( 
    input$file,{
      updateSelectInput(session,
                        'genecol',
                        choices=colnames(data_file()))
      updateSelectInput(session,
                        'treatcol',
                        choices=colnames(data_file()))
      updateSelectInput(session,
                        'ctcol',
                        choices=c('',colnames(data_file())))
      updateSelectInput(session,
                        'repcol',
                        choices=colnames(data_file()))
      })
  observeEvent(
    input$genecol,{
      updateSelectInput(session,
                      'refgene',
                      choices=c(unique(as.character(data_file()[[input$genecol]]))))
  })
  observeEvent(
    input$genecol,{
      updateSelectInput(session,
                        'refgene2',
                        choices=c(unique(as.character(data_file()[[input$genecol]]))))
    })
      observeEvent(
        input$genecol,{
          updateSelectInput(session,
                            'refgene3',
                            choices=c(unique(as.character(data_file()[[input$genecol]]))))
    })
  observeEvent(
    input$genecol,{
      updateSelectInput(session,
                        'genestoshow',
                        choices=c(unique(as.character(data_file()[[input$genecol]]))))
    })
  observeEvent(
    input$treatcol,{
      updateSelectInput(session,
                        'controlsample',
                        choices=c(unique(as.character(data_file()[[input$treatcol]]))))
    })

## delta CT data ----  
  deltasummary <- reactive({
    req(input$ctcol)
    ref1 <- input$refgene
    ref2 <- input$refgene2
    ref3 <- input$refgene3
    df <- data_file()
    df <- df %>% rename('CT'=input$ctcol,
                        'Treatment'=input$treatcol,
                        'Gene'=input$genecol,
                        'Replicate'=input$repcol)
    df$CT<- as.numeric(as.character(df$CT))
    df <- df[!(is.na(df$Treatment) | df$Treatment==""), ]
    df <- df %>% filter(!CT == 'Undetermined')
    
    if(input$refcount == 1) {
    df1 <- df %>%
    group_by(Replicate, Treatment, Gene) %>% 
    summarise(delta.CT = mean(CT)) %>% 
    mutate_at(vars(matches("delta")), funs(.- .[Gene == ref1])) %>% 
    filter(Gene != ref1)
    } else if(input$refcount == 2) {
      df2 <- df %>% 
        filter(Gene == ref1 | Gene == ref2) %>% 
        group_by(Replicate, Treatment, Gene) %>%
        summarise(delta.CT = mean(CT))
      df3 <- df2 %>% 
        group_by(Replicate, Treatment) %>% 
        summarise(delta.CT = sqrt(prod(delta.CT)))
      df3$Gene <- "GeoMean_REF"
      df1 <- df %>%
        group_by(Replicate, Treatment, Gene) %>% 
        summarise(delta.CT = mean(CT))
      df3 <- df3[c("Replicate", "Treatment", "Gene", "delta.CT")]
      df1 <- df1[c("Replicate", "Treatment", "Gene", "delta.CT")]
      df1 <- rbind(df1,df3)
      df1 <- df1 %>%
        mutate_at(vars(matches("delta")), funs(.- .[Gene == 'GeoMean_REF'])) %>% 
        filter(Gene != "GeoMean_REF", Gene != ref1, Gene != ref2)
    } else {
      df2 <- df %>% 
        filter(Gene == ref1 | Gene == ref2 | Gene == ref3) %>% 
        group_by(Replicate, Treatment, Gene) %>%
        summarise(delta.CT = mean(CT))
      df3 <- df2 %>% 
        group_by(Replicate, Treatment) %>% 
        summarise(delta.CT = sqrt(prod(delta.CT)))
      df3$Gene <- "GeoMean_REF"
      df1 <- df %>%
        group_by(Replicate, Treatment, Gene) %>% 
        summarise(delta.CT = mean(CT))
      df3 <- df3[c("Replicate", "Treatment", "Gene", "delta.CT")]
      df1 <- df1[c("Replicate", "Treatment", "Gene", "delta.CT")]
      df1 <- rbind(df1,df3)
      df1 <- df1 %>%
        mutate_at(vars(matches("delta")), funs(.- .[Gene == 'GeoMean_REF'])) %>% 
        filter(Gene != "GeoMean_REF", Gene != ref1, Gene != ref2, Gene != ref3)
    }

    df1 <- df1 %>% group_by(Treatment, Gene) %>% 
      summarise_each(funs(mean(., na.rm = T), 
                          n = sum(!is.na(.)), 
                          se = sd(., na.rm = T)/sqrt(sum(!is.na(.)))), 
                     delta.CT) %>% 
      mutate(TrA = 2^-mean, # transcript abundance (TrA)
             perc.TrA = TrA * 100, # percentage TrA
             lw.TrA = 2^ - (mean + se),  # lower error position
             up.TrA = 2^ - (mean - se)) %>%  # upper error position
      rename(dCT = mean) %>% 
      mutate_if(is.numeric, round, 5)
})
  
## ddct data ----  
  ddCT_file <- reactive({
    req(input$ctcol)
    df <- deltasummary()
    cont_samp <- input$controlsample
  df2 <- df %>% #####
    group_by(Gene) %>%
    mutate_at(vars(matches("dCT")), funs(.- .[Treatment == cont_samp])) %>% #specify control treatment
    select(1:5) %>% 
    rename(ddCT = dCT) %>% 
    mutate(fc = 2^-ddCT, # fold change ddCT
           fc.se = se^2) %>% # SE fold change
    mutate_at(vars(matches("fc.se")), funs(sqrt(.+. [Treatment == cont_samp]))) %>% 
    mutate(lw.fc = 2^ - (ddCT + fc.se),  # lower error position
           up.fc = 2^ - (ddCT - fc.se)) %>% # upper error position
    mutate_if(is.numeric, round, 5)
  df2
  })
  
## Create dCt Plot ----
  dctplot <- reactive({
    req(input$ctcol)
    data <- deltasummary()
    ref1 <- input$refgene
    ref2 <- input$refgene2
    ref3 <- input$refgene3
    if(input$geneshow == TRUE) {
      df=data
    } else {
      df = subset(data, Gene %in% input$genestoshow)
      }
    p <- ggplot(data = df, 
         aes(x = Gene, 
             y = TrA, 
             fill = Treatment)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.7), 
             colour = "black", width=0.6) +
    geom_errorbar(aes(ymin = lw.TrA, ymax = up.TrA, width = 0.1), # add standard error bars
                  position = position_dodge(width = 0.7)) +
      ggtitle('Transcript Abundance plot') +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.position="right",
            panel.background = element_blank(),
            axis.line = element_line(colour = "black")) +
      scale_y_continuous(expand = c(0, 0))
    if(input$refcount == 1) {
     p <- p + ylab(paste('Transcript Abundance relative to',ref1))
    } else if(input$refcount == 2) {
     p <- p +ylab(paste('Transcript Abundance relative to \n Geometric mean of', ref1, "and", ref2))
    } else {
    p <- p +ylab(paste('Transcript Abundance relative to \n Geometric mean of', ref1, "and", ref2, "and", ref3))
    }
  })
## Output dCt Plot ----
  output$dctplot <- renderPlot({
    req(dctplot)
    dctplotPrint <- dctplot()
    dctplotPrint
  })
## Create ddCt Plot ---- 
  ddctplot <- reactive({
    req(input$ctcol)
    data <- ddCT_file()
    
    if(input$geneshow == TRUE) {
      df=data
    } else {
      df = subset(data, Gene %in% input$genestoshow)
    }
    p <- ggplot(data = df, 
                aes(x = Gene, 
                    y = fc, 
                    fill = Treatment)) +
      geom_bar(stat = "identity", position = position_dodge(width = 0.7), 
               colour = "black", width=0.6) +
      geom_errorbar(aes(ymin = lw.fc, ymax = up.fc, width = 0.1), # add standard error bars
                    position = position_dodge(width = 0.7)) +
      ylab(paste('Fold Change relative to treatment',input$controlsample)) +
      ggtitle('Fold Change Plot') +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.position="right",
            panel.background = element_blank(),
            axis.line = element_line(colour = "black")) +
      scale_y_continuous(expand = c(0, 0))
  })
## Output ddCt Plot ----
  output$ddctplot <- renderPlot({
    req(ddctplot)
    ddctplotPrint <- ddctplot()
    ddctplotPrint
  })
## Download DCt Plot ----  
  output$downloadPlot1 <- downloadHandler(
    filename = function(){paste('DeltaCtPlot.png')},
    content = function(file) {
      ggsave(file, plot = dctplot(), device = "png", width=8, height=6, units='in')
    }
  )
## Download DDCt Plot ----  
  output$downloadPlot2 <- downloadHandler(
    filename = function(){paste('DeltaDeltaCtPlot.png')},
    content = function(file) {
      ggsave(file, plot = ddctplot(), device = "png", width=8, height=6, units='in')
    }
  )
## Download DCt Table ----
  output$dct_table_DL <- downloadHandler(
    filename = function(){paste('DeltaCtSummary.csv')},
    content = function(file) {
      write.csv(deltasummary(), file, row.names = FALSE)
    }
  )
## Download DDCt Table ----
  output$ddct_table_DL <- downloadHandler(
    filename = function(){paste('DeltaDeltaCtSummary.csv')},
    content = function(file) {
      write.csv(ddCT_file(), file, row.names = FALSE)
    }
  )
} # end server

##---- RUN APP ----
shinyApp(ui = ui, server = server)