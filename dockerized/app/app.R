library("BiocManager")
options(repos = BiocManager::repositories())
getOption("repos")
library(shiny)
library(plotly)
library(DT)
library(dplyr)
options(shiny.maxRequestSize = 100*1024^2)



# Define UI
ui <- shinyUI(navbarPage(title = "RNA-seq: Unsupervised Exploration",
                         
                         tabPanel(title = "Main",
                                  sidebarLayout(
                                    sidebarPanel(
                                      fileInput("file", label = ""),
                                      actionButton(inputId="run","RUN"),
                                      # plotOutput("hist"), 
                                      width = 3), 
                                    
                                    mainPanel(
                                      DT::dataTableOutput("mytable1")
                                    )
                                  )
                         ),
                         
                         tabPanel(title = "Clust", plotOutput(outputId = "clust", width = 1080, height = 1080)),
                         
                         tabPanel("3D-PCA", 
                                  sidebarLayout(
                                    sidebarPanel(
                                      uiOutput("groups"),width = 3), 
                                    
                                    mainPanel(
                                      plotlyOutput("plot"))
                                  )
                         ),
                         
                         tabPanel("2D-PCA", 
                                  sidebarLayout(
                                    sidebarPanel(
                                      uiOutput("groups2"),width = 4,
                                      plotOutput(outputId = "scree")), 
                                    
                                    mainPanel(
                                      plotlyOutput("plot2"))
                                  )
                         )
                         
                         
                         
                         
)
)


# Define server logic
server <- shinyServer(function(input, output, session) {
  
  
  observeEvent(input$run,{
    if ( is.null(input$file)) return(NULL)
    inFile <- input$file
    file <- inFile$datapath
    # load the file into new environment and get it from there
    e = new.env()
    name <- load(file, envir = .GlobalEnv)
    
    progress <- Progress$new(session, min=1, max=20)
    on.exit(progress$close())
    
    progress$set(message = 'Execution in progress',
                 detail = 'This may take a while...')
    for (i in 1:4) {
      progress$set(value = i)
      Sys.sleep(0.5)
    }
    
    # Plot the data
    target2 <- target %>% select(label, group)
    output$mytable1 <- DT::renderDataTable({
      DT::datatable(target2)
    })
    
    # output$groups <- renderUI({
    #   tagList(
    #     selectInput(inputId = "choice", label = "Select Groupings", 
    #                 choices = colnames(target), selected = "group"))
    #   
    # })
    
    
   
    # ################################
    # ################################
    # # DESEQ2
    #
    library("DESeq2")
    library("dplyr")
    library("tidyverse")
    
    #
    ## -------------------------------------------------------------------------------------------------------------------
    dds <- DESeqDataSetFromMatrix(countData = counts,
                                  colData = target,
                                  design = ~ group)
    
    #
    #
    #
    # ## -------------------------------------------------------------------------------------------------------------------
    for (i in 4:12) {
      progress$set(value = i)
      Sys.sleep(0.1)
    }
    
    dds <- DESeq(dds)
    resultsNames(dds)
    #
    # ## -------------------------------------------------------------------------------------------------------------------
    vsd <- varianceStabilizingTransformation(dds, blind=T)
    
    for (i in 13:17) {
      progress$set(value = i)
      Sys.sleep(0.5)
    }
    #
    # ## -------------------------------------------------------------------------------------------------------------------
    rv <- rowVars(assay(dds))
    # select 500 top genes by variance
    select <- order(rv, decreasing=TRUE)[seq_len(min(500, length(rv)))]
    # perform a PCA on the data in assay(x) for the selected genes
    pca <- prcomp(t(assay(dds)[select,]))
    #
    # # the contribution to the total variance for each component
    percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
    cond <- target$group
    intgroup <- "group"
    intgroup.df <- as.data.frame(colData(dds)[, intgroup, drop=FALSE])
    
    #
    # # assembly the data for the plot
    d <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], group=cond, intgroup.df, name=colnames(dds))
    d2 <- data.frame(PC1=pca$x[,1], PC3=pca$x[,3], group=cond, intgroup.df, name=colnames(dds))
    
    # ## -------------------------------------------------------------------------------------------------------------------
    #
    PC1= as.data.frame(pca$x[,1])
    PC2= as.data.frame(pca$x[,2])
    PC3= as.data.frame(pca$x[,3])
    #
    pComp.df <- data.frame(PC1, PC2 , PC3)
    
    # pComp.df <- rownames_to_column(pComp.df, var = "samples")
    # pComp.df <- left_join(pComp.df, target, by = c("samples" = "label"))
    # pComp.df <- column_to_rownames(pComp.df, var = "samples")
    # pComp.df <- pComp.df %>% mutate_at(vars(matches("pca")), ~ ./100000)
    
    pComp.df <- rownames_to_column(pComp.df, var = "samples")
    pComp.df <- left_join(target, pComp.df, by = c("label" = "samples"))
    pComp.df <- pComp.df %>% mutate_at(vars(matches("pca")), ~ ./100000)
    
    hc2 <- hclust(dist(t(assay(vsd))), method="ward.D") 
    
    
    
    
    output$clust <- renderPlot(plot(hc2, hang=-1, ylab="Height", las=2, 
                                    xlab="Method: Euclidean distance - Ward criterion", 
                                    main="Cluster Dendrogram"))
    
    output$groups <- renderUI({
#      if(is.null(colnames(pComp.df))){return()}
      
      tagList(
        selectInput(inputId = "choice", label = "Select Groupings",
                    choices = colnames(pComp.df), selected = "group")
        )
      
      
      
    })
    
    output$groups2 <- renderUI({
      #      if(is.null(colnames(pComp.df))){return()}
      
      tagList(
        selectInput(inputId = "choice2", label = "Select Groupings",
                    choices = colnames(pComp.df), selected = "group")
      )
      
      
      
    })
    
    output$plot <- renderPlotly({
      
      if(is.null(input$choice)){return()}
      plot_ly(data = pComp.df, x = ~ pca.x...1., y = ~ pca.x...2., z = ~ pca.x...3., 
              color = ~ get(input$choice),
              marker = list(size = 8,
                            line = list(color = ~ group, width = 1))) %>%
        layout(autosize = F, width = 800, height = 800, 
               scene = list(xaxis = list(title = 'PC1'),
                            yaxis = list(title = 'PC2'),
                            zaxis = list(title = 'PC3')))
    })
    
    
    output$plot2 <- renderPlotly({
      
      if(is.null(input$choice2)){return()}
      plot_ly(data = pComp.df, x = ~ pca.x...1., y = ~ pca.x...2., 
              color = ~ get(input$choice2),
              marker = list(size = 8,
                            line = list(color = ~ group, width = 1))) %>%
        layout(autosize = F, width = 800, height = 800,
               xaxis = list(title = 'PC1'),
               yaxis = list(title = 'PC2'))
      
    })
    
    output$scree <- renderPlot({
      screeplot(pca, npcs = ncol(pca$x), type = "lines", main = "Var Explained by PC's")

    })

    
    for (i in 18:20) {
      progress$set(value = i)
      Sys.sleep(0.1)
    }
    
    # output$groups <- renderUI({
    #   
    #   if(is.null(input$choice)){return()}
    #   
    #   tagList(
    #     selectInput(inputId = "choice", label = "Select Groupings",
    #                 choices = colnames(pComp.df), selected = "group"))
    #   
    # 
    #   
    # })
    
    
  })
})

# Run the application 
shinyApp(ui = ui, server = server)
