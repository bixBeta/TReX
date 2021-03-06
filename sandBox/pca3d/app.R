library(plotly)
library(shiny)

source("global.R")
ui <- fluidPage(
  sidebarLayout(
    sidebarPanel(
      
      # verbatimTextOutput("hover"),
      # verbatimTextOutput("click"),width = 3,
      
      selectInput(inputId = "choice", label = "Select Groupings", 
                  choices = colnames(pComp.df), selected = "group"), width = 2),
      
      mainPanel(
        tabsetPanel(
          id = 'TABS',
          tabPanel("plotPCA", plotlyOutput("plot"), width = 4),
          
          tabPanel("Clust", plotOutput(outputId = "clust", width = 1080, height = 1080))
        )
        
    #     ),
    # 
    # 
    # mainPanel(plotlyOutput("plot"), width = 4)
    
  )))


server <- function(input, output, session) {
  
  output$plot <- renderPlotly({
    plot_ly(pComp.df, x = ~ pca.x...1., y = ~ pca.x...2., z = ~ pca.x...3., 
            color = ~ get(input$choice)) %>%
      layout(autosize = F, width = 800, height = 800, 
             scene = list(xaxis = list(title = 'PC1'),
                          yaxis = list(title = 'PC2'),
                          zaxis = list(title = 'PC3')))
  })
  
  output$clust <- renderPlot(plot(hc, hang=-1, ylab="Height", las=2, 
                                  xlab="Method: Euclidean distance - Ward criterion", 
                                  main="Cluster dendrogram")
  )
  
  # output$hover <- renderPrint({
  #   d <- event_data("plotly_hover")
  #   if (is.null(d)) "Hover events appear here (unhover to clear)" else d
  # })
  # 
  # output$click <- renderPrint({
  #   d <- event_data("plotly_click")
  #   if (is.null(d)) "Click events appear here (double-click to clear)" else d
  # })
  
}

shinyApp(ui, server)
