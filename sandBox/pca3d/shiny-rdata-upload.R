library(shiny)
library(plotly)
library(DT)
library(dplyr)
options(shiny.maxRequestSize = 50*1024^2)
# # Define two datasets and store them to disk
# x <- rnorm(100)
# save(x, file = "x.RData")
# rm(x)
# y <- rnorm(100, mean = 2)
# save(y, file = "y.RData")
# rm(y)

# Define UI
ui <- shinyUI(navbarPage(title = "RNA-seq: Unsupervised Discovery",
  
                tabPanel(title = "Main",
                         sidebarLayout(
                           sidebarPanel(
                              fileInput("file", label = ""),
                              actionButton(inputId="run","RUN"),
                              plotOutput("hist"), width = 3), 
                         
                        mainPanel(
                         DT::dataTableOutput("mytable1"),uiOutput("groups")
                         )
                )
                ),
                
                tabPanel(title = "Clust", plotOutput(outputId = "clust", width = 1080, height = 1080))
            
    )
    )


# Define server logic
server <- shinyServer(function(input, output) {
  
  observeEvent(input$run,{
    if ( is.null(input$file)) return(NULL)
    inFile <- input$file
    file <- inFile$datapath
    # load the file into new environment and get it from there
    e = new.env()
    name <- load(file, envir = .GlobalEnv)

    # Plot the data
    output$mytable1 <- DT::renderDataTable({
      DT::datatable(target)
      })
    
    output$groups <- renderUI({
      tagList(
        selectInput(inputId = "choice", label = "Select Groupings", 
                    choices = colnames(target), selected = "group"))
  
      })
    
    })
})

# Run the application 
shinyApp(ui = ui, server = server)
