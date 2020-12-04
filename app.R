# Define UI for app that draws a histogram ----
library("shiny")
library("tidyverse")
source("preprocessing.R")
library("shinyMatrix")


countries = read.csv("countries_names.csv")
countries_list = list()
for (c in unique(countries$country_name)){
  if (sum(countries$country_name == c ) >1 ){
    countries_list[[c]]= as.list(sapply(countries$region[countries$country_name == c], function(x){toString(x)}))
  }else{
    countries_list[[c]]= list(c("Main territory"))
  }
}

equivalences <- list("H" = 3,
                   "W" = 1.5,
                   "S" =	4,
                    "C" = 1.3)

ui <- fluidPage(
  
  # App title ----
  titlePanel("Compute the sensitivity of the pool",
             windowTitle="Pool Testing under Correlations"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      numericInput(inputId = "p",
                   label = "How many people per pool?",
                   value = 90,
                   min=0),
      radioButtons(inputId = "existence_Correlations",
                   label = "Do the people in the pool know each other?",
                   choices =  c(  "Not that I know of"= 0,
                                  "Yes" = 1),
                   selected = 0,
                   inline = FALSE),
      tags$hr(),
      conditionalPanel(
        condition = "input.existence_Correlations == 1",
      selectInput(
        "where",
        "Select the participants that know someone else in the pool",
        choices = 1:100,
        selected = NULL,
        multiple = TRUE,
        selectize = TRUE,
        width = NULL,
        size = NULL
      )
      ),
      conditionalPanel(
         condition = "input.where.length > 0",
         matrixInput(
           "Correlation",
           label = "Fill in the type of relationships between your selection of participants. S=spouse, H=same household, C=same suspected cluster, W=coworkers",
           value = diag(rep(1,3)),
           class = "character",
           rows = list(
             extend = FALSE
           ),
           cols = list(
             names = FALSE
           )
         )
         
       ),
      tags$hr(),
      selectInput(inputId = "country",
                  label = "Which country do you live in?",
                  choices = unique(countries$country_name),
                  selected = "United Kingdom"),
      selectInput(inputId = "region",
                  label = "Which region do you live in?",
                  choices = countries_list,
                  selected="Main territory"),
      dateInput(inputId = "date_event",
                label ="When will you be pooling samples?", 
                value ="2020-11-15",
                min = NULL,
                max = "2020-11-29",
                format = "yyyy-mm-dd",
                startview = "month",
                weekstart = 0,
                language = "en",
                width = NULL,
                autoclose = TRUE,
                datesdisabled = NULL,
                daysofweekdisabled = NULL),
      tags$hr(),
      # Input: Select a file ----
      fileInput("file1", "Upload participants' information/risk factors (choose CSV File)",
                multiple = FALSE,
                accept = c("text/csv",
                           "text/comma-separated-values,text/plain",
                           ".csv")),
      
      # Horizontal line ----
      uiOutput("region")
  ),



    # Main panel for displaying outputs ----
  # Main panel for displaying outputs ----
  mainPanel(
    
    # Output: Histogram ----
    tabsetPanel(
      tabPanel("Introduction", htmlOutput("disclaimer")),
      tabPanel("Correlation Struture in the pool", plotOutput("matPlot")),
      #tabPanel("Table", tableOutput("probs")),
      tabPanel("Report", h3(htmlOutput("Report"))))
    )
  )
)



# Define server logic required to draw a histogram ----
server <- function(input, output, session) {
  
  # Histogram of the Old Faithful Geyser Data ----
  # with requested number of bins
  # This expression that generates a histogram is wrapped in a call
  # to renderPlot to indicate that:
  #
  # 1. It is "reactive" and therefore should be automatically
  #    re-executed when inputs (input$bins) change
  # 2. Its output type is a plot
  # ------------------ App virtualenv setup (Do not edit) ------------------- #
  
  observe({
    x <- input$p
    updateSelectInput(session, "where",
                       label = "Select the participants that know someone else in the pool",
                        choices = 1:x,
                        selected = NULL
                      )
    
  })
  
  observe({
    x <- input$country
    
    # Can use character(0) to remove all choices
    if (is.null(x))
      x <- character(0)
    
    # Can also set the label and select items
    country_regions = dplyr::filter(countries, country_name == input$country)
    updateSelectInput(session, "region",
                      label = "What region do you live in?",
                      choices =  country_regions$region,
                      selected = "Main territory"
    )
  })
  
  observe({
    x <- input$where
    print(x)
    # Can use character(0) to remove all choices
    if (is.null(x) || length(x)<2){
      x <- character(0)
    } else{
    # Can also set the label and select items
    
    mat = data.frame(diag(rep(1,length(x))), row.names = x)
    names(mat) = x
    updateMatrixInput(session, "Correlation", value = as.matrix(mat))
    }
  })
  
  #output$probs <- renderDataTable({
  #output$region = renderUI({
  #  country_regions = dplyr::filter(regions, country == input$country)
  #  selectInput('region2', 'What region do you live in?', country_regions$region)
  #})
  
  
  dataInput <- reactive({
    #### Convert matrix into values
    m = matrix(0, nrow(input$Correlation),nrow(input$Correlation))
    diag(m) <- 1
    for (i in names(equivalences) ){
      if (sum(input$Correlation == i)>0){
        ind = which(input$Correlation == i, arr.ind = TRUE)
        for(u in 1:nrow(ind)){
          m[ind[u,1], ind[u,2]] = equivalences[[i]]
          m[ind[u,2], ind[u,1]] = equivalences[[i]]
        }
      }
    }
    print(m)
    return(m)
    
   })
  
  output$matPlot <- renderPlot({
    x =  dataInput()
    if (is.null(x) == FALSE){
    #### Convert matrix into values
    heatmap(x, Rowv = FALSE, Colv = FALSE, keep.dendro = FALSE)
    }
  })
  
  output$distPlot <- renderPlot({
    x =  dataInput()
    
    #hist(100, breaks = seq(from=0, to=100, by=2.5), col = "#75AADB", border = "white",
    #     xlab = "Probability (in %)",
    #     main ="Your probability distribution" )
  })
  
  output$disclaimer = renderUI({
    #tags$head(
      # Facebook OpenGraph tags - image to share when social sharing app
      #tags$meta(property = "og:title", content = share$title),
      #tags$meta(property = "og:type", content = "website"),
      #tags$meta(property = "og:url", content = share$url),
      #tags$meta(property = "og:description", content = share$description))
    
    tags$div(
      tags$p("The purpose of this calculator is to assess the efficiency of adding correlations on the expected sensitivity of the pooled test."), 
    )
  })
  
}

shinyApp(ui, server)
