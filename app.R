# Define UI for app that draws a histogram ----
library("shiny")
library("tidyverse")
library("shinyMatrix")

source("preprocessing.R")
source("MCMC.R")
above.lod = 0.05  #### Need to update this depending on the LOD of the test
prev_tot=0.01  ### Need to compute this based on the date and geographical location

v0 = c(1,2,3,4,1,1,1,2,3,1,2,3,1,1,1,2,3,1,3,1,2,3,1)
countries = read.csv("countries_names.csv")
countries_list = list()
for (c in unique(countries$country_name)){
  if (sum(countries$country_name == c ) >1 ){
    countries_list[[c]]= as.list(sapply(countries$region[countries$country_name == c], function(x){toString(x)}))
  }else{
    countries_list[[c]]= list(c("Main territory"))
  }
}

equivalences <- list("H" = 0.3,
                     "W" = 0.1,
                     "S" =	0.4,
                     "C" = 0.12)
data = read_csv( "lookup_sensitivity.csv")
ui <- fluidPage(
  
  # App title ----
  titlePanel("Compute the sensitivity of the pool",
             windowTitle="Pool Testing under Correlations"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      numericInput(inputId = "p",
                   label = "What is the maximum number of people per pool?",
                   value = 20,
                   min=0),
      radioButtons(inputId = "above.llod",
                   label = "What is the limit of detection threshold for your PCR machine?",
                   choices = c(35, 40, 50),
                   selected= 35),
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
        choices = 1:20,
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
      tabPanel("Sensitivity", plotOutput("distPlot")),
      #tabPanel("Table", tableOutput("probs")),
      tabPanel("Summary of the results", h3(htmlOutput("Report"))))
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
    if (is.null(input$where) == FALSE){
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
      v = v0[1:length(input$where)]
    }else{
      m= NULL
      v = c()
    }
    
    #### return the appropriate entry in the data table
    prev= rep(prev_tot,input$p)
    data_temp = data %>% dplyr::filter((above.llod == 0.05) & (pool >= G))
    data_temp$probability = 0
    data_temp$probability_null = 0
    #### Do the simulations to compute the probabilities
    for (g in max(2, length(v)):input$p){
      if((g-length(v))>0){
        proba = 1-exp(sum(log(1-v*prev_tot)) + sum(log(1-prev[1:(g-length(v))])))
        dist<-simulate_infections(v, m, prev[1:(g-length(v))], B=1000)
      }else{
        print("Esta quai")
        proba = 1-exp(sum(log(1-v*prev_tot)) )
        dist<-simulate_infections(v, m, c(), B=1000)

      }
      dist$Freq <- (proba) * dist$Freq
      dist = rbind(data.frame(temp=0, Freq = 1-proba), dist)
      dist$Freq = dist$Freq/sum(dist$Freq)
      data_temp[which(data_temp$pool == g), 'probability'] = unlist(apply(data_temp[which(data_temp$pool == g), ],
                                                                          1,
                                                                          function(x){as.numeric(dist$Freq[which(dist$temp == x['pos'])])}))
      data_temp[which(data_temp$pool == g), 'probability_null'] = sapply(0:g, function(x){dbinom(x,g, prev_tot)})
    }
    

    #### Compute the sensitivity
    group_by(data_temp, pool, above.llod, limit) %>% 
      summarize(pos=weighted.mean(pos, w=probability),
                total.tests1=weighted.mean(total.tests, w=probability), 
                tests.per.sample1=weighted.mean(tests.per.sample, w=probability)/mean(pool), # calculate tests per sample
                tn1=weighted.mean(tn, w=probability), 
                tp1=weighted.mean(tp, w=probability), 
                fn1=weighted.mean(fn, w=probability), 
                fp1=weighted.mean(fp, w=probability),
                total.tests0=weighted.mean(total.tests, w=probability_null), 
                tests.per.sample0=weighted.mean(tests.per.sample, w=probability_null)/mean(pool), # calculate tests per sample
                tn0=weighted.mean(tn, w=probability_null), 
                tp0=weighted.mean(tp, w=probability_null), 
                fn0=weighted.mean(fn, w=probability_null), 
                fp0=weighted.mean(fp, w=probability_null)) -> apw
    
    return(list(m=m, x=apw, pool = input$p))
   })
  
  output$matPlot <- renderPlot({
    print("data input m")
    print(dataInput()$m)
    x =  dataInput()$m
    if (is.null(x) == FALSE || (dataInput()$pool == 1)){
    #### Convert matrix into values
    heatmap(x, Rowv = FALSE, Colv = FALSE, keep.dendro = FALSE)
    }
  })
  
  output$distPlot <- renderPlot({
    x =  dataInput()$x
    
    ggplot(x) +
      geom_line(aes(x=pool, y=tp1/(tp1 + fn1), colour="red"),size=2) + 
      geom_line(aes(x=pool, y=tp0/(tp0 + fn0), colour="black"),linetype = "dashed", size=1) + 
      theme_bw()
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
