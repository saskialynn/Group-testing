# Define UI for app that draws a histogram ----
library("shiny")
library("shinydashboard")
library("tidyverse")
library("ggpubr")

source("prob_laws_528.R")
source("sensitivity_fxn_shinyapp.R")
source("beta_params.R")
source("covid_case_predictions.R")
source("under_ascertainment_bias.R")

above.lod = 0.25  #### Need to update this depending on the LOD of the test


# change
COUNTRY_DATA <- read.csv(file="https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/owid-covid-data.csv", header=T, sep=",")
COUNTRY_DATA$date <- (as.Date(COUNTRY_DATA$date, "%Y-%m-%d"))
ifr_data=read.csv("country_ifr_data.csv", header=T)

COUNTRY_LIST =  intersect(unique(COUNTRY_DATA$location),ifr_data$location)



countries = read.csv("countries_names.csv")
countries_list = list()
for (c in unique(countries$country_name)){
  if (sum(countries$country_name == c ) >1 ){
    countries_list[[c]]= as.list(sapply(countries$region[countries$country_name == c], function(x){toString(x)}))
  }else{
    countries_list[[c]]= list(c("Main territory"))
  }
}
tau_alphas <- list("From Child Index Case" = 8.38, 
                   "Healthcare Setting" = 8.3, 
                   "Spouses" = 21.78,
                   "Household (Asymptomatic Index Case)" = 0.74, 
                   "Household (Symptomatic Index Case)"= 64.95)
tau_betas <- list("From Child Index Case" = 59.43, 
                  "Healthcare Setting" = 359.61, 
                  "Spouses" = 35.92,
                  "Household (Asymptomatic Index Case)" = 62.23, 
                  "Household (Symptomatic Index Case)"= 296.26)

ui <- dashboardPage(
  dashboardHeader(title = "Pooled Testing under Correlations",
                  tags$li(class = "dropdown",
                          tags$style(".main-header {max-height: 50px}"),
                          tags$style(".main-header .logo {height: 50px;}"),
                          tags$style(".sidebar-toggle {height: 50px; width=50px;padding-top: 1px !important;}"),
                          tags$style(".navbar {min-height:50px !important}"),
                          tags$style(".main-header .navbar { margin-left: 400px;}"),
                          tags$style(".main-header .logo { width: 400px;}")
                          
                  )),
  dashboardSidebar(
    actionButton("go", "Launch Calculations"),
    # sliderInput( inputId = "N",
    #              label = "What is the maximum number of people per pool?",
    #              value = 30,
    #              min=0, max=50),
     numericInput(inputId = "N",
                       label = "What is the maximum number of people per pool?",
                       value = 30,
                       min=0),
      radioButtons(inputId = "above.llod",
                  label = "What proportion of your pool do you estimate to be above the LOD of your PCR machine?",
                  choices = mat2[,2] * 100,
                  selected= 20),
     sliderInput("prev", label = "What is the estimated 95%-confidence range for the number of active COVID-19 infections per million in the community?\n (See right tab if unknown.)" , min = 0, 
                max = 100000, value = c(50, 150),
                step = 50),
      ##numericInput(inputId = "prev",
          #              label="What is the estimated prevalence of active COVID-19 infections in the community?",
          #              value=1,
          #              min=0),
          # numericInput(inputId = "prev_lw_ci",
          #              label="What is the lower 95% CI on the prevalence?",
          #              value=1,
          #              min=0),
          # numericInput(inputId = "prev_up_ci",
          #              label="What is the upper 95% CI on the prevalence?",
          #              value=1,
          #              min=0),
          radioButtons(inputId="tau_setting",
                        label="What is the network transmission scenario within the pool?",
                        choices=c("From Child Index Case", "Healthcare Setting", "Spouses",
                                 "Household (Asymptomatic Index Case)", "Household (Symptomatic Index Case)")),
    width = 300
        
  ),
  dashboardBody(
    fluidRow(
      box(
        # Title can include an icon
        title = "Results: Sensitivity of the pool",
        plotOutput("distPlot")
        ),
      box(
        side = "right",
        title = "Predicting the prevalence",
        actionButton("go_prev", "Launch local prevalence inference"),
        selectInput(inputId = "country",
                    label = "Which country do you live in?",
                    choices = unique(countries$country_name),
                    selected = "United Kingdom"),
        # selectInput(inputId = "region",
        #             label = "Which region do you live in?",
        #             choices = countries_list,
        #             selected="Main territory"),
        dateInput(inputId = "date_event",
                  label ="When will you be pooling samples?", 
                  value =paste0(Sys.Date()),
                  min = NULL,
                  max = paste0(Sys.Date() + 40) ,
                  format = "yyyy-mm-dd",
                  startview = "month",
                  weekstart = 0,
                  language = "en",
                  width = NULL,
                  autoclose = TRUE,
                  datesdisabled = NULL,
                  daysofweekdisabled = NULL),
        tags$hr(),
        plotOutput("plotprev")
      )
    ),
    fluidRow(
      tabBox(
        title = "Our method",
        # The id lets us use input$tabset1 on the server to find the current tab
        id = "tabset1",
        tabPanel("About this calculator", htmlOutput("about")),
        tabPanel("Disclaimer", htmlOutput("disclaimer"))
      )
    )
  )
)

# ui <- fluidPage(
#   
#   # App title ----
#   titlePanel("Compute the sensitivity of the pool",
#              windowTitle="Pool Testing under Correlations"),
#   
#   
#   
#   # Sidebar layout with input and output definitions ----
#   sidebarLayout(
#     
#     # Sidebar panel for inputs ----
#     sidebarPanel(
#       actionButton("go", "Launch Calculations"),
#       numericInput(inputId = "N",
#                    label = "What is the maximum number of people per pool?",
#                    value = 30,
#                    min=0),
#       radioButtons(inputId = "above.llod",
#                    label = "What is the limit of detection threshold for your PCR machine?",
#                    choices = c(35, 40, 50),
#                    selected= 35),
#       numericInput(inputId = "prev",
#                    label="What is the estimated prevalence of active COVID-19 infections in the community?",
#                    value=1, 
#                    min=0),
#       numericInput(inputId = "prev_lw_ci",
#                    label="What is the lower 95% CI on the prevalence?",
#                    value=1, 
#                    min=0),
#       numericInput(inputId = "prev_up_ci",
#                    label="What is the upper 95% CI on the prevalence?",
#                    value=1, 
#                    min=0),
#       radioButtons(inputId="tau_setting", 
#                    label="What is the network structure of the pool?", 
#                    choices=c("From Child Index Case", "Healthcare Setting", "Spouses", 
#                              "Household (Asymptomatic Index Case)", "Household (Symptomatic Index Case)")),
#     )),
#     
#     
#     
#     # Main panel for displaying outputs ----
#     # Main panel for displaying outputs ----
# 
#     mainPanel(
#       
#       # Output: Histogram ----
#       
#     #   fluidRow(
#     #     tabBox(
#     #       title = "First tabBox",
#     #       # The id lets us use input$tabset1 on the server to find the current tab
#     #       id = "tabset1", height = "250px",
#     #       tabPanel("Tab1", "First tab content"),
#     #       tabPanel("Tab2", "Tab content 2")
#     #     ),
#     #     tabBox(
#     #       side = "right", height = "250px",
#     #       selected = "Tab3",
#     #       tabPanel("Tab1", "Tab content 1"),
#     #       tabPanel("Tab2", "Tab content 2"),
#     #       tabPanel("Tab3", "Note that when side=right, the tab order is reversed.")
#     #     )
#     #   ),
#     #   fluidRow(
#     #     tabBox(
#     #       # Title can include an icon
#     #       title = tagList(shiny::icon("gear"), "tabBox status"),
#     #       tabPanel("Tab1",
#     #                "Currently selected tab from first box:",
#     #                verbatimTextOutput("tabset1Selected")
#     #       ),
#     #       tabPanel("Tab2", "Tab content 2")
#     #     )
#     #   )
#     # )
#       
#       
#       # tabsetPanel(
#       #   tabPanel("Introduction", htmlOutput("disclaimer")),
#       #   tabPanel("Correlation Struture in the pool", plotOutput("matPlot")),
#       #   tabPanel("Sensitivity", plotOutput("distPlot")),
#       #   #tabPanel("Table", tableOutput("probs")),
#       #   tabPanel("Summary of the results", h3(htmlOutput("Report"))))
#     
#   )
# )



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
  
  
  dataInputPrev <- eventReactive(input$go_prev, {
  
    bias<- compute_underascertainment_bias(min(as.Date(as.character(input$date_event)),  Sys.Date()) - (21+14),
                                           input$country, COUNTRY_DATA, 
                                           Case_to_death_delay= 21, 
                                           date_max=min(as.Date(as.character(input$date_event)) +14, Sys.Date()),
                                           plot=FALSE)
    
    bias_corr  = mean(bias$value)
    bias_corr = ifelse(is.na(bias_corr), 1, bias_corr )
    bias_sd  = sd(bias$value)
    
    
    
    ##### A BUNCH OF PRELIMINARY PARAMETERS
    MAX_DATE = as.Date(Sys.Date() -1, fmt="%Y-%m-%d")
    #PERIOD_FOR_PREDICTING = as.numeric(as.Date(input$date_event)- as.Date(max(british_prev$date, na.rm = TRUE), fmt="%Y-%m-%d"))
    PERIOD_FOR_PREDICTING = max(28, (as.numeric(as.Date(input$date_event)- MAX_DATE-1)))
    PERIOD_FOR_FITTING = 28
    
    
    if(as.Date(input$date_event) >=Sys.Date()+ 2){
      data2fit = COUNTRY_DATA %>% 
        dplyr::filter(date >= Sys.Date() - PERIOD_FOR_FITTING , location == input$country) %>% 
        dplyr::select(new_cases_smoothed_per_million)
      if (nrow(data2fit) == 0){
        data2fit = COUNTRY_DATA %>%
          dplyr::filter(date >= Sys.Date() - PERIOD_FOR_FITTING , location == "World") %>% 
          dplyr::select(new_cases_smoothed_per_million)
      }
      res_df = compute_prevalence(as.Date(as.character(input$date_event)), data2fit,
                                  country_data =COUNTRY_DATA, nb_curves=20, 
                                  distance=Difference_function, 
                                  period4predicting=PERIOD_FOR_PREDICTING + 1, period4fitting = min(PERIOD_FOR_FITTING, nrow(data2fit)), 
                                  div=2, distance_type="MSE")
      prevalence_df = res_df$res
      samples_prev = res_df$output
      
      #### Correct for all of these under-ascertainment issues
      future_prevalence_df = prevalence_df %>%
        dplyr::filter(time > 0) %>%
        dplyr::mutate(prevalence = 1/bias_corr * prevalence,
               sd_prevalence = 1/bias_corr * sd_prevalence)
    }else{
      data2fit = COUNTRY_DATA %>% 
        dplyr::filter(date >= as.Date(as.character(input$date_event)) - PERIOD_FOR_PREDICTING - 2,
                      date <= as.Date(as.character(input$date_event)), 
                      location == input$country) %>%  
        dplyr::select(new_cases_smoothed_per_million)
      #### Extrapolate if the data is 
      if (nrow(data2fit) < PERIOD_FOR_PREDICTING +1){
        data2fit= rbind(data2fit, data.frame("new_cases_smoothed_per_million" = rep(data2fit$new_cases_smoothed_per_million[nrow(data2fit)], PERIOD_FOR_PREDICTING +1 -nrow(data2fit) ) ))
      }
      ##### Compute SD using k-nn (here n=10 should suffice)
      uu = sapply(3: (PERIOD_FOR_PREDICTING),function(i){sd(data2fit$new_cases_smoothed[(i-2): (i+2)])})
      uu = c(uu, rep(uu[length(uu)], PERIOD_FOR_PREDICTING + 1  - length(uu)))
    
  
      future_prevalence_df  = data.frame("time" = 1:(PERIOD_FOR_PREDICTING+1),
                                         "prevalence" =  1/bias_corr * 1e-6 * data2fit$new_cases_smoothed_per_million[2: (PERIOD_FOR_PREDICTING+2) ],
                                         "sd_prevalence" = 1e-6/bias_corr *uu,
                                         "Date_of_cases" =seq(from=as.Date(as.character(input$date_event)) - PERIOD_FOR_PREDICTING, to = as.Date(as.character(input$date_event)), by="day"))
    }
    
    future_prevalence_df = future_prevalence_df %>% 
      dplyr::mutate(ymin = prevalence -2*sd_prevalence,ymax = prevalence +2*sd_prevalence)
    future_prevalence_df["ymin"] = sapply(future_prevalence_df["ymin"], function(x){ifelse(x<0,1e-6,x)})
    future_prevalence_df["ymax"] = sapply(future_prevalence_df["ymax"], function(x){ifelse(x>1,1,x)})
    return(list(prevalence_df=future_prevalence_df))
  })
  
  
  
  
  
  dataInput <- eventReactive(input$go, {
    
    #### Step 0: load parameters from inputs
    B=100
    pool.max = input$N ### size of the group
    prev = mean(input$prev)/1e6
    params  = beta.parms.from.quantiles(q=c(input$prev[1]/1e6, input$prev[2]/1e6))
    prev_alpha = params$a
    prev_beta = params$b
    tau_setting = input$tau_setting
    tau_alpha = tau_alphas[[tau_setting]]
    tau_beta = tau_betas[[tau_setting]]
    tau = tau_alpha/(tau_alpha + tau_beta)
    #### Step 1: 
    #### Do the simulations to compute the probabilities
    withProgress(message = paste0('Running simulations across different pool sizes'), max = 1,value = 0, {
      it=1
      #for (g in max(1, G):input$N){
      for (N in 1:input$N){
        #dist<-simulate_infections(G, g - G, m, prev, B=1000)
        a = list(proba_laws(N, prev, prev_alpha, prev_beta, 
                            tau, tau_alpha, tau_beta,
                            tau_graph_effect=NULL,
                            prev_graph_effect=NULL,
                            null_mod = NULL,
                            null_mod_prev_graph=NULL,
                            B=B),
                 proba_laws(N, prev, prev_alpha, prev_beta, 
                            tau, tau_alpha, tau_beta,
                            tau_graph_effect=tau_graph_effect,
                            prev_graph_effect=NULL,
                            null_mod = NULL,
                            null_mod_prev_graph=NULL,
                            B=B),
                 proba_laws(N, prev, prev_alpha, prev_beta, 
                            tau, tau_alpha, tau_beta,
                            tau_graph_effect=NULL,
                            prev_graph_effect=prev_graph_effect,
                            null_mod = NULL,
                            null_mod_prev_graph=NULL,
                            B=B),
                 proba_laws(N, prev, prev_alpha, prev_beta, 
                            tau, tau_alpha, tau_beta,
                            tau_graph_effect=tau_graph_effect,
                            prev_graph_effect=prev_graph_effect,
                            null_mod = NULL,
                            null_mod_prev_graph=NULL,
                            B=B), 
                 proba_laws(N, prev, prev_alpha, prev_beta, 
                            tau, tau_alpha, tau_beta,
                            tau_graph_effect=NULL,
                            prev_graph_effect=NULL,
                            null_mod = TRUE,
                            null_mod_prev_graph=NULL,
                            B=B),
                 proba_laws(N, prev, prev_alpha, prev_beta, 
                            tau, tau_alpha, tau_beta,
                            tau_graph_effect=NULL,
                            prev_graph_effect=NULL,
                            null_mod = NULL,
                            null_mod_prev_graph=TRUE,
                            B=B)
        )
        
        names = c("Fixed Prevalence/Fixed Network Transmission", "Fixed Prevalence/Random Network Transmission",
                  "Random Prevalence/Fixed Network Transmission","Random Prevalence/Random Network Transmission", 
                  "Null Model (No Correlations, Fixed Prevalence)", "Null Model (No Correlations, Random Prevalence)")
        
        #incProgress(1/(input$p - max(1, G) + 1), detail = paste("*"))
        
        for (n in 1:length(names)){
          sens <- data.frame(sens_fxn(N, above.llod=as.numeric(input$above.llod))) # prob(test positive | sum(Y_i) = k)
          sens <- sens$mean 
          sens_indiv <- data.frame(sens_fxn(1, above.llod=as.numeric(input$above.llod)))$mean
          prob_pos <- sapply(1:B, function(b){sum(a[[n]][2:(N+1), b] * sens)}) 
          # a[[n]][2:(N+1), b] is the vector of probabilities (for the bth simulation)
          # of having 1:N positives in a pool of size N
          
          if (it ==1){
            res = data.frame(
              sensit = sapply(1:B, function(b){
                sum(a[[n]][2:(N+1), b] * sens)/sum(a[[n]][2:(N+1), b])}), #sensitivity
              
              ppa = sapply(1:B, function(b){
                (sum(a[[n]][2:(N+1), b] * sens)/sum(a[[n]][2:(N+1), b]))/sens_indiv}),
              
              num_tests = (1/N) + prob_pos,
              
              caught_cases = sapply(1:B, function(b){
                (sum(c(1:N)*a[[n]][2:(N+1), b] * sens)/sum(a[[n]][2:(N+1), b]))})/((1/N) + prob_pos),
              
              missed_cases_persample = sapply(1:B, function(b){ # missed cases per sample
                (sum(c(1:N)*a[[n]][2:(N+1), b] * (1-sens)))})/((1/N) + prob_pos),
              
              type = names[n], 
              pool_size = N, 
              prev = prev,
              prev_param = paste(prev_alpha, ",", prev_beta),
              #prev_setting = pi_dist$state_date[p],
              tau = tau,
              tau_param = paste(tau_alpha, ",", tau_beta), 
              tau_setting = tau_setting)
          }else{
            res = rbind(res,
                        data.frame(sensit = sapply(1:B, function(b){
                          sum(a[[n]][2:(N+1), b] * sens)/sum(a[[n]][2:(N+1), b])}),
                          
                          ppa = sapply(1:B, function(b){
                            (sum(a[[n]][2:(N+1), b] * sens)/sum(a[[n]][2:(N+1), b]))/sens_indiv}),
                          
                          num_tests = (1/N) + prob_pos,
                          
                          caught_cases = sapply(1:B, function(b){
                            (sum(c(1:N)*a[[n]][2:(N+1), b] * sens)/sum(a[[n]][2:(N+1), b]))})/((1/N) + prob_pos),
                          
                          missed_cases_persample = sapply(1:B, function(b){ # missed cases per sample
                            (sum(c(1:N)*a[[n]][2:(N+1), b] * (1-sens)))})/((1/N) + prob_pos),
                          type = names[n], 
                          pool_size = N, 
                          prev = prev,
                          prev_param = paste(prev_alpha, ",", prev_beta),
                          #prev_setting = pi_dist$state_date[p],
                          tau = tau,
                          tau_param = paste(tau_alpha, ",", tau_beta), 
                          tau_setting = tau_setting))
          }
          
          it = it + 1
          
        }
        incProgress(1/(input$N+2), detail = paste0("Simulation for pool size ", N + 1, "/", input$N))
      }
   
    incProgress(1 /(input$N+2), detail = "Summarizing results")
    conf_int <- res %>% 
      group_by(type, pool_size, tau, prev, tau_param, prev_param, tau_setting) %>% #prev_setting, 
      summarise(sd_sens= sd(sensit), 
                sd_ppa = sd(ppa), 
                sd_tests = sd(num_tests), 
                sd_missed = sd(missed_cases_persample),
                mean_sens= mean(sensit), 
                mean_ppa = mean(ppa), 
                mean_tests = mean(num_tests), 
                mean_missed = mean(missed_cases_persample), 
                sens_q025 = quantile(sensit, probs=0.025), 
                sens_q975 = quantile(sensit, probs = 0.975), 
                ppa_q025 = quantile(ppa, probs=0.025), 
                ppa_q975 = quantile(ppa, probs = 0.975), 
                tests_q025 = quantile(num_tests, probs=0.025), 
                tests_q975 = quantile(num_tests, probs = 0.975), 
                missed_q025 = quantile(missed_cases_persample, probs=0.025), 
                missed_q975 = quantile(missed_cases_persample, probs = 0.975), 
                sens_fda = mean(sensit >= 0.85), 
                sens_fda_q025 = mean(sens_q025 >=0.85), 
                sens_fda_q975 = mean(sens_q975>=0.85))
    conf_int$prev <- round(conf_int$prev, 5)
    conf_int$tau <- round(conf_int$tau, 5)
    })
    #return(list(m=m, x=res, pool = input$p))
    return(list(x=conf_int))
  })
  

  
  
  output$distPlot <- renderPlot({
    conf_int =  dataInput()$x
    indiv_sens <- data.frame(sens_fxn(1, above.llod = as.numeric(input$above.llod) ))$mean
    #conf_int = data.frame(conf_int)
    # labels and factors for plotting
    types <- factor(c("Null Model (No Correlations, Fixed Prevalence)", "Null Model (No Correlations, Random Prevalence)", 
                      "Random Prevalence/Random Network Transmission", 
                      "Fixed Prevalence/Fixed Network Transmission", 
                      "Random Prevalence/Fixed Network Transmission",
                      "Fixed Prevalence/Random Network Transmission"), 
                    levels = c("Null Model (No Correlations, Fixed Prevalence)",
                               "Null Model (No Correlations, Random Prevalence)",
                               "Random Prevalence/Random Network Transmission", 
                               "Fixed Prevalence/Fixed Network Transmission",
                               "Random Prevalence/Fixed Network Transmission",
                               "Fixed Prevalence/Random Network Transmission"))
    conf_int$type <- factor(conf_int$type, 
                            levels = c("Null Model (No Correlations, Fixed Prevalence)",
                                       "Null Model (No Correlations, Random Prevalence)",
                                       "Random Prevalence/Random Network Transmission", 
                                       "Fixed Prevalence/Fixed Network Transmission",
                                       "Random Prevalence/Fixed Network Transmission",
                                       "Fixed Prevalence/Random Network Transmission"))
    lines <- c("Network Transmission" = "solid", "No Correlation" = "dashed")
    lines_sensplots <- c("Network Transmission" = "solid", "No Correlation" = "dashed", 
                         "Individual Testing" = "dotdash", "FDA Threshold" = "dotted")
    
    conf_int$prev_labels <- paste(100*conf_int$prev, "%")
    #require(graphics)
    conf_int$prev_labels <- factor(conf_int$prev_labels)
    
    # clean up tau_setting column
    conf_int$tau_labels <- factor(conf_int$tau_setting, 
                                  levels = c("From Child Index Case", 
                                             "Healthcare Setting", 
                                             "Spouses", 
                                             "Household (Asymptomatic Index Case)", 
                                             "Household (Symptomatic Index Case)"),
                                  labels = c(expression("From Child Index Case"),
                                             expression("Healthcare Setting"),
                                             expression("Spouses"),
                                             expression("Household \n (Asymptomatic Index Case)"), 
                                             expression("Household \n (Symptomatic Index Case)")
                                             ))
    # Fixed and Tau Graph Effect
    sens_plots_pifixed <- lapply(types[c(4, 6)], function(i){ # Fixed and Tau Graph Effect
      ggplot(conf_int %>% filter(type == i), aes(x=pool_size, y= mean_sens)) +
        geom_line(aes(color=tau_labels, linetype = "Network Transmission")) +
        geom_line(data=conf_int[conf_int$type == "Null Model (No Correlations, Fixed Prevalence)",], aes(linetype = "No Correlation"))+
        geom_hline(aes(yintercept=indiv_sens, linetype = "Individual Testing"), color = "grey35")+
        geom_hline(aes(yintercept=0.85, linetype = "FDA Threshold"), color = "grey35")+
        facet_wrap(~ prev_labels) +
        geom_ribbon(data = . %>% filter(type==i), 
                    aes(x=pool_size, ymin = sens_q025, ymax = sens_q975, 
                        fill = tau_labels), show.legend = F, alpha = 0.2, colour = NA)+
        scale_linetype_manual(values = lines_sensplots) +
        theme_bw() + 
        labs(x = "Pool Size", y = "Sensitivity", linetype = lines) +
        ggtitle(i) +
        guides(color=guide_legend(title=expression(paste("Transmission (", tau, ")"))), linetype=guide_legend(title="Model Type")) 
    })
    #sens_plots_pifixed
    
    # Pi Graph Effect and All Graph Effect
    sens_plots_pigraph <- lapply(types[c(3, 5)], function(i){
      ggplot(conf_int %>% filter(type == i), aes(x=pool_size, y= mean_sens)) +
        geom_line(aes(color=tau_labels, linetype = "Network Transmission")) +
        geom_line(data=conf_int[conf_int$type == "Null Model (No Correlations, Random Prevalence)",], aes(linetype = "No Correlation"))+
        geom_hline(aes(yintercept=indiv_sens, linetype = "Individual Testing"), color = "grey35")+
        geom_hline(aes(yintercept=0.85, linetype = "FDA Threshold"), color = "grey35")+
        facet_wrap(~ prev_labels) +
        geom_ribbon(data = . %>% filter(type==i), aes(x=pool_size, ymin = sens_q025, ymax = sens_q975, 
                                                      fill = tau_labels), show.legend = F, alpha = 0.2, colour = NA)+
        geom_ribbon(data = conf_int[conf_int$type == "Null Model (No Correlations, Random Prevalence)",], 
                    aes(x=pool_size, ymin = sens_q025, ymax = sens_q975, 
                        fill = tau_labels), show.legend = F, alpha = 0.07, colour = NA)+
        scale_linetype_manual(values = lines_sensplots) +
        theme_bw() + 
        labs(x = "Pool Size", y = "Sensitivity", linetype = lines) +
        ggtitle(i) +
        guides(color=guide_legend(title=expression(paste("Transmission (", tau, ")"))), linetype=guide_legend(title="Model Type")) 
    })
    #sens_plots_pigraph
    ggarrange(sens_plots_pifixed[[1]], sens_plots_pifixed[[2]],
              sens_plots_pigraph[[1]], sens_plots_pigraph[[2]], 
              nrow = 2, ncol =2, 
              common.legend = TRUE, legend = "bottom")
  })
  
  output$plotprev =renderPlot({
    x =  dataInputPrev()
    dd = x$prevalence_df %>% 
      dplyr::filter(Date_of_cases == max(x$prevalence_df$Date_of_cases))
    ggplot(x$prevalence_df)+
      geom_line(aes(x=Date_of_cases, y=1e6 * prevalence), colour="red")+
      geom_ribbon(aes(x=Date_of_cases, ymin=1e6 * ymin, ymax=1e6*ymax), colour="grey", alpha=0.5)+
      theme_bw()+ xlab("Date") + ylab("Incidence (per Million)")+
      ggtitle(paste0("Predicted Incidence (per Million) Per Region\n 95% Confidence Region for Day of Pooling:\n [",
                     round(1e6 * dd$ymin,0), ";", round(1e6*dd$ymax,0), "]" )
              )+scale_y_log10() + 
      theme(plot.title = element_text(hjust = 0.5,
                                      face = "bold",
                                      size = 14))
    
  })
  
  output$disclaimer = renderUI({
    #tags$head(
    # Facebook OpenGraph tags - image to share when social sharing app
    #tags$meta(property = "og:title", content = share$title),
    #tags$meta(property = "og:type", content = "website"),
    #tags$meta(property = "og:url", content = share$url),
    #tags$meta(property = "og:description", content = share$description))
    
    tags$div(
      tags$p("The purpose of this calculator is to assess the efficiency of adding correlations on the expected sensitivity of the pooled sample.")
    )
  })
  
  output$about = renderUI({
    #tags$head(
    # Facebook OpenGraph tags - image to share when social sharing app
    #tags$meta(property = "og:title", content = share$title),
    #tags$meta(property = "og:type", content = "website"),
    #tags$meta(property = "og:url", content = share$url),
    #tags$meta(property = "og:description", content = share$description))
    
    tags$div(
      tags$p("The purpose of this calculator is to assess the efficiency of adding correlations on the expected sensitivity of the pooled sample.")
    )
  })
  
}

shinyApp(ui, server)
