# Define UI for app that draws a histogram ----
library("shiny")
library("tidyverse")
library("shinyMatrix")
library("ggpubr")

#source("preprocessing.R") # change this?
#source("MCMC.R")
source("prob_laws_528.R")
source("sensitivity_fxn_shinyapp.R")
source("beta_params.R")
above.lod = 0.25  #### Need to update this depending on the LOD of the test


# change
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

ui <- fluidPage(
  
  # App title ----
  titlePanel("Compute the sensitivity of the pool",
             windowTitle="Pool Testing under Correlations"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      actionButton("go", "Launch Calculations"),
      numericInput(inputId = "N",
                   label = "What is the maximum number of people per pool?",
                   value = 30,
                   min=0),
      radioButtons(inputId = "above.llod",
                   label = "What is the limit of detection threshold for your PCR machine?",
                   choices = c(35, 40, 50),
                   selected= 35),
      numericInput(inputId = "prev",
                   label="What is the estimated prevalence of active COVID-19 infections in the community?",
                   value=1, 
                   min=0),
      numericInput(inputId = "prev_lw_ci",
                   label="What is the lower 95% CI on the prevalence?",
                   value=1, 
                   min=0),
      numericInput(inputId = "prev_up_ci",
                   label="What is the upper 95% CI on the prevalence?",
                   value=1, 
                   min=0),
      radioButtons(inputId="tau_setting", 
                   label="What is the network structure of the pool?", 
                   choices=c("From Child Index Case", "Healthcare Setting", "Spouses", 
                             "Household (Asymptomatic Index Case)", "Household (Symptomatic Index Case)")),
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
  
  dataInput <- eventReactive(input$go, {
    
    #### Step 0: load parameters from inputs
    B=100
    pool.max = input$N ### size of the group
    prev = input$prev
    prev_alpha = beta.parms.from.quantiles(q=c(input$prev_lw_ci, input$prev_up_ci))$a
    prev_beta = beta.parms.from.quantiles(q=c(input$prev_lw_ci, input$prev_up_ci))$b
    tau_setting = input$tau_setting
    tau_alpha = tau_alphas[[tau_setting]]
    tau_beta = tau_betas[[tau_setting]]
    tau = tau_alpha/(tau_alpha + tau_beta)
    #### Step 1: 
    #### Do the simulations to compute the probabilities
    withProgress(message = paste0('Running simulations across different pool sizes'), value = 0, {
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
        
        names = c("Fixed", "Tau Graph Effect","Pi Graph Effect","All Graph Effect", 
                  "Null Model (Fixed)", "Null Model (Pi Graph Effect)")
        
        #incProgress(1/(input$p - max(1, G) + 1), detail = paste("*"))
        
        for (n in 1:length(names)){
          sens <- data.frame(sens_fxn(N)) # prob(test positive | sum(Y_i) = k)
          sens <- sens$mean 
          sens_indiv <- data.frame(sens_fxn(1))$mean
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
      }
    })
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
    conf_int$prev <- round(conf_int$prev, 3)
    conf_int$tau <- round(conf_int$tau, 3)
    #return(list(m=m, x=res, pool = input$p))
    return(list(x=conf_int))
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
    conf_int =  dataInput()$x
    indiv_sens <- data.frame(sens_fxn(1))$mean
    #conf_int = data.frame(conf_int)
    # labels and factors for plotting
    types <- factor(c("Null Model (Fixed)", "Null Model (Pi Graph Effect)", 
                      "All Graph Effect", "Fixed", "Pi Graph Effect", "Tau Graph Effect"), 
                    levels = c("Null Model (Fixed)","Null Model (Pi Graph Effect)",
                               "All Graph Effect", "Fixed", "Pi Graph Effect", "Tau Graph Effect"))
    conf_int$type <- factor(conf_int$type, 
                            levels = c("Null Model (Fixed)","Null Model (Pi Graph Effect)",
                                       "All Graph Effect", "Fixed", "Pi Graph Effect", "Tau Graph Effect"))
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
        geom_line(data=conf_int[conf_int$type == "Null Model (Fixed)",], aes(linetype = "No Correlation"))+
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
        geom_line(data=conf_int[conf_int$type == "Null Model (Pi Graph Effect)",], aes(linetype = "No Correlation"))+
        geom_hline(aes(yintercept=indiv_sens, linetype = "Individual Testing"), color = "grey35")+
        geom_hline(aes(yintercept=0.85, linetype = "FDA Threshold"), color = "grey35")+
        facet_wrap(~ prev_labels) +
        geom_ribbon(data = . %>% filter(type==i), aes(x=pool_size, ymin = sens_q025, ymax = sens_q975, 
                                                      fill = tau_labels), show.legend = F, alpha = 0.2, colour = NA)+
        geom_ribbon(data = conf_int[conf_int$type == "Null Model (Pi Graph Effect)",], 
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
