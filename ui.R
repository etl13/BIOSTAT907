ui <- fluidPage(
  titlePanel("Two-Stage UMVUE and MLE Explorer (Jung & Kim, 2004)"),
  
  sidebarLayout(
    sidebarPanel(
      numericInput("n1", "Stage 1 sample size (n1):", 13, min = 1),
      numericInput("n2", "Stage 2 sample size (n2):", 30, min = 1),
      numericInput("a1", "Stage 1 futility cutoff (a1):", 3, min = 0),
      numericInput("a", "Final cutoff (a):", 12, min = 0),
      numericInput("alpha", "Alpha (1 - confidence level):", 0.05, min = 0.001, max = 0.5, step = 0.01),
      textInput("pvals", "P values (comma-separated):", "0.1,0.2,0.3,0.4,0.5"),
      actionButton("run", "Run Analysis", class = "btn-primary"),
      downloadButton("downloadData", "Download Table (CSV)"),
      helpText("Outputs UMVUE, MLE, confidence intervals, and f(m,s|p) for specified p values.")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Table", DTOutput("results_table")),
        tabPanel("PMF Plots", fluidRow(column(8, plotOutput("results_plot1",width="500px",height="300px"),
                          column(8, plotOutput("results_plot2",width="500px",height="300px")))

      )))
  )
))
    
