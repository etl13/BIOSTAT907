
library(DT)
# UI
# ========================================================
ui <- fluidPage(
  titlePanel("Two-Stage UMVUE & MLE Explorer (Updated)"),
  
  sidebarLayout(
    sidebarPanel(
      numericInput("n1", "Stage 1 sample size (n1):", 13, min=1),
      numericInput("n2", "Stage 2 sample size (n2):", 30, min=1),
      numericInput("a1", "Stage 1 futility cutoff (a1):", 3, min=0),
      numericInput("a",  "Final cutoff (a):", 12, min=0),
      numericInput("alpha","Alpha (1 - CI level):", 0.05, min=0.001, max=0.5),
      textInput("pvals", "Probability mass for true p at each observation:", "0.1,0.2,0.3,0.4,0.5"),
      numericInput("p0", "Null p0 for p-value test:", 0.2, min=0, max=1, step=0.01),
      
      checkboxInput("show_pmf", "Show PMF Columns", value=TRUE),
      
      actionButton("run", "Run Analysis", class="btn-primary"),
      downloadButton("downloadData", "Download Table (CSV)")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Estimates", DTOutput("estimates_table")),
        tabPanel("Confidence Intervals", DTOutput("ci_table")),
        tabPanel("P-Values", DTOutput("pvalue_table")),
        tabPanel("PMF Plots",
                 fluidRow(
                   column(6, plotOutput("results_plot1")),
                   column(6, plotOutput("results_plot2"))
                 ))
      )
    )
  )
)

shinyApp(ui, server) 
