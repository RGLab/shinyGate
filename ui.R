library(rCharts)

data <- list.dirs("gs", full.names=TRUE)[-1]
shinyUI( pageWithSidebar(
  
  headerPanel("iGate"),
  
  sidebarPanel(
    
    ## these are generated through JavaScript; we just need placeholders here
    uiOutput("gating_controls"),
    uiOutput("sample_controls"),
    uiOutput("apply_controls"),
    
    selectInput("data_selection", choices=data, label='Select Data', selected="gs/LoveLab"),
    selectInput("sample_selection", choices='', label='Select Samples to Work With', multiple=TRUE),
    actionButton("select_all_samples", label="Select All Samples"),
    selectInput("parent_selection", choices="root", label='Select Parent Population'),
    selectInput("channel_1", choices="Live", label='Select Channel 1'),
    selectInput("channel_2", choices="Dead", label='Select Channel 2'),
    selectInput("gating_method",
      choices=list(
        `Minimum Density`="mindensity", 
        Cytokine="cytokine",
        flowClust.1d="flowClust.1d",
        flowClust.2d="flowClust.2d"
      ),
      label="Gating Method"
    ),
    conditionalPanel("input.gating_method != 'flowClust.2d'",
      checkboxInput("gate_positive_and_negative", "Gate both positive and negative sides?"),
      radioButtons("select_channels", 
        "Select Channels to be used in Gating", 
        choices=list(
          "Channel 1"="channel_1", 
          "Channel 2"="channel_2"
        )
      )
    ),
    textInput("gate_alias", label="Gate Alias", value=''),
    actionButton("apply_gate", "Apply Gate"),
    br(),
    br(),
    tags$div( style="overflow: auto;",
      downloadButton("save_gs", "Save GatingSet"),
      checkboxInput("autosave", "Autosave?", TRUE)
    )
  ),
  
  mainPanel(
    tabsetPanel(
      tabPanel("1D Density Plots",
        numericInput("1d_densityplot_n", "Maximum Number of Cells to Plot", value=10000),
        plotOutput("1d_densityplot"),
        htmlOutput("1d_densityplot_diag")
      ),
      tabPanel("2D Density Plots",
        tags$div( style="overflow: auto;",
          selectInput("2d_overlay", "Overlay", ""),
          numericInput("2d_densityplot_n", "Maximum Number of Cells to Plot", value=10000)
          
        ),
        plotOutput("2d_densityplot"),
        htmlOutput("2d_densityplot_diag")
        
      ),
      tabPanel("Gate Plots",
        plotOutput("gate_plot"),
        actionButton("remove_gate", label="Remove Gates"),
        uiOutput("gate_plot_diag")
      ),
      tabPanel("Gating Hierarchy",
        plotOutput("gating_hierarchy")
      ),
      tabPanel("Debug",
        tags$div(
          textInput("R_input", "Enter an R Command", ''),
          actionButton("R_send", label="Send Command")
        ),
        tags$div(
          verbatimTextOutput("debug")
        )
      )
    ),
    
    includeCSS("www/css/jquery-ui-1.10.3.custom.css"),
    includeCSS("www/css/styles.css"),
    
    includeScript("www/js/jquery-1.9.1.js"),
    includeScript("www/js/jquery-ui-1.10.3.custom.min.js"),
    includeScript("www/js/jquery.dialogextend.min.js"),
    includeScript("www/js/gating_controls.js"),
    includeScript("www/js/apply_gate.js"),
    includeScript("www/js/R_input.js")
    
  )
  
) )
