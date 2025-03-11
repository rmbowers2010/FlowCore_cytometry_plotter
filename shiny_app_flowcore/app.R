options(shiny.maxRequestSize = 10000 * 1024^2)

library(shiny)
library(flowCore)
library(ggplot2)
library(rlang)
library(tools)      # for file_path_sans_ext()
library(viridis)    # for color scales
library(MASS)       # for simulation (if used)

# Helper function: reads the FCS file from the given path and returns a data frame
read_fcs_df <- function(fp) {
  message("Inside helper, using file path: ", fp)
  fcs_obj <- read.FCS(fp, transformation = FALSE)
  as.data.frame(flowCore::exprs(fcs_obj))
}

ui <- fluidPage(
  titlePanel("Professional Flow Cytometry Data Viewer"),
  sidebarLayout(
    sidebarPanel(
      checkboxInput("simulate", "Simulate Data", value = FALSE),
      conditionalPanel(
        condition = "input.simulate == true",
        sliderInput("numPopulations", "Number of Populations", min = 1, max = 4, value = 1, step = 1),
        actionButton("simulateButton", "Generate Simulated Data")
      ),
      conditionalPanel(
        condition = "input.simulate == false",
        fileInput("fcs_file", "Choose FCS File", accept = ".fcs")
      ),
      radioButtons("plot_type", "Plot Type",
                   choices = c("Histogram", "2D Plot"), inline = TRUE),
      conditionalPanel(
        condition = "input.plot_type == 'Histogram'",
        uiOutput("channel_select")
      ),
      conditionalPanel(
        condition = "input.plot_type == '2D Plot'",
        uiOutput("x_channel_select"),
        uiOutput("x_slider"),
        uiOutput("y_channel_select"),
        uiOutput("y_slider"),
        radioButtons("plot2d_style", "2D Plot Style",
                     choices = c("Scatter", "Density", "Contour"), inline = TRUE)
      ),
      br(),
      downloadButton("downloadPlot", "Download Plot (PDF)"),
      br(), br(),
      downloadButton("downloadSummary", "Download Overall Summary (CSV)"),
      br(), br(),
      conditionalPanel(
        condition = "input.plot_type == '2D Plot'",
        downloadButton("downloadGated", "Download Gated Summary (CSV)")
      )
    ),
    mainPanel(
      plotOutput("plot", brush = brushOpts(id = "plot_brush")),
      br(),
      h3("Overall Data Summary"),
      tableOutput("summary_table"),
      conditionalPanel(
        condition = "input.plot_type == '2D Plot'",
        h3("Gated Data Summary"),
        tableOutput("gated_summary_table")
      )
    )
  )
)

server <- function(input, output, session) {
  
  fileName <- reactive({
    if (isTRUE(input$simulate)) {
      "SimulatedData"
    } else {
      req(input$fcs_file)
      file_path_sans_ext(input$fcs_file$name)
    }
  })
  
  # Use withProgress() to wrap the file reading/processing; also check file size
  dataFCS <- reactive({
    if (isTRUE(input$simulate)) {
      req(simData())
      simData()
    } else {
      req(input$fcs_file)
      # Check file size and warn if it's very large
      if (input$fcs_file$size > 100 * 1024^2) {
        showNotification("Warning: Large file upload may be slow on mobile devices.", type = "warning")
      }
      withProgress(message = "Uploading and processing file...", value = 0, {
        fp <- as.character(input$fcs_file$datapath)
        incProgress(0.3, detail = "Reading file...")
        dfFCS <- tryCatch({
          read_fcs_df(fp)
        }, error = function(e) {
          message("Error reading FCS: ", e)
          return(NULL)
        })
        req(!is.null(dfFCS))
        incProgress(0.6, detail = "Processing data...")
        colnames(dfFCS) <- make.names(as.character(colnames(dfFCS)), unique = TRUE)
        incProgress(1, detail = "Complete")
        dfFCS
      })
    }
  })
  
  simData <- eventReactive(input$simulateButton, {
    req(input$numPopulations)
    nPop <- input$numPopulations
    dataList <- lapply(1:nPop, function(i) {
      mu <- c(runif(1, 0, 100), runif(1, 0, 100))
      sigma <- matrix(c(runif(1, 5, 20), 0, 0, runif(1, 5, 20)), nrow = 2)
      MASS::mvrnorm(100, mu = mu, Sigma = sigma)
    })
    simMatrix <- do.call(rbind, dataList)
    dfSim <- as.data.frame(simMatrix)
    names(dfSim) <- c("SimX", "SimY")
    dfSim
  })
  
  output$channel_select <- renderUI({
    req(dataFCS())
    selectInput("channel", "Select Channel", choices = names(dataFCS()))
  })
  output$x_channel_select <- renderUI({
    req(dataFCS())
    selectInput("x_channel", "Select X Channel", choices = names(dataFCS()))
  })
  output$y_channel_select <- renderUI({
    req(dataFCS())
    selectInput("y_channel", "Select Y Channel", choices = names(dataFCS()))
  })
  
  output$x_slider <- renderUI({
    req(dataFCS(), input$x_channel)
    x_data <- dataFCS()[[input$x_channel]]
    sliderInput("x_range", "X Range",
                min = min(x_data, na.rm = TRUE),
                max = max(x_data, na.rm = TRUE),
                value = c(min(x_data, na.rm = TRUE), max(x_data, na.rm = TRUE)),
                step = (max(x_data, na.rm = TRUE) - min(x_data, na.rm = TRUE)) / 100)
  })
  output$y_slider <- renderUI({
    req(dataFCS(), input$y_channel)
    y_data <- dataFCS()[[input$y_channel]]
    sliderInput("y_range", "Y Range",
                min = min(y_data, na.rm = TRUE),
                max = max(y_data, na.rm = TRUE),
                value = c(min(y_data, na.rm = TRUE), max(y_data, na.rm = TRUE)),
                step = (max(y_data, na.rm = TRUE) - min(y_data, na.rm = TRUE)) / 100)
  })
  
  currentPlot <- reactive({
    req(dataFCS())
    df <- dataFCS()
    title_text <- paste(fileName(), "-", input$plot_type)
    
    if (input$plot_type == "Histogram") {
      req(input$channel)
      x_data <- df[[input$channel]]
      m <- mean(x_data, na.rm = TRUE)
      med <- median(x_data, na.rm = TRUE)
      total_count <- nrow(df)
      ggplot(df, aes(x = !!sym(input$channel))) +
        geom_histogram(bins = 30, fill = "blue", color = "black", alpha = 0.7) +
        geom_vline(xintercept = m, color = "red", linetype = "dashed") +
        geom_vline(xintercept = med, color = "green", linetype = "dotted") +
        annotate("text", x = Inf, y = Inf, label = paste("Total =", total_count),
                 hjust = 1.1, vjust = 1.1, size = 5, color = "black") +
        annotate("text", x = m, y = Inf, label = paste("Mean =", round(m, 2)),
                 vjust = -0.5, color = "red", size = 5) +
        annotate("text", x = med, y = Inf, label = paste("Median =", round(med, 2)),
                 vjust = -1.5, color = "green", size = 5) +
        labs(x = paste("Channel:", input$channel),
             y = "Count",
             title = title_text) +
        theme(axis.title = element_text(size = 16),
              axis.text = element_text(size = 14),
              plot.title = element_text(size = 18, face = "bold"))
    } else if (input$plot_type == "2D Plot") {
      req(input$x_channel, input$y_channel, input$x_range, input$y_range, input$plot2d_style)
      filtered_df <- df[df[[input$x_channel]] >= input$x_range[1] &
                          df[[input$x_channel]] <= input$x_range[2] &
                          df[[input$y_channel]] >= input$y_range[1] &
                          df[[input$y_channel]] <= input$y_range[2], ]
      total_count <- nrow(filtered_df)
      
      if (input$plot2d_style == "Scatter") {
        if(nrow(filtered_df) > 1) {
          corr <- round(cor(filtered_df[[input$x_channel]], filtered_df[[input$y_channel]], use = "complete.obs"), 2)
        } else { 
          corr <- NA 
        }
        brushed <- brushedPoints(filtered_df, input$plot_brush)
        gated_count <- ifelse(is.null(brushed), 0, nrow(brushed))
        ggplot(filtered_df, aes(x = !!sym(input$x_channel), y = !!sym(input$y_channel))) +
          geom_point(color = "grey") +
          geom_point(data = brushed, color = "red", size = 2) +
          annotate("text", x = Inf, y = Inf, 
                   label = paste("Total =", total_count, "\nGated =", gated_count),
                   hjust = 1.1, vjust = 1.1, size = 5, color = "black") +
          labs(x = paste("X Channel:", input$x_channel),
               y = paste("Y Channel:", input$y_channel),
               title = paste(title_text, "\nCorrelation =", corr)) +
          theme(axis.title = element_text(size = 16),
                axis.text = element_text(size = 14),
                plot.title = element_text(size = 18, face = "bold"))
      } else if (input$plot2d_style == "Density") {
        ggplot(filtered_df, aes(x = !!sym(input$x_channel), y = !!sym(input$y_channel))) +
          geom_point(color = "grey", alpha = 0.5) +
          geom_density_2d(color = "blue", size = 1) +
          annotate("text", x = Inf, y = Inf, 
                   label = paste("Total =", total_count),
                   hjust = 1.1, vjust = 1.1, size = 5, color = "black") +
          labs(x = paste("X Channel:", input$x_channel),
               y = paste("Y Channel:", input$y_channel),
               title = title_text) +
          theme(axis.title = element_text(size = 16),
                axis.text = element_text(size = 14),
                plot.title = element_text(size = 18, face = "bold"))
      } else if (input$plot2d_style == "Contour") {
        ggplot(filtered_df, aes(x = !!sym(input$x_channel), y = !!sym(input$y_channel))) +
          geom_point(color = "grey", alpha = 0.5) +
          stat_density_2d(color = "blue", size = 1) +
          annotate("text", x = Inf, y = Inf, 
                   label = paste("Total =", total_count),
                   hjust = 1.1, vjust = 1.1, size = 5, color = "black") +
          labs(x = paste("X Channel:", input$x_channel),
               y = paste("Y Channel:", input$y_channel),
               title = title_text) +
          theme(axis.title = element_text(size = 16),
                axis.text = element_text(size = 14),
                plot.title = element_text(size = 18, face = "bold"))
      }
    }
  })
  
  output$plot <- renderPlot({
    currentPlot()
  })
  
  output$downloadPlot <- downloadHandler(
    filename = function() {
      paste0(fileName(), "_", input$plot_type, ".pdf")
    },
    content = function(file) {
      ggsave(file, plot = currentPlot(), device = "pdf", width = 8, height = 6)
    }
  )
  
  overallSummary <- reactive({
    req(dataFCS())
    df <- dataFCS()
    t(sapply(df, function(x) {
      if (is.numeric(x)) {
        c(n = length(x),
          min = min(x, na.rm = TRUE),
          median = median(x, na.rm = TRUE),
          mean = mean(x, na.rm = TRUE),
          max = max(x, na.rm = TRUE))
      } else {
        rep(NA, 5)
      }
    }))
  })
  
  output$summary_table <- renderTable({
    overallSummary()
  }, rownames = TRUE)
  
  output$downloadSummary <- downloadHandler(
    filename = function() {
      paste0(fileName(), "_overall_summary.csv")
    },
    content = function(file) {
      write.csv(overallSummary(), file)
    }
  )
  
  gatedSummary <- reactive({
    req(input$plot_brush, dataFCS(), input$x_channel, input$y_channel)
    df <- dataFCS()
    filtered_df <- df[df[[input$x_channel]] >= input$x_range[1] &
                        df[[input$x_channel]] <= input$x_range[2] &
                        df[[input$y_channel]] >= input$y_range[1] &
                        df[[input$y_channel]] <= input$y_range[2], ]
    gated <- brushedPoints(filtered_df, input$plot_brush)
    if (nrow(gated) == 0) return(NULL)
    t(sapply(gated[, c(input$x_channel, input$y_channel)], function(x) {
      c(n = length(x),
        min = min(x, na.rm = TRUE),
        median = median(x, na.rm = TRUE),
        mean = mean(x, na.rm = TRUE),
        max = max(x, na.rm = TRUE))
    }))
  })
  
  output$gated_summary_table <- renderTable({
    gatedSummary()
  }, rownames = TRUE)
  
  output$downloadGated <- downloadHandler(
    filename = function() {
      paste0(fileName(), "_gated_summary.csv")
    },
    content = function(file) {
      gs <- gatedSummary()
      if (!is.null(gs)) write.csv(gs, file) else write.csv(data.frame(), file)
    }
  )
}

shinyApp(ui, server)