options(shiny.maxRequestSize = 100000 * 1024^2)

library(shiny)
library(flowCore)
library(ggplot2)
library(rlang)
library(tools)      # for file_path_sans_ext()
library(viridis)    # for color scales
library(MASS)       # for simulation (if used)
library(Rtsne)      # for t-SNE analysis
library(scales)

# Helper function: reads the FCS file and returns a data frame
# Improve the read_fcs_df function with better error handling
read_fcs_df <- function(fp) {
  message("Inside helper, using file path: ", fp)
  tryCatch({
    if (!file.exists(fp)) {
      stop("File does not exist: ", fp)
    }
    fcs_obj <- read.FCS(fp, transformation = FALSE)
    df <- as.data.frame(flowCore::exprs(fcs_obj))
    # Remove rows with zero values
    original_rows <- nrow(df)
    df <- df[!apply(df == 0, 1, any), ]
    removed_rows <- original_rows - nrow(df)
    if (removed_rows > 0) {
      message("Removed ", removed_rows, " rows containing zero values")
    }
    if (nrow(df) == 0) {
      warning("No data remaining after removing zeros")
      return(NULL)
    }
    return(df)
  }, error = function(e) {
    message("Error reading FCS file: ", e$message)
    return(NULL)
  })
}

ui <- fluidPage(
  titlePanel("Flow Cytometry Data Viewer"),
  # Add to UI sidebar
  conditionalPanel(
    condition = "input.simulate == false",
    fileInput("fcs_files", "Choose FCS Files (max 8)", multiple = TRUE, accept = ".fcs"),
    # Add this checkbox for zero removal option
    checkboxInput("remove_zeros", "Remove Zero Values", value = TRUE)
  ),
  sidebarLayout(
    sidebarPanel(
      # Data source controls
      checkboxInput("simulate", "Simulate Data", value = FALSE),
      conditionalPanel(
        condition = "input.simulate == true",
        sliderInput("numPopulations", "Number of Populations", min = 1, max = 8, value = 1, step = 1),
        actionButton("simulateButton", "Generate Simulated Data")
      ),
      conditionalPanel(
        condition = "input.simulate == false",
        fileInput("fcs_files", "Choose FCS Files (max 8)", multiple = TRUE, accept = ".fcs")
      ),
      # Plot type selection
      radioButtons("plot_type", "Plot Type",
                   choices = c("2D Plot", "Histogram"), selected = "2D Plot", inline = TRUE),
      # Reset button
      actionButton("reset_button", "Reset Plot Ranges", 
                   style = "color: #fff; background-color: #337ab7; border-color: #2e6da4"),
      # Histogram controls: channel and a slider for x-axis range
      conditionalPanel(
        condition = "input.plot_type == 'Histogram'",
        uiOutput("channel_select"),
        uiOutput("hist_slider"),
        checkboxInput("log_scale_x_hist", "Log Scale X-axis", value = FALSE)
      ),
      # Add transparency control to the UI in the 2D Plot conditionalPanel
      conditionalPanel(
        condition = "input.plot_type == '2D Plot'",
        uiOutput("x_channel_select"),
        uiOutput("y_channel_select"),
        uiOutput("x_slider"),
        uiOutput("y_slider"),
        radioButtons("plot2d_style", "2D Plot Style",
                    choices = c("Scatter", "Density", "Contour"), inline = TRUE),
        # Add transparency slider for scatter plots
        conditionalPanel(
          condition = "input.plot2d_style == 'Scatter'",
          sliderInput("scatter_alpha", "Point Transparency", 
                      min = 0.1, max = 1.0, value = 0.6, step = 0.1)
        ),
        checkboxInput("sync_brushing", "Synchronize Selection Across Plots", value = TRUE),
        checkboxInput("persistent_brush", "Keep Brush Selections Visible", value = TRUE),
        # Additive Gating option
        checkboxInput("additive_gating", "Progressive Filtering (Additive Gating)", value = FALSE),
        # NEW GATING CONTROLS GO HERE:
        radioButtons("gating_mode", "Gating Mode",
                    choices = c("Filter Data" = "filter", "Visual Gates Only" = "visual"), 
                    selected = "filter", inline = TRUE),
        checkboxInput("individual_gates", "Individual Plot Gates", value = FALSE),
        conditionalPanel(
          condition = "input.individual_gates == true",
          h5("Individual Gate Controls:"),
          fluidRow(
            column(4, actionButton("clear_gate_1", "Clear Gate 1", style = "font-size: 10px; padding: 2px 8px;")),
            column(4, actionButton("clear_gate_2", "Clear Gate 2", style = "font-size: 10px; padding: 2px 8px;")),
            column(4, actionButton("clear_gate_3", "Clear Gate 3", style = "font-size: 10px; padding: 2px 8px;"))
          ),
          fluidRow(
            column(4, actionButton("clear_gate_4", "Clear Gate 4", style = "font-size: 10px; padding: 2px 8px;")),
            column(4, actionButton("clear_gate_5", "Clear Gate 5", style = "font-size: 10px; padding: 2px 8px;")),
            column(4, actionButton("clear_gate_6", "Clear Gate 6", style = "font-size: 10px; padding: 2px 8px;"))
          ),
          fluidRow(
            column(4, actionButton("clear_gate_7", "Clear Gate 7", style = "font-size: 10px; padding: 2px 8px;")),
            column(4, actionButton("clear_gate_8", "Clear Gate 8", style = "font-size: 10px; padding: 2px 8px;")),
            column(4, actionButton("clear_all_gates", "Clear All", style = "font-size: 10px; padding: 2px 8px; background-color: #d9534f; color: white;"))
          )
        ),
        br(),
        # Log scale options for 2D plots
        h5("Log Scale Options:"),
        checkboxInput("log_scale_x", "Log Scale X-axis", value = FALSE),
        checkboxInput("log_scale_y", "Log Scale Y-axis", value = FALSE)
      ),
      br(),
      # t-SNE parameters
      selectInput("tsne_dims", "t-SNE Dimensions", 
                  choices = c("All" = "all", "PCA-reduced (faster)" = "pca"), 
                  selected = "pca"),
      sliderInput("tsne_perplexity", "t-SNE Perplexity", min = 5, max = 50, value = 30),
      sliderInput("tsne_iter", "t-SNE Iterations", min = 250, max = 1000, value = 500),
      br(),
      # Reset All button
      actionButton("reset_all_button", "Reset All", 
                   style = "color: #fff; background-color: #d9534f; border-color: #d43f3d"),
      br(), br(),
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
      tabsetPanel(
        tabPanel("Diagnostics", 
                 verbatimTextOutput("diagnosticInfo")),
        tabPanel("Plots",
                 fluidRow(
                   column(width = 6, 
                          plotOutput("plot1", brush = brushOpts(id = "plot1_brush"), height = "300px")
                   ),
                   column(width = 6, 
                          plotOutput("plot2", brush = brushOpts(id = "plot2_brush"), height = "300px")
                   )
                 ),
                 fluidRow(
                   column(width = 6, 
                          plotOutput("plot3", brush = brushOpts(id = "plot3_brush"), height = "300px")
                   ),
                   column(width = 6, 
                          plotOutput("plot4", brush = brushOpts(id = "plot4_brush"), height = "300px")
                   )
                 ),
                 fluidRow(
                   column(width = 6, 
                          plotOutput("plot5", brush = brushOpts(id = "plot5_brush"), height = "300px")
                   ),
                   column(width = 6, 
                          plotOutput("plot6", brush = brushOpts(id = "plot6_brush"), height = "300px")
                   )
                 ),
                 fluidRow(
                   column(width = 6, 
                          plotOutput("plot7", brush = brushOpts(id = "plot7_brush"), height = "300px")
                   ),
                   column(width = 6, 
                          plotOutput("plot8", brush = brushOpts(id = "plot8_brush"), height = "300px")
                   )
                 )
        ),
        tabPanel("t-SNE",
                 selectInput("tsne_file_select", "Select File for t-SNE", choices = c("No files loaded")),
                 plotOutput("tsne_plot", brush = brushOpts(id = "tsne_brush"))
        ),
        tabPanel("Summary",
                 selectInput("summary_file_select", "Select File for Summary", choices = c("No files loaded")),
                 h3("Overall Data Summary"),
                 tableOutput("summary_table"),
                 conditionalPanel(
                   condition = "input.plot_type == '2D Plot' && input.additive_gating",
                   h3("Gated Data Summary"),
                   tableOutput("gated_summary_table")
                 )
        )
      )
    )
  )
)

server <- function(input, output, session) {
  
  # Store all loaded FCS data
  loaded_data <- reactiveVal(list())
  file_names <- reactiveVal(character(0))
  # Store individual plot gates for visual gating
  individual_gates <- reactiveVal(list())
  
  # Store original data ranges
  original_ranges <- reactiveVal(NULL)
  
  # Reactive range values for all plots
  range_values <- reactiveValues(
    x_min = 0,
    x_max = 100,
    y_min = 0,
    y_max = 100
  )
  
  # Store brushed selections for persistence
  brushed_selections <- reactiveVal(list())
  
  # Store progressively filtered data for additive gating
  filtered_gates <- reactiveVal(list())
  
  # Diagnostic output
  output$diagnosticInfo <- renderPrint({
    data_list <- loaded_data()
    names_list <- file_names()
    
    cat("Number of loaded files:", length(data_list), "\n")
    cat("File names:", paste(names_list, collapse=", "), "\n\n")
    
    if (length(data_list) > 0) {
      cat("Row counts before filtering:\n")
      for (i in seq_along(data_list)) {
        cat("File", i, ":", names_list[i], "- Rows:", nrow(data_list[[i]]), "\n")
      }
      
      filtered <- filteredData()
      cat("\nRow counts after filtering:\n")
      for (i in seq_along(filtered)) {
        cat("File", i, ":", names_list[i], "- Rows:", nrow(filtered[[i]]), "\n")
      }
      
      if (input$additive_gating && input$plot_type == "2D Plot") {
        gates <- filtered_gates()
        if (length(gates) > 0) {
          cat("\nRow counts after additive gating:\n")
          for (i in seq_along(gates)) {
            if (!is.null(gates[[i]]) && nrow(gates[[i]]) > 0) {
              cat("File", i, ":", names_list[i], "- Rows:", nrow(gates[[i]]), "\n")
            }
          }
        }
      }
      
      bs <- brushed_selections()
      if (length(bs) > 0) {
        cat("\nPersistent Brush Selections:\n")
        for (i in seq_along(bs)) {
          if (!is.null(bs[[i]])) {
            cat("Plot", i, ": Active\n")
          }
        }
      }
    }
  })
  
  # Overall Data Summary
  output$summary_table <- renderTable({
    req(loaded_data(), input$summary_file_select)
    names_list <- file_names()
    idx <- match(input$summary_file_select, names_list)
    data_list <- loaded_data()
    
    if (!is.na(idx) && idx <= length(data_list)) {
      df <- data_list[[idx]]
      summary_stats <- data.frame(
        Channel = names(df),
        Mean = sapply(df, function(x) mean(x, na.rm = TRUE)),
        Median = sapply(df, function(x) median(x, na.rm = TRUE)),
        SD = sapply(df, function(x) sd(x, na.rm = TRUE)),
        Min = sapply(df, function(x) min(x, na.rm = TRUE)),
        Max = sapply(df, function(x) max(x, na.rm = TRUE)),
        N = sapply(df, function(x) sum(!is.na(x)))
      )
      summary_stats
    } else {
      NULL
    }
  })
  
  # Gated Data Summary (only for 2D plots with gating)
  output$gated_summary_table <- renderTable({
    req(input$plot_type == "2D Plot", input$summary_file_select, input$additive_gating)
    gates <- filtered_gates()
    names_list <- file_names()
    idx <- match(input$summary_file_select, names_list)
    data_list <- loaded_data()
    
    if (!is.na(idx) && idx <= length(gates) && idx <= length(data_list)) {
      gated_df <- gates[[idx]]
      full_df <- data_list[[idx]]
      
      if (!is.null(gated_df) && nrow(gated_df) > 0) {
        gated_stats <- data.frame(
          Channel = names(gated_df),
          Mean = sapply(gated_df, function(x) mean(x, na.rm = TRUE)),
          Median = sapply(gated_df, function(x) median(x, na.rm = TRUE)),
          SD = sapply(gated_df, function(x) sd(x, na.rm = TRUE)),
          Min = sapply(gated_df, function(x) min(x, na.rm = TRUE)),
          Max = sapply(gated_df, function(x) max(x, na.rm = TRUE)),
          N = sapply(gated_df, function(x) sum(!is.na(x))),
          PercentGated = round(100 * nrow(gated_df) / nrow(full_df), 2)
        )
        gated_stats
      } else {
        NULL
      }
    } else {
      NULL
    }
  })
  
  output$tsne_plot <- renderPlot({
    req(tsne_data(), input$tsne_dims)
    df <- tsne_data()
    
    numeric_df <- df[, sapply(df, is.numeric)]
    
    if (input$tsne_dims == "pca" && ncol(numeric_df) > 10) {
      pca_res <- prcomp(numeric_df, scale. = TRUE)
      tsne_input <- pca_res$x[, 1:min(10, ncol(pca_res$x))]
    } else {
      tsne_input <- numeric_df
    }
    
    set.seed(42)
    tsne_out <- Rtsne(tsne_input, perplexity = input$tsne_perplexity, max_iter = input$tsne_iter)
    
    tsne_df <- data.frame(
      TSNE1 = tsne_out$Y[, 1],
      TSNE2 = tsne_out$Y[, 2]
    )
    
    ggplot(tsne_df, aes(x = TSNE1, y = TSNE2)) +
      geom_point(alpha = 0.6, color = "darkblue") +
      labs(title = paste("t-SNE Plot:", input$tsne_file_select)) +
      theme_classic() +
      theme(axis.title = element_text(size = 14),
            axis.text = element_text(size = 12),
            plot.title = element_text(size = 16, face = "bold"))
  })
  
  observeEvent(input$reset_all_button, {
    data_ranges <- original_ranges()
    have_data <- length(loaded_data()) > 0
    
    if (have_data && !is.null(data_ranges)) {
      if (input$plot_type == "Histogram" && !is.null(input$channel)) {
        chan <- input$channel
        if (chan %in% names(data_ranges)) {
          x_min <- floor(data_ranges[[chan]][1])
          x_max <- ceiling(data_ranges[[chan]][2])
          updateSliderInput(session, "hist_xrange", min = x_min, max = x_max, value = c(x_min, x_max))
        }
      }
      if (!is.null(input$x_channel) && input$x_channel %in% names(data_ranges)) {
        x_min <- floor(data_ranges[[input$x_channel]][1])
        x_max <- ceiling(data_ranges[[input$x_channel]][2])
        range_values$x_min <- x_min
        range_values$x_max <- x_max
        updateSliderInput(session, "x_range", min = x_min, max = x_max, value = c(x_min, x_max))
        updateNumericInput(session, "x_min_input", value = x_min)
        updateNumericInput(session, "x_max_input", value = x_max)
      }
      if (!is.null(input$y_channel) && input$y_channel %in% names(data_ranges)) {
        y_min <- floor(data_ranges[[input$y_channel]][1])
        y_max <- ceiling(data_ranges[[input$y_channel]][2])
        range_values$y_min <- y_min
        range_values$y_max <- y_max
        updateSliderInput(session, "y_range", min = y_min, max = y_max, value = c(y_min, y_max))
      }
    } else {
      updateSliderInput(session, "hist_xrange", value = c(0, 100))
      updateSliderInput(session, "x_range", value = c(0, 100))
      updateSliderInput(session, "y_range", value = c(0, 100))
      updateNumericInput(session, "x_min_input", value = 0)
      updateNumericInput(session, "x_max_input", value = 100)
      range_values$x_min <- 0
      range_values$x_max <- 100
      range_values$y_min <- 0
      range_values$y_max <- 100
    }
    
    brushed_selections(list())
    filtered_gates(list())
    
    for (i in 1:8) {
      session$resetBrush(paste0("plot", i, "_brush"))
    }
    session$resetBrush("tsne_brush")
    updateCheckboxInput(session, "additive_gating", value = FALSE)
  })
  
  # Individual gate clearing observers
  lapply(1:8, function(i) {
    observeEvent(input[[paste0("clear_gate_", i)]], {
      gates <- individual_gates()
      if (length(gates) >= i) {
        gates[[i]] <- NULL
      }
      individual_gates(gates)
    })
  })

  observeEvent(input$clear_all_gates, {
    individual_gates(list())
    for (i in 1:8) {
      session$resetBrush(paste0("plot", i, "_brush"))
    }
  })


  observeEvent(input$simulateButton, {
    req(input$numPopulations)
    nPop <- input$numPopulations
    simData_list <- list()
    sim_names <- character(0)
    
    for (i in 1:min(8, nPop)) {
      dataList <- lapply(1:nPop, function(j) {
        mu <- c(runif(1, 0, 100), runif(1, 0, 100))
        sigma <- matrix(c(runif(1, 5, 20), 0, 0, runif(1, 5, 20)), nrow = 2)
        MASS::mvrnorm(100, mu = mu, Sigma = sigma)
      })
      simMatrix <- do.call(rbind, dataList)
      dfSim <- as.data.frame(simMatrix)
      names(dfSim) <- c("SimX", "SimY")
      
      simData_list[[i]] <- dfSim
      sim_names[i] <- paste0("SimulatedData_", i)
    }
    
    loaded_data(simData_list)
    file_names(sim_names)
    brushed_selections(list())
    filtered_gates(list())
    
    updateSelectInput(session, "tsne_file_select", choices = sim_names, selected = sim_names[1])
    updateSelectInput(session, "summary_file_select", choices = sim_names, selected = sim_names[1])
    
    calculate_original_ranges()
    
    for (i in 1:8) {
      session$resetBrush(paste0("plot", i, "_brush"))
    }
    session$resetBrush("tsne_brush")
  })
  
  observeEvent(input$fcs_files, {
    req(input$fcs_files)
    
    files_to_process <- head(input$fcs_files, 8)
    data_list <- list()
    names_list <- character(0)
    
    print("Files to process:")
    for (i in 1:nrow(files_to_process)) {
      print(paste("File", i, "Name:", files_to_process$name[i], "Path:", files_to_process$datapath[i]))
    }
    
    withProgress(message = "Loading FCS files...", value = 0, {
      successful_idx <- 1
      
      for (i in 1:nrow(files_to_process)) {
        fp <- files_to_process$datapath[i]
        name <- file_path_sans_ext(files_to_process$name[i])
        
        if (is.na(fp) || !file.exists(fp)) {
          message("Invalid file path for file: ", name, " (", fp, ")")
          next
        }
        
        incProgress(0.1 * i, detail = paste("Reading file", i, "of", nrow(files_to_process)))
        
        result <- tryCatch({
          df <- read_fcs_df(fp)
          if (!is.null(df)) {
            colnames(df) <- make.names(as.character(colnames(df)), unique = TRUE)
            list(success = TRUE, data = df)
          } else {
            list(success = FALSE, data = NULL)
          }
        }, error = function(e) {
          message("Error reading FCS file ", name, ": ", e$message)
          list(success = FALSE, data = NULL)
        })
        
        if (result$success) {
          data_list[[successful_idx]] <- result$data
          names_list[successful_idx] <- name
          successful_idx <- successful_idx + 1
        }
      }
    })
    
    if (length(data_list) > 0) {
      loaded_data(data_list)
      file_names(names_list)
      brushed_selections(list())
      filtered_gates(list())
      
      updateSelectInput(session, "tsne_file_select", choices = names_list, selected = names_list[1])
      updateSelectInput(session, "summary_file_select", choices = names_list, selected = names_list[1])
      
      calculate_original_ranges()
      
      for (i in 1:8) {
        session$resetBrush(paste0("plot", i, "_brush"))
      }
      session$resetBrush("tsne_brush")
    } else {
      showNotification("No valid FCS files could be loaded", type = "error")
    }
  })
  
  calculate_original_ranges <- function() {
    data_list <- loaded_data()
    if (length(data_list) == 0) return()
    
    channels <- names(data_list[[1]])
    ranges <- list()
    
    for (chan in channels) {
      chan_min <- Inf
      chan_max <- -Inf
      
      for (df in data_list) {
        if (chan %in% names(df)) {
          values <- df[[chan]]
          chan_min <- min(chan_min, min(values, na.rm = TRUE))
          chan_max <- max(chan_max, max(values, na.rm = TRUE))
        }
      }
      
      if (is.finite(chan_min) && is.finite(chan_max)) {
        ranges[[chan]] <- c(chan_min, chan_max)
      }
    }
    
    original_ranges(ranges)
  }
  
  observeEvent(input$reset_button, {
    req(loaded_data(), original_ranges())
    ranges <- original_ranges()
    
    if (input$plot_type == "Histogram" && !is.null(input$channel)) {
      chan <- input$channel
      if (chan %in% names(ranges)) {
        x_min <- floor(ranges[[chan]][1])
        x_max <- ceiling(ranges[[chan]][2])
        updateSliderInput(session, "hist_xrange", min = x_min, max = x_max, value = c(x_min, x_max))
      }
    } else if (input$plot_type == "2D Plot" && !is.null(input$x_channel) && !is.null(input$y_channel)) {
      x_chan <- input$x_channel
      y_chan <- input$y_channel
      
      if (x_chan %in% names(ranges)) {
        x_min <- floor(ranges[[x_chan]][1])
        x_max <- ceiling(ranges[[x_chan]][2])
        range_values$x_min <- x_min
        range_values$x_max <- x_max
        updateSliderInput(session, "x_range", min = x_min, max = x_max, value = c(x_min, x_max))
        updateNumericInput(session, "x_min_input", value = x_min)
        updateNumericInput(session, "x_max_input", value = x_max)
      }
      
      if (y_chan %in% names(ranges)) {
        y_min <- floor(ranges[[y_chan]][1])
        y_max <- ceiling(ranges[[y_chan]][2])
        range_values$y_min <- y_min
        range_values$y_max <- y_max
        updateSliderInput(session, "y_range", min = y_min, max = y_max, value = c(y_min, y_max))
      }
    }
    
    brushed_selections(list())
    filtered_gates(list())
    
    for (i in 1:8) {
      session$resetBrush(paste0("plot", i, "_brush"))
    }
  })
  
  output$channel_select <- renderUI({
    req(loaded_data())
    data_list <- loaded_data()
    if (length(data_list) == 0) return(NULL)
    selectInput("channel", "Select Channel", choices = names(data_list[[1]]))
  })
  
  output$x_channel_select <- renderUI({
    req(loaded_data())
    data_list <- loaded_data()
    if (length(data_list) == 0) return(NULL)
    selectInput("x_channel", "Select X Channel", choices = names(data_list[[1]]))
  })
  
  output$y_channel_select <- renderUI({
    req(loaded_data())
    data_list <- loaded_data()
    if (length(data_list) == 0) return(NULL)
    channels <- names(data_list[[1]])
    default_y <- ifelse("SSC.A" %in% channels, "SSC.A", channels[min(2, length(channels))])
    selectInput("y_channel", "Select Y Channel", choices = channels, selected = default_y)
  })
  
  observeEvent(input$channel, {
    req(loaded_data(), input$channel)
    data_list <- loaded_data()
    if (length(data_list) == 0) return()
    
    ranges <- original_ranges()
    if (!is.null(ranges) && input$channel %in% names(ranges)) {
      x_min <- floor(ranges[[input$channel]][1])
      x_max <- ceiling(ranges[[input$channel]][2])
      updateSliderInput(session, "hist_xrange", min = x_min, max = x_max, value = c(x_min, x_max))
    } else {
      x_min <- Inf; x_max <- -Inf
      for (df in data_list) {
        if (input$channel %in% names(df)) {
          x_values <- df[[input$channel]]
          x_min <- min(x_min, min(x_values, na.rm = TRUE))
          x_max <- max(x_max, max(x_values, na.rm = TRUE))
        }
      }
      if (is.finite(x_min) && is.finite(x_max)) {
        x_min <- floor(x_min)
        x_max <- ceiling(x_max)
        updateSliderInput(session, "hist_xrange", min = x_min, max = x_max, value = c(x_min, x_max))
      }
    }
  })
  
  observeEvent(list(input$x_channel, input$y_channel), {
    req(loaded_data(), input$x_channel, input$y_channel)
    data_list <- loaded_data()
    if (length(data_list) == 0) return()
    
    filtered_gates(list())
    
    ranges <- original_ranges()
    if (!is.null(ranges) && input$x_channel %in% names(ranges) && input$y_channel %in% names(ranges)) {
      x_min <- floor(ranges[[input$x_channel]][1])
      x_max <- ceiling(ranges[[input$x_channel]][2])
      y_min <- floor(ranges[[input$y_channel]][1])
      y_max <- ceiling(ranges[[input$y_channel]][2])
      
      range_values$x_min <- x_min
      range_values$x_max <- x_max
      range_values$y_min <- y_min
      range_values$y_max <- y_max
      
      updateSliderInput(session, "x_range", min = x_min, max = x_max, value = c(x_min, x_max))
      updateSliderInput(session, "y_range", min = y_min, max = y_max, value = c(y_min, y_max))
      updateNumericInput(session, "x_min_input", value = x_min, min = x_min, max = x_max)
      updateNumericInput(session, "x_max_input", value = x_max, min = x_min, max = x_max)
    } else {
      x_min <- Inf; x_max <- -Inf; y_min <- Inf; y_max <- -Inf
      for (df in data_list) {
        if (input$x_channel %in% names(df) && input$y_channel %in% names(df)) {
          x_values <- df[[input$x_channel]]
          y_values <- df[[input$y_channel]]
          x_min <- min(x_min, min(x_values, na.rm = TRUE))
          x_max <- max(x_max, max(x_values, na.rm = TRUE))
          y_min <- min(y_min, min(y_values, na.rm = TRUE))
          y_max <- max(y_max, max(y_values, na.rm = TRUE))
        }
      }
      if (is.finite(x_min) && is.finite(x_max) && is.finite(y_min) && is.finite(y_max)) {
        x_min <- floor(x_min)
        x_max <- ceiling(x_max)
        y_min <- floor(y_min)
        y_max <- ceiling(y_max)
        range_values$x_min <- x_min
        range_values$x_max <- x_max
        range_values$y_min <- y_min
        range_values$y_max <- y_max
        updateSliderInput(session, "x_range", min = x_min, max = x_max, value = c(x_min, x_max))
        updateSliderInput(session, "y_range", min = y_min, max = y_max, value = c(y_min, y_max))
        updateNumericInput(session, "x_min_input", value = x_min, min = x_min, max = x_max)
        updateNumericInput(session, "x_max_input", value = x_max, min = x_min, max = x_max)
      }
    }
    
    brushed_selections(list())
    for (i in 1:8) {
      session$resetBrush(paste0("plot", i, "_brush"))
    }
  })
  
  output$hist_slider <- renderUI({
    req(loaded_data(), input$channel)
    data_list <- loaded_data()
    if (length(data_list) == 0) return(NULL)
    
    ranges <- original_ranges()
    if (!is.null(ranges) && input$channel %in% names(ranges)) {
      x_min <- floor(ranges[[input$channel]][1])
      x_max <- ceiling(ranges[[input$channel]][2])
    } else {
      x_min <- 0; x_max <- 100
      for (df in data_list) {
        if (input$channel %in% names(df)) {
          x_values <- df[[input$channel]]
          x_min <- min(x_min, floor(min(x_values, na.rm = TRUE)))
          x_max <- max(x_max, ceiling(max(x_values, na.rm = TRUE)))
        }
      }
    }
    
    range_diff <- x_max - x_min
    step_size <- max(1, round(range_diff / 100))
    
    sliderInput("hist_xrange", "Histogram X Range", min = x_min, max = x_max, value = c(x_min, x_max), step = step_size, ticks = FALSE)
  })
  
  # Modify the x_slider output to include numeric inputs
  output$x_slider <- renderUI({
    req(range_values$x_min, range_values$x_max)
    range_diff <- range_values$x_max - range_values$x_min
    step_size <- max(1, round(range_diff / 100))
    x_min <- floor(range_values$x_min)
    x_max <- ceiling(range_values$x_max)
    
    tagList(
      sliderInput("x_range", "X Range", min = x_min, max = x_max, 
                  value = c(x_min, x_max), step = step_size, ticks = FALSE),
      fluidRow(
        column(6, numericInput("x_min_input", "X Min:", value = x_min, 
                               min = x_min, max = x_max, step = step_size)),
        column(6, numericInput("x_max_input", "X Max:", value = x_max, 
                               min = x_min, max = x_max, step = step_size))
      )
    )
  })
  
  output$y_slider <- renderUI({
    req(range_values$y_min, range_values$y_max)
    range_diff <- range_values$y_max - range_values$y_min
    step_size <- max(1, round(range_diff / 100))
    y_min <- floor(range_values$y_min)
    y_max <- ceiling(range_values$y_max)
    sliderInput("y_range", "Y Range", min = y_min, max = y_max, value = c(y_min, y_max), step = step_size, ticks = FALSE)
  })
  
  # Add observers for numeric inputs to sync with sliders
  observeEvent(input$x_min_input, {
    req(input$x_min_input, input$x_range)
    if (!is.na(input$x_min_input)) {
      new_range <- c(input$x_min_input, input$x_range[2])
      updateSliderInput(session, "x_range", value = new_range)
    }
  })
  
  observeEvent(input$x_max_input, {
    req(input$x_max_input, input$x_range)
    if (!is.na(input$x_max_input)) {
      new_range <- c(input$x_range[1], input$x_max_input)
      updateSliderInput(session, "x_range", value = new_range)
    }
  })
  
  observe({
    if (input$plot_type != "2D Plot") return()
    
    if (input$individual_gates) {
      # Individual plot gating
      gates <- individual_gates()
      updated <- FALSE
      
      for (i in 1:8) {
        brush_name <- paste0("plot", i, "_brush")
        current_brush <- input[[brush_name]]
        
        if (!is.null(current_brush)) {
          if (length(gates) < i) {
            gates <- c(gates, vector("list", i - length(gates)))
          }
          gates[[i]] <- current_brush
          updated <- TRUE
        }
      }
      
      if (updated) {
        individual_gates(gates)
      }
    } else {
      # Existing synchronized brushing logic
      if (input$persistent_brush) {
        bs <- list()
        for (i in 1:8) {
          brush_name <- paste0("plot", i, "_brush")
          bs[[i]] <- input[[brush_name]]
        }
        brushed_selections(bs)
      } else {
        active_brush <- NULL
        active_idx <- NULL
        for (i in 1:8) {
          brush_name <- paste0("plot", i, "_brush")
          if (!is.null(input[[brush_name]])) {
            active_brush <- input[[brush_name]]
            active_idx <- i
            break
          }
        }
        if (!is.null(active_idx)) {
          bs <- vector("list", 8)
          for (i in 1:8) {
            bs[[i]] <- if (i == active_idx) active_brush else NULL
          }
          brushed_selections(bs)
        } else {
          brushed_selections(list())
        }
      }
    }
  })
  
  observeEvent(input$additive_gating, {
    if (!input$additive_gating) {
      filtered_gates(list())
      for (i in 1:8) {
        session$resetBrush(paste0("plot", i, "_brush"))
      }
    }
  }, ignoreInit = TRUE)
  
  observe({
    req(input$additive_gating, input$plot_type == "2D Plot")
    
    active_brush <- NULL
    active_idx <- NULL
    for (i in 1:8) {
      brush_name <- paste0("plot", i, "_brush")
      if (!is.null(input[[brush_name]])) {
        active_brush <- input[[brush_name]]
        active_idx <- i
        break
      }
    }
    if (is.null(active_brush) || is.null(active_idx)) return()
    
    req(input$x_channel, input$y_channel)
    current_filtered <- filteredData()
    if (active_idx > length(current_filtered)) return()
    df <- current_filtered[[active_idx]]
    
    gates <- filtered_gates()
    if (length(gates) > 0 && active_idx <= length(gates) &&
        !is.null(gates[[active_idx]]) && nrow(gates[[active_idx]]) > 0) {
      df <- gates[[active_idx]]
    }
    if (!all(c(input$x_channel, input$y_channel) %in% names(df))) return()
    
    gated_df <- tryCatch({
      result <- brushedPoints(df, active_brush)
      if (is.data.frame(result) && nrow(result) > 0) result else NULL
    }, error = function(e) {
      message("Error in brushing: ", e$message)
      NULL
    })
    
    if (!is.null(gated_df) && nrow(gated_df) > 0) {
      new_gates <- gates
      if (length(new_gates) == 0) {
        new_gates <- vector("list", length(current_filtered))
      }
      new_gates[[active_idx]] <- gated_df
      filtered_gates(new_gates)
    }
  })
  
  observe({
    req(input$sync_brushing, input$plot_type == "2D Plot")
    active_brush <- NULL
    active_idx <- NULL
    for (i in 1:8) {
      brush_name <- paste0("plot", i, "_brush")
      if (!is.null(input[[brush_name]])) {
        active_brush <- input[[brush_name]]
        active_idx <- i
        break
      }
    }
    if (!is.null(active_brush) && !is.null(active_idx)) {
      updateSliderInput(session, "x_range", value = c(active_brush$xmin, active_brush$xmax))
      updateSliderInput(session, "y_range", value = c(active_brush$ymin, active_brush$ymax))
      updateNumericInput(session, "x_min_input", value = active_brush$xmin)
      updateNumericInput(session, "x_max_input", value = active_brush$xmax)
      
      range_values$x_min <- active_brush$xmin
      range_values$x_max <- active_brush$xmax
      range_values$y_min <- active_brush$ymin
      range_values$y_max <- active_brush$ymax
      
      if (!input$persistent_brush) {
        for (i in 1:8) {
          if (i != active_idx) {
            session$resetBrush(paste0("plot", i, "_brush"))
          }
        }
      }
    }
  })
  
  filteredData <- reactive({
    req(loaded_data())
    data_list <- loaded_data()
    if (length(data_list) == 0) return(list())
    
    if (input$plot_type == "Histogram") {
      req(input$hist_xrange, input$channel)
      lapply(data_list, function(df) {
        if (input$channel %in% names(df)) {
          df[df[[input$channel]] >= input$hist_xrange[1] & df[[input$channel]] <= input$hist_xrange[2], ]
        } else {
          df
        }
      })
    } else {
      req(input$x_range, input$y_range, input$x_channel, input$y_channel)
      lapply(data_list, function(df) {
        if (input$x_channel %in% names(df) && input$y_channel %in% names(df)) {
          df[df[[input$x_channel]] >= input$x_range[1] & df[[input$x_channel]] <= input$x_range[2] &
             df[[input$y_channel]] >= input$y_range[1] & df[[input$y_channel]] <= input$y_range[2], ]
        } else {
          df
        }
      })
    }
  })
  
  plotData <- reactive({
    req(filteredData())
    filtered <- filteredData()
    if (input$additive_gating && input$plot_type == "2D Plot") {
      gates <- filtered_gates()
      if (length(gates) > 0) {
        result <- list()
        for (i in seq_along(filtered)) {
          if (i <= length(gates) && !is.null(gates[[i]]) && nrow(gates[[i]]) > 0) {
            result[[i]] <- gates[[i]]
          } else {
            result[[i]] <- filtered[[i]]
          }
        }
        return(result)
      }
    }
    filtered
  })
  
  generate_plot <- function(df, file_idx, file_name, brush) {
    if (is.null(df) || nrow(df) == 0) {
      return(ggplot() + 
               annotate("text", x = 0.5, y = 0.5, label = "No data available") +
               theme_void())
    }
    
    # Handle zero/negative values for log scales by converting them to 1
    if (input$plot_type == "Histogram" && input$log_scale_x_hist) {
      req(input$channel)
      if (input$channel %in% names(df)) {
        df[[input$channel]][df[[input$channel]] <= 0 | is.na(df[[input$channel]])] <- 1
      }
    } else if (input$plot_type == "2D Plot") {
      if (input$log_scale_x) {
        req(input$x_channel)
        if (input$x_channel %in% names(df)) {
          df[[input$x_channel]][df[[input$x_channel]] <= 0 | is.na(df[[input$x_channel]])] <- 1
        }
      }
      if (input$log_scale_y) {
        req(input$y_channel)
        if (input$y_channel %in% names(df)) {
          df[[input$y_channel]][df[[input$y_channel]] <= 0 | is.na(df[[input$y_channel]])] <- 1
        }
      }
    }
    
    if (input$plot_type == "Histogram") {
      req(input$channel)
      if (!(input$channel %in% names(df))) {
        return(ggplot() + 
                 annotate("text", x = 0.5, y = 0.5, label = paste("Channel", input$channel, "not found in file", file_idx)) +
                 theme_void())
      }
      x_data <- df[[input$channel]]
      m <- mean(x_data, na.rm = TRUE)
      med <- median(x_data, na.rm = TRUE)
      total_count <- nrow(df)
      
      p <- ggplot(df, aes(x = !!sym(input$channel))) +
        geom_histogram(bins = 30, fill = "blue", color = "black", alpha = 0.7) +
        geom_vline(xintercept = m, color = "red", linetype = "dashed") +
        geom_vline(xintercept = med, color = "green", linetype = "dotted") +
        annotate("text", x = Inf, y = Inf, label = paste("Total =", total_count),
                 hjust = 1.1, vjust = 1.1, size = 3, color = "black") +
        labs(x = paste("Channel:", input$channel),
             y = "Count",
             title = paste("File", file_idx, "-", file_name)) +
        theme_classic() +
        theme(axis.title = element_text(size = 14),
              axis.text = element_text(size = 12),
              plot.title = element_text(size = 16, face = "bold"))
      
      # Add xlim only if not using log scale (log scale handles its own range)
      if (!input$log_scale_x_hist) {
        p <- p + xlim(input$hist_xrange[1], input$hist_xrange[2])
      }
      
      # Add log scale for x-axis if selected
      if (input$log_scale_x_hist) {
        p <- p + scale_x_log10(labels = scales::label_number(accuracy = 1))
      } else {
        p <- p + scale_x_continuous(labels = scales::label_number(accuracy = 1))
      }
      
      return(p)
      
    } else {
      req(input$x_channel, input$y_channel)
      if (!(input$x_channel %in% names(df)) || !(input$y_channel %in% names(df))) {
        return(ggplot() + 
                 annotate("text", x = 0.5, y = 0.5, label = paste("Channels not found in file", file_idx)) +
                 theme_void())
      }
      total_count <- nrow(df)
      is_gated <- FALSE
      if (input$additive_gating) {
        gates <- filtered_gates()
        if (length(gates) >= file_idx && !is.null(gates[[file_idx]]) && nrow(gates[[file_idx]]) > 0) {
          is_gated <- TRUE
        }
      }
      if (input$persistent_brush && is.null(brush)) {
        saved_brushes <- brushed_selections()
        if (length(saved_brushes) >= file_idx) {
          brush <- saved_brushes[[file_idx]]
        }
      }
      # Replace the scatter plot section in generate_plot function
# Modify the generate_plot function - replace the scatter plot section (around line 808)

if (input$plot2d_style == "Scatter") {
  corr <- if (nrow(df) > 1) round(cor(df[[input$x_channel]], df[[input$y_channel]], use = "complete.obs"), 2) else NA
  brushed <- NULL
  current_brush <- NULL
  
  # Determine which brush to use
  if (input$individual_gates) {
    gates <- individual_gates()
    if (length(gates) >= file_idx && !is.null(gates[[file_idx]])) {
      current_brush <- gates[[file_idx]]
    }
  } else {
    current_brush <- brush
    if (input$persistent_brush && is.null(current_brush)) {
      saved_brushes <- brushed_selections()
      if (length(saved_brushes) >= file_idx) {
        current_brush <- saved_brushes[[file_idx]]
      }
    }
  }
  
  # Get brushed points
  if (!is.null(current_brush)) {
    brushed <- tryCatch({
      result <- brushedPoints(df, current_brush)
      if (is.data.frame(result) && nrow(result) > 0) result else NULL
    }, error = function(e) {
      NULL
    })
  }
  
  alpha_value <- if(is.null(input$scatter_alpha)) 0.6 else input$scatter_alpha
  
  # Create base plot
  p <- ggplot(df, aes(x = !!sym(input$x_channel), y = !!sym(input$y_channel))) +
    geom_point(color = ifelse(is_gated, "darkgreen", "black"), 
              size = 0.5, 
              alpha = alpha_value)
  
  # Add gated points overlay and statistics
  gate_percent <- 0
  if (!is.null(brushed) && nrow(brushed) > 0) {
    gate_percent <- round(100 * nrow(brushed) / nrow(df), 1)
    
    if (input$gating_mode == "visual" || input$individual_gates) {
      # Visual gating - highlight points without filtering
      p <- p + geom_point(data = brushed, 
                        color = "red", 
                        size = 1, 
                        alpha = min(1.0, alpha_value + 0.3)) +
        annotate("text", x = Inf, y = Inf,
                label = paste("Total =", total_count, 
                            "\nGated =", nrow(brushed), 
                            paste0("(", gate_percent, "%)")),
                hjust = 1.1, vjust = 1.1, size = 3, color = "black")
    } else {
      # Filter mode - existing behavior
      p <- p + geom_point(data = brushed, color = "red", size = 1, alpha = alpha_value) +
        annotate("text", x = Inf, y = Inf,
                label = paste("Total =", total_count, "\nGated =", nrow(brushed)),
                hjust = 1.1, vjust = 1.1, size = 3, color = "black")
    }
  } else {
    p <- p + annotate("text", x = Inf, y = Inf,
                    label = paste("Total =", total_count),
                    hjust = 1.1, vjust = 1.1, size = 3, color = "black")
  }
  
  # Add gate rectangle for ANY valid brush (not just individual gates)
  if (!is.null(current_brush)) {
    p <- p + 
      # Gate rectangle
      annotate("rect", 
              xmin = current_brush$xmin, xmax = current_brush$xmax,
              ymin = current_brush$ymin, ymax = current_brush$ymax,
              fill = NA, color = "red", size = 1, linetype = "dashed", alpha = 0.7) +
      # Gate percentage annotation
      annotate("text", 
              x = current_brush$xmax, 
              y = current_brush$ymax + (current_brush$ymax - current_brush$ymin) * 0.1,
              label = paste("Gate =", gate_percent, "%"), 
              hjust = 1, vjust = 0, 
              size = 4, color = "red", fontface = "bold",
              bg.colour = "white", bg.r = 0.1)
  }
  
  title_text <- paste("File", file_idx, "-", file_name, "\nCorr =", corr)
  if (is_gated) title_text <- paste0(title_text, " (Gated)")
  if (gate_percent > 0) title_text <- paste0(title_text, " | ", gate_percent, "% selected")
  
  # Create the result plot with existing coordinate and scale logic
  result_plot <- p + labs(x = paste("X:", input$x_channel),
                        y = paste("Y:", input$y_channel),
                        title = title_text) +
                  theme_classic() +
                  theme(axis.title = element_text(size = 14),
                        axis.text = element_text(size = 12),
                        plot.title = element_text(size = 16, face = "bold"))
                  
  # Apply coordinate limits only if not using log scales (log scales manage their own ranges)
  if (!input$log_scale_x && !input$log_scale_y) {
    result_plot <- result_plot + coord_cartesian(xlim = c(input$x_range[1], input$x_range[2]),
                                                ylim = c(input$y_range[1], input$y_range[2]))
  } else if (!input$log_scale_x && input$log_scale_y) {
    result_plot <- result_plot + coord_cartesian(xlim = c(input$x_range[1], input$x_range[2]))
  } else if (input$log_scale_x && !input$log_scale_y) {
    result_plot <- result_plot + coord_cartesian(ylim = c(input$y_range[1], input$y_range[2]))
  }
  
  # Apply scale transformations
  result_plot + 
    {if (input$log_scale_x) scale_x_log10(labels = scales::label_number(accuracy = 1)) else scale_x_continuous(labels = scales::label_number(accuracy = 1))} +
    {if (input$log_scale_y) scale_y_log10(labels = scales::label_number(accuracy = 1)) else scale_y_continuous(labels = scales::label_number(accuracy = 1))}
} else if (input$plot2d_style == "Density") {
  # Get current brush for density plots too
  current_brush <- NULL
  if (input$individual_gates) {
    gates <- individual_gates()
    if (length(gates) >= file_idx && !is.null(gates[[file_idx]])) {
      current_brush <- gates[[file_idx]]
    }
  } else {
    current_brush <- brush
    if (input$persistent_brush && is.null(current_brush)) {
      saved_brushes <- brushed_selections()
      if (length(saved_brushes) >= file_idx) {
        current_brush <- saved_brushes[[file_idx]]
      }
    }
  }
  
  p <- ggplot(df, aes(x = !!sym(input$x_channel), y = !!sym(input$y_channel))) +
    geom_point(color = ifelse(is_gated, "darkgreen", "black"), alpha = 0.5, size = 0.5) +
    stat_density_2d_filled(alpha = 0.5, contour_var = "density") +
    geom_density_2d(color = "blue", size = 0.5) +
    annotate("text", x = Inf, y = Inf,
             label = paste("Total =", total_count),
             hjust = 1.1, vjust = 1.1, size = 3, color = "black")
  
  # Add gate rectangle for density plots
  if (!is.null(current_brush)) {
    # Calculate gate percentage
    brushed <- tryCatch({
      result <- brushedPoints(df, current_brush)
      if (is.data.frame(result) && nrow(result) > 0) result else NULL
    }, error = function(e) {
      NULL
    })
    
    gate_percent <- if (!is.null(brushed) && nrow(brushed) > 0) {
      round(100 * nrow(brushed) / nrow(df), 1)
    } else {
      0
    }
    
    p <- p + 
      # Gate rectangle
      annotate("rect", 
              xmin = current_brush$xmin, xmax = current_brush$xmax,
              ymin = current_brush$ymin, ymax = current_brush$ymax,
              fill = NA, color = "red", size = 1, linetype = "dashed", alpha = 0.7) +
      # Gate percentage annotation
      annotate("text", 
              x = current_brush$xmax, 
              y = current_brush$ymax + (current_brush$ymax - current_brush$ymin) * 0.1,
              label = paste("Gate =", gate_percent, "%"), 
              hjust = 1, vjust = 0, 
              size = 4, color = "red", fontface = "bold",
              bg.colour = "white", bg.r = 0.1)
  }
  
  title_text <- paste("File", file_idx, "-", file_name)
  if (is_gated) title_text <- paste0(title_text, " (Gated)")
  
  # Create the base plot
  result_plot <- p + labs(x = paste("X:", input$x_channel),
                         y = paste("Y:", input$y_channel),
                         title = title_text) +
                theme_classic() +
                theme(axis.title = element_text(size = 14),
                      axis.text = element_text(size = 12),
                      plot.title = element_text(size = 16, face = "bold"),
                      legend.position = "none")
                        
  # Apply coordinate limits only if not using log scales
  if (!input$log_scale_x && !input$log_scale_y) {
    result_plot <- result_plot + coord_cartesian(xlim = c(input$x_range[1], input$x_range[2]),
                                                ylim = c(input$y_range[1], input$y_range[2]))
  } else if (!input$log_scale_x && input$log_scale_y) {
    result_plot <- result_plot + coord_cartesian(xlim = c(input$x_range[1], input$x_range[2]))
  } else if (input$log_scale_x && !input$log_scale_y) {
    result_plot <- result_plot + coord_cartesian(ylim = c(input$y_range[1], input$y_range[2]))
  }
  
  # Apply scale transformations
  result_plot + 
    {if (input$log_scale_x) scale_x_log10(labels = scales::label_number(accuracy = 1)) else scale_x_continuous(labels = scales::label_number(accuracy = 1))} +
    {if (input$log_scale_y) scale_y_log10(labels = scales::label_number(accuracy = 1)) else scale_y_continuous(labels = scales::label_number(accuracy = 1))}
  
} else if (input$plot2d_style == "Contour") {
  # Get current brush for contour plots too
  current_brush <- NULL
  if (input$individual_gates) {
    gates <- individual_gates()
    if (length(gates) >= file_idx && !is.null(gates[[file_idx]])) {
      current_brush <- gates[[file_idx]]
    }
  } else {
    current_brush <- brush
    if (input$persistent_brush && is.null(current_brush)) {
      saved_brushes <- brushed_selections()
      if (length(saved_brushes) >= file_idx) {
        current_brush <- saved_brushes[[file_idx]]
      }
    }
  }
  
  p <- ggplot(df, aes(x = !!sym(input$x_channel), y = !!sym(input$y_channel))) +
    geom_point(color = ifelse(is_gated, "darkgreen", "black"), alpha = 0.5, size = 0.5) +
    stat_density_2d(aes(color = ..level..), size = 0.7, bins = 10) +
    scale_color_viridis_c() +
    annotate("text", x = Inf, y = Inf,
             label = paste("Total =", total_count),
             hjust = 1.1, vjust = 1.1, size = 3, color = "black")
  
  # Add gate rectangle for contour plots
  if (!is.null(current_brush)) {
    # Calculate gate percentage
    brushed <- tryCatch({
      result <- brushedPoints(df, current_brush)
      if (is.data.frame(result) && nrow(result) > 0) result else NULL
    }, error = function(e) {
      NULL
    })
    
    gate_percent <- if (!is.null(brushed) && nrow(brushed) > 0) {
      round(100 * nrow(brushed) / nrow(df), 1)
    } else {
      0
    }
    
    p <- p + 
      # Gate rectangle
      annotate("rect", 
              xmin = current_brush$xmin, xmax = current_brush$xmax,
              ymin = current_brush$ymin, ymax = current_brush$ymax,
              fill = NA, color = "red", size = 1, linetype = "dashed", alpha = 0.7) +
      # Gate percentage annotation
      annotate("text", 
              x = current_brush$xmax, 
              y = current_brush$ymax + (current_brush$ymax - current_brush$ymin) * 0.1,
              label = paste("Gate =", gate_percent, "%"), 
              hjust = 1, vjust = 0, 
              size = 4, color = "red", fontface = "bold",
              bg.colour = "white", bg.r = 0.1)
  }
  
  title_text <- paste("File", file_idx, "-", file_name)
  if (is_gated) title_text <- paste0(title_text, " (Gated)")
  
  # Create the base plot
  result_plot <- p + labs(x = paste("X:", input$x_channel),
                         y = paste("Y:", input$y_channel),
                         title = title_text) +
                theme_classic() +
                theme(axis.title = element_text(size = 14),
                      axis.text = element_text(size = 12),
                      plot.title = element_text(size = 16, face = "bold"),
                      legend.position = "none")
  
  # Apply coordinate limits only if not using log scales
  if (!input$log_scale_x && !input$log_scale_y) {
    result_plot <- result_plot + coord_cartesian(xlim = c(input$x_range[1], input$x_range[2]),
                                                ylim = c(input$y_range[1], input$y_range[2]))
  } else if (!input$log_scale_x && input$log_scale_y) {
    result_plot <- result_plot + coord_cartesian(xlim = c(input$x_range[1], input$x_range[2]))
  } else if (input$log_scale_x && !input$log_scale_y) {
    result_plot <- result_plot + coord_cartesian(ylim = c(input$y_range[1], input$y_range[2]))
  }
  
  # Apply scale transformations
  result_plot + 
    {if (input$log_scale_x) scale_x_log10(labels = scales::label_number(accuracy = 1)) else scale_x_continuous(labels = scales::label_number(accuracy = 1))} +
    {if (input$log_scale_y) scale_y_log10(labels = scales::label_number(accuracy = 1)) else scale_y_continuous(labels = scales::label_number(accuracy = 1))}
}
    }
  }
  
  output$plot1 <- renderPlot({
    req(loaded_data(), file_names())
    data_list <- plotData()
    names_list <- file_names()
    if (length(data_list) >= 1) {
      tryCatch({
        generate_plot(data_list[[1]], 1, names_list[1], input$plot1_brush)
      }, error = function(e) {
        message("Error in plot 1: ", e$message)
        ggplot() + annotate("text", x = 0.5, y = 0.5, label = paste("Error:", e$message)) + theme_void()
      })
    } else {
      ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No data loaded for plot 1") + theme_void()
    }
  })
  
  output$plot2 <- renderPlot({
    req(loaded_data(), file_names())
    data_list <- plotData()
    names_list <- file_names()
    if (length(data_list) >= 2) {
      tryCatch({
        generate_plot(data_list[[2]], 2, names_list[2], input$plot2_brush)
      }, error = function(e) {
        message("Error in plot 2: ", e$message)
        ggplot() + annotate("text", x = 0.5, y = 0.5, label = paste("Error:", e$message)) + theme_void()
      })
    } else {
      ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No data loaded for plot 2") + theme_void()
    }
  })
  
  output$plot3 <- renderPlot({
    req(loaded_data(), file_names())
    data_list <- plotData()
    names_list <- file_names()
    if (length(data_list) >= 3) {
      tryCatch({
        generate_plot(data_list[[3]], 3, names_list[3], input$plot3_brush)
      }, error = function(e) {
        message("Error in plot 3: ", e$message)
        ggplot() + annotate("text", x = 0.5, y = 0.5, label = paste("Error:", e$message)) + theme_void()
      })
    } else {
      ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No data loaded for plot 3") + theme_void()
    }
  })
  
  output$plot4 <- renderPlot({
    req(loaded_data(), file_names())
    data_list <- plotData()
    names_list <- file_names()
    if (length(data_list) >= 4) {
      tryCatch({
        generate_plot(data_list[[4]], 4, names_list[4], input$plot4_brush)
      }, error = function(e) {
        message("Error in plot 4: ", e$message)
        ggplot() + annotate("text", x = 0.5, y = 0.5, label = paste("Error:", e$message)) + theme_void()
      })
    } else {
      ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No data loaded for plot 4") + theme_void()
    }
  })
  
  output$plot5 <- renderPlot({
    req(loaded_data(), file_names())
    data_list <- plotData()
    names_list <- file_names()
    if (length(data_list) >= 5) {
      tryCatch({
        generate_plot(data_list[[5]], 5, names_list[5], input$plot5_brush)
      }, error = function(e) {
        message("Error in plot 5: ", e$message)
        ggplot() + annotate("text", x = 0.5, y = 0.5, label = paste("Error:", e$message)) + theme_void()
      })
    } else {
      ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No data loaded for plot 5") + theme_void()
    }
  })
  
  output$plot6 <- renderPlot({
    req(loaded_data(), file_names())
    data_list <- plotData()
    names_list <- file_names()
    if (length(data_list) >= 6) {
      tryCatch({
        generate_plot(data_list[[6]], 6, names_list[6], input$plot6_brush)
      }, error = function(e) {
        message("Error in plot 6: ", e$message)
        ggplot() + annotate("text", x = 0.5, y = 0.5, label = paste("Error:", e$message)) + theme_void()
      })
    } else {
      ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No data loaded for plot 6") + theme_void()
    }
  })
  
  output$plot7 <- renderPlot({
    req(loaded_data(), file_names())
    data_list <- plotData()
    names_list <- file_names()
    if (length(data_list) >= 7) {
      tryCatch({
        generate_plot(data_list[[7]], 7, names_list[7], input$plot7_brush)
      }, error = function(e) {
        message("Error in plot 7: ", e$message)
        ggplot() + annotate("text", x = 0.5, y = 0.5, label = paste("Error:", e$message)) + theme_void()
      })
    } else {
      ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No data loaded for plot 7") + theme_void()
    }
  })
  
  output$plot8 <- renderPlot({
    req(loaded_data(), file_names())
    data_list <- plotData()
    names_list <- file_names()
    if (length(data_list) >= 8) {
      tryCatch({
        generate_plot(data_list[[8]], 8, names_list[8], input$plot8_brush)
      }, error = function(e) {
        message("Error in plot 8: ", e$message)
        ggplot() + annotate("text", x = 0.5, y = 0.5, label = paste("Error:", e$message)) + theme_void()
      })
    } else {
      ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No data loaded for plot 8") + theme_void()
    }
  })
  
  tsne_data <- reactive({
    req(loaded_data(), input$tsne_file_select)
    data_list <- loaded_data()
    names_list <- file_names()
    idx <- match(input$tsne_file_select, names_list)
    if (!is.na(idx) && idx <= length(data_list)) {
      df <- data_list[[idx]]
      if (input$additive_gating && input$plot_type == "2D Plot") {
        gates <- filtered_gates()
        if (length(gates) >= idx && !is.null(gates[[idx]]) && nrow(gates[[idx]]) > 0) {
          df <- gates[[idx]]
        } else {
          if (input$x_channel %in% names(df) && input$y_channel %in% names(df)) {
            df <- df[df[[input$x_channel]] >= input$x_range[1] & df[[input$x_channel]] <= input$x_range[2] &
                     df[[input$y_channel]] >= input$y_range[1] & df[[input$y_channel]] <= input$y_range[2], ]
          }
        }
      } else if (input$plot_type == "2D Plot" && input$x_channel %in% names(df) && input$y_channel %in% names(df)) {
        df <- df[df[[input$x_channel]] >= input$x_range[1] & df[[input$x_channel]] <= input$x_range[2] &
                 df[[input$y_channel]] >= input$y_range[1] & df[[input$y_channel]] <= input$y_range[2], ]
      }
      df
    } else {
      NULL
    }
  })
  
  # Download handlers

output$downloadPlot <- downloadHandler(
  filename = function() {
    paste0("flow_cytometry_plots_", Sys.Date(), ".pdf")
  },
  content = function(file) {
    req(loaded_data(), file_names())
    
    pdf(file, width = 12, height = 8)
    
    data_list <- plotData()
    names_list <- file_names()
    
    for (i in seq_along(data_list)) {
      if (i <= 8) { # Only save first 8 plots
        # Get the appropriate brush/gate information for this plot
        current_brush <- NULL
        
        if (input$plot_type == "2D Plot") {
          if (input$individual_gates) {
            # Use individual gate information
            gates <- individual_gates()
            if (length(gates) >= i && !is.null(gates[[i]])) {
              current_brush <- gates[[i]]
            }
          } else if (input$persistent_brush) {
            # Use persistent brush selections
            saved_brushes <- brushed_selections()
            if (length(saved_brushes) >= i && !is.null(saved_brushes[[i]])) {
              current_brush <- saved_brushes[[i]]
            }
          }
        }
        
        plot_obj <- generate_plot(data_list[[i]], i, names_list[i], current_brush)
        print(plot_obj)
      }
    }
    
    dev.off()
  }
)
  
  output$downloadSummary <- downloadHandler(
    filename = function() {
      paste0("flow_cytometry_summary_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(loaded_data(), input$summary_file_select)
      names_list <- file_names()
      idx <- match(input$summary_file_select, names_list)
      data_list <- loaded_data()
      
      if (!is.na(idx) && idx <= length(data_list)) {
        df <- data_list[[idx]]
        summary_stats <- data.frame(
          Channel = names(df),
          Mean = sapply(df, function(x) mean(x, na.rm = TRUE)),
          Median = sapply(df, function(x) median(x, na.rm = TRUE)),
          SD = sapply(df, function(x) sd(x, na.rm = TRUE)),
          Min = sapply(df, function(x) min(x, na.rm = TRUE)),
          Max = sapply(df, function(x) max(x, na.rm = TRUE)),
          N = sapply(df, function(x) sum(!is.na(x)))
        )
        write.csv(summary_stats, file, row.names = FALSE)
      }
    }
  )
  
  output$downloadGated <- downloadHandler(
    filename = function() {
      paste0("gated_summary_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(input$plot_type == "2D Plot", input$summary_file_select, input$additive_gating)
      gates <- filtered_gates()
      names_list <- file_names()
      idx <- match(input$summary_file_select, names_list)
      data_list <- loaded_data()
      
      if (!is.na(idx) && idx <= length(gates) && idx <= length(data_list)) {
        gated_df <- gates[[idx]]
        full_df <- data_list[[idx]]
        
        if (!is.null(gated_df) && nrow(gated_df) > 0) {
          gated_stats <- data.frame(
            Channel = names(gated_df),
            Mean = sapply(gated_df, function(x) mean(x, na.rm = TRUE)),
            Median = sapply(gated_df, function(x) median(x, na.rm = TRUE)),
            SD = sapply(gated_df, function(x) sd(x, na.rm = TRUE)),
            Min = sapply(gated_df, function(x) min(x, na.rm = TRUE)),
            Max = sapply(gated_df, function(x) max(x, na.rm = TRUE)),
            N = sapply(gated_df, function(x) sum(!is.na(x))),
            PercentGated = round(100 * nrow(gated_df) / nrow(full_df), 2)
          )
          write.csv(gated_stats, file, row.names = FALSE)
        }
      }
    }
  )
}

shinyApp(ui, server)