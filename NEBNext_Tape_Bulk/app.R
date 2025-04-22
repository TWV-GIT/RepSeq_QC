# NEBNext TapeStation QC Shiny App - Deployment Version
# Modified to accept multiple input files simultaneously

# Load standard libraries
library(shiny)
library(ggplot2)
library(gridExtra)
library(grid)
library(dplyr)
library(shinydashboard)
library(DT)
library(reshape2)
library(plyr)
library(scales)
library(XML)
library(methods)
library(tools) # For file basename and extension functions

# Source local bioanalyzeR functions instead of loading the package
source("bioanalyzeR_local.R")

# Define UI
ui <- dashboardPage(
  dashboardHeader(title = "NEBNext TapeStation QC"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Input Files", tabName = "input", icon = icon("file-upload")),
      menuItem("Results", tabName = "results", icon = icon("chart-line")),
      menuItem("Settings", tabName = "settings", icon = icon("cog"))
    )
  ),
  dashboardBody(
    tabItems(
      tabItem(tabName = "input",
              fluidRow(
                box(
                  title = "Upload Files", width = 12, status = "primary", solidHeader = TRUE,
                  fileInput("xmlFiles", "Select XML Files", accept = ".xml", multiple = TRUE),
                  fileInput("csvFiles", "Select Electropherogram CSV Files", accept = ".csv", multiple = TRUE),
                  fileInput("annotationsFile", "Select the Annotations file (optional)", accept = ".csv"),
                  actionButton("matchFiles", "Match Files", class = "btn-info"),
                  hr(),
                  DTOutput("matchedFilesTable"),
                  hr(),
                  actionButton("process", "Process Files", class = "btn-success")
                )
              )
      ),
      tabItem(tabName = "settings",
              fluidRow(
                box(
                  title = "Parameters", width = 6, status = "warning", solidHeader = TRUE,
                  numericInput("dilution", "Dilution Factor:", 10, min = 1, max = 100),
                  numericInput("lowerBound", "Lower Bound (bp):", 600, min = 0, max = 1000),
                  numericInput("upperBound", "Upper Bound (bp):", 800, min = 0, max = 1000),
                  textInput("excludeWells", "Exclude Wells (comma separated):", ""),
                  numericInput("bcrRefline", "BCR Reference Line:", 635, min = 0, max = 1000),
                  numericInput("tcrRefline", "TCR Reference Line:", 660, min = 0, max = 1000),
                  numericInput("ighRefline", "IGH Reference Line:", 700, min = 0, max = 1000)
                ),
                box(
                  title = "Sample Annotations", width = 6, status = "warning", solidHeader = TRUE,
                  radioButtons("annotationSource", "Annotation Source:",
                               choices = list("Use headers from CSV" = "csv_headers",
                                              "Use custom annotations file" = "custom_file",
                                              "Use standard annotations" = "standard")),
                  downloadButton("downloadTemplate", "Download Annotation Template")
                )
              )
      ),
      tabItem(tabName = "results",
              fluidRow(
                box(
                  title = "Sample Selection", width = 12, status = "primary", solidHeader = TRUE,
                  uiOutput("sampleSelectionUI")
                )
              ),
              fluidRow(
                box(
                  title = "Sample Information", width = 12, status = "info", solidHeader = TRUE,
                  DTOutput("sampleTable")
                )
              ),
              fluidRow(
                box(
                  title = "Concentration Plot", width = 12, status = "success", solidHeader = TRUE,
                  plotOutput("concentrationPlot")
                )
              ),
              fluidRow(
                box(
                  title = "Molarity Plot", width = 12, status = "success", solidHeader = TRUE,
                  plotOutput("molarityPlot")
                )
              ),
              fluidRow(
                box(
                  width = 12,
                  downloadButton("downloadPDF", "Download PDF Report"),
                  downloadButton("downloadCSV", "Download Data as CSV")
                )
              )
      )
    )
  )
)

# Define server logic
server <- function(input, output, session) {
  
  # Add this near the top of your server function
  session$onSessionEnded(function() {
    cat("Session ended\n")
  })
  
  # Define standard annotations
  standard_annotations <- data.frame(
    well.number = c("B1", "C1", "D1", "E1", "F1", "G1", "H1", "A2", "B2", "C2", "D2", "E2", "F2", "G2", "H2"),
    sample.type = c("IGH", "BCR", "TCR", "IGH", "BCR", "TCR", "IGH", "BCR", "TCR", "IGH", "BCR", "TCR", "IGH", "BCR", "TCR")
  )
  
  # Reactive values to store data
  values <- reactiveValues(
    matchedFiles = data.frame(XML = character(), CSV = character(), RunID = character(), stringsAsFactors = FALSE),
    processedData = list(),
    combinedTable = NULL,
    selectedSamples = NULL,
    filteredElectrophoresis = NULL,
    conc_plot = NULL,
    mol_plot = NULL
  )
  
  # Match files when button is clicked
  observeEvent(input$matchFiles, {
    req(input$xmlFiles, input$csvFiles)
    
    # Extract file names and paths
    xml_files <- input$xmlFiles
    csv_files <- input$csvFiles
    
    xml_names <- basename(xml_files$name)
    csv_names <- basename(csv_files$name)
    
    # Prepare a data frame to store matched files
    matched_df <- data.frame(XML = character(), 
                             CSV = character(), 
                             XMLPath = character(),
                             CSVPath = character(),
                             RunID = character(),
                             stringsAsFactors = FALSE)
    
    # Match files by name
    for (i in 1:length(xml_names)) {
      xml_name <- xml_names[i]
      xml_path <- xml_files$datapath[i]
      run_id <- tools::file_path_sans_ext(xml_name)
      
      # Look for matching CSV files
      for (j in 1:length(csv_names)) {
        csv_name <- csv_names[j]
        csv_path <- csv_files$datapath[j]
        
        # Check if this CSV matches the XML (naming pattern: runID*Electropherogram.csv)
        if (grepl(paste0("^", run_id, ".*Electropherogram\\.csv$"), csv_name, ignore.case = TRUE)) {
          matched_df <- rbind(matched_df, data.frame(
            XML = xml_name,
            CSV = csv_name,
            XMLPath = xml_path,
            CSVPath = csv_path,
            RunID = run_id,
            stringsAsFactors = FALSE
          ))
          break
        }
      }
    }
    
    values$matchedFiles <- matched_df
    
    # Show notification about matching results
    if (nrow(matched_df) == 0) {
      showNotification("No matching file pairs found. Please check your file names.", type = "error")
    } else {
      showNotification(paste(nrow(matched_df), "file pairs matched successfully!"), type = "message")
    }
  })
  
  # Display matched files in a table
  output$matchedFilesTable <- renderDT({
    req(values$matchedFiles)
    if (nrow(values$matchedFiles) == 0) {
      return(NULL)
    }
    
    # Create a display version of the matched files table without the paths
    display_df <- values$matchedFiles[, c("RunID", "XML", "CSV")]
    datatable(display_df, options = list(pageLength = 10))
  })
  
  # Download annotation template
  output$downloadTemplate <- downloadHandler(
    filename = function() {
      "annotation_template.csv"
    },
    content = function(file) {
      write.csv(standard_annotations, file, row.names = FALSE)
    }
  )
  
  # Process files when button is clicked
  observeEvent(input$process, {
    req(values$matchedFiles)
    
    if (nrow(values$matchedFiles) == 0) {
      showNotification("No matched files to process. Please match files first.", type = "error")
      return()
    }
    
    # Display progress
    withProgress(message = 'Processing files...', value = 0, {
      
      # Clear previous processed data
      values$processedData <- list()
      values$combinedTable <- NULL
      
      # Calculate the progress step based on number of files
      files_count <- nrow(values$matchedFiles)
      step_size <- 1 / (files_count + 1)
      
      # Process each matched pair
      for (i in 1:files_count) {
        incProgress(step_size, detail = paste("Processing file", i, "of", files_count))
        
        xml_path <- values$matchedFiles$XMLPath[i]
        csv_path <- values$matchedFiles$CSVPath[i]
        run_id <- values$matchedFiles$RunID[i]
        
        cat("Processing", run_id, "\n")
        
        # Read TapeStation data
        tryCatch({
          electrophoresis <- read.tapestation(
            xml_path,
            csv_path,
            method = "hyman",
            extrapolate = FALSE
          )
          
          # Determine annotations
          if (input$annotationSource == "standard") {
            annotations <- standard_annotations
          } else if (input$annotationSource == "custom_file" && !is.null(input$annotationsFile)) {
            annotations <- read.csv(input$annotationsFile$datapath)
          } else if (input$annotationSource == "csv_headers") {
            # Extract sample names from CSV headers
            csv_data <- read.csv(csv_path, check.names = FALSE)
            
            # grab just the sample‑columns
            orig <- colnames(csv_data)[-1]
            
            # define a single regex with two capture groups:
            #  1) the well (letter+digits)
            #  2) whatever comes after the ": "
            pat <- "^([A-H][0-9]{1,2}):\\s*(.*)$"
            
            # apply it once to all names 
            m <- regexec(pat, orig)
            vals <- regmatches(orig, m)
            
            # extract into two vectors, falling back to "Unknown#" when there was no match
            well_numbers <- ifelse(
              lengths(vals) == 3,
              vapply(vals, `[`, character(1), 2),
              paste0("Unknown", seq_along(orig))
            )
            
            clean_sample_names <- ifelse(
              lengths(vals) == 3,
              vapply(vals, `[`, character(1), 3),
              orig
            )
            
            # Handle duplicate wells
            duplicate_wells <- well_numbers[duplicated(well_numbers)]
            if (length(duplicate_wells) > 0) {
              cat("Duplicate well numbers detected:", paste(duplicate_wells, collapse=", "), "\n")
              
              for (i in which(duplicated(well_numbers))) {
                well_numbers[i] <- paste0(well_numbers[i], "_dup", sum(well_numbers[1:i] == well_numbers[i]))
              }
            }
            
            # Determine sample types from sample names
            if (any(nchar(clean_sample_names) < 3)) {
              cat("Warning: Some sample names are less than 3 characters long\n")
              sample_types <- sapply(clean_sample_names, function(name) {
                if (nchar(name) < 3) {
                  return("UNK")  # Default type for short names
                } else {
                  return(substr(name, nchar(name) - 2, nchar(name)))
                }
              })
            } else {
              sample_types <- substr(clean_sample_names, nchar(clean_sample_names) - 2, nchar(clean_sample_names))
              clean_sample_names <- substr(clean_sample_names, 1, nchar(clean_sample_names) - 4)
            }
            
            annotations <- data.frame(
              well.number = well_numbers,
              sample.name = clean_sample_names,
              sample.type = sample_types,
              stringsAsFactors = FALSE
            )
          }
          
          # Create a unique run identifier for each well to avoid duplicates across runs
          # Add run ID prefix to well numbers in annotations to make them unique across runs
          annotations$unique_well <- paste(run_id, annotations$well.number, sep = "_")
          
          # Annotate electrophoresis
          an_electrophoresis <- annotate.electrophoresis(
            electrophoresis,
            annotations,
            header = TRUE,
            row.names = NULL,
            sep = "\t",
            stringsAsFactors = TRUE
          )
          
          # Add unique well identifier to annotated electrophoresis data
          an_electrophoresis$unique_well <- paste(run_id, an_electrophoresis$well.number, sep = "_")
          
          # Apply well exclusions
          if (input$excludeWells != "") {
            exclude_wells <- unlist(strsplit(input$excludeWells, ",\\s*"))
            f_electrophoresis <- subset(an_electrophoresis, !(well.number %in% exclude_wells))
          } else {
            f_electrophoresis <- an_electrophoresis
          }
          
          # Create table data
          table_data <- subset(an_electrophoresis[["samples"]], 
                               select = c('well.number', 'sample.name', 'sample.type'))
          
          # Add run ID to the table data and sample names
          table_data$run_id <- run_id
          table_data$unique_well <- paste(run_id, table_data$well.number, sep = "_")
          table_data$full_name <- paste(run_id, table_data$sample.name, sep = " - ")
          
          # Calculate molarity
          table_data$molarity <- input$dilution * round(integrate.custom(
            an_electrophoresis,
            lower.bound = input$lowerBound,
            upper.bound = input$upperBound,
            bound.variable = "length",
            sum.variable = "molarity"
          ), 0)
          
          # Calculate concentration
          table_data$concentration <- input$dilution * round(integrate.custom(
            an_electrophoresis,
            lower.bound = input$lowerBound,
            upper.bound = input$upperBound,
            bound.variable = "length",
            sum.variable = "concentration"
          ), 0)
          
          # Store all relevant data for this run
          values$processedData[[run_id]] <- list(
            electrophoresis = electrophoresis,
            annotations = annotations,
            an_electrophoresis = an_electrophoresis,
            f_electrophoresis = f_electrophoresis,
            table_data = table_data
          )
          
        }, error = function(e) {
          showNotification(paste("Error processing", run_id, ":", e$message), type = "error")
        })
      }
      
      # Combine all table data
      if (length(values$processedData) > 0) {
        combined_table <- do.call(rbind, lapply(values$processedData, function(x) x$table_data))
        values$combinedTable <- combined_table
        
        # Initialize selected samples to all samples
        values$selectedSamples <- combined_table$full_name
      }
      
      # Complete progress
      incProgress(step_size, detail = "Complete")
      
      # Notify user
      showNotification("Processing complete! Go to Results tab to select samples for plotting.", type = "message")
    })
  })
  
  # Generate sample selection UI dynamically
  output$sampleSelectionUI <- renderUI({
    req(values$combinedTable)
    
    # Group samples by run ID
    runs <- unique(values$combinedTable$run_id)
    selection_inputs <- list()
    
    # Overall select/deselect all button
    selection_inputs[[1]] <- fluidRow(
      column(6, actionButton("selectAllSamples", "Select All", class = "btn-info")),
      column(6, actionButton("selectNoneSamples", "Deselect All", class = "btn-info"))
    )
    
    # Add a horizontal rule
    selection_inputs[[2]] <- tags$hr()
    
    # Create checkboxes grouped by run
    for (i in seq_along(runs)) {
      run_id <- runs[i]
      samples <- subset(values$combinedTable, run_id == run_id)
      
      selection_inputs[[i+2]] <- box(
        title = paste("Run:", run_id), 
        status = "primary", 
        width = 12,
        collapsible = TRUE,
        checkboxGroupInput(
          inputId = paste0("samples_", run_id),
          label = NULL,
          choices = setNames(samples$full_name, samples$full_name),
          selected = samples$full_name
        )
      )
    }
    
    do.call(tagList, selection_inputs)
  })
  
  # Handle select all button
  observeEvent(input$selectAllSamples, {
    req(values$combinedTable)
    runs <- unique(values$combinedTable$run_id)
    
    for (run_id in runs) {
      samples <- subset(values$combinedTable, run_id == run_id)
      updateCheckboxGroupInput(
        session,
        inputId = paste0("samples_", run_id),
        selected = samples$full_name
      )
    }
  })
  
  # Handle deselect all button
  observeEvent(input$selectNoneSamples, {
    req(values$combinedTable)
    runs <- unique(values$combinedTable$run_id)
    
    for (run_id in runs) {
      updateCheckboxGroupInput(
        session,
        inputId = paste0("samples_", run_id),
        selected = character(0)
      )
    }
  })
  
  # Update selected samples when checkboxes change
  observe({
    req(values$combinedTable)
    runs <- unique(values$combinedTable$run_id)
    
    selected <- character(0)
    for (run_id in runs) {
      input_id <- paste0("samples_", run_id)
      if (!is.null(input[[input_id]])) {
        selected <- c(selected, input[[input_id]])
      }
    }
    
    values$selectedSamples <- selected
    
    # After updating selected samples, filter and combine electrophoresis data for plotting
    updateFilteredElectrophoresis()
  })
  
  # Function to update filtered electrophoresis based on selected samples
  updateFilteredElectrophoresis <- function() {
    req(values$selectedSamples, values$processedData)
    
    # Get selected sample info
    selected_table <- values$combinedTable[values$combinedTable$full_name %in% values$selectedSamples, ]
    
    # Create an empty combined electrophoresis dataset
    combined_electrophoresis <- NULL
    
    for (run_id in unique(selected_table$run_id)) {
      # Get the processed data for this run
      run_data <- values$processedData[[run_id]]
      if (is.null(run_data)) {
        next
      }
      
      # Get selected samples for this run
      run_selected <- selected_table[selected_table$run_id == run_id, ]
      
      # Get unique well identifiers for this run's selected samples
      selected_unique_wells <- run_selected$unique_well
      
      # Filter electrophoresis data for this run to include only selected samples
      # Use unique_well to ensure we're matching correctly across runs
      f_electro <- subset(run_data$an_electrophoresis, unique_well %in% selected_unique_wells)
      
      # Add run ID to sample names
      f_electro$original_sample.name <- f_electro$sample.name
      f_electro$sample.name <- paste(run_id, f_electro$sample.name, sep = " - ")
      
      # Combine with previous data
      if (is.null(combined_electrophoresis)) {
        combined_electrophoresis <- f_electro
      } else {
        # Combine correctly - this depends on the structure of your electrophoresis object
        combined_electrophoresis <- rbind(combined_electrophoresis, f_electro)
      }
    }
    
    values$filteredElectrophoresis <- combined_electrophoresis
    
    # Update plots
    createPlots()
  }
  
  # Function to create plots
  createPlots <- function() {
    req(values$filteredElectrophoresis)
    
    # Set up vertical reference lines
    vlines <- data.frame(
      sample.type = c("IGH", "BCR", "TCR"),
      xintercept = c(input$ighRefline, input$bcrRefline, input$tcrRefline),
      linetype = c("dashed", "dashed", "dashed"),
      color = c("black", "black", "black"),
      alpha = c(0.3, 0.3, 0.3)
    )
    vlines$text_label <- vlines$xintercept
    
    # Concentration plot
    p_conc <- qplot.electrophoresis(
      values$filteredElectrophoresis,
      x = "length",
      y = "concentration",
      color = sample.name,
      log = "",
      normalize = "none",
      facets = ~sample.type,
      margins = FALSE,
      scales = "fixed",
      geom = c("line", "area"),
      include.ladder = FALSE,
      include.markers = FALSE,
      lower.marker.spread = 10,
      xlim = c(200, 1000),
      ylim = c(NA, NA),
      show.peaks = "none",
      region.alpha = 0.2,
      area.alpha = 0.2,
      title = 'NEBNext IS (human) library size distribution for undiluted samples',
      xlab = NULL,
      ylab = NULL
    ) + theme_minimal() +
      theme(panel.spacing = unit(2, "lines")) +
      scale_color_discrete(name = "Sample")
    
    # Scale based on dilution
    conc_df <- p_conc[["data"]]
    conc_df <- conc_df %>% mutate(y.scaled = y.scaled * input$dilution)
    p_conc[["data"]] <- conc_df
    
    values$conc_plot <- p_conc + 
      geom_vline(data = vlines, aes(xintercept = xintercept), 
                 linetype = vlines$linetype, color = vlines$color, alpha = vlines$alpha) +
      geom_text(data = vlines, aes(x = xintercept+20, y = 0, label = text_label),
                vjust = 1.5, hjust = 0, size = 3.5)
    
    # Molarity plot
    p_mol <- qplot.electrophoresis(
      values$filteredElectrophoresis,
      x = "length",
      y = "molarity",
      color = sample.name,
      log = "",
      normalize = "none",
      facets = ~sample.type,
      margins = FALSE,
      scales = "fixed",
      geom = c("line", "area"),
      include.ladder = FALSE,
      include.markers = FALSE,
      lower.marker.spread = 10,
      xlim = c(200, 1000),
      ylim = c(NA, NA),
      show.peaks = "none",
      region.alpha = 0.2,
      area.alpha = 0.2,
      title = 'NEBNext IS (human) library size distribution for undiluted samples',
      xlab = NULL,
      ylab = NULL
    ) + theme_minimal() +
      theme(panel.spacing = unit(2, "lines")) +
      scale_color_discrete(name = "Sample")
    
    # Scale based on dilution
    mol_df <- p_mol[["data"]]
    mol_df <- mol_df %>% mutate(y.scaled = y.scaled * input$dilution)
    p_mol[["data"]] <- mol_df
    
    values$mol_plot <- p_mol + 
      geom_vline(data = vlines, aes(xintercept = xintercept), 
                 linetype = vlines$linetype, color = vlines$color, alpha = vlines$alpha) +
      geom_text(data = vlines, aes(x = xintercept+20, y = 0, label = text_label),
                vjust = 1.5, hjust = 0, size = 3.5)
  }
  
  # Render the sample table
  output$sampleTable <- renderDT({
    req(values$combinedTable, values$selectedSamples)
    
    # Filter the table to show only selected samples
    filtered_table <- values$combinedTable[values$combinedTable$full_name %in% values$selectedSamples, ]
    
    # Display columns
    display_cols <- c('run_id', 'well.number', 'sample.name', 'sample.type', 'concentration', 'molarity')
    datatable(filtered_table[, display_cols], options = list(pageLength = 15))
  })
  
  # Render the concentration plot
  output$concentrationPlot <- renderPlot({
    req(values$conc_plot)
    values$conc_plot
  }, height = 300)
  
  # Render the molarity plot
  output$molarityPlot <- renderPlot({
    req(values$mol_plot)
    values$mol_plot
  }, height = 300)
  
  # Download PDF report
  output$downloadPDF <- downloadHandler(
    filename = function() {
      "NEBNext_library_size_plot.pdf"
    },
    content = function(file) {
      req(values$combinedTable, values$selectedSamples, values$conc_plot, values$mol_plot)
      
      # Filter the table to show only selected samples
      filtered_table <- values$combinedTable[values$combinedTable$full_name %in% values$selectedSamples, ]
      
      # Create PDF
      pdf(file, width = 11, height = 8)
      
      # Summary page
      plot(1, type="n", xlab="", ylab="", xlim=c(0,10), ylim=c(0,10), axes=FALSE)
      lines <- c(
        paste0("NEBNext TapeStation QC Analysis"),
        paste0("Analysis date: ", format(Sys.time(), "%Y-%m-%d at %H:%M:%S")),
        "",
        paste0("Number of runs: ", length(unique(filtered_table$run_id))),
        paste0("Number of samples: ", nrow(filtered_table)),
        paste0("Wells excluded from plotting: ", input$excludeWells),
        paste0("Dilution factor: ", input$dilution),
        paste0("Molarity and concentration values are for fragments of size ", input$lowerBound, " to ", input$upperBound, " bp."),
        paste0("Molarity is expressed in pM. Concentration is expressed in pg/µl."),
        paste0("Values shown are for undiluted samples (considering ", input$dilution, "X dilution before analysis).")
      )
      text(x=5, y=seq(9, 5, by=-0.5), labels=lines, cex=1)
      
      # Data table
      grid.newpage()
      display_cols <- c('run_id', 'well.number', 'sample.name', 'sample.type', 'concentration', 'molarity')
      grid.table(filtered_table[, display_cols], rows = NULL)
      
      # Plots
      print(values$conc_plot)
      print(values$mol_plot)
      
      dev.off()
    }
  )
  
  # Download CSV data
  output$downloadCSV <- downloadHandler(
    filename = function() {
      "NEBNext_analysis_data.csv"
    },
    content = function(file) {
      req(values$combinedTable, values$selectedSamples)
      
      # Filter the table to show only selected samples
      filtered_table <- values$combinedTable[values$combinedTable$full_name %in% values$selectedSamples, ]
      
      # Remove the full_name and unique_well columns before saving
      filtered_table$full_name <- NULL
      filtered_table$unique_well <- NULL
      
      write.csv(filtered_table, file, row.names = FALSE)
    }
  )
}

# Run the application
shinyApp(ui = ui, server = server)