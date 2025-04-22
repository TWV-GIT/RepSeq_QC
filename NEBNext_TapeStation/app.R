# NEBNext TapeStation QC Shiny App - Deployment Version
# Modified to use local bioanalyzeR functions

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
                  fileInput("csvFile", "Select the Electropherogram .csv file", accept = ".csv"),
                  fileInput("xmlFile", "Select the .xml file", accept = ".xml"),
                  fileInput("annotationsFile", "Select the Annotations file (optional)", accept = ".csv"),
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
  
  # For file upload debugging
  observeEvent(input$csvFile, {
    cat("CSV file uploaded:", input$csvFile$name, "\n")
  })
  
  observeEvent(input$xmlFile, {
    cat("XML file uploaded:", input$xmlFile$name, "\n")
  })
  
  # Define standard annotations
  standard_annotations <- data.frame(
    well.number = c("B1", "C1", "D1", "E1", "F1", "G1", "H1", "A2", "B2", "C2", "D2", "E2", "F2", "G2", "H2"),
    sample.type = c("IGH", "BCR", "TCR", "IGH", "BCR", "TCR", "IGH", "BCR", "TCR", "IGH", "BCR", "TCR", "IGH", "BCR", "TCR")
  )
  
  # Empty reactive values to store processed data
  values <- reactiveValues(
    electrophoresis = NULL,
    annotations = NULL,
    an_electrophoresis = NULL,
    f_electrophoresis = NULL,
    table_data = NULL,
    conc_plot = NULL,
    mol_plot = NULL,
    run_id = NULL
  )
  
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
    req(input$xmlFile, input$csvFile)
    
    # Display progress
    withProgress(message = 'Processing files...', value = 0, {
      
      # Step 1: Read files
      incProgress(0.2, detail = "Reading files")
      
      # Get run ID from XML filename
      values$run_id <- gsub(".*/|\\.xml", "", input$xmlFile$name)
      
      # Read TapeStation data
      tryCatch({
        values$electrophoresis <- read.tapestation(
          input$xmlFile$datapath,
          input$csvFile$datapath,
          method = "hyman",
          extrapolate = FALSE
        )
      }, error = function(e) {
        showNotification(paste("Error reading TapeStation files:", e$message), type = "error")
        return(NULL)
      })
      
      # Step 2: Handle annotations
      incProgress(0.2, detail = "Processing annotations")
      
      if (input$annotationSource == "standard") {
        values$annotations <- standard_annotations
      } else if (input$annotationSource == "custom_file" && !is.null(input$annotationsFile)) {
        values$annotations <- read.csv(input$annotationsFile$datapath)
      } else if (input$annotationSource == "csv_headers") {
        # Extract sample names from CSV headers
        # read in
        csv_data <- read.csv(input$csvFile$datapath, check.names = FALSE)
        
        # grab just the sample‑columns
        orig <- colnames(csv_data)[-1]
        
        # define a single regex with two capture groups:
        #  1) the well (letter+digits)
        #  2) whatever comes after the “: ”
        pat <- "^([A-H][0-9]{1,2}):\\s*(.*)$"
        
        # apply it once to all names 
        m <- regexec(pat, orig)
        vals <- regmatches(orig, m)
        
        # extract into two vectors, falling back to “Unknown#” when there was no match
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
        
        # inspect
        well_numbers
        clean_sample_names
        
        
        # Before creating the annotations dataframe, check for duplicates
        duplicate_wells <- well_numbers[duplicated(well_numbers)]
        if (length(duplicate_wells) > 0) {
          # Log the duplicates
          cat("Duplicate well numbers detected:", paste(duplicate_wells, collapse=", "), "\n")
          
          # Option 1: Make them unique by appending a suffix
          for (i in which(duplicated(well_numbers))) {
            well_numbers[i] <- paste0(well_numbers[i], "_dup", sum(well_numbers[1:i] == well_numbers[i]))
          }
        }
        
        # Create annotations dataframe with corrected well numbers and clean sample names
        # Check if all sample names are at least 3 characters long
        if (any(nchar(clean_sample_names) < 3)) {
          cat("Warning: Some sample names are less than 3 characters long\n")
          # Handle this case - maybe use a default type
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
        
        values$annotations <- data.frame(
          well.number = well_numbers,
          sample.name = clean_sample_names,
          sample.type = sample_types,
          stringsAsFactors = FALSE
        )
      }
      
      # Step 3: Annotate electrophoresis
      incProgress(0.2, detail = "Annotating data")
      
      if (!is.null(values$electrophoresis) && !is.null(values$annotations)) {
        values$an_electrophoresis <- annotate.electrophoresis(
          values$electrophoresis,
          values$annotations,
          header = TRUE,
          row.names = NULL,
          sep = "\t",
          stringsAsFactors = TRUE
        )
        
        # Apply well exclusions
        # Apply well exclusions
        if (input$excludeWells != "") {
          exclude_wells <- unlist(strsplit(input$excludeWells, ",\\s*"))
          values$f_electrophoresis <- subset(values$an_electrophoresis, !(well.number %in% exclude_wells))
        } else {
          values$f_electrophoresis <- values$an_electrophoresis
        }
        
        # Create table data
        values$table_data <- subset(values$an_electrophoresis[["samples"]], 
                                    select = c('well.number', 'sample.name', 'sample.type'))
        
        # Calculate molarity
        values$table_data$molarity <- input$dilution * round(integrate.custom(
          values$an_electrophoresis,
          lower.bound = input$lowerBound,
          upper.bound = input$upperBound,
          bound.variable = "length",
          sum.variable = "molarity"
        ), 0)
        
        # Calculate concentration
        values$table_data$concentration <- input$dilution * round(integrate.custom(
          values$an_electrophoresis,
          lower.bound = input$lowerBound,
          upper.bound = input$upperBound,
          bound.variable = "length",
          sum.variable = "concentration"
        ), 0)
      }
      
      # Step 4: Create plots
      incProgress(0.2, detail = "Creating plots")
      
      if (!is.null(values$f_electrophoresis)) {
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
          values$f_electrophoresis,
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
          values$f_electrophoresis,
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
      
      # Complete progress
      incProgress(0.2, detail = "Complete")
      
      # Notify user
      showNotification("Processing complete! View results in the Results tab.", type = "message")
    })
  })
  
  # Render the sample table
  output$sampleTable <- renderDT({
    req(values$table_data)
    datatable(values$table_data, options = list(pageLength = 15))
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
      paste0(values$run_id, "_library_size_plot.pdf")
    },
    content = function(file) {
      req(values$table_data, values$conc_plot, values$mol_plot)
      
      # Create PDF
      pdf(file, width = 11, height = 8)
      
      # Summary page
      plot(1, type="n", xlab="", ylab="", xlim=c(0,10), ylim=c(0,10), axes=FALSE)
      lines <- c(
        paste0("NEBNext TapeStation QC Analysis"),
        paste0("Run ID: ", values$run_id),
        paste0("Analysis date: ", format(Sys.time(), "%Y-%m-%d at %H:%M:%S")),
        "",
        paste0("Wells excluded from plotting: ", input$excludeWells),
        paste0("Dilution factor: ", input$dilution),
        paste0("Molarity and concentration values are for fragments of size ", input$lowerBound, " to ", input$upperBound, " bp."),
        paste0("Molarity is expressed in pM. Concentration is expressed in pg/µl."),
        paste0("Values shown are for undiluted samples (considering ", input$dilution, "X dilution before analysis).")
      )
      text(x=5, y=seq(9, 5, by=-0.5), labels=lines, cex=1)
      
      # Data table
      grid.newpage()
      grid.table(values$table_data, rows = NULL)
      
      # Plots
      print(values$conc_plot)
      print(values$mol_plot)
      
      dev.off()
    }
  )
  
  # Download CSV data
  output$downloadCSV <- downloadHandler(
    filename = function() {
      paste0(values$run_id, "_analysis_data.csv")
    },
    content = function(file) {
      req(values$table_data)
      write.csv(values$table_data, file, row.names = FALSE)
    }
  )
}

# Run the application
shinyApp(ui = ui, server = server)