library(shiny)
library(reshape2)
library(dplyr)
library(ggplot2)
library(pracma)
library(DT)
library(shinyjs)

# Define UI --------------------------------------------------------------------
ui <- fluidPage(
  useShinyjs(),
  titlePanel("NEBNext IS qPCR Absolute Quantification"),
  
  sidebarLayout(
    sidebarPanel(
      h4("Settings"),
      fileInput("amp_results", "Select Quantification Amplification Results (CSV)", accept = ".csv"),
      fileInput("quant_summary", "Select Quantification Summary (CSV)", accept = ".csv"),
      hr(),
      numericInput("cycle_cutoff", "Cycle cutoff (upper):", value = 30),
      numericInput("lower_cycle_cutoff", "Cycle cutoff (lower):", value = 0),
      numericInput("sample_detect_threshold", "Sample detection threshold:", value = 500),
      textInput("wavelength", "Wavelength:", value = "SYBR"),
      numericInput("abs_threshold", "Absolute threshold:", value = 188),
      uiOutput("baseline_range"),
      hr(),
      textInput("plot_title", "Plot title:", value = "NEBNext IS library quantification plot"),
      hr(),
      textAreaInput("filter_names", "Samples to filter (comma separated):", value = "NTC,"),
      textAreaInput("ladder_names", "Standard sample names (comma separated):", value = "100,10,1,0.1"),
      textInput("ladder_unit", "Standard unit:", value = "nM"),
      hr(),
      numericInput("ladder_amplicon_size", "Standard amplicon size:", value = 399),
      numericInput("igh_amplicon_size", "IGH amplicon size:", value = 700),
      numericInput("bcr_amplicon_size", "BCR amplicon size:", value = 635),
      numericInput("tcr_amplicon_size", "TCR amplicon size:", value = 660),
      numericInput("sample_dilution", "Sample dilution factor:", value = 2000),
      hr(),
      textAreaInput("highlight_names", "Samples to highlight (comma separated):", value = "A"),
      hr(),
      textAreaInput("well_filter", "Wells to filter (comma separated):", value = ""),
      hr(),
      actionButton("process_btn", "Process Data", class = "btn-primary"),
      hr(),
      downloadButton("download_plots", "Download Plots (PDF)"),
      downloadButton("download_tables", "Download Tables (TSV)")
    ),
    
    mainPanel(
      tabsetPanel(
        id = "main_tabs",
        tabPanel("Sample Renaming", 
                 h4("Rename Samples"),
                 p("Here you can rename your samples before plotting. The original names are shown on the left."),
                 DT::dataTableOutput("rename_table"),
                 actionButton("apply_names_btn", "Apply New Names", class = "btn-success")
        ),
        tabPanel("Amplification Curves",
                 plotOutput("plot_log", height = "500px"),
                 plotOutput("plot_linear", height = "500px"),
                 plotOutput("plot_selected", height = "500px")
        ),
        tabPanel("Standard Regression",
                 plotOutput("plot_regression", height = "500px")
        ),
        tabPanel("Concentrations",
                 plotOutput("plot_concentrations", height = "500px"),
                 DT::dataTableOutput("concentration_table")
        ),
        tabPanel("Raw Data",
                 DT::dataTableOutput("raw_data_table")
        )
      )
    )
  )
)

# Define server ----------------------------------------------------------------
server <- function(input, output, session) {
  
  # Initialize reactive values
  rv <- reactiveValues(
    data_filtered = NULL,
    average_sd = NULL,
    df_summary = NULL,
    data_filtered_meta = NULL,
    sample_names = NULL,
    original_to_new_names = NULL,
    has_processed = FALSE
  )
  
  # Dynamic UI for baseline range
  output$baseline_range <- renderUI({
    sliderInput("baseline", "Baseline range:", 
                min = input$lower_cycle_cutoff, 
                max = 10, 
                value = c(input$lower_cycle_cutoff, 3))
  })
  
  # Process data when button is clicked
  observeEvent(input$process_btn, {
    req(input$amp_results, input$quant_summary)
    
    # Load datasets
    RFU <- read.csv(input$amp_results$datapath, header = TRUE, check.names = FALSE)
    RFU <- RFU[, names(RFU) != ""]
    RFU$Cycle <- as.numeric(rownames(RFU))
    RFU_long <- melt(RFU, id.vars = "Cycle", variable.name = "Well", value.name = "fluorescence")
    RFU_long$Well <- sprintf("%s%02d", substr(RFU_long$Well, 1, 1), as.integer(substr(RFU_long$Well, 2, nchar(as.character(RFU_long$Well)))))
    
    Cq <- read.csv(input$quant_summary$datapath, header = TRUE, check.names = FALSE)
    Cq <- Cq[, !(names(Cq) == "")]
    
    merged_table <- merge(RFU_long, Cq[, c("Well", "Sample", "Target")], by = "Well")
    merged_table <- merged_table[order(as.numeric(merged_table$Sample), merged_table$Cycle), ]
    
    # Parse inputs
    filternames <- unlist(strsplit(input$filter_names, ",\\s*"))
    laddernames <- unlist(strsplit(input$ladder_names, ",\\s*"))
    highlightnames <- unlist(strsplit(input$highlight_names, ",\\s*"))
    well_filter <- unlist(strsplit(input$well_filter, ",\\s*"))
    
    # Store unique sample names for renaming
    rv$sample_names <- unique(merged_table$Sample)
    rv$original_to_new_names <- data.frame(
      Original = rv$sample_names,
      New = rv$sample_names,
      stringsAsFactors = FALSE
    )
    
    # Calculations
    data_filtered <- merged_table %>%
      group_by(Sample, Target) %>%
      filter(max(fluorescence) > input$sample_detect_threshold)
    
    if (length(well_filter) > 0 && well_filter[1] != "") {
      data_filtered <- data_filtered[!data_filtered$Well %in% well_filter, ]
    }
    
    data_filtered <- data_filtered[!data_filtered$Sample %in% filternames, ]
    data_filtered <- subset(data_filtered, Cycle < input$cycle_cutoff)
    data_filtered <- subset(data_filtered, Cycle > input$lower_cycle_cutoff)
    
    # Baseline calculation
    baseline <- seq(input$baseline[1], input$baseline[2])
    df_baseline <- data_filtered[data_filtered$Cycle %in% baseline, ]
    baseline_fluorescence <- aggregate(df_baseline$fluorescence, by = list(Well = df_baseline$Well), FUN = mean)
    colnames(baseline_fluorescence) <- c("Well", "BaselineFluorescence")
    data_filtered <- merge(data_filtered, baseline_fluorescence, by = "Well")
    data_filtered$AdjustedFluorescence <- data_filtered$fluorescence - data_filtered$BaselineFluorescence
    
    # Calculate Ct values
    data_filtered <- data_filtered %>%
      group_by(Well, Target) %>%
      mutate(Ct = approx(AdjustedFluorescence, Cycle, xout = input$abs_threshold)$y)
    
    average_ct_values <- data_filtered %>%
      group_by(Sample, Target) %>%
      summarize(Avg_Ct = mean(Ct, na.rm = TRUE))
    
    sd_ct_values <- data_filtered %>%
      group_by(Sample, Target) %>%
      summarize(sd = sd(Ct))
    
    average_sd <- merge(average_ct_values, sd_ct_values[, c("Sample", "Target", "sd")], by = c("Sample", "Target"))
    average_sd <- average_sd %>%
      mutate(ModifiedSampleName = paste0(Sample, " ", Target, " (", round(Avg_Ct, 2), "±", round(sd, 3), ")"))
    
    data_filtered <- left_join(data_filtered, average_sd, by = c("Sample", "Target"))
    data_filtered_select <- data_filtered[data_filtered$Sample %in% highlightnames, ]
    average_ct_values_select <- average_ct_values[average_ct_values$Sample %in% highlightnames, ]
    
    # Standard regression
    data_filtered_meta <- data_filtered %>% 
      select(-fluorescence, -Cycle, -AdjustedFluorescence) %>% distinct()
    
    data_ladder <- data_filtered_meta[data_filtered_meta$Sample %in% laddernames, ]
    data_ladder$Sample <- as.numeric(data_ladder$Sample)
    
    model <- lm(Ct ~ log(Sample, 10), data = data_ladder)
    r_squared <- summary(model)$r.squared
    slope <- summary(model)$coefficients[2, 1]
    intercept <- summary(model)$coefficients[1, 1]
    efficiency <- (10^(-1/slope)) - 1
    
    data_filtered_meta$avg_ct_fit <- 10^((data_filtered_meta$Avg_Ct - intercept) / slope)
    data_filtered_meta$undiluted <- (10^((data_filtered_meta$Ct - intercept) / slope)) * input$sample_dilution
    data_filtered_meta <- data_filtered_meta %>% mutate(concentration = 
                                                          ifelse(Target == "IGH", undiluted * (input$ladder_amplicon_size / input$igh_amplicon_size),
                                                                 ifelse(Target == "BCR", undiluted * (input$ladder_amplicon_size / input$bcr_amplicon_size),
                                                                        ifelse(Target == "TCR", undiluted * (input$ladder_amplicon_size / input$tcr_amplicon_size), undiluted))))
    
    df_summary <- data_filtered_meta[!data_filtered_meta$Sample %in% laddernames, ] %>%
      group_by(Sample, Target) %>%
      summarise(
        mean_concentration = mean(concentration / 1000, na.rm = TRUE),
        se_concentration = sd(concentration / 1000, na.rm = TRUE) / sqrt(n())
      )
    
    # Store processed data
    rv$data_filtered <- data_filtered
    rv$average_sd <- average_sd
    rv$data_filtered_meta <- data_filtered_meta
    rv$df_summary <- df_summary
    rv$has_processed <- TRUE
    
    updateTabsetPanel(session, "main_tabs", selected = "Sample Renaming")
  })
  
  # Sample renaming table
  output$rename_table <- DT::renderDataTable({
    req(rv$original_to_new_names)
    DT::datatable(rv$original_to_new_names, 
                  editable = TRUE,
                  options = list(pageLength = 25, searching = FALSE, ordering = FALSE))
  })
  
  # Handle edited cells in renaming table
  observeEvent(input$rename_table_cell_edit, {
    info <- input$rename_table_cell_edit
    if(info$col == 2){
      rv$original_to_new_names[info$row, info$col] <- info$value
    }
  })
  
  # Apply new names button
  observeEvent(input$apply_names_btn, {
    req(rv$has_processed, rv$original_to_new_names)
    map_names <- function(old_name) {
      idx <- match(old_name, rv$original_to_new_names$Original)
      if (!is.na(idx)) return(rv$original_to_new_names$New[idx])
      else return(old_name)
    }
    if (!is.null(rv$data_filtered)) {
      rv$data_filtered$Sample_Original <- rv$data_filtered$Sample
      rv$data_filtered$Sample <- sapply(rv$data_filtered$Sample, map_names)
      
      rv$average_sd$Sample_Original <- rv$average_sd$Sample
      rv$average_sd$Sample <- sapply(rv$average_sd$Sample, map_names)
      rv$average_sd$ModifiedSampleName <- paste0(
        rv$average_sd$Sample, " ", 
        rv$average_sd$Target, " (", 
        round(rv$average_sd$Avg_Ct, 2), "±", 
        round(rv$average_sd$sd, 3), ")"
      )
      
      rv$data_filtered <- rv$data_filtered %>% 
        select(-ModifiedSampleName) %>% 
        left_join(rv$average_sd[, c("Sample", "Target", "ModifiedSampleName")], by = c("Sample", "Target"))
      
      rv$data_filtered_meta$Sample_Original <- rv$data_filtered_meta$Sample
      rv$data_filtered_meta$Sample <- sapply(rv$data_filtered_meta$Sample, map_names)
      
      rv$df_summary$Sample_Original <- rv$df_summary$Sample
      rv$df_summary$Sample <- sapply(rv$df_summary$Sample, map_names)
    }
    updateTabsetPanel(session, "main_tabs", selected = "Amplification Curves")
  })
  
  # Generate plots in main panel
  output$plot_log <- renderPlot({
    req(rv$data_filtered, rv$has_processed)
    data_filtered <- rv$data_filtered
    average_ct_values <- data_filtered %>%
      group_by(Sample, Target) %>%
      summarize(Avg_Ct = mean(Ct, na.rm = TRUE))
    
    ggplot(data_filtered, aes(x = Cycle, y = AdjustedFluorescence, group = Well, color = ModifiedSampleName)) +
      geom_line() +
      geom_hline(yintercept = input$abs_threshold, linetype = "dashed", color = "black") +
      geom_vline(data = average_ct_values, aes(xintercept = Avg_Ct), linetype = "dashed", color = "grey") +
      scale_x_continuous(breaks = seq(min(data_filtered$Cycle), max(data_filtered$Cycle), by = 1)) +
      scale_y_log10(breaks = round(logspace(log10(1), log10(max(data_filtered$fluorescence)), 20), 0)) +
      labs(x = "Cycle", y = paste0("Fluorescence (", input$wavelength, ")"),
           title = input$plot_title,
           subtitle = paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " - Log-linear plot")) +
      theme_minimal() +
      theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) +
      annotate("text", x = 1, y = input$abs_threshold, label = "Ct threshold", 
               hjust = 0.1, vjust = -0.2, color = "black") +
      scale_color_discrete(name = "Samples (Ct val.±SD)")
  })
  
  output$plot_linear <- renderPlot({
    req(rv$data_filtered, rv$has_processed)
    data_filtered <- rv$data_filtered
    average_ct_values <- data_filtered %>%
      group_by(Sample, Target) %>%
      summarize(Avg_Ct = mean(Ct, na.rm = TRUE))
    
    ggplot(data_filtered, aes(x = Cycle, y = AdjustedFluorescence, group = Well, color = ModifiedSampleName)) +
      geom_line() +
      geom_hline(yintercept = input$abs_threshold, linetype = "dashed", color = "black") +
      geom_vline(data = average_ct_values, aes(xintercept = Avg_Ct), linetype = "dashed", color = "grey") +
      scale_x_continuous(breaks = seq(min(data_filtered$Cycle), max(data_filtered$Cycle), by = 1)) +
      scale_y_continuous(breaks = seq(0, max(data_filtered$fluorescence), 
                                      by = ceiling(max(data_filtered$fluorescence)/10))) +
      labs(x = "Cycle", y = paste0("Fluorescence (", input$wavelength, ")"),
           title = input$plot_title,
           subtitle = paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " - Linear-linear plot")) +
      theme_minimal() +
      theme(panel.grid.minor.x = element_blank()) +
      annotate("text", x = 1, y = input$abs_threshold, label = "Ct threshold", 
               hjust = 0.1, vjust = -0.2, color = "black") +
      scale_color_discrete(name = "Samples (Ct val.±SD)")
  })
  
  output$plot_selected <- renderPlot({
    req(rv$data_filtered, rv$has_processed)
    highlightnames <- unlist(strsplit(input$highlight_names, ",\\s*"))
    mapped_highlightnames <- sapply(highlightnames, function(name) {
      idx <- match(name, rv$original_to_new_names$Original)
      if (!is.na(idx)) return(rv$original_to_new_names$New[idx])
      else return(name)
    })
    data_filtered <- rv$data_filtered
    data_filtered_select <- data_filtered[data_filtered$Sample %in% mapped_highlightnames, ]
    
    if (nrow(data_filtered_select) == 0) {
      return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No samples selected for highlighting") + theme_void())
    }
    
    average_ct_values_select <- data_filtered_select %>%
      group_by(Sample, Target) %>%
      summarize(Avg_Ct = mean(Ct, na.rm = TRUE))
    
    ggplot(data_filtered_select, aes(x = Cycle, y = AdjustedFluorescence, group = Well, color = ModifiedSampleName)) +
      geom_line() +
      geom_hline(yintercept = input$abs_threshold, linetype = "dashed", color = "black") +
      geom_vline(data = average_ct_values_select, aes(xintercept = Avg_Ct), linetype = "dashed", color = "grey") +
      scale_x_continuous(breaks = seq(min(data_filtered$Cycle), max(data_filtered$Cycle), by = 1)) +
      scale_y_continuous(breaks = seq(0, max(data_filtered$fluorescence), 
                                      by = ceiling(max(data_filtered$fluorescence)/10))) +
      labs(x = "Cycle", y = paste0("Fluorescence (", input$wavelength, ")"),
           title = paste(input$plot_title, 'of selected samples'),
           subtitle = paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " - Linear-linear plot")) +
      theme_minimal() +
      theme(panel.grid.minor.x = element_blank()) +
      annotate("text", x = 1, y = input$abs_threshold, label = "Ct threshold", 
               hjust = 0.1, vjust = -0.2, color = "black") +
      scale_color_discrete(name = "Samples (Ct val.±SD)")
  })
  
  output$plot_regression <- renderPlot({
    req(rv$data_filtered_meta, rv$has_processed)
    laddernames <- unlist(strsplit(input$ladder_names, ",\\s*"))
    data_filtered_meta <- rv$data_filtered_meta
    data_ladder <- data_filtered_meta[data_filtered_meta$Sample %in% laddernames, ]
    data_ladder$Sample <- as.numeric(data_ladder$Sample_Original)
    
    model <- lm(Ct ~ log(Sample, 10), data = data_ladder)
    r_squared <- summary(model)$r.squared
    slope <- summary(model)$coefficients[2, 1]
    intercept <- summary(model)$coefficients[1, 1]
    efficiency <- (10^(-1/slope)) - 1
    
    ggplot(data_ladder, aes(x = Sample, y = Ct)) +
      geom_smooth(method = "lm", col = "grey") +
      geom_point() +
      theme_minimal() +
      scale_x_log10() +
      labs(title = "Linear regression of standards",
           subtitle = paste("R-squared: ", round(r_squared, 3), "; Slope: ", round(slope, 2),
                            "; Efficiency: ", round(efficiency, 3) * 100, "%"),
           x = paste0("Concentration (", input$ladder_unit, ")"), y = "Ct-value") +
      geom_point(data = data_filtered_meta %>% filter(Target != "ladder"), 
                 aes(x = avg_ct_fit, y = Avg_Ct, color = Target), position = "identity")
  })
  
  output$plot_concentrations <- renderPlot({
    req(rv$df_summary, rv$has_processed)
    df_summary <- rv$df_summary
    ggplot(df_summary, aes(x = Sample, y = mean_concentration, fill = Target)) +
      geom_bar(stat = "identity", position = position_dodge()) +
      geom_errorbar(aes(ymin = mean_concentration - se_concentration, ymax = mean_concentration + se_concentration),
                    width = 0.2, position = position_dodge(0.9)) +
      labs(title = 'Absolute concentration of libraries ± SE',
           subtitle = paste0('Sample dilution = ', input$sample_dilution, "X;\nLibrary amplicon size of IGH / BCR / TCR: ", 
                             input$igh_amplicon_size, " / ", input$bcr_amplicon_size, " / ", input$tcr_amplicon_size,
                             "\nStandard amplicon size: ", input$ladder_amplicon_size),
           x = "Sample", y = paste0("Concentration (", input$ladder_unit, ")")) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })
  
  # Data tables
  output$concentration_table <- DT::renderDataTable({
    req(rv$df_summary)
    DT::datatable(rv$df_summary, options = list(pageLength = 25, scrollX = TRUE))
  })
  
  output$raw_data_table <- DT::renderDataTable({
    req(rv$data_filtered_meta)
    DT::datatable(rv$data_filtered_meta, options = list(pageLength = 25, scrollX = TRUE))
  })
  
  # Download handler for plots (PDF)
  output$download_plots <- downloadHandler(
    filename = function() {
      paste0("library_quantification_plots_", format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), ".pdf")
    },
    content = function(file) {
      laddernames <- unlist(strsplit(input$ladder_names, ",\\s*"))
      filternames <- unlist(strsplit(input$filter_names, ",\\s*"))
      well_filter <- unlist(strsplit(input$well_filter, ",\\s*"))
      highlightnames <- unlist(strsplit(input$highlight_names, ",\\s*"))
      
      mapped_highlightnames <- sapply(highlightnames, function(name) {
        idx <- match(name, rv$original_to_new_names$Original)
        if (!is.na(idx)) return(rv$original_to_new_names$New[idx])
        else return(name)
      })
      
      data_filtered <- rv$data_filtered
      data_filtered_meta <- rv$data_filtered_meta
      df_summary <- rv$df_summary
      
      average_ct_values <- data_filtered %>% group_by(Sample, Target) %>% 
        summarize(Avg_Ct = mean(Ct, na.rm = TRUE))
      
      data_filtered_select <- data_filtered[data_filtered$Sample %in% mapped_highlightnames, ]
      average_ct_values_select <- data_filtered_select %>% group_by(Sample, Target) %>% 
        summarize(Avg_Ct = mean(Ct, na.rm = TRUE))
      
      data_ladder <- data_filtered_meta[data_filtered_meta$Sample_Original %in% laddernames, ]
      data_ladder$Sample <- as.numeric(data_ladder$Sample_Original)
      
      model <- lm(Ct ~ log(Sample, 10), data = data_ladder)
      r_squared <- summary(model)$r.squared
      slope <- summary(model)$coefficients[2, 1]
      intercept <- summary(model)$coefficients[1, 1]
      efficiency <- (10^(-1/slope)) - 1
      
      # Build plots
      p_log <- ggplot(data_filtered, aes(x = Cycle, y = AdjustedFluorescence, group = Well, color = ModifiedSampleName)) +
        geom_line() +
        geom_hline(yintercept = input$abs_threshold, linetype = "dashed", color = "black") +
        geom_vline(data = average_ct_values, aes(xintercept = Avg_Ct), linetype = "dashed", color = "grey") +
        scale_x_continuous(breaks = seq(min(data_filtered$Cycle), max(data_filtered$Cycle), by = 1)) +
        scale_y_log10(breaks = round(logspace(log10(1), log10(max(data_filtered$fluorescence)), 20), 0)) +
        labs(x = "Cycle", y = paste0("Fluorescence (", input$wavelength, ")"), title = input$plot_title,
             subtitle = paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " - Log-linear plot")) +
        theme_minimal() +
        theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) +
        annotate("text", x = 1, y = input$abs_threshold, label = "Ct threshold", hjust = 0.1, vjust = -0.2, color = "black") +
        scale_color_discrete(name = "Samples (Ct val.±SD)")
      
      p_linear <- ggplot(data_filtered, aes(x = Cycle, y = AdjustedFluorescence, group = Well, color = ModifiedSampleName)) +
        geom_line() +
        geom_hline(yintercept = input$abs_threshold, linetype = "dashed", color = "black") +
        geom_vline(data = average_ct_values, aes(xintercept = Avg_Ct), linetype = "dashed", color = "grey") +
        scale_x_continuous(breaks = seq(min(data_filtered$Cycle), max(data_filtered$Cycle), by = 1)) +
        scale_y_continuous(breaks = seq(0, max(data_filtered$fluorescence), by = ceiling(max(data_filtered$fluorescence)/10))) +
        labs(x = "Cycle", y = paste0("Fluorescence (", input$wavelength, ")"), title = input$plot_title,
             subtitle = paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " - Linear-linear plot")) +
        theme_minimal() +
        theme(panel.grid.minor.x = element_blank()) +
        annotate("text", x = 1, y = input$abs_threshold, label = "Ct threshold", hjust = 0.1, vjust = -0.2, color = "black") +
        scale_color_discrete(name = "Samples (Ct val.±SD)")
      
      p_s_linear <- if(nrow(data_filtered_select) > 0) {
        ggplot(data_filtered_select, aes(x = Cycle, y = AdjustedFluorescence, group = Well, color = ModifiedSampleName)) +
          geom_line() +
          geom_hline(yintercept = input$abs_threshold, linetype = "dashed", color = "black") +
          geom_vline(data = average_ct_values_select, aes(xintercept = Avg_Ct), linetype = "dashed", color = "grey") +
          scale_x_continuous(breaks = seq(min(data_filtered$Cycle), max(data_filtered$Cycle), by = 1)) +
          scale_y_continuous(breaks = seq(0, max(data_filtered$fluorescence), by = ceiling(max(data_filtered$fluorescence)/10))) +
          labs(x = "Cycle", y = paste0("Fluorescence (", input$wavelength, ")"),
               title = paste(input$plot_title, 'of selected samples'),
               subtitle = paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " - Linear-linear plot")) +
          theme_minimal() +
          theme(panel.grid.minor.x = element_blank()) +
          annotate("text", x = 1, y = input$abs_threshold, label = "Ct threshold", hjust = 0.1, vjust = -0.2, color = "black") +
          scale_color_discrete(name = "Samples (Ct val.±SD)")
      } else {
        NULL
      }
      
      reg <- ggplot(data_ladder, aes(x = Sample, y = Ct)) +
        geom_smooth(method = "lm", col = "grey") +
        geom_point() +
        theme_minimal() +
        scale_x_log10() +
        labs(title = "Linear regression of standards",
             subtitle = paste("R-squared: ", round(r_squared, 3), "; Slope: ", round(slope, 2),
                              "; Efficiency: ", round(efficiency, 3) * 100, "%"),
             x = paste0("Concentration (", input$ladder_unit, ")"), y = "Ct-value") +
        geom_point(data = data_filtered_meta %>% filter(Target != "ladder"), 
                   aes(x = avg_ct_fit, y = Avg_Ct, color = Target), position = "identity")
      
      con <- ggplot(df_summary, aes(x = Sample, y = mean_concentration, fill = Target)) +
        geom_bar(stat = "identity", position = position_dodge()) +
        geom_errorbar(aes(ymin = mean_concentration - se_concentration, ymax = mean_concentration + se_concentration),
                      width = 0.2, position = position_dodge(0.9)) +
        labs(title = 'Absolute concentration of libraries ± SE',
             subtitle = paste0('Sample dilution = ', input$sample_dilution, "X;\nLibrary amplicon size of IGH / BCR / TCR: ", 
                               input$igh_amplicon_size, " / ", input$bcr_amplicon_size, " / ", input$tcr_amplicon_size,
                               "\nStandard amplicon size: ", input$ladder_amplicon_size),
             x = "Sample", y = paste0("Concentration (", input$ladder_unit, ")")) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      
      # Create PDF with header info and all plots
      pdf(file, width = 8, height = 6)
      plot(1, type = "n", xlab = "", ylab = "", xlim = c(0, 10), ylim = c(0, 10), axes = FALSE)
      header_text <- c(
        "PDF generated with the NEBNext IS qPCR Absolute Quantification Shiny App.",
        paste0("Amplification input file: ", input$amp_results$name),
        paste0("Quantification summary file: ", input$quant_summary$name),
        paste0("Filtered samples: ", paste(filternames, collapse = ", ")),
        paste0("Filtered wells: ", paste(well_filter, collapse = ", ")),
        paste0("Threshold value: ", input$abs_threshold),
        paste0("Baseline range: ", input$lower_cycle_cutoff, " to ", input$baseline[2], " cycles"),
        paste0("PDF generated on ", format(Sys.time(), "%Y-%m-%d at %H:%M:%S"))
      )
      text(x = 5, y = seq(10, 1, by = -1.2), labels = header_text, cex = 1)
      print(p_log)
      print(p_linear)
      if (!is.null(p_s_linear)) print(p_s_linear)
      print(reg)
      print(con)
      dev.off()
    }
  )
  
  # Download handler for tables (ZIP file containing TSVs)
  output$download_tables <- downloadHandler(
    filename = function() {
      paste0("library_quantification_tables_", format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), ".zip")
    },
    content = function(file) {
      tmpdir <- tempdir()
      file1 <- file.path(tmpdir, paste0(format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), "_Ct_values_", input$abs_threshold, "_threshold.tsv"))
      file2 <- file.path(tmpdir, paste0(format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), "_Library-concentrations_", input$abs_threshold, "_threshold.tsv"))
      write.table(rv$average_sd, file = file1, sep = "\t", row.names = FALSE)
      write.table(rv$df_summary, file = file2, sep = "\t", row.names = FALSE)
      utils::zip(zipfile = file, files = c(file1, file2))
    },
    contentType = "application/zip"
  )
}

shinyApp(ui, server)
