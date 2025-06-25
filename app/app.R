# ==============================================================================
#
# Shiny App for Displaying iCLIP Data from "Imp Syp iCLIP" Sci Adv Publication
#
# Author: Jeffrey Y Lee
# Date: June 2025
#
# Description:
# This Shiny application provides an interactive platform to explore the datasets
# associated with the "Imp Syp iCLIP" Sci Adv publication. It includes an
# introduction to the study, a searchable table of the full dataset, and a
# tool to explore specific gene binding sites.
#
#
# --- DIRECTORY STRUCTURE (for shinylive) ---
#
# /app_root/
# |-- app.R         (this file)
# |-- data/
# |   |-- full_data.rds
# |   |-- target_gene_names.RDS  
# |   |-- plots/                 
# |   |   |-- pros.RDS
# |   |   |-- GeneB.RDS
# |   |   |-- ...
# |
# |-- www/          
# |   |-- science_advances_logo.png 
# |   |-- intro_image_1.png         
# |   |-- intro_image_2.png         
# |   |-- intro_image_3.png         
#
# ==============================================================================


# 1. LOAD LIBRARIES
# ------------------------------------------------------------------------------
library(shiny)
library(DT)
library(dplyr)
library(ggplot2)
library(grid)

# 2. LOAD DATA
# ------------------------------------------------------------------------------
# Load the main dataset
full_data <- tryCatch({
  readRDS("data/full_data.rds")
}, error = function(e) {
  # Create a placeholder dataframe if the file is missing
  data.frame(
    gene_name = c("pros", "GeneB", "GeneC"), # Ensure this column exists
    chromosome = c("chr1", "chr2", "chrX"),
    log2FC = c(1.5, -0.5, 2.1),
    p_value = c(0.001, 0.2, 0.0005),
    peak_sequence = c("AUCG", "GCUA", "UUUU"),
    notes = c("Note 1", "Note 2", "Note 3"),
    stringsAsFactors = FALSE
  )
})

# Load the character vector of gene names for the dropdown menu.
# This is more efficient for shinylive than scanning a directory.
gene_choices <- tryCatch({
  readRDS("data/target_gene_names.RDS")
}, error = function(e) {
  # Create a placeholder vector if the file is missing
  c("pros", "GeneB", "GeneC")
})


# 3. DEFINE UI (USER INTERFACE)
# ------------------------------------------------------------------------------
ui <- navbarPage(
  title = "Imp Syp iCLIP targets",
  collapsible = TRUE,
  
  # --- Main Page Tab ---
  tabPanel("Data",
           # --- Custom CSS Styling ---
           tags$head(
             tags$style(HTML("
        /* Target the full data table and its components */
        #full_data_table .table td, 
        #full_data_table .table th {
          font-size: 90%; /* Adjust data cell and header font size */
        }

        #full_data_table .form-control {
           font-size: 90%; /* Adjust column filter input font size */
        }
        
        /* Target the wrapper for global search, pagination, etc. */
        .dataTables_wrapper {
          font-size: 90%; /* Adjust font size for table controls */
        }
      "))
           ),
           
           sidebarLayout(
             # --- Sidebar Panel ---
             sidebarPanel(
               width = 3,
               
               # Publication Information Panel
               wellPanel(
                 h4("Publication"),
                 tags$p("This web resource accompanies the publication:"),
                 tags$p(
                   tags$strong("Imp/IGF2BP and Syp/SYNCRIP temporal RNA interactomes uncover combinatorial networks of regulators of Drosophila brain development")
                 ),
                 tags$p("Lee JY, et al. (2025)"),
                 tags$a(
                   href = "https://www.science.org/doi/10.1126/sciadv.adr6682", 
                   target = "_blank",
                   tags$img(src = "science_advances_logo.png", height = "30px", style = "margin-right: 10px;"),
                   tags$br(),
                   "Link to the Article"
                 ),
                 tags$br(),
                 tags$br(),
                 tags$p(
                   tags$a(href = "https://github.com/jefflee1103/Lee2024_Imp-Syp-iCLIP", icon("github"), "GitHub source data", target = "_blank")
                 )
               ),
               
               # Contact Details Panel
               wellPanel(
                 h4("Contact"),
                 tags$p("Jeffrey Y Lee", tags$a(href = "mailto:jeff.lee@glasgow.ac.uk", icon("envelope"), "Email", target = "_blank")),
                 tags$p("Ilan Davis", tags$a(href = "mailto:ilan.davis@glasgow.ac.uk", icon("envelope"), "Email", target = "_blank")),
                 tags$br(),
                 tags$p(
                   tags$a(href = "https://x.com/jefflee1103", icon("x-twitter"), "Twitter", target = "_blank")
                 ),
                 tags$p(
                   tags$a(href = "https://orcid.org/0000-0002-5146-0037", icon("orcid"), "ORCID", target = "_blank")
                 ),
                 tags$p(
                   tags$a(href = "https://bsky.app/profile/jeffylee.bsky.social", icon("bluesky"), "Bluesky", target = "_blank")
                 ),
                 tags$p(
                   tags$a(href = "https://scholar.google.com/citations?user=CxsNxLsAAAAJ&hl=en", icon("google"), "Google Scholar", target = "_blank")
                 )
               )
             ),
             
             # --- Main Panel ---
             mainPanel(
               width = 9,
               tabsetPanel(
                 id = "main_tabs",
                 
                 # Tab 1: Introduction
                 tabPanel("Introduction",
                          div(style = "max-width: 85%;",
                              h3("Imp Syp iCLIP Data Explorer"),
                              p(HTML("This web application provides a user-friendly interface to explore the data from our recent study on the RNA-binding proteins (RBPs) IGFP2BP (Imp) and Syncrip (Syp). Here, you can browse the complete dataset, search for specific genes of interest, and visualise the binding profiles. The goal of this resource is to enhance the reproducibility and accessibility of our findings.")),
                              hr(),
                              h4("Question."),
                              p(HTML("Understanding how the immense complexity of the brain arises from a finite pool of neural stem cells is a central question in developmental neuroscience. A key, evolutionarily conserved strategy for generating this neuronal diversity is temporal patterning, a process where neural stem cells produce different types of neurons in a stereotyped birth order.")),
                              p(HTML("In the fruit fly, <em>Drosophila melanogaster</em>, post-embryonic brain development critically depends on the opposing expression gradients of two conserved RNA-binding proteins (RBPs): <strong>Imp (IGF2BP)</strong> and <strong>Syp (SYNCRIP)</strong> (<em>Fig.1</em>). As development proceeds, Imp levels decrease while Syp levels rise within neural stem cells (<em>Fig.2</em>). This dynamic interplay is known to regulate the generation of diverse neuronal fates, the timing of cell cycle exit, and even the growth of the brain. However, a comprehensive, system-wide understanding of how Imp and Syp achieve this regulation has been missing. The complete set of messenger RNA (mRNA) molecules they bind to directly <em>in vivo</em>, and how these interactions change over time, remained largely unexplored. This study aimed to fill that critical knowledge gap by mapping the temporal RNA interactomes of both Imp and Syp.")),
                              tags$br(),
                              fluidRow(
                                column(6,
                                       tags$figure(
                                         tags$img(src = "intro_image_1.png", width = "100%", style="border: 1px solid #ddd; border-radius: 4px; padding: 5px;"),
                                         tags$figcaption("Fig. 1: Temporal patterning of RBPs regulates neural diversity generation.")
                                       )),
                                column(6,
                                       tags$figure(
                                         tags$img(src = "intro_image_2.png", width = "100%", style="border: 1px solid #ddd; border-radius: 4px; padding: 5px;"),
                                         tags$figcaption("Fig. 2: Imp and Syp expression dynamics across larval brain development.")
                                       )
                                )
                              ),
                              tags$br(),
                              hr(),
                              h4("Experimental approach."),
                              p(HTML("To identify the direct RNA targets of Imp and Syp, we employed <strong>individual-nucleotide resolution UV cross-linking and immunoprecipitation (iCLIP)</strong> (<em>Fig.3</em>). This powerful technique allows for the precise mapping of protein-RNA interaction sites on a transcriptome-wide scale. This was performed at three key developmental time points—late L1, L2, and wandering L3 larval stages—to capture the changing regulatory landscape as the Imp/Syp gradient progresses.")),
                              p(HTML("The resulting dataset is a high-resolution, temporal map of the Imp and Syp RNA interactomes throughout a crucial phase of brain development. It reveals a highly overlapping set of target mRNAs that are dynamically regulated over time, providing a rich resource for dissecting the post-transcriptional networks that drive the specification of neuronal identity.")),
                              tags$br(),
                              fluidRow(
                                column(6,
                                       tags$figure(
                                         tags$img(src = "intro_image_3.png", width = "100%", style="border: 1px solid #ddd; border-radius: 4px; padding: 5px;"),
                                         tags$figcaption("Fig. 3: Schematic of the iCLIP experiment workflow")
                                       )
                                )
                              ),
                              tags$br(),
                              tags$br(),
                              tags$br()
                          )
                 ),
                 
                 # Tab 2: Explore Full Data Table
                 tabPanel("Explore Imp/Syp target table",
                          h3("Complete Imp and Syp iCLIP target dataset"),
                          p("The table below contains the full set of identified Imp and Syp RNA targets. Use the search boxes at the top of each column to filter the data. You can also sort columns by clicking on the headers. For detailed table with further analyses, please see suplementary files in the original publication."),
                          hr(),
                          div(style="background-color: #f5f5f5; padding: 15px; border-radius: 5px; margin-bottom: 20px;",
                              h4("Column Definitions"),
                              p(HTML("<strong>stage_RBP:</strong> Larval developmental stage and RNA-binding protein (RBP) used for iCLIP (e.g. 'L1_Imp' means Imp iCLIP targets at the L1 stage.)")),
                              p(HTML("<strong>gene_id:</strong> Target gene ID.")),
                              p(HTML("<strong>gene_name:</strong> Target gene name.")),
                              p(HTML("<strong>avg_iCLIP_tpm:</strong> iCLIP crosslinks normalised to transcripts per million, average of 3 replicates.")),
                              p(HTML("<strong>max_log2foldchange:</strong> The log2 fold-change of the most enriched binding site in iCLIP compared to size-matched input.")),
                              p(HTML("<strong>n_binding_site:</strong> The number of RBP binding sites on target RNA.")),
                              p(HTML("<strong>binding_feature:</strong> The location of RBP binding sites on target RNA feature.")),
                              p(HTML("<strong>human_homologs:</strong> High confidence human homologs (DIOPT score > 7).")),
                              p(HTML("<strong>mammalian_imp_target:</strong> Conserved target with mammalian IMP1-3 in human pluripotent stem cells.")),
                              p(HTML("<strong>mammalian_syp_target:</strong> Conserved target with mammalian SYNCRIP or HNRNPR in rodent neurons."))
                          ),
                          hr(),
                          DT::dataTableOutput("full_data_table"),
                          downloadButton("download_full_data", "Download Full Table (.csv)", class = "btn-primary", style = "margin-top: 15px;"),
                          tags$br(),
                          tags$br(),
                          tags$br()
                 ),
                 
                 # Tab 3: Explore Binding Sites
                 tabPanel("Explore Imp/Syp binding sites",
                          h3("Visualise Binding Sites for a Specific Gene"),
                          p("Select a gene from the dropdown menu (or type to search) and click 'Search' to visualise binding profile of Imp/Syp."),
                          
                          fluidRow(
                            column(4, 
                                   selectizeInput("gene_search_input", 
                                                  "Select or type a Gene Symbol:", 
                                                  choices = gene_choices,
                                                  selected = "chinmo",
                                                  # options = list(placeholder = 'e.g., pros')
                                                  )
                            ),
                            column(2, 
                                   actionButton("search_button", "Search", icon = icon("search"), style="margin-top: 25px;")
                            )
                          ),
                          
                          hr(),
                          
                          h4(textOutput("plot_title")),
                          plotOutput("gene_plot", height = "800px"),
                          
                          hr(),
                          
                          h4(textOutput("table_title")),
                          DT::dataTableOutput("filtered_gene_table"),
                          downloadButton("download_filtered_data", "Download Filtered Table (.csv)", class = "btn-primary", style = "margin-top: 15px;")
                 )
               )
             )
           )
  )
)


# 4. DEFINE SERVER LOGIC
# ------------------------------------------------------------------------------
server <- function(input, output, session) {
  
  # --- Server logic for Tab 2: Full Data Table ---
  output$full_data_table <- DT::renderDataTable({
    DT::datatable(
      full_data,
      options = list(
        pageLength = 15,
        autoWidth = TRUE,
        scrollX = TRUE
      ),
      filter = 'top',
      rownames = FALSE,
      class = 'cell-border stripe'
    )
  })
  
  # Download handler for the full dataset
  output$download_full_data <- downloadHandler(
    filename = function() {
      paste("ImpSyp_iCLIP_full_dataset_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(full_data, file, row.names = FALSE)
    }
  )
  
  # --- Server logic for Tab 3: Gene-specific Exploration ---
  
  selected_gene <- eventReactive(input$search_button, {
    # Ensure input is not empty before allowing the event to fire
    if (input$gene_search_input == "") {
      return(NULL)
    }
    input$gene_search_input
  }, ignoreNULL = FALSE)
  
  output$plot_title <- renderText({
    req(selected_gene())
    paste("iCLIP binding sites of Imp/Syp on", selected_gene(), "transcript.")
  })
  
  output$table_title <- renderText({
    req(selected_gene())
    paste("Binding Sites in", selected_gene())
  })
  
  # Render the plot on-demand based on the selected gene
  output$gene_plot <- renderPlot({
    gene_name <- selected_gene()
    req(gene_name)
    
    validate(
      need(gene_name %in% gene_choices, "Plot not available. Please select a valid gene and click Search.")
    )
    
    # Construct the file path and load the specific RDS file
    plot_file_path <- file.path("data", "plots", paste0(gene_name, ".RDS"))
    
    validate(
      need(file.exists(plot_file_path), "Plot data file not found. Please contact the administrator.")
    )
    
    grob_to_plot <- readRDS(plot_file_path)
    grid::grid.draw(grob_to_plot)
  })
  
  # Reactive expression to filter data based on the gene selected via the button
  filtered_data_reactive <- reactive({
    gene_name <- selected_gene()
    req(gene_name)
    
    # Correctly filter the data frame using the 'gene_name' column
    filter(full_data, .data$gene_name == gene_name)
  })
  
  output$filtered_gene_table <- renderDataTable({
    df <- filtered_data_reactive()
    
    validate(
      need(nrow(df) > 0, "No data to display. Select a gene and click 'Search'.")
    )
    
    DT::datatable(
      df,
      options = list(pageLength = 10, autoWidth = TRUE),
      rownames = FALSE,
      class = 'cell-border stripe'
    )
  })
  
  # Download handler for the filtered dataset
  output$download_filtered_data <- downloadHandler(
    filename = function() {
      paste("ImpSyp_iCLIP_filtered_", selected_gene(), "_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(filtered_data_reactive(), file, row.names = FALSE)
    }
  )
  
}


# 5. RUN THE APPLICATION
# ------------------------------------------------------------------------------
shinyApp(ui = ui, server = server)
