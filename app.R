#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
#library(checkpoint)
#checkpoint(snapshotDate ='2019-12-17')
#library(AMR)
library(data.table)
library(DT)
library(ggridges)
library(lubridate)
library(plotly)
library(rintrojs)
library(shiny)
library(shinyBS)
library(shinycssloaders)
library(shinydashboard)
library(shinyjs)
library(shinyWidgets)
library(survival)
library(ggplot2)
#library(tidyverse)
library(viridis)
library(readr)
library(readxl)


# Read data
Escherichia_coli_UPEC_536_As <- read.csv("UPEC_As_1As_2As_3Oss_1_UP_Sp_1Sp_2Sp_3_DOWN.csv", sep=",")
Escherichia_coli_EPEC_0127_H6_E2348_69_As <- read.csv("EPEC_As_1As_2As_3_UP.csv", sep=",")
Enterococcus_faecalis_As <- read.csv("ENTFA_As_1As_2As_3_UP.csv", sep=",")

###########Salmonella_As############
Salmonella_As =read.csv("SALMT_As_1As_2As_3_UP.csv",  sep=",")
Salmonella_expr_As =as.data.frame(Salmonella_As[,2:ncol(Salmonella_As)])
rownames(Salmonella_expr_As) = Salmonella_As[,1]

################SALMT_Sp1_up_createdbyCansu###############
Salmonella_enterica_Sp1_up=read.csv("SALMT_Sp_1Sp_2Sp_3_UP.csv", sep=",")
Salmonella_enterica_expr_Sp1_up=as.data.frame(Salmonella_enterica_Sp1_up[,2:ncol(Salmonella_enterica_Sp1_up)])
rownames(Salmonella_enterica_expr_Sp1_up) =Salmonella_enterica_Sp1_up[,1]

################SALMTSp2_up_createdbyCansu###############
Salmonella_enterica_Sp2_up=read.csv("SALMT_Sp_1Sp_2Sp_3Vic_3_UP.csv", sep=",")
Salmonella_enterica_expr_Sp2_up=as.data.frame(Salmonella_enterica_Sp2_up[,2:ncol(Salmonella_enterica_Sp2_up)])
rownames(Salmonella_enterica_expr_Sp2_up) =Salmonella_enterica_Sp2_up[,1]

##################PGFams##########################
PGFams_locustags13df64a3fa583c4875d73be0d98473ba3c3d01cbe802357f34ddfce7e49976a1 <- read_excel("PGFams_locustags13df64a3fa583c4875d73be0d98473ba3c3d01cbe802357f34ddfce7e49976a1.xlsx")

Annotations_all_PGFAM <- PGFams_locustags13df64a3fa583c4875d73be0d98473ba3c3d01cbe802357f34ddfce7e49976a1

#Read again export mathes GO#####

matches_df <- read_csv("matches_df.csv")
subset_GO_df_full_BP <- subset(matches_df, matches_df$`GO.ASPECT` == "P")
annotations = subset_GO_df_full_BP[,c(8,3,4)]
colnames(annotations) <- c("SYMBOL", "GO_TERM", "GO_NAME")

#####Read full dataset########
#########################Salmonella#############################
Salmonella = read.csv("Salmonella_enterica_R_input.csv", sep=";")
dim(Salmonella) #analising input

#getting TPM values
Salmonella_expr =as.data.frame(Salmonella[,2:ncol(Salmonella)])
dim(Salmonella_expr) #analysing dataset

#transforming dataset to start running WGCNA functions
rownames(Salmonella_expr) = Salmonella[,1]
all_genes <- rownames(Salmonella_expr)

# UI -----------------------------------------------------------------------
ui <- dashboardPage(
  skin = "black",
  title = "G0-Pathogenex",
  
  # HEADER ------------------------------------------------------------------
  dashboardHeader(
    title = tags$span(tags$img(src = "Go_pathogenex.png", height = 35), ""),
    titleWidth = 300,
    
    tags$li(
      class = "dropdown",
      a(
        strong("ABOUT GO-PATHOgenex"),
        href = "https://github.com/microbioinformatic/Co-PATHOgenex",
        title = "",
        target = "_blank"
      )
    )
  ),
  
  # SIDEBAR -----------------------------------------------------------------
  dashboardSidebar(
    sidebarMenu(
      menuItem("Tutorial", tabName = "tutorial"),
      menuItem(
        "Selected stimulons",
        menuSubItem("Acid stress", tabName = "acid_stress"),
        menuSubItem("Stationary phase", tabName = "stationary_phase"),
        menuSubItem("Nutritional downwshift", tabName = "nutritional_downwshift")
      ),
      menuItem("General info", tabName = "general_info")
    )
  ),
  
  # BODY --------------------------------------------------------------------
  dashboardBody(
    tabItems(
      tabItem(
        tabName = "tutorial",
        fluidRow(
          box(
            title = "Tutorial",
            "Put your tutorial content here."
          )
        )
      ),
      tabItem(
        tabName = "acid_stress",
        fluidRow(
          box(
            title = "Acid stress stimulons per species",
            radioButtons(
              inputId = "select_acid_stress",
              label = "Select Species",
              choices = c("Salmonella enterica As (+)")
            ),
         #   actionButton("submit_acid_stress", "Submit")
          ),
          box(
            title = textOutput("box_title"),
            plotOutput("Acidstress_boxplot")
          ),
          box(
            title = "Stimulon gene annotation data",
            dataTableOutput("mytableAs")
          )
        )
      ),
      tabItem(
        tabName = "stationary_phase",
        fluidRow(
          box(
            title = "Stationary phase stimulons per species",
            radioButtons(
              inputId = "select_sp",
              label = "Select Species",
              choices = c("Salmonella enterica Sp (+) (I)", "Salmonella enterica Sp (+) (II)")
            ),
          #  actionButton("submit_sp_stress", "Submit")
          ),
          box(
            title = textOutput("box_title_sp"),
            plotOutput("SPstress_boxplot")
          ),
          box(
            title = "Stimulon gene annotation data",
            dataTableOutput("mytableSp")
          )
        )
      ),
      tabItem(
        tabName = "nutritional_downwshift",
        fluidRow(
          # Similar layout structure as "acid_stress" tab
          # Adjust as needed
        )
      ),
      tabItem(
        tabName = "general_info",
        fluidRow(
          box(
            title = "General info",
            "Put your general info content here."
          )
        )
      )
    )
  )
)

    
  
  
    
 


# Server  ---------------------------------------------------------------
server <- function(input, output, session) {
  
  ###Define reactive data based on user input#####
  ###########AS_stimulon############# 
  datasetInputStimulonAcid <- reactive({
    switch(input$select_acid_stress,
           "Salmonella enterica As (+)"= (Salmonella_expr_As)
           # Add more cases for other species...
    )
  })
  
  ###########SP_stimulon#############
  datasetInputStimulonSp <- reactive({
    switch(input$select_sp,
           "Salmonella enterica Sp (+) (I)" = (Salmonella_enterica_expr_Sp1_up),
           "Salmonella enterica Sp (+) (II)" = (Salmonella_enterica_expr_Sp2_up)
           
           
    )
  })
 
  # Initialize reactiveValues to store selected options
  selectedOptions <- reactiveValues(acid_stress = NULL, stationary_phase = NULL)
  
  #####Function boxplot########
  generateBoxPlot <- function(datasetInput) {
    return(reactive({
      AcidStress <- datasetInput()
      if (nrow(AcidStress) == 36 ) {
        generic_stress_sample_names <- c("As_1","As_2","As_3","Bs_1","Bs_2","Bs_3","Ctrl_1","Ctrl_2","Ctrl_3","Li_1","Li_2","Li_3","Mig_1","Mig_2","Mig_3","Nd_1","Nd_2","Nd_3","Ns_1","Ns_2","Ns_3","Oss_1","Oss_2","Oss_3","Oxs_1","Oxs_2","Oxs_3","Sp_1","Sp_2","Sp_3","Tm_1","Tm_2","Tm_3","Vic_1","Vic_2","Vic_3")
        ME_lengthOfVector <- rep(NA, ncol(AcidStress) * 36)
        Module_names <- rep(colnames(AcidStress), 36)
        lengthOfVector_col <- ncol(AcidStress)
        position_ME <- 1
        for(i in 1:lengthOfVector_col) {
          for(j in 1:36) {
            ME_lengthOfVector[position_ME] <- AcidStress[j,i]
            Module_names[position_ME] <- colnames(AcidStress)[i]
            position_ME <- position_ME + 1
          }
          position_ME <- i * 36 + 1
        }
        stress <- rep(generic_stress_sample_names, ncol(AcidStress))
        stress_1 <- gsub("_1","", stress)
        stress_2 <- gsub("_2","", stress_1)
        stress_3 <- gsub("_3","", stress_2)
        stress <- stress_3
        
        ME_Kruskall_dataset <- as.data.frame(cbind(ME_lengthOfVector, stress, Module_names))
        colnames(ME_Kruskall_dataset) <- c("MEigengene", "Stress", "Module_names")
        Plot <- ggplot(ME_Kruskall_dataset, aes(x=Stress, y=as.numeric(MEigengene), fill = Stress)) + labs(y="Standard deviation of gene expression", x = "Stress conditions")
        Plot + geom_violin(width=1, color="white") + geom_boxplot(width = .1, color = "black", alpha = .5) + scale_fill_brewer(palette = "Paired", aesthetics = "fill", guide = "none") + theme_classic(base_size = 16)
      } else if (nrow(AcidStress) == 33) {
        generic_stress_sample_names <- c("As_1","As_2","As_3","Ctrl_1","Ctrl_2","Ctrl_3","Li_1","Li_2","Li_3","Mig_1","Mig_2","Mig_3","Nd_1","Nd_2","Nd_3","Ns_1","Ns_2","Ns_3","Oss_1","Oss_2","Oss_3","Oxs_1","Oxs_2","Oxs_3","Sp_1","Sp_2","Sp_3","Tm_1","Tm_2","Tm_3","Vic_1","Vic_2","Vic_3")
        ME_lengthOfVector <- rep(NA, ncol(AcidStress) * 33)
        Module_names <- rep(colnames(AcidStress), 33)
        lengthOfVector_col <- ncol(AcidStress)
        position_ME <- 1
        for(i in 1:lengthOfVector_col) {
          for(j in 1:33) {
            ME_lengthOfVector[position_ME] <- AcidStress[j,i]
            Module_names[position_ME] <- colnames(AcidStress)[i]
            position_ME <- position_ME + 1
          }
          position_ME <- i * 33 + 1
        }
        stress <- rep(generic_stress_sample_names, ncol(AcidStress))
        stress_1 <- gsub("_1","", stress)
        stress_2 <- gsub("_2","", stress_1)
        stress_3 <- gsub("_3","", stress_2)
        stress <- stress_3
        
        ME_Kruskall_dataset <- as.data.frame(cbind(ME_lengthOfVector, stress, Module_names))
        colnames(ME_Kruskall_dataset) <- c("MEigengene", "Stress", "Module_names")
        Plot <- ggplot(ME_Kruskall_dataset, aes(x=Stress, y=as.numeric(MEigengene), fill = Stress)) + labs(y="Standard deviation of gene expression", x = "Stress conditions")
        Plot + geom_violin(width=1, color="white") + geom_boxplot(width = .1, color = "black", alpha = .5) + scale_fill_brewer(palette = "Paired", aesthetics = "fill", guide = "none") + theme_classic(base_size = 16)
      } 
    }))
  }
  ###end boxplot function#####
  # Observe changes in select_acid_stress input and update selected option
  observeEvent(input$select_acid_stress, {
    selectedOptions$acid_stress <- switch(input$select_acid_stress,
                                          "Salmonella enterica As (+)" = (Salmonella_expr_As)
                                          # Add more cases for other species...
    )
  })
  
  # Observe changes in select_sp input and update selected option
  observeEvent(input$select_sp, {
    selectedOptions$stationary_phase <- switch(input$select_sp,
                                               "Salmonella enterica Sp (+) (I)" = (Salmonella_enterica_expr_Sp1_up),
                                               "Salmonella enterica Sp (+) (II)" = (Salmonella_enterica_expr_Sp2_up)
    )
  })
  
  
  ######Define function for AS######## 
  # Define reactive expression to dynamically generate the box title
  box_title <- reactive({
    paste("Stimulon expression for:", input$select_acid_stress)
  })
  
  # box_plot_AS <- generateBoxPlot(datasetInputStimulonAcid) 
  
  # Generate box plots using the function
  output$Acidstress_boxplot <- renderPlot({
    req(selectedOptions$acid_stress)
    box_plot_Acidstress <- generateBoxPlot(datasetInputStimulonAcid)
    box_plot_Acidstress()
  })
  
  
  # Update box title dynamically
  output$box_title <- renderText({
    box_title()
  })
  
  
  
  
  #####data table AS#########
  
  datAs <- reactive({
    dataset_As <- datasetInputStimulonAcid()
    Locus_Tag_As= colnames(dataset_As)
    dataset_transpose_As=t(dataset_As)
    rownames(dataset_transpose_As) <- NULL
    colnames(dataset_transpose_As) <- NULL
    
    probes = Locus_Tag_As
    Genes_As=Annotations_all_PGFAM$ref_genome_gene_name[match(probes,Annotations_all_PGFAM$comp_genome_1_patric_id)]
    PGFam_As=Annotations_all_PGFAM$ref_genome_function[match(probes,Annotations_all_PGFAM$comp_genome_1_patric_id)]
    table_2_inspect_As= as.data.frame(cbind(Locus_Tag_As,Genes_As,PGFam_As,dataset_transpose_As))
    
    if (nrow(dataset_As) == 36 ) {
      colnames(table_2_inspect_As) <- c("Locus Tag","Genes", "PGFam description","As_1","As_2","As_3","Bs_1","Bs_2","Bs_3","Ctrl_1","Ctrl_2","Ctrl_3","Li_1","Li_2","Li_3","Mig_1","Mig_2","Mig_3","Nd_1","Nd_2","Nd_3","Ns_1","Ns_2","Ns_3","Oss_1","Oss_2","Oss_3","Oxs_1","Oxs_2","Oxs_3","Sp_1","Sp_2","Sp_3","Tm_1","Tm_2","Tm_3","Vic_1","Vic_2","Vic_3")
      return(as.data.frame(table_2_inspect_As))
    } else if (nrow(dataset_As) == 33 ) {
      colnames(table_2_inspect_As) <- c("Locus Tag","Genes", "PGFam description","As_1","As_2","As_3","Ctrl_1","Ctrl_2","Ctrl_3","Li_1","Li_2","Li_3","Mig_1","Mig_2","Mig_3","Nd_1","Nd_2","Nd_3","Ns_1","Ns_2","Ns_3","Oss_1","Oss_2","Oss_3","Oxs_1","Oxs_2","Oxs_3","Sp_1","Sp_2","Sp_3","Tm_1","Tm_2","Tm_3","Vic_1","Vic_2","Vic_3")
      return(as.data.frame(table_2_inspect_As))
    } 
  })
  
  output$table_box <- renderUI({
    box(
      title = "Your Table Title",
      status = "info",
      solidHeader = TRUE,
      collapsible = TRUE,
      width = 10,
      DTOutput("mytableAs")
    )
  })
  
  output$mytableAs <- DT::renderDataTable({
    # Your table rendering logic here
    dat <- datAs()  # Assuming you have a function datAs() to get your data
    datatable(dat, extensions = c('Buttons', 'FixedColumns'),  # Include FixedColumns extension
              options = list(
                paging = TRUE, 
                lengthMenu = c(50, 100, 150), 
                scrollX = TRUE,  # Enable horizontal scrolling
                lengthChange = TRUE, 
                searching = TRUE,
                initComplete = JS(js),
                dom = 'Bfrtip',
                buttons = c('copy', 'csv', 'excel'),
                fixedColumns = list(leftColumns = 3)  # Freeze the first column
              ),
              rownames = FALSE, 
              class = "display"
    ) %>% formatStyle(names(dat[-1:-3]), backgroundColor = styleInterval(brks, clrs))
  })
  
  

  
  
  
  

  
  
  
  
  ####Define function for SP#######
  
  # Define reactive expression to dynamically generate the box title
  box_title_sp <- reactive({
    paste("Stimulon expression for:", input$select_sp)
  })
  
  # box_plot_AS <- generateBoxPlot(datasetInputStimulonAcid) 
  
  # Generate box plots using the function
  output$SPstress_boxplot <- renderPlot({
    req(selectedOptions$stationary_phase)
    box_plot_SPstress <- generateBoxPlot(datasetInputStimulonSp)
    box_plot_SPstress()
  })
  
  # Update box title dynamically for SP
  output$box_title_sp <- renderText({
    box_title_sp()
  })
  
  
  
  #####data table SP#########
  
  datSp <- reactive({
    dataset_As <- datasetInputStimulonSp()
    Locus_Tag_As= colnames(dataset_As)
    dataset_transpose_As=t(dataset_As)
    rownames(dataset_transpose_As) <- NULL
    colnames(dataset_transpose_As) <- NULL
    
    probes = Locus_Tag_As
    Genes_As=Annotations_all_PGFAM$ref_genome_gene_name[match(probes,Annotations_all_PGFAM$comp_genome_1_patric_id)]
    PGFam_As=Annotations_all_PGFAM$ref_genome_function[match(probes,Annotations_all_PGFAM$comp_genome_1_patric_id)]
    table_2_inspect_As= as.data.frame(cbind(Locus_Tag_As,Genes_As,PGFam_As,dataset_transpose_As))
    
    if (nrow(dataset_As) == 36 ) {
      colnames(table_2_inspect_As) <- c("Locus Tag","Genes", "PGFam description","As_1","As_2","As_3","Bs_1","Bs_2","Bs_3","Ctrl_1","Ctrl_2","Ctrl_3","Li_1","Li_2","Li_3","Mig_1","Mig_2","Mig_3","Nd_1","Nd_2","Nd_3","Ns_1","Ns_2","Ns_3","Oss_1","Oss_2","Oss_3","Oxs_1","Oxs_2","Oxs_3","Sp_1","Sp_2","Sp_3","Tm_1","Tm_2","Tm_3","Vic_1","Vic_2","Vic_3")
      return(as.data.frame(table_2_inspect_As))
    } else if (nrow(dataset_As) == 33 ) {
      colnames(table_2_inspect_As) <- c("Locus Tag","Genes", "PGFam description","As_1","As_2","As_3","Ctrl_1","Ctrl_2","Ctrl_3","Li_1","Li_2","Li_3","Mig_1","Mig_2","Mig_3","Nd_1","Nd_2","Nd_3","Ns_1","Ns_2","Ns_3","Oss_1","Oss_2","Oss_3","Oxs_1","Oxs_2","Oxs_3","Sp_1","Sp_2","Sp_3","Tm_1","Tm_2","Tm_3","Vic_1","Vic_2","Vic_3")
      return(as.data.frame(table_2_inspect_As))
    } 
  })
  
  
  
  output$mytableSp <- DT::renderDataTable({
    # Your table rendering logic here
    dat <- datSp()  # Assuming you have a function datAs() to get your data
    datatable(dat, extensions = c('Buttons', 'FixedColumns'),  # Include FixedColumns extension
              options = list(
                paging = TRUE, 
                lengthMenu = c(50, 100, 150), 
                scrollX = TRUE,  # Enable horizontal scrolling
                lengthChange = TRUE, 
                searching = TRUE,
                initComplete = JS(js),
                dom = 'Bfrtip',
                buttons = c('copy', 'csv', 'excel'),
                fixedColumns = list(leftColumns = 3)  # Freeze the first column
              ),
              rownames = FALSE, 
              class = "display"
    ) %>% formatStyle(names(dat[-1:-3]), backgroundColor = styleInterval(brks, clrs))
  })
  
  
  
  js <- c(
    "function(settings){",
    "  var instance = settings.oInstance;",
    "  var table = instance.api();",
    "  var input = instance.parent().find('.dataTables_filter input');",
    "  input.off('keyup search input').on('keyup', function(){",
    "    var keyword = '\\\\b' + input.val() + '\\\\b';",
    "    table.search(keyword, true, false).draw();",
    "  });",
    "}"
  )
  color_from_middle <- function (data,color1,color2) 
  {
    max_val=max(abs(data))
    JS(sprintf("isNaN(parseFloat(value)) || value < 0 ? 'linear-gradient(90deg, transparent, transparent ' + (50 + value/%s * 50) + '%%, %s ' + (50 + value/%s * 50) + '%%,%s  50%%,transparent 50%%)': 'linear-gradient(90deg, transparent, transparent 50%%, %s 50%%, %s ' + (50 + value/%s * 50) + '%%, transparent ' + (50 + value/%s * 50) + '%%)'",
               max_val,color1,max_val,color1,color2,color2,max_val,max_val))
  } 
  
  brks = c(-2,-1.5, -1.25, -1, -0.5,-0.25, 0, 0.25, 0.5, 1, 1.25, 1.5, 2)
  clrs = colorRampPalette(c("blue", "white", "red"))(length(brks) + 0.25)
  
 

  
  
}

# Run the application 
shinyApp(ui = ui, server = server)