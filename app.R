
### Load in required libraries
library(BiocManager)
options(repos = BiocManager::repositories())
options(rsconnect.http = "curl")
library(aws.s3)
library(shiny)
library(shinyWidgets)
library(ComplexHeatmap)
library(ggplot2)
library(ggExtra)
library(circlize)
library(Matrix)
library(shinycssloaders)
library(shinyhelper)
library(magrittr)
library(dplyr)
library(shinyjs)
library(viridis)
library(egg)
library(readxl)
library(cowplot)
library(patchwork)
library(RColorBrewer)
library(Seurat)

options(Seurat.object.assay.version = 'v5')

get_bucket("s3:////hmsrc-cbdm-data.s3.amazonaws.com", region = "us-east-2", check_region = F,verbose = TRUE)


### Comment this line out when deploying
# setwd(dirname(rstudioapi::getSourceEditorContext()$path))

integrated_data<-read.csv('./integrated_umap.csv')
integrated_data$cell_name<-integrated_data$X
integrated_data$cell_name<-sub("_[^_]+$", "", integrated_data$cell_name)
integrated_data<-integrated_data[!duplicated(integrated_data$cell_name),] 
integrated_data$Integrated_Clusters<-as.character(integrated_data$cluster)
rownames(integrated_data)<-integrated_data$cell_name

#integrated_data<-integrated_data[-c(1,2,3,4)]

### Read in the ADT Information table
ADT_table <-
  as.data.frame(read.csv("./adt_hash_seq_file_2022_panel.csv"))

### Master file for experiments
sample_metadata <- read.csv("./sample_metadata.csv")

### Displaying webpage
ui <- fluidPage(
  tags$head(includeHTML("google-analytics.html")),
  shinyjs::useShinyjs(),
  
  
  ### Header
  titlePanel(fluidRow(column(
    4, img(src = "ImmGenTIcon.png", height="60", width="50%")
  ),
  column(
    4,
    h1(
      id = "heading",
      "Rosetta",
      align = "center",
      style = 'padding:0px;'
    )
  )),
  windowTitle = "Protein/RNA Mosaic Viewer"),
  tags$head(tags$style(HTML(
    "#heading{color: darkred}"
  )),
  tags$title("Protein/RNA Mosaic Viewer")),
  
  sidebarLayout(
    sidebarPanel(
      id = "sidebar",
      actionButton(
        'help',
        img(src = "What.png",
            style = "padding-bottom:0px;border-radius: 0px;border-width: 0px")
      ),
      HTML("<br> <br>"),
      selectInput(
        "selected_igt_id",
        "Dataset",
        choices = unique(sample_metadata$igt_description),
        selected = '2_Spleen BM Thymus Blood at steadystate'
      ),

      conditionalPanel(condition = "input.tabs == 'Heatmap'",
                       selectizeInput("cluster", "Cluster color", as.list(c("RNA_clusters","Protein_clusters", 'sample_name','rna/protein expression','list of cells')),
                                      selected = "RNA_clusters")),
      conditionalPanel(condition = "input.tabs == 'Heatmap_GC'",
                       selectizeInput("cluster2", "Cluster color", as.list(c("RNA_clusters","Protein_clusters", 'sample_name')),
                                      selected = "sample_name")),
      conditionalPanel(condition = "input.tabs == 'Heatmap_GC' ",
                       selectizeInput("datatype", "Data Type", as.list(c("RNA Plot","Protein Plot")),
                                      selected = "RNA Plot")),
      conditionalPanel(condition = "input.tabs == 'Heatmap' | input.tabs == 'facs'",
                       selectizeInput("selected_igt_samples", "Samples",        
                                      choices = NULL,
                                      multiple = TRUE)),
      
      conditionalPanel(condition = "input.tabs == 'Heatmap_GC' ",
                       selectizeInput("datatype", "Data Type", as.list(c("RNA Plot","Protein Plot")),
                                      selected = "RNA Plot")),
      
      conditionalPanel(condition = "input.tabs == 'Integrated_Data' ",
                       selectizeInput("datatype2", "IGT Data Type", as.list(c("RNA Plot","Protein Plot")),
                                      selected = "RNA Plot")),
      conditionalPanel(condition = "input.tabs == 'Integrated_Data'" ,
                       selectizeInput("cluster_int", "Cluster color", choices=NULL
                       )),

      
      
      conditionalPanel(condition = "input.tabs == 'Heatmap_GC'",
                       selectizeInput("selected_samples_a", "Samples_a",        
                                      choices = NULL,
                                      multiple = TRUE)),
      conditionalPanel(condition = "input.tabs == 'Heatmap_GC'",
                       selectizeInput("selected_samples_b", "Samples_b",        
                                      choices = NULL,
                                      multiple = TRUE)),
      
      conditionalPanel(condition = "input.tabs == 'Heatmap' | input.tabs == 'facs'",
                       downloadLink('downloadCells', '
                   Download Selected Cells')
                       
      ),
      conditionalPanel(condition = "input.tabs == 'facs'",
                       downloadLink('downloadData', 'Ab list is here')
                       
      ),
      conditionalPanel('output.fileUploaded & input.tabs == "Heatmap_GC"',
                       downloadLink('download_DE_genes', '
                   Download full list DE Genes')
                       
      ),
      conditionalPanel(  'output.fileUploaded & input.tabs == "Heatmap_GC"',
                       downloadLink('download_DE_proteins', '
                   Download full list DE Proteins')
                       
      ),
      conditionalPanel(
        condition = "input.cluster == 'rna/protein expression' && input.tabs == 'Heatmap'",
        selectizeInput("gene", "Select Gene",choices = NULL,options = NULL)
      ),
      conditionalPanel(condition = "input.tabs == 'Heatmap'",
                       uiOutput("ProteinSelector")),
      conditionalPanel(
        condition = "input.cluster == 'list of cells' && input.tabs == 'Heatmap'",
        fileInput("file1", "Choose CSV File",
                  accept = c(
                    "text/csv",
                    "text/comma-separated-values,text/plain",
                    ".csv")
        )),
      textOutput("description"),
    
      width = 2
    ),
    mainPanel(tabsetPanel(
      id = 'tabs',
      tabPanel(
        'Heatmap',
        value = 'Heatmap',
        splitLayout(withSpinner(
          plotOutput(
            "RNA_plot",
            height = 400,
            width = 500,
            brush = brushOpts(
              id = "plot4_brush",
              resetOnNew =
                FALSE,
              delay = 10000,
              delayType =
                "debounce"
            )
          )
        ),
        withSpinner(
          plotOutput(
            "ADT_plot",
            height = 400,
            width = 500,
            brush = brushOpts(
              id = "ADT_brush",
              resetOnNew =
                FALSE,
              delay = 10000,
              delayType =
                "debounce"
            )
          )
        ),
        align = "center"),
        withSpinner(plotOutput("heatmap"))
      ),
      tabPanel(
        'Flow Cytometry Style Plots',
        value = 'facs',
        fluidRow(
          align = "center",
          column(
            id = "first",
            3,
            style = 'padding:10px;',
            withSpinner(
              plotOutput(
                "prot_plot1",
                height = 200,
                width = 200,
                brush = brushOpts(id =
                                    "facs1_brush",
                                  resetOnNew =
                                    TRUE)
              )
            ),
            fixedRow(
              column(
                6,
                align = "center",
                selectizeInput("plot1_X",
                               "X axis:",
                               as.list(c("")),
                               selected =
                                 "TCRB")
              ),
              column(
                6,
                align = "center",
                selectizeInput("plot1_Y", "Y axis:", as.list(c("")),
                               selected =
                                 "CD19")
              ),
              column(6, align =
                       "center")
              
              
            ),
            align = "center"
          ),
          column(
            id = "second",
            3,
            fixedRow(
              column(
                1,
                radioButtons(
                  "highlight_filter",
                  "",
                  choices = c("Highlight", "Filter"),
                  width = "10px"
                )
              ),
              column(
                style = 'padding:10px;',
                3,
                offset = 2,
                withSpinner(
                  plotOutput(
                    "prot_plot2",
                    height = 200,
                    width = 200,
                    brush = brushOpts(id =
                                        "facs2_brush", resetOnNew = FALSE)
                  )
                )
              )
            ),
            fixedRow(
              column(
                6,
                align = "center",
                selectizeInput("plot2_X", "X axis:", as.list(c("")),
                               selected =
                                 "CD4")
              ),
              column(
                6,
                align = "center",
                selectizeInput("plot2_Y", "Y axis:", as.list(c("")),
                               selected =
                                 "CD8a")
              ),
              column(6, align =
                       "center", selectizeInput("plot2_gate", "Gate 2:", c("A", "none")))
            )
          ),
          column(
            id = "third",
            3,
            fixedRow(
              column(
                1,
                radioButtons(
                  "highlight_filter2",
                  "",
                  choices = c("Highlight", "Filter"),
                  width = "10px"
                )
              ),
              column(
                style = 'padding:10px;',
                3,
                offset = 2,
                withSpinner(
                  plotOutput(
                    "prot_plot3",
                    height = 200,
                    width = 200,
                    brush = brushOpts(id =
                                        "facs3_brush", resetOnNew = FALSE)
                  )
                )
              )
            ),
            fixedRow(
              column(
                6,
                align = "center",
                selectizeInput("plot3_X", "X axis:", as.list(c("")),
                               selected =
                                 "CD44")
              ),
              column(
                6,
                align = "center",
                selectizeInput("plot3_Y", "Y axis:", as.list(c("")),
                               selected =
                                 "CD62L")
              ),
              column(
                6,
                align = "center",
                selectizeInput("plot3_gate", "Gate 3:", c("A", "B", "A+B", "none"), selected =
                                 "A+B")
              )
            ),
            align = "right"
          ),
          
          column(
            id = "fourth",
            3,
            fixedRow(
              column(
                1,
                radioButtons(
                  "highlight_filter3",
                  "",
                  choices = c("Highlight", "Filter"),
                  width = "10px"
                )
              ),
              column(
                style = 'padding:10px;',
                3,
                offset = 2,
                withSpinner(plotOutput(
                  "prot_plot4", height = 200, width = 200
                ))
              )
            ),
            fixedRow(
              column(
                6,
                align = "center",
                selectizeInput("plot4_X", "X axis:", as.list(c("")),
                               selected =
                                 "CD44")
              ),
              column(
                6,
                align = "center",
                selectizeInput("plot4_Y", "Y axis:", as.list(c("")),
                               selected =
                                 "CD62L")
              ),
              column(
                6,
                align = "center",
                selectizeInput("plot4_gate", "Gate 4:", c("A", "B", "C","A+B+C", "none"), selected =
                                 "A+B+C")
              )
            ),
            align = "right"
          )
          
          
          
          
          
        ),
        column(
          style = 'padding-top:70px;padding-bottom:0px',
          width = 12,
          selectizeInput(
            "UMAP_gate",
            "UMAP gate:",
            c("A", "B", "C","A+B+C", "none"),
            selected = "A+B+C",
            width = "200px"
          ),
          align = "center"
        ),
        tags$style(
          type = 'text/css',
          ".selectize-input { font-size: 12px} .selectize-dropdown { text-align:center !important; font-size: 12px } .container-fluid {  min-width: 1800px; }",
          HTML(
            "
      #first {
          border: 1px solid black;height=350px;width=350px;
      }
      #second {
          border: 1px solid black;height=350px;width=350px;
      }
      #third {
          border: 1px solid black;height=350px;width=350px;
      }
    "
          )
        ),
        withSpinner(plotOutput("rnaumap"))
      ),
      tabPanel(
        'Heatmap Group Comparisons',
        value = 'Heatmap_GC',
        fluidRow(
          align = "center",
          column(
            id = "Plot1",
            6,
            style = 'padding:10px;',
            withSpinner(
              plotOutput(
                "RNA_plot_GC_1",
                height = 400,
                width = 500,
                brush = brushOpts(
                  id = "plot5_brush",
                  resetOnNew =
                    FALSE,
                  delay = 10000,
                  delayType =
                    "debounce"
                )
              )
            )),
          column(
            id = "Plot2",
            6,
            style = 'padding:10px;',
            withSpinner(
              plotOutput(
                "RNA_plot_GC_2",
                height = 400,
                width = 500,
                brush = brushOpts(
                  id = "plot6_brush",
                  resetOnNew =
                    FALSE,
                  delay = 10000,
                  delayType =
                    "debounce"
                )
              )
            )),
          column(
            id = "HM1",
            6,
            style = 'padding:10px;',
            withSpinner((plotOutput("heatmap_GC1")))),
          column(
            id = "HM2",
            6,
            style = 'padding:10px;',
            withSpinner(plotOutput("heatmap_GC2")))
        )
        
      ),
      tabPanel(
        'Integrated_Data',
        value = 'Integrated_Data',
        splitLayout(
          withSpinner(
            plotOutput(
              "ADT_plot_3",
              height = 380,
              width = 450,
              brush = brushOpts(
                id = "ADT_brush_3",
                resetOnNew =
                  FALSE,
                delay = 10000,
                delayType =
                  "debounce"
              )
            )
          ),withSpinner(
            plotOutput(
              "Integrated_plot",
              height = 500,
              width = 580
            )
          ),
          align = "center")    )
    )
    )
  )
)




server <- function(input, output, session) {
  observeEvent(input$help, {
    showModal(
      modalDialog(
        align = "center",
        title = "Help",
        HTML(
          "This is an interactive tool to help compare single cell Protein and RNA landscapes.
      <br> <br>
      On the Heatmap tab (default), select cells on either the Protein or RNA UMAP.
           The web tool will output the most differential proteins and most differential genes for your selected
           cell population as two heatmaps. It will also highlight the cells you selected on the other unselected UMAP.
           The outputed heatmaps show the same cell orders between the two heatmaps. <br> <br> Note: It may take
           10-30 seconds for the sample list to appear for subsetting when a dataset is selected.
           <br> <br> To clear,
           click anywhere on the plot you highlighted. <br> <br>"
        ),
        img(
          src = "HeatmapStep1.png",
          height = "50%",
          width = "50%"
        ),
        HTML(
          "<br> <br>
      Select the flow cytometry style plots tab to plot protein markers on XY axes. Each subsequent plot will output
      the selected cells according to the selected proteins for the XY axes. Gates can be selected and
      RNA and ADT umaps will be outputed after each plot depending on the highlighting or filtering used. <br> <br>
      On the flow cytometry style tab, select 2 protein markers for the first XY axis and select the population
      you are interested in.<br> <br>"
        ),
        img(
          src = "FlowStep1.png",
          height = "50%",
          width = "50%"
        ),
        HTML(
          "<br> <br> The selected populations will show on the plot plotted as the X axis and Y axis. <br> <br>"
        ),
        img(
          src = "Gatechange.png",
          height = "25%",
          width = "25%"
        ),
        HTML(
          "<br> <br> Next, select the cells of interest in the second plot in order to plot in the third and final plot.
      This final plot will show the cells plotted on the XY axis selected."
        ),
        img(
          src = "Flowrevisited.png",
          height = "50%",
          width = "50%"
        ),
        HTML(
          "<br> <br> To change what cells are shown in the umap (either Gate A, Gate B, Gate A and B,
      or all cells in dataset), select the appropriate option under UMAP gate."
        )
      )
    )
  })
  
  observeEvent(selectedIGT(), {
    
    if (input$selected_igt_id == 'Integrated_Data'){
      sc2 <- readRDS('Integrated/totalvi_igt1_56_20231030_allgenes_seurat_sketch.Rds')
      sc <- sc2
      Idents(sc) <- 'RNA_clusters'
      choices<- levels(sc)
      meta_int<- c('none',"RNA_clusters")
      updateSelectizeInput(
        session,
        inputId = "cluster_int",
        label = "IGT Data: Cluster color",
        choices = as.list(meta_int),
        selected = "none"
      )
      meta_int_2<- c('Integrated_Clusters')
      
      
      
      updateSelectInput(inputId = "selected_igt_samples", choices = choices, selected = NULL)
      updateSelectInput(inputId = "selected_samples_a", choices = choices, selected = NULL)
      updateSelectInput(inputId = "selected_samples_b", choices = choices, selected = NULL)
      
    
    }
    else{
    
    sc2 <- aws.s3::s3read_using(readRDS,object = paste0("s3:////hmsrc-cbdm-data/hsmrc-cbdm-data/IGT", strsplit(input$selected_igt_id, "_")[[1]][1], "/dataset_clean.Rds"))
    
     sc <- sc2
    Idents(sc) <- 'sample_name'
    choices<- levels(sc)
    
    
    
    meta_int<- c('none',"RNA_clusters", 'sample_name')
    updateSelectizeInput(
      session,
      inputId = "cluster_int",
      label = "IGT Data: Cluster color",
      choices = as.list(meta_int),
      selected = "none"
    )
    
    meta_int_2<- c('Integrated_Clusters')

    
    
    updateSelectInput(inputId = "selected_igt_samples", choices = choices, selected = NULL)
    updateSelectInput(inputId = "selected_samples_a", choices = choices, selected = NULL)
    updateSelectInput(inputId = "selected_samples_b", choices = choices, selected = NULL)
    
    }
    
  })
  
  #according to the dataset selected, update the adt protein list
  datasetInput <- reactive({
    if (input$selected_igt_id == 'Integrated_Data'){
      sc2 <- readRDS('Integrated/totalvi_igt1_56_20231030_allgenes_seurat_sketch.Rds')
      sc <- sc2
      
      

      sc_a<-sc

      adtumap1_a <- sc_a@reductions$umap_adt@cell.embeddings[, 1]
      adtumap2_a <- sc_a@reductions$umap_adt@cell.embeddings[, 2]
      rnaumap1_a <- sc_a@reductions$umap_rna@cell.embeddings[, 1]
      rnaumap2_a <- sc_a@reductions$umap_rna@cell.embeddings[, 2]
      
      data_a <- as.data.frame(t(as.matrix(sc_a@assays$ADT@counts)))
      data_a <- log2(data_a + 1)
      data_a <- data_a[, order(colnames(data_a))]
      
      rna_data_a <- sc_a@assays$RNA@data
      
      
      
      sc_b<-sc

      
      adtumap1_b <- sc_b@reductions$umap_adt@cell.embeddings[, 1]
      adtumap2_b <- sc_b@reductions$umap_adt@cell.embeddings[, 2]
      rnaumap1_b <- sc_b@reductions$umap_rna@cell.embeddings[, 1]
      rnaumap2_b <- sc_b@reductions$umap_rna@cell.embeddings[, 2]
      
      data_b <- as.data.frame(t(as.matrix(sc_b@assays$ADT@counts)))
      data_b <- log2(data_b + 1)
      data_b <- data_b[, order(colnames(data_b))]
      
      rna_data_b <- sc_b@assays$RNA@data
      

      adtumap1 <- sc@reductions$umap_adt@cell.embeddings[, 1]
      adtumap2 <- sc@reductions$umap_adt@cell.embeddings[, 2]
      rnaumap1 <- sc@reductions$umap_rna@cell.embeddings[, 1]
      rnaumap2 <- sc@reductions$umap_rna@cell.embeddings[, 2]
      
      data <- as.data.frame(t(as.matrix(sc@assays$ADT@counts)))
      data <- log2(data + 1)
      data <- data[, order(colnames(data))]
      
      rna_data <- sc@assays$RNA@data
      proteinList <- rownames(sc@assays$ADT)
      geneList <- sc@assays$RNA@data@Dimnames[1]
      
    }
    else{
    
    
    sc2 <- aws.s3::s3read_using(readRDS,object = paste0("s3:////hmsrc-cbdm-data/hsmrc-cbdm-data/IGT", strsplit(input$selected_igt_id, "_")[[1]][1], "/dataset_clean.Rds"))
    
    
    
    sc <- sc2
    
    
    samples_to_subset <- input$selected_igt_samples
 
    sc_a<-sc
    samples_to_subset_a <- input$selected_samples_a
    
    Idents(sc_a) <- sc_a$sample_name
    sc_a <- subset(sc_a, idents = samples_to_subset_a)
    
    adtumap1_a <- sc_a@reductions$umap_adt@cell.embeddings[, 1]
    adtumap2_a <- sc_a@reductions$umap_adt@cell.embeddings[, 2]
    rnaumap1_a <- sc_a@reductions$umap_rna@cell.embeddings[, 1]
    rnaumap2_a <- sc_a@reductions$umap_rna@cell.embeddings[, 2]
    
    data_a <- as.data.frame(t(as.matrix(sc_a@assays$ADT@counts)))
    data_a <- log2(data_a + 1)
    data_a <- data_a[, order(colnames(data_a))]
    
    rna_data_a <- sc_a@assays$RNA@data
    
    
    
    sc_b<-sc
    samples_to_subset_b <- input$selected_samples_b
    
    Idents(sc_b) <- sc_b$sample_name
    sc_b <- subset(sc_b, idents = samples_to_subset_b)
    
    adtumap1_b <- sc_b@reductions$umap_adt@cell.embeddings[, 1]
    adtumap2_b <- sc_b@reductions$umap_adt@cell.embeddings[, 2]
    rnaumap1_b <- sc_b@reductions$umap_rna@cell.embeddings[, 1]
    rnaumap2_b <- sc_b@reductions$umap_rna@cell.embeddings[, 2]
    
    data_b <- as.data.frame(t(as.matrix(sc_b@assays$ADT@counts)))
    data_b <- log2(data_b + 1)
    data_b <- data_b[, order(colnames(data_b))]
    
    rna_data_b <- sc_b@assays$RNA@data
    
    Idents(sc) <- sc$sample_name
    sc <- subset(sc, idents = samples_to_subset)
    
    adtumap1 <- sc@reductions$umap_adt@cell.embeddings[, 1]
    adtumap2 <- sc@reductions$umap_adt@cell.embeddings[, 2]
    rnaumap1 <- sc@reductions$umap_rna@cell.embeddings[, 1]
    rnaumap2 <- sc@reductions$umap_rna@cell.embeddings[, 2]
    
    data <- as.data.frame(t(as.matrix(sc@assays$ADT@counts)))
    data <- log2(data + 1)
    data <- data[, order(colnames(data))]
    
    rna_data <- sc@assays$RNA@data
    proteinList <- rownames(sc@assays$ADT)
    geneList <- sc@assays$RNA@data@Dimnames[1]

    }
    return(list(
      data,
      adtumap1,
      adtumap2,
      rnaumap1,
      rnaumap2,
      sc@meta.data,
      rna_data,
      sc,
      geneList,
      proteinList,
      data_a,
      adtumap1_a,
      adtumap2_a,
      rnaumap1_a,
      rnaumap2_a,
      sc_a@meta.data,
      rna_data_a,
      data_b,
      adtumap1_b,
      adtumap2_b,
      rnaumap1_b,
      rnaumap2_b,
      sc_b@meta.data,
      rna_data_b,
      sc_a,
      sc_b
      
    ))
  })
  

    
    
  datasetInput_d <- datasetInput %>% debounce(100)
  
  
  
  
  
  
  
  
  
  
  
  
  observeEvent(input$selected_igt_id, {
    reset("ProteinSelector",asis = FALSE)
    
    session$resetBrush("ADT_brush")
    session$resetBrush("plot4_brush")
    session$resetBrush("facs1_brush")
    session$resetBrush("facs2_brush")
    session$resetBrush("facs3_brush")
    session$resetBrush("plot1_brush")
    session$resetBrush("plot4FC_brush")
    session$resetBrush("plot5_brush")
    session$resetBrush("plot6_brush")
    
    
    data = datasetInput_d()
    
    geneList = data[[9]]
    proteinList = data[[10]]
    
    updateSelectizeInput(
      session,
      inputId = "cluster",
      label = "Cluster color",
      choices = as.list(c("RNA_clusters","Protein_clusters", 'sample_name','rna/protein expression','list of cells')),
      selected = "RNA_clusters"
    )
    updateSelectizeInput(
      session,
      inputId = "cluster2",
      label = "Cluster color",
      choices = as.list(c("RNA_clusters","Protein_clusters", 'sample_name')),
      selected = "RNA_clusters"
    )
    updateSelectizeInput(
      session,
      inputId = "plot1_X",
      label = "X axis",
      choices = as.list(proteinList),
      selected = "TCRB"
    )
    updateSelectizeInput(
      session,
      inputId = "plot1_Y",
      label = "Y axis",
      choices = as.list(proteinList),
      selected = "TCRGD"
    )
    updateSelectizeInput(
      session,
      inputId = "plot2_X",
      label = "X axis",
      choices = as.list(proteinList),
      selected = "CD4"
    )
    updateSelectizeInput(
      session,
      inputId = "plot2_Y",
      label = "Y axis",
      choices = as.list(proteinList),
      selected = "CD8A"
    )
    updateSelectizeInput(
      session,
      inputId = "plot3_X",
      label = "X axis",
      choices = as.list(proteinList),
      selected = "CD44"
    )
    updateSelectizeInput(
      session,
      inputId = "plot3_Y",
      label = "Y axis",
      choices = as.list(proteinList),
      selected = "CD62L"
    )
    updateSelectizeInput(
      session,
      inputId = "plot4_X",
      label = "X axis",
      choices = as.list(proteinList),
      selected = "CD44"
    )
    updateSelectizeInput(
      session,
      inputId = "plot4_Y",
      label = "Y axis",
      choices = as.list(proteinList),
      selected = "CD62L"
    )
    
  })
  
  observeEvent(input$datatype, {

    session$resetBrush("plot5_brush")
    session$resetBrush("plot6_brush")
    session$resetBrush("ADT_brush_3")
  })
  
  observeEvent(input$cluster2, {
    
    session$resetBrush("plot5_brush")
    session$resetBrush("plot6_brush")
  })
  
  observeEvent(input$cluster_int, {
    session$resetBrush("ADT_brush_3")
    
  })
  
  observeEvent(input$selected_igt_id, {
    reset("ProteinSelector",asis = FALSE)
    
    session$resetBrush("ADT_brush")
    session$resetBrush("ADT_brush_3")
  })
  
  observeEvent(input$tabs, {
    reset("ProteinSelector",asis = FALSE)
    
    session$resetBrush("ADT_brush")
    session$resetBrush("ADT_brush_3")

    
  })
  
  observeEvent(input$cluster_int, {
    session$resetBrush("ADT_brush_3")
    
  })
  observeEvent(input$datatype2, {
    session$resetBrush("ADT_brush_3")
    
  })
  observeEvent(input$cluster, {
    
    if (input$cluster == 'rna/protein expression'){
      data = datasetInput_d()
      
      geneList = data[[9]]
      
      
      updateSelectizeInput(
        session,
        inputId = "gene",
        label = "Select Gene",
        selected = "Cd8a",
        choices = sort(unlist(geneList, recursive = TRUE, use.names = TRUE)),
        server= TRUE
      )
      
      
    }else {
    }
    
  })
  output$ProteinSelector <- renderUI({
    
    if (input$cluster == 'rna/protein expression'){
      data = datasetInput_d()
      
      proteinList = data[[10]]
      
      reset("protein",asis = FALSE)
      reset("ProteinSelector",asis = FALSE)
      
      
      #updateSelectizeInput(
      #session,
      selectizeInput("protein", "Select Protein", as.list(sort(c(proteinList))),
                     selected = "CD8A")
    }else {
    }
    
  })
  
  
  
  #for the link to download the antibody list used for the experiment
  output$downloadData <- downloadHandler(
    filename <- function() {
      "adt_hash_seq_file_2022_panel.csv"
    },
    content <- function(file) {
      write.csv(ADT_table, file)
    })
  
  #depending on the dataset, select the clusters needed to show on the main screen
  
  
  
  output$description <- renderText({
    description <- read.csv("descriptions.csv")
    desc <-
      description$Description[description$Dataset == input$selected_igt_id]
  })
  
  # When a double-click happens, check if there's a brush on the plot.
  # If so, zoom to the brush bounds; if not, reset the zoom.
  #OBSERVE EVENT ----
  observeEvent(input$plot1_dblclick, {
    brush <- input$plot1_brush
    if (!is.null(brush)) {
      ranges$x <- c(brush$xmin, brush$xmax)
      ranges$y <- c(brush$ymin, brush$ymax)
    } else {
      ranges$x <- NULL
      ranges$y <- NULL
    }
  })
  
  # -------------------------------------------------------------------
  # Linked plots (middle and right)
  facs1 <- reactiveValues(x = NULL, y = NULL)
  ranges4 <- reactiveValues(x = NULL, y = NULL)
  ADTrange <- reactiveValues(x = NULL, y = NULL)
  facs2 <- reactiveValues(x = NULL, y = NULL)
  facs3 <- reactiveValues(x = NULL, y = NULL)
  ranges5 <- reactiveValues(x = NULL, y = NULL)
  ranges6 <- reactiveValues(x = NULL, y = NULL)
  ADTrange_3 <- reactiveValues(x = NULL, y = NULL)
  
  selectedIGT <- reactive({
    filter(sample_metadata, igt_id == strsplit(input$selected_igt_id, "_")[[1]][1])
    
  })
  


  
  
  
  output$prot_plot1 <- renderPlot({
    
    data = datasetInput_d()
    
    data = data[[1]]
    
    ggplot(data, aes(data[, input$plot1_X], data[, input$plot1_Y])) + geom_point(col =
                                                                                   rgb(0, 0, 0, 0.1)) +
      xlab(input$plot1_X) + ggtitle(paste0(input$plot1_X, " vs ", input$plot1_Y, " (select gate A)")) +
      ylab(input$plot1_Y) + xlim(0, 10) + ylim(0, 10) + theme_bw()
  })
  
  output$prot_plot2 <- renderPlot({
    data <- datasetInput_d()
    data = data[[1]]
    
    if (!is.null(facs1$x)) {
      if (input$plot2_gate == "none") {
        data %>% ggplot(aes(data[, input$plot2_X], data[, input$plot2_Y])) + geom_point(col =
                                                                                          rgb(0, 0, 0, 0.1)) +
          xlab(input$plot2_X) + ggtitle(paste0(input$plot2_X, " vs ", input$plot2_Y, " (select gate B)")) +
          ylab(input$plot2_Y) +
          xlim(0, 10) + ylim(0, 10) + theme_bw()
      } else{
        if (input$plot2_gate == "A") {
          data2 <-
            data[which(data[, input$plot1_X] < facs1$x[2] &
                         data[, input$plot1_X] > facs1$x[1]),]
          data2 <-
            data2[which(data2[, input$plot1_Y] < facs1$y[2] &
                          data2[, input$plot1_Y] > facs1$y[1]),]
          
          if (input$highlight_filter == "Filter") {
            data2 %>% ggplot(aes(data2[, input$plot2_X], data2[, input$plot2_Y])) +
              geom_point(data = data2,
                         aes(data2[, input$plot2_X], data2[, input$plot2_Y]),
                         col = rgb(1, 0, 0, 0.1)) +
              xlab(input$plot2_X) + ggtitle(paste0(
                input$plot2_X,
                " vs ",
                input$plot2_Y,
                " (select gate B)"
              )) + ylab(input$plot2_Y) +
              xlim(0, 10) + ylim(0, 10) + theme_bw()
          } else{
            data2 %>% ggplot(aes(data2[, input$plot2_X], data2[, input$plot2_Y])) + geom_point(data =
                                                                                                 data,
                                                                                               aes(data[, input$plot2_X], data[, input$plot2_Y]),
                                                                                               col = "grey89") +
              geom_point(data = data2,
                         aes(data2[, input$plot2_X], data2[, input$plot2_Y]),
                         col = rgb(1, 0, 0, 0.1)) +
              xlab(input$plot2_X) + ggtitle(paste0(
                input$plot2_X,
                " vs ",
                input$plot2_Y,
                " (select gate B)"
              )) + ylab(input$plot2_Y) +
              xlim(0, 10) + ylim(0, 10) + theme_bw()
          }
        } else{
          if (input$plot2_gate == "B" & !is.null(facs2$x)) {
            data2 <-
              data[which(data[, input$plot2_X] < facs2$x[2] &
                           data[, input$plot2_X] > facs2$x[1]),]
            data2 <-
              data2[which(data2[, input$plot2_Y] < facs2$y[2] &
                            data2[, input$plot2_Y] > facs2$y[1]),]
            
            data2 %>% ggplot(aes(data2[, input$plot2_X], data2[, input$plot2_Y])) + geom_point(col =
                                                                                                 rgb(0, 0, 0, 0.1)) +
              #geom_point(data=data2,aes(data2[,input$ADTmarker],data2[,input$ADTmarker2[2]]),col="red")+
              xlab(input$plot2_X) + ggtitle(paste0(
                input$plot2_X,
                " vs ",
                input$plot2_Y,
                " (select gate B)"
              )) + ylab(input$plot2_Y) +
              xlim(0, 10) + ylim(0, 10) + theme_bw()
            
            
          } else{
            if (input$plot2_gate == "A+B" & !is.null(facs2$x)) {
              data2 <-
                data[which(data[, input$plot1_X] < facs1$x[2] &
                             data[, input$plot1_X] > facs1$x[1]),]
              data2 <-
                data2[which(data2[, input$plot1_Y] < facs1$y[2] &
                              data2[, input$plot1_Y] > facs1$y[1]),]
              
              data2 <-
                data2[which(data2[, input$plot2_X] < facs2$x[2] &
                              data2[, input$plot2_X] > facs2$x[1]),]
              data2 <-
                data2[which(data2[, input$plot2_Y] < facs2$y[2] &
                              data2[, input$plot2_Y] > facs2$y[1]),]
              
              data2 %>% ggplot(aes(data2[, input$plot2_X], data2[, input$plot2_Y])) + geom_point(col =
                                                                                                   rgb(0, 0, 0, 0.1)) +
                #geom_point(data=data2,aes(data2[,input$ADTmarker],data2[,input$ADTmarker2[2]]),col="red")+
                xlab(input$plot2_X) + ggtitle(paste0(
                  input$plot2_X,
                  " vs ",
                  input$plot2_Y,
                  " (select gate B)"
                )) + ylab(input$plot2_Y) +
                xlim(0, 10) + ylim(0, 10) + theme_bw()
              
            }
          }
        }
      }
    } else{
      if (input$plot2_X == "") {
        ggplot() + geom_point() + theme_bw()
      } else{
        ggplot(data, aes(data[, input$plot2_X], data[, input$plot2_Y])) + geom_point(col =
                                                                                       rgb(0, 0, 0, 0.1)) +
          xlab(input$plot2_X) + ggtitle(paste0(input$plot2_X, " vs ", input$plot2_Y, " (select gate B)")) +
          ylab(input$plot2_Y) + xlim(0, 10) + ylim(0, 10) + theme_bw()
      }
    }
    
  })
  
  output$prot_plot3 <- renderPlot({
    data <- datasetInput_d() 
    data = data[[1]]
    
    if (!is.null(facs2$x)) {
      if (input$plot3_gate == "none") {
        ggplot(data, aes(data[, input$plot3_X], data[, input$plot3_Y])) + geom_point(col =
                                                                                       rgb(0, 0, 0, 0.1)) +
          xlab(input$plot3_X) + ggtitle(paste0(input$plot3_X, " vs ", input$plot3_Y)) +
          ylab(input$plot3_Y) + xlim(0, 10) + ylim(0, 10) + theme_bw()
      } else{
        if (input$plot3_gate == "A") {
          #subset data to gate A
          data2 <-
            data[which(data[, input$plot1_X] < facs1$x[2] &
                         data[, input$plot1_X] > facs1$x[1]),]
          data2 <-
            data2[which(data2[, input$plot1_Y] < facs1$y[2] &
                          data2[input$plot1_Y] > facs1$y[1]),]
          if (input$highlight_filter2 == "Filter") {
            ggplot() +
              geom_point(data = data[rownames(data2),],
                         aes(data[rownames(data2), input$plot3_X], data[rownames(data2), input$plot3_Y]),
                         col = rgb(0, 0, 1, 0.1)) +
              xlab(input$plot3_X) + ggtitle(paste0(input$plot3_X, " vs ", input$plot3_Y)) +
              ylab(input$plot3_Y) + xlim(0, 10) + ylim(0, 10) + theme_bw()
          } else{
            ggplot() + geom_point(data = data,
                                  aes(data[, input$plot3_X], data[, input$plot3_Y]),
                                  col = "grey89") +
              geom_point(data = data[rownames(data2),],
                         aes(data[rownames(data2), input$plot3_X], data[rownames(data2), input$plot3_Y]),
                         col = rgb(0, 0, 1, 0.1)) +
              xlab(input$plot3_X) + ggtitle(paste0(input$plot3_X, " vs ", input$plot3_Y)) +
              ylab(input$plot3_Y) + xlim(0, 10) + ylim(0, 10) + theme_bw()
          }
        } else{
          if (input$plot3_gate == "B" && input$plot2_gate != "A") {
            #subset data to gate B
            data2 <-
              data[which(data[, input$plot2_X] < facs2$x[2] &
                           data[, input$plot2_X] > facs2$x[1]),]
            data2 <-
              data2[which(data2[, input$plot2_Y] < facs2$y[2] &
                            data2[input$plot2_Y] > facs2$y[1]),]
            if (input$highlight_filter2 == "Filter") {
              ggplot() +
                geom_point(data = data[rownames(data2),],
                           aes(data[rownames(data2), input$plot3_X], data[rownames(data2), input$plot3_Y]),
                           col = rgb(0, 0, 1, 0.1)) +
                xlab(input$plot3_X) + ggtitle(paste0(input$plot3_X, " vs ", input$plot3_Y)) +
                ylab(input$plot3_Y) + xlim(0, 10) + ylim(0, 10) + theme_bw()
            } else{
              ggplot() +  geom_point(data = data,
                                     aes(data[, input$plot3_X], data[, input$plot3_Y]),
                                     col = "grey89") +
                geom_point(data = data[rownames(data2),],
                           aes(data[rownames(data2), input$plot3_X], data[rownames(data2), input$plot3_Y]),
                           col = rgb(0, 0, 1, 0.1)) +
                xlab(input$plot3_X) + ggtitle(paste0(input$plot3_X, " vs ", input$plot3_Y)) +
                ylab(input$plot3_Y) + xlim(0, 10) + ylim(0, 10) + theme_bw()
            }
          } else{
            if ((input$plot3_gate == "A+B") ||
                (input$plot3_gate == "B" &&
                 input$plot2_gate == "A")) {
              if (!is.null(facs1$x)) {
                data2 <-
                  data[which(data[, input$plot1_X] < facs1$x[2] &
                               data[, input$plot1_X] > facs1$x[1]),]
                data2 <-
                  data2[which(data2[, input$plot1_Y] < facs1$y[2] &
                                data2[input$plot1_Y] > facs1$y[1]),]
                
                data3 <-
                  data2[which(data2[, input$plot2_X] < facs2$x[2] &
                                data2[, input$plot2_X] > facs2$x[1]),]
                data3 <-
                  data3[which(data3[, input$plot2_Y] < facs2$y[2] &
                                data3[input$plot2_Y] > facs2$y[1]),]
              } else{
                data3 <-
                  data[which(data[, input$plot2_X] < facs2$x[2] &
                               data[, input$plot2_X] > facs2$x[1]),]
                data3 <-
                  data3[which(data3[, input$plot2_Y] < facs2$y[2] &
                                data3[input$plot2_Y] > facs2$y[1]),]
              }
              
              if (input$highlight_filter2 == "Filter") {
                ggplot() +
                  geom_point(data = data[rownames(data3),],
                             aes(data[rownames(data3), input$plot3_X],
                                 data[rownames(data3), input$plot3_Y]),
                             col = rgb(1, 0, 1, 0.1)) +
                  geom_point(col = rgb(1, 0, 1, 0.5)) +
                  xlab(input$plot3_X) + ggtitle(paste0(input$plot3_X, " vs ", input$plot3_Y)) +
                  ylab(input$plot3_Y) + xlim(0, 10) + ylim(0, 10) +
                  theme_bw()
              } else{
                ggplot() + geom_point(data = data,
                                      aes(data[, input$plot3_X], data[, input$plot3_Y]),
                                      col = "grey89") +
                  geom_point(data = data[rownames(data3),],
                             aes(data[rownames(data3), input$plot3_X],
                                 data[rownames(data3), input$plot3_Y]),
                             col = rgb(1, 0, 1, 0.1)) +
                  geom_point(col = rgb(1, 0, 1, 0.5)) +
                  xlab(input$plot3_X) + ggtitle(paste0(input$plot3_X, " vs ", input$plot3_Y)) +
                  ylab(input$plot3_Y) + xlim(0, 10) + ylim(0, 10) +
                  theme_bw()
              }
            }
          }
        }
      }
    } else{
      if (input$plot3_X == "") {
        ggplot() + geom_point() + theme_bw()
      } else{
        ggplot(data, aes(data[, input$plot3_X], data[, input$plot3_Y])) + geom_point(col =
                                                                                       rgb(0, 0, 0, 0.1)) +
          xlab(input$plot3_X) + ggtitle(paste0(input$plot3_X, " vs ", input$plot3_Y)) +
          ylab(input$plot3_Y) + xlim(0, 10) + ylim(0, 10) + theme_bw()
      }
    }
    
    
  })
  output$prot_plot4 <- renderPlot({
    data <- datasetInput_d() 
    data = data[[1]]
    
    if (!is.null(facs3$x)) {
      if (input$plot4_gate == "none") {
        ggplot(data, aes(data[, input$plot4_X], data[, input$plot4_Y])) + geom_point(col =
                                                                                       rgb(0, 0, 0, 0.1)) +
          xlab(input$plot4_X) + ggtitle(paste0(input$plot4_X, " vs ", input$plot4_Y)) +
          ylab(input$plot4_Y) + xlim(0, 10) + ylim(0, 10) + theme_bw()
      } else{
        if (input$plot4_gate == "A") {
          #subset data to gate A
          data2 <-
            data[which(data[, input$plot1_X] < facs1$x[2] &
                         data[, input$plot1_X] > facs1$x[1]),]
          data2 <-
            data2[which(data2[, input$plot1_Y] < facs1$y[2] &
                          data2[input$plot1_Y] > facs1$y[1]),]
          if (input$highlight_filter3 == "Filter") {
            ggplot() +
              geom_point(data = data[rownames(data2),],
                         aes(data[rownames(data2), input$plot4_X], data[rownames(data2), input$plot4_Y]),
                         col = rgb(0, 0, 1, 0.1)) +
              xlab(input$plot4_X) + ggtitle(paste0(input$plot4_X, " vs ", input$plot4_Y)) +
              ylab(input$plot4_Y) + xlim(0, 10) + ylim(0, 10) + theme_bw()
          } else{
            ggplot() + geom_point(data = data,
                                  aes(data[, input$plot4_X], data[, input$plot4_Y]),
                                  col = "grey89") +
              geom_point(data = data[rownames(data2),],
                         aes(data[rownames(data2), input$plot4_X], data[rownames(data2), input$plot4_Y]),
                         col = rgb(0, 0, 1, 0.1)) +
              xlab(input$plot4_X) + ggtitle(paste0(input$plot4_X, " vs ", input$plot4_Y)) +
              ylab(input$plot4_Y) + xlim(0, 10) + ylim(0, 10) + theme_bw()
          }
        } else{
          if (input$plot4_gate == "B" && input$plot2_gate != "A" ) {
            #subset data to gate B
            data2 <-
              data[which(data[, input$plot2_X] < facs2$x[2] &
                           data[, input$plot2_X] > facs2$x[1]),]
            data2 <-
              data2[which(data2[, input$plot2_Y] < facs2$y[2] &
                            data2[input$plot2_Y] > facs2$y[1]),]
            if (input$highlight_filter3 == "Filter") {
              ggplot() +
                geom_point(data = data[rownames(data2),],
                           aes(data[rownames(data2), input$plot4_X], data[rownames(data2), input$plot4_Y]),
                           col = rgb(0, 0, 1, 0.1)) +
                xlab(input$plot4_X) + ggtitle(paste0(input$plot4_X, " vs ", input$plot4_Y)) +
                ylab(input$plot4_Y) + xlim(0, 10) + ylim(0, 10) + theme_bw()
            } else{
              ggplot() +  geom_point(data = data,
                                     aes(data[, input$plot4_X], data[, input$plot4_Y]),
                                     col = "grey89") +
                geom_point(data = data[rownames(data2),],
                           aes(data[rownames(data2), input$plot4_X], data[rownames(data2), input$plot4_Y]),
                           col = rgb(0, 0, 1, 0.1)) +
                xlab(input$plot4_X) + ggtitle(paste0(input$plot4_X, " vs ", input$plot4_Y)) +
                ylab(input$plot4_Y) + xlim(0, 10) + ylim(0, 10) + theme_bw()
            }
          } else{
            if (input$plot4_gate == "B" && input$plot2_gate == "A" ) {
              
              data2 <-
                data[which(data[, input$plot1_X] < facs1$x[2] &
                             data[, input$plot1_X] > facs1$x[1]),]
              data2 <-
                data2[which(data2[, input$plot1_Y] < facs1$y[2] &
                              data2[input$plot1_Y] > facs1$y[1]),]
              
              data3 <-
                data2[which(data2[, input$plot2_X] < facs2$x[2] &
                              data2[, input$plot2_X] > facs2$x[1]),]
              data3 <-
                data3[which(data3[, input$plot2_Y] < facs2$y[2] &
                              data3[input$plot2_Y] > facs2$y[1]),]
              
              
              if (input$highlight_filter3 == "Filter") {
                ggplot() +
                  geom_point(data = data[rownames(data3),],
                             aes(data[rownames(data3), input$plot4_X], data[rownames(data3), input$plot4_Y]),
                             col = rgb(0, 0, 1, 0.1)) +
                  xlab(input$plot4_X) + ggtitle(paste0(input$plot4_X, " vs ", input$plot4_Y)) +
                  ylab(input$plot4_Y) + xlim(0, 10) + ylim(0, 10) + theme_bw()
              } else{
                ggplot() +  geom_point(data = data,
                                       aes(data[, input$plot4_X], data[, input$plot4_Y]),
                                       col = "grey89") +
                  geom_point(data = data[rownames(data3),],
                             aes(data[rownames(data3), input$plot4_X], data[rownames(data3), input$plot4_Y]),
                             col = rgb(0, 0, 1, 0.1)) +
                  xlab(input$plot4_X) + ggtitle(paste0(input$plot4_X, " vs ", input$plot4_Y)) +
                  ylab(input$plot4_Y) + xlim(0, 10) + ylim(0, 10) + theme_bw()
              }
              
            } else{
              
              if (input$plot4_gate == "C" && input$plot3_gate == "none") {
                #subset data to gate C
                data2 <-
                  data[which(data[, input$plot3_X] < facs3$x[2] &
                               data[, input$plot3_X] > facs3$x[1]),]
                data2 <-
                  data2[which(data2[, input$plot3_Y] < facs3$y[2] &
                                data2[input$plot3_Y] > facs3$y[1]),]
                if (input$highlight_filter3 == "Filter") {
                  ggplot() +
                    geom_point(data = data[rownames(data2),],
                               aes(data[rownames(data2), input$plot4_X], data[rownames(data2), input$plot4_Y]),
                               col = rgb(0, 0, 1, 0.1)) +
                    xlab(input$plot4_X) + ggtitle(paste0(input$plot4_X, " vs ", input$plot4_Y)) +
                    ylab(input$plot4_Y) + xlim(0, 10) + ylim(0, 10) + theme_bw()
                } else{
                  ggplot() +  geom_point(data = data,
                                         aes(data[, input$plot4_X], data[, input$plot4_Y]),
                                         col = "grey89") +
                    geom_point(data = data[rownames(data2),],
                               aes(data[rownames(data2), input$plot4_X], data[rownames(data2), input$plot4_Y]),
                               col = rgb(0, 0, 1, 0.1)) +
                    xlab(input$plot4_X) + ggtitle(paste0(input$plot4_X, " vs ", input$plot4_Y)) +
                    ylab(input$plot4_Y) + xlim(0, 10) + ylim(0, 10) + theme_bw()
                }
              }
              
              
              else{
                if (input$plot4_gate == "C" && input$plot3_gate == "A") {
                  
                  data2 <-
                    data[which(data[, input$plot1_X] < facs1$x[2] &
                                 data[, input$plot1_X] > facs1$x[1]),]
                  data2 <-
                    data2[which(data2[, input$plot1_Y] < facs1$y[2] &
                                  data2[input$plot1_Y] > facs1$y[1]),]
                  
                  data3 <-
                    data2[which(data2[, input$plot3_X] < facs3$x[2] &
                                  data2[, input$plot3_X] > facs3$x[1]),]
                  data3 <-
                    data3[which(data3[, input$plot3_Y] < facs3$y[2] &
                                  data3[input$plot3_Y] > facs3$y[1]),]
                  
                  
                  if (input$highlight_filter3 == "Filter") {
                    ggplot() +
                      geom_point(data = data[rownames(data3),],
                                 aes(data[rownames(data3), input$plot4_X], data[rownames(data3), input$plot4_Y]),
                                 col = rgb(0, 0, 1, 0.1)) +
                      xlab(input$plot4_X) + ggtitle(paste0(input$plot4_X, " vs ", input$plot4_Y)) +
                      ylab(input$plot4_Y) + xlim(0, 10) + ylim(0, 10) + theme_bw()
                  } else{
                    ggplot() +  geom_point(data = data,
                                           aes(data[, input$plot4_X], data[, input$plot4_Y]),
                                           col = "grey89") +
                      geom_point(data = data[rownames(data3),],
                                 aes(data[rownames(data3), input$plot4_X], data[rownames(data3), input$plot4_Y]),
                                 col = rgb(0, 0, 1, 0.1)) +
                      xlab(input$plot4_X) + ggtitle(paste0(input$plot4_X, " vs ", input$plot4_Y)) +
                      ylab(input$plot4_Y) + xlim(0, 10) + ylim(0, 10) + theme_bw()
                  }
                  
                }
                
                
                
                
                
                
                else{
                  
                  
                  if ((input$plot4_gate == "A+B+C") ||
                      ( input$plot4_gate == "C" &&
                        input$plot3_gate == "B" &&
                        input$plot2_gate == "A") ||
                      ( input$plot4_gate == "C" &&
                        input$plot3_gate == "A+B" &&
                        input$plot2_gate == "A")) 
                  {
                    if (!is.null(facs1$x)) {
                      data2 <-
                        data[which(data[, input$plot1_X] < facs1$x[2] &
                                     data[, input$plot1_X] > facs1$x[1]),]
                      data2 <-
                        data2[which(data2[, input$plot1_Y] < facs1$y[2] &
                                      data2[input$plot1_Y] > facs1$y[1]),]
                      
                      data3 <-
                        data2[which(data2[, input$plot2_X] < facs2$x[2] &
                                      data2[, input$plot2_X] > facs2$x[1]),]
                      data3 <-
                        data3[which(data3[, input$plot2_Y] < facs2$y[2] &
                                      data3[input$plot2_Y] > facs2$y[1]),]
                      data4 <-
                        data3[which(data3[, input$plot3_X] < facs3$x[2] &
                                      data3[, input$plot3_X] > facs3$x[1]),]
                      data4 <-
                        data4[which(data4[, input$plot3_Y] < facs3$y[2] &
                                      data4[input$plot3_Y] > facs3$y[1]),]
                    } else{
                      data3 <-
                        data[which(data[, input$plot2_X] < facs2$x[2] &
                                     data[, input$plot2_X] > facs2$x[1]),]
                      data3 <-
                        data3[which(data3[, input$plot2_Y] < facs2$y[2] &
                                      data3[input$plot2_Y] > facs2$y[1]),]
                      data4 <-
                        data3[which(data3[, input$plot3_X] < facs3$x[2] &
                                      data3[, input$plot3_X] > facs3$x[1]),]
                      data4 <-
                        data4[which(data4[, input$plot3_Y] < facs3$y[2] &
                                      data4[input$plot3_Y] > facs3$y[1]),]
                    }
                    
                    if (input$highlight_filter3 == "Filter") {
                      ggplot() +
                        geom_point(data = data[rownames(data4),],
                                   aes(data[rownames(data4), input$plot4_X],
                                       data[rownames(data4), input$plot4_Y]),
                                   col = rgb(0, 1, 0.5, 0.3)) +
                        geom_point(col = rgb(1, 0, 1, 0.5)) +
                        xlab(input$plot4_X) + ggtitle(paste0(input$plot4_X, " vs ", input$plot4_Y)) +
                        ylab(input$plot4_Y) + xlim(0, 10) + ylim(0, 10) +
                        theme_bw()
                    } else{
                      ggplot() + geom_point(data = data,
                                            aes(data[, input$plot4_X], data[, input$plot4_Y]),
                                            col = "grey89") +
                        geom_point(data = data[rownames(data4),],
                                   aes(data[rownames(data4), input$plot4_X],
                                       data[rownames(data4), input$plot4_Y]),
                                   col = rgb(0, 1, 0.5, 0.3)) +
                        geom_point(col = rgb(1, 0, 1, 0.5)) +
                        xlab(input$plot4_X) + ggtitle(paste0(input$plot4_X, " vs ", input$plot4_Y)) +
                        ylab(input$plot4_Y) + xlim(0, 10) + ylim(0, 10) +
                        theme_bw()
                    }
                  }
                }
              }
            }
          }
        }
      }
      
    } else{
      if (input$plot4_X == "") {
        ggplot() + geom_point() + theme_bw()
      } else{
        ggplot(data, aes(data[, input$plot4_X], data[, input$plot4_Y])) + geom_point(col =
                                                                                       rgb(0, 0, 0, 0.1)) +
          xlab(input$plot4_X) + ggtitle(paste0(input$plot4_X, " vs ", input$plot4_Y)) +
          ylab(input$plot4_Y) + xlim(0, 10) + ylim(0, 10) + theme_bw()
      }
    }
    
    
  })
  
  output$rnaumap <- renderPlot({
    data_list <- datasetInput_d() 
    data = data_list[[1]]
    umaps = data.frame(
      adtumap_1 = data_list[[2]],
      adtumap_2 = data_list[[3]],
      rnaumap_1 = data_list[[4]],
      rnaumap_2 = data_list[[5]]
    )
    
    
    if (input$tabs == "facs") {
      if (is.null(facs1$x) & is.null(facs2$x) & is.null(facs3$x) ) {
        output$downloadCells <- downloadHandler(
          filename = function() {
            paste("Cell_Selection_Warning.csv")
          },
          content = function(file) {
            write.csv(as.data.frame(c('No Cells Selected')), file,row.names = F)
          }
        )
        g1 <-
          ggplot(umaps, aes(rnaumap_1, rnaumap_2)) + geom_point(size = 0.5, col =
                                                                  "grey89") +
          ggtitle("RNA UMAP") + theme_bw() + theme(plot.title = element_text(hjust = 0.5, size =
                                                                               20)) + removeGrid()
        g2 <-
          ggplot(umaps, aes(adtumap_1, adtumap_2)) + geom_point(size = 0.5, col =
                                                                  "grey89") +
          ggtitle("Protein UMAP") + theme_bw() + theme(plot.title = element_text(hjust = 0.5, size =
                                                                                   20)) + removeGrid()
        g1 + g2
        
      } 
      else{
        if (((!is.null(facs2$x) && is.null(facs1$x) && is.null(facs3$x))) || ((!is.null(facs1$x) && is.null(facs2$x) && is.null(facs3$x))) || ((!is.null(facs3$x) && is.null(facs1$x) && is.null(facs2$x))))  {
          if (is.null(facs1$x) && is.null(facs3$x) ) {
            data2 <-
              data[which(data[, input$plot2_X] < facs2$x[2] &
                           data[, input$plot2_X] > facs2$x[1]),]
            data2 <-
              data2[which(data2[, input$plot2_Y] < facs2$y[2] &
                            data2[, input$plot2_Y] > facs2$y[1]),]
          } else{
            if (is.null(facs2$x) && is.null(facs3$x) ) {
              data2 <-
                data[which(data[, input$plot1_X] < facs1$x[2] &
                             data[, input$plot1_X] > facs1$x[1]),]
              data2 <-
                data2[which(data2[, input$plot1_Y] < facs1$y[2] &
                              data2[, input$plot1_Y] > facs1$y[1]),]
            }
            else{
              if (is.null(facs1$x) && is.null(facs2$x) ) {
                data2 <-
                  data[which(data[, input$plot3_X] < facs3$x[2] &
                               data[, input$plot3_X] > facs3$x[1]),]
                data2 <-
                  data2[which(data2[, input$plot3_Y] < facs3$y[2] &
                                data2[, input$plot3_Y] > facs3$y[1]),]
              }
            }
          }
          
          output$downloadCells <- downloadHandler(
            filename = function() {
              paste("SelectedCells_FC_plots.csv")
            },
            content = function(file) {
              write.csv(as.data.frame(rownames(data2)), file,row.names = F)
            }
          )
          g1 <-
            ggplot(umaps, aes(rnaumap_1, rnaumap_2)) + geom_point(size = 0.5, col =
                                                                    "grey89") +
            geom_point(
              data = umaps[rownames(data2),],
              aes(rnaumap_1, rnaumap_2),
              col = rgb(1, 0, 0),
              size = 0.5
            ) +
            ggtitle("RNA UMAP") + theme_bw() + theme(plot.title = element_text(hjust = 0.5, size =
                                                                                 20)) + removeGrid()
          g2 <-
            ggplot(umaps, aes(adtumap_1, adtumap_2)) + geom_point(size = 0.5, col =
                                                                    "grey89") +
            geom_point(
              data = umaps[rownames(data2),],
              aes(adtumap_1, adtumap_2),
              col = rgb(1, 0, 0),
              size = 0.5
            ) +
            ggtitle("Protein UMAP") + theme_bw() + theme(plot.title = element_text(hjust = 0.5, size =
                                                                                     20)) + removeGrid()
          g1 + g2
          
          
          
        } 
        else{
          if (((!is.null(facs2$x) && !is.null(facs1$x) && is.null(facs3$x))) || ((!is.null(facs1$x) && is.null(facs2$x) && is.null(!facs3$x))) || ((!is.null(facs3$x) && is.null(facs1$x) && !is.null(facs2$x))))  {
            if (!is.null(facs2$x) && !is.null(facs1$x) ) {
              data2 <-
                data[which(data[, input$plot1_X] < facs1$x[2] &
                             data[, input$plot1_X] > facs1$x[1]),]
              data2 <-
                data2[which(data2[, input$plot1_Y] < facs1$y[2] &
                              data2[, input$plot1_Y] > facs1$y[1]),]
              
              data3 <-
                data2[which(data2[, input$plot2_X] < facs2$x[2] &
                              data2[, input$plot2_X] > facs2$x[1]),]
              data3 <-
                data3[which(data3[, input$plot2_Y] < facs2$y[2] &
                              data3[input$plot2_Y] > facs2$y[1]),]
              
            } else{
              if (!is.null(facs2$x) && !is.null(facs3$x) ) {
                data2 <-
                  data[which(data[, input$plot2_X] < facs2$x[2] &
                               data[, input$plot2_X] > facs2$x[1]),]
                data2 <-
                  data2[which(data2[, input$plot2_Y] < facs2$y[2] &
                                data2[, input$plot2_Y] > facs2$y[1]),]
                
                data3 <-
                  data2[which(data2[, input$plot3_X] < facs3$x[2] &
                                data2[, input$plot3_X] > facs3$x[1]),]
                data3 <-
                  data3[which(data3[, input$plot3_Y] < facs3$y[2] &
                                data3[input$plot3_Y] > facs3$y[1]),]
                
              }
              else{
                if (!is.null(facs1$x) && !is.null(facs3$x) ) {
                  data2 <-
                    data[which(data[, input$plot1_X] < facs1$x[2] &
                                 data[, input$plot1_X] > facs1$x[1]),]
                  data2 <-
                    data2[which(data2[, input$plot1_Y] < facs1$y[2] &
                                  data2[, input$plot1_Y] > facs1$y[1]),]
                  
                  data3 <-
                    data2[which(data2[, input$plot3_X] < facs3$x[2] &
                                  data2[, input$plot3_X] > facs3$x[1]),]
                  data3 <-
                    data3[which(data3[, input$plot3_Y] < facs3$y[2] &
                                  data3[input$plot3_Y] > facs3$y[1]),]
                  
                }
              }
            }
            
            if (input$UMAP_gate == "none") {
              output$downloadCells <- downloadHandler(
                filename = function() {
                  paste("Cell_Selection_Warning.csv")
                },
                content = function(file) {
                  write.csv(as.data.frame(c('No Cells Selected')), file,row.names = F)
                }
              )
              g1 <-
                ggplot(umaps, aes(rnaumap_1, rnaumap_2)) + geom_point(size = 0.5, col =
                                                                        "grey89") +
                ggtitle("RNA UMAP") + theme_bw() + theme(plot.title = element_text(hjust = 0.5, size =
                                                                                     20)) + removeGrid()
              g2 <-
                ggplot(umaps, aes(adtumap_1, adtumap_2)) + geom_point(size = 0.5, col =
                                                                        "grey89") +
                ggtitle("Protein UMAP") + theme_bw() + theme(plot.title = element_text(hjust = 0.5, size =
                                                                                         20)) + removeGrid()
              g1 + g2
              
            } else{
              if (input$UMAP_gate == "A") {
                data2 <-
                  data[which(data[, input$plot1_X] < facs1$x[2] &
                               data[, input$plot1_X] > facs1$x[1]),]
                data2 <-
                  data2[which(data2[, input$plot1_Y] < facs1$y[2] &
                                data2[, input$plot1_Y] > facs1$y[1]),]
                
                output$downloadCells <- downloadHandler(
                  filename = function() {
                    paste("SelectedCells_FC_plots.csv")
                  },
                  content = function(file) {
                    write.csv(as.data.frame(rownames(data2)), file,row.names = F)
                  }
                )
                g1 <-
                  ggplot(umaps, aes(rnaumap_1, rnaumap_2)) + geom_point(size = 0.5, col =
                                                                          "grey89") +
                  geom_point(
                    data = umaps[rownames(data2),],
                    aes(rnaumap_1, rnaumap_2),
                    col = rgb(0, 0, 1),
                    size = 0.5
                  ) +
                  ggtitle("RNA UMAP") + theme_bw() + theme(plot.title = element_text(hjust = 0.5, size =
                                                                                       20)) + removeGrid()
                g2 <-
                  ggplot(umaps, aes(adtumap_1, adtumap_2)) + geom_point(size = 0.5, col =
                                                                          "grey89") +
                  geom_point(
                    data = umaps[rownames(data2),],
                    aes(adtumap_1, adtumap_2),
                    col = rgb(0, 0, 1),
                    size = 0.5
                  ) +
                  ggtitle("Protein UMAP") + theme_bw() + theme(plot.title = element_text(hjust = 0.5, size =
                                                                                           20)) + removeGrid()
                g1 + g2
                
              } else{
                if (input$UMAP_gate == "B") {
                  data2 <-
                    data[which(data[, input$plot2_X] < facs2$x[2] &
                                 data[, input$plot2_X] > facs2$x[1]),]
                  data2 <-
                    data2[which(data2[, input$plot2_Y] < facs2$y[2] &
                                  data2[input$plot2_Y] > facs2$y[1]),]
                  output$downloadCells <- downloadHandler(
                    filename = function() {
                      paste("SelectedCells_FC_plots.csv")
                    },
                    content = function(file) {
                      write.csv(as.data.frame(rownames(data2)), file,row.names = F)
                    }
                  )
                  
                  g1 <-
                    ggplot(umaps, aes(rnaumap_1, rnaumap_2)) + geom_point(size = 0.5, col =
                                                                            "grey89") +
                    geom_point(
                      data = umaps[rownames(data2),],
                      aes(rnaumap_1, rnaumap_2),
                      col = rgb(1, 0, 0),
                      size = 0.5
                    ) +
                    ggtitle("RNA UMAP") + theme_bw() + theme(plot.title = element_text(hjust = 0.5, size =
                                                                                         20)) + removeGrid()
                  g2 <-
                    ggplot(umaps, aes(adtumap_1, adtumap_2)) + geom_point(size = 0.5, col =
                                                                            "grey89") +
                    geom_point(
                      data = umaps[rownames(data2),],
                      aes(adtumap_1, adtumap_2),
                      col = rgb(1, 0, 0),
                      size = 0.5
                    ) +
                    ggtitle("Protein UMAP") + theme_bw() + theme(plot.title = element_text(hjust = 0.5, size =
                                                                                             20)) + removeGrid()
                  g1 + g2
                  
                } 
                else{
                  if (input$UMAP_gate == "C") {
                    data2 <-
                      data[which(data[, input$plot3_X] < facs3$x[2] &
                                   data[, input$plot3_X] > facs3$x[1]),]
                    data2 <-
                      data2[which(data2[, input$plot3_Y] < facs3$y[2] &
                                    data2[input$plot3_Y] > facs3$y[1]),]
                    output$downloadCells <- downloadHandler(
                      filename = function() {
                        paste("SelectedCells_FC_plots.csv")
                      },
                      content = function(file) {
                        write.csv(as.data.frame(rownames(data2)), file,row.names = F)
                      }
                    )
                    
                    g1 <-
                      ggplot(umaps, aes(rnaumap_1, rnaumap_2)) + geom_point(size = 0.5, col =
                                                                              "grey89") +
                      geom_point(
                        data = umaps[rownames(data2),],
                        aes(rnaumap_1, rnaumap_2),
                        col = rgb(1, 0, 1, 0.1),
                        size = 0.5
                      ) +
                      ggtitle("RNA UMAP") + theme_bw() + theme(plot.title = element_text(hjust = 0.5, size =
                                                                                           20)) + removeGrid()
                    g2 <-
                      ggplot(umaps, aes(adtumap_1, adtumap_2)) + geom_point(size = 0.5, col =
                                                                              "grey89") +
                      geom_point(
                        data = umaps[rownames(data2),],
                        aes(adtumap_1, adtumap_2),
                        col = rgb(1, 0, 1, 0.1),
                        size = 0.5
                      ) +
                      ggtitle("Protein UMAP") + theme_bw() + theme(plot.title = element_text(hjust = 0.5, size =
                                                                                               20)) + removeGrid()
                    g1 + g2
                    
                    
                  } 
                  else{
                    
                    if (input$UMAP_gate == "A+B+C") {
                      
                      output$downloadCells <- downloadHandler(
                        filename = function() {
                          paste("SelectedCells_FC_plots.csv")
                        },
                        content = function(file) {
                          write.csv(as.data.frame(rownames(data3)), file,row.names = F)
                        }
                      )
                      g1 <-
                        ggplot(umaps, aes(rnaumap_1, rnaumap_2)) + geom_point(size = 0.5, col =
                                                                                "grey89") +
                        geom_point(
                          data = umaps[rownames(data3),],
                          aes(rnaumap_1, rnaumap_2),
                          col = rgb(1, 0, 1, 0.5),
                          size = 0.5
                        ) +
                        ggtitle("RNA UMAP") + theme_bw() + theme(plot.title = element_text(hjust = 0.5, size =
                                                                                             20)) + removeGrid()
                      g2 <-
                        ggplot(umaps, aes(adtumap_1, adtumap_2)) + geom_point(size = 0.5, col =
                                                                                "grey89") +
                        geom_point(
                          data = umaps[rownames(data3),],
                          aes(adtumap_1, adtumap_2),
                          col = rgb(1, 0, 1, 0.5),
                          size = 0.5
                        ) +
                        ggtitle("Protein UMAP") + theme_bw() + theme(plot.title = element_text(hjust = 0.5, size =
                                                                                                 20)) + removeGrid()
                      g1 + g2
                      
                      
                    }
                  }
                }
              }
            }
          }
          else{
            
            if (!is.null(facs1$x) & !is.null(facs2$x) & !is.null(facs3$x)) {
              data2 <-
                data[which(data[, input$plot1_X] < facs1$x[2] &
                             data[, input$plot1_X] > facs1$x[1]),]
              data2 <-
                data2[which(data2[, input$plot1_Y] < facs1$y[2] &
                              data2[, input$plot1_Y] > facs1$y[1]),]
              
              data3 <-
                data2[which(data2[, input$plot2_X] < facs2$x[2] &
                              data2[, input$plot2_X] > facs2$x[1]),]
              data3 <-
                data3[which(data3[, input$plot2_Y] < facs2$y[2] &
                              data3[input$plot2_Y] > facs2$y[1]),]
              
              data4 <-
                data3[which(data3[, input$plot3_X] < facs3$x[2] &
                              data3[, input$plot3_X] > facs3$x[1]),]
              data4 <-
                data4[which(data4[, input$plot3_Y] < facs3$y[2] &
                              data4[input$plot3_Y] > facs3$y[1]),]
              
              
              if (input$UMAP_gate == "none") {
                output$downloadCells <- downloadHandler(
                  filename = function() {
                    paste("Cell_Selection_Warning.csv")
                  },
                  content = function(file) {
                    write.csv(as.data.frame(c('No Cells Selected')), file,row.names = F)
                  }
                )
                g1 <-
                  ggplot(umaps, aes(rnaumap_1, rnaumap_2)) + geom_point(size = 0.5, col =
                                                                          "grey89") +
                  ggtitle("RNA UMAP") + theme_bw() + theme(plot.title = element_text(hjust = 0.5, size =
                                                                                       20)) + removeGrid()
                g2 <-
                  ggplot(umaps, aes(adtumap_1, adtumap_2)) + geom_point(size = 0.5, col =
                                                                          "grey89") +
                  ggtitle("Protein UMAP") + theme_bw() + theme(plot.title = element_text(hjust = 0.5, size =
                                                                                           20)) + removeGrid()
                g1 + g2
                
              } else{
                if (input$UMAP_gate == "A") {
                  data2 <-
                    data[which(data[, input$plot1_X] < facs1$x[2] &
                                 data[, input$plot1_X] > facs1$x[1]),]
                  data2 <-
                    data2[which(data2[, input$plot1_Y] < facs1$y[2] &
                                  data2[, input$plot1_Y] > facs1$y[1]),]
                  
                  output$downloadCells <- downloadHandler(
                    filename = function() {
                      paste("SelectedCells_FC_plots.csv")
                    },
                    content = function(file) {
                      write.csv(as.data.frame(rownames(data2)), file,row.names = F)
                    }
                  )
                  g1 <-
                    ggplot(umaps, aes(rnaumap_1, rnaumap_2)) + geom_point(size = 0.5, col =
                                                                            "grey89") +
                    geom_point(
                      data = umaps[rownames(data2),],
                      aes(rnaumap_1, rnaumap_2),
                      col = rgb(0, 0, 1),
                      size = 0.5
                    ) +
                    ggtitle("RNA UMAP") + theme_bw() + theme(plot.title = element_text(hjust = 0.5, size =
                                                                                         20)) + removeGrid()
                  g2 <-
                    ggplot(umaps, aes(adtumap_1, adtumap_2)) + geom_point(size = 0.5, col =
                                                                            "grey89") +
                    geom_point(
                      data = umaps[rownames(data2),],
                      aes(adtumap_1, adtumap_2),
                      col = rgb(0, 0, 1),
                      size = 0.5
                    ) +
                    ggtitle("Protein UMAP") + theme_bw() + theme(plot.title = element_text(hjust = 0.5, size =
                                                                                             20)) + removeGrid()
                  g1 + g2
                  
                } else{
                  if (input$UMAP_gate == "B") {
                    data2 <-
                      data[which(data[, input$plot2_X] < facs2$x[2] &
                                   data[, input$plot2_X] > facs2$x[1]),]
                    data2 <-
                      data2[which(data2[, input$plot2_Y] < facs2$y[2] &
                                    data2[input$plot2_Y] > facs2$y[1]),]
                    output$downloadCells <- downloadHandler(
                      filename = function() {
                        paste("SelectedCells_FC_plots.csv")
                      },
                      content = function(file) {
                        write.csv(as.data.frame(rownames(data2)), file,row.names = F)
                      }
                    )
                    
                    g1 <-
                      ggplot(umaps, aes(rnaumap_1, rnaumap_2)) + geom_point(size = 0.5, col =
                                                                              "grey89") +
                      geom_point(
                        data = umaps[rownames(data2),],
                        aes(rnaumap_1, rnaumap_2),
                        col = rgb(1, 0, 0),
                        size = 0.5
                      ) +
                      ggtitle("RNA UMAP") + theme_bw() + theme(plot.title = element_text(hjust = 0.5, size =
                                                                                           20)) + removeGrid()
                    g2 <-
                      ggplot(umaps, aes(adtumap_1, adtumap_2)) + geom_point(size = 0.5, col =
                                                                              "grey89") +
                      geom_point(
                        data = umaps[rownames(data2),],
                        aes(adtumap_1, adtumap_2),
                        col = rgb(1, 0, 0),
                        size = 0.5
                      ) +
                      ggtitle("Protein UMAP") + theme_bw() + theme(plot.title = element_text(hjust = 0.5, size =
                                                                                               20)) + removeGrid()
                    g1 + g2
                    
                  } 
                  else{
                    if (input$UMAP_gate == "C") {
                      data2 <-
                        data[which(data[, input$plot3_X] < facs3$x[2] &
                                     data[, input$plot3_X] > facs3$x[1]),]
                      data2 <-
                        data2[which(data2[, input$plot3_Y] < facs3$y[2] &
                                      data2[input$plot3_Y] > facs3$y[1]),]
                      output$downloadCells <- downloadHandler(
                        filename = function() {
                          paste("SelectedCells_FC_plots.csv")
                        },
                        content = function(file) {
                          write.csv(as.data.frame(rownames(data2)), file,row.names = F)
                        }
                      )
                      
                      g1 <-
                        ggplot(umaps, aes(rnaumap_1, rnaumap_2)) + geom_point(size = 0.5, col =
                                                                                "grey89") +
                        geom_point(
                          data = umaps[rownames(data2),],
                          aes(rnaumap_1, rnaumap_2),
                          col = rgb(1, 0, 1, 0.5),
                          size = 0.5
                        ) +
                        ggtitle("RNA UMAP") + theme_bw() + theme(plot.title = element_text(hjust = 0.5, size =
                                                                                             20)) + removeGrid()
                      g2 <-
                        ggplot(umaps, aes(adtumap_1, adtumap_2)) + geom_point(size = 0.5, col =
                                                                                "grey89") +
                        geom_point(
                          data = umaps[rownames(data2),],
                          aes(adtumap_1, adtumap_2),
                          col = rgb(1, 0, 1, 0.5),
                          size = 0.5
                        ) +
                        ggtitle("Protein UMAP") + theme_bw() + theme(plot.title = element_text(hjust = 0.5, size =
                                                                                                 20)) + removeGrid()
                      g1 + g2
                      
                      
                      
                      
                      
                      
                    } else{
                      if (input$UMAP_gate == "A+B") {
                        if (!is.null(facs1$x)) {
                          data2 <-
                            data[which(data[, input$plot1_X] < facs1$x[2] &
                                         data[, input$plot1_X] > facs1$x[1]),]
                          data2 <-
                            data2[which(data2[, input$plot1_Y] < facs1$y[2] &
                                          data2[input$plot1_Y] > facs1$y[1]),]
                          
                          data3 <-
                            data2[which(data2[, input$plot2_X] < facs2$x[2] &
                                          data2[, input$plot2_X] > facs2$x[1]),]
                          data3 <-
                            data3[which(data3[, input$plot2_Y] < facs2$y[2] &
                                          data3[input$plot2_Y] > facs2$y[1]),]
                        } else{
                          data3 <-
                            data[which(data[, input$plot2_X] < facs2$x[2] &
                                         data[, input$plot2_X] > facs2$x[1]),]
                          data3 <-
                            data3[which(data3[, input$plot2_Y] < facs2$y[2] &
                                          data3[input$plot2_Y] > facs2$y[1]),]
                        }
                        output$downloadCells <- downloadHandler(
                          filename = function() {
                            paste("SelectedCells_FC_plots.csv")
                          },
                          content = function(file) {
                            write.csv(as.data.frame(rownames(data3)), file,row.names = F)
                          }
                        )
                        
                        g1 <-
                          ggplot(umaps, aes(rnaumap_1, rnaumap_2)) + geom_point(size = 0.5, col =
                                                                                  "grey89") +
                          geom_point(
                            data = umaps[rownames(data3),],
                            aes(rnaumap_1, rnaumap_2),
                            col = rgb(1, 0, 1, 0.5),
                            size = 0.5
                          ) +
                          theme_bw() + removeGrid() + ggtitle("RNA UMAP") +
                          theme(plot.title = element_text(hjust = 0.5, size = 20))
                        g2 <-
                          ggplot(umaps, aes(adtumap_1, adtumap_2)) + geom_point(size = 0.5, col =
                                                                                  "grey89") +
                          geom_point(
                            data = umaps[rownames(data3),],
                            aes(adtumap_1, adtumap_2),
                            col = rgb(1, 0, 1, 0.5),
                            size = 0.5
                          ) +
                          theme_bw() + removeGrid() + ggtitle("Protein UMAP") +
                          theme(plot.title = element_text(hjust = 0.5, size = 20))
                        g1 + g2
                        
                      } else{
                        if (input$UMAP_gate == "A+B+C") {
                          if (!is.null(facs1$x) & !is.null(facs2$x) & !is.null(facs3$x)) {
                            data2 <-
                              data[which(data[, input$plot1_X] < facs1$x[2] &
                                           data[, input$plot1_X] > facs1$x[1]),]
                            data2 <-
                              data2[which(data2[, input$plot1_Y] < facs1$y[2] &
                                            data2[input$plot1_Y] > facs1$y[1]),]
                            
                            data3 <-
                              data2[which(data2[, input$plot2_X] < facs2$x[2] &
                                            data2[, input$plot2_X] > facs2$x[1]),]
                            data3 <-
                              data3[which(data3[, input$plot2_Y] < facs2$y[2] &
                                            data3[input$plot2_Y] > facs2$y[1]),]
                            
                            data4 <-
                              data3[which(data3[, input$plot3_X] < facs3$x[2] &
                                            data3[, input$plot3_X] > facs3$x[1]),]
                            data4 <-
                              data4[which(data4[, input$plot3_Y] < facs3$y[2] &
                                            data4[input$plot3_Y] > facs3$y[1]),]
                            
                            
                          } else{
                            if (is.null(facs1$x) & !is.null(facs2$x) & is.null(facs3$x)) {
                              data4 <-
                                data[which(data[, input$plot2_X] < facs2$x[2] &
                                             data[, input$plot2_X] > facs2$x[1]),]
                              data4 <-
                                data4[which(data4[, input$plot2_Y] < facs2$y[2] &
                                              data4[input$plot2_Y] > facs2$y[1]),]
                              
                            } else{
                              if (!is.null(facs1$x) & is.null(facs2$x) & is.null(facs3$x)) {
                                data4 <-
                                  data[which(data[, input$plot1_X] < facs1$x[2] &
                                               data[, input$plot1_X] > facs1$x[1]),]
                                data4 <-
                                  data4[which(data4[, input$plot1_Y] < facs1$y[2] &
                                                data4[input$plot1_Y] > facs1$y[1]),]
                              } else{
                                if (is.null(facs1$x) & is.null(facs2$x) & !is.null(facs3$x)) {
                                  data4 <-
                                    data[which(data[, input$plot3_X] < facs3$x[2] &
                                                 data[, input$plot3_X] > facs3$x[1]),]
                                  data4 <-
                                    data4[which(data4[, input$plot3_Y] < facs3$y[2] &
                                                  data4[input$plot3_Y] > facs3$y[1]),]
                                }
                              }
                            }
                          }
                          
                          output$downloadCells <- downloadHandler(
                            filename = function() {
                              paste("SelectedCells_FC_plots.csv")
                            },
                            content = function(file) {
                              write.csv(as.data.frame(rownames(data4)), file,row.names = F)
                            }
                          )
                          
                          g1 <-
                            ggplot(umaps, aes(rnaumap_1, rnaumap_2)) + geom_point(size = 0.5, col =
                                                                                    "grey89") +
                            geom_point(
                              data = umaps[rownames(data4),],
                              aes(rnaumap_1, rnaumap_2),
                              col = rgb(0, 1, 0.5,1),
                              size = 0.5
                            ) +
                            theme_bw() + removeGrid() + ggtitle("RNA UMAP") +
                            theme(plot.title = element_text(hjust = 0.5, size = 20))
                          g2 <-
                            ggplot(umaps, aes(adtumap_1, adtumap_2)) + geom_point(size = 0.5, col =
                                                                                    "grey89") +
                            geom_point(
                              data = umaps[rownames(data4),],
                              aes(adtumap_1, adtumap_2),
                              col = rgb(0, 1, 0.5, 1),
                              size = 0.5
                            ) +
                            theme_bw() + removeGrid() + ggtitle("Protein UMAP") +
                            theme(plot.title = element_text(hjust = 0.5, size = 20))
                          g1 + g2
                          
                        }
                        
                        
                      }
                    }
                  }
                }
              }
            }
          }
        }
        
        
      }
    }
  })
  
  
  
  
  autoInvalidate <- reactiveTimer(2000)
  
  output$RNA_plot <- renderPlot({
    if (input$cluster == 'rna/protein expression'){
      
      
      data_list <- datasetInput_d()
      
      sc <- data_list[[8]]
      
      adtumap1 <- sc@reductions$umap_adt@cell.embeddings[, 1]
      adtumap2 <- sc@reductions$umap_adt@cell.embeddings[, 2]
      rnaumap1 <- sc@reductions$umap_rna@cell.embeddings[, 1]
      rnaumap2 <- sc@reductions$umap_rna@cell.embeddings[, 2]
      DefaultAssay(sc)<- 'RNA'
      gene <- input$gene
      data_list <- datasetInput_d()
      rnaumap_1 = rnaumap1
      rnaumap_2 = rnaumap2
      DefaultAssay(sc)<- 'RNA'
      
      feature<- FeaturePlot(sc,features = c(gene))+xlim(min(rnaumap_1), max(rnaumap_1)) + ylim(min(rnaumap_2), max(rnaumap_2))
      feature$labels$x <- 'UMAP_1'
      feature$labels$y <- 'UMAP_2'
      feature$data$UMAP_1 <- rnaumap_1
      feature$data$UMAP_2 <- rnaumap_2
      feature$data$umapadt_1 <- rnaumap_1
      feature$data$umapadt_2 <- rnaumap_2
      feature
      
    }
    
    
    else {
      
      
      
      if (input$cluster == 'list of cells'){
        data_list <- datasetInput_d()
        data = data_list[[1]]
        umaps = data.frame(
          adtumap_1 = data_list[[2]],
          adtumap_2 = data_list[[3]],
          rnaumap_1 = data_list[[4]],
          rnaumap_2 = data_list[[5]]
        )
        metadata = data_list[[6]]
        sc <- data_list[[8]]
        
        inFile <- input$file1
        
        if (is.null(inFile))
          return(NULL)
        
        cells<- read.csv(inFile$datapath)
        rownames(cells)<- cells[,1]
        sc_sub<- subset(sc,cells=rownames(cells))
        adtumap_1 <- sc_sub@reductions$umap_adt@cell.embeddings[, 1]
        adtumap_2 <- sc_sub@reductions$umap_adt@cell.embeddings[, 2]
        rnaumap_1 <- sc_sub@reductions$umap_rna@cell.embeddings[, 1]
        rnaumap_2 <- sc_sub@reductions$umap_rna@cell.embeddings[, 2]
        rnaumap_2 <- sc_sub@reductions$umap_rna@cell.embeddings[, 2]
        umaps2 = data.frame(
          adtumap_1=adtumap_1,
          adtumap_2=adtumap_2,
          rnaumap_1=rnaumap_1,
          rnaumap_2=rnaumap_2
        )
        
        umaps %>% ggplot(aes(rnaumap_1, rnaumap_2)) + geom_point(col = "grey89") +
          ggtitle("Highlighted select cells") +
          geom_point(
            data = umaps2,
            aes(rnaumap_1, rnaumap_2),
            col = "red",
            size = 0.5
          ) +
          theme_bw() + removeGrid()
        
        
      }
      else{
        
        data_list <- datasetInput_d()
        
        data = data_list[[1]]
        umaps = data.frame(
          adtumap_1 = data_list[[2]],
          adtumap_2 = data_list[[3]],
          rnaumap_1 = data_list[[4]],
          rnaumap_2 = data_list[[5]]
        )
        metadata = data_list[[6]]
        
        
        umaps2 <-
          umaps[which(
            umaps$adtumap_1 < ADTrange$x[2]  & umaps$adtumap_1 > ADTrange$x[1] &
              umaps$adtumap_2 < ADTrange$y[2] &
              umaps$adtumap_2 > ADTrange$y[1]
          ),]
        
        
        
        cluster <- input$cluster
        umaps$cluster <- metadata[, cluster]
        
        if (!is.null(input$ADT_brush)) {
          umaps %>% ggplot(aes(rnaumap_1, rnaumap_2)) + geom_point(col = "grey89") +
            ggtitle("RNA UMAP (select cells to show top differential heatmaps)") +
            geom_point(
              data = umaps2,
              aes(rnaumap_1, rnaumap_2),
              col = "red",
              size = 0.5
            ) +
            labs(colour = cluster) +
            theme_bw() + removeGrid()
        } else{
          umaps %>% ggplot(aes(rnaumap_1, rnaumap_2, col = cluster)) + geom_point(size =
                                                                                    0.5) +
            ggtitle("RNA UMAP (select cells to show top differential heatmaps)") +
            geom_point(
              data = umaps2,
              aes(rnaumap_1, rnaumap_2),
              col = "red",
              size = 0.5
            ) +
            labs(colour = cluster) +
            theme_bw() + removeGrid()

        }
      }
    }
  })
  
  output$ADT_plot <- renderPlot({
    if (input$selected_igt_id == 'Integrated_Data'){
    }
    else{
    if (input$cluster == 'rna/protein expression'){
      
      data_list <- datasetInput_d()
      
      sc <- data_list[[8]]
      
      adtumap1 <- sc@reductions$umap_adt@cell.embeddings[, 1]
      adtumap2 <- sc@reductions$umap_adt@cell.embeddings[, 2]
      rnaumap1 <- sc@reductions$umap_rna@cell.embeddings[, 1]
      rnaumap2 <- sc@reductions$umap_rna@cell.embeddings[, 2]
      DefaultAssay(sc)<- 'ADT'
      
      protein <- input$protein
      data_list <- datasetInput_d()
      adtumap_1 = adtumap1
      adtumap_2 = adtumap2
      DefaultAssay(sc)<- 'ADT'
      feature2<- FeaturePlot(sc,features = c(protein))+xlim(min(adtumap_1), max(adtumap_1)) + ylim(min(adtumap2), max(adtumap2))
      feature2$labels$x <- 'UMAP_1'
      feature2$labels$y <- 'UMAP_2'
      feature2$data$UMAP_1 <- adtumap_1
      feature2$data$UMAP_2 <- adtumap2
      feature2$data$umapadt_1 <- adtumap_1
      feature2$data$umapadt_2 <- adtumap2
      feature2
      
    }
    else {
      
      
      
      if (input$cluster == 'list of cells'){
        data_list <- datasetInput_d()
        data = data_list[[1]]
        umaps = data.frame(
          adtumap_1 = data_list[[2]],
          adtumap_2 = data_list[[3]],
          rnaumap_1 = data_list[[4]],
          rnaumap_2 = data_list[[5]]
        )
        metadata = data_list[[6]]
        sc <- data_list[[8]]
        
        inFile <- input$file1
        
        if (is.null(inFile))
          return(NULL)
        
        cells<- read.csv(inFile$datapath)
        rownames(cells)<- cells[,1]
        sc_sub<- subset(sc,cells=rownames(cells))
        adtumap_1 <- sc_sub@reductions$umap_adt@cell.embeddings[, 1]
        adtumap_2 <- sc_sub@reductions$umap_adt@cell.embeddings[, 2]
        rnaumap_1 <- sc_sub@reductions$umap_rna@cell.embeddings[, 1]
        rnaumap_2 <- sc_sub@reductions$umap_rna@cell.embeddings[, 2]
        rnaumap_2 <- sc_sub@reductions$umap_rna@cell.embeddings[, 2]
        umaps2 = data.frame(
          adtumap_1=adtumap_1,
          adtumap_2=adtumap_2,
          rnaumap_1=rnaumap_1,
          rnaumap_2=rnaumap_2
        )
        
        umaps %>% ggplot(aes(adtumap_1, adtumap_2)) + geom_point(col = "grey89") +
          ggtitle("Highlighted select cells") +
          geom_point(
            data = umaps2,
            aes(adtumap_1, adtumap_2),
            col = "red",
            size = 0.5
          ) +
          theme_bw() + removeGrid()
        
        
      }
      else{
        data_list <- datasetInput_d() 
        data = data_list[[1]]
        umaps = data.frame(
          adtumap_1 = data_list[[2]],
          adtumap_2 = data_list[[3]],
          rnaumap_1 = data_list[[4]],
          rnaumap_2 = data_list[[5]]
        )
        metadata = data_list[[6]]
        
        
        umaps2 <-
          umaps[which(
            umaps$rnaumap_1 < ranges4$x[2]  & umaps$rnaumap_1 > ranges4$x[1] &
              umaps$rnaumap_2 < ranges4$y[2] &
              umaps$rnaumap_2 > ranges4$y[1]
          ),]
        
        cluster <- input$cluster
        umaps$cluster <- metadata[, cluster]
        
        if (!is.null(input$plot4_brush)) {
          umaps %>% ggplot(aes(adtumap_1, adtumap_2)) + geom_point(col = "grey89") +
            ggtitle("Protein UMAP (select cells to show top differential heatmaps)") +
            geom_point(
              data = umaps2,
              aes(adtumap_1, adtumap_2),
              col = "red",
              size = 0.5
            ) +
            labs(colour = cluster) +
            theme_bw() + removeGrid()
        } else{
          umaps %>% ggplot(aes(adtumap_1, adtumap_2, col = cluster)) + geom_point(size =
                                                                                    0.5) +
            ggtitle("Protein UMAP (select cells to show top differential heatmaps)") +
            geom_point(
              data = umaps2,
              aes(adtumap_1, adtumap_2),
              col = "red",
              size = 0.5
            ) +
            labs(colour = cluster) +
            theme_bw() + removeGrid()
        }
      }
    }
    }
  })
  
  output$heatmap <- renderPlot({
    if (is.null(input$plot4_brush) & is.null(input$ADT_brush)) {
      output$downloadCells <- downloadHandler(
        filename = function() {
          paste("Cell_Selection_Warning.csv")
        },
        content = function(file) {
          write.csv(as.data.frame(c('No Cells Selected')), file,row.names = F)
        }
      )
      
    }
    #  if (input$cluster == 'rna/protein expression'){
    #   
    #}
    else {
      if (input$cluster == 'list of cells'){
        
      }
      else{
        data_list <- datasetInput_d() 
        data = data_list[[1]]
        umaps = data.frame(
          adtumap_1 = data_list[[2]],
          adtumap_2 = data_list[[3]],
          rnaumap_1 = data_list[[4]],
          rnaumap_2 = data_list[[5]]
        )
        metadata = data_list[[6]]
        
        rna_data = data_list[[7]]
        
        cluster <- input$cluster
        umaps$cluster <- metadata[, cluster]
        
        if (!is.null(input$plot4_brush) & is.null(input$ADT_brush)) {
          umaps2 <-
            umaps[which(
              umaps$rnaumap_1 < ranges4$x[2]  & umaps$rnaumap_1 > ranges4$x[1] &
                umaps$rnaumap_2 < ranges4$y[2] &
                umaps$rnaumap_2 > ranges4$y[1]
            ),]
          data2 <-
            data[which(
              umaps$rnaumap_1 < ranges4$x[2]  & umaps$rnaumap_1 > ranges4$x[1] &
                umaps$rnaumap_2 < ranges4$y[2] &
                umaps$rnaumap_2 > ranges4$y[1]
            ),]
          
          
          output$downloadCells <- downloadHandler(
            filename = function() {
              paste("SelectedCells_RNA.csv")
            },
            content = function(file) {
              write.csv(as.data.frame(colnames(data2)), file,row.names = F)
            }
          )
          
          
          
          
        } else{
          if (!is.null(input$ADT_brush) & is.null(input$plot4_brush)) {
            umaps2 <-
              umaps[which(
                umaps$adtumap_1 < ADTrange$x[2]  & umaps$adtumap_1 > ADTrange$x[1] &
                  umaps$adtumap_2 < ADTrange$y[2] &
                  umaps$adtumap_2 > ADTrange$y[1]
              ),]
            data2 <-
              data[which(
                umaps$adtumap_1 < ADTrange$x[2]  & umaps$adtumap_1 > ADTrange$x[1] &
                  umaps$adtumap_2 < ADTrange$y[2] &
                  umaps$adtumap_2 > ADTrange$y[1]
              ),]
            
            output$downloadCells <- downloadHandler(
              filename = function() {
                paste("SelectedCells_Protein.csv")
              },
              content = function(file) {
                write.csv(as.data.frame(colnames(data2)), file,row.names = F)
              }
            )
            
            
          }
        }
        
        if (!is.null(input$plot4_brush) |
            !is.null(input$ADT_brush)) {
          FC_ADT <-
            apply(data2, 2, mean) / apply(data[-which(rownames(data) %in% rownames(data2)),], 2, mean)
          FC_ADT <- names(FC_ADT[order(FC_ADT, decreasing = T)])[1:15]
          
          print(rna_data)
          FC_RNA <-
            rowMeans(rna_data[, rownames(data2)]) / rowMeans(rna_data[,-which(colnames(rna_data) %in%
                                                                                rownames(data2))])
          FC_RNA <-
            FC_RNA[names(which((
              rowSums(rna_data > 0) / dim(data2)[2] > 0.05 * dim(data2)[2]
            ) == T))]
          FC_RNA <- names(FC_RNA[order(FC_RNA, decreasing = T)])[1:15]
          FC_RNA <- rowMeans(rna_data[FC_RNA,])
          FC_RNA <- names(FC_RNA)
          
          if (dim(data2)[1] > 2000) {
            s <- sample(1:dim(data2)[1], 2000, replace = F)
            data2 <- data2[s,]
          }
          rnadata <-
            as.data.frame(as.matrix(rna_data[rownames(rna_data) %in% FC_RNA,
                                             rownames(data2)]))
          
          
          column_order = names(colSums(data2[, FC_ADT])[order(colSums(data2[, FC_ADT]), decreasing =
                                                                T)])
          data2 <- t(data2[, column_order])
          h <-
            Heatmap(
              as.matrix(data2),
              name = "Protein(log1p)",
              cluster_rows = F,
              cluster_columns = T,
              show_row_names = T,
              row_names_gp = gpar(fontsize = 15),
              show_heatmap_legend= T,
              # heatmap_legend_param = list(title = ""),
              show_column_dend = F,
              show_column_names = F,
              column_title = "Protein Heatmap (top protein markers differentiating selected cells from rest of cells)",
              col = c('black','red','yellow')
              
            )
          
          ht <- draw(h)
          h2 <-
            Heatmap(
              as.matrix(rnadata[, column_order(ht)]),
              name = "RNA(log1p)",
              cluster_rows = T,
              cluster_columns = F,
              row_names_gp = gpar(fontsize = 15),
              show_heatmap_legend= T,
              #heatmap_legend_param = list(title = ""),
              show_column_names = F,
              column_title = "RNA Heatmap (top gene markers differentiating selected cells from rest of cells)",
              col = c('black','red','yellow')
            )
          
          draw(h2 + h, auto_adjust = F)
        }
      }
    }
  })
  
  
  output$RNA_plot_GC_1 <- renderPlot({
    
    if (input$datatype == 'RNA Plot'){
      
      
      if (input$cluster2 == 'rna/protein expression'){
        
        
        data_list <- datasetInput_d()
        
        sc <- data_list[[8]]
        
        adtumap1 <- sc@reductions$umap_rna@cell.embeddings[, 1]
        adtumap2 <- sc@reductions$umap_rna@cell.embeddings[, 2]
        rnaumap1 <- sc@reductions$umap_rna@cell.embeddings[, 1]
        rnaumap2 <- sc@reductions$umap_rna@cell.embeddings[, 2]
        
        gene <- input$gene
        data_list <- datasetInput_d()
        rnaumap_1 = rnaumap1
        rnaumap_2 = rnaumap2
        feature<- FeaturePlot(sc,features = c(gene))+xlim(min(rnaumap_1), max(rnaumap_1)) + ylim(min(rnaumap_2), max(rnaumap_2))
        feature$labels$x <- 'UMAP_1'
        feature$labels$y <- 'UMAP_2'
        feature$data$UMAP_1 <- rnaumap_1
        feature$data$UMAP_2 <- rnaumap_2
        feature$data$umapadt_1 <- rnaumap_1
        feature$data$umapadt_2 <- rnaumap_2
        feature
        
      }
      
      
      else{

        data_list <- datasetInput_d()
        data = data_list[[11]]
        umaps = data.frame(
          adtumap_1 = data_list[[14]],
          adtumap_2 = data_list[[15]],
          rnaumap_1 = data_list[[14]],
          rnaumap_2 = data_list[[15]]
        )
        metadata = data_list[[16]]
        
        
        
        
        
        

        
        umaps2 <-
          umaps[which(
            umaps$adtumap_1 < ranges5$x[2]  & umaps$adtumap_1 > ranges5$x[1] &
              umaps$adtumap_2 < ranges5$y[2] &
              umaps$adtumap_2 > ranges5$y[1]
          ),]
        
        cluster2 <- input$cluster2
        umaps$cluster2 <- metadata[, cluster2]
        

          
          umaps %>% ggplot(aes(adtumap_1, adtumap_2, col = cluster2)) + geom_point(size =
                                                                                     0.5) +
            ggtitle("A") +
            theme(plot.title = element_text(size = 40, face = "bold"),
                  axis.title.x=element_blank(),
                  axis.title.y=element_blank())+
            geom_point(
              data = umaps2,
              aes(adtumap_1, adtumap_2),
              col = "green",
              size = 0.5
            ) +
            labs(colour = cluster2) +
            removeGrid()
        

        
      }
    }
    
    
    else{
      
      if (input$cluster2 == 'rna/protein expression'){
        
        
        data_list <- datasetInput_d()
        
        sc <- data_list[[8]]
        
        adtumap1 <- sc@reductions$umap_adt@cell.embeddings[, 1]
        adtumap2 <- sc@reductions$umap_adt@cell.embeddings[, 2]
        rnaumap1 <- sc@reductions$umap_adt@cell.embeddings[, 1]
        rnaumap2 <- sc@reductions$umap_adt@cell.embeddings[, 2]
        
        gene <- input$gene
        data_list <- datasetInput_d()
        rnaumap_1 = rnaumap1
        rnaumap_2 = rnaumap2
        feature<- FeaturePlot(sc,features = c(gene))+xlim(min(rnaumap_1), max(rnaumap_1)) + ylim(min(rnaumap_2), max(rnaumap_2))
        feature$labels$x <- 'UMAP_1'
        feature$labels$y <- 'UMAP_2'
        feature$data$UMAP_1 <- rnaumap_1
        feature$data$UMAP_2 <- rnaumap_2
        feature$data$umapadt_1 <- rnaumap_1
        feature$data$umapadt_2 <- rnaumap_2
        feature
        
      }
      
      
      else{
        
        data_list <- datasetInput_d()
        
        data = data_list[[11]]
        umaps = data.frame(
          adtumap_1 = data_list[[12]],
          adtumap_2 = data_list[[13]],
          rnaumap_1 = data_list[[12]],
          rnaumap_2 = data_list[[13]]
        )
        metadata = data_list[[16]]
        
        
        
        
        
        
        umaps2 <-
          umaps[which(
            umaps$adtumap_1 < ranges5$x[2]  & umaps$adtumap_1 > ranges5$x[1] &
              umaps$adtumap_2 < ranges5$y[2] &
              umaps$adtumap_2 > ranges5$y[1]
          ),]
        
        cluster2 <- input$cluster2
        umaps$cluster2 <- metadata[, cluster2]
        
        
        
        umaps %>% ggplot(aes(adtumap_1, adtumap_2, col = cluster2)) + geom_point(size =
                                                                                   0.5) +
          ggtitle("A") +
          theme(plot.title = element_text(size = 40, face = "bold"),
                axis.title.x=element_blank(),
                axis.title.y=element_blank())+
          geom_point(
            data = umaps2,
            aes(adtumap_1, adtumap_2),
            col = "green",
            size = 0.5
          ) +
          labs(colour = cluster2) +
          removeGrid()
        
      }
      
      
      
    }
  })
  
  output$RNA_plot_GC_2 <- renderPlot({
    
    if (input$datatype == 'Protein Plot'){
      data = datasetInput_d()
      
      geneList = data[[9]]
      
      
      updateSelectizeInput(
        session,
        inputId = "gene",
        label = "Select Gene",
        selected = "Cd8a",
        choices = sort(unlist(geneList, recursive = TRUE, use.names = TRUE)),
        server= TRUE
      )
      
      
      
      if (input$cluster2 == 'rna/protein expression'){
        
        data_list <- datasetInput_d()
        
        sc <- data_list[[8]]
        
        adtumap1 <- sc@reductions$umap_adt@cell.embeddings[, 1]
        adtumap2 <- sc@reductions$umap_adt@cell.embeddings[, 2]
        rnaumap1 <- sc@reductions$umap_adt@cell.embeddings[, 1]
        rnaumap2 <- sc@reductions$umap_adt@cell.embeddings[, 2]
        
        protein <- input$protein
        data_list <- datasetInput_d()
        adtumap_1 = adtumap1
        adtumap_2 = adtumap2
        feature2<- FeaturePlot(sc,features = c(protein))+xlim(min(adtumap_1), max(adtumap_1)) + ylim(min(adtumap2), max(adtumap2))
        feature2$labels$x <- 'UMAP_1'
        feature2$labels$y <- 'UMAP_2'
        feature2$data$UMAP_1 <- adtumap_1
        feature2$data$UMAP_2 <- adtumap2
        feature2$data$umapadt_1 <- adtumap_1
        feature2$data$umapadt_2 <- adtumap2
        feature2
        
      }
      else {

        data_list <- datasetInput_d()
        
        data_list <- datasetInput_d() 
        data = data_list[[18]]
        umaps = data.frame(
          adtumap_1 = data_list[[19]],
          adtumap_2 = data_list[[20]],
          rnaumap_1 = data_list[[19]],
          rnaumap_2 = data_list[[20]]
        )
        metadata = data_list[[23]]
        
        
        
        
        umaps2 <-
          umaps[which(
            umaps$adtumap_1 < ranges6$x[2]  & umaps$adtumap_1 > ranges6$x[1] &
              umaps$adtumap_2 < ranges6$y[2] &
              umaps$adtumap_2 > ranges6$y[1]
          ),]
        
        cluster2 <- input$cluster2
        umaps$cluster2 <- metadata[, cluster2]
        
        
        
        umaps %>% ggplot(aes(adtumap_1, adtumap_2, col = cluster2)) + geom_point(size =
                                                                                   0.5) +
          ggtitle("B") +
          theme(plot.title = element_text(size = 40, face = "bold"),
                axis.title.x=element_blank(),
                axis.title.y=element_blank())+
          geom_point(
            data = umaps2,
            aes(adtumap_1, adtumap_2),
            col = "purple",
            size = 0.5
          ) +
          labs(colour = cluster2) +
          removeGrid()
        
      }
    }
    else{
      data = datasetInput_d()
      
      geneList = data[[9]]
      
      

      
      
      
      if (input$cluster2 == 'rna/protein expression'){
        
        data_list <- datasetInput_d()

        sc <- data_list[[18]]
        
        adtumap1 <- sc@reductions$umap_rna@cell.embeddings[, 21]
        adtumap2 <- sc@reductions$umap_rna@cell.embeddings[, 22]
        rnaumap1 <- sc@reductions$umap_rna@cell.embeddings[, 21]
        rnaumap2 <- sc@reductions$umap_rna@cell.embeddings[, 22]
        
        protein <- input$protein
        data_list <- datasetInput_d()
        adtumap_1 = adtumap1
        adtumap_2 = adtumap2
        feature2<- FeaturePlot(sc,features = c(protein))+xlim(min(adtumap_1), max(adtumap_1)) + ylim(min(adtumap2), max(adtumap2))
        feature2$labels$x <- 'UMAP_1'
        feature2$labels$y <- 'UMAP_2'
        feature2$data$UMAP_1 <- adtumap_1
        feature2$data$UMAP_2 <- adtumap2
        feature2$data$umapadt_1 <- adtumap_1
        feature2$data$umapadt_2 <- adtumap2
        feature2
        
      }
      else {
        
        
        data_list <- datasetInput_d() 
        data = data_list[[18]]
        umaps = data.frame(
          adtumap_1 = data_list[[21]],
          adtumap_2 = data_list[[22]],
          rnaumap_1 = data_list[[21]],
          rnaumap_2 = data_list[[22]]
        )
        metadata = data_list[[23]]
        
        
        
        umaps2 <-
          umaps[which(
            umaps$adtumap_1 < ranges6$x[2]  & umaps$adtumap_1 > ranges6$x[1] &
              umaps$adtumap_2 < ranges6$y[2] &
              umaps$adtumap_2 > ranges6$y[1]
          ),]
        
        cluster2 <- input$cluster2
        umaps$cluster2 <- metadata[, cluster2]
        
        
   
        umaps %>% ggplot(aes(adtumap_1, adtumap_2, col = cluster2)) + geom_point(size =
                                                                                   0.5) +
          ggtitle("B") +
          theme(plot.title = element_text(size = 40, face = "bold"),
                axis.title.x=element_blank(),
                axis.title.y=element_blank())+
          geom_point(
            data = umaps2,
            aes(adtumap_1, adtumap_2),
            col = "purple",
            size = 0.5
          ) +
          labs(colour = cluster2) +
          removeGrid()
        
        
        
      }
      
      
      
      
      
      
      
    }
  })
  
  
  output$fileUploaded <- reactive({
    val <- !(is.null(input$plot5_brush) | is.null(input$plot6_brush))
    
    
    
  })
  outputOptions(output, 'fileUploaded', suspendWhenHidden=FALSE)
  
  
  
  
  
 
  

  

  
  
  

  output$heatmap_GC1 <- renderPlot({
    
    if (is.null(input$plot5_brush) |
        is.null(input$plot6_brush)) {
    }else{

    if (input$datatype == 'Protein Plot'){

      
      data_list <- datasetInput_d() 
      data = data_list[[11]]
      umaps = data.frame(
        adtumap_1 = data_list[[12]],
        adtumap_2 = data_list[[13]],
        rnaumap_1 = data_list[[12]],
        rnaumap_2 = data_list[[13]]
      )
      metadata = data_list[[16]]
      
      rna_data = data_list[[17]]
      sc <- data_list[[25]]
      
      
      
      data_list_b <- datasetInput_d() 
      data_b = data_list_b[[18]]
      umaps_b = data.frame(
        adtumap_1_b = data_list_b[[19]],
        adtumap_2_b = data_list_b[[20]],
        rnaumap_1_b = data_list_b[[19]],
        rnaumap_2_b = data_list_b[[20]]
      )
      metadata_b = data_list_b[[23]]
      
      rna_data_b = data_list_b[[24]]
      sc_b <- data_list_b[[26]]
      
      cluster2 <- input$cluster2
      umaps$cluster2 <- metadata[, cluster2]
      umaps_b$cluster2 <- metadata_b[, cluster2]
      
      rna_data<- cbind(rna_data,rna_data_b)
     
    } else{
      

      data_list <- datasetInput_d() 
      data = data_list[[11]]
      umaps = data.frame(
        adtumap_1 = data_list[[14]],
        adtumap_2 = data_list[[15]],
        rnaumap_1 = data_list[[14]],
        rnaumap_2 = data_list[[15]]
      )
      metadata = data_list[[16]]
      
      rna_data = data_list[[17]]
      sc <- data_list[[25]]
      
      
      
      data_list_b <- datasetInput_d() 
      data_b = data_list_b[[18]]
      umaps_b = data.frame(
        adtumap_1_b = data_list_b[[21]],
        adtumap_2_b = data_list_b[[22]],
        rnaumap_1_b = data_list_b[[21]],
        rnaumap_2_b = data_list_b[[22]]
      )
      metadata_b = data_list_b[[23]]
      
      rna_data_b = data_list_b[[24]]
      sc_b <- data_list_b[[26]]
        
        cluster2 <- input$cluster2
        umaps$cluster2 <- metadata[, cluster2]
        umaps_b$cluster2 <- metadata_b[, cluster2]
        rna_data<- cbind(rna_data,rna_data_b)
       
    }
      
      
      
      if (!is.null(input$plot5_brush) &&
          !is.null(input$plot6_brush)) {
        
        
        
        umaps2 <-
          umaps[which(
            umaps$adtumap_1 < ranges5$x[2]  & umaps$adtumap_1 > ranges5$x[1] &
              umaps$adtumap_2 < ranges5$y[2] &
              umaps$adtumap_2 > ranges5$y[1]
          ),]
        data2 <-
          data[which(
            umaps$adtumap_1 < ranges5$x[2]  & umaps$adtumap_1 > ranges5$x[1] &
              umaps$adtumap_2 < ranges5$y[2] &
              umaps$adtumap_2 > ranges5$y[1]
          ),]
        
        umaps3 <-
          umaps_b[which(
            umaps_b$adtumap_1_b < ranges6$x[2]  & umaps_b$adtumap_1_b > ranges6$x[1] &
              umaps_b$adtumap_2_b < ranges6$y[2] &
              umaps_b$adtumap_2_b > ranges6$y[1]
          ),]
        data3 <-
          data_b[which(
            umaps_b$adtumap_1_b < ranges6$x[2]  & umaps_b$adtumap_1_b > ranges6$x[1] &
              umaps_b$adtumap_2_b < ranges6$y[2] &
              umaps_b$adtumap_2_b > ranges6$y[1]
          ),]
        
     
        
        
        overlap<-rownames(data2) %in% rownames(data3)
        if (sum(overlap) == 0) {
        
          
          cells1 <- colnames(sc)
          cells2 <- colnames(sc_b)
          
          # Identify shared cells between both objects
          shared_cells <- intersect(x = cells1, y = cells2)

          
          # Create Obj1 minus shared cells
          if (length(shared_cells) < length(cells1)){
          obj1_minus <- subset(sc, cells = shared_cells, invert = TRUE)
          sc <- merge(obj1_minus, sc_b)
        
          }else{
            sc<- sc

          }
          
     
          

    
          sc_selected <- subset(sc, cells=c(rownames(data3),rownames(data2)) )
          sc_selected$selection<- ifelse(colnames(sc_selected) %in% rownames(data2), "One", "Two")
          Idents(sc_selected) <- "selection"
     
          
          
        } else{
          data4<- data3[-which(rownames(data3) %in% rownames(data2)),]
          data2<- data2[-which(rownames(data2) %in% rownames(data3)),]
          data3<- data4
          
          cells1 <- colnames(sc)
          cells2 <- colnames(sc_b)
          
          # Identify shared cells between both objects
          shared_cells <- intersect(x = cells1, y = cells2)
          
          # Create Obj1 minus shared cells
          if (length(shared_cells) < length(cells1)){
            obj1_minus <- subset(sc, cells = shared_cells, invert = TRUE)
            sc <- merge(obj1_minus, sc_b)
          }else{
            sc<- sc
          }
          
          
          sc_selected <- subset(sc, cells=c(rownames(data4),rownames(data2)) )
          sc_selected$selection<- ifelse(colnames(sc_selected) %in% rownames(data2), "One", "Two")
          Idents(sc_selected) <- "selection"
          
          }
          comp_marks <- FindMarkers(sc_selected, ident.1 = "One", ident.2 = "Two", assay= "RNA",min.diff.pct = 0.1, logfc.threshold = 0)
          download_markers_rna<-comp_marks
          comp_marks<- comp_marks[comp_marks[,1] < 0.05 ,]
          comp_marks_down <- comp_marks[comp_marks[,2] < 0 ,]
          comp_marks_up<- comp_marks[comp_marks[,2] >0,]
          
          output$download_DE_genes <- downloadHandler(
            filename = function() {
              paste("Stats_Differentially_Expressed_Genes.csv")
            },
            content = function(file) {
              write.csv(as.data.frame(download_markers_rna), file,row.names = T)
            }
          )
          


          if (nrow(comp_marks_down)<15 && nrow(comp_marks_up)>15 ){
            comp_marks<- rbind(comp_marks_down[1:nrow(comp_marks_down),],comp_marks_up[1:15,])
            split_up <- rep("A>B",15)
            split_down <- rep("A<B",nrow(comp_marks_down))
          }else{
            if (nrow(comp_marks_up)<15 && nrow(comp_marks_down)>15  ){
              comp_marks<- rbind(comp_marks_down[1:15,],comp_marks_up[1:nrow(comp_marks_up),])
              split_down <- rep("A<B",15)
              split_up <- rep("A>B",nrow(comp_marks_up))
            }else{
              if (nrow(comp_marks_up)<15 && nrow(comp_marks_down)<15  ){
                comp_marks<- rbind(comp_marks_down[1:nrow(comp_marks_down),],comp_marks_up[1:nrow(comp_marks_up),])
                split_down <- rep("A<B",nrow(comp_marks_down))
                split_up <- rep("A>B",nrow(comp_marks_up))
                
              }else {
                if (nrow(comp_marks_up)>15 && nrow(comp_marks_down)>15  ){
                  comp_marks<- rbind(comp_marks_down[1:15,],comp_marks_up[1:15,])
                  split_up <- rep("A>B",15)
                  split_down <- rep("A<B",15)
                }
              }
            }
          }

   
          
          
          FC_RNA <- rownames(comp_marks)

          
          if (dim(data2)[1] > 2000) {
            s <- sample(1:dim(data2)[1], 2000, replace = F)
            data2 <- data2[s,]
          }
          rnadata <-
            as.data.frame(as.matrix(rna_data[FC_RNA,
                                             c(rownames(data2),rownames(data3))]))
          rnadata<- rnadata[order(match(rnadata[,1],FC_RNA)),]
          

          split_A <- rep("A",nrow(data2))
          split_B <- rep("B",nrow(data3))
          
      } 


        h2 <-
          Heatmap(
            as.matrix(rnadata),
            name = "RNA(log1p)",
            cluster_row_slices = FALSE,
            row_split = c(split_down,split_up),
            column_split = c(split_A,split_B),
            column_title_gp = gpar(fill = c("green", "purple"), font = 1:3),
            cluster_rows = F,
            cluster_columns = F,
            show_row_names = T,
            row_names_gp = gpar(fontsize = 15),
            show_heatmap_legend= T,
            show_column_names = F,
            col = c('black','red','yellow')
            #column_title = "RNA Heatmap of top differentially expressed genes (A vs B)"
          )
        
        draw(h2,auto_adjust = F)
    }
      })
  
  output$heatmap_GC2 <- renderPlot({ 
    if (is.null(input$plot5_brush) |
        is.null(input$plot6_brush)) {
    }
    else{
    
    if (input$datatype == 'Protein Plot'){
    

    
        
        
        data_list <- datasetInput_d() 
        data = data_list[[11]]
        umaps = data.frame(
          adtumap_1 = data_list[[12]],
          adtumap_2 = data_list[[13]],
          rnaumap_1 = data_list[[12]],
          rnaumap_2 = data_list[[13]]
        )
        metadata = data_list[[16]]
        
        rna_data = data_list[[17]]
        sc <- data_list[[25]]
        
        
        
        data_list_b <- datasetInput_d() 
        data_b = data_list_b[[18]]
        umaps_b = data.frame(
          adtumap_1_b = data_list_b[[19]],
          adtumap_2_b = data_list_b[[20]],
          rnaumap_1_b = data_list_b[[19]],
          rnaumap_2_b = data_list_b[[20]]
        )
        metadata_b = data_list_b[[23]]
        
        rna_data_b = data_list_b[[24]]
        sc_b <- data_list_b[[26]]
        
        cluster2 <- input$cluster2
        umaps$cluster2 <- metadata[, cluster2]
        umaps_b$cluster2 <- metadata_b[, cluster2]
        
        rna_data<- cbind(rna_data,rna_data_b)
   
      } else{
        data_list <- datasetInput_d() 
        data = data_list[[11]]
        umaps = data.frame(
          adtumap_1 = data_list[[14]],
          adtumap_2 = data_list[[15]],
          rnaumap_1 = data_list[[14]],
          rnaumap_2 = data_list[[15]]
        )
        metadata = data_list[[16]]
        
        rna_data = data_list[[17]]
        sc <- data_list[[25]]
        
        
        
        data_list_b <- datasetInput_d() 
        data_b = data_list_b[[18]]
        umaps_b = data.frame(
          adtumap_1_b = data_list_b[[21]],
          adtumap_2_b = data_list_b[[22]],
          rnaumap_1_b = data_list_b[[21]],
          rnaumap_2_b = data_list_b[[22]]
        )
        metadata_b = data_list_b[[23]]
        
        rna_data_b = data_list_b[[24]]
        sc_b <- data_list_b[[26]]
    
    cluster2 <- input$cluster2
    umaps$cluster2 <- metadata[, cluster2]
    umaps_b$cluster2 <- metadata_b[, cluster2]
    rna_data<- cbind(rna_data,rna_data_b)
    
  }
    
    
    
    if (!is.null(input$plot5_brush) &&
        !is.null(input$plot6_brush)) {
      
      
      
      umaps2 <-
        umaps[which(
          umaps$adtumap_1 < ranges5$x[2]  & umaps$adtumap_1 > ranges5$x[1] &
            umaps$adtumap_2 < ranges5$y[2] &
            umaps$adtumap_2 > ranges5$y[1]
        ),]
      data2 <-
        data[which(
          umaps$adtumap_1 < ranges5$x[2]  & umaps$adtumap_1 > ranges5$x[1] &
            umaps$adtumap_2 < ranges5$y[2] &
            umaps$adtumap_2 > ranges5$y[1]
        ),]
      
      umaps3 <-
        umaps_b[which(
          umaps_b$adtumap_1_b < ranges6$x[2]  & umaps_b$adtumap_1_b > ranges6$x[1] &
            umaps_b$adtumap_2_b < ranges6$y[2] &
            umaps_b$adtumap_2_b > ranges6$y[1]
        ),]
      data3 <-
        data_b[which(
          umaps_b$adtumap_1_b < ranges6$x[2]  & umaps_b$adtumap_1_b > ranges6$x[1] &
            umaps_b$adtumap_2_b < ranges6$y[2] &
            umaps_b$adtumap_2_b > ranges6$y[1]
        ),]
      
      
      
      
      
      overlap<-rownames(data2) %in% rownames(data3)
      if (sum(overlap) == 0) {
      
        cells1 <- colnames(sc)
        cells2 <- colnames(sc_b)
        
        shared_cells <- intersect(x = cells1, y = cells2)
        
        # Create Obj1 minus shared cells
        if (length(shared_cells) < length(cells1)){
          obj1_minus <- subset(sc, cells = shared_cells, invert = TRUE)
          sc <- merge(obj1_minus, sc_b)
        }else{
          sc<- sc
        }
        
        sc_selected <- subset(sc, cells=c(rownames(data3),rownames(data2)) )
        sc_selected$selection<- ifelse(colnames(sc_selected) %in% rownames(data2), "One", "Two")
        Idents(sc_selected) <- "selection"
        
      } else{
        
        data4<- data3[-which(rownames(data3) %in% rownames(data2)),]
        data2<- data2[-which(rownames(data2) %in% rownames(data3)),]
        data3<- data4
        cells1 <- colnames(sc)
        cells2 <- colnames(sc_b)
        
        shared_cells <- intersect(x = cells1, y = cells2)
        
        # Create Obj1 minus shared cells
        if (length(shared_cells) < length(cells1)){
          obj1_minus <- subset(sc, cells = shared_cells, invert = TRUE)
          sc <- merge(obj1_minus, sc_b)
        }else{
          sc<- sc
        }
        
        sc_selected <- subset(sc, cells=c(rownames(data4),rownames(data2)) )
        sc_selected$selection<- ifelse(colnames(sc_selected) %in% rownames(data2), "One", "Two")
        Idents(sc_selected) <- "selection"
        
      }

      
      comp_marks_ADT <- FindMarkers(sc_selected, ident.1 = "One", ident.2 = "Two", assay= "ADT", logfc.threshold = 0)
      download_markers_protein<-comp_marks_ADT
      comp_marks_ADT<- comp_marks_ADT[comp_marks_ADT[,1] < 0.05 ,]
      comp_marks_ADT_down <- comp_marks_ADT[comp_marks_ADT[,2] < 0 ,]
      comp_marks_ADT_up<- comp_marks_ADT[comp_marks_ADT[,2] >0,]
      
     
      output$download_DE_proteins <- downloadHandler(
        filename = function() {
          paste("Stats_Differentially_Expressed_Proteins.csv")
        },
        content = function(file) {
          write.csv(as.data.frame(download_markers_protein), file,row.names = T)
        }
      )
      
      
      if (nrow(comp_marks_ADT_down)<15 && nrow(comp_marks_ADT_up)>15 ){
        comp_marks_ADT<- rbind(comp_marks_ADT_down[1:nrow(comp_marks_ADT_down),],comp_marks_ADT_up[1:15,])
        split_up <- rep("A>B",15)
        split_down <- rep("A<B",nrow(comp_marks_ADT_down))
      }else{
        if (nrow(comp_marks_ADT_up)<15 && nrow(comp_marks_ADT_down)>15  ){
          comp_marks_ADT<- rbind(comp_marks_ADT_down[1:15,],comp_marks_ADT_up[1:nrow(comp_marks_ADT_up),])
          split_down <- rep("A<B",15)
          split_up <- rep("A>B",nrow(comp_marks_ADT_up))
        }else{
          if (nrow(comp_marks_ADT_up)<15 && nrow(comp_marks_ADT_down)<15  ){
            comp_marks_ADT<- rbind(comp_marks_ADT_down[1:nrow(comp_marks_ADT_down),],comp_marks_ADT_up[1:nrow(comp_marks_ADT_up),])
            split_down <- rep("A<B",nrow(comp_marks_ADT_down))
            split_up <- rep("A>B",nrow(comp_marks_ADT_up))
            
          }else {
            if (nrow(comp_marks_ADT_up)>15 && nrow(comp_marks_ADT_down)>15  ){
              comp_marks_ADT<- rbind(comp_marks_ADT_down[1:15,],comp_marks_ADT_up[1:15,])
              split_up <- rep("A>B",15)
              split_down <- rep("A<B",15)
            }
          }
        }
      }
      
      
      FC_ADT<-  rownames(comp_marks_ADT)
      
      
      if (dim(data2)[1] > 2000) {
        s <- sample(1:dim(data2)[1], 2000, replace = F)
        data2 <- data2[s,]
      }
      
      split_A <- rep("A",nrow(data2))
      split_B <- rep("B",nrow(data3))
      
      
      data2<- t(data2)
      data3<- t(data3)
      dataf1 <-as.data.frame(as.matrix(data2[FC_ADT,]))
      dataf2 <-as.data.frame(as.matrix(data3[FC_ADT,]))
      dataf<- cbind(dataf1,dataf2)
    } 
    h <-
      Heatmap(
        as.matrix(dataf),
        name = "Protein(log1p)",
        cluster_rows = F,
        cluster_columns = F,
        row_split = c(split_down,split_up),
        column_split = c(split_A,split_B),
        column_title_gp = gpar(fill = c("green", "purple"), font = 1:3),
        #column_split = column_split,
        show_row_names = T,
        row_names_gp = gpar(fontsize = 15),
        show_heatmap_legend= T,
        # heatmap_legend_param = list(title = ""),
        show_row_dend = T,
        show_column_names = F,
        col = c('black','red','yellow')
        #column_title = "Protein Heatmap of top differentially expressed proteins (A vs B)"
        
      )
 
   
        draw(h,auto_adjust = F)
    }

  })
  
  
  
  
  
  output$ADT_plot_3 <- renderPlot({
    
    if (input$cluster_int == 'rna/protein expression'){
      
      data_list <- datasetInput_d()
      
      sc <- data_list[[8]]
      
      adtumap1 <- sc@reductions$umap_adt@cell.embeddings[, 1]
      adtumap2 <- sc@reductions$umap_adt@cell.embeddings[, 2]
      rnaumap1 <- sc@reductions$umap_rna@cell.embeddings[, 1]
      rnaumap2 <- sc@reductions$umap_rna@cell.embeddings[, 2]
      
      protein <- input$protein
      data_list <- datasetInput_d()
      adtumap_1 = adtumap1
      adtumap_2 = adtumap2
      feature2<- FeaturePlot(sc,features = c(protein))+xlim(min(adtumap_1), max(adtumap_1)) + ylim(min(adtumap2), max(adtumap2))
      feature2$labels$x <- 'UMAP_1'
      feature2$labels$y <- 'UMAP_2'
      feature2$data$UMAP_1 <- adtumap_1
      feature2$data$UMAP_2 <- adtumap2
      feature2$data$umapadt_1 <- adtumap_1
      feature2$data$umapadt_2 <- adtumap2
      feature2
      
    }
    else {
    if (input$cluster_int == 'none'){
      data_list <- datasetInput_d() 
      data = data_list[[1]]
      umaps = data.frame(
        adtumap_1 = data_list[[2]],
        adtumap_2 = data_list[[3]],
        rnaumap_1 = data_list[[4]],
        rnaumap_2 = data_list[[5]]
      )
      metadata = data_list[[6]]
      
      if (input$datatype2 == 'RNA Plot'){
      
      umaps %>% ggplot(aes(rnaumap_1, rnaumap_2, col = 'grey94')) + geom_point(size =
                                                                                0.5,col='grey86') +
        ggtitle("IGT UMAP (select clusters from sidebar to color on Integrated UMAP)") +
        
        theme_bw() + removeGrid()+ theme(legend.position="none")
      }
      else{
        if (input$datatype2 == 'Protein Plot'){
        umaps %>% ggplot(aes(adtumap_1, adtumap_2, col = 'grey94')) + geom_point(size =
                                                                                   0.5,col='grey86') +
          ggtitle("IGT UMAP (select clusters from sidebar to color on Integrated UMAP)") +
          
          theme_bw() + removeGrid()+ theme(legend.position="none")
        
        }
      }
      
    }
    else {
      
      
      
      data_list <- datasetInput_d() 
      data = data_list[[1]]
      umaps = data.frame(
        adtumap_1 = data_list[[2]],
        adtumap_2 = data_list[[3]],
        rnaumap_1 = data_list[[4]],
        rnaumap_2 = data_list[[5]]
      )
      metadata = data_list[[6]]
      
      
      
      cluster <- input$cluster_int
      umaps$cluster <- metadata[, cluster]
      
      
      if (input$datatype2 == 'RNA Plot'){
        
        
        umaps_int2 <- umaps %>% group_by(cluster) %>% select(rnaumap_1, rnaumap_2) %>% summarize_all(mean)
        #ggplot(data, aes(tSNE_1,tSNE_2, col = factor(seurat_clusters, labels= c("..",".."))),label= TRUE)+ geom_point(size=0.3,alpha=0.3) +   labs(color= "Clusters") +  ggrepel::geom_label_repel(data2,aes(label = seurat_clusters))
        
 
       # umaps %>% ggplot(aes(rnaumap_1, rnaumap_2,col = factor(cluster, labels= c(""))),label= TRUE) + geom_point(size =
        umaps %>% ggplot(aes(rnaumap_1, rnaumap_2,col = cluster),label= TRUE) + geom_point(size =
        
                                                                                  0.5) +
          ggtitle("IGT UMAP (select clusters from sidebar to color on Integrated UMAP") +
          
          theme_bw() + removeGrid()+ ggrepel::geom_label_repel(data = umaps_int2,aes(label = cluster))
        
      }
      else{
        if (input$datatype2 == 'Protein Plot'){
        
        
        umaps_int2 <- umaps %>% group_by(cluster) %>% select(adtumap_1, adtumap_2) %>% summarize_all(mean)
        #ggplot(data, aes(tSNE_1,tSNE_2, col = factor(seurat_clusters, labels= c("..",".."))),label= TRUE)+ geom_point(size=0.3,alpha=0.3) +   labs(color= "Clusters") +  ggrepel::geom_label_repel(data2,aes(label = seurat_clusters))
        
        
        # umaps %>% ggplot(aes(rnaumap_1, rnaumap_2,col = factor(cluster, labels= c(""))),label= TRUE) + geom_point(size =
        umaps %>% ggplot(aes(adtumap_1, adtumap_2,col = cluster),label= TRUE) + geom_point(size =
                                                                                             
                                                                                             0.5) +
          ggtitle("IGT UMAP (select clusters from sidebar to color on Integrated UMAP") +
          
          theme_bw() + removeGrid()+ ggrepel::geom_label_repel(data = umaps_int2,aes(label = cluster))
        }
      }
      
    }
    }
  })
  
  
  
  output$Integrated_plot <- renderPlot({
    
  
    
    data_list <- datasetInput_d() 
    data = data_list[[1]]
    umaps = data.frame(
      adtumap_1 = data_list[[2]],
      adtumap_2 = data_list[[3]],
      rnaumap_1 = data_list[[4]],
      rnaumap_2 = data_list[[5]]
    )
    metadata = data_list[[6]]
    
    
    
    
    umaps2 <-
      umaps[which(
        umaps$adtumap_1 < ADTrange_3$x[2]  & umaps$adtumap_1 > ADTrange_3$x[1] &
          umaps$adtumap_2 < ADTrange_3$y[2] &
          umaps$adtumap_2 > ADTrange_3$y[1]
      ),]
    
    
    rownames(umaps2) %in% rownames(integrated_data)
    
    
    
    
    total_umap1 <- integrated_data$umaptotalvi_1
    total_umap2 <- integrated_data$umaptotalvi_2
    #metadata<-integrated_data@meta.data
    #cluster <- input$cluster_int_2
    #umaps$cluster <- metadata[, cluster]
    
    metadata = data_list[[6]]
    

    
    umaps_int = data.frame(total_umap1,total_umap2)
    
    
    
    
    
    if (input$cluster_int == 'none'){
    umaps_int$cluster <- integrated_data$Integrated_Clusters
    
    
    umaps2 <-
      umaps[which(
        umaps$adtumap_1 < ADTrange_3$x[2]  & umaps$adtumap_1 > ADTrange_3$x[1] &
          umaps$adtumap_2 < ADTrange_3$y[2] &
          umaps$adtumap_2 > ADTrange_3$y[1]
      ),]
    
    
    
    integrated_data_highlight <- integrated_data[rownames(integrated_data) %in% rownames(umaps2),]
    
    
    
    #data$seurat_clusters <- factor(data$seurat_clusters,labels=c("..",".."))
    
    umaps_int2 <- umaps_int %>% group_by(cluster) %>% select(total_umap1, total_umap2) %>% summarize_all(mean)
    #ggplot(data, aes(tSNE_1,tSNE_2, col = factor(seurat_clusters, labels= c("..",".."))),label= TRUE)+ geom_point(size=0.3,alpha=0.3) +   labs(color= "Clusters") +  ggrepel::geom_label_repel(data2,aes(label = seurat_clusters))
    
    
    
    
    
    #umaps_int %>% ggplot(aes(total_umap1, total_umap2, col = cluster)) + geom_point(size =0.5)  +ggtitle("Integrated Data UMAP)") + theme_bw() + removeGrid()
    umaps_int %>% ggplot(aes(total_umap1, total_umap2, col = cluster),label= TRUE) + geom_point(size =0.3)  +ggtitle("Integrated Data UMAP") + theme_bw() + removeGrid()+ ggrepel::geom_label_repel(data = umaps_int2,aes(label = cluster))+ theme(legend.position="none")
    
    #umaps_int %>% ggplot(aes(total_umap1, total_umap2, col = cluster)) + geom_point(size =0.5)  +ggtitle("Integrated Data UMAP)") + theme_bw() + removeGrid()
    
    }
    
    else{

      
      cluster <- input$cluster_int
      
      
      rownames(umaps_int)<-rownames(integrated_data)
      
      
      igt<-as.data.frame(metadata[, cluster])
      rownames(igt)<-rownames(metadata)
      colnames(igt)<-'color'
      igt$cells<-rownames(igt)
      
      umaps_int$cells<-rownames(umaps_int)
      umaps_int$cells<-rownames(umaps_int)
      umaps_int_small<-umaps_int
      igt<-igt[rownames(igt) %in% rownames(umaps_int),]
      
      
      
      #color_final = left_join(x=umaps_int_color,y=igt,by="cells")
      #color_final$color[is.na(color_final$color)] <- 'NA'
      #umaps_int$cluster <-color_final$color
      #umaps_int2 <- umaps_int %>% group_by(cluster) %>% select(total_umap1, total_umap2) %>% summarize_all(mean)
      #umaps_int %>% ggplot(aes(total_umap1, total_umap2, col = factor(cluster, labels= c(""))),label= TRUE) + geom_point(size =0.5)  +ggtitle("Integrated Data UMAP)") + theme_bw() + removeGrid()+ ggrepel::geom_label_repel(data = umaps_int2,aes(label = cluster))+ theme(legend.position="none")
      
      color_final = right_join(x=umaps_int,y=igt,by="cells")
      #color_final$color[is.na(color_final$color)] <- 'NA'
      umaps_int$cluster<-'NA'
      
      umaps_int_small<-color_final
      umaps_int_small$cluster<-color_final$color
      umaps_int2 <- umaps_int_small %>% group_by(cluster) %>% select(total_umap1, total_umap2) %>% summarize_all(mean)
     # umaps_int_small %>% ggplot(aes(total_umap1, total_umap2, col = cluster),label= TRUE) + geom_point(size =0.3)  +ggtitle("Integrated Data UMAP)") + theme_bw() + removeGrid()
      #+ ggrepel::geom_label_repel(data = umaps_int2,aes(label = cluster))
      #+ theme(legend.position="none")
      
      umaps_int %>% ggplot(aes(total_umap1, total_umap2, col = cluster)) + geom_point(size =
                                                                                        0.3,col='grey86') +
        geom_point(
          data = umaps_int_small,
          aes(total_umap1, total_umap2, col = cluster),
          size = 0.3) +ggtitle("Integrated Data UMAP") + theme_bw() + removeGrid()+ ggrepel::geom_label_repel(data = umaps_int2,aes(label = cluster))
      #+ theme(legend.position="none")
      
      
      
      
      
    }
  })
  
  
  
  
  
  
  output$sample_table <- renderTable({
    data_list <- datasetInput_d() 
    metadata = data_list[[6]] 
    
    samples <- factor(metadata$sample_name)
    samples
  })
  
  # OBSERVE ----
  
  observe({
    brush_facs1 <- input$facs1_brush
    if (!is.null(brush_facs1)) {
      facs1$x <- c(brush_facs1$xmin, brush_facs1$xmax)
      facs1$y <- c(brush_facs1$ymin, brush_facs1$ymax)
      
    } else {
      facs1$x <- NULL
      facs1$y <- NULL
    }
  })
  
  observe({
    brush <- input$plot4_brush
    if (!is.null(brush)) {
      ranges4$x <- c(brush$xmin, brush$xmax)
      ranges4$y <- c(brush$ymin, brush$ymax)
      
    } else {
      ranges4$x <- NULL
      ranges4$y <- NULL
    }
  })
  
  observe({
    brush_5 <- input$plot5_brush
    if (!is.null(brush_5)) {
      ranges5$x <- c(brush_5$xmin, brush_5$xmax)
      ranges5$y <- c(brush_5$ymin, brush_5$ymax)
      
    } else {
      ranges5$x <- NULL
      ranges5$y <- NULL
    }
  })
  
  observe({
    brush_6 <- input$plot6_brush
    if (!is.null(brush_6)) {
      ranges6$x <- c(brush_6$xmin, brush_6$xmax)
      ranges6$y <- c(brush_6$ymin, brush_6$ymax)
      
    } else {
      ranges6$x <- NULL
      ranges6$y <- NULL
    }
  })
  
  
  observe({
    brush_ADT <- input$ADT_brush
    if (!is.null(brush_ADT)) {
      ADTrange$x <- c(brush_ADT$xmin, brush_ADT$xmax)
      ADTrange$y <- c(brush_ADT$ymin, brush_ADT$ymax)
      
    } else {
      ADTrange$x <- NULL
      ADTrange$y <- NULL
    }
  })
  
  observe({
    brush_facs2 <- input$facs2_brush
    if (!is.null(brush_facs2)) {
      facs2$x <- c(brush_facs2$xmin, brush_facs2$xmax)
      facs2$y <- c(brush_facs2$ymin, brush_facs2$ymax)
      
    } else {
      facs2$x <- NULL
      facs2$y <- NULL
    }
  })
  observe({
    brush_facs3 <- input$facs3_brush
    if (!is.null(brush_facs3)) {
      facs3$x <- c(brush_facs3$xmin, brush_facs3$xmax)
      facs3$y <- c(brush_facs3$ymin, brush_facs3$ymax)
      
    } else {
      facs3$x <- NULL
      facs3$y <- NULL
    }
  })
  observe({
    brush_ADT_3 <- input$ADT_brush_3
    if (!is.null(brush_ADT_3)) {
      ADTrange_3$x <- c(brush_ADT_3$xmin, brush_ADT_3$xmax)
      ADTrange_3$y <- c(brush_ADT_3$ymin, brush_ADT_3$ymax)
      
    } else {
      ADTrange_3$x <- NULL
      ADTrange_3$y <- NULL
    }
  })
  
}

shinyApp(ui, server)
