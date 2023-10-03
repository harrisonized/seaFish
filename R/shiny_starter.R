## Plots the integrated dataset
## Includes a dropdown menu to select the dataset and a text box to type in the gene

wd = dirname(this.path::here())  # wd = '~/github/R/harrisonRTools'
library('shiny')
suppressMessages(library('shinydashboard'))
suppressMessages(library('Seurat'))
library('optparse')
source(file.path(wd, 'R', 'functions', 'file_io.R'))  # load_rdata

options(shiny.host = '127.0.0.1')
options(shiny.port = 3838)


# args
option_list = list(
    make_option(c("-i", "--input-dir"),
                default='data/scrnaseq-ballesteros/integrated',
                metavar='data/scrnaseq-ballesteros/integrated',
                type="character",
                help="output of scrnaseq_integrated_analysis.R"),
    
    make_option(c("-f", "--figures-dir"),
                default="figures/scrnaseq-ballesteros/integrated",
                metavar="figures/scrnaseq-ballesteros/integrated", type="character",
                help="set the figure directory"),

    make_option(c("-g", "--gene-of-interest"), default="Dnase1l1",
                metavar="Dnase1l1", type="character",
                help="choose a default gene")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)


# user interface elements and layout
ui <- fluidPage(
    titlePanel("Ballesteros scRNAseq"),
    sidebarLayout(

        sidebarPanel(
            sidebarMenu(
                menuItem("Select the data",
                         tabName = "Select the data",
                         icon = icon("bar-chart")
                ),
                selectInput("dataset",
                            "Dataset",
                            c("bm", "lung", "pbzt", 'spleen')
                ),
                textInput("gene",
                          "Gene Name",
                          value = "Dnase1l1"
                )
            )
        ),

        mainPanel(
            plotOutput(outputId = "fig"),  # gene expression
            plotOutput(outputId = "img")  # clusters
        )
    )
)


# server-side computations
server <- function(input, output, session) {

    # plot gene expression
    output$fig <- renderPlot({
        seurat_obj <- load_rdata(
            file.path(wd, opt[['input-dir']],
            input$dataset, 'integrated_seurat.RData')
        )
        FeaturePlot(
            seurat_obj,
            reduction = "umap", 
            features = input$gene,
            pt.size = 0.4, 
            order = TRUE,
            split.by = "treatment",
            min.cutoff = 'q10',
            label = FALSE
        )
    })

    # plot clusters
    output$img <- renderImage({
        filepath = file.path(
            wd, opt[['figures-dir']], input$dataset,
            paste0('umap-', input$dataset, '-integrated-immgen_labeled.png')
        )
        list(
            src = filepath,
            width = "650px", 
            height = "400px"
        )
    }, deleteFile = FALSE)

}

# run app
shinyApp(ui = ui, server = server)
