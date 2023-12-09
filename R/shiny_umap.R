## Plots the integrated dataset
## Includes a dropdown menu to select the dataset and a text box to type in the gene

wd = dirname(this.path::here())  # wd = '~/github/R/seaFish'
library('shiny')
library('ggplot2')
suppressMessages(library('shinydashboard'))
suppressMessages(library('Seurat'))
library('optparse')
source(file.path(wd, 'R', 'functions', 'file_io.R'))  # load_rdata

options(shiny.host = '127.0.0.1')
options(shiny.port = 3838)


# args
option_list = list(
    make_option(c("-i", "--input-dir"),
                default='data/ballesteros-2020/output/integrated',
                metavar='data/ballesteros-2020/output/integrated',
                type="character",
                help="output of scrnaseq_integrated_analysis.R"),
    
    make_option(c("-f", "--figures-dir"),
                default="figures/ballesteros-2020/output/integrated",
                metavar="figures/ballesteros-2020/output/integrated",
                type="character",
                help="set the figure directory"),

    make_option(c("-g", "--gene-of-interest"), default="Dnase1l1",
                metavar="Dnase1l1", type="character",
                help="choose a default gene")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)


# user interface elements and layout
ui <- fluidPage(
    titlePanel("scRNAseq Data"),
    inputPanel(
        selectInput("dataset",
                "Dataset",
                c("bm", "lung", "pb", 'spleen')
        ),
        textInput("gene",
                "Gene Name",
                value = "Dnase1l1"
        ),
        align="center"
    ),
    fluidRow(
        splitLayout(
            cellWidths = c("50%", "50%"),
            plotOutput("clusters_img", height=600),
            plotOutput("gene_fig", height=600)
        ),
        align="center",
    )
)


# server-side computations
server <- function(input, output, session) {

    # clusters
    output$clusters_img <- renderImage({
        filepath = file.path(
            wd, opt[['figures-dir']],
            paste0('umap-integrated-', input$dataset, '.png')
        )
        list(
            src = filepath,
            height = "600px", width = "750px"
        )
    }, deleteFile = FALSE)

    # gene of interest
    output$gene_fig <- renderPlot({
        seurat_obj <- load_rdata(
            file.path(wd, opt[['input-dir']],
            paste0(input$dataset, '.RData'))
        )
        FeaturePlot(
            seurat_obj,
            reduction = "umap", 
            features = input$gene,
            pt.size = 0.8, 
            order = TRUE,
            # split.by = "treatment",
            min.cutoff = 'q10',
            label = FALSE
        ) + ggtitle(input$gene)
    }, height=600, width=750, scaling=1.5)

}

# run app
shinyApp(ui = ui, server = server)
