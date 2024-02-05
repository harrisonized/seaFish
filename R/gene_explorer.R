## Plots UMAP of labeled on the left, gene-of-interest expression on the right
## Useful for exploratory data analysis

wd = dirname(this.path::here())  # wd = '~/github/R/seaFish'
library('shiny')
library('optparse')
import::from(ggplot2, ggtitle)
import::from(Seurat, 'DimPlot', 'FeaturePlot')
import::from(file.path(wd, 'R', 'tools', 'file_io.R'),
    'load_rdata', .character_only=TRUE)

options(shiny.host = '127.0.0.1')
options(shiny.port = 3838)


# args
option_list = list(
    make_option(c("-i", "--input-dir"),
                default='data/ballesteros-2020/output/rdata',
                metavar='data/ballesteros-2020/output/rdata',
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
ui <- shinyUI(fluidPage(

    titlePanel("scRNAseq Data"),

    inputPanel(
        selectInput(
            inputId="dataset",
            label="Dataset",
            choices=c("bm", "lung", "pb", 'spleen')
        ),
        selectInput(
            inputId="gene",
            label="Gene Name",
            choices=""
        ),
        align="center"
    ),

    fluidRow(
        splitLayout(
            cellWidths = c("50%", "50%"),
            plotOutput("clusters_img", height=600),
            plotOutput("gene_fig", height=600)
        ), align="center",
    )
))


# server-side computations
server <- function(input, output, session) {

    # UMAP with cell_type labels
    output$clusters_img <- renderPlot({

        seurat_obj <- load_rdata(
            file.path(wd, opt[['input-dir']], paste0('integrated-', input$dataset, '.RData'))
        )

        num_cells_per_label <- as.data.frame(table(seurat_obj$cell_type))  # value counts
        populations_to_keep <- num_cells_per_label[num_cells_per_label['Freq'] > 50, 'Var1']
        seurat_obj_subset <- seurat_obj[, seurat_obj$cell_type %in% populations_to_keep]

        DimPlot(
            seurat_obj_subset,
            reduction = "umap", 
            group.by = "cell_type",
            split.by = NULL,
            label = TRUE) +
            ggtitle(input$dataset)

    }, height=480, width=600, scaling=1.5)

    # UMAP with gene of interest
    output$gene_fig <- renderPlot({

        seurat_obj <- load_rdata(
            file.path(wd, opt[['input-dir']], paste0('integrated-', input$dataset, '.RData'))
        )

        FeaturePlot(
            seurat_obj,
            reduction = "umap", features = input$gene,
            pt.size = 0.4, min.cutoff = 'q10', order = TRUE, label = FALSE) +
            ggtitle(input$gene)

    }, height=480, width=600, scaling=1.5)

    # Gene dropdown menu
    gene_names <- reactive({
        seurat_obj <- load_rdata(
            file.path(wd, opt[['input-dir']], paste0('integrated-', input$dataset, '.RData'))
        )
        sort(unique(seurat_obj@assays[['RNA']]@data@Dimnames[[1]]))
    })

    observe({
        updateSelectInput(
            session,
            inputId="gene",
            choices = gene_names()
        )
    })
}

# run app
shinyApp(ui = ui, server = server)
