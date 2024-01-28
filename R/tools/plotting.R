# library('Matrix')
wd = dirname(dirname(this.path::here()))
import::here(ggplot2, 'ggsave')
import::here(grid, 'grid.newpage', 'grid.draw')
import::here(file.path(wd, 'R', 'tools', 'list_tools.R'),
    'filter_list_for_match', .character_only=TRUE)

## Functions
## save_fig

#' Save Figure
#' 
#' @description Switch case to reduce the number of lines in the main script
#' 
savefig <- function(
    filepath,
    fig=NULL,
    height=800, width=1200, dpi=300, units="px", scaling=0.5,
    makedir=FALSE,
    troubleshooting=FALSE,
    lib='ggplot'  # choose: ggplot, grid
) {
    if (!troubleshooting) {

        # make directory
        dirpath <- dirname(filepath)
        if (makedir && !dir.exists(dirpath)) {
            dir.create(dirpath, recursive=TRUE)
        }

        if (lib=='ggplot') {
            ggsave(
                filepath,
                height=height, width=width, dpi=dpi, units=units, scaling=scaling
            )
        } else if (lib=='grid') {
            png(filepath,
                height=height, width=width, res=dpi, units=units
            )
            grid.newpage()
            grid.draw(fig$gtable)
            dev.off()
        } else {
            warning(paste0("lib='", lib, "' not found"))
        }
    }
}
