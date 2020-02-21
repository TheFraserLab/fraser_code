library(gplots)

#' Wrapper function to easily make a good looking heatmap
#'
#' @param x numeric  matrix - data use for heatmap, must have column names and row names
#' @param heatmap_gradient colorRampPalette - a color gradient used for coloring values in the heatmap
#' @param dendrogram character - same as in dendrogram of gplots::heatmap2
#' @param Rowv logical - same as in Rowv of gplots::heatmap2
#' @param Colv logical - same as in Colv of gplots::heatmap2
#' @param symbreaks logical - same as in symbreaks of gplots::heatmap2
#' @param trace character - same as in trace of gplots::heatmap2
#' @param density character - same as in density of gplots::heatmap2
#' @param ... passed to gplots::heatmap2
#'
#' @return draws a heatmap
#'
#' @examples
#' mat <- matrix(rnorm(2500),  ncol=50, nrow=50)
#' Heatmap(mat)

Heatmap <- function(x, heatmap_gradient = NULL, dendrogram = "none", Rowv = FALSE, Colv = FALSE, symbreaks = F, trace = "none", density = "none", ...) {
    
    library(gplots)
    
    if(is.null(heatmap_gradient)) 
        heatmap_gradient <- colorRampPalette(c("#FFF5F0", "#3399ff"))(20)
    
    x <- as.matrix(x)
    
    #----
    # Gathering pars for heatmap2
    heatPars <- list(x = x,
                     col = heatmap_gradient, 
                     dendrogram = dendrogram, Rowv = Rowv, Colv = Colv,
                     symbreaks = symbreaks,
                     #scale = "column",
                     trace = trace, density = density,
                     ...
                     )
    
    #----
    
    return(do.call(heatmap.2, heatPars))
    
    
}
