ggplotCombinedSameLegend <- function(listPlots, layout)
{
  ## Require ----
  require(gridExtra)
  require(grid)
  require(ggplot2)
  
  ## Main ----
  numPlots <- length(listPlots)
  stopifnot(numPlots > 1) ## You don't need this function

  g <- ggplotGrob(listPlots[[1]] + theme(legend.position = "bottom"))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  ggGrid <- grid.arrange(
    do.call(arrangeGrob, lapply(listPlots, function(x)
      x + theme(legend.position="none"))),
    legend,
    ncol = 1,
    heights = unit.c(unit(1, "npc") - lheight, lheight))
  
  return(ggGrid)
}

ggplotCombined <- function(listPlots, layout, main, ...)
{
  ## Require ----
  require(ggplot2)
  require(ggpubr)
  
  ggGrid <- ggarrange(plotlist = listPlots, nrow = layout[1], ncol = layout[2]) + theme(...)
  ggGrid <- annotate_figure(ggGrid, top = text_grob(main, face = 'bold', size = 14))
  
  return(ggGrid)
}