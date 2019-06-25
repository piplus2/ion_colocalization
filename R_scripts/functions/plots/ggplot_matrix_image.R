ggplotMatrixImageContinuous <- function(xMat, RGB = F, theme = 'minimal',
                                        xlabel = 'X', ylabel = 'Y', main = NULL,
                                        addROI = NULL)
{
  ## Require ----
  require(ggplot2)
  require(reshape)
  require(viridis)

  ## Main ----
  df <- melt(xMat)

  gg <- ggplot(df, aes(x = X1, y = X2, fill = value)) +
    geom_raster() +
    coord_fixed() +
    xlab(xlabel) +
    ylab(ylabel)

  if (theme == 'minimal')
    gg <- gg + theme_minimal()
  if (theme == 'void')
    gg <- gg + theme_void()

  gg <- gg + theme(axis.title.y = element_text(face = 'bold', size = 10), #14
                   axis.title.x = element_text(face = 'bold', size = 10))

  if (RGB)
  {
    gg <- gg + scale_fill_identity()
  } else
  {
    gg <- gg + scale_fill_viridis(name = 'intensity', option = 'E')
  }

  if (!is.null(main))
    gg <- gg + ggtitle(main) +
    theme(plot.title = element_text(face = 'bold', hjust = 0.5, size = 14))

  if (!is.null(addROI))
  {
    df_roi <- melt(addROI)
    df_roi <- df_roi[df_roi$value == 1, ]
    gg <- gg + geom_point(data = df_roi, aes(x = X1, y = X2),
                         inherit.aes = F, color = 'white', size = 0.1)
  }

  return(gg)
}


ggplotMatrixImageDiscrete <- function(xMat, title = NULL, showLegend = TRUE)
{
  ## Require ----
  require(ggplot2)
  require(reshape)
  require(RColorBrewer)

  ## Main ----
  df <- melt(xMat)
  df$value <- factor(df$value)

  gg <- ggplot(df, aes(x = X1, y = X2, fill = value)) +
    geom_raster() +
    # scale_fill_brewer(name = 'cluster', type = 'qual', palette = 'Set1') +
    scale_fill_viridis_d(name = 'cluster') +
    coord_fixed() +
    xlab('X') + ylab('Y')

  if (!showLegend)
    gg <- gg + guides(fill = FALSE)

  if (!is.null(title))
    gg <- gg + ggtitle(title)

  gg <- gg + theme_bw() +
    theme(panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(colour = "black", fill = NA, size = 1),
          axis.text = element_blank(),
          axis.title = element_blank(),
          plot.title = element_text(face = 'bold', hjust = 0.5, size = 14))

  return(gg)
}
