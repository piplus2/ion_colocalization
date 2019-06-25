## Collection of functions used for plotting
##

plotModuleGraph <- function(adj,
                            modules,
                            network = NULL,
                            selectedModule,
                            outputDir,
                            onlyTest = FALSE,
                            colorNodes = NULL,
                            plotSignificant = TRUE,
                            showOnlyPositive = TRUE,
                            signifThreshold = 0.7,
                            quantileThreshold = 0.99,
                            baseSize = 2,
                            nodeSize = NULL,
                            title = c('control', 'test'))
{
  .qgraphTemplate <- function(color = NULL, ...)
  {
    stopifnot(require(qgraph))
    
    if (!is.null(color))
    {
      qgraph::qgraph(layout = "spring", diag = FALSE, edge.label.bg = TRUE,
                     aspect = TRUE, edge.label.cex = 1, color = color, title.cex = 3, esize = 5, 
                     maximum = 1, ...)
    } else
    {
      qgraph::qgraph(layout = "spring", diag = FALSE, edge.label.bg = TRUE,
                     aspect = TRUE, edge.label.cex = 1, title.cex = 3, esize = 5, ...)
    }
  }
  
  stopifnot(length(selectedModule) == 1)
  
  if (is.null(network))
    network <- adj
  
  mz <- as.numeric(colnames(adj[[1]]))
  numVars <- length(mz)
  ## Match the ions ----
  modIons <- as.numeric(names(modules)[modules == selectedModule])
  stopifnot(length(modIons) > 0)
  modIdx <- which(round(mz, 4) %in% round(modIons, 4))
  stopifnot(length(modIdx) == length(modIons))
  ## Check whether directed or undirected graph
  directed <- array(F, 2)
  for (set in 1:2)
  {
    if (!isSymmetric(adj[[set]]))
    {
      directed[set] <- T
    }
  }
  ## Edge significance ----
  r1 <- adj[[1]][modIdx, modIdx]
  r2 <- adj[[2]][modIdx, modIdx]
  diag(r1) <- 0
  diag(r2) <- 0
  
  n1 <- network[[1]][modIdx, modIdx]
  n2 <- network[[2]][modIdx, modIdx]
  
  rm(network)
  
  if (showOnlyPositive)
  {
    m1 <- r1 < 0
    m2 <- r2 < 0
    
    r1[m1] <- 0
    r2[m2] <- 0
    
    n1[m1] <- 0
    n2[m2] <- 0
    
    diag(n1) <- 0
    diag(n2) <- 0
    
    rm(m1, m2)
  }
  if (!is.null(colorNodes))
  {
    colorNodes <- colorNodes[modIdx]
  }
  ## Plot the module ----
  outControlFile <- paste0(outputDir, "/module_", selectedModule, "_control.pdf")
  outTestFile <- paste0(outputDir, "/module_", selectedModule, "_test.pdf")
  s <- NULL
  if (plotSignificant)
  {
    ## Plot the significant edges ----
    adjLowerTri <- c(adj[[1]][lower.tri(adj[[1]], diag=FALSE)])
    q <- quantile(adjLowerTri, probs = quantileThreshold)
    cat(sprintf('edge weight quantile (%.2f%%): %.4f\n', quantileThreshold*100, q))
    s <- (n1 / mean(n1)) / (n1 / mean(n1) + n2 / mean(n2) + 1e-5)
    sMask <- (s >= signifThreshold & r1 >= q) * 1
    sMask[sMask == 1] <- "s"
    sMask[sMask == 0] <- NA
    
    rm(n1, n2)
    
    if (is.null(colorNodes))
    {
      colorNodes <- rep('white', nrow(r1))
      rc <- which(sMask == 's', arr.ind = T)
      rc <- unique(rc)
      colorNodes[rc] <- 'yellow'
    }
    
    if (!onlyTest)
    {
      pdf(file = outControlFile, w = 8*3, h = 6*3)
      .qgraphTemplate(input=r1, labels=round(modIons, 2), #edge.labels=sMask,
                      directed=directed[1], color=colorNodes, vsize = exp(nodeSize[[1]]) * baseSize, width = 8*2, height = 6*2,
                      mar = c(1, 1, 1, 1), title = title[1])
      dev.off()
    }
    pdf(file = outTestFile, w = 8*3, h = 6*3)
    .qgraphTemplate(input=r2, labels=round(modIons, 2), #edge.labels=sMask,
                   directed=directed[2], color=colorNodes, vsize = exp(nodeSize[[2]]) * baseSize, width = 8*2, height = 6*2,
                   mar = c(1, 1, 1, 1), title = title[2])
    dev.off()
  } else
  {
    ## Plot without significant edges ----
    if (!onlyTest)
    {
      pdf(file = outControlFile, w = 8*3, h = 6*3)
      .qgraphTemplate(input=r1, labels=round(modIons, 4), directed=directed[1],
                      color=colorNodes, vsize = exp(nodeSize[[1]]) * baseSize, width = 8*2, height = 6*2,
                      mar = c(1, 1, 1, 1), title = title[1])
      dev.off()
    }
    pdf(file = outTestFile, w = 8*3, h = 6*3)
    .qgraphTemplate(input=r2, labels=round(modIons, 4), directed=directed[2], color=colorNodes,
                    color=colorNodes, vsize = exp(nodeSize[[2]]) * baseSize, width = 8*2, height = 6*2,
                    mar = c(1, 1, 1, 1), title = title[2])
    dev.off()
  }
  return(s)
}
