insertMinorTicks <- function(xStart, xEnd, stepSize, stepTick)
{
  axisBreaks <- seq(xStart, xEnd, by = stepSize)

  tickLabels <- array("", length(axisBreaks))

  idxTicks <- seq(1, length(tickLabels), by = stepTick)
  tickLabels[idxTicks] <- axisBreaks[idxTicks]

  if (tickLabels[length(tickLabels)] == "")
    warning('Select a different stepSize or stepTick')

  return(tickLabels)
}
