# Re-calibrate using LOESS fit on matched reference mass
#
# Author: Paolo Inglese, Imperial College London


recalibratePeaks <- function(listPeaks, referenceMZ, imageShape,
                             tolerancePPM = 10,
                             genPlot = TRUE, saveRecal = TRUE,
                             saveFilename) {
  require(MALDIquant)
  require(MALDIquantForeign)

  numspectra <- length(listPeaks)

  if (genPlot) {
    par(mfrow = c(2, 2))
  }

  # Find the peaks closer than tolerancePPM to the referenceMZ
  matched <- lapply(listPeaks, function(z) {
    d <- sapply(z@mass, function(w) abs(w - referenceMZ) / referenceMZ * 1e6)

    if (sum(d <= tolerancePPM) > 0) {
      matchidx <- which(d <= tolerancePPM)
      if (length(matchidx) > 1) {
        matchidx <- matchidx[which.max(z@intensity[matchidx])]
      }
      matchdist <- z@mass[matchidx] - referenceMZ
    } else {
      matchidx <- NA
      matchdist <- NA
    }
    return(c(matchidx, matchdist))
  })
  matched_mat <- matrix(NA, numspectra, 2,
    dimnames = list(NULL, c("matchedIdx", "distFromRef"))
  )

  for (j in 1:numspectra) {
    matched_mat[j, ] <- matched[[j]]
  }

  rm(matched)

  # Convert to dataframe
  matched_mat <- data.frame(matched_mat)
  matched_pixels_idx <- which(!is.na(matched_mat$matchedIdx))

  # Plot the intensity image of the matched peaks
  matched_intensity <- sapply(1:numspectra, function(z) {
    listPeaks[[z]]@intensity[matched_mat$matchedIdx[z]]
  })

  if (genPlot) {
    image(matrix(log2(matched_intensity + 1), imageShape[1], imageShape[2]))
  }

  print(min(matched_mat$distFromRef, na.rm = TRUE))
  print(max(matched_mat$distFromRef, na.rm = TRUE))

  if (genPlot) {
    plot(1:numspectra, matched_mat$distFromRef,
      ylim = c(
        -max(abs(matched_mat$distFromRef), na.rm = TRUE),
        max(abs(matched_mat$distFromRef), na.rm = TRUE)
      ),
      pch = 20
    )
  }

  # Fit LOESS
  df <- data.frame(
    x = matched_pixels_idx,
    y = matched_mat$distFromRef[matched_pixels_idx]
  )
  fit <- loess(y ~ x, data = df)

  # Update the peaks position
  fitteddist <- predict(fit, seq(1, numspectra))

  if (genPlot) {
    points(1:numspectra, fitteddist, col = "red", pch = 20)
  }

  if (any(is.na(fitteddist))) {
    # If some predictions fail, take the average of the shifts in the interval
    # [-2, -1, x, +1, +2]
    fittedNA <- which(is.na(fitteddist))
    for (i in 1:length(fittedNA))
    {
      avginterval <- seq(fittedNA[i] - 2, fittedNA[i] + 2)
      avginterval <- avginterval[avginterval >= 1 & avginterval <= numspectra]
      fitteddist[fittedNA[i]] <- mean(fitteddist[avginterval], na.rm = T)
    }
  }

  # Shift the m/z values
  for (j in 1:numspectra)
  {
    if (is.na(fitteddist[j])) {
      next()
    }
    listPeaks[[j]]@mass <- listPeaks[[j]]@mass - fitteddist[j]
    listPeaks[[j]] <- createMassPeaks(
      mass = listPeaks[[j]]@mass,
      intensity = listPeaks[[j]]@intensity,
      metaData = listPeaks[[j]]@metaData
    )
  }

  # Fit LOESS after recalibrating
  recalDist <- sapply(1:numspectra, function(z) {
    listPeaks[[z]]@mass[matched_mat$matchedIdx[[z]]] - referenceMZ
  })

  if (genPlot) {
    plot(1:numspectra, recalDist, col = "blue", pch = 20)
  }

  if (saveRecal) {
    cat("saving re-calibrated data...\n")
    ## Convert into spectra
    spectra <- lapply(listPeaks, function(z) {
      createMassSpectrum(mass = z@mass, intensity = z@intensity, metaData = z@metaData)
    })
    exportImzMl(x = spectra, file = saveFilename, processed = TRUE, force = TRUE)

    rm(spectra)
    gc()
  }

  return(listPeaks)
}
