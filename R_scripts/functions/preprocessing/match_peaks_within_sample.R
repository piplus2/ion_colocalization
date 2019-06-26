matchPeaksWithinSample <- function(imzMLPath,
                                   minSignalPixelsFrac = 0.005,
                                   deiso = TRUE,
                                   tolPPM = 2000,
                                   verbose = TRUE) {
  # REQUIRE ----
  require(MALDIquant)
  require(MALDIquantForeign)
  require(magrittr)

  .unlist <- function(x) {
    unlist(x, recursive = FALSE, use.names = FALSE)
  }

  # Read the centroided data
  if (verbose) {
    cat("reading peak list...\n")
  }

  rawDir <- strsplit(imzMLPath, split = "/")[[1]] %>%
    .[-length(.)] %>%
    paste0(., collapse = "/")

  pixelPeaks <- importImzMl(imzMLPath, centroided = T, verbose = F)
  numPixels <- length(pixelPeaks)

  # Determine the spatial dimensions
  sz <- pixelPeaks[[numPixels]]@metaData$imaging$pos

  if (numPixels != prod(sz)) {
    stop("image size and number of pixels do not match.")
  }

  # Match pixelPeaks -----------------------------------------------------------

  numPixels <- length(pixelPeaks)

  # customized binPeaks from MALDIquant
  if (verbose) {
    cat("matching peaks...\n")
  }

  mass <- unname(.unlist(lapply(pixelPeaks, function(x) x@mass)))
  intensities <- .unlist(lapply(pixelPeaks, function(x) x@intensity))
  samples <- c(rep.int(seq_along(pixelPeaks), lengths(pixelPeaks)))
  s <- sort.int(mass, method = "quick", index.return = TRUE)
  mass <- s$x
  intensities <- intensities[s$ix]
  samples <- samples[s$ix]

  # binning
  mass <- .binPeaks_(
    mass = mass, intensities = intensities, samples = samples,
    tolerance = tolPPM / 1e6, grouper = .grouperStrict_
  )
  mass <- round(mass, 6)
  tableMasses <- table(mass)
  # Remove peaks present in less than minSignalPixelsFrac
  tableMasses <- as.numeric(names(tableMasses[tableMasses >=
    minSignalPixelsFrac * numPixels]))
  stopifnot(!any(duplicated(tableMasses)))
  commonMasses <- sort(unique(tableMasses))

  rm(pixelPeaks)
  gc()

  if (deiso) {
    cat("de-isotoping...\n")
    deiso_cmz_list <- deisotope(commonMasses)
    commonMasses <- as.numeric(names(deiso_cmz_list))
    cat(paste0("length m/z vector = ", length(commonMasses), "\n"))
  }

  # Filter also the other vectors
  cat("selecting elements...\n")
  keep_idx <- which(mass %in% commonMasses)
  stopifnot(sort(unique(mass[keep_idx])) == commonMasses)

  mass <- mass[keep_idx]
  intensities <- intensities[keep_idx]
  samples <- samples[keep_idx]

  rm(keep_idx)

  # Generate intensity matrix ----
  if (verbose) {
    cat("Generating feature matrix...\n")
  }

  in_ <- findInterval(mass, commonMasses)

  X <- matrix(NA_real_,
    nrow = numPixels, ncol = length(commonMasses),
    dimnames = list(NULL, commonMasses)
  )
  X[cbind(samples, in_)] <- intensities
  attr(X, "mass") <- commonMasses

  if (verbose) {
    cat("dim X:", dim(X), "\n")
  }

  return(list(X = X, shape = sz, mz = as.numeric(colnames(X))))
}

.binPeaks_ <- function(mass, intensities, samples, tolerance,
                       grouper = .grouperStrict_, ...) {
  n <- length(mass)

  # calculate difference
  d <- diff(mass)

  # grouper function
  grouper <- match.fun(grouper)

  ## stack based implementation taken from
  ## caMassClass 1.9 R/msc.peaks.clust.R written by
  ## Jarek Tuszynski <jaroslaw.w.tuszynski@saic.com>
  ## it is a lot of faster than recursion

  ## store boundaries in a stack
  nBoundaries <- max(20L, floor(3L * log(n)))
  boundary <- list(left = double(nBoundaries), right = double(nBoundaries))

  currentBoundary <- 1L
  boundary$left[currentBoundary] <- 1L
  boundary$right[currentBoundary] <- n

  ## workhorse loop
  while (currentBoundary > 0L) {
    ## find largest gap
    left <- boundary$left[currentBoundary]
    right <- boundary$right[currentBoundary]
    currentBoundary <- currentBoundary - 1L
    gaps <- d[left:(right - 1L)]

    gapIdx <- which.max(gaps) + left - 1L

    ## left side
    l <- grouper(
      mass = mass[left:gapIdx],
      intensities = intensities[left:gapIdx],
      samples = samples[left:gapIdx],
      tolerance = tolerance, ...
    )
    ## further splitting needed?
    if (is.na(l[1L])) {
      currentBoundary <- currentBoundary + 1L
      boundary$left[currentBoundary] <- left
      boundary$right[currentBoundary] <- gapIdx
    } else {
      mass[left:gapIdx] <- l
    }

    ## right side
    r <- grouper(
      mass = mass[(gapIdx + 1L):right],
      intensities = intensities[(gapIdx + 1L):right],
      samples = samples[(gapIdx + 1L):right],
      tolerance = tolerance, ...
    )
    ## further splitting needed?
    if (is.na(r[1L])) {
      currentBoundary <- currentBoundary + 1L
      boundary$left[currentBoundary] <- gapIdx + 1L
      boundary$right[currentBoundary] <- right
    } else {
      mass[(gapIdx + 1L):right] <- r
    }

    ## stack size have to be increased?
    ## (should rarely happen because recursion deep is mostly < 20)
    if (currentBoundary == nBoundaries) {
      nBoundaries <- floor(nBoundaries * 1.5)
      boundary$left <- c(
        boundary$left,
        double(nBoundaries - currentBoundary)
      )
      boundary$right <- c(
        boundary$right,
        double(nBoundaries - currentBoundary)
      )
    }
  }
  mass
}

# Extracted from MALDIquant
.grouperStrict_ <- function(mass, intensities, samples, tolerance) {
  # don't accept two or more peaks of the same sample
  if (anyDuplicated(samples)) {
    return(NA)
  }
  meanMass <- mean(mass)
  # all peaks in range?
  if (any(abs(mass - meanMass) / meanMass > tolerance)) {
    return(NA)
  }
  return(meanMass)
}


deisotope <- function(mzVector, isoDelta = c(1.002, 1.0045)) {
  stopifnot(length(isoDelta) == 2)
  stopifnot(is.numeric(isoDelta))

  listDeiso <- list()
  listDeiso[[1]] <- mzVector[1]
  names(listDeiso)[1] <- mzVector[1]
  mzVector <- mzVector[-1]

  while (length(mzVector) > 0) {
    currMZ <- mzVector[1]
    mzVector <- mzVector[-1]

    foundIso <- F
    skip <- F

    for (i in 1:length(listDeiso))
    {
      if (foundIso || skip) {
        break
      }

      for (j in 1:length(listDeiso[[i]]))
      {
        deltaMZ <- abs(currMZ - listDeiso[[i]][j])
        # If the current distance is larger than the second extremal, it means
        # that we are aleady too far from a possible isotope. Thus, we can skip.
        if (deltaMZ > isoDelta[2]) {
          skip <- T
          break
        }
        # If the current m/z difference falls into the isotope distance
        # interval, then assign it as an isotope of current m/z value.
        if (deltaMZ >= isoDelta[1] & deltaMZ <= isoDelta[2]) {
          # Append the mass to the listDeiso and remove the current m/z from
          # the vector
          listDeiso[[j]] <- c(listDeiso[[j]], currMZ)
          foundIso <- T
          break
        }
      }
    }

    # If no isotope is found, then add the current m/z as a new element of the
    # de-isotoped list
    if (!foundIso) {
      listDeiso[[length(listDeiso) + 1]] <- currMZ
      names(listDeiso)[length(listDeiso)] <- currMZ
    }
  }

  return(listDeiso)
}
