.ppmDist <- function(m1, m2, method="average")
{
  validMethods <- c("average", "mutual", "fixed")
  if (!(method %in% validMethods))
  {
    stop("Accepted values for 'ref' are:", paste0(validMethods, collapse = ","), "\n")
  }
  if (method == 'average')
  {
    refm <- mean(c(m1, m2))
    d <- mean(c(abs(m1 - refm) / refm, abs(m2 - refm) / refm)) * 1e6
  } else if (method == 'mutual')
  {
    d <- mean(c(abs(m1 - m2) / m1, abs(m1 - m2) / m2)) * 1e6
  } else if (method == 'fixed')
  {
    d <- abs(m2 - m1) / m1 * 1e6
  }
  return(d)
}

.ppmDistArray <- function(mArray, mRef, method = "mutual")
{
  validMethods <- c("mutual", "fixed")
  if (!(method %in% validMethods))
  {
    stop("Accepted values for 'method' are:", paste0(validMethods, collapse = ","), "\n")
  }
  d <- sapply(mArray, function(z) {
    .ppmDist(z, mRef, method)
  })
  return(d)
}

loadReferenceSpectra <- function(listDataDirs,
                                 filename = "avg_spectrum_within.RData")
{
  numSamples <- length(listDataDirs)

  refSpectra <- vector(mode = "list", length = numSamples)
  for (i in 1:numSamples)
  {
    cat(i, "/", numSamples, ": ", listDataDirs[i], "\n")

    ## Load the average spectrum
    data_env <- new.env()
    load(paste0(listDataDirs[i], "/", filename), envir = data_env)
    refSpectra[[i]] <- data_env$avg_spectrum

    rm(data_env)
  }
  rm(i)

  return(refSpectra)
}


matchPeaksBetweenSamples <- function(samplesPath,
                                     refPeaksList,
                                     freqThreshold=1.0,
                                     tolerance=20,
                                     deiso = TRUE,
                                     noisePeaks = NULL,
                                     inFile=NULL,
                                     outFile=NULL,
                                     verbose=TRUE)
{
  require(MALDIquant)

  .unlist <- function(x) { unlist(x, recursive = FALSE, use.names = FALSE) }

  ## Check if MALDIquant peaks

  numSamples <- length(samplesPath)
  stopifnot(length(refPeaksList) == numSamples)

  if (!is.null(noisePeaks))
  {
    stopifnot(length(noisePeaks) == numSamples)
    refPeaksList <- lapply(1:numSamples, function(z) {
      if (length(noisePeaks[[z]]) > 0)
        return(createMassPeaks(mass = refPeaksList[[z]]@mass[-noisePeaks[[z]]],
                               intensity = refPeaksList[[z]]@intensity[-noisePeaks[[z]]]))
      return(refPeaksList[[z]])
    })
  }

  ## customized binPeaks from MALDIquant
  mass <- unname(.unlist(lapply(refPeaksList, function(x)x@mass)))
  intensities <- .unlist(lapply(refPeaksList, function(x)x@intensity))
  samples <- c(rep.int(seq_along(refPeaksList), lengths(refPeaksList)))
  s <- sort.int(mass, method = "quick", index.return = TRUE)
  mass <- s$x
  intensities <- intensities[s$ix]
  samples <- samples[s$ix]

  # binning
  mass <- .binPeaks_(mass = mass, intensities = intensities, samples = samples,
                     tolerance = tolerance/1e6, grouper = .grouperStrict_)
  tableMasses <- table(mass)
  ## Group mass/intensities by sample ids
  lIdx <- split(seq_along(mass), samples)

  tableMasses <- as.numeric(names(tableMasses[tableMasses >=
                                                freqThreshold * numSamples]))
  stopifnot(!any(duplicated(tableMasses)))
  commonMasses <- sort(unique(tableMasses))

  if (verbose)
    cat("num. common masses: ", length(commonMasses), "\n")

  if (deiso)
  {
    cat('de-isotoping...\n')
    deiso_cmz_list <- deisotope(commonMasses)
    commonMasses <- as.numeric(names(deiso_cmz_list))
    cat('length m/z vector =', length(commonMasses), '\n')
  }

  ## Match intensity matrices
  .assignNewIntensityMat(samplesPath, commonMasses, mass, lIdx,
                         inFile = inFile, outFile = outFile, verbose = verbose)

  return(commonMasses)
}


# Extracted from MALDIquant
.grouperStrict_ <- function(mass, intensities, samples, tolerance)
{
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

.assignNewIntensityMat <- function(samplesPath, commonMasses, masses, indices,
                                   inFile=NULL, outFile=NULL, verbose=TRUE)
{
  numSamples <- length(samplesPath)
  for (s in 1:numSamples)
  {
    if (verbose)
      message(sprintf("%d/%d: %s", s, numSamples, samplesPath[s]))

    load(paste0(samplesPath[s], "/", inFile))

    shape <- sz
    sampleMasses <- mz
    stopifnot(length(mz) == ncol(X))
    imageShape <- shape
    stopifnot(nrow(X) == prod(shape))
    rm(shape)

    if (verbose)
      cat("dim. sample matrix:", dim(X), "\n")

    ## Assign NA to all the zeros
    X[X == 0] <- NA

    ## Match the sample masses with the common masses
    if (verbose)
      cat("num. masses sample:", length(sampleMasses), ".\n")

    matchIdx1 <- which(!is.na(match(round(commonMasses, 4),
                                    round(masses[indices[[s]]], 4))))
    matchIdx2 <- which(!is.na(match(round(masses[indices[[s]]], 4),
                                    round(commonMasses, 4))))

    if (verbose)
      cat("matched ", length(matchIdx1), " masses.\n")

    ## Fill the new intensity matrix
    Xmatched <- matrix(NA, nrow(X), length(commonMasses),
                       dimnames = list(NULL, commonMasses))
    Xmatched[, matchIdx1] <- X[, matchIdx2]
    X <- Xmatched
    rm(Xmatched)
    gc()

    ## Save

    cat("saving the matched intensity matrix...\n")

    if (is.null(outFile))
      outFile <- "X_match_between.RData"

    tryCatch({
      if (file.exists(paste0(samplesPath[s], "/", outFile)))
        file.remove(paste0(samplesPath[s], "/", outFile))
      mz <- commonMasses
      shape <- imageShape
      save(list = c('X', 'mz', 'shape'), file = paste0(samplesPath[s], "/", outFile))
    }, error=function(e) {
      stop("Error saving ", paste0(samplesPath[s], "/", outFile), "\n")
    })
  }
  return(invisible(NULL))
}


.binPeaks_ <- function(mass, intensities, samples, tolerance,
                       grouper=.grouperStrict_, ...) {
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
  boundary <- list(left=double(nBoundaries), right=double(nBoundaries))

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
    l <- grouper(mass=mass[left:gapIdx],
                 intensities=intensities[left:gapIdx],
                 samples=samples[left:gapIdx],
                 tolerance=tolerance, ...)
    ## further splitting needed?
    if (is.na(l[1L])) {
      currentBoundary <- currentBoundary + 1L
      boundary$left[currentBoundary] <- left
      boundary$right[currentBoundary] <- gapIdx
    } else {
      mass[left:gapIdx] <- l
    }

    ## right side
    r <- grouper(mass=mass[(gapIdx + 1L):right],
                 intensities=intensities[(gapIdx + 1L):right],
                 samples=samples[(gapIdx + 1L):right],
                 tolerance=tolerance, ...)
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
      boundary$left <- c(boundary$left,
                         double(nBoundaries - currentBoundary))
      boundary$right <- c(boundary$right,
                          double(nBoundaries - currentBoundary))
    }
  }
  mass
}


deisotope <- function(mzVector, isoDelta = c(1.002, 1.0045))
{
  stopifnot(length(isoDelta) == 2)
  stopifnot(is.numeric(isoDelta))

  listDeiso <- list()
  listDeiso[[1]] <- mzVector[1]
  names(listDeiso)[1] <- mzVector[1]
  mzVector <- mzVector[-1]

  while (length(mzVector) > 0)
  {
    currMZ <- mzVector[1]
    mzVector <- mzVector[-1]

    foundIso <- F
    skip <- F

    for (i in 1:length(listDeiso))
    {
      if (foundIso || skip)
      {
        break
      }

      for (j in 1:length(listDeiso[[i]]))
      {
        deltaMZ <- abs(currMZ - listDeiso[[i]][j])
        # If the current distance is larger than the second extremal, it means
        # that we are aleady too far from a possible isotope. Thus, we can skip.
        if (listDeiso[[i]][j] - currMZ > isoDelta[2])
        {
          skip <- T
          break
        }
        # If the current m/z difference falls into the isotope distance
        # interval, then assign it as an isotope of current m/z value.
        if (deltaMZ >= isoDelta[1] & deltaMZ <= isoDelta[2])
        {
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
    if (!foundIso)
    {
      listDeiso[[length(listDeiso) + 1]] <- currMZ
      names(listDeiso)[length(listDeiso)] <- currMZ
    }
  }

  return(listDeiso)
}

matchPeaksWithCMZ <- function(listDataDirs,
                              commonMZ,
                              tolPPM = 20,
                              inputFilename = 'X_matched_within_SPUTNIK.RData',
                              outFilename = 'X_matched_with_cmz.RData',
                              verbose = TRUE)
{
  for (i in 1:length(listDataDirs))
  {
    if (verbose)
      cat(sprintf('%d/%d: %s\n', i, length(listDataDirs), listDataDirs[i]))

    temp <- new.env()
    load(paste0(listDataDirs[i], '/', inputFilename), envir = temp)  # Load X

    ## Match the m/z values with the common m/z
    matched_ix <- array(NA, length(commonMZ))
    for (j in 1:length(commonMZ))
    {
      candidate <- c()
      for (k in 1:length(temp$mz)) {
        if (.ppmDist(temp$mz[k], commonMZ[j]) <= tolPPM) {
          candidate <- c(candidate, k)
        }
      }
      if (length(candidate) == 1) {
        matched_ix[j] <- k
      } else if (length(candidate) > 1) {
        matched_ix[j] <- candidate[which.max(apply(temp$X[, k], 1, mean, na.rm = TRUE))]
      }
    }
    
    # Check that there are no duplicated assignments. In that case, reduce the tolerance.
    stopifnot(length(matched_ix[!is.na(matched_ix)]) == length(unique(matched_ix[!is.na(matched_ix)])))
    
    if (verbose)
      cat('unmatched peaks:', sum(is.na(matched_ix)), '\n')
    stopifnot(all(diff(matched_ix) > 0, na.rm = T))

    X <- temp$X[, matched_ix]
    colnames(X) <- commonMZ
    attr(X, 'mass') <- commonMZ

    save(X, file = paste0(list_samples_dirs[idx_test[i]], '/', outFilename))
    rm(temp)
    gc()
  }
  return(invisible(NULL))
}
