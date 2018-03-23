#' Feature construction
#'
#' Helper function to make a data frame ready for training/validation/testing
#' from a bsseq object that contains the CpG/DNA features for the model
#'
#' @param bs a \code{BSseq} object
#' @param sample which sample to construct features for; default is 1, the
#' first sample.
#' @param minCov the minimum coverage, below which values will be set to NA.
#' Default is 10.
#' @param sampleAvg boolean of whether to not to include the sample average
#' as a feature. Default is TRUE.
#' @param neighbMeth boolean of whether or not to include nearest non-missing
#' neighboring CpG methylation values. Default is TRUE.
#' @param neighbDist boolean of whether or not to include nearest non-missing
#' neighboring CpG distances. Default is TRUE.
#' @param featureBEDs optional named list of paths to BED files to be included
#' as features in the model. Names are used as the feature name;
#' e.g. list(chromState = "chromatinStates.bed")
#' @return a data.frame, each row is a CpG; first column is the actual beta
#' value, subsequent columns are features
#'
#'
#' @importClassesFrom bsseq BSseq
#'
#' @importFrom dplyr bind_cols
#'
#' @import bsseq
#' @import GenomicRanges
#'
#' @export

constructFeatures <- function(bs,
                              sample = 1,
                              minCov = 10,
                              sampleAvg = TRUE,
                              neighbMeth = TRUE,
                              neighbDist = TRUE,
                              featureBEDs = NULL) {
  stopifnot(class(bs) == "BSseq")
  if (nrow(pData(bs)) < 3 & sampleAvg == TRUE)
    stop("At least 3 samples are needed to use the sample average feature")
  y <- as.vector(getMeth(bs[, sample], type = "raw"))
  coverage <- as.vector(getCoverage(bs[, sample], type = "Cov"))
  y[which(coverage < minCov)] <- NA
  ranges <- granges(bs)
  features <- as.data.frame(ranges)
  if (sampleAvg) {
    # compute the sample average feature;
    # need to make a copy of the data 'cause bsseq doesn't allow NAs
    otherMeths <- getMeth(bs[, -sample], type = "raw")
    otherCovs <- getCoverage(bs[, -sample], type = "Cov")
    if (class(otherMeths) == "DelayedMatrix") {
      otherMeths <- as.matrix(otherMeths)
      otherCovs <- as.matrix(otherCovs)
    }
    otherMeths[which(otherCovs < minCov)] <- NA
    sampleAverage <- rowMeans(otherMeths, na.rm = T)
    features$sampleAvg <- sampleAverage
  }
  if (neighbMeth | neighbDist) {
    neighb <- neighbors(as.matrix(data.frame(chr = as.integer(as.character(
      gsub("chr", "",
           features$seqnames))),
      pos = features$start,
      meth = y)))
    if (neighbMeth) {
      features$upBeta <- neighb$upstreamBetas
      features$downBeta <- rev(neighb$downstreamBetas)
      # TODO: make it so don't need to reverse after running neighbors
    }
    if (neighbDist) {
      features$upDist <- neighb$upstreamDistances
      features$downDist <- rev(neighb$downstreamDistances)
    }
  }
  if (!is.null(featureBEDs)) {
    for (i in 1:length(featureBEDs)) {
      features <- addFeatureFromBED(features, ranges,
                                    featureName = names(featureBEDs)[i],
                                    file = featureBEDs[[i]])
    }
  }
  dat <- bind_cols(data.frame(y = y), features)
  dat <- dat[, -which(names(dat) %in% c("seqnames", "start", "end", "width",
                                        "strand"))]
  dat
}
