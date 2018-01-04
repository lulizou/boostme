#' Add feature to the model from a BED file
#'
#' Helper function to add features from BED files. All features will become
#' binary features, 0 indicating no overlap, 1 indicating overlap.
#' If the feature has more than one factor, each factor will become its own
#' binary feature.
#'
#' @param features a data.frame to add the features to
#' @param ranges a \code{GRanges} object with the ranges to overlap features
#' (must have same length as the number of rows in the features data frame)
#' @param featureName what the feature should be called
#' @param file path to the BED File with 4 columns (chr, start, end, feature)
#'
#' @return the original features data.frame with the new feature(s) added from
#' the BED file as the last column(s)
#'
#' @importClassesFrom bsseq BSseq
#' @importMethodsFrom bsseq pData seqnames sampleNames start width
#'
#' @import bsseq
#' @import GenomicRanges
#'
#' @export


addFeatureFromBED <- function(features,
                              ranges,
                              featureName = NULL,
                              file = NULL) {
  if (is.null(file)) {
    stop("Must specify the BED file to add using the file param")
  }
  if (length(ranges) != nrow(features)) {
    stop("Ranges are not the same length as the features data frame = bad")
  }
  dat <- read.table(file, header = F, sep = "\t")
  colnames(dat)[1:3] <- c("chr", "start", "end")
  colnames(dat)[4] <- featureName
  dat$start <- dat$start + 1 # make sure in same coordinates
  grl <- split(makeGRangesFromDataFrame(dat), dat[,4])
  for (i in 1:length(grl)) {
    overlaps <- findOverlaps(ranges, grl[i])
    eval(parse(text = paste0("features$`", names(grl)[i], "` <- 0")))
    eval(parse(text = paste0("features$`", names(grl)[i],
                            "`[queryHits(overlaps)] <- 1")))
  }
  features
}
