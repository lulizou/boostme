#' Add feature to the model from a BED file
#'
#' Helper function to add features from BED files. All features will become
#' binary features, 0 indicating no overlap, 1 indicating overlap.
#' If the feature has more than one factor, each factor will become its own
#' binary feature.
#'
#' @param features a data.frame to add the features to.
#' @param ranges a \code{GRanges} object with the ranges to overlap features
#' (must have same length as the number of rows in the features data frame).
#' @param featureName what the feature should be called.
#' @param file path to the BED File with >=3 columns (chr, start, end,
#' feature(s)). If only 3 columns, the name of the file will be used as the
#' feature name. If 4 columns, and the 4th column has multiple factors,
#' uses the factor names of the 4th column as the feature names. Does not
#' support files with features in > the 4th column (please separate into
#' multiple BED files).
#'
#' @return the original features data.frame with the new feature(s) added from
#' the BED file as the last column(s)
#'
#' @importClassesFrom bsseq BSseq
#' @importClassesFrom S4Vectors Hits
#' @importMethodsFrom bsseq pData seqnames sampleNames start width
#'
#' @import bsseq
#' @import GenomicRanges
#' @import S4Vectors
#'
#' @export


addFeatureFromBED <- function(features,
                              ranges,
                              featureName = NULL,
                              file = NULL) {
  # checks
  if (is.null(file)) {
    stop("Must specify the BED file to add using the file param")
  }
  if (length(ranges) != nrow(features)) {
    stop("Ranges are not the same length as the features data frame = bad")
  }
  dat <- read.table(file, header = F, sep = "\t")
  if (ncol(dat) < 3) {
    stop(paste("BED file", file, "does not have >= 3 columns."))
  }
  if (ncol(dat) > 4) {
    message(paste("Warning: there are more than 4 columns in", file,
                  "so ignoring everything after the 4th column."))
  }
  colnames(dat)[1:3] <- c("chr", "start", "end")
  dat$start <- dat$start + 1 # make sure in same coordinates
  if (ncol(dat) == 4) {
    for (i in 4:ncol(dat)) {
      grl <- split(makeGRangesFromDataFrame(dat), dat[, i])
      for (j in 1:length(grl)) {
        overlaps <- findOverlaps(ranges, grl[j])
        eval(parse(text = paste0("features$`", names(grl)[j], "` <- 0")))
        eval(parse(text = paste0("features$`", names(grl)[j],
                                 "`[queryHits(overlaps)] <- 1")))
      }
    }
  } else { # if the file only has 3 columns, assume that it's one feature
    # without multiple factors and use the name of the file as the feature name
    grl <- makeGRangesFromDataFrame(dat)
    overlaps <- findOverlaps(ranges, grl)
    eval(parse(text = paste0("features$`", basename(file), "` <- 0")))
    eval(parse(text = paste0("features$`", basename(file),
                             "`[queryHits(overlaps)] <- 1")))
  }

  features
}
