#' Function for training and and imputing with a boostme model.
#'
#'
#' Uses the \code{xgboost} framework (C) Tianqi Chen, Tong He, Michael Benesty,
#' Vadim Khotilovich, Yuan Tang. Sample average feature requires at
#' least 3 samples in the bsseq object.
#'
#' @param bs a bsseq object containing the methylation & coverage values
#' as well as the features loaded into \code{pData(bs)}. If no features
#' are loaded into \code{pData(bs)}, the model will simply use neighboring
#' CpGs and the sample average of the other CpGs.
#' @param imputeAndReplace boolean of whether or not to impute and replace
#' CpG methylation values below the minCov. Default is TRUE. Set to FALSE if
#' want to do a dry run and see the RMSE for each sample.
#' @param randomCpGs boolean of whether or not to select a simple random
#' sample of CpGs genome-wide or not. Default is FALSE. If TRUE, will ignore
#' the trainChr, validateChr, and testChr parameters and select CpGs for the
#' training, validation, and test sets at random. Can modify how large
#' each of these sets will be individually using the trainSize, validateSize,
#' and testSize parameters. Defaults are 1 million CpGs each. NOTE: this takes
#' way longer to do than simply dividing by chromosome, and achieves similar
#' accuracy.
#' @param trainChr which chromosome(s) to use for training.
#' default = chr3 (approximately 1.5 million CpGs). Note that the more CpGs
#' used for training, the more memory required to train and store the model.
#' @param validateChr which chromosome(s) to use for validation.
#' default = chr22 (approximately 550,000 CpGs).
#' @param testChr which chromosome(s) to use for testing.
#' default = chr4 (approximately 1.4 million CpGs).
#' @param trainSize integer of how many CpGs to use for train set. Default is
#' 1 million. NOTE: only kicks into effect when randomCpGs = TRUE.
#' @param validateSize integer of how many CpGs to use for validation set.
#' Default is 1 million. NOTE: only kicks into effect when randomCpGs = TRUE.
#' @param testSize integer of how many CpGs to use for test set. Default is 1
#' million. NOTE: only kicks into effect when randomCpGs = TRUE.
#' @param minCov the minimum coverage required to consider a methylation
#' value trainable, default is 10 (i.e. 10 total reads at a CpG). Also used
#' as the cutoff below which to impute and replace the methylation value, given
#' that none of the features used for that CpG are NA. E.g. if a CpG has
#' coverage of 2, but sampleAvg = TRUE and < 2 samples have coverage >=10 for
#' that CpG, then that CpG's value will not be imputed and replaced.
#' @param sampleAvg boolean of whether to not to include the sample average
#' as a feature. Default is TRUE.
#' @param neighbMeth boolean of whether or not to include nearest non-missing
#' neighboring CpG methylation values. Default is TRUE.
#' @param neighbDist boolean of whether or not to include nearest non-missing
#' neighboring CpG distances. Default is TRUE.
#' @param featureBEDs optional named list of paths to BED files to be included
#' as features in the model. Names are used as the feature name;
#' e.g. list(chromState = "chromatinStates.bed")
#' @param threads (optional) number of threads to use for training. default = 2
#' @return a data frame that has the imputed values (if imputeAndReplace
#' is TRUE). Otherwise doesn't return anything; just prints RMSE for each
#' sample (dry run).
#'
#' @importClassesFrom bsseq BSseq
#' @importMethodsFrom bsseq pData seqnames sampleNames start width
#' @importFrom dplyr bind_rows bind_cols
#'
#' @import bsseq
#' @import GenomicRanges
#' @import xgboost
#'
#' @export

boostme <- function(bs,
                    imputeAndReplace = TRUE,
                    randomCpGs = FALSE,
                    trainChr = "chr1",
                    validateChr = "chr22",
                    testChr = "chr2",
                    trainSize = 1000000,
                    validateSize = 1000000,
                    testSize = 1000000,
                    minCov = 10,
                    sampleAvg = TRUE,
                    neighbMeth = TRUE,
                    neighbDist = TRUE,
                    featureBEDs = NULL,
                    threads = 2) {
  # checks
  stopifnot(class(bs) == "BSseq")
  if (nrow(pData(bs)) < 3 & sampleAvg == TRUE)
    stop("At least 3 samples are needed to use the sample average feature")

  # only use the autosome
  bs <- chrSelectBSseq(bs, seqnames = paste("chr", 1:22, sep=""))

  imputed <- getMeth(bs, type = "raw")
  message(paste(Sys.time(), "Extracting positions from bs file"))
  rownames(imputed) <- as.character(granges(bs))
  for (i in 1:nrow(pData(bs))) { # Train a model for each sample
    # TODO: add in parallel option for this instead of loop
    message(paste(Sys.time(), "Training on", sampleNames(bs)[i]))
    message(paste(Sys.time(), "Building features..."))
    if (randomCpGs) { # need to randomly sample and also make sure have
      # complete cases for the required amount of CpGs.
      datSize <- 0
      targetSize <- trainSize + validateSize + testSize
      alreadySampled <- NULL # for making sure no overlapping sampling
      myAll <- NULL
      l <- c(1:length(bs))
      while (datSize < targetSize) {
        sampleRows <- sample(l[!l %in% alreadySampled],
                             size = (targetSize - datSize))
        alreadySampled <- append(alreadySampled, sampleRows)
        bigBS <- bs[sampleRows, ]
        bigBS <- constructFeatures(bigBS, sample = i, minCov = minCov,
                                   sampleAvg = sampleAvg,
                                   neighbMeth = neighbMeth,
                                   neighbDist = neighbDist,
                                   featureBEDs = featureBEDs)
        myAll <- bind_rows(myAll, complete.cases(bigBS))
        datSize <- nrow(myAll)
      }
      rm(bigBS)
      myTrain <- myAll[1:trainSize, ]
      myValidate <- myAll[(trainSize + 1):(trainSize + validateSize), ]
      myTest <- myAll[(trainSize + validateSize + 1):targetSize, ]
    } else {
      message(paste(Sys.time(), "Using", trainChr, "for training,",
                    validateChr, "for validation, and", testChr, "for testing..."))
      train <- chrSelectBSseq(bs, seqnames = trainChr)
      validate <- chrSelectBSseq(bs, seqnames = validateChr)
      test <- chrSelectBSseq(bs, seqnames = testChr)
      myTrain <- constructFeatures(train, sample = i, minCov = minCov,
                                   sampleAvg = sampleAvg,
                                   neighbMeth = neighbMeth,
                                   neighbDist = neighbDist,
                                   featureBEDs = featureBEDs)
      myValidate <- constructFeatures(validate, sample = i, minCov = minCov,
                                      sampleAvg = sampleAvg,
                                      neighbMeth = neighbMeth,
                                      neighbDist = neighbDist,
                                      featureBEDs = featureBEDs)
      myTest <- constructFeatures(test, sample = i, minCov = minCov,
                                  sampleAvg = sampleAvg,
                                  neighbMeth = neighbMeth,
                                  neighbDist = neighbDist,
                                  featureBEDs = featureBEDs)

      # only take complete cases (i.e. no NA features)
      myTrain <- myTrain[complete.cases(myTrain), ]
      myValidate <- myValidate[complete.cases(myValidate), ]
      myTest <- myTest[complete.cases(myTest), ]
    }

    # convert data frames to xgboost-ready matrices
    myTrain[] <- lapply(myTrain, as.numeric)
    myValidate[] <- lapply(myValidate, as.numeric)
    myTest[] <- lapply(myTest, as.numeric)

    dtrain <- xgb.DMatrix(data = data.matrix(myTrain[, -1]),
                          label = myTrain[, 1])
    dvalidate <- xgb.DMatrix(data = data.matrix(myValidate[, -1]),
                             label = myValidate[, 1])
    dtest <- xgb.DMatrix(data = data.matrix(myTest[, -1]),
                         label = myTest[, 1])

    # train the model
    watchlist <- list(train = dtrain, validate = dvalidate)
    message(paste(Sys.time(), "Training the model..."))
    my_model <- xgb.train(data = dtrain,
                          nthread = threads,
                          nrounds = 500,
                          watchlist = watchlist,
                          objective = 'reg:linear',
                          early_stopping_round = 10,
                          verbose = 0,
                          silent = 1)
    testPreds <- predict(my_model, data.matrix(myTest[, -1]))
    testRMSE <- sqrt(mean((myTest[, 1]-testPreds)^2))
    message(paste(Sys.time(), "Overall Testing RMSE for", sampleNames(test)[i],
                  ":", testRMSE))

    if (imputeAndReplace) {
      yCov <- as.vector(getCoverage(bs[, i]))
      replaceThese <- which(yCov < minCov)
      message(paste(Sys.time(), length(replaceThese), "CpGs",
                    "out of", length(yCov), "(",
                    length(replaceThese)/length(yCov), "%) in",
                    sampleNames(bs)[i], "CpGs have minCov <", minCov))

      # build features for all sites in the sample
      message(paste(Sys.time(), "Constructing features for all sites in the sample
                    (may take a while)"))
      dat <- constructFeatures(bs, sample = i, minCov = minCov,
                               sampleAvg = sampleAvg,
                               neighbMeth = neighbMeth,
                               neighbDist = neighbDist,
                               featureBEDs = featureBEDs)
      enoughInfoToImpute <- replaceThese[which(
        complete.cases(dat[replaceThese, -1]))]
      message(paste(Sys.time(), "Able to impute", length(enoughInfoToImpute),
                    "out of", length(replaceThese)))
      dat <- dat[enoughInfoToImpute,]
      dat[] <- lapply(dat, as.numeric)

      # impute
      message(paste(Sys.time(), "Imputing..."))
      imputedValues <- predict(my_model, data.matrix(dat[, -1]))

      # replace
      message(paste(Sys.time(), "Replacing..."))
      newY <- getMeth(bs[, i], type = "raw")
      newY[enoughInfoToImpute] <- imputedValues
      imputed[, i] <- newY
    }
  }
  imputed
}
