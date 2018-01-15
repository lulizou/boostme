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
#' CpG methylation values below the minCov. Default is TRUE.
#' @param trainChr which chromosome(s) to use for training.
#' default = chr3 (approximately 1.5 million CpGs). Note that the more CpGs
#' used for training, the more memory required to train and store the model.
#' @param validateChr which chromosome(s) to use for validation.
#' default = chr22 (approximately 550,000 CpGs).
#' @param testChr which chromosome(s) to use for testing.
#' default = chr4 (approximately 1.4 million CpGs).
#' @param minCov the minimum coverage required to consider a methylation
#' value trainable, default is 10.
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
#' @return a new bsseq object that has the imputed values (if imputeAndReplace
#' is TRUE). Otherwise doesn't return anything; just prints RMSE (dry run).
#'
#' @importClassesFrom bsseq BSseq
#' @importMethodsFrom bsseq pData seqnames sampleNames start width
#'
#'
#' @import bsseq
#' @import GenomicRanges
#' @import xgboost
#'
#' @export

boostme <- function(bs,
                    imputeAndReplace = TRUE,
                    trainChr = "chr1",
                    validateChr = "chr22",
                    testChr = "chr2",
                    minCov = 10,
                    sampleAvg = TRUE,
                    neighbMeth = TRUE,
                    neighbDist = TRUE,
                    featureBEDs = NULL,
                    threads = 2) {
  # checks
  stopifnot(class(bs) == "BSseq")
  if (nrow(pData(bs)) < 3)
    stop("At least 3 samples are needed to use BoostMe")

  # divide into train, validate, and test
  train <- chrSelectBSseq(bs, seqnames = trainChr)
  validate <- chrSelectBSseq(bs, seqnames = validateChr)
  test <- chrSelectBSseq(bs, seqnames = testChr)
  message(paste(Sys.time(), "Using", trainChr, "for training,",
                validateChr, "for validation, and", testChr, "for testing..."))

  for (i in 1:nrow(pData(train))) { # Train a model for each sample
    message(paste(Sys.time(), "Training on", sampleNames(bs)[i]))

    message(paste(Sys.time(), "Building features..."))
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

    # only take complete cases
    myTrain <- myTrain[complete.cases(myTrain),]
    myValidate <- myValidate[complete.cases(myValidate),]
    myTest <- myTest[complete.cases(myTest),]

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
                          silent = 1)
    testPreds <- predict(my_model, data.matrix(myTest[,-1]))
    testRMSE <- sqrt(mean((myTest[,1]-testPreds)^2))
    message(paste(Sys.time(), "Overall Testing RMSE for", sampleNames(test)[i],
                  ":", testRMSE))

    if (imputeAndReplace) {
      yCov <- as.vector(getCoverage(bs[,i]))
      replaceThese <- which(yCov < minCov)
      message(paste(Sys.time(), length(replaceThese), "CpGs",
                    "out of", length(yCov), "in", sampleNames(bs)[i],
                    "CpGs have minCov <", minCov))

      # build features for all sites in the sample
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
      imputedValues <- predict(my_model, data.matrix(dat[,-1]))

      # replace
      message(paste(Sys.time(), "Replacing..."))
      newY <- y
      newY[enoughInfoToImpute] <- imputedValues

      bs[,i] <- BSseq(M = matrix(newY), Cov = getCoverage(bs[,i]),
                      pData = pData(bs[,i]), gr = granges(bs[,i]))
    }

  }
  bs
}
