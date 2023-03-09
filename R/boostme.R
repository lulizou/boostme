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
#' @param impute boolean of whether or not to impute CpG methylation
#' values below the minCov. Default is TRUE. Set to FALSE if want to do a dry
#' run and see the RMSE for each sample.
#' @param randomCpGs logical value of whether or not to select a simple random
#' sample of CpGs genome-wide or not. If TRUE, will ignore the trainChr,
#' validateChr, and testChr parameters and select CpGs for the
#' training, validation, and test sets at random. Can modify how large
#' each of these sets will be individually using the trainSize, validateSize,
#' and testSize parameters. Defaults are 1 million CpGs each. Default is TRUE.
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
#' @param sampleAvg logical value of whether to not to include the sample
#' average as a feature. Default is TRUE.
#' @param neighbMeth logical value of whether or not to include nearest non-
#' missing neighboring CpG methylation values. Default is TRUE.
#' @param neighbDist logical value of whether or not to include nearest non-
#' missing neighboring CpG distances. Default is TRUE.
#' @param featureBEDs optional vector of paths to BED files to be included
#' as features in the model. If the file only has 3 columns, assumes that the
#' file is just one feature. If the file has 4 columns, 4th column is considered
#' to be the feature(s). If the 4th column has multiple factors (i.e. multiple
#' different strings) then each factor is converted to its own binary feature
#' (1 if present, else 0). Ignores all columns past the 4th. Ranges must be
#' non-overlapping (can be solved with \code{bedtools merge} or separating into
#' multiple BED files).
#' @param mask logical value of whether to mask all unimputed values as NA
#' in resulting data frame.
#' @param seed optional integer random seed used before selecting random CpGs
#' for training/validation/testing.
#' @param threads (optional) number of threads to use for training. default = 2
#' @param save (optional) file path to save metrics to (e.g. results.txt)
#' @param verbose logical value of whether to print status messages. Default is
#' TRUE.
#' @return a matrix that has the imputed values (if impute is TRUE). Otherwise
#' doesn't return anything; just prints RMSE for each sample (dry run).
#'
#' @importClassesFrom bsseq BSseq
#' @importClassesFrom DelayedArray DelayedArray DelayedMatrix
#' @importMethodsFrom bsseq pData seqnames sampleNames start width
#' @importMethodsFrom S4Vectors as.data.frame
#'
#' @importFrom PRROC pr.curve roc.curve
#' @importFrom dplyr bind_rows bind_cols sample_n
#' @importFrom DelayedArray realize
#'
#' @import DelayedArray
#' @import bsseq
#' @import GenomicRanges
#' @import xgboost
#' @import PRROC
#' @import S4Vectors
#'
#'
#' @export

boostme <- function(bs,
                    impute = TRUE,
                    randomCpGs = TRUE,
                    trainChr = "chr3",
                    validateChr = "chr22",
                    testChr = "chr4",
                    trainSize = 1000000,
                    validateSize = 1000000,
                    testSize = 1000000,
                    minCov = 10,
                    sampleAvg = TRUE,
                    neighbMeth = TRUE,
                    neighbDist = TRUE,
                    featureBEDs = NULL,
                    mask = FALSE,
                    threads = 2,
                    seed = 1,
                    save = NULL,
                    verbose = TRUE,
                    n_autosomes=22) {
  # checks
  stopifnot(class(bs) == "BSseq")
  if (trainSize + validateSize + testSize > nrow(bs)) {
    stop('trainSize + validateSize + testSize must be <= num of CpGs in bsseq')
  }
  stopifnot(trainSize + validateSize + testSize <= nrow(bs))
  if (nrow(pData(bs)) < 3 & sampleAvg == TRUE) {
    stop("At least 3 samples are needed to use the sample average feature")
  }
  if (is.null(save)) {
    if (verbose) {
      message(paste("Warning: Performance metrics will be reported but not",
                    "saved; to save metrics, specify the save parameter, e.g.",
                    "save = 'results.txt'"))
    }
  }

  metrics <- data.frame(matrix(nrow = ncol(bs), ncol = 9))
  colnames(metrics) <- c("sample",
                         "rmse_val", "auroc_val", "auprc_val", "acc_val",
                         "rmse_test", "auroc_test", "auprc_test", "acc_test")

  # only use the autosome
  if (grepl('chr', trainChr)) { # if using chr in the names
    bs <- chrSelectBSseq(bs, seqnames = paste("chr", 1:n_autosomes, sep=""))
  } else { # else don't use chr in the names
    bs <- chrSelectBSseq(bs, seqnames = seq(1,n_autosomes,1))
  }
  if (mask) {
    message(paste(Sys.time(),
                  "Masking unimputed values in resulting matrix"))
    imputed <- matrix(NA, nrow = nrow(bs), ncol = ncol(bs))
  } else {
    imputed <- getMeth(bs, type = "raw")
  }
  if (typeof(imputed)=='S4') {
    imputed <- as.data.frame(realize(imputed)@seed)
  }
  if (verbose) {
    message(paste(Sys.time(),
                  "Extracting positions from bs object (takes a bit)"))
  }
  rownames(imputed) <- as.character(granges(bs))
  for (i in 1:nrow(pData(bs))) { # Train a model for each sample
    # TODO: add in parallel option for this instead of loop (?)
    if (verbose) {
      message(paste(Sys.time(), sampleNames(bs)[i]))
      message(paste(Sys.time(), "... Building features"))
    }
    if (randomCpGs) { # need to randomly sample and also make sure have
      # complete cases for the required amount of CpGs.
      if (verbose) {
        message(paste(Sys.time(), "... Randomly selecting CpGs for",
                      "training, validation, and testing"))
      }
      targetSize <- trainSize + validateSize + testSize
      bigBS <- constructFeatures(bs, sample = i, minCov = minCov,
                               sampleAvg = sampleAvg,
                               neighbMeth = neighbMeth,
                               neighbDist = neighbDist,
                               featureBEDs = featureBEDs)
      bigBS <- bigBS[complete.cases(bigBS), ]
      set.seed(seed)
      bigBS <- sample_n(bigBS, targetSize)
      myTrain <- bigBS[1:trainSize, ]
      myValidate <- bigBS[(trainSize + 1):(trainSize + validateSize), ]
      myTest <- bigBS[(trainSize + validateSize + 1):targetSize, ]
      rm(bigBS)
    } else {
      if (verbose) {
        message(paste(Sys.time(), "... Using", trainChr, "for training,",
                      validateChr, "for validation, and", testChr,
                      "for testing"))
      }
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
    if (verbose) {
      print(str(myTrain))
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
    if (verbose) {
      message(paste(Sys.time(), "... Training the model"))
    }
    my_model <- xgb.train(data = dtrain,
                          nthread = threads,
                          nrounds = 500,
                          watchlist = watchlist,
                          objective = 'reg:squarederror',
                          early_stopping_round = 10,
                          verbose = 0)
    valPreds <- predict(my_model, data.matrix(myValidate[, -1]))
    testPreds <- predict(my_model, data.matrix(myTest[, -1]))

    valMetrics <- performanceMetrics(valPreds, myValidate[, 1])
    testMetrics <- performanceMetrics(testPreds, myTest[, 1])

    if (verbose) {
      message(paste(Sys.time(), "...... Validation RMSE:",
                    round(valMetrics$rmse, digits = 4)))
      message(paste(Sys.time(), "...... Validation AUROC:",
                    round(valMetrics$auroc, digits = 4)))
      message(paste(Sys.time(), "...... Validation AUPRC:",
                    round(valMetrics$auprc, digits = 4)))
      message(paste(Sys.time(), "...... Validation Accuracy:",
                    round(valMetrics$acc, digits = 4)))
      message(paste(Sys.time(), "...... Testing RMSE:",
                    round(testMetrics$rmse, digits = 4)))
      message(paste(Sys.time(), "...... Testing AUROC:",
                    round(testMetrics$auroc, digits = 4)))
      message(paste(Sys.time(), "...... Testing AUPRC:",
                    round(testMetrics$auprc, digits = 4)))
      message(paste(Sys.time(), "...... Testing Accuracy:",
                    round(testMetrics$acc, digits = 4)))
    }
    metrics[i, 1] <- sampleNames(bs)[i]
    metrics[i, 2:ncol(metrics)] <-
      c(valMetrics$rmse, valMetrics$auroc, valMetrics$auprc, valMetrics$acc,
        testMetrics$rmse, testMetrics$auroc, testMetrics$auprc, testMetrics$acc)

    if (impute) {
      yCov <- as.vector(getCoverage(bs[, i]))
      replaceThese <- which(yCov < minCov)
      if (verbose) {
        message(paste(Sys.time(), "...", length(replaceThese), "CpGs",
                      "out of", length(yCov), "(",
                      round(length(replaceThese)/length(yCov)*100, digits = 4),
                      "% ) in", sampleNames(bs)[i], "have minCov <", minCov))
        message(paste(Sys.time(),
                      "...... Constructing features for all sites in the",
                      "sample (may take a while)"))
      }
      # build features for all sites in the sample
      dat <- constructFeatures(bs, sample = i, minCov = minCov,
                               sampleAvg = sampleAvg,
                               neighbMeth = neighbMeth,
                               neighbDist = neighbDist,
                               featureBEDs = featureBEDs)
      enoughInfoToImpute <- replaceThese[which(
        complete.cases(dat[replaceThese, -1]))]
      if (verbose) {
        message(paste(Sys.time(), "...... Able to impute",
                      length(enoughInfoToImpute),
                      "out of", length(replaceThese), "(",
                      round(length(enoughInfoToImpute)/length(replaceThese)*100,
                            digits = 4),
                      "% )"))
      }
      dat <- dat[enoughInfoToImpute,]
      dat[] <- lapply(dat, as.numeric)

      # impute
      if (verbose) {
        message(paste(Sys.time(), "...... Imputing"))
      }
      imputedValues <- predict(my_model, data.matrix(dat[, -1]))
      imputedValues[imputedValues < 0] <- 0
      newY <- getMeth(bs[, i], type = "raw") # returns DelayedArray for large bs
      if (typeof(newY)=='S4') {
        newY <- as.data.frame(realize(newY)@seed)[, 1]
      }
      newY[enoughInfoToImpute] <- imputedValues
      if (mask) {
        newY[-enoughInfoToImpute] <- NA
      }
      imputed[, i] <- newY
    }
  }
  if (!is.null(save)) {
    write.table(metrics, file = save, quote = F, sep = "\t", row.names = F)
    if (verbose) {
      message(paste(Sys.time(), "... Saved results to", save))
    }
  }
  imputed
}
