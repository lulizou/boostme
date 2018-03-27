#' Performance metrics
#'
#' Helper function to calculate performance metrics from labelled predictions.
#'
#' @param pred a vector containing the predicted methylation proportions.
#' @param actual a vector containing the actual methylated proportions.
#' @return a list of RMSE, AUROC, AUPRC, and accuracy
#'
#'
#' @importFrom PRROC pr.curve roc.curve
#'
#' @import PRROC
#'
#' @export

performanceMetrics <- function(pred, actual) {
  rmse <- sqrt(mean((actual - pred)^2))
  # binarize for everything else
  actualBinary <- ifelse(actual >= 0.5, 1, 0)
  fg <- pred[actualBinary == 1]
  bg <- pred[actualBinary == 0]
  auroc <- roc.curve(scores.class0 = fg, scores.class1 = bg)
  auprc <- pr.curve(scores.class0 = fg, scores.class1 = bg)
  acc <- length(which(ifelse(pred >= 0.5, 1, 0) == actualBinary))/
    length(pred)
  list(rmse = rmse, auroc = auroc$auc, auprc = auprc$auc.integral, acc = acc)
}
