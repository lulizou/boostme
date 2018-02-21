#include <Rcpp.h>
#include <vector>
using namespace Rcpp;

//' Get the nearest non-missing neighbor betas and distances
//'
//' This function finds the nearest non-missing neighboring CpG methylation
//' values upstream and downstream for each CpG.
//'
//' @param posAndBetas A numeric matrix that looks like:
//' chr  pos  beta
//' 1    5000 0.1
//' 1    5012 0.9
//' ...
//' 3    10020 0.01
//' 3    10112 0.09
//' (with some NAs in there)
//' @export


// [[Rcpp::export]]
List neighbors(NumericMatrix posAndBetas) {
  int nrow = posAndBetas.nrow();
  std::vector<double> upstreamBetas;
  std::vector<double> upstreamDistances;
  std::vector<double> downstreamBetas;
  std::vector<double> downstreamDistances;
  upstreamBetas.reserve(nrow);
  upstreamDistances.reserve(nrow);
  downstreamBetas.reserve(nrow);
  downstreamDistances.reserve(nrow);
  double previousBeta;
  double previousStart;
  bool everythingBeforeIsNA = true;
  for (int i = 0; i < nrow; ++i) { // get upstream betas
    if (i != nrow - 1) { // if not the last value
      if (posAndBetas(i, 0) == posAndBetas(i+1, 0)) {  // if same chromosome
        if (everythingBeforeIsNA) {
          upstreamBetas.push_back(NA_REAL);
          upstreamDistances.push_back(NA_REAL);
          if (!(NumericVector::is_na(posAndBetas(i, 2)))) {
            everythingBeforeIsNA = false;
          }
        } else {
          upstreamBetas.push_back(previousBeta);
          upstreamDistances.push_back(posAndBetas(i, 1) - previousStart);
        }
        if (!(NumericVector::is_na(posAndBetas(i, 2)))) {
          previousBeta = posAndBetas(i, 2);
          previousStart = posAndBetas(i, 1);
        }
      } else { // if not the same chromosome, start over
        everythingBeforeIsNA = true;
        upstreamBetas.push_back(previousBeta);
        upstreamDistances.push_back(posAndBetas(i, 1) - previousStart);
      }
    } else { // if the last value
      upstreamBetas.push_back(previousBeta);
      upstreamDistances.push_back(posAndBetas(i, 1) - previousStart);
    }
  }

  everythingBeforeIsNA = true;
  // TODO: maybe looping twice isn't the *most* efficient, but reserving space
  // for the vectors seems to make it speedy enough
  for (int i = nrow - 1; i >= 0; --i) { // get downstream betas
    if (i != 0) { // if not the first value
      if (posAndBetas(i, 0) == posAndBetas(i-1, 0)) { // if the same chromosome
        if (everythingBeforeIsNA) {
          downstreamBetas.push_back(NA_REAL);
          downstreamDistances.push_back(NA_REAL);
          if (!(NumericVector::is_na(posAndBetas(i, 2)))) {
            everythingBeforeIsNA = false;
          }
        } else {
          downstreamBetas.push_back(previousBeta);
          downstreamDistances.push_back(previousStart - posAndBetas(i, 1));
        }
        if (!(NumericVector::is_na(posAndBetas(i, 2)))) {
          previousBeta = posAndBetas(i, 2);
          previousStart = posAndBetas(i, 1);
        }
      } else { // if not the same chromosome, start over
        everythingBeforeIsNA = true;
        downstreamBetas.push_back(previousBeta);
        downstreamDistances.push_back(previousStart - posAndBetas(i, 1));
      }
    } else { // if the first value
      downstreamBetas.push_back(previousBeta);
      downstreamDistances.push_back(previousStart - posAndBetas(i, 1));
    }

  }
  return Rcpp::List::create(Rcpp::Named("upstreamBetas") = upstreamBetas,
                            Rcpp::Named("upstreamDistances") = upstreamDistances,
                            Rcpp::Named("downstreamBetas") = downstreamBetas,
                            Rcpp::Named("downstreamDistances") = downstreamDistances);
}
