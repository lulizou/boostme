# BoostMe: DNA methylation prediction within whole-genome bisulfite sequencing

BoostMe is a machine learning method for imputing the continuous methylation
values of CpGs sequenced at low coverage within whole-genome bisulfite
sequencing data (WGBS). BoostMe relies on [XGBoost](https://github.com/dmlc/xgboost), a previously
developed gradient boosting machine learning algorithm, and the availability of
multiple samples to achieve both higher accuracy and faster runtimes than
previously reported methods.

## Getting started
Installation (requires R >= 3.4):
```
devtools::install_github("lulizou/boostme")
```
BoostMe accepts data input WGBS data as a `bsseq` object, which you can learn more about [here](https://bioconductor.org/packages/release/bioc/html/bsseq.html). Highest accuracy is achieved when multiple samples are used, but if you want, imputation can still be done using only neighboring CpG information by setting `sampleAvg = FALSE`.


## More information

To learn more about BoostMe, see the manuscript:

> Zou, L.S., Erdos, M.R., Taylor, D.L., Chines, P.S., Varshney, A., The
> McDonnell Genome Institute, Parker, S.C.J., Collins, F.S., and Didion, J.P.
> BoostMe accurately predicts DNA methylation values in whole-genome bisulfite
> sequencing of multiple human tissues. *bioRxiv* 207506, 2017.
> [10.1101/207056](https://www.biorxiv.org/content/early/2018/01/12/207506)
