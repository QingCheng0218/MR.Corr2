MR.Corr2
=======

**MR.Corr2** is a package for two-sample Mendelian randomization  using correlated instrumental variants  accounting correlated horizontal pleiotropy.

Installation
============
Install the development version of **MR.Corr2** by use of the 'devtools' package. Note that **MR.Corr2** depends on the 'Rcpp' and 'RcppArmadillo' package, which also requires appropriate setting of Rtools and Xcode for Windows and Mac OS/X, respectively.
```
library(devtools)
install_github("QingCheng0218/MR.Corr2")
```

If you have errors installing this package on Linux, try the following commands in R:
```
Sys.setlocale(category = "LC_ALL", locale = "English_United States.1252") 
Sys.setenv(LANG="en_US.UTF-8")
Sys.setenv(LC_ALL="en_US.UTF-8")

library(devtools)
install_github("QingCheng0218/MR.Corr2")
```

Usage
=========
The ['MR.Corr2' vignette](https://github.com/QingCheng0218/MR.Corr2/blob/master/vignettes/MR.Corr2.pdf) will provide a good start point for two-sample Mendelian randomization analysis using **MR.Corr2** package. 

References
==========
Qing Cheng, Tingting Qiu, Xiaoran Chai, Baoluo Sun, Yingcun Xia, Xingjie<sup>+</sup>, ShiJin Liu<sup>+</sup>. (2021) [MR-Corr$^$2: a two-sample Mendelian randomization method that accounts for correlated horizontal pleiotropy using correlated instrumental variants.](https://academic.oup.com/bioinformatics/advance-article-abstract/doi/10.1093/bioinformatics/btab646/6367765)

Development
===========

This package is developed and maintained by Qing Cheng (qingcheng0218@gmail.com). 
