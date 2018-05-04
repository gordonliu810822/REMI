REMI
===
REMI is a regression with marginal information approach with applications in genome-wide association studies (GWAS).

Installation 
===========

To install the development version of REMI, it's easiest to use the 'devtools' package. Note that REMI depends on the 'Rcpp' package, which also requires appropriate setting of Rtools and Xcode for Windows and Mac OS/X, respectively.

```
#install.packages("devtools")
library(devtools)
install_github("gordonliu810822/REMI")
```

Usage
===========

The ['REMI' vignette](https://github.com/gordonliu810822/REMI/blob/master/vignettes/REMI_packages.pdf) will provide a good start point for the genetic analysis using REMI package. The following help page will also provide quick references for REMI package and the example command lines:

```
library(REMI)
package?REMI
```

Reproducing Results of Huang et al. (2018)
===========
All the simulation results can be reproduced by using the code at [simulation](https://drive.google.com/drive/folders/1ic9Q7Onq0iDNSkex-S4V_YRv-_aL4UxL?usp=sharing). Before running simulation to reproduce the results, please familarize yourself with REMI using [demo.R](https://drive.google.com/file/d/1LY6W7nEZROd7dofxClDagXp5W_2U2DWK/view?usp=sharing) and ['REMI' vignette](https://github.com/gordonliu810822/REMI/blob/master/vignettes/REMI_packages.pdf). Simulation results can be reproduced using [simulation.R](https://drive.google.com/file/d/1VDlborxyv7Lm3X6HhtZ6Ap-7nKHhuDUC/view?usp=sharing) with a batch script [batch_submit.txt](https://drive.google.com/file/d/1ggsW7Xc8VXxD7WxrcKlmI9uqPAzpz6sW/view?usp=sharing). Then Figures 1 and 2 in Huang et al. (2018) can be reproduced using [plotsInPaper.R](https://drive.google.com/file/d/1KoH9_uCb_QtpC8EK3yzUSH18oBdr9a0v/view?usp=sharing).

In addition, we provide summary statistics data to produce results from the real data analysis. All summary statistics are stored in [link](https://drive.google.com/drive/folders/1TgU4M9k8gwbgwNwd7P9IL-XjvsCA2QHT?usp=sharing). The solution path of NFBC data can be reproduced by running [analysis_NFBC.R](https://drive.google.com/file/d/1oqHbbMRnMz0td-oUSJy0nimNbyyenfW1/view?usp=sharing). For the analysis of ten traits in Figures 4 and 5, the results can be reproduced by running [analysis_trait.R](https://drive.google.com/file/d/13vXiljTLjff9tZMnhJCtrMv6DE5eNK_-/view?usp=sharing).


References
==========
J. Huang, Y. Jiao, J. Liu and C. Yang. (2018) [REMI: Regression with marginal information with applications in genome-wide association studies](https://arxiv.org/pdf/1805.01284.pdf).


Development
===========

This package is developed and maintained by Jin Liu (jin.liu@duke-nus.edu.sg).
