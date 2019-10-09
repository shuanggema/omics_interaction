Maintainer
-------
Yaqing Xu <yaqing.xu@yale.edu>

Publication
-------
Identifying gene-environment interactions using integrative omics data for cancer outcomes (2019, In preparation)

Description
-------
1. functions.R:   + testbic : algorithm for identify regulatory modules via sparse clustering  + update_bic : to update estimated regulatory relationship by subtracting identified modules    + cv.multi_elnet: to estimate regulatory relationship theta suing multivariate regression with lasso penalization
  + joint_model: Joint model for hierarchical G-E interactions using lasso-based penalization 
2. example.R 
Example program of the proposed method using integrative omics data for identifying hierarchical interactions 