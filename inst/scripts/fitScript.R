#Here is how the fits in inst/extdata/ZhangFits.RData were generated:
data(Zhang)
#Unconstrained
 microMetaboInt = combi(
 list("microbiome" = zhangMicrobio, "metabolomics" = zhangMetabo),
 distributions = c("quasi", "gaussian"), compositional = c(TRUE, FALSE),
 logTransformGaussian = FALSE)
#Constrained
 microMetaboIntConstr = combi(
     list("microbiome" = zhangMicrobio, "metabolomics" = zhangMetabo),
     distributions = c("quasi", "gaussian"), compositional = c(TRUE, FALSE),
     logTransformGaussian = FALSE, covariates = zhangMetavars)
