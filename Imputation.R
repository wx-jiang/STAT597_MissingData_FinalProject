require(DirichletReg)
require(mice)

#### Load Data ----
load("Data/MetaData_and_Missing_Data")

MetaData = dataStore[[1]]
arch3mis2pct10 = dataStore[[2]]
arch3mis2pct30 = dataStore[[3]]
arch3mis2pct90 = dataStore[[4]]

arch6mis2pct10 = dataStore[[5]]
arch6mis2pct30 = dataStore[[6]]
arch6mis2pct90 = dataStore[[7]]

arch6mis4pct10 = dataStore[[8]]
arch6mis4pct30 = dataStore[[9]]
arch6mis4pct90 = dataStore[[10]]

arch12mis2pct10 = dataStore[[11]]
arch12mis2pct30 = dataStore[[12]]
arch12mis2pct90 = dataStore[[13]]

arch12mis4pct10 = dataStore[[14]]
arch12mis4pct30 = dataStore[[15]]
arch12mis4pct90 = dataStore[[16]]

arch12mis11pct10 = dataStore[[17]]
arch12mis11pct30 = dataStore[[18]]
arch12mis11pct90 = dataStore[[19]]

#### Complete case data ----
# output a list of complete case arch3 matrices
# can be accessed by names
arch3_Complete_Case = lapply(
  mget(paste0("arch3mis2pct", c(10, 30, 90))),
  function(X) X[complete.cases(X),]
)

names(arch3_Complete_Case)

# output a list of complete case arch6 matrices
# can be accessed by names
arch6_Complete_Case = lapply(
  mget(paste0("arch6mis", c(rep(2, 3), rep(4, 3)), "pct", c(10, 30, 90))),
  function(X) X[complete.cases(X),]
)

names(arch6_Complete_Case)

# output a list of complete case arch12 matrices
# can be accessed by names
arch12_Complete_Case = lapply(
  mget(paste0("arch12mis", c(rep(2, 3), rep(4, 3), rep(11, 3)), "pct", c(10, 30, 90))),
  function(X) X[complete.cases(X),]
)

names(arch12_Complete_Case)

## save RData
save(list = c("arch3_Complete_Case", "arch6_Complete_Case", "arch12_Complete_Case"), 
     file = "imputed_data/Complete_case_data.RData")

# load("imputed_data/Complete_case_data.RData")


#### Regression imputation ----

# confirm no missing data in MetaData
sum(complete.cases(MetaData))
nrow(MetaData)

reg_imp <- function(d, ncol_missing, covariate_data) {
  d_with_covariate = cbind(d, covariate_data)
  complete_idx = complete.cases(d_with_covariate)
  # d_with_covariate_complete = d_with_covariate[complete_idx,]
  
  if (ncol_missing == 2) {
    # fit regression on the second last variable using complete column and covariates
    lm_model = lm(d_with_covariate[, ncol(d) - 1] ~ .,
                  data = d_with_covariate[, -c(ncol(d) - 1, ncol(d))])
    
    # impute second last variable
    d[!complete_idx, ncol(d) - 1] = predict.lm(lm_model, newdata = d_with_covariate[!complete_idx, -c(ncol(d) - 1, ncol(d))])
    # impute last column
    d[!complete_idx, ncol(d)] = 1 - rowSums(d[!complete_idx, ], na.rm = T)
    # fix the issue that the row sum might go above 1
    d = ifelse(d < 0, 0, d)
    d = d / rowSums(d)
    return(d)
    
  } else if (ncol_missing == 4) {
    # fit regression on the second last variable using complete column and covariates
    lm_model1 = lm(d_with_covariate[, ncol(d) - 1] ~ .,
                   data = d_with_covariate[, -c(ncol(d) - c(0, 1, 2, 3))])
    # fit regression on the third last variable using complete column and covariates
    lm_model2 = lm(d_with_covariate[, ncol(d) - 2] ~ .,
                   data = d_with_covariate[, -c(ncol(d) - c(0, 1, 2, 3))])
    # fit regression on the fourth last variable using complete column and covariates
    lm_model3 = lm(d_with_covariate[, ncol(d) - 3] ~ .,
                   data = d_with_covariate[, -c(ncol(d) - c(0, 1, 2, 3))])
    
    # impute second last variable
    d[!complete_idx, ncol(d) - 1] = predict.lm(lm_model1, newdata = d_with_covariate[!complete_idx, -c(ncol(d) - c(0, 1, 2, 3))])
    # impute third last variable
    d[!complete_idx, ncol(d) - 2] = predict.lm(lm_model2, newdata = d_with_covariate[!complete_idx, -c(ncol(d) - c(0, 1, 2, 3))])
    # impute fourth last variable
    d[!complete_idx, ncol(d) - 3] = predict.lm(lm_model3, newdata = d_with_covariate[!complete_idx, -c(ncol(d) - c(0, 1, 2, 3))])
    # impute last column
    d[!complete_idx, ncol(d)] = 1 - rowSums(d[!complete_idx, ], na.rm = T)
    # fix the issue that the row sum might go above 1
    d = ifelse(d < 0, 0, d)
    d = d / rowSums(d)
    return(d)
    
  } else if (ncol_missing == 11) {
    # fit regression on the second last variable using complete column and covariates
    lm_model1 = lm(d_with_covariate[, ncol(d) - 1] ~ .,
                   data = d_with_covariate[, -c(ncol(d) - c(0:10))])
    # fit regression on the third last variable using complete column and covariates
    lm_model2 = lm(d_with_covariate[, ncol(d) - 2] ~ .,
                   data = d_with_covariate[, -c(ncol(d) - c(0:10))])
    # fit regression on the fourth last variable using complete column and covariates
    lm_model3 = lm(d_with_covariate[, ncol(d) - 3] ~ .,
                   data = d_with_covariate[, -c(ncol(d) - c(0:10))])
    # fit regression on the fifth last variable using complete column and covariates
    lm_model4 = lm(d_with_covariate[, ncol(d) - 4] ~ .,
                   data = d_with_covariate[, -c(ncol(d) - c(0:10))])
    # fit regression on the sixth last variable using complete column and covariates
    lm_model5 = lm(d_with_covariate[, ncol(d) - 5] ~ .,
                   data = d_with_covariate[, -c(ncol(d) - c(0:10))])
    # fit regression on the seventh last variable using complete column and covariates
    lm_model6 = lm(d_with_covariate[, ncol(d) - 6] ~ .,
                   data = d_with_covariate[, -c(ncol(d) - c(0:10))])
    # fit regression on the eighth last variable using complete column and covariates
    lm_model7 = lm(d_with_covariate[, ncol(d) - 7] ~ .,
                   data = d_with_covariate[, -c(ncol(d) - c(0:10))])
    # fit regression on the nineth last variable using complete column and covariates
    lm_model8 = lm(d_with_covariate[, ncol(d) - 8] ~ .,
                   data = d_with_covariate[, -c(ncol(d) - c(0:10))])
    # fit regression on the tenth last variable using complete column and covariates
    lm_model9 = lm(d_with_covariate[, ncol(d) - 9] ~ .,
                   data = d_with_covariate[, -c(ncol(d) - c(0:10))])
    # fit regression on the eleventh last variable using complete column and covariates
    lm_model10 = lm(d_with_covariate[, ncol(d) - 10] ~ .,
                    data = d_with_covariate[, -c(ncol(d) - c(0:10))])
    
    # impute second last variable
    d[!complete_idx, ncol(d) - 1] = predict.lm(lm_model1, newdata = d_with_covariate[!complete_idx, -c(ncol(d) - c(0:10))])
    # impute third last variable
    d[!complete_idx, ncol(d) - 2] = predict.lm(lm_model2, newdata = d_with_covariate[!complete_idx, -c(ncol(d) - c(0:10))])
    # impute fourth last variable
    d[!complete_idx, ncol(d) - 3] = predict.lm(lm_model3, newdata = d_with_covariate[!complete_idx, -c(ncol(d) - c(0:10))])
    # impute 5 last variable
    d[!complete_idx, ncol(d) - 4] = predict.lm(lm_model4, newdata = d_with_covariate[!complete_idx, -c(ncol(d) - c(0:10))])
    # impute 6 last variable
    d[!complete_idx, ncol(d) - 5] = predict.lm(lm_model5, newdata = d_with_covariate[!complete_idx, -c(ncol(d) - c(0:10))])
    # impute 7 last variable
    d[!complete_idx, ncol(d) - 6] = predict.lm(lm_model6, newdata = d_with_covariate[!complete_idx, -c(ncol(d) - c(0:10))])
    # impute 8 last variable
    d[!complete_idx, ncol(d) - 7] = predict.lm(lm_model7, newdata = d_with_covariate[!complete_idx, -c(ncol(d) - c(0:10))])
    # impute 9 last variable
    d[!complete_idx, ncol(d) - 8] = predict.lm(lm_model8, newdata = d_with_covariate[!complete_idx, -c(ncol(d) - c(0:10))])
    # impute 10 last variable
    d[!complete_idx, ncol(d) - 9] = predict.lm(lm_model9, newdata = d_with_covariate[!complete_idx, -c(ncol(d) - c(0:10))])
    # impute 11 last variable
    d[!complete_idx, ncol(d) - 10] = predict.lm(lm_model10, newdata = d_with_covariate[!complete_idx, -c(ncol(d) - c(0:10))])
    # impute last column
    d[!complete_idx, ncol(d)] = 1 - rowSums(d[!complete_idx, ], na.rm = T)
    # fix the issue that the row sum might go above 1
    d = ifelse(d < 0, 0, d)
    d = d / rowSums(d)
    return(d)
  }
}


# output a list of regression_imp case arch3 matrices
# can be accessed by names
arch3_Regression_Imp = lapply(
  mget(paste0("arch3mis2pct", c(10, 30, 90))),
  function(X) reg_imp(X, 2, MetaData)
)

names(arch3_Regression_Imp)

# output a list of complete case arch6 matrices
# can be accessed by names
arch6_Regression_Imp = lapply(
  mget(paste0("arch6mis", c(rep(2, 3), rep(4, 3)), "pct", c(10, 30, 90))),
  function(X) {
    if (sum(colSums(is.na(X)) != 0) == 4) {
      reg_imp(X, 4, MetaData)
    } else {reg_imp(X, 2, MetaData)}
  }
)

names(arch6_Regression_Imp)

# output a list of complete case arch12 matrices
# can be accessed by names
arch12_Regression_Imp = lapply(
  mget(paste0("arch12mis", c(rep(2, 3), rep(4, 3), rep(11, 3)), "pct", c(10, 30, 90))), {
    function(X) if (sum(colSums(is.na(X)) != 0) == 11) { # number of missing columns
      reg_imp(X, 11, MetaData)
    } else if (sum(colSums(is.na(X)) != 0) == 4) {
      reg_imp(X, 4, MetaData)
    } else if (sum(colSums(is.na(X)) != 0) == 2) {
      reg_imp(X, 2, MetaData)
    }
  }
)

names(arch12_Regression_Imp)

## save RData
save(list = c("arch3_Regression_Imp", "arch6_Regression_Imp", "arch12_Regression_Imp"), 
     file = "imputed_data/Linear_regression_Imp_data.RData")

# load("imputed_data/Linear_regression_Imp_data.RData")


#### Dirichlet Regression Imputation using DirichletReg R Package -----

dirichlet_reg_imp <- function(d, n_archtype, covariate_data) {
  d_with_covariate = cbind(d, covariate_data)
  complete_idx = complete.cases(d_with_covariate)
  # d_with_covariate_complete = d_with_covariate[complete_idx,]
  n_missing = sum(colSums(is.na(d)) != 0)
  
  if (n_archtype == 3) {
    # save a copy of non-missing columns
    observed_data = d[,1, drop = F]
    
    d_with_covariate$Y = DR_data(d_with_covariate[, 1:3])
    # fit regression on the second last variable using complete column and covariates
    dirichelet_model = DirichReg(Y ~ Type + Trial + Duration, data = d_with_covariate)
    
    # impute
    d[!complete_idx, ] = predict(dirichelet_model, 
                                 newdata = d_with_covariate[!complete_idx,])
    
    # fix the issue that observed data is changed in prediction
    residual = 1 - (observed_data[!complete_idx,])
    d[!complete_idx, 2:3] = (d[!complete_idx, 2:3] / rowSums(d[!complete_idx, 2:3])) * residual
    d[!complete_idx, 1] = observed_data[!complete_idx, ]
    
    return(d)
    
  } else if (n_archtype == 6) {
    if (n_missing == 2) {
      # save a copy of non-missing columns
      observed_data = d[,1:4, drop = F]
      
      d_with_covariate$Y = DR_data(d_with_covariate[, 1:6])
      # fit regression on the second last variable using complete column and covariates
      dirichelet_model = DirichReg(Y ~ Type + Trial + Duration, data = d_with_covariate)
      
      # impute
      d[!complete_idx, ] = predict(dirichelet_model, 
                                   newdata = d_with_covariate[!complete_idx,])
      
      # fix the issue that observed data is changed in prediction
      residual = 1 - rowSums(observed_data[!complete_idx, , drop = F])
      d[!complete_idx, 5:6] = (d[!complete_idx, 5:6] / rowSums(d[!complete_idx, 5:6])) * residual
      d[!complete_idx, 1:4] = observed_data[!complete_idx, ]
      
      return(d)
      
    } else if (n_missing == 4) {
      # save a copy of non-missing columns
      observed_data = d[,1:2, drop = F]
      
      d_with_covariate$Y = DR_data(d_with_covariate[, 1:6])
      # fit regression on the second last variable using complete column and covariates
      dirichelet_model = DirichReg(Y ~ Type + Trial + Duration, data = d_with_covariate)
      
      # impute
      d[!complete_idx, ] = predict(dirichelet_model, 
                                   newdata = d_with_covariate[!complete_idx,])
      
      # fix the issue that observed data is changed in prediction
      residual = 1 - rowSums(observed_data[!complete_idx, , drop = F])
      d[!complete_idx, 3:6] = (d[!complete_idx, 3:6] / rowSums(d[!complete_idx, 3:6])) * residual
      d[!complete_idx, 1:2] = observed_data[!complete_idx, ]
      
      return(d)
    }
    
  } else if (n_archtype == 12) {
    if (n_missing == 2) {
      # save a copy of non-missing columns
      observed_data = d[,1:10, drop = F]
      
      d_with_covariate$Y = DR_data(d_with_covariate[, 1:12])
      # fit regression on the second last variable using complete column and covariates
      dirichelet_model = DirichReg(Y ~ Type + Trial + Duration, data = d_with_covariate)
      
      # impute
      d[!complete_idx, ] = predict(dirichelet_model, 
                                   newdata = d_with_covariate[!complete_idx,])
      
      # fix the issue that observed data is changed in prediction
      residual = 1 - rowSums(observed_data[!complete_idx, , drop = F])
      d[!complete_idx, 11:12] = (d[!complete_idx, 11:12] / rowSums(d[!complete_idx, 11:12])) * residual
      d[!complete_idx, 1:10] = observed_data[!complete_idx, ]
      
      return(d)
      
    } else if (n_missing == 4) {
      # save a copy of non-missing columns
      observed_data = d[,1:8, drop = F]
      
      d_with_covariate$Y = DR_data(d_with_covariate[, 1:12])
      # fit regression on the second last variable using complete column and covariates
      dirichelet_model = DirichReg(Y ~ Type + Trial + Duration, data = d_with_covariate)
      
      # impute
      d[!complete_idx, ] = predict(dirichelet_model, 
                                   newdata = d_with_covariate[!complete_idx,])
      
      # fix the issue that observed data is changed in prediction
      residual = 1 - rowSums(observed_data[!complete_idx, , drop = F])
      d[!complete_idx, 9:12] = (d[!complete_idx, 9:12] / rowSums(d[!complete_idx, 9:12])) * residual
      d[!complete_idx, 1:8] = observed_data[!complete_idx, ]
      
      return(d)
      
    } else if (n_missing == 11) {
      # save a copy of non-missing columns
      observed_data = d[,1, drop = F]
      
      d_with_covariate$Y = DR_data(d_with_covariate[, 1:12])
      # fit regression on the second last variable using complete column and covariates
      dirichelet_model = DirichReg(Y ~ Type + Trial + Duration, data = d_with_covariate)
      
      # impute
      d[!complete_idx, ] = predict(dirichelet_model, 
                                   newdata = d_with_covariate[!complete_idx,])
      
      # fix the issue that observed data is changed in prediction
      residual = 1 - rowSums(observed_data[!complete_idx, , drop = F])
      d[!complete_idx, 2:12] = (d[!complete_idx, 2:12] / rowSums(d[!complete_idx, 2:12])) * residual
      d[!complete_idx, 1:11] = observed_data[!complete_idx, ]
      
      return(d)
    }
  }
}


# output a list of regression_imp case arch3 matrices
# can be accessed by names
arch3_Dirichlet_Imp = lapply(
  mget(paste0("arch3mis2pct", c(10, 30, 90))), {
    function(X) {
      dirichlet_reg_imp(X, 3, MetaData)
    }
  }
)

names(arch3_Dirichlet_Imp)

# output a list of complete case arch6 matrices
# can be accessed by names
arch6_Dirichlet_Imp = lapply(
  mget(paste0("arch6mis", c(rep(2, 3), rep(4, 3)), "pct", c(10, 30, 90))), {
    function(X) {
      dirichlet_reg_imp(X, 6, MetaData)
    }
  }
)

names(arch6_Dirichlet_Imp)

# output a list of complete case arch12 matrices
# can be accessed by names
arch12_Dirichlet_Imp = lapply(
  mget(paste0("arch12mis", c(rep(2, 3), rep(4, 3), rep(11, 3)), "pct", c(10, 30, 90))), {
    function(X) {
      dirichlet_reg_imp(X, 12, MetaData)
    }
  }
)

names(arch12_Dirichlet_Imp)

## save RData
save(list = c("arch3_Dirichlet_Imp", "arch6_Dirichlet_Imp", "arch12_Dirichlet_Imp"), 
     file = "imputed_data/Dirichlet_regression_Imp_data.RData")

# load("imputed_data/Dirichlet_regression_Imp_data.RData")


#### Multiple imputation ----

# output a list of regression_imp case arch3 matrices
# can be accessed by names
arch3_Mice = lapply(
  mget(paste0("arch3mis2pct", c(10, 30, 90))), {
    function(X) {
      X_joint = cbind(X, MetaData)
      colnames(X_joint) = paste0("V", colnames(X_joint)) # handle mice error when variable starts with an integer
      mice(X_joint)
    }
  }
)

arch3_multiple_imp_data = lapply(
  mget(paste0("arch3mis2pct", c(10, 30, 90))), {
    function(X) {
      X_joint = cbind(X, MetaData)
      colnames(X_joint) = paste0("V", colnames(X_joint)) # handle mice error when variable starts with an integer
      complete(mice(X_joint))[, 1:3] / rowSums(complete(mice(X_joint))[, 1:3]) # remove the covariates and normalize
    }
  }
)

names(arch3_Mice)
names(arch3_multiple_imp_data)

# output a list of complete case arch6 matrices
# can be accessed by names
arch6_Mice = lapply(
  mget(paste0("arch6mis", c(rep(2, 3), rep(4, 3)), "pct", c(10, 30, 90))), {
    function(X) {
      X_joint = cbind(X, MetaData)
      colnames(X_joint) = paste0("V", colnames(X_joint)) # handle mice error when variable starts with an integer
      mice(X_joint)
    }
  }
)

arch6_multiple_imp_data = lapply(
  mget(paste0("arch6mis", c(rep(2, 3), rep(4, 3)), "pct", c(10, 30, 90))), {
    function(X) {
      X_joint = cbind(X, MetaData)
      colnames(X_joint) = paste0("V", colnames(X_joint)) # handle mice error when variable starts with an integer
      complete(mice(X_joint))[, 1:6] / rowSums(complete(mice(X_joint))[, 1:6]) # remove the covariates and normalize
    }
  }
)

names(arch6_Mice)
names(arch6_multiple_imp_data)

# output a list of complete case arch12 matrices
# can be accessed by names
arch12_Mice = lapply(
  mget(paste0("arch12mis", c(rep(2, 3), rep(4, 3), rep(11, 3)), "pct", c(10, 30, 90))), {
    function(X) {
      X_joint = cbind(X, MetaData)
      colnames(X_joint) = paste0("V", colnames(X_joint)) # handle mice error when variable starts with an integer
      mice(X_joint)
    }
  }
)

arch12_multiple_imp_data = lapply(
  mget(paste0("arch12mis", c(rep(2, 3), rep(4, 3), rep(11, 3)), "pct", c(10, 30, 90))), {
    function(X) {
      X_joint = cbind(X, MetaData)
      colnames(X_joint) = paste0("V", colnames(X_joint)) # handle mice error when variable starts with an integer
      complete(mice(X_joint))[, 1:12] / rowSums(complete(mice(X_joint))[, 1:12]) # remove the covariates and normalize
    }
  }
)

names(arch12_Mice)
names(arch12_multiple_imp_data)

## save RData
save(list = c("arch3_Mice", "arch6_Mice", "arch12_Mice"), 
     file = "imputed_data/MICE_objects.RData")

save(list = c("arch3_multiple_imp_data", "arch6_multiple_imp_data", "arch12_multiple_imp_data"), 
     file = "imputed_data/Multiple_Imp_data.RData")

# load("imputed_data/MICE_objects.RData")
# load("imputed_data/Multiple_Imp_data.RData")


