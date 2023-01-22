

# packages
{
  require(drc)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
}

# functions
{

  #Tanimoto matrix for existing compounds with existing compounds
  tanimoto <- function(x, similarity=FALSE) {
    n = dim(x)[1] # nrow x
    k = dim(x)[2]
    x = matrix(as.logical(x), n)
    res = apply(x,2,function(x1) {i=colSums(x1 & x) / colSums(x1 | x); ifelse(is.na(i), 1, i)})
    res = new("dspMatrix", Dim=c(k,k), x = res[upper.tri(res,diag=TRUE)],uplo="U")
    if(similarity) return(res)
    else return(1-res)
  }

  #cholesky inverse function, quicker than chol2inv(col)
  cholinv <- function(A) {
    ## invert a +ve definite matrix via pivoted Cholesky with diagonal
    ## pre-conditioning.
    d <- 1/sqrt(diag(A))
    R <- chol(d*t(A*d),pivot=TRUE)
    ipiv <- piv <- attr(R,"pivot")
    ipiv[piv] <- 1:ncol(A)
    d*t(chol2inv(R)[ipiv,ipiv]*d)
  } # cholinv
  Woodbury_determinant <- function(A,B,C){
    Ua2 <- chol(A) # A=t(U)xU, take U from this
    Ub2 <- chol(B)
    A_inv <- chol2inv(Ua2)# A^-1
    B_inv <- chol2inv(Ub2) # B^-1
    Ua2 <- chol(A_inv) # A^-1=U^tU
    UC <- Ua2 %*% C # UxC
    CA <- crossprod(UC,Ua2) # t(UxC)xU = t(C)xt(U)xU = t(C)xA^-1
    BpCAC <- B_inv + crossprod(UC) # C^TxA^-1xC
    BpCAC_ch = chol(BpCAC)
    BpCAC_inv <- chol2inv(BpCAC_ch) # (B^-1 + t(C)xAxC)^-1
    logdetA <- -2*sum(log(diag(Ua2)))
    logdetB <- -2*sum(log(diag(Ub2)))
    logdetbr <- 2*sum(log(diag(BpCAC_ch)))
    logdetin <- logdetbr+logdetB+logdetA # final determinant calculation
    V <- chol(BpCAC_inv)
    VCA <- V %*% CA
    ACYCA <- crossprod(VCA)
    #ACYCA <- as(ACYCA, "dspMatrix")
    X <-  A_inv- ACYCA
    list(X,logdetin)
  }
  likelihood_eval <- function(K,lambda,n,P,H,Ky_inv,y){ 
    diag(K) <- diag(K) + 1e-6# correlation matrix for compound
    a <- Diagonal(n,lambda) # diagonal matrix with variance terms
    R <- Woodbury_determinant(a,K,P)
    Ky_inv <- as.matrix(R[[1]])
    logdetKy <- R[[2]]
    cp <- crossprod(H,Ky_inv)
    b1 <- solve(cp%*%H)
    b2 <- cp%*%y
    beta_hat <- as.double(b1%*%b2)
    mu <- rowSums(t(beta_hat*t(H)))
    inner <- as.double(t(as.matrix((y-mu)))%*%Ky_inv%*%as.matrix(y-mu))
    sgsq <- inner/n
    ll <- - 0.5*n*log(sgsq) - 0.5*logdetKy - 0.5*inner/sgsq
    -ll
  }
  ltan <- function(par,H,y,n,dmc,P){
    lambda <- exp(par[1])
    K <- 1 - dmc # correlation matrix for compound
    likelihood_eval(K,lambda,n,P,H,Ky_inv,y)
  }
  lme <- function(par,H,y,n,dmc,P){
    lambda <- exp(par[1])
    K <- diag(1,nrow(dmc))
    likelihood_eval(K,lambda,n,P,H,Ky_inv,y)
  }
  lgau <- function(par,H,y,n,dmc,P){
    lambda <- exp(par[1])
    l <- exp(par[2])
    K <- exp(-dmc/l)
    likelihood_eval(K,lambda,n,P,H,Ky_inv,y)
  }
  lexp <- function(par,H,y,n,dmc,P){
    lambda <- exp(par[1])
    l <- exp(par[2])
    K <- exp(-(dmc/l)^0.5)
    likelihood_eval(K,lambda,n,P,H,Ky_inv,y)
  }
  
}

options(warn=-1)

#### Read and process fingerprint data ####

fp <- read.csv('compounds.csv',header = T)

compound_levels <- sort(fp$C)

fp$C <- factor(fp$C, levels = compound_levels)

fp <- fp[order(fp$C),] # reorder according to factor levels

#### Read and process fictitious data ####

dat <- read.csv('data.csv',header = T)

n <- nrow(dat)

factors <- c("c1", "c2", "c3", "c4")
dat[factors] <- lapply(dat[factors], factor) # make categorical
dat$C <- factor(dat$C, levels = compound_levels) # define factor levels of the compounds

dat$log_n1 <- log(dat$n1/max(dat$n1)) # take the log of the scaled rates
dat$yt <- qlogis(pmin(99, pmax(1, dat$y))/100) # logit transform the damages

Fixed <- ~ log_n1 + c1 + c2 + c3 + c4 # regressor variables
Random <- ~ C - 1 # GP without intercept

gp_regression <- function(Fixed, Random, 
                          fp, dat, 
                          gpcov,
                          tts, params,
                          returnData = TRUE){
  
  #' fp: matrix, fingerprints with label column
  #' dat: matrix, data 
  #' gpcov: string,covariance function for GP, choice of c("mixed","tanimoto","exponential","gaussian") 
  #' tts: float between 0 and 1, fraction for train/test split, e.g., tts = 0.8 is 80% training data 20% test.
  #' returnData: logical, return data used for model fitting/testing
  
  out <- list()
  out$Model <- out$Data <- out$Fit <- list()  
  
  #### Train test split ####
  
  set.seed(1)
  train_prop <- tts
  train_samp <- sample(1:n,train_prop*n)
  test_samp <- (1:n)[-train_samp]
  
  train_set <- dat[train_samp,]
  test_set <- dat[test_samp,]
  
  y_train <- train_set$yt
  y_test <- test_set$yt
  
  #### Design matrices for regressors and random effects ####
  
  H <- model.matrix(Fixed, data = train_set) # regressor design matrix for training set
  H_star <- model.matrix(Fixed, data = test_set) # regressor design matrix for test set
  
  compound_train <- model.matrix(Random, data = train_set) # indicator matrix for compounds in train set
  compound_test <- model.matrix(Random, data = test_set) # indicator matrix for compounds in test set
  
  distance_compound <- as.matrix(dist(fp[!names(fp)%in%labels(terms(Random))], method = "binary")) # distance matrix based on Tanimoto (Jaccard) metric
  
  if (gpcov=="mixed"){
    ll <- lme
  }
  if (gpcov=="tanimoto"){
    ll <- ltan
  }
  if (gpcov=="exponential"){
    ll <- lexp
  }
  if (gpcov=="gaussian"){
    ll <- lgau
  }
  
  #### optimise the profile-likelihood ####
  
  op_time <- system.time({op_par <- optim(params, # starting values, log(variance) and log(scale) for exponential and Gaussian K
                                          ll,
                                          H = H,
                                          y = y_train,
                                          n = nrow(H),
                                          dmc = distance_compound,
                                          P = compound_train,
                                          hessian = T)})
  
  run_time = as.numeric(op_time[3]/60) # optimisation time in minutes
  
  ssq = round(exp(op_par$par[1]), 3) # GP variance
  
  if (!gpcov %in% c("mixed","tanimoto")){
    l = round(exp(op_par$par[2]), 3)
  }
  out$Fit$Startingvals <- params
  out$Fit$Value <- op_par$value
  out$Fit$Convergence <- op_par$convergence
  out$Fit$Hessian <- out$Fit$Hessian
  
  if (gpcov=="mixed"){
    K <- diag(1,nrow(distance_compound))
    
  }
  if (gpcov=="tanimoto"){
    K <- 1 - distance_compound
  }
  if (gpcov=="exponential"){
    K <- exp(-(distance_compound/l)^0.5)
    
  }
  if (gpcov=="gaussian"){
    K <- exp(-distance_compound/l)
  }
  
  ss <- Diagonal(nrow(train_set), ssq) # diagonal matrix with variance terms
  
  R <- Woodbury_determinant(ss, K, compound_train) # get inverse and determinant using lemmas
  logdetKy <- R[[2]] # log of determinant of K_y
  
  #### covariance matrices####
  
  Ky_inv  <- R[[1]] # inverse of K_y, where K_y = (K + ssq*I)
  K_star <- compound_train%*%tcrossprod(K, compound_test) # K(X,X*)
  K_starstar <- compound_test%*%tcrossprod(K, compound_test) # K(X*,X*)
  
  
  beta_hat <- solve(crossprod(H, Ky_inv)%*%H)%*%crossprod(H, Ky_inv)%*%y_train # regressor effects
  
  #### Predictive mean ####
  
  mu <- colSums(beta_hat*t(H)) # predicted value based on regressors variables for training set
  mu_star <- colSums(beta_hat*t(H_star)) # predicted value based on regressors variables for test set
  
  mu_gstar <-  mu_star + crossprod(K_star,Ky_inv)%*%(y_train-mu) # predicted value factoring spatial effect for test set
  
  #### Predictive variance ####
  
  r <- t(H_star) - t(H)%*%Ky_inv%*%K_star # variance component
  
  cov_fstar <- K_starstar - t(K_star)%*%Ky_inv%*%K_star # conditional predictive variance
  
  var_gstar <- cov_fstar + t(r)%*%solve(t(H)%*%Ky_inv%*%H)%*%r # variances due to GP and regressor variables
  
  #### Parameters ####
  
  error <- as.double(t(y_train-mu)%*%Ky_inv%*%(y_train-mu))/length(y_train) # noise variance
  
  se_op_par <- sqrt(diag(solve(op_par$hessian))) # GP variance
  
  var_beta_hat <- matrix(error, length(beta_hat), 1)*solve(crossprod(H, Ky_inv)%*%H) # var-cov matrix of ols estimators
  se_beta_hat <- sqrt(diag(var_beta_hat)) # standard errors of ols estimators
  
  parameters <- data.frame(Parameters = c("Variance", "Error", colnames(H)),
                           Estimate = round(c(ssq, 0, as.numeric(beta_hat)), 3),
                           SE = round(c(se_op_par[1], error, se_beta_hat), 3))
  
  if (!gpcov %in% c("mixed", "tanimoto")){
    parameters <- rbind(c("Scale",round(c(l, se_op_par[2]),3)), parameters)
  }
  
  out$Parameters <- parameters
  
  #### Predictions summary ####
  
  LC <- mu_gstar - 1.96*sqrt(diag(var_gstar)) # lower 95% CI for damage predictions
  UC <- mu_gstar + 1.96*sqrt(diag(var_gstar)) # Upper 95% CI for damage predictions
  
  CI <- data.frame(Lower_CI = as.numeric(LC),
                   y_hat = as.numeric(mu_gstar),
                   Upper_CI = as.numeric(UC),
                   y_true = y_test) %>% round(2)
  
  #### Untransformed predictions summary ####
  
  tmu_gstar <- 100/(1 + exp(-mu_gstar))
  tLC <- 100/(1 + exp(-LC))
  tUC <- 100/(1 + exp(-UC))
  
  tCI <- data.frame(Lower_CI = as.numeric(tLC),
                    y_hat = as.numeric(tmu_gstar),
                    Upper_CI = as.numeric(tUC),
                    y_true = dat$y[test_samp]) %>% round(2)
  
  
  #### Model metrics ####
  
  RMSE <- sqrt(mean((tCI$y_hat-tCI$y_true)^2)) # root mean squared error
  ABSE <- mean(abs(tCI$y_hat-tCI$y_true)) # mean absolute error
  AIC <- 2*(length(op_par$par) - log(op_par$value)) # AIC
  
  Score <- mean(dnorm(tCI$y_true, # observed damages in test set
                      tCI$y_hat, # fitted values in test set
                      as.numeric(sqrt(diag(var_gstar)))*sqrt(length(tCI$y_true)), # standard deviations of fitted values, SD = SE*sqrt(n)
                      log=TRUE)) # mean score for model, accounts for uncertainty of predictions
  
  
  model_metrics <- data.frame(RMSE = RMSE, ABSE = ABSE, AIC = AIC, Score = Score, op_min = run_time) %>% round(3)
  
  #### outputs ####
  
  out$Metrics <- model_metrics
  
  if (isTRUE(returnData)) {
    out$Data$ResponseTrain <- y_train
    out$Data$Predictions <- CI
    out$Data$Predictions_transformed <- tCI
    out$Data$DesignFixedTrain <- H
    out$Data$DesignFixedTest <- H_star
    out$Data$DesignGPTrain <- compound_train
    out$Data$DesignGPTest <- compound_test
    out$Data$Fingerprints <- fp
    out$Data$DistMat <- distance_compound
    
  } else {
    out$Data$ResponseTrain <- NULL
    out$Data$ResponseTrain <- NULL
    out$Data$Predictions <- NULL
    out$Data$Predictions_transformed <- NULL
    out$Data$DesignFixedTrain <- NULL
    out$Data$DesignFixedTest <- NULL
    out$Data$DesignGPTrain <- NULL
    out$Data$DesignGPTest <- NULL
    out$Data$Fingerprints <- NULL
    out$Data$DistMat <- NULL
  }
  
  out$Model$formula_fixed <- Fixed
  out$Model$formula_random <- Random
  out$Model$fcor <- gpcov
  
  return(out)
  
}

covs <- c("mixed", "tanimoto", "exponential", "gaussian")


#### Single Model Fitting ####

tts <- 0.8 #  80% training data 
gpmix <- gp_regression(Fixed, Random, fp, dat, "mixed", tts, -1)
gptan <- gp_regression(Fixed, Random, fp, dat, "tanimoto", tts, -1)
gpexp <- gp_regression(Fixed, Random, fp, dat, "exponential", tts, c(-1,-1))
gpgau <- gp_regression(Fixed, Random, fp, dat, "gaussian", tts, c(-1,-1))

df_metric <- rbind(gpmix$Metrics,
                   gptan$Metrics,
                   gpexp$Metrics,
                   gpgau$Metrics)

rownames(df_metric) <- covs

#### Cross Validation ####

nfolds <- 2 # number of folds for cross validation
models <- expand.grid(1:nfolds, covs)
names(models) <- c("run", "K")

train_prop <- 0.8
model_metrics <- data.frame()
samples <- matrix(, nrow(dat), nfolds)

for (mod in 1:nrow(models)){  # cross validation
  
  print(paste(mod, 'of', nrow(models)))
  
  set.seed(models[mod, 1])
  inds <- sample(1:nrow(dat), nrow(dat))
  samples[, models[mod, 1]] <- inds
  dat <- dat[inds,]
  
  if (models[mod, 2] %in% c("mixed", "tanimoto")) pars <- -1 else pars <- c(-1, -1)
  
  gpmod <- gp_regression(Fixed, Random,
                         fp, dat,
                         models[mod, 2],
                         train_prop, pars,
                         FALSE)
  
  model_metrics <- rbind(cbind(gpmod$Model$fcor, gpmod$Metrics), model_metrics)

}

metrics_grouped = rename(model_metrics, K = 1) %>%
  group_by(K) %>%
  summarise(across(everything(), mean)) %>%
  data.frame()

metrics_grouped
