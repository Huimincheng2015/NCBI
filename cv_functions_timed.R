##############################################################
# Modified BM and CV functions with separate timing tracking
##############################################################

BM_timed = function(task="par_sel", Phat=NULL, model, A, K, tau=NULL,
                    theta=NULL, par1=NULL, par2=NULL, par_method=NULL,
                    index_random=T) {

  A_ori = A
  n = dim(A)[1]
  A_up = A_ori[upper.tri(A_ori)]
  rd = floor(length(A_up)/K)

  if(index_random == T){
    index = sample(1:length(A_up), length(A_up), replace=F)
  } else {
    index = 1:length(A_up)
  }

  # Timing: CV overhead - Step 1: construct K mirror networks
  time_cv_overhead_start = Sys.time()

  lap_jj = lapply(1:K, function(iiii){
    A_up = A_ori[upper.tri(A_ori, diag = FALSE)]
    set.seed(iiii)

    # K-fold cross validation
    sap = (index)[(1+rd*(iiii-1)): min(rd*iiii, length(A_up))]
    if(length(theta) == 0){
      theta = sum(A_up[-sap])/length(A_up[-sap])
    }
    A_up[sap] = rbinom(length(sap), 1, theta)
    mat_temp = matrix(0, nrow(A), ncol(A))
    mat_temp[upper.tri(mat_temp)] = A_up
    A1 = mat_temp + t(mat_temp)
    list(A1=A1, sap=sap, theta=theta)
  })

  time_cv_overhead_fold = as.numeric(difftime(Sys.time(), time_cv_overhead_start, units="secs"))

  # Timing: Estimation and validation
  time_estimation_total = 0
  time_cv_overhead_validation = 0

  cv_score_k = lapply(lap_jj, function(x){
    theta = x$theta

    # TIME: Model estimation
    time_est_start = Sys.time()
    if(task == "par_sel"){
      res2 = model(x$A1, par1)
      phat = res2$P
      diag(phat) = 0
    } else {
      if(length(par_method) == 0){
        res2 = model(x$A1)
      } else {
        res2 = model(x$A1, par_method)
      }
      phat = res2$P
      diag(phat) = 0
    }
    time_estimation = as.numeric(difftime(Sys.time(), time_est_start, units="secs"))

    # TIME: CV overhead - validation (rescaling and score computation)
    time_val_start = Sys.time()

    sap = x$sap
    phat = (K*phat - theta)/(K-1)  # rescale
    phat[phat >= 1] = 1
    phat[phat <= 0] = 0

    # MSE on validation set
    diag(phat) = 0
    mse_sap = mean((phat[upper.tri(phat)][sap] - A_ori[upper.tri(A_ori)][sap])^2)

    time_validation = as.numeric(difftime(Sys.time(), time_val_start, units="secs"))

    list(mse=mse_sap, time_est=time_estimation, time_val=time_validation)
  })

  # Aggregate timing across folds
  time_estimation_total = sum(sapply(cv_score_k, function(x) x$time_est))
  time_cv_overhead_validation = sum(sapply(cv_score_k, function(x) x$time_val))
  time_cv_overhead_total = time_cv_overhead_fold + time_cv_overhead_validation

  mse_stat = mean(sapply(cv_score_k, function(x) x$mse))

  return(list(
    cv_score = mse_stat,
    time_estimation = time_estimation_total,
    time_cv_overhead = time_cv_overhead_total,
    time_total = time_estimation_total + time_cv_overhead_total
  ))
}


CV_timed = function(task="par_sel", Phat=NULL, model, A, K, tau=0.9,
                    par1=NULL, par2=NULL, index_random=T, rank_=NULL) {

  A_ori = A
  n = dim(A)[1]
  A_up = A[upper.tri(A)]
  rd = round(length(A_up)/K)

  if(index_random == T){
    index = sample(1:length(A_up), length(A_up))
  } else {
    index = 1:length(A_up)
  }

  # Timing: CV overhead - Step 1: construct K mirror networks with imputation
  time_cv_overhead_start = Sys.time()

  lap_jj = lapply(1:K, function(iiii){
    set.seed(iiii)
    A_up = A_ori[upper.tri(A_ori)]

    # K-fold cross validation
    sap = (index)[(1+rd*(iiii-1)): min(rd*iiii, length(index))]
    tau = (K-1)/K

    # Mask entries
    A_up[sap] = 0
    mat_temp = matrix(0, nrow(A), ncol(A))
    mat_temp[upper.tri(mat_temp)] = A_up
    A1 = mat_temp + t(mat_temp)

    # SVD imputation
    if(length(rank_) == 0){
      rank_ = rankMatrix(A1)
    }
    svd_ = svd(A1 / tau)
    A_hat = svd_$u[, 1:rank_] %*% diag(svd_$d[1:rank_]) %*% t(svd_$v[, 1:rank_])
    diag(A_hat) = 0
    A_hat[A_hat <= 0.5] = 0
    A_hat[A_hat > 0.5] = 1

    A1 = A_hat
    list(A1=A1, sap=sap)
  })

  time_cv_overhead_fold = as.numeric(difftime(Sys.time(), time_cv_overhead_start, units="secs"))

  # Timing: Estimation and validation
  time_estimation_total = 0
  time_cv_overhead_validation = 0

  cv_score_k = lapply(lap_jj, function(x){
    # TIME: Model estimation
    time_est_start = Sys.time()
    if(task == "par_sel"){
      res2 = model(x$A1, par1)
      phat = res2$P
    } else {
      res2 = model(x$A1)
      phat = res2$P
    }
    time_estimation = as.numeric(difftime(Sys.time(), time_est_start, units="secs"))

    # TIME: CV overhead - validation (score computation)
    time_val_start = Sys.time()

    sap = x$sap
    diag(phat) = 0
    mse_sap = mean((phat[upper.tri(phat)][sap] - A_ori[upper.tri(A_ori)][sap])^2)

    time_validation = as.numeric(difftime(Sys.time(), time_val_start, units="secs"))

    list(mse=mse_sap, time_est=time_estimation, time_val=time_validation)
  })

  # Aggregate timing across folds
  time_estimation_total = sum(sapply(cv_score_k, function(x) x$time_est))
  time_cv_overhead_validation = sum(sapply(cv_score_k, function(x) x$time_val))
  time_cv_overhead_total = time_cv_overhead_fold + time_cv_overhead_validation

  mse_stat = mean(sapply(cv_score_k, function(x) x$mse))

  return(list(
    cv_score = mse_stat,
    time_estimation = time_estimation_total,
    time_cv_overhead = time_cv_overhead_total,
    time_total = time_estimation_total + time_cv_overhead_total
  ))
}
