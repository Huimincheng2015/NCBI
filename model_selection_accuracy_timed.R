#####model selection accuracy - WITH SEPARATE TIMING###########
# Modified to separately track estimation time vs. CV overhead time

rm(list=ls())

source("ICE.R")
library(graphon)
library(dplyr)
library(CVXR)
library(transport)
library(igraph)
library(ggplot2)
library(reshape2)
library(plyr)
library(viridis)
library(Matrix)
library(pheatmap)
library(gridExtra)
library(pheatmap)
library(foreach)
library(doParallel)
library(reshape2)
proc=proc.time()

source("function.R")
source("network_generate.R")
source("graphon-master/R/auxiliary.R")
source("graphon-master/R/est.nbdsmooth.R")
source("cv_functions_timed.R")  # NEW: Load timing-aware CV functions

result = NULL
mean_dat = g_graphonid = g_graphonid_temp = NULL
sample_size_lis = c(50,100,150,200)
core = 9
repeatn = 5
tau_lis = c(0.7)
graphon_id_lis = c(5,16,19,4)
rank_ = 3
K_cv = 3
tau_cv = 0.8
run_ecv = FALSE

for(graphon_id in graphon_id_lis){
  sample_id_temp = 1
  for(sample_size in sample_size_lis){
    for(tau in tau_lis){
      gp = gp_generate(sample_size, graphon_id, type="seq")
      P = gp$P
      diag(P) = 0

      cat(graphon_id, sample_size, "\n")
      cl = makeCluster(core, outfile = "")
      registerDoParallel(cl)
      clusterExport(cl=cl, varlist=c("P","repeatn","sample_size","graphon_id",
                                      "rank_","K_cv","tau_cv",
                                      "EL","likeli",
                                      "run_ecv",
                                      "est.nbdsmooth2","BM_timed","CV_timed","tau","KL",
                                      "gp_generate","asSymmetric"), envir=environment())
      clusterEvalQ(cl, {
        library(graphon)
        library(matlabr)
        library(Matrix)
        library(Rfast)
        library(reshape2)
        library(ggplot2)
        library(plyr)
        library(dplyr)
        library(gridExtra)
        library(pheatmap)
        library(seewave)
        library(irlba)
        source("graphon-master/R/auxiliary.R")
        source("graphon-master/R/est.nbdsmooth.R")
        source("function.R")
        source("cv_functions_timed.R")  # NEW: Load in workers
      })

      lap = parLapply(cl, 1:repeatn, function(repeat_i){
        set.seed(repeat_i)
        gp = gp_generate(sample_size, graphon_id, type="unif")
        P = gp$P
        diag(P) = 0
        A = gmodel.P(P, rep=1, symmetric.out=TRUE)
        diag(A) = 0

        ###NB###
        model = est.nbdsmooth2

        h_lis = c(1, seq(0.1, 5, 0.2))
        A_up = A[upper.tri(A)]
        theta = sum(A_up)/length(A_up)

        mse_lis = BM_lis = EL_lis = CV_lis = kl_lis = rep(NA, length(h_lis))

        # NEW: Separate timing arrays
        time_BM_total = time_BM_estimation = time_BM_cv_overhead = rep(NA, length(h_lis))
        time_CV_total = time_CV_estimation = time_CV_cv_overhead = rep(NA, length(h_lis))

        for(i in 1:length(h_lis)){
          tryCatch({
            par = h_lis[i]
            ice_res = model(A, par)
            temp = ice_res$P
            diag(temp) = 0
            phat = temp

            #----kl------
            kl_lis[i] = KL(phat, P)
            #----mse------
            err_ICE = mean((phat - P)^2)
            mse_lis[i] = err_ICE

            #----GGCV (BM) with detailed timing------
            K = ceiling(0.2*nrow(A))
            index_random = T

            # NEW: Use timed version
            BM_result = BM_timed(task="par_sel", phat, model, A, K=K, tau=0.9,
                                 theta=NULL, par1=h_lis[i], par2=NULL,
                                 index_random=index_random)

            BM_lis[i] = BM_result$cv_score
            time_BM_total[i] = BM_result$time_total
            time_BM_estimation[i] = BM_result$time_estimation
            time_BM_cv_overhead[i] = BM_result$time_cv_overhead

            #----EL------
            EL_lis[i] = EL(phat, P, K=K)

            #----ECV with detailed timing------
            K_cv = K
            tau_cv = (K-1)/K

            if(run_ecv == TRUE){
              # NEW: Use timed version
              CV_result = CV_timed(task="par_sel", phat, model, A, K_cv, tau_cv,
                                   par1=h_lis[i], par2=NULL,
                                   index_random=index_random, rank_)

              CV_lis[i] = CV_result$cv_score
              time_CV_total[i] = CV_result$time_total
              time_CV_estimation[i] = CV_result$time_estimation
              time_CV_cv_overhead[i] = CV_result$time_cv_overhead
            } else {
              CV_lis[i] = 1
              time_CV_total[i] = 0
              time_CV_estimation[i] = 0
              time_CV_cv_overhead[i] = 0
            }

          }, error=function(e){cat("ERROR:", conditionMessage(e), "\n")})
        }

        # Aggregate times (sum across all hyperparameter values)
        time_BM_total_sum = sum(time_BM_total, na.rm=TRUE)
        time_BM_estimation_sum = sum(time_BM_estimation, na.rm=TRUE)
        time_BM_cv_overhead_sum = sum(time_BM_cv_overhead, na.rm=TRUE)

        time_CV_total_sum = sum(time_CV_total, na.rm=TRUE)
        time_CV_estimation_sum = sum(time_CV_estimation, na.rm=TRUE)
        time_CV_cv_overhead_sum = sum(time_CV_cv_overhead, na.rm=TRUE)

        BM_MSE = mse_lis[which(BM_lis %in% min(BM_lis))[length(which(BM_lis %in% min(BM_lis)))]]
        CV_MSE = mse_lis[which(CV_lis %in% min(CV_lis))[1]]

        dat = data.frame(
          graphon_id = graphon_id,
          sample_size = sample_size,
          repeat_i = repeat_i,
          tau = tau,
          h_lis = h_lis,
          mse_lis = mse_lis,
          kl_lis = kl_lis,
          BM = BM_MSE,
          CV = CV_MSE,
          NAIVE = mse_lis[which.max(EL_lis)],
          NB = mse_lis[1],
          BM_index = BM_lis,
          CV_index = CV_lis,

          # NEW: Detailed timing columns
          time_BM_total = time_BM_total_sum,
          time_BM_estimation = time_BM_estimation_sum,
          time_BM_cv_overhead = time_BM_cv_overhead_sum,
          time_CV_total = time_CV_total_sum,
          time_CV_estimation = time_CV_estimation_sum,
          time_CV_cv_overhead = time_CV_cv_overhead_sum
        )
        dat = na.omit(dat)
        dat
      })

      stopCluster(cl)
      lap2 = do.call("rbind", lap)

      # Create visualization for parameter selection
      dap = ddply(lap2, .(h_lis), function(x){
        data.frame(
          mean = c(median(x$BM_index), median(x$CV_index),
                   median(x$mse_lis), median(x$kl_lis)),
          sd = c(sd(x$BM_index), sd(x$CV_index), sd(x$mse_lis), sd(x$kl_lis)),
          method = c("GGCV","ECV","MSE","KL")
        )
      })

      xx = dap[dap$method %in% "GGCV","mean"]
      dap[dap$method %in% "GGCV","mean"] = (xx - min(xx))/(max(xx) - min(xx))
      xx = dap[dap$method %in% "ECV","mean"]
      dap[dap$method %in% "ECV","mean"] = (xx - min(xx))/(max(xx) - min(xx))
      xx = dap[dap$method %in% "MSE","mean"]
      dap[dap$method %in% "MSE","mean"] = (xx - min(xx))/(max(xx) - min(xx))
      xx = dap[dap$method %in% "KL","mean"]
      dap[dap$method %in% "KL","mean"] = (xx - min(xx))/(max(xx) - min(xx))

      dap$method = factor(dap$method, levels=c("GGCV","KL","ECV","MSE"))
      dap = filter(dap, method %in% c("GGCV","MSE"))

      dap$mean = round(dap$mean, 3)
      g1 = ggplot(dap, aes(x=h_lis, y=mean, group=method, color=method, shape=method)) +
        geom_line(size=0.5) +
        geom_point(size=1) +
        theme_bw() +
        scale_color_manual(values=c("salmon","black")) +
        labs(x="Neighborhood Size M", y="Normalized Score")

      pdf(paste0("../box_results_ns/h",graphon_id,"_",sample_size,".pdf"))
      print(g1)
      dev.off()

      g_graphonid_temp[[1]][[sample_id_temp]] = g1
      sample_id_temp = sample_id_temp + 1

      result = rbind(result, lap2)
      write.csv(result, "par_sel_ns_timed.csv")

      dat0 = result[result$sample_size %in% sample_size,]
      dat0 = dat0[dat0$graphon_id %in% graphon_id,]
      dat_our = ddply(dat0, .(repeat_i), function(x){mean(x$BM)})
      dat_ECV = ddply(dat0, .(repeat_i), function(x){mean(x$CV)})
      print(median(dat_our$V1))
      print(median(dat_ECV$V1))
    }
  }

  if(length(sample_size_lis) == 4){
    pdf(paste0("../box_results_ns/h_paper",graphon_id,".pdf"), width=14, height=10)
    grid.arrange(g_graphonid_temp[[1]][[1]], g_graphonid_temp[[1]][[2]],
                 g_graphonid_temp[[1]][[3]], g_graphonid_temp[[1]][[4]],
                 nrow = 6, ncol=4)
    dev.off()
  }
  g_graphonid = c(g_graphonid, g_graphonid_temp)
}

pdf(paste0("../box_results_ns/h_paper",".pdf"), width=12, height=12)
grid.arrange(g_graphonid[[1]][[1]], g_graphonid[[1]][[2]], g_graphonid[[1]][[3]], g_graphonid[[1]][[4]],
             g_graphonid[[2]][[1]], g_graphonid[[2]][[2]], g_graphonid[[2]][[3]], g_graphonid[[2]][[4]],
             g_graphonid[[3]][[1]], g_graphonid[[3]][[2]], g_graphonid[[3]][[3]], g_graphonid[[3]][[4]],
             g_graphonid[[4]][[1]], g_graphonid[[4]][[2]], g_graphonid[[4]][[3]], g_graphonid[[4]][[4]],
             nrow = 6, ncol=4)
dev.off()


##### Accuracy #################
dat = read.csv("par_sel_ns_timed.csv")
library(dplyr)
library(plyr)
data_lis = NULL
graphon_id_lis = c(5,16,19,4)

for(graphon_id in graphon_id_lis){
  for(sample_size in 200){
    dat0 = dat[dat$sample_size %in% sample_size,]
    dat0 = dat0[dat0$graphon_id %in% graphon_id,]
    dat_default = dat0[dat0$h_lis %in% 1,]

    dat_our = ddply(dat0, .(repeat_i), function(x){
      BM_index = x$BM_index
      BM_index = round(BM_index, 3)
      wh = which(BM_index %in% min(BM_index))
      x$mse_lis[wh[length(wh)]]
    })
    dat_ECV = ddply(dat0, .(repeat_i), function(x){
      CV_index = x$CV_index
      CV_index = round(x$CV_index, 3)
      wh = which(CV_index %in% min(CV_index))
      x$mse_lis[wh[1]]
    })

    cat(graphon_id, sample_size, "\n")
    print(round(c(mean(dat_our$V1), sd(dat_our$V1))*100, 2))
    print(round(c(mean(dat_ECV$V1), sd(dat_ECV$V1))*100, 2))
    print(round(c(mean(dat_default$mse_lis), sd(dat_default$mse_lis))*100, 2))

    data_lis = rbind(data_lis, data.frame(graphon_id = graphon_id, sample_size=sample_size,
                                          mean = mean(dat_our$V1), sd = sd(dat_our$V1), method= "BM"))
    data_lis = rbind(data_lis, data.frame(graphon_id = graphon_id, sample_size=sample_size,
                                          mean = mean(dat_ECV$V1), sd = sd(dat_ECV$V1), method= "CV"))
  }
}

g1_lis = NULL
for(i in 1:length(graphon_id_lis)){
  graphon_id1 = graphon_id_lis[i]
  dap = dplyr::filter(data_lis, graphon_id %in% graphon_id1)
  g1 = ggplot(dap, aes(x=sample_size, y = mean, fill=method, color = method, shape=method)) +
    geom_line(size=1) +
    geom_point(size=3) +
    theme_bw() +
    theme(legend.position="none") +
    labs(x="n", y="MSE") +
    scale_color_manual(values=c("salmon","#AD84C1")) +
    scale_fill_manual(values=c("salmon","#AD84C1"))
  g1_lis[[i]] = g1
}
pdf(paste0("../box_results_ns/accuracy",".pdf"))
grid.arrange(g1_lis[[1]], g1_lis[[2]],
             g1_lis[[3]], g1_lis[[4]],
             nrow = 3, ncol=2)
dev.off()


#####TIME - TOTAL (original)#################
dat = read.csv("par_sel_ns_timed.csv")
library(dplyr)
library(plyr)
data_lis = NULL

for(i in 1:length(graphon_id_lis)){
  graphon_id = graphon_id_lis[i]
  for(sample_size in sample_size_lis){
    dat0 = dat[dat$sample_size %in% sample_size,]
    dat0 = dat0[dat0$graphon_id %in% graphon_id,]

    # Use new timing columns
    dat_our = ddply(dat0, .(repeat_i), function(x){unique(x$time_BM_total)})
    dat_ECV = ddply(dat0, .(repeat_i), function(x){unique(x$time_CV_total)})

    data_lis = rbind(data_lis, data.frame(graphon_id = graphon_id, sample_size=sample_size,
                                          mean = mean(dat_our$V1), sd = sd(dat_our$V1), method= "BM"))
    data_lis = rbind(data_lis, data.frame(graphon_id = graphon_id, sample_size=sample_size,
                                          mean = mean(dat_ECV$V1), sd = sd(dat_ECV$V1), method= "CV"))
  }
}

g1_lis = NULL
for(i in 1:length(graphon_id_lis)){
  graphon_id1 = graphon_id_lis[i]
  dap = dplyr::filter(data_lis, graphon_id %in% graphon_id1)
  g1 = ggplot(dap, aes(x=sample_size, y = mean, fill=method, color = method, shape=method)) +
    geom_line(size=1) +
    geom_point(size=3) +
    theme_bw() +
    theme(legend.position="none") +
    labs(x="n", y="Total Time (secs)") +
    scale_color_manual(values=c("salmon","#AD84C1")) +
    scale_fill_manual(values=c("salmon","#AD84C1"))
  g1_lis[[i]] = g1
}
pdf(paste0("../time_results/time_NS_total",".pdf"), width=9)
grid.arrange(g1_lis[[1]], g1_lis[[2]],
             g1_lis[[3]], g1_lis[[4]],
             nrow = 4, ncol=4)
dev.off()


##### NEW: TIME - CV OVERHEAD ONLY #################
dat = read.csv("par_sel_ns_timed.csv")
library(dplyr)
library(plyr)
data_lis_cv_only = NULL

for(i in 1:length(graphon_id_lis)){
  graphon_id = graphon_id_lis[i]
  for(sample_size in sample_size_lis){
    dat0 = dat[dat$sample_size %in% sample_size,]
    dat0 = dat0[dat0$graphon_id %in% graphon_id,]

    # Extract CV overhead time only
    dat_our = ddply(dat0, .(repeat_i), function(x){unique(x$time_BM_cv_overhead)})
    dat_ECV = ddply(dat0, .(repeat_i), function(x){unique(x$time_CV_cv_overhead)})

    data_lis_cv_only = rbind(data_lis_cv_only,
                             data.frame(graphon_id = graphon_id, sample_size=sample_size,
                                       mean = mean(dat_our$V1), sd = sd(dat_our$V1),
                                       method= "GGCV"))
    data_lis_cv_only = rbind(data_lis_cv_only,
                             data.frame(graphon_id = graphon_id, sample_size=sample_size,
                                       mean = mean(dat_ECV$V1), sd = sd(dat_ECV$V1),
                                       method= "ECV"))
  }
}

g1_lis_cv_only = NULL
for(i in 1:length(graphon_id_lis)){
  graphon_id1 = graphon_id_lis[i]
  dap = dplyr::filter(data_lis_cv_only, graphon_id %in% graphon_id1)
  g1 = ggplot(dap, aes(x=sample_size, y = mean, fill=method, color = method, shape=method)) +
    geom_line(size=1) +
    geom_point(size=3) +
    theme_bw() +
    theme(legend.position="none") +
    labs(x="n", y="CV Overhead Time (secs)") +
    scale_color_manual(values=c("salmon","#AD84C1")) +
    scale_fill_manual(values=c("salmon","#AD84C1"))
  g1_lis_cv_only[[i]] = g1
}

pdf(paste0("../time_results/time_NS_cv_overhead_only",".pdf"), width=9, height=6)
grid.arrange(g1_lis_cv_only[[1]], g1_lis_cv_only[[2]],
             g1_lis_cv_only[[3]], g1_lis_cv_only[[4]],
             nrow = 2, ncol=2)
dev.off()


##### NEW: TIME BREAKDOWN ANALYSIS #################
# Create detailed breakdown showing estimation vs CV overhead

dat = read.csv("par_sel_ns_timed.csv")
timing_breakdown = NULL

for(graphon_id in graphon_id_lis){
  for(sample_size in sample_size_lis){
    dat0 = dat[dat$sample_size %in% sample_size,]
    dat0 = dat0[dat0$graphon_id %in% graphon_id,]

    # BM method breakdown
    timing_breakdown = rbind(timing_breakdown,
                            data.frame(
                              graphon_id = graphon_id,
                              sample_size = sample_size,
                              method = "GGCV",
                              time_type = "Estimation",
                              mean_time = mean(unique(dat0$time_BM_estimation)),
                              sd_time = sd(unique(dat0$time_BM_estimation))
                            ))
    timing_breakdown = rbind(timing_breakdown,
                            data.frame(
                              graphon_id = graphon_id,
                              sample_size = sample_size,
                              method = "GGCV",
                              time_type = "CV Overhead",
                              mean_time = mean(unique(dat0$time_BM_cv_overhead)),
                              sd_time = sd(unique(dat0$time_BM_cv_overhead))
                            ))

    # ECV method breakdown
    timing_breakdown = rbind(timing_breakdown,
                            data.frame(
                              graphon_id = graphon_id,
                              sample_size = sample_size,
                              method = "ECV",
                              time_type = "Estimation",
                              mean_time = mean(unique(dat0$time_CV_estimation)),
                              sd_time = sd(unique(dat0$time_CV_estimation))
                            ))
    timing_breakdown = rbind(timing_breakdown,
                            data.frame(
                              graphon_id = graphon_id,
                              sample_size = sample_size,
                              method = "ECV",
                              time_type = "CV Overhead",
                              mean_time = mean(unique(dat0$time_CV_cv_overhead)),
                              sd_time = sd(unique(dat0$time_CV_cv_overhead))
                            ))
  }
}

# Save detailed timing breakdown
write.csv(timing_breakdown, "../time_results/timing_breakdown_detailed.csv", row.names=FALSE)

# Create stacked bar plot showing breakdown
timing_breakdown$time_type = factor(timing_breakdown$time_type,
                                   levels = c("Estimation", "CV Overhead"))

g_breakdown_lis = NULL
for(i in 1:length(graphon_id_lis)){
  graphon_id1 = graphon_id_lis[i]
  dap = dplyr::filter(timing_breakdown, graphon_id %in% graphon_id1)

  g1 = ggplot(dap, aes(x=factor(sample_size), y=mean_time, fill=time_type)) +
    geom_bar(stat="identity", position="stack") +
    facet_wrap(~method, ncol=2) +
    theme_bw() +
    labs(x="Sample Size (n)", y="Time (secs)", fill="Time Component",
         title=paste("Graphon", graphon_id1)) +
    scale_fill_manual(values=c("Estimation"="#4DAF4A", "CV Overhead"="#E41A1C"))

  g_breakdown_lis[[i]] = g1
}

pdf(paste0("../time_results/time_breakdown_stacked",".pdf"), width=10, height=8)
grid.arrange(g_breakdown_lis[[1]], g_breakdown_lis[[2]],
             g_breakdown_lis[[3]], g_breakdown_lis[[4]],
             nrow = 2, ncol=2)
dev.off()

cat("\n=== TIMING ANALYSIS COMPLETE ===\n")
cat("Results saved:\n")
cat("  - CV overhead only plot: ../time_results/time_NS_cv_overhead_only.pdf\n")
cat("  - Total time plot: ../time_results/time_NS_total.pdf\n")
cat("  - Timing breakdown (stacked): ../time_results/time_breakdown_stacked.pdf\n")
cat("  - Detailed timing data: ../time_results/timing_breakdown_detailed.csv\n")
