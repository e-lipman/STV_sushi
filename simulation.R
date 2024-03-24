library(tidyverse)
source("functions.R")

# hyperparams
N <- 1000
a <- 2^{3:0}/2

k <- exp(-a)
J <- length(a)+1

# read in number of distinct votes by delta vector
del_info <- readRDS(file.path("delta_counts",paste0("delta_",J,".RDS")))
del_mat <- reduce(del_info$del, rbind)
num_del <- del_info$n_del

v_info <- readRDS(file.path("delta_counts",paste0("v_",J,".RDS")))

# simulate data and estimate MLE - fixed order
del_counts_sim <- sim_del_counts(a,N,del_info)

MLE <- est_MLE_a(del_counts_sim,del_mat,num_del)
res <- tibble(a=a,est=MLE$par)
ui <- make_constraint_matrix(J)
res$const_a=(ui%*%res$a)[,1]
res$const_est=(ui%*%res$est)[,1]

plot(res$a,res$est)
abline(a=0,b=1)
res

# simulate data and estimate MLE - estimate order
v_counts_sim <- sim_v_counts(a,N,v_info)

MLE <- est_MLE_greedy(v_counts_sim,del_mat,num_del,
                      perm_init=c(1,5:2))
MLE <- est_MLE_enumerate(v_counts_sim,del_mat,num_del)

res <- tibble(a=a,est=MLE$a)
ui <- make_constraint_matrix(J)
res$const_a=(ui%*%res$a)[,1]
res$const_est=(ui%*%res$est)[,1]

MLE$order
res
