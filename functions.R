library(tidyverse)
library(combinat)

# utility functions
get_del_one <- function(v, s=NULL){
  if (!is.null(s)){
    perm <- map_dbl(1:10, ~which(s==.x))
    v_perm=perm_v(v,perm)
  } else {v_perm=v}
  
  j <- length(v_perm):2
  place <- map_dbl(j, ~which(v_perm==.x)-1)
  
  map2_int(j,place,~as.numeric(all(v_perm[min(1,.y):.y]>.x)))
}

perm_v <- function(v,perm){
  J=length(v)
  perm_df <- tibble(orig=1:J, new=perm)
  left_join(tibble(v=v),perm_df,by=c("v"="orig")) %>% pull(new)
}

get_del_perm <- function(perm, v_info, v_counts){
  xwalk <- select(v_info,v_str, del_str) 
  v_dat <- select(v_info, v) %>% mutate(v_count=v_counts)
  
  v_dat_perm <- v_dat %>% 
    mutate(v_perm=map(v, ~paste0(perm_v(.x,perm=perm), collapse=""))) %>%
    unnest(v_perm)
  
  v_dat_perm %>%
    left_join(xwalk, by=c("v_perm"="v_str")) %>%
    group_by(del_str) %>%
    summarise(del_count=sum(v_count)) %>%
    ungroup() %>% pull(del_count)
}

# simulate data
sim_del_counts <- function(a, N=10000, del_info){
  stopifnot(length(a)==length(del_info$del[[1]]))
  
  del_info %>%
    mutate(prob_v=map(del, ~exp(-sum(a*.x)))) %>%
    unnest(prob_v) %>%
    mutate(p_del=n_del*prob_v,
           p_del=p_del/sum(p_del),
           count=rmultinom(1,N, p_del)[,1]) %>%
    pull(count)
}

sim_v_counts <- function(a,N,v_info){
  v_info %>%
    mutate(prob_v=map(del, ~exp(-sum(a*.x)))) %>%
    unnest(prob_v) %>%
    mutate(prob_v=prob_v/sum(prob_v),
           count=rmultinom(1,N, prob_v)[,1]) %>%
    pull(count)
}

# MLE for a
make_constraint_matrix <- function(J){
  ui <- diag(J-1)
  for (i in 1:(J-2)){
    for (j in (i+1):(J-1)){
      ui[i,j] <- (-1)
    }
  }
  return(ui)
}

compute_llik <- function(a, del_counts, del_mat, num_del){
  del_a <- del_mat%*%a  
  l_denom <- log(sum(del_info$n_del * exp(-del_a)))
  l_nums <- -del_a
  llik <- sum(l_nums*del_counts) - N*l_denom
  return(llik)
}

# llik1 <- sum(-del1*del_a) - N*l_denom
# llik2 <- sum(-del2*del_a) - N*l_denom
# 
# a <- res1$a
# del_a <- del_mat%*%a  
# sum(del1*del_a)
# sum(del2*del_a)
 
# nj <- map_df(1:ncol(del_mat),
#               ~tibble(nj_1=sum(del1[del_mat[,.x]==1]),
#                       nj_2=sum(del2[del_mat[,.x]==1]))) 
# sum((nj$nj_1*a))
# sum((nj$nj_2*a))

est_MLE_a <- function(del_counts, del_mat, num_del, eps=0){
  ui <- make_constraint_matrix(J)
  ci <- rep(eps,J-1)
  
  init <- solve(ui)%*%rep(2*max(eps,.001),J-1)
  MLE_constr <- constrOptim(init,
                            function(a){-compute_llik(a,del_counts,
                                                      del_mat,num_del)},
                            grad=NULL,
                            ui=ui,ci=ci,
                            control = list(reltol = 10^{-6}))
  return(MLE_constr)
}

est_MLE_enumerate <- function(v_counts, del_mat, num_del){
  perms <- permn(J)
  del_counts <- map(perms, get_del_perm, v_info, v_counts)
  
  MLES <- map(del_counts, est_MLE_a,del_mat,num_del)
  
  neg_llik <- map_dbl(MLES, ~(.$value))
  order_opt <- perms[which.min(neg_llik)][[1]]
  a_opt <- MLES[[which.min(neg_llik)]]$par
  
  return(list(order=order_opt, a=a_opt))
}

est_MLE_greedy <- function(v_counts, del_mat, num_del,
                           eps=0, s_init=NULL){
  if (is.null(s_init)){s_init=1:J}
  s_curr <- s_init
  del_counts <- get_del_perm(perm_curr, v_info, v_counts)
  neg_llik_curr <- est_MLE_a(del_counts, del_mat, num_del, 
                             eps=eps)$value
  
  any_change=T
  while (any_change){
    any_change=F
    for (j in 2:J){
      print(s_curr)
      print(neg_llik_curr)
      # make new perm by flipping entries j and j-1
      s_prop <- s_curr
      s_prop[j] <- s[j-1]
      s_prop[j-1] <- s_curr[j]
      
      # calc new llik
      del_counts <- get_del_perm(s_prop, v_info, v_counts)
      neg_llik_prop <- est_MLE_a(del_counts, del_mat, 
                                 num_del,eps=eps)$value
      print(neg_llik_prop)
      
      if (neg_llik_prop<neg_llik_curr){
        s_curr<-s_prop
        neg_llik_curr<-neg_llik_prop
        any_change=T
      }
    }
    print(s_curr)  
    print(neg_llik_curr)
  }
  
  order_opt <- s_curr
  del_counts <- get_del_perm(order_opt, v_info, v_counts)
  a_opt <- est_MLE_a(del_counts, del_mat, num_del, eps=eps)$par
  return(list(order=order_opt, a=a_opt))
}

