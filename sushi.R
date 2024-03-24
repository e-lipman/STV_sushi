library(tidyverse)
library(STV)
library(xtable)

source("functions.R")

# read sushi data
item_names <- c("shrimp","sea eel","tuna","squid",
                "sea urchin","salmon roe","egg","fatty tuna",
                "tuna roll","cucumber roll")
items <- tibble(item=1:10, name=item_names)

dat_mat <- read_delim(file.path("data","sushi3a.5000.10.order"), 
                      skip=1, col_names=F) %>% .[,3:12] + 1

dat <- dat_mat %>%
  set_names(1:10) %>%
  mutate(id=1:n()) %>%
  pivot_longer(-id, names_to="rank", values_to="item") %>%
  mutate(rank=as.numeric(rank)) %>%
  left_join(items, by="item") 

# STV
do_stv_rank <- function(dat){
  M <- length(unique(dat$item))
  ranking <- rep(NA, M)
  d <- dat
  
  for (i in 1:M){
    counts <- filter(d, rank==1) %>% 
      count(item) %>%
      arrange(n) 
    print(counts)
      
    min_item <- pull(counts, item) %>% .[1]
    ranking[i] <- min_item
    
    d <- filter(d, item!=min_item) %>%
      group_by(id) %>%
      mutate(rank=order(rank)) %>%
      ungroup()
  }
  return(ranking[M:1])
}

STV_rank <- do_stv_rank(dat)
STV_rank

# estimation
N <- nrow(dat_mat)
J <- max(dat_mat)
stopifnot(min(dat_mat)==1)

del_info <- readRDS(file.path("delta_counts",paste0("delta_",J,".RDS")))
del_mat <- reduce(del_info$del, rbind)
num_del <- del_info$n_del

make_del_counts <- function(s){
  v_dat <- 
    tibble(v=map(1:nrow(dat_mat), ~(unlist(dat_mat[.x,])))) %>%
    mutate(del=map(v, get_del_one, s=s),
           del_str=map(del, paste0, collapse="")) %>%
    unnest(del_str)
  
  del_counts <- left_join(select(del_info, del_str), 
                          count(v_dat, del_str),
                          by="del_str")  %>%
    mutate(n=ifelse(is.na(n),0,n)) %>%
    pull(n)  
  
}


get_MLE_greedy <- function(dat_mat, del_mat, num_del,
                           s_init=STV_rank, eps=0){
  if (is.null(s_init)){s_init=1:J}
  s_curr <- s_init
  del_counts <- make_del_counts(s_init)
  neg_llik_curr <- est_MLE_a(del_counts, del_mat, num_del,eps=eps)$value
  
  any_change=T
  while (any_change){
    any_change=F
    for (j in 2:J){
      print(s_curr)
      print(neg_llik_curr)
      # make new perm by flipping entries j and j-1
      s_prop <- s_curr
      s_prop[j] <- s_curr[j-1]
      s_prop[j-1] <- s_curr[j]
      
      # calc new llik
      del_counts <- make_del_counts(s_prop)
      neg_llik_prop <- est_MLE_a(del_counts, del_mat, num_del,eps=eps)$value
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
  del_counts <- make_del_counts(s_curr)
  a_opt <- est_MLE_a(del_counts, del_mat, num_del,eps=eps)$par
  return(list(order=order_opt, a=a_opt))
}

# MLE for s fixing a restrictions to have STV as MLE
v_dat <- 
  tibble(v=map(1:nrow(dat_mat), ~(unlist(dat_mat[.x,])))) %>%
  mutate(del=map(v, get_del_one, s=STV_rank),
         del_str=map(del, paste0, collapse="")) %>%
  unnest(del_str)

del_counts <- left_join(select(del_info, del_str), 
                        count(v_dat, del_str),
                        by="del_str")  %>%
  mutate(n=ifelse(is.na(n),0,n)) %>%
  pull(n) 

MLE_STV <- est_MLE_a(del_counts, del_mat, num_del, eps=0)

ui <- make_constraint_matrix(J)
res1 <- tibble(a=MLE_STV$par, k=exp(-a))
res1$a_margin=(ui%*%res1$a)[,1]
res1$k_margin=exp(-(ui%*%res1$a)[,1])

tab1 <- bind_rows(res1, tibble(a=0, k=1)) %>%
  arrange(a) %>%
  mutate(item=item_names[STV_rank])  %>%
  select(item, a, k, a_margin, k_margin)
MLE_STV$convergence
tab1

# MLE for s over all sets of a restrictions
rerun=F
eps=0
filename <- paste0("sushi_MLE_overall_eps",eps,".RDS")
if (rerun){
  MLE_overall <- get_MLE_greedy(dat_mat, del_mat, num_del, s_init = STV_rank)
  saveRDS(MLE_overall,file.path("results",filename))
} else {
  MLE_overall <- readRDS(file.path("results",filename))  
}

ui <- make_constraint_matrix(J)
res2 <- tibble(a=MLE_overall$a, k=exp(-a))
res2$a_margin=(ui%*%res2$a)[,1]
res2$k_margin=exp(-(ui%*%res2$a)[,1])

tab2 <- bind_rows(res2, tibble(a=0, k=1)) %>%
  arrange(a) %>%
  mutate(item=item_names[MLE_overall$order])  %>%
  select(item, a, k, a_margin, k_margin)
tab2

xtable(tab1[,c(1,2,3,5)], digits=4) %>% print(include.rownames=F)

# chk
del1 <- make_del_counts(STV_rank)
del2 <- make_del_counts(MLE_overall$order)

MLE1 <- est_MLE_a(del1, del_mat, num_del, eps=0)
MLE2 <- est_MLE_a(del2, del_mat, num_del, eps=0)

MLE1$convergence
MLE2$convergence

a1 <- MLE1$par
a2 <- MLE2$par

## want to confirm that del2 actually does better than del1 
## (min over a1 and a2). If not it is issues with optimizer

compute_llik(a1, del1, del_mat, num_del)
compute_llik(a2, del2, del_mat, num_del)

# max del for a1 (del2)
compute_llik(a1, del1, del_mat, num_del)
compute_llik(a1, del2, del_mat, num_del)

# max del for a2 (del2)
compute_llik(a2, del1, del_mat, num_del)
compute_llik(a2, del2, del_mat, num_del)

# max a for del1 (a2)
compute_llik(a1, del1, del_mat, num_del) 
compute_llik(a2, del1, del_mat, num_del) 

# max a for del2 (a2)
compute_llik(a1, del2, del_mat, num_del) 
compute_llik(a2, del2, del_mat, num_del)


# Goodness of fit

## expected vs observed delta
counts <- 
  del_info %>%
  mutate(prob_v=map(del, ~exp(-sum(a1*.x))),
         first1=map(del, ~(which(cumsum(.x)==1))[1]),
         first1=map(first1,  ~ifelse(is.na(.x), 1, 11-.x))) %>%
  unnest(c(prob_v, first1)) %>%
  mutate(p_del=n_del*prob_v,
         p_del=p_del/sum(p_del),
         expected=p_del*N,
         observed = del1,
         diff=observed-expected)

maxdiff <- arrange(counts, desc(abs(diff))) %>% 
  head(4)

ggplot(counts, aes(x=expected, y=observed)) +
  geom_point() +
  geom_abline() +
  geom_label(data=maxdiff, 
             aes(label=del_str, 
                 x=expected+c(300,300,-100,300),
                 y=observed+c(0,0,150,0))) +
  xlim(c(0, max(counts$expected,counts$observed))) +
  ylim(c(0, max(counts$expected,counts$observed))) +
  theme_bw() +
  ggtitle("Observed versus expected")
ggsave(file.path("figs","observed_expected.jpeg"),
       width=4)

ggplot(counts, aes(x=expected, y=observed)) +
  geom_point() +
  geom_abline() +
  scale_x_continuous(trans='log1p') +
  scale_y_continuous(trans='log1p') +
  theme_bw() +
  ggtitle("Observed versus expected - log scale")
ggsave(file.path("figs","observed_expected_log.jpeg"),
       width=4)

## E vs O of first place rankings
counts_first <- counts %>%
  group_by(first1) %>%
  summarise(expected=sum(expected), observed=sum(observed),
            diff=observed-expected, 
            percent_diff=(observed-expected)/expected) %>%
  mutate(item=items$name[STV_rank[1:J]]) %>%
  select(item, expected, observed, diff)
counts_first %>% xtable(digits=0) %>% print(include.rownames=F)

## E vs O of first second place given fatty tuna (top item) first
ft_first <- 
  filter(dat, rank==2) %>%
  semi_join(filter(dat, rank==1 & item==8), by="id") %>%
  count(name) %>%
  mutate(expected=sum(n)/(J-1),
         diff=n-expected,
         name=fct_relevel(name, item_names[STV_rank][2:10])) %>%
  select(item=name, expected, observed=n,diff) %>%
  arrange(item)
ft_first  %>% xtable(digits=0) %>% print(include.rownames=F)

