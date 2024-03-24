library(tidyverse)
library(combinat)
library(tictoc)

source("functions.R")

J <- 4
save_v <- T

tic()
v_info <- tibble(v=permn(J)) %>%
  mutate(del=map(v,get_del_one),
         v_str=map(v,paste, collapse=""),
         del_str=map(del,paste, collapse="")) %>%
  unnest(c(v_str,del_str)) %>%
  arrange(v_str)

if (save_v){
  saveRDS(v_info, 
          file.path("delta_counts",paste0("v_",J,".RDS")))
  
}
  
del_info <- v_info %>%
  group_by(del_str) %>%
  mutate(n_del=n()) %>%
  ungroup() %>%
  select(del,del_str,n_del) %>%
  distinct() %>% 
  arrange(del_str)
  
saveRDS(del_info,
        file.path("delta_counts",paste0("delta_",J,".RDS")))

toc()

print(nrow(del_info)) # 2^{J-1}



