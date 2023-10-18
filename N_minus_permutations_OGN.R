library(tidyverse)
library(combinat)
library(gtools)
library(data.table)
10-> full_length_OGN
5-> n_minus_OGN
77.93 ->PS_linkage

combinations(n= 10, r = n_minus_OGN, v = c(1:10))->test
print(combinations(n= 10, r = n_minus_OGN, v = c(1:10)))
t(test) %>%as_tibble() %>% map(sum)%>% as_vector()
#no repeats

### EXAMPLE
#c(329.03470,320.02320,345.02970,329.03480,329.034700000001,403.06960,394.06200,393.07590,393.07600,315.14300)-> mz_diff_10mer

#combinations(n= 10, r = 5, v = mz_diff_10mer)%>%as_tibble() ->fivemer_comb_from10mer
#t(fivemer_comb_from10mer) %>%as_tibble() %>% map(sum)%>% as_vector()->sum_fivemer_comb_from10mer
#sum_fivemer_comb_from10mer%>%round(digits=2) %>%unique() ->sum_fivemer_comb_from10mer_rounded_digit2
#no repeat fivemers, only uniques 


c(329.03470,320.02320,345.02970,329.03480,329.034700000001,403.06960,394.06200,393.07590,393.07600,393.07601)-> mz_diff_10mer
c( "MOErA", "T", "G","A","A", "2MOErA", "2MOErT", "2MOErC", "2MOErC", "2MOErC") -> lbls

temp <- sort(mz_diff_10mer, index.return=TRUE)
mz_diff_10mer <- temp$x
lbls <- lbls[temp$ix] 
  
combinations(n= full_length_OGN, r = n_minus_OGN, v = mz_diff_10mer)%>%as_tibble() ->sixmer_comb_from10mer
t(sixmer_comb_from10mer) %>%as_tibble() %>% map(sum)%>% as_vector()->sum_sixmer_comb_from10mer
sum_sixmer_comb_from10mer%>%round(digits=2) %>%unique() ->sum_sixmer_comb_from10mer_rounded_digit2
sum_sixmer_comb_from10mer_rounded_digit2%>%as_vector()-PS_linkage %>% sort()->permutation_N_minus_failed_seq
permutation_N_minus_failed_seq%>% range()
print(permutation_N_minus_failed_seq)
#no repeat fivemers, only uniques 


combinations(n= full_length_OGN, r = n_minus_OGN, v = mz_diff_10mer)%>%as_tibble() ->sixmer_comb_from10mer
combinations(n= full_length_OGN, r = n_minus_OGN, v = 1:length(mz_diff_10mer))%>%as_tibble() -> indices_comb
rowSums(sixmer_comb_from10mer) -> sums_sixmer_comb_from10mer
sums_sixmer_comb_from10mer %>% round(digits=2) %>% as_vector()-PS_linkage -> sums_adjusted
unique(sums_adjusted) -> unique_sums
sums_adjusted %>% as.data.table() -> test
test[, list(list(.I)), by = test] -> big_giant_list

combos <- list()
combos_sums <- list()
for (i in 1:length(big_giant_list[[1]])) {
  indices_comb[big_giant_list[[i,2]],] -> almost_there
  lbls[t(almost_there) %>% as_vector()] -> almost_almost_there
  t(array(almost_almost_there,dim=c(n_minus_OGN,length(almost_almost_there) / n_minus_OGN))) -> voila
  combos[[i]] <- voila
  
  mz_diff_10mer[t(almost_there) %>% as_vector()] -> almost_almost_there2
  t(array(almost_almost_there2,dim=c(n_minus_OGN,length(almost_almost_there2) / n_minus_OGN))) -> voila2
  combos_sums[[i]] <- rowSums(voila2) - PS_linkage 
}

for (i in 1:length(big_giant_list[[1]])) {
  print(combos_sums[[i]])
  print(combos[[i]])
}
