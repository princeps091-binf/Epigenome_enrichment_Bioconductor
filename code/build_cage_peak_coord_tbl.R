#Build the CAGE-I table for different cell lines across both enhancers and TSS
library(tidyverse)
library(vroom)
library(parallel)
#---------------------------------------------------------------------------
## Utility functions
cage_tbl_subset<-function(cage_tbl,ID_col,col_set){
  return(cage_tbl%>%dplyr::select(ID_col)%>%bind_cols(cage_tbl%>%dplyr::select(contains(col_set))))
}
cage_sub_fn<-function(cage_tbl){
  print("compute m")
  cl<-makeCluster(5)
  tmp_m<-parallel::parApply(cl,X = as.matrix(cage_tbl[,-1]),MARGIN = 1,function(x){
    mean(x)
  })
  stopCluster(cl)
  rm(cl)
  cage_tbl<-cage_tbl%>%mutate(m=tmp_m)%>%filter(m>0)
  return(cage_tbl)
}
data_tbl_load_fn<-function(file){
  out_tbl<-get(base::load(file))
  tmp_obj<-names(mget(base::load(file)))
  rm(list=tmp_obj)
  rm(tmp_obj)
  return(out_tbl)
}
#---------------------------------------------------------------------------

H1<-c("CNhs14067","CNhs14068","CNhs13964")
GM12878<-c('CNhs12331','CNhs12332','CNhs12333')
HMEC<-c('CNhs11077','CNhs11382','CNhs12032')
col_set<-GM12878
#---------------------------------------------------------------------------
# CAGE-peak for TSS with tpm normalisation
cage<-vroom("~/Documents/multires_bhicect/data/epi_data/hg19.cage_peak_phase1and2combined_tpm_ann.osc.txt",comment = '#')
colnames(cage)[1]<-"Id"

cage_tss_tbl<-cage %>% 
  cage_tbl_subset(.,"Id",col_set) %>%
  cage_sub_fn()

tmp_coord<-str_split_fixed(cage_tss_tbl$Id,pattern = ":|\\.\\.|,",4)
cage_tss_tbl<-cage_tss_tbl %>% 
  mutate(chr=tmp_coord[,1],
         start=as.numeric(tmp_coord[,2]),
         end=as.numeric(tmp_coord[,3])) %>% 
  filter(!(is.na(start)))

save(cage_tss_tbl,file = "./data/CAGE_tbl/GM12878_CAGE_TSS_tbl.Rda")
#---------------------------------------------------------------------------
# CAGE-peak for enhancers
cage_enh_tbl<-vroom("~/Documents/multires_bhicect/data/epi_data/human_permissive_enhancers_phase_1_and_2_expression_tpm_matrix.txt",comment = '#')
colnames(cage_enh_tbl)[1]<-"Id"

cage_enh_tbl<-cage_enh_tbl %>% 
  cage_tbl_subset(.,"Id",col_set) %>%
  cage_sub_fn() %>% 
  mutate(chr=str_split_fixed(Id,pattern = ":|-",3)[,1],
         start=as.numeric(str_split_fixed(Id,pattern = ":|-",3)[,2]),
         end=as.numeric(str_split_fixed(Id,pattern = ":|-",3)[,3])) %>% 
  filter(!(is.na(start)))

save(cage_enh_tbl,file = "./data/CAGE_tbl/GM12878_CAGE_enh_tbl.Rda")
#---------------------------------------------------------------------------
#Union of peaks
tss_file<-"./data/CAGE_tbl/GM12878_CAGE_TSS_tbl.Rda"
enh_file<-"./data/CAGE_tbl/GM12878_CAGE_enh_tbl.Rda"
tss_tbl<-data_tbl_load_fn(tss_file)
enh_tbl<-data_tbl_load_fn(enh_file)

tss_tbl %>% 
  mutate(set='tss') %>% 
  bind_rows(.,enh_tbl %>% 
              mutate(set='enh')) %>% 
  ggplot(.,aes(m,color=set))+geom_density()+scale_x_log10()
