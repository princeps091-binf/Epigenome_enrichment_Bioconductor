library(tidyverse)
library(GenomicRanges)
library(furrr)
options(scipen = 999999999)
res_set <- c('1Mb','500kb','100kb','50kb','10kb','5kb')
res_num <- c(1e6,5e5,1e5,5e4,1e4,5e3)
names(res_num)<-res_set

#-----------------------------------------
#Utils. Fn
data_tbl_load_fn<-function(file){
  out_tbl<-get(base::load(file))
  tmp_obj<-names(mget(base::load(file)))
  rm(list=tmp_obj)
  rm(tmp_obj)
  return(out_tbl)
}

#-----------------------------------------
tbl_file<-"./data/CAGE_tbl/HMEC_CAGE_TSS_tbl.Rda"
peak_tbl<-data_tbl_load_fn(tbl_file)
plan(multisession,workers=5)
peak_tbl<-peak_tbl %>% 
  mutate(GRange=future_pmap(list(chr,start,end),function(chr,start,end){
    return(GRanges(seqnames=chr,
                                ranges = IRanges(start=start,
                                                 end=end
                                )))
    
  }))
plan(sequential)

save(peak_tbl,file=tbl_file)
