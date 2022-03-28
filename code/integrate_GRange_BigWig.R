library(GenomicRanges)
library(rtracklayer)
library(tidyverse)

data_tbl_load_fn<-function(file){
  out_tbl<-get(base::load(file))
  tmp_obj<-names(mget(base::load(file)))
  rm(list=tmp_obj)
  rm(tmp_obj)
  return(out_tbl)
}

CAGE_tss_GRange_file<-"~/Documents/multires_bhicect/Bootstrapp_fn/data/GRanges/CAGE_tss_HMEC_Grange.Rda"
CAGE_enh_GRange_file<-"~/Documents/multires_bhicect/Bootstrapp_fn/data/GRanges/CAGE_enh_HMEC_Grange.Rda"

CAGE_enh_tbl_file<-"~/Documents/multires_bhicect/HiC_enrichment/HiC_zscore_explore/data/HMEC_enh_H3K27ac_pval_tbl.Rda"

cage_enh_GRange<-data_tbl_load_fn(CAGE_enh_GRange_file)
cage_tss_GRange<-data_tbl_load_fn(CAGE_tss_GRange_file)

bwf_manual <-BigWigFile("~/Documents/multires_bhicect/data/epi_data/HMEC/CTCF/ENCODE/BigWig/ENCFF101RHP_HMEC_FC.bigWig")

test_manual<-import(bwf_manual,which=cage_tss_GRange,as="NumericList")
cage_enh_H3K27ac_tbl<-data_tbl_load_fn(CAGE_enh_tbl_file)
comp_tbl<-tibble(as.data.frame(cage_enh_GRange)) %>% 
  mutate(H3K27ac.vec=as.list(test_manual)) %>% 
  left_join(.,cage_enh_H3K27ac_tbl %>% 
              dplyr::select(enh,fc,pval) %>% 
              mutate(seqnames=str_split_fixed(enh,"_",3)[,1],start=as.integer(str_split_fixed(enh,"_",3)[,2]),end=as.integer(str_split_fixed(enh,"_",3)[,3]))
  )

comp_tbl %>% 
  mutate(m=map_dbl(H3K27ac.vec,function(x)mean(x))) %>% 
  ggplot(.,aes(m,fc))+geom_point(alpha=0.1)+scale_x_log10()+scale_y_log10()

enh.n<-134
plot(seq(start(cage_tss_GRange[enh.n]),end(cage_tss_GRange[enh.n])),test_manual[[enh.n]],type='l')


