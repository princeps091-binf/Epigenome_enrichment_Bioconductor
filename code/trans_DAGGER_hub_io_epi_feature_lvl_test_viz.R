library(tidyverse)
library(GenomicRanges)
library(rtracklayer)
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
Build_coord_fn<-function(top_compound_hub_5kb_tbl,spec_res_file){
  coord_tbl<-do.call(bind_rows,map(unique(top_compound_hub_5kb_tbl$chr),function(chromo){
    message(chromo)
    base::load(paste0(spec_res_file,chromo,"_spec_res.Rda"))
    tmp_tbl<-top_compound_hub_5kb_tbl %>% 
      filter(chr==chromo) %>% 
      mutate(bins=chr_spec_res$cl_member[node]) %>% 
      mutate(bins=map(bins,as.numeric)) 
    
  }))
  return(coord_tbl)
}

Build_GRange_fn<-function(chromo,res,bins,res_num){
  inter_cl_Grange<-   GRanges(seqnames=chromo,
                              ranges = IRanges(start=bins,
                                               end=bins + res_num[res]-1
                              ))
  inter_cl_Grange<-GenomicRanges::reduce(inter_cl_Grange)
  return(inter_cl_Grange)
  
}
#-----------------------------------------
#candidate_hub_file<-"~/Documents/multires_bhicect/Bootstrapp_fn/data/DAGGER_tbl/GM12878_union_trans_res_dagger_tbl.Rda"
candidate_hub_file<-"~/Documents/multires_bhicect/Bootstrapp_fn/data/DAGGER_tbl/trans_res/HMEC_union_top_trans_res_dagger_tbl.Rda"

spec_res_file<-"~/Documents/multires_bhicect/data/HMEC/spec_res/"

CAGE_enh_GRange_file<-"~/Documents/multires_bhicect/Bootstrapp_fn/data/GRanges/CAGE_enh_HMEC_Grange.Rda"
feature_bigWig_file<-"~/Documents/multires_bhicect/data/epi_data/HMEC/ENCODE/RNAP2/ENCFF728WWJ_HMEC_RNAP2_FC.bigWig"
#-----------------------------------------
cage_enh_GRange<-data_tbl_load_fn(CAGE_enh_GRange_file)
bwf_manual <-BigWigFile(feature_bigWig_file)
tmp_max<-summary(bwf_manual,cage_enh_GRange,type="max")


feature_tbl<-tibble(as.data.frame(unlist(tmp_max))) %>% 
  dplyr::rename(max.lvl=score) %>% 
  mutate(enh=paste(seqnames,start,end,sep="_"))

trans_hub_tbl<-data_tbl_load_fn(candidate_hub_file) %>% 
  mutate(res=str_split_fixed(node,"_",2)[,1])

plan(multisession,workers=3)
trans_hub_tbl<-Build_coord_fn(trans_hub_tbl,spec_res_file) %>% 
  mutate(GRange=future_pmap(list(chr,res,bins),function(chromo,res,bins){
    Build_GRange_fn(chromo,res,bins,res_num)
  }))
plan(sequential)

plan(multisession,workers=3)
trans_hub_tbl<-trans_hub_tbl %>% 
  mutate(peak.content=future_map(GRange,function(x){
    tmp_hit<-unique(queryHits(findOverlaps(cage_enh_GRange,x)))
    paste(seqnames(cage_enh_GRange)[tmp_hit],start(cage_enh_GRange)[tmp_hit],end(cage_enh_GRange)[tmp_hit],sep="_")
  }))
plan(sequential)

in_set<-trans_hub_tbl %>% 
  filter(res=="50kb") %>% 
  dplyr::select(peak.content) %>% 
  unnest(cols=c(peak.content)) %>% 
  distinct() %>% unlist


zero_tresh<-10**(floor(min(log10(feature_tbl$max.lvl[feature_tbl$max.lvl>0]),na.rm=T)) -1)
feature_tbl %>% 
  mutate(hub.io=ifelse(enh %in% in_set,"in","out")) %>% 
  ggplot(.,aes(max.lvl+zero_tresh,color=hub.io))+geom_density()+
  scale_x_log10(breaks=c(zero_tresh,0.1,1,10,100),labels=c(0,0.1,1,10,100))+ 
  scale_color_brewer(palette="Set1")+
  xlab("fold-change")
ggsave("~/Documents/multires_bhicect/Poster/img/F4/RNAP2_HMEC_enh_hub_io.svg")

in_vec<-feature_tbl %>% 
  mutate(hub.io=ifelse(enh %in% in_set,"in","out")) %>% 
  filter(hub.io=="in") %>% 
  dplyr::select(max.lvl) %>% 
  unlist

out_vec<-feature_tbl %>% 
  mutate(hub.io=ifelse(enh %in% in_set,"in","out")) %>% 
  filter(hub.io=="out") %>% 
  dplyr::select(max.lvl) %>% 
  unlist

wilcox.test(in_vec,out_vec,alternative = "greater")$p.value
