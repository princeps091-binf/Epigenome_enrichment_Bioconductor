library(tidyverse)
library(GenomicRanges)
library(rtracklayer)
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
      mutate(bins=chr_spec_res$cl_member[parent.hub]) %>% 
      mutate(bins=map(bins,as.numeric)) 
    
  }))
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
candidate_hub_file<-"~/Documents/multires_bhicect/Bootstrapp_fn/data/candidate_compound_hub/H1_5kb_tss_compound_hub.Rda"
spec_res_file<-"~/Documents/multires_bhicect/data/H1/Dekker/spec_res/"

CAGE_GRange_file<-"./data/CAGE_tbl/H1_CAGE_enh_tbl.Rda"
#-----------------------------------------
compound_hub_5kb_tbl<-data_tbl_load_fn(candidate_hub_file)

top_compound_hub_5kb_tbl<-do.call(bind_rows,map(unique(compound_hub_5kb_tbl$chr),function(chromo){
  tmp_tbl<-compound_hub_5kb_tbl %>% 
    filter(chr==chromo)
  return(tmp_tbl %>% 
           filter(parent.hub %in% tmp_tbl$parent.hub[which(!(tmp_tbl$parent.hub %in% unique(unlist(tmp_tbl$ch.hub))))])
  )
  
}))

top_compound_hub_5kb_tbl<-Build_coord_fn(top_compound_hub_5kb_tbl,spec_res_file) %>% 
  mutate(GRange=pmap(list(chr,res,bins),function(chromo,res,bins){
    Build_GRange_fn(chromo,res,bins,res_num)
  }))

feature_tbl<-data_tbl_load_fn(CAGE_GRange_file)

cage_GRange<-GRanges(seqnames=feature_tbl$chr,
                         ranges = IRanges(start=feature_tbl$start,
                                          end=feature_tbl$end
                         ))
mcols(cage_GRange)<-tibble(Id=feature_tbl$Id)

top_compound_hub_5kb_tbl<-top_compound_hub_5kb_tbl %>% 
  mutate(peak.content=map(GRange,function(x){
    tmp_hit<-unique(queryHits(findOverlaps(cage_GRange,x)))
    return(mcols(cage_GRange)$Id[tmp_hit])
  }))


in_set<-top_compound_hub_5kb_tbl %>% 
  #  filter(res=="50kb") %>% 
  dplyr::select(peak.content) %>% 
  unnest(cols=c(peak.content)) %>% 
  distinct() %>% unlist


#zero_tresh<-10**(floor(min(log10(feature_tbl$m[feature_tbl$m>0]),na.rm=T)) -1)
feature_tbl %>% 
  mutate(hub.io=ifelse(Id %in% in_set,"in","out")) %>% 
  ggplot(.,aes(m,color=hub.io))+geom_density()+
  scale_x_log10()+ 
  scale_color_brewer(palette="Set1")+
  xlab("CAGE Intensity (tpm)")
ggsave("~/Documents/multires_bhicect/weeklies/weekly55/img/H1_enh_io_cage_I.png")

in_vec<-feature_tbl %>% 
  mutate(hub.io=ifelse(Id %in% in_set,"in","out")) %>% 
  filter(hub.io=="in") %>% 
  dplyr::select(m) %>% 
  unlist

out_vec<-feature_tbl %>% 
  mutate(hub.io=ifelse(Id %in% in_set,"in","out")) %>% 
  filter(hub.io=="out") %>% 
  dplyr::select(m) %>% 
  unlist

wilcox.test(in_vec,out_vec,alternative = "greater")$p.value

