library(tidyverse)
library(GenomicRanges)
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
CAGE_tss_file<-"./data/CAGE_tbl/HMEC_CAGE_TSS_tbl.Rda"
CAGE_enh_file<-"./data/CAGE_tbl/HMEC_CAGE_enh_tbl.Rda"
candidate_hub_file<-"~/Documents/multires_bhicect/Bootstrapp_fn/data/candidate_compound_hub/HMEC_5kb_tss_compound_hub.Rda"
spec_res_file<-"~/Documents/multires_bhicect/data/HMEC/spec_res/"
#-----------------------------------------
enh_tbl<-data_tbl_load_fn(CAGE_enh_file)
tss_tbl<-data_tbl_load_fn(CAGE_tss_file)

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
hub_GRange<-GRangesList(top_compound_hub_5kb_tbl$GRange)

enh_GRange<-GRanges(seqnames=enh_tbl$chr,
                    ranges = IRanges(start=enh_tbl$start,
                                     end=enh_tbl$end
                    ))
tss_GRange<-GRanges(seqnames=tss_tbl$chr,
                    ranges = IRanges(start=tss_tbl$start,
                                     end=tss_tbl$end
                    ))
countOverlaps(hub_GRange,enh_GRange)
countOverlaps(hub_GRange,tss_GRange)



hub_order<-tibble(chr=top_compound_hub_5kb_tbl$chr,hub=top_compound_hub_5kb_tbl$parent.hub,enh=countOverlaps(hub_GRange,enh_GRange),tss=countOverlaps(hub_GRange,tss_GRange)) %>% 
  mutate(f=(enh)/(tss+enh)) %>% arrange(desc(f)) %>% 
  mutate(ID=paste(chr,hub,sep="_")) %>% 
  dplyr::select(ID)

hub_re_count_tbl<-tibble(hub=rep(paste(top_compound_hub_5kb_tbl$chr,top_compound_hub_5kb_tbl$parent.hub,sep="_"),2),count=c(countOverlaps(hub_GRange,enh_GRange),countOverlaps(hub_GRange,tss_GRange)),set=rep(c('enh','tss'),each=nrow(top_compound_hub_5kb_tbl)))

hub_re_count_tbl %>% 
  mutate(chr=str_split_fixed(hub,pattern = "_",n=2)[,1]) %>% 
  left_join(.,chr_ratio_tbl) %>% 
  mutate(hub=fct_relevel(hub,hub_order$ID)) %>% 
  ggplot(.,aes(hub,count,fill=set))+
  geom_bar(stat='identity',position='fill')+
  theme( axis.text.x=element_blank(),
         axis.ticks.x=element_blank())

enh_tbl %>% 
  group_by(chr) %>% 
  summarise(n=n(),set="enh") %>% 
  bind_rows(.,tss_tbl %>% 
              group_by(chr) %>% 
              summarise(n=n(),set='tss')) %>%
  mutate(chr=fct_relevel(chr,paste0('chr',c(1:22,'X','M')))) %>% 
  ggplot(.,aes(chr,n,fill=set))+
  geom_bar(stat='identity',position='fill')

enh_tss_count_tbl<-tibble(chr=top_compound_hub_5kb_tbl$chr,hub=top_compound_hub_5kb_tbl$parent.hub,enh=countOverlaps(hub_GRange,enh_GRange),tss=countOverlaps(hub_GRange,tss_GRange))
chr_ratio_tbl<-enh_tss_count_tbl %>% 
  group_by(chr) %>% 
  summarise(f=sum(enh)/(sum(tss)+sum(enh)))

dt <- enh_tss_count_tbl %>% 
  left_join(.,chr_ratio_tbl) %>% 
  mutate(f.dev=(enh/(tss+enh))/f)
dens <- density(dt$f.dev)
df <- data.frame(x=dens$x, y=dens$y)
probs <- c(0, 0.25, 0.5, 0.75, 1)
quantiles <- quantile(dt$f.dev, prob=probs)
df$quant <- factor(findInterval(df$x,quantiles))
ggplot(df, aes(x,y)) + geom_line() + 
  geom_ribbon(aes(ymin=0, ymax=y, fill=quant)) +
  geom_vline(xintercept = 1,lty=2)+
  scale_x_continuous(breaks=quantiles) + 
  scale_fill_brewer(guide="none")
