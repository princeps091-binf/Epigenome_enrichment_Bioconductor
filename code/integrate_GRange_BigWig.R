library(GenomicRanges)
library(tidyverse)

data_tbl_load_fn<-function(file){
  out_tbl<-get(base::load(file))
  tmp_obj<-names(mget(base::load(file)))
  rm(list=tmp_obj)
  rm(tmp_obj)
  return(out_tbl)
}

CAGE_enh_GRange_file<-"~/Documents/multires_bhicect/Bootstrapp_fn/data/GRanges/CAGE_enh_GM12878_Grange.Rda"

cage_enh_GRange<-data_tbl_load_fn(CAGE_enh_GRange_file)

bwf_bioc <-BigWigFile("~/Documents/multires_bhicect/data/epi_data/HMEC/Bioconductor/H3K27ac/hmec_H3K27ac_FC.bigWig")
bwf_manual <-BigWigFile("~/Documents/multires_bhicect/data/epi_data/HMEC/ENCODE/H3K27ac/ENCFF981WTU_FC_HMEC_H3K27ac.bigWig")

test_manual<-import(bwf_manual,which=cage_enh_GRange,as="NumericList")
test_bioc<-import(bwf_bioc,which=cage_enh_GRange,as="NumericList")

enh.n<-1234
plot(seq(start(cage_enh_GRange[enh.n]),end(cage_enh_GRange[enh.n])),test[[enh.n]],type='l')


plot(unlist(test_bioc),unlist(test_manual))
