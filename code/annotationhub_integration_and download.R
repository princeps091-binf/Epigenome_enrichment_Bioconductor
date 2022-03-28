library(AnnotationHub)
library(rtracklayer)

ahub <- AnnotationHub()

ahub.bw <- subset(ahub, rdataclass == "BigWigFile" & species == "Homo sapiens" & genome == "hg19")
ahub.gr <- subset(ahub, rdataclass == "GRanges" & species == "Homo sapiens" & genome == "hg19")

# Tissue code: https://www.ncbi.nlm.nih.gov/research/snpdelscore/tissues/
# HMEC = E119
# H1 = E003
# GM12878 = E116

grep("",mcols(ahub.bw)$title,value=T)

grep("TfbsH1hescCtcf",mcols(ahub.gr)$title,value=T)
grep("TfbsHmecCtcf",mcols(ahub.gr)$title,value=T)
grep("TfbsGm12878Ctcf",mcols(ahub.gr)$title,value=T)

grep("TfbsH1hesc",mcols(ahub.gr)$title,value=T)
ahub.gr[grep("TfbsHmec",mcols(ahub.gr)$title)]
grep("TfbsGm12878",mcols(ahub.gr)$title,value=T)


grep("Hmec|HMEC|hmec",mcols(ahub.gr)$title,value=T)
grep("H1hesc|\\.H1\\.",mcols(ahub.gr)$title,value=T)
grep("Gm12878|gm12878",mcols(ahub.gr)$title,value=T)

hmec_H3K27ac_FC_BW<-ahub.bw[[grep("E119-H3K27ac.fc.signal",mcols(ahub.bw)$title)]]
seqinfo(hmec_H3K27ac_FC_BW)
export.bw(hmec_H3K27ac_FC_BW,"~/Documents/multires_bhicect/data/epi_data/HMEC/Bioconductor/H3K27ac/hmec_H3K27ac_FC.bigWig")
ahub.bw[grep("E119-H3K27ac.imputed.pval.signal",mcols(ahub.bw)$title)]

ahub.bw[grep("E003-H3K27ac.fc.signal",mcols(ahub.bw)$title)]
ahub.bw[grep("E003-H3K27ac.imputed.pval.signal",mcols(ahub.bw)$title)]

ahub.bw[grep("E116-H3K27ac.fc.signal",mcols(ahub.bw)$title)]
ahub.bw[grep("E116-H3K27ac.imputed.pval.signal",mcols(ahub.bw)$title)]



