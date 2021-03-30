#---------------------------------------#
####  Get TSS from Ensembl gtf file  ####
#---------------------------------------#



library(GenomicFeatures)
library(tidyverse)

SsalTxDB <- makeTxDbFromGFF(file = "ref_genome/Salmo_salar.ICSASG_v2.103.chr.gtf")


TSS <- promoters(transcripts(SsalTxDB), upstream=0, downstream=1)
TSS <- unique(TSS)

TSSbed <- data.frame(seqnames=seqnames(TSS),
                     start=start(TSS)-1,
                     end=end(TSS)) 

write.table(TSSbed,file = "ref_genome/SsalTSS_ENSEMBL.bed",sep="\t",col.names = F,row.names = F,quote = F, append = F)




#-------------------------------#
####  Using NCBI annotation  ####
#-------------------------------#



# generate TxDB from gff
SsalTxDB <- makeTxDbFromGFF(file = "ref_genome/GCF_000233375.1_ICSASG_v2_genomic.gff")


TSS <- promoters(transcripts(SsalTxDB), upstream=0, downstream=1)
TSS <- unique(TSS)

TSSbed <- data.frame(seqnames=seqnames(TSS),
                     start=start(TSS)-1,
                     end=end(TSS)) %>% 
  # Only chromosomes
  filter(grepl("NC_0273...1",seqnames)) %>% 
  # convert seqnames to ssaXX 
  mutate( seqnames = sprintf("ssa%02i", 1 + as.integer(sub("NC_0273(..).1","\\1",seqnames)) ) )

write.table(TSSbed,file = "SsalTSS_NCBI.bed",sep="\t",col.names = F,row.names = F,quote = F, append = F)

