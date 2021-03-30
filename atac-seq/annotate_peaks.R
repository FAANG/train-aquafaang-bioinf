
# Annotating peaks with ChIPseeker

library(ChIPseeker)
library(GenomicFeatures)

SsalTxDB <- makeTxDbFromGFF(file = "ref_genome/Salmo_salar.ICSASG_v2.103.chr.gtf")

peakAnno <- annotatePeak("results/Genrich/B-1.narrowPeak", tssRegion=c(-3000, 3000),
                         TxDb=SsalTxDB)

peakAnnoTbl <- as.data.frame(peakAnno)
