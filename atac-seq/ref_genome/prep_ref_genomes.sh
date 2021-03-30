
#  Note: requires samtools to index fasta file

#------------------#
####  Download  ####
#------------------#


wget http://ftp.ensembl.org/pub/release-103/fasta/salmo_salar/dna/Salmo_salar.ICSASG_v2.dna.toplevel.fa.gz

wget http://ftp.ensembl.org/pub/release-103/gtf/salmo_salar/Salmo_salar.ICSASG_v2.103.chr.gtf.gz

# get the NCBI annotations to get the TSS's
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/233/375/GCF_000233375.1_ICSASG_v2/GCF_000233375.1_ICSASG_v2_genomic.gff.gz

# unizip
gunzip *.gz

# index fasta
samtools faidx Salmo_salar.ICSASG_v2.dna.toplevel.fa

# get list of chromosomes (required for ataqv)
grep ^ssa Salmo_salar.ICSASG_v2.dna.toplevel.fa.fai | cut -f 1 > Ssal_chrm.txt

# get list of chromosome sizes (required for bedGraphToBigWig)
grep ^ssa Salmo_salar.ICSASG_v2.dna.toplevel.fa.fai | cut -f 1,2 > Ssal_chrm_sizes.txt
