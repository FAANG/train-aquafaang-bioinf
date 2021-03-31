# ATAC-seq tutorial (Under development)

> This is an ATAC-seq tutorial currently under development for the EBI aqua-faang course.

Read mapping and peak calling covered on day 1 with chip data. Since read mapping is the same for both we can start the ATAC-seq session with already mapped reads (.bam files).

Rough overview:

* Peak calling with Genrich
  + convert to BigWig
* Quality control with [ATAQV](https://github.com/ParkerLab/ataqv/)
* Peak annotation - Connect peaks to genes with [ChIPseeker](https://bioconductor.org/packages/release/bioc/html/ChIPseeker.html)
* [DiffBind](https://bioconductor.org/packages/release/bioc/html/DiffBind.html) to find differential open chromatin
* Visualise with [IGV](https://software.broadinstitute.org/software/igv/) (ATAC-seq pile-ups, peaks and DiffBind results)

## TODO:

> Currently this repo only contains an example snakemake workflow and not an actual tutorial...

* Generate a tutorial from this (document + slides + pre-recorded video)
* Install sofware and data on the VM
  + List required software (e.g. IGV, docker/singularity?, snakemake, R, R packages, RStudio?)
  + Add data to VM
  + Test the workflow in the VM

## Files in this repo:

* data/ - directory for input bam files. (needs to be downloaded)
* ref_genome/ - directory for reference genome files (download files using the script "prep_ref_genomes.sh")
* get_TSS.R - Script that extracts TSS from gtf file
* workflow/Snakefile - Snakemake workflow that runs Genrich and ATAQV.
* config/workflow_config.yaml - Configuration for workflow.
* config/samples.tsv - Sample table (paths to bam files and sample metadata)
* diffBindAnalyses.Rmd - Example Rmarkdown script using diffBind to detect differtially accessible peaks.
* annotate_peaks.R - Example R script using ChIPseeker to annotate peaks.



## Tutorial steps

### Input Data

The input data are bam files for 8 ATAC-seq samples mapped to the Atlantic salmon genome. There are four samples from brain (B) and four from liver (L). (These files can be downloaded from http://arken.nmbu.no/~lagr/atac-seq-training/data and placed in the "data" folder.)

The bam files are are sorted by mapping position. Duplicates/multimapped reads and reads mapping to mitochondrial chromosome or unplaced scaffolds have been removed.

### Peak calling

> Note: Genrich requires that the bam files are sorted by read names, not by mapping position.

* sort bam file by reads

  samtools sort -n -o {output} -@ {threads} {input}

* run genrich

  Genrich -j -q 0.05 -t {sampleID}.sortedByRead.bam -o {sampleID}.narrowPeak -k {sampleID}.k.bedGraph

> See https://github.com/jsh58/Genrich for more details on parameters

(What about blacklists?)

At this point we could already visualise results in IGV. However, the bedgraph(ish) files are not suitable for viewing in a genome browser so they should be converted to BigWig (.bw) file. 

* convert bedGraph to .bw

> The output .bedGraph from Genrich contains some extra columns. We just want to keep column 1-4 and drop the header before running bedGraphToBigWig

  cat B-1.k.bedGraph | grep ^ssa | cut -f 1-4 > B-1.bedGraph
  
  bedGraphToBigWig B-1.bedGraph ~/genomes/Salmo_salar/ENSEMBL/chrmSizes.txt B-1.bw

### Setting up IGV

Genomes -> Load genome from File...
ref_genome/Salmo_salar.ICSASG_v2.dna.toplevel.fa

Add tracks:

Genes:
ref_genome/Salmo_salar.ICSASG_v2.103.chr.gtf

Pileup and peaks:
results/Genrich/B-1.bw 
results/Genrich/B-1.narrowPeak
etc..

### ATAQV

* bam file needs to sorted by position and indexed
* Get TSS positions and list of chromosome names (see get_TSS.R and ref_genome/prep_ref_genomes.sh)
  + Note that there is a clear difference in TSS enrichment depending on whether you use the NCBI or Ensembl gene annotations for defining TSS loci
* run ATAQV 

ataqv Species_name {input.bam} \
      --metrics-file {output.ataqv.json} \
      --peak-file {peaks} \
      --tss-file {tssfile} \
      --autosomal-reference-file {chrmList} 

* generate html report and serve?

mkarv {outdir} {*.ataqv.json}

### diffBind

See diffBindAnalysis.Rmd

diffBind peaks with differential binding results:

results/diffBind/diffBind.bed

### Annotating peaks

See "annotate_peaks.R" for an example.
