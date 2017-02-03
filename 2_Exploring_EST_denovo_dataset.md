# 2_Exploring_EST_Denovo_Data
Anna A. Hippmann  
February 3, 2017  



### Exploring EST de novo data

What I would like to do here is explore the denovo data a bit and getting some feel of the subversiveness of redundancies and how we might need to deal with them. The information that we do have about each transcript is if it:

* hits a predicted reference gene on the genome of _T. oceanica_ (CCMP 1005) ( __gene__)

* hits the chloroplast genome ( __chl__)

* hits any non-coding part of the genome ( __genome__)

* hits nothing on the genome ( __not__)

What I would like to explore in more detail

1) What is the proportion of each of the above defined 4 pools (i.e. ORF, chl, gene, not)

2) how much redundencies are there (i.e. in ORF and chl pool, how many transcripts are hitting the same reading frame?) 
  - Those transcripts that are hitting the same genes will be dubbed to belong the same __pool__ of transcripts
  
3) are the redundant transcripts of one pool expressed in a similar way or are there __sub-pools__ of expression

####Why do I want to know/deal with the redundancies?

With my proteomics data I would like to investigate the proteomic response of _T. oceanica_ (CCMP 1003) to low Cu. At the moment I am particularly looking into the carbon metabolism. Earlier work, especially by Allen et al, Gruber et al and Smith et al, has shown that many metabolic pathways include several isozymes that are often times targeted to different cellular subcompartments. Furthermore some are used within the whole cellular metabolism in differnet ways than we would find in plants or animals. All of this makes it crucial to know where a  protein is really targeted to.

_T. oceanica's_ publicly available genome is helpful. However, it is in a rather rudimentary _Scaffold_ state when compared to the golden standard of genomes such as _T. pseudonanas_ or _P. tricornutums_. Oftentimes, gene models in _T. oceanica_ are missing their target sequence. THe EST library could help enhance the current gene models, especially upstream of the genes of interest to me to allow a beeter cellular targeting prediction.


##Libraries used

```r
suppressPackageStartupMessages(library(visreg))
suppressPackageStartupMessages(library(dplyr))
library(tidyr)
library(ggplot2)
library(knitr)
#library (gridExtra)
#library(cowplot)
```



```r
EST_all <- read.delim("Input_Data/T.oceanica_annotation_denovo_2015.txt")

EST_small <- EST_all %>% 
  select(1:9) %>% 
  rename(transcript = orf_id, ref_gene=ref_gene_id, ref_genome = ref_genome_id)


names(EST_small)
```

```
## [1] "transcript"         "orig.orf_id_row."   "hitting"           
## [4] "ref_gene"           "ref_genome"         "ref_genome_start"  
## [7] "ref_genome_end"     "orf_contam_type"    "AnnotationCombined"
```
####Now I owuld like to have a summary of each column




