# 2017_B_allEST_TO03_Master filtering contigs
Anna A. Hippmann  
February 24, 2017  

I am exploring the 2nd new proteomics data I got from Jenny now (March 2017). What she did:
using my old proteomic peptide files from the LC-MS/MS runs, she let's them being re-searched again, this time against a database comprised of all protein sequences from both __and using now the second time around the older MASCOT version__:

1) gene ORFs of TO 1005 (publicly available, ~32K)
2) contigs of my ESTs of TO 1003 (in collaboration with Andy and John at JGVI, ~145K)

__NOTE__: this time, we used the complete set of contig sequences. However, the contig sequences can be divided into  three subgroups:

a) mapping to a TO05 gene ORF
b) mapping to a TO05 scaffold (aka genome)
c) not mapping to TO05 at all

Many peptides mapped to both contigs and or TO05 PRFs. This makes sense, as many many contigs map to these TO05 ORfs.  While exploring the data, I would like to adress a few questions:

1) to how many proteins/contigs could we map the proteomics peptides?
2) of these, how many created useable ratios?
3) of these, how many were part of mapping to 
  + TO05 ORFs (identifiers such as THAOC_#####)
  + contigs and TO05 ORFs
  + only contigs
  
Then, I would like to double-check, if those that hit only contigs, if these are indeed those that JCVI deemed not to have mapped to any region of TO05.

And last but not least looking into sig UP and DOWN regulated proteins and how all of this compares to my old/original  data (where the peptides have only been mapped to TO05 genome ORFs)

<a id="BackUP"></a>

And I do this for both

[TO 1003](#TO1003)

[TO 1005](#TO1005)




```r
suppressPackageStartupMessages(library(dplyr)) # to help manipulating the data
suppressPackageStartupMessages(library(readr))      # to help loading the data 
suppressPackageStartupMessages(library(ggplot2))    # to help plotting the data
suppressPackageStartupMessages(library(gridExtra))  # to help plotting multiplot plots
suppressPackageStartupMessages(library(knitr))
suppressPackageStartupMessages(library(VennDiagram))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(stringr))
```

and this is the summary of what I found out:





<a id="TO1003"></a>

# TO 1003

[Back Up](#BackUP)


```r
TO03 <- read.delim("Input_Data/2017_B_allEST_TO03_Master.txt") #2017_01_TO03_Master using all ESTs and TO05 genome to be mapped

ESTs <- read.delim("Input_Data/T.oceanica_annotation_denovo_2015_04.txt")

EST_narrow <- ESTs %>% 
  select(orf_id, hitting)

EST_NOTmapped <- EST_narrow %>% 
  filter(hitting == "NOTmapped" )

EST_onGenome <-  EST_narrow %>% 
  filter(hitting == "onGenome")


TO03_ratio <- TO03 %>% 
  filter(!is.na (Ratio.M.L.normalized)|!is.na (Ratio.H.L.normalized)|!is.na (Ratio.H.M.normalized)) #filters out only tohse that create ratio somewhere

TO03_ratio_narrow <- TO03_ratio %>% 
  select(original_row, Protein.IDs, Ratio.M.L.normalized, Ratio.H.L.normalized, Ratio.H.M.normalized)


TO03_ratio_narrow_info <- TO03_ratio_narrow %>% 
  mutate(THAOC = str_extract(Protein.IDs, "THAOC_\\d{5}")) %>% #uses RegEx to find first THAOC# in string and deposit result in new variable
  mutate (PS = str_extract(Protein.IDs, ";[a-z]{3}[A-Z]")) %>%  #deposits those on chloroplast genome e.g. PSI, PSII, ATP
  mutate (Rib = str_extract(Protein.IDs, ";rp[l,s]\\d")) %>% #RegEx to find ribosonmal proteins
  mutate (contig = str_extract(Protein.IDs,"^contig_")) %>% #getting all that START with a contig (as best hits)
  mutate (orf_id = str_extract(Protein.IDs, "^contig_\\d{1,6}_\\d{1,5}_\\d{2,5}_\\-|^contig_\\d{1,6}_\\d{1,5}_\\d{2,5}_\\+")) #getting all that START with a contig (as best hits) and depositing the actual contig identifyer into new variable


TO03_only_contigs <-TO03_ratio_narrow_info %>% 
  filter(contig == "contig_") %>%  #taking only those that START with a contig, 
  filter (is.na(Rib) & is.na(THAOC) & is.na(PS))# and then that do not have any other possible hits in TO05 chloroplast or nuclear genome

TO03_only_contigs_map_info <- left_join(TO03_only_contigs,EST_narrow, by = "orf_id" )
```

```
## Warning in left_join_impl(x, y, by$x, by$y, suffix$x, suffix$y): joining
## factor and character vector, coercing into character vector
```

```r
TO03_only_contigs_map_info_annot <- left_join(TO03_only_contigs, ESTs, by = "orf_id")
```

```
## Warning in left_join_impl(x, y, by$x, by$y, suffix$x, suffix$y): joining
## factor and character vector, coercing into character vector
```

```r
write.table(TO03_only_contigs_map_info_annot, "Output_Data/2017_02_TO03_proteomics_ProteinID_contigs_only.txt", sep="\t", row.names = FALSE, col.names = TRUE)
```

As to answering some of my questions:

1) How many protein identification all together hits? 1168
2) How many resulted in a ratio? 1005
3) Of these, how many were only contigs? 184

and here is a table with information on these contigs and where they mapped to: 

|hitting     |   n|
|:-----------|---:|
|Chloroplast |   2|
|NOTmapped   |  18|
|onGenome    |  32|
|RefGene     | 132|

<a id="TO1005"></a>

# TO 1005

[Back Up](#BackUP)


```r
TO05 <- read.delim("Input_Data/2017_B_allEST_TO05_Master.txt")



TO05_ratio <- TO05 %>% 
  filter(!is.na (Ratio.M.L.normalized)|!is.na (Ratio.H.L.normalized)|!is.na (Ratio.H.M.normalized)) #filters out only tohse that create ratio somewhere

TO05_ratio_narrow <- TO05_ratio %>% 
  select(original_row, Protein.IDs, Ratio.M.L.normalized, Ratio.H.L.normalized, Ratio.H.M.normalized)


TO05_ratio_narrow_info <- TO05_ratio_narrow %>% 
  mutate(THAOC = str_extract(Protein.IDs, "THAOC_\\d{5}")) %>% #uses RegEx to find first THAOC# in string and deposit result in new variable
  mutate (PS = str_extract(Protein.IDs, ";[a-z]{3}[A-Z]")) %>%  #deposits those on chloroplast genome e.g. PSI, PSII, ATP
  mutate (Rib = str_extract(Protein.IDs, ";rp[l,s]\\d")) %>% #RegEx to find ribosonmal proteins
  mutate (contig = str_extract(Protein.IDs,"^contig_")) %>% #getting all that START with a contig (as best hits)
  mutate (orf_id = str_extract(Protein.IDs, "^contig_\\d{1,6}_\\d{1,5}_\\d{2,5}_\\-|^contig_\\d{1,6}_\\d{1,5}_\\d{2,5}_\\+"))


TO05_only_contigs <-TO05_ratio_narrow_info %>% 
  filter(contig == "contig_") %>%  #taking only those that START with a contig, 
  filter (is.na(Rib) & is.na(THAOC) & is.na(PS))# and then that do not have any other possible hits in TO05 chloroplast or nuclear genome

TO05_only_contigs_map_info <- left_join(TO05_only_contigs,EST_narrow, by = "orf_id" )
```

```
## Warning in left_join_impl(x, y, by$x, by$y, suffix$x, suffix$y): joining
## factor and character vector, coercing into character vector
```

```r
TO05_only_contigs_map_info_annot <- left_join(TO05_only_contigs, ESTs, by = "orf_id")
```

```
## Warning in left_join_impl(x, y, by$x, by$y, suffix$x, suffix$y): joining
## factor and character vector, coercing into character vector
```

```r
write.table(TO05_only_contigs_map_info_annot, "Output_Data/2017_02_TO05_proteomics_ProteinID_contigs_only.txt", sep="\t", row.names = FALSE, col.names = TRUE)
```

As to answering some of my questions:

1) How many protein identification all together hits? 2448
2) How many resulted in a ratio? 2049
3) Of these, how many were only contigs? 382

and here is a table with information on these contigs and where they mapped to: 

|hitting     |   n|
|:-----------|---:|
|Chloroplast |   1|
|NOTmapped   |  31|
|onGenome    |  92|
|RefGene     | 258|
