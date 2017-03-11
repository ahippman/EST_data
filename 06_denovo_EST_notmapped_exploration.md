# 06_denovo_EST_notmapped
Anna A. Hippmann  
March 10, 2017  

I would like to explore a bit more the ~45K EST contig ORFs that DID NOT MAP to any of the TO05 genome (neither scaffold, called gene ORFs or Chloroplast genome). Andy suggest in his email (March 7th, 2017) that I should look into the best hit species and look at the percent cover of the orfs on the contigs.

Hence, I will do the latter here.

First I will read in the fasta files of both contig ORFs and contigs, then calculate the respective lengths in nt space and finally calculate the percentage. I will do this for all contig_ORFs, so I can compare the subset of mapped vs notmapped contig ORFs.



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

As done in "04_Creating_new_fasta_subfiles" I will now import the fasta file from the FULL contig sequence (NOT the contig_ORF sequence!) into a data.frame. This time I will use it to determine the length of the whole contig to be able to draw conclusions on how large the percentage of the OORF is on that contig.


```r
#input <- read_lines(file = "Input_Data/Tocean.assembly_contigSeq.fa") 
#output <- file("Input_Data/Tocean.assembly_contigSeq.csv","w")

#currentSeq <- 0
#newLine <- 0

#for(i in 1:length(input)) {
#  if(strtrim(input[i], 1) == ">") {
#    if(currentSeq == 0) {
#      writeLines(paste(input[i],"\t"), output, sep="")
#      currentSeq <- currentSeq + 1
#    } else {
#      writeLines(paste("\n",input[i],"\t", sep =""), output, sep="")
#    }
#  } else {
#    writeLines(paste(input[i]), output, sep="")
#  }
#}

#close(output)
```


```r
EST_info <- read.delim("Input_Data/T.oceanica_annotation_denovo_2015_04.txt")

EST_narrow <- EST_info %>% 
  select(orf_id, hitting) %>% 
  mutate(sign = ">") 

EST_narrow <- unite(EST_narrow, orf_id, c(sign, orf_id), sep = "")

contig_ORF <- read.delim("Input_Data/contig_fasta_all.csv", header = FALSE)

full_contigs <- read.delim("Input_Data/Tocean.assembly_contigSeq.csv", header = FALSE)
```

Now I would like to count the nts of the contig ORFs and the full contigs and store that infomation in a respective new variable

```r
#contig_ORF_DATA
contig_ORF <- contig_ORF %>% 
  rename(orf_id = V1, AA_seq = V2) 

contig_ORF$AA_seq <- as.character(contig_ORF$AA_seq)

contig_ORF_length <- contig_ORF %>% 
  mutate(ORF_AA_length = nchar(AA_seq)) %>% 
  mutate(ORF_nt_length = ORF_AA_length * 3)

contig_ORF_length <- contig_ORF_length %>% 
  mutate(contig_id = str_extract(orf_id, "contig_\\d{1,6}"))  #uses RegEx to find first THAOC# in string and deposit 

#FULL contig DATA
full_contigs$V2 <- as.character(full_contigs$V2)

full_contigs_length <- full_contigs %>% 
  mutate(contig_id = str_extract(V1, "contig_\\d{1,6}")) %>% 
  mutate(full_contig_nt = nchar(V2)) %>% 
  select(contig_id, full_contig_nt )


contig_ORF_and_full_contig_length <- left_join(contig_ORF_length, full_contigs_length, by = "contig_id")

contig_ORF_and_full_contig_length <- contig_ORF_and_full_contig_length %>% 
  mutate(diff_full.vs.ORF = full_contig_nt-ORF_nt_length) %>% 
  mutate(ORF_percent_cover = 100*ORF_nt_length/full_contig_nt)

#Adding hit information
contig_ORF_and_full_contig_length <- left_join(contig_ORF_and_full_contig_length, EST_narrow, by = "orf_id")
```

```
## Warning in left_join_impl(x, y, by$x, by$y, suffix$x, suffix$y): joining
## character vector and factor, coercing into character vector
```

```r
fs <- c("min", "max", "mean")
vars <- c("ORF_nt_length", "ORF_percent_cover")

test <- contig_ORF_and_full_contig_length %>% 
  group_by(hitting) %>% 
  summarise_each_(funs_(fs ), var = vars)

kable(test, format = "markdown")
```



|hitting     | ORF_nt_length_min| ORF_percent_cover_min| ORF_nt_length_max| ORF_percent_cover_max| ORF_nt_length_mean| ORF_percent_cover_mean|
|:-----------|-----------------:|---------------------:|-----------------:|---------------------:|------------------:|----------------------:|
|Chloroplast |                87|             0.5423048|              4089|              99.68847|           724.3729|               24.30610|
|NOTmapped   |                60|             0.9398076|              7158|             100.49261|           335.6065|               77.25193|
|onGenome    |                60|             0.3916646|             13584|             100.42373|           289.0090|               62.22049|
|RefGene     |                63|             0.6863811|             25455|             100.47170|           478.9328|               82.34433|
|NA          |               219|            98.2062780|               219|              98.20628|           219.0000|               98.20628|

```r
write.table(test, "Output_Data/Overview_contig_length_stats_per_hit_group.txt", sep="\t", row.names = FALSE, col.names = TRUE)
```

#Best Hit SPecies

```r
EST_info <- read.delim("Input_Data/T.oceanica_annotation_denovo_2015_04.txt")

EST_species <- EST_info %>% 
  select(orf_id, best_hit_species) %>% 
  mutate(sign = ">") %>% 
  mutate(Species = str_extract(best_hit_species, "^([A-Z])\\w+"))

EST_species <- unite(EST_species, orf_id, c(sign, orf_id), sep = "")

contig_ORF_hitting <- contig_ORF_and_full_contig_length %>% 
  select(orf_id, hitting)

EST_species <- left_join(EST_species, contig_ORF_hitting, by = "orf_id")



  
test_summary<- EST_species %>% 
  group_by(hitting) %>% 
  count(Species) %>%
  rename(occur = n) %>% 
  mutate(percent = 100*occur/sum(occur))

Summary_NOTmapped <- test_summary %>% 
  filter(hitting == "NOTmapped") %>% 
  rename(NOTmapped = occur, NOTmapped.percent = percent)

Summary_Chloroplast <- test_summary %>% 
  filter(hitting == "Chloroplast") %>% 
  rename(Chloroplast = occur, Chloroplast.percent = percent)

Summary_onGenome <- test_summary %>% 
  filter(hitting == "onGenome") %>% 
  rename(onGenome = occur, onGenome.percent = percent)


Summary_RefGene <- test_summary %>% 
  filter(hitting == "RefGene") %>% 
  rename(RefGene = occur, RefGene.percent = percent)

Species <- test_summary %>% 
  ungroup(hitting) %>% 
  select(Species) 

Species <- unique(Species)

Table_Species <- full_join(Species, Summary_onGenome,by = "Species")
Table_Species <- full_join(Table_Species, Summary_RefGene, by = "Species")
Table_Species <- full_join(Table_Species, Summary_Chloroplast, by = "Species")
Table_Species <- full_join(Table_Species, Summary_NOTmapped, by = "Species")

kable(head(arrange(Table_Species, desc(RefGene.percent))), format = "markdown") 
```



|Species       |hitting.x | onGenome| onGenome.percent|hitting.y | RefGene| RefGene.percent|hitting.x.x | Chloroplast| Chloroplast.percent|hitting.y.y | NOTmapped| NOTmapped.percent|
|:-------------|:---------|--------:|----------------:|:---------|-------:|---------------:|:-----------|-----------:|-------------------:|:-----------|---------:|-----------------:|
|Thalassiosira |onGenome  |     2873|       19.0782921|RefGene   |   37817|      45.9233983|Chloroplast |         109|          92.3728814|NOTmapped   |       535|         1.1253208|
|Pelagomonas   |onGenome  |     1280|        8.4999004|RefGene   |   23779|      28.8762326|Chloroplast |           5|           4.2372881|NOTmapped   |      1752|         3.6851626|
|NA            |onGenome  |     9323|       61.9098214|RefGene   |   17083|      20.7448876|NA          |          NA|                  NA|NOTmapped   |     23004|        48.3866897|
|Skeletonema   |onGenome  |       76|        0.5046816|RefGene   |     760|       0.9229125|NA          |          NA|                  NA|NOTmapped   |        38|         0.0799293|
|Detonula      |onGenome  |       67|        0.4449167|RefGene   |     595|       0.7225434|NA          |          NA|                  NA|NOTmapped   |        13|         0.0273442|
|Chaetoceros   |onGenome  |      182|        1.2085796|RefGene   |     249|       0.3023753|Chloroplast |           1|           0.8474576|NOTmapped   |       156|         0.3281309|

```r
write.table(Table_Species, "Output_Data/Summary_Table_contig_ORfs_Species_per_HITgroup.txt", sep="\t", row.names = FALSE, col.names = TRUE)
```
