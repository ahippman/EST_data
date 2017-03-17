# 5_Adding gene_ORF_id column to EST matched Proteomic Data
Anna A. Hippmann  
March 1, 2017  

I got the new dataset (2017_A_TO03_master_...). Old peptides have been re-searched against the combined protein sequences of TO05 genome (~32K) + my own generates EST library (~ 145K contig ORFs). THe contigs hit either the Chloroplast, Gene ORFs, somewhere else on the scaffolds or nothing. The protein hit column includes all proteins that have been hit by the LC-MS/MS peptides. This results in often including identifiers for many contigs as well as THAOC_ and/or Chloroplast encoded proteins.

I would like to compare the results with the original ones. FOr that, I would like to create a "gene_ORF_id" column where each row is associated with only one identifier. I will decide on the identifyer based on the following criteria:

1) if a chloroplast protein identifier is part of the list, use THAT ONE
2) if THAOC identifyer(s) are part, use the first one
3) if only contig identifiers are used, use the first one

Meanwhile, I also got the second re-searched proteomics dataset. In this one, Everything is the same, BUT the Mascot software version! The second re-searched dataset is __"2017_B_TO03_Master..."__



<a id="BackUP"></a>

And I do this for both

[TO 1003 2017_A](#TO1003)

[TO 1005 2017_A](#TO1005)

[TO 1003 2017_B](#TO1003_B)

[TO 1005 2017_B](#TO1005_B)

##TO 1003


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



```r
TO03 <- read.delim("Input_Data/2017_allEST_TO03_Master.txt") #2017_01_TO03_Master using all ESTs and TO05 genome to be mapped
TO05 <- read.delim("Input_Data/2017_allEST_TO05_Master.txt") #2017_01_TO05_Master using all ESTs and TO05 genome to be mapped

TO03_B <- read.delim("Input_Data/2017_B_allEST_TO03_Master.txt") #2017_01_TO03_Master using all ESTs and TO05 genome to be mapped
TO05_B <- read.delim("Input_Data/2017_B_allEST_TO05_Master.txt") #2017_01_TO05_Master using all ESTs and TO05 genome to be mapped
```



```r
TO03_ProteinIDs <- TO03 %>% 
  select(original_row, Protein.IDs)


TO03_ProteinIDs_info <- TO03_ProteinIDs %>% 
  mutate(THAOC = str_extract(Protein.IDs, "THAOC_\\d{5}")) %>% #uses RegEx to find first THAOC# in string and deposit result in new variable
  mutate (PS_start = str_extract(Protein.IDs, "^[a-z]{3}[A-Z]")) %>%  #deposits those on chloroplast genome e.g. PSI, PSII, ATP
  mutate (PS_middle = str_extract(Protein.IDs, ";[a-z]{3}[A-Z]")) %>%  #deposits those on chloroplast mutate (Rib = str_extract(Protein.IDs, "rp[l,s]\\d")) %>% #RegEx to find ribosonmal proteins %>% 
  mutate(PS_middle = str_extract(PS_middle, "[a-z]{3}[A-Z]")) %>%
  mutate (Rib = str_extract(Protein.IDs, "rp[l,s]\\d{1,2}")) %>% #RegEx to find ribosonmal proteins
  mutate (contig = str_extract(Protein.IDs,"^contig_\\d{1,6}_\\d{1,5}_\\d{2,5}_\\-|^contig_\\d{1,6}_\\d{1,5}_\\d{2,5}_\\+")) #getting all that START with a contig (as best hits)


TO03_ProteinIDs_info_comb <- unite(TO03_ProteinIDs_info, gene_ORF_comb, c(PS_start,PS_middle, Rib, THAOC, contig), sep = "; ", remove = FALSE)

TO03_ProteinIDs_info_comb <- TO03_ProteinIDs_info_comb %>% 
  mutate(colon = ";" )

TO03_ProteinIDs_info_comb <- unite(TO03_ProteinIDs_info_comb, together, c(gene_ORF_comb, colon), sep = "")

#now I would like to remove NAs
TO03_ProteinIDs_info_comb$together <- gsub("NA;", "", TO03_ProteinIDs_info_comb$together)
TO03_ProteinIDs_info_comb$together <- trimws(TO03_ProteinIDs_info_comb$together, which = "both")

#TO03_ProteinIDs_info_comb <- TO03_ProteinIDs_info_comb %>% 
#  mutate(gene_ORF_id = str_extract(together, "^THAOC_\\d{5}"))

TO03_ProteinIDs_info_comb <- TO03_ProteinIDs_info_comb %>% 
  mutate(gene_ORF_id_try = str_extract(together, "^.{4,27};"))
#this will lead to psbA etc and THAOC together
#so I will need to separate this column with using Sep=";"
TO03_ProteinIDs_info_comb <- separate(TO03_ProteinIDs_info_comb, gene_ORF_id_try, c("gene_ORF_id", "extra_1", "extra_2"), sep = ";")
```

```
## Warning: Too few values at 1236 locations: 1, 2, 3, 4, 5, 6, 7, 8, 9, 11,
## 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, ...
```

```r
#now I would like to remove ";" - not needed anymore
#TO03_ProteinIDs_info_comb$gene_ORF_id <- gsub(";", "", TO03_ProteinIDs_info_comb$gene_ORF_id)

#and removing whitespaces
TO03_ProteinIDs_info_comb$gene_ORF_id <- trimws(TO03_ProteinIDs_info_comb$gene_ORF_id, which = "both") #leading white space "left" trailing white space "right"



write.table(TO03_ProteinIDs_info_comb, "Output_Data/2017_A_TO03_proteomics_gene_ORF_ids_forMaster.txt", sep="\t", row.names = FALSE, col.names = TRUE)

kable(head(TO03_ProteinIDs_info_comb), format = "markdown")
```



| original_row|Protein.IDs                                                                         |together               |THAOC |PS_start |PS_middle |Rib |contig                |gene_ORF_id           |extra_1 |extra_2 |
|------------:|:-----------------------------------------------------------------------------------|:----------------------|:-----|:--------|:---------|:---|:---------------------|:---------------------|:-------|:-------|
|            1|atpAADB27547.1ATPsynthaseCF1subunitalpha;contig_46002_1_294_+;contig_119588_1_338_+ |atpA;                  |NA    |atpA     |NA        |NA  |NA                    |atpA                  |        |NA      |
|            2|atpFADB27545.1ATPsynthaseCF0subunitIBchain                                          |atpF;                  |NA    |atpF     |NA        |NA  |NA                    |atpF                  |        |NA      |
|            7|contig_100473_1_208_+;contig_110443_1_287_-                                         |contig_100473_1_208_+; |NA    |NA       |NA        |NA  |contig_100473_1_208_+ |contig_100473_1_208_+ |        |NA      |
|            8|contig_101299_1_202_+                                                               |contig_101299_1_202_+; |NA    |NA       |NA        |NA  |contig_101299_1_202_+ |contig_101299_1_202_+ |        |NA      |
|            9|contig_102546_1_216_-                                                               |contig_102546_1_216_-; |NA    |NA       |NA        |NA  |contig_102546_1_216_- |contig_102546_1_216_- |        |NA      |
|           10|contig_103488_1_207_+;contig_120062_1_241_+                                         |contig_103488_1_207_+; |NA    |NA       |NA        |NA  |contig_103488_1_207_+ |contig_103488_1_207_+ |        |NA      |

<a id="TO1005"></a>

##TO 1005

[Back Up](#BackUP)



```r
TO05_ProteinIDs <- TO05 %>% 
  select(original_row, Protein.IDs)


TO05_ProteinIDs_info <- TO05_ProteinIDs %>% 
  mutate(THAOC = str_extract(Protein.IDs, "THAOC_\\d{5}")) %>% #uses RegEx to find first THAOC# in string and deposit result in new variable
  mutate (PS_start = str_extract(Protein.IDs, "^[a-z]{3}[A-Z]")) %>%  #deposits those on chloroplast genome e.g. PSI, PSII, ATP
  mutate (PS_middle = str_extract(Protein.IDs, ";[a-z]{3}[A-Z]")) %>%  #deposits those on chloroplast mutate (Rib = str_extract(Protein.IDs, "rp[l,s]\\d")) %>% #RegEx to find ribosonmal proteins %>% 
  mutate(PS_middle = str_extract(PS_middle, "[a-z]{3}[A-Z]")) %>%
  mutate (Rib = str_extract(Protein.IDs, "rp[l,s]\\d{1,2}")) %>% #RegEx to find ribosonmal proteins
  mutate (contig = str_extract(Protein.IDs,"^contig_\\d{1,6}_\\d{1,5}_\\d{2,5}_\\-|^contig_\\d{1,6}_\\d{1,5}_\\d{2,5}_\\+")) #getting all that START with a contig (as best hits)


TO05_ProteinIDs_info_comb <- unite(TO05_ProteinIDs_info, gene_ORF_comb, c(PS_start,PS_middle, Rib, THAOC, contig), sep = "; ", remove = FALSE)

TO05_ProteinIDs_info_comb <- TO05_ProteinIDs_info_comb %>% 
  mutate(colon = ";" )

TO05_ProteinIDs_info_comb <- unite(TO05_ProteinIDs_info_comb, together, c(gene_ORF_comb, colon), sep = "")

#now I would like to remove NAs
TO05_ProteinIDs_info_comb$together <- gsub("NA;", "", TO05_ProteinIDs_info_comb$together)
TO05_ProteinIDs_info_comb$together <- trimws(TO05_ProteinIDs_info_comb$together, which = "both")

#TO05_ProteinIDs_info_comb <- TO05_ProteinIDs_info_comb %>% 
#  mutate(gene_ORF_id = str_extract(together, "^THAOC_\\d{5}"))

TO05_ProteinIDs_info_comb <- TO05_ProteinIDs_info_comb %>% 
  mutate(gene_ORF_id_try = str_extract(together, "^.{4,27};"))
#this will lead to psbA etc and THAOC together
#so I will need to separate this column with using Sep=";"
TO05_ProteinIDs_info_comb <- separate(TO05_ProteinIDs_info_comb, gene_ORF_id_try, c("gene_ORF_id", "extra_1", "extra_2"), sep = ";")
```

```
## Warning: Too few values at 2472 locations: 1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
## 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, ...
```

```r
#now I would like to remove ";" - not needed anymore
#TO05_ProteinIDs_info_comb$gene_ORF_id <- gsub(";", "", TO05_ProteinIDs_info_comb$gene_ORF_id)

#and removing whitespaces
TO05_ProteinIDs_info_comb$gene_ORF_id <- trimws(TO05_ProteinIDs_info_comb$gene_ORF_id, which = "both") #leading white space "left" trailing white space "right"



write.table(TO05_ProteinIDs_info_comb, "Output_Data/2017_A_TO05_proteomics_gene_ORF_ids_forMaster.txt", sep="\t", row.names = FALSE, col.names = TRUE)

kable(head(TO05_ProteinIDs_info_comb), format = "markdown")
```



| original_row|Protein.IDs                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              |together                             |THAOC       |PS_start |PS_middle |Rib |contig                 |gene_ORF_id           |extra_1 |extra_2 |
|------------:|:------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|:------------------------------------|:-----------|:--------|:---------|:---|:----------------------|:---------------------|:-------|:-------|
|            1|atpAADB27547.1ATPsynthaseCF1subunitalpha                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 |atpA;                                |NA          |atpA     |NA        |NA  |NA                     |atpA                  |        |NA      |
|            2|atpFADB27545.1ATPsynthaseCF0subunitIBchain                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               |atpF;                                |NA          |atpF     |NA        |NA  |NA                     |atpF                  |        |NA      |
|            9|contig_100071_60_204_-;contig_83969_1_154_+;contig_99757_91_200_-;contig_90515_98_207_-;contig_96649_1_131_+;contig_85924_210_344_-;contig_107383_1_132_+;contig_95379_79_219_-;contig_85911_51_200_-;contig_87982_50_203_-;contig_99688_63_216_-;contig_86910_62_215_-;contig_80773_63_216_-;contig_78821_50_216_-;contig_94062_28_205_-;contig_76920_31_207_-;contig_74862_27_207_-;contig_96797_38_225_-;contig_86128_38_225_-;contig_80567_38_225_-;contig_106363_1_200_-;contig_89761_1_201_-;contig_115461_1_200_+;contig_104579_1_202_+;contig_86514_1_203_-;contig_81818_1_203_+;contig_121266_1_217_+;contig_113607_1_222_-;contig_130909_1_231_-;contig_70370_8_407_-;THAOC_17785EJK61682.1hypotheticalprotein |THAOC_17785; contig_100071_60_204_-; |THAOC_17785 |NA       |NA        |NA  |contig_100071_60_204_- |THAOC_17785           |        |NA      |
|           10|contig_100119_1_202_-;contig_106652_1_211_+;THAOC_25974EJK54403.1hypotheticalprotein                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     |THAOC_25974; contig_100119_1_202_-;  |THAOC_25974 |NA       |NA        |NA  |contig_100119_1_202_-  |THAOC_25974           |        |NA      |
|           11|contig_100206_1_201_+                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    |contig_100206_1_201_+;               |NA          |NA       |NA        |NA  |contig_100206_1_201_+  |contig_100206_1_201_+ |        |NA      |
|           12|contig_102056_1_211_-                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    |contig_102056_1_211_-;               |NA          |NA       |NA        |NA  |contig_102056_1_211_-  |contig_102056_1_211_- |        |NA      |


```r
TO03_table <- TO03_ProteinIDs_info_comb %>% 
  summarise_each(funs(n_distinct)) %>% 
  mutate(chloroplast = PS_start + PS_middle) %>% 
  rename(ribosome = Rib) %>% 
  mutate(strain = "TO03") %>% 
  select(one_of(c("strain", "original_row", "together", "chloroplast", "ribosome", "THAOC", "contig", "gene_ORF_id"))) 
TO05_table <- TO05_ProteinIDs_info_comb %>% 
  summarise_each(funs(n_distinct)) %>% 
  mutate(chloroplast = PS_start + PS_middle) %>% 
  rename(ribosome = Rib) %>% 
  mutate(strain = "TO05") %>% 
  select(one_of(c("strain", "original_row", "together", "chloroplast", "ribosome", "THAOC", "contig", "gene_ORF_id"))) 
Sum_Table <- bind_rows(TO03_table, TO05_table) 
kable(Sum_Table, format="markdown")
```



|strain | original_row| together| chloroplast| ribosome| THAOC| contig| gene_ORF_id|
|:------|------------:|--------:|-----------:|--------:|-----:|------:|-----------:|
|TO03   |         1244|     1244|          31|       21|   936|    727|        1244|
|TO05   |         2480|     2480|          38|       31|  1883|   1380|        2480|

```r
write.table(Sum_Table,"Output_Data/2017_A_TO03_TO05_geneORF_id_Count_Table.txt", sep="\t", row.names = FALSE, col.names = TRUE)
```
<a id="TO1003_B"></a>

##TO 1003_B

[Back Up](#BackUP)



```r
TO03_B_ProteinIDs <- TO03_B %>% 
  select(original_row, Protein.IDs)


TO03_B_ProteinIDs_info <- TO03_B_ProteinIDs %>% 
  mutate(THAOC = str_extract(Protein.IDs, "THAOC_\\d{5}")) %>% #uses RegEx to find first THAOC# in string and deposit result in new variable
  mutate (PS_start = str_extract(Protein.IDs, "^[a-z]{3}[A-Z]")) %>%  #deposits those on chloroplast genome e.g. PSI, PSII, ATP
  mutate (PS_middle = str_extract(Protein.IDs, ";[a-z]{3}[A-Z]")) %>%  #deposits those on chloroplast mutate (Rib = str_extract(Protein.IDs, "rp[l,s]\\d")) %>% #RegEx to find ribosonmal proteins %>% 
  mutate(PS_middle = str_extract(PS_middle, "[a-z]{3}[A-Z]")) %>%
  mutate (Rib = str_extract(Protein.IDs, "rp[l,s]\\d{1,2}")) %>% #RegEx to find ribosonmal proteins
  mutate (contig = str_extract(Protein.IDs,"^contig_\\d{1,6}_\\d{1,5}_\\d{2,5}_\\-|^contig_\\d{1,6}_\\d{1,5}_\\d{2,5}_\\+")) #getting all that START with a contig (as best hits)


TO03_B_ProteinIDs_info_comb <- unite(TO03_B_ProteinIDs_info, gene_ORF_comb, c(PS_start,PS_middle, Rib, THAOC, contig), sep = "; ", remove = FALSE)

TO03_B_ProteinIDs_info_comb <- TO03_B_ProteinIDs_info_comb %>% 
  mutate(colon = ";" )

TO03_B_ProteinIDs_info_comb <- unite(TO03_B_ProteinIDs_info_comb, together, c(gene_ORF_comb, colon), sep = "")

#now I would like to remove NAs
TO03_B_ProteinIDs_info_comb$together <- gsub("NA;", "", TO03_B_ProteinIDs_info_comb$together)
TO03_B_ProteinIDs_info_comb$together <- trimws(TO03_B_ProteinIDs_info_comb$together, which = "both")

#TO03_B_ProteinIDs_info_comb <- TO03_B_ProteinIDs_info_comb %>% 
#  mutate(gene_ORF_id = str_extract(together, "^THAOC_\\d{5}"))

TO03_B_ProteinIDs_info_comb <- TO03_B_ProteinIDs_info_comb %>% 
  mutate(gene_ORF_id_try = str_extract(together, "^.{4,27};"))
#this will lead to psbA etc and THAOC together
#so I will need to separate this column with using Sep=";"
TO03_B_ProteinIDs_info_comb <- separate(TO03_B_ProteinIDs_info_comb, gene_ORF_id_try, c("gene_ORF_id", "extra_1", "extra_2"), sep = ";")
```

```
## Warning: Too few values at 1161 locations: 1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
## 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, ...
```

```r
#now I would like to remove ";" - not needed anymore
#TO03_B_ProteinIDs_info_comb$gene_ORF_id <- gsub(";", "", TO03_B_ProteinIDs_info_comb$gene_ORF_id)

#and removing whitespaces
TO03_B_ProteinIDs_info_comb$gene_ORF_id <- trimws(TO03_B_ProteinIDs_info_comb$gene_ORF_id, which = "both") #leading white space "left" trailing white space "right"



write.table(TO03_B_ProteinIDs_info_comb, "Output_Data/2017_B_TO03_proteomics_gene_ORF_ids_forMaster.txt", sep="\t", row.names = FALSE, col.names = TRUE)

kable(head(TO03_B_ProteinIDs_info_comb), format = "markdown")
```



| original_row|Protein.IDs                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      |together                           |THAOC       |PS_start |PS_middle |Rib |contig                |gene_ORF_id           |extra_1 |extra_2 |
|------------:|:--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|:----------------------------------|:-----------|:--------|:---------|:---|:---------------------|:---------------------|:-------|:-------|
|            1|atpA ADB27547.1 ATP synthase CF1 subunit alpha                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   |atpA;                              |NA          |atpA     |NA        |NA  |NA                    |atpA                  |        |NA      |
|            2|atpF ADB27545.1 ATP synthase CF0 subunit I B chain                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               |atpF;                              |NA          |atpF     |NA        |NA  |NA                    |atpF                  |        |NA      |
|            7|contig_100917_1_82_-;contig_87723_1_200_-;contig_101998_1_201_+;contig_95211_1_202_-;contig_89498_1_201_+;contig_88580_1_201_+;contig_86622_1_201_+;contig_86568_1_202_+;contig_83277_1_201_-;contig_82806_1_201_-;contig_82389_1_201_-;contig_80753_1_201_-;contig_79360_1_202_+;contig_76851_1_201_+;contig_114990_1_200_-;contig_107188_1_201_-;contig_102731_1_200_+;contig_90425_1_205_+;contig_87738_1_205_-;contig_86964_1_205_-;contig_85529_1_204_-;contig_85239_1_204_+;contig_82749_1_205_-;contig_82313_1_205_-;contig_76915_1_205_+;contig_75354_1_205_-;contig_74852_1_205_-;contig_74228_1_205_-;contig_102865_1_205_+;contig_93455_1_206_-;contig_81506_1_206_+;contig_77405_1_208_-;contig_76031_1_208_-;contig_73026_1_207_-;contig_106255_1_208_-;contig_83252_1_209_-;contig_81913_1_209_-;contig_87730_1_213_-;contig_73109_1_213_-;contig_106541_1_213_+;contig_76155_1_223_+;contig_73747_1_221_+;contig_110265_1_221_+;contig_109983_1_226_-;contig_70127_1_308_+;contig_121352_1_397_-;contig_65372_1_453_+;THAOC_20341 EJK59439.1 hypothetical protein |THAOC_20341; contig_100917_1_82_-; |THAOC_20341 |NA       |NA        |NA  |contig_100917_1_82_-  |THAOC_20341           |        |NA      |
|            8|contig_101299_1_202_+                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            |contig_101299_1_202_+;             |NA          |NA       |NA        |NA  |contig_101299_1_202_+ |contig_101299_1_202_+ |        |NA      |
|            9|contig_102546_1_216_-                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            |contig_102546_1_216_-;             |NA          |NA       |NA        |NA  |contig_102546_1_216_- |contig_102546_1_216_- |        |NA      |
|           10|contig_103488_1_207_+;contig_120062_1_241_+                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      |contig_103488_1_207_+;             |NA          |NA       |NA        |NA  |contig_103488_1_207_+ |contig_103488_1_207_+ |        |NA      |

<a id="TO1005_B"></a>

##TO 1005_B

[Back Up](#BackUP)



```r
TO05_B_ProteinIDs <- TO05_B %>% 
  select(original_row, Protein.IDs)


TO05_B_ProteinIDs_info <- TO05_B_ProteinIDs %>% 
  mutate(THAOC = str_extract(Protein.IDs, "THAOC_\\d{5}")) %>% #uses RegEx to find first THAOC# in string and deposit result in new variable
  mutate (PS_start = str_extract(Protein.IDs, "^[a-z]{3}[A-Z]")) %>%  #deposits those on chloroplast genome e.g. PSI, PSII, ATP
  mutate (PS_middle = str_extract(Protein.IDs, ";[a-z]{3}[A-Z]")) %>%  #deposits those on chloroplast mutate (Rib = str_extract(Protein.IDs, "rp[l,s]\\d")) %>% #RegEx to find ribosonmal proteins %>% 
  mutate(PS_middle = str_extract(PS_middle, "[a-z]{3}[A-Z]")) %>%
  mutate (Rib = str_extract(Protein.IDs, "rp[l,s]\\d{1,2}")) %>% #RegEx to find ribosonmal proteins
  mutate (contig = str_extract(Protein.IDs,"^contig_\\d{1,6}_\\d{1,5}_\\d{2,5}_\\-|^contig_\\d{1,6}_\\d{1,5}_\\d{2,5}_\\+")) #getting all that START with a contig (as best hits)


TO05_B_ProteinIDs_info_comb <- unite(TO05_B_ProteinIDs_info, gene_ORF_comb, c(PS_start,PS_middle, Rib, THAOC, contig), sep = "; ", remove = FALSE)

TO05_B_ProteinIDs_info_comb <- TO05_B_ProteinIDs_info_comb %>% 
  mutate(colon = ";" )

TO05_B_ProteinIDs_info_comb <- unite(TO05_B_ProteinIDs_info_comb, together, c(gene_ORF_comb, colon), sep = "")

#now I would like to remove NAs
TO05_B_ProteinIDs_info_comb$together <- gsub("NA;", "", TO05_B_ProteinIDs_info_comb$together)
TO05_B_ProteinIDs_info_comb$together <- trimws(TO05_B_ProteinIDs_info_comb$together, which = "both")

#TO05_B_ProteinIDs_info_comb <- TO05_B_ProteinIDs_info_comb %>% 
#  mutate(gene_ORF_id = str_extract(together, "^THAOC_\\d{5}"))

TO05_B_ProteinIDs_info_comb <- TO05_B_ProteinIDs_info_comb %>% 
  mutate(gene_ORF_id_try = str_extract(together, "^.{4,27};"))
#this will lead to psbA etc and THAOC together
#so I will need to separate this column with using Sep=";"
TO05_B_ProteinIDs_info_comb <- separate(TO05_B_ProteinIDs_info_comb, gene_ORF_id_try, c("gene_ORF_id", "extra_1", "extra_2"), sep = ";")
```

```
## Warning: Too few values at 2439 locations: 1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
## 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, ...
```

```r
#now I would like to remove ";" - not needed anymore
#TO05_B_ProteinIDs_info_comb$gene_ORF_id <- gsub(";", "", TO05_B_ProteinIDs_info_comb$gene_ORF_id)

#and removing whitespaces
TO05_B_ProteinIDs_info_comb$gene_ORF_id <- trimws(TO05_B_ProteinIDs_info_comb$gene_ORF_id, which = "both") #leading white space "left" trailing white space "right"



write.table(TO05_B_ProteinIDs_info_comb, "Output_Data/2017_B_TO05_proteomics_gene_ORF_ids_forMaster.txt", sep="\t", row.names = FALSE, col.names = TRUE)

kable(head(TO05_B_ProteinIDs_info_comb), format = "markdown")
```



| original_row|Protein.IDs                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 |together                             |THAOC       |PS_start |PS_middle |Rib |contig                 |gene_ORF_id           |extra_1 |extra_2 |
|------------:|:---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|:------------------------------------|:-----------|:--------|:---------|:---|:----------------------|:---------------------|:-------|:-------|
|            1|atpA ADB27547.1 ATP synthase CF1 subunit alpha                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              |atpA;                                |NA          |atpA     |NA        |NA  |NA                     |atpA                  |        |NA      |
|            2|atpF ADB27545.1 ATP synthase CF0 subunit I B chain                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          |atpF;                                |NA          |atpF     |NA        |NA  |NA                     |atpF                  |        |NA      |
|            9|contig_100071_60_204_-;contig_83969_1_154_+;contig_99757_91_200_-;contig_90515_98_207_-;contig_96649_1_131_+;contig_85924_210_344_-;contig_107383_1_132_+;contig_95379_79_219_-;contig_85911_51_200_-;contig_87982_50_203_-;contig_99688_63_216_-;contig_86910_62_215_-;contig_80773_63_216_-;contig_78821_50_216_-;contig_94062_28_205_-;contig_76920_31_207_-;contig_74862_27_207_-;contig_96797_38_225_-;contig_86128_38_225_-;contig_80567_38_225_-;contig_106363_1_200_-;contig_89761_1_201_-;contig_115461_1_200_+;contig_104579_1_202_+;contig_86514_1_203_-;contig_81818_1_203_+;contig_121266_1_217_+;contig_113607_1_222_-;contig_130909_1_231_-;contig_70370_8_407_-;THAOC_17785 EJK61682.1 hypothetical protein |THAOC_17785; contig_100071_60_204_-; |THAOC_17785 |NA       |NA        |NA  |contig_100071_60_204_- |THAOC_17785           |        |NA      |
|           10|contig_100119_1_202_-;contig_106652_1_211_+;THAOC_25974 EJK54403.1 hypothetical protein                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     |THAOC_25974; contig_100119_1_202_-;  |THAOC_25974 |NA       |NA        |NA  |contig_100119_1_202_-  |THAOC_25974           |        |NA      |
|           11|contig_100206_1_201_+                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       |contig_100206_1_201_+;               |NA          |NA       |NA        |NA  |contig_100206_1_201_+  |contig_100206_1_201_+ |        |NA      |
|           12|contig_102056_1_211_-                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       |contig_102056_1_211_-;               |NA          |NA       |NA        |NA  |contig_102056_1_211_-  |contig_102056_1_211_- |        |NA      |
