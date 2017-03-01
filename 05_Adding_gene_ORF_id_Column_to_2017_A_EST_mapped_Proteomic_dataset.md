# 5_Adding gene_ORF_id column to EST matched Proteomic Data
Anna A. Hippmann  
March 1, 2017  

I got the new dataset (2017_A_TO03_master_...). Old peptides have been re-searched against the combined protein sequences of TO05 genome (~32K) + my own generates EST library (~ 145K contig ORFs). THe contigs hit either the Chloroplast, Gene ORFs, somewhere else on the scaffolds or nothing. The protein hit column includes all proteins that have been hit by the LC-MS/MS peptides. This results in often including identifiers for many contigs as well as THAOC_ and/or Chloroplast encoded proteins.

I would like to compare the results with the original ones. FOr that, I would like to create a "gene_ORF_id" column where each row is associated with only one identifier. I will decide on the identifyer based on the following criteria:

1) if a chloroplast protein identifier is part of the list, use THAT ONE
2) if THAOC identifyer(s) are part, use the first one
3) if only contig identifiers are used, use the first one





<a id="BackUP"></a>

And I do this for both

[TO 1003](#TO1003)

[TO 1005](#TO1005)


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

TO03_ProteinIDs_info_comb <- TO03_ProteinIDs_info_comb %>% 
  mutate(gene_ORF_id = str_extract(together, "^THAOC_\\d{5}"))

TO03_ProteinIDs_info_comb <- TO03_ProteinIDs_info_comb %>% 
  mutate(gene_ORF_id = str_extract(together, "^.{4,27};"))

#now I would like to remove ;
TO03_ProteinIDs_info_comb$gene_ORF_id <- gsub(";", "", TO03_ProteinIDs_info_comb$gene_ORF_id)


write.table(TO03_ProteinIDs_info_comb, "Output_Data/2017_A_TO03_proteomics_gene_ORF_ids_forMaster.txt", sep="\t", row.names = FALSE, col.names = TRUE)
```

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

TO05_ProteinIDs_info_comb <- TO05_ProteinIDs_info_comb %>% 
  mutate(gene_ORF_id = str_extract(together, "^THAOC_\\d{5}"))

TO05_ProteinIDs_info_comb <- TO05_ProteinIDs_info_comb %>% 
  mutate(gene_ORF_id = str_extract(together, "^.{4,27};"))

#now I would like to remove ;
TO05_ProteinIDs_info_comb$gene_ORF_id <- gsub(";", "", TO05_ProteinIDs_info_comb$gene_ORF_id)


write.table(TO05_ProteinIDs_info_comb, "Output_Data/2017_A_TO05_proteomics_gene_ORF_ids_forMaster.txt", sep="\t", row.names = FALSE, col.names = TRUE)
```


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
|TO03   |         1244|     1244|          31|       21|   936|    727|        1224|
|TO05   |         2480|     2480|          38|       31|  1883|   1380|        2430|

```r
write.table(Sum_Table,"Output_Data/2017_A_TO03_TO05_geneORF_id_Count_Table.txt", sep="\t", row.names = FALSE, col.names = TRUE)
```

