# README - EST
Anna A. Hippmann  
January 18, 2016  

# EST_data
EST_Data_Analysis

This Repository will be the home of my EST analysis and all my related files!

Cheers
Anna

### 01)[Comparing EST expression with Proteomics Expression](https://github.com/ahippman/EST_data/blob/master/01_Comparing_EST_expression_vs_Proteomic_Expression.md) has lots of plots comparing

* edgeR results for EST vs TO03 or TO05
  + for both low Fe and low Cu
  + for all proteins, only soluble, only insoluble

* TO03 proteins vs TO05 proteins
  + for both low Fe and low Cu
  + for all proteins, only soluble, only insoluble



# Re-searched proteomics data

here, Jenny has used the original peptides and mapped them against a cobvined database of all 32K TO05 predicted proteins and the 145K TO03 EST predicted proteins

We have two sets for each TO03 and TO05

* set __2017_A__ was mapped using an updated software of MASCOT
* set __2017_B__ was mapped using the original software of MASCOT

### 05)[Adding gene_ORF_id column to EST matched Proteomic Data](https://github.com/ahippman/EST_data/blob/master/05_Adding_gene_ORF_id_Column_to_2017_A_EST_mapped_Proteomic_dataset.md)

the protein_ID column of my respective Master files gives the short proteinID of ALL contigs/ORFs that have been mapped by the peptides. SOmetimes it is only one protein ID, sometimes there are two or even up to >40!!! In order to getting a better grasp, I want to have only ONE identifier. I use the following ranking to make a desicion as to which identifyer of many I will use:

1) if htere is a chloroplast id, I take that one
2) if there is a THAOC id, I will take that one (if there are more than one, I will take the first one)
3) if there is only contig (EST) identifiers, I will take the first one of these

### 3) [3_Exploring_2017_01_proteomics_matched_against_EST_dataset.md](https://github.com/ahippman/EST_data/blob/master/3_Exploring_2017_01_proteomics_matched_against_EST_dataset.md)
here, I create a table in which I count the occurances of the differnet subsets of chloropalst / THAOC /contig hitting peptides and save also those that are only hitting contigs and the subsets of the contigs (i.e. chloorplast hitting, Ref gene hitting, Genome hitting, NOT MAPPING)
