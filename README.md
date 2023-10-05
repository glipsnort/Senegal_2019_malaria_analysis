# Senegal_2019_malaria_analysis
Key analysis scripts for paper analyzing 2019 falciparum data from Senegal ("Malaria surveillance reveals parasite relatedness, signatures of selection, and correlates of transmission across Senegal"). Also included are two files for handling sample names. One file, all_sampdupes.txt, lists samples that are replicate sequencing of existing samples; these were dropped from analysis. The other,all_sampmap.txt, enables the user to translate the sample names as uploaded to the SRA to the sample names used for the paper; the latter names contain the three-letter abbreviations for the sample sites, which are important for analysis.

The scripts are quite specific to this dataset and the analysis in the paper and not at all pretty. They work from a vcf file containing the dataset and assume a directory structure with the following subdirectories under the directory containing the script: /data, /output, /seq, /results. They assume that whole genome sequence data is in a file named data/all.vcf.gz and that barcode data is in three files: data/Senegal_AllBarcodes_20200701.txt, data/mono_barcodes_filtered_2019.csv, and data/poly_barcodes_filtered_2019.csv. (These last three are included in this repository).                                                                                                                          
- `extract_het.pl all`    Calculates the number of het calls for each sample.                                                   
- `hetrate.py all yes`    Classifies samples as good/bad based on coverage and as mono/poly based on het rate.                   
- `extract_vcf.pl all 5`  Extracts genotypes for good mono samples for input to hmmIBD.                                        
- `hmmIBD -i seq/all_seq.txt -o output/all -m 20 -n 40`  Runs hmmIBD (which must be downloaded and compiled separately). -m and -n values can be lowered to, say, 10 with little loss of accuracy.
- `filter_ibd.py all`     Extracts hmmIBD results for related pairs only                                                          
- `list_sequenced_clones.py`  Creates list of clones to skip for some analyses                                                                               
- `count_samples.py`                                              
- `anal_ibd.py` Calculate partial and clonal relatedness within and between sites.
              
