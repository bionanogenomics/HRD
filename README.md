# Homologous Recombination Deficiency Script #
This script is designed to calculate a homologous recombination deficiency (HRD) score. The HRD score is the summation of the three HRD signatures:

1. Loss of heterozygosity (HRD-LOH): the number of regions with a loss  >15 Mbp but shorter than the whole chromosome
2. Telomeric allelic imbalance (TAI): the number of regions of gain and loss >10Mbp that extend to a subtelomere but do not cross the centromere
3. Large-scale state transition LST: the number of chromosomal breakpoints whose SV size >10Mb but not the whole chromosome


### Running Perl Script ###
---
Current execution command
```
perl calculate_HRD_LST_LOH_TAI_score_smap_cnv_3.7.pl --smapFile --cnvFile --aneuploidyFile --outFile --centromereFile --chromLengthFile --cnvMaskFile
```

** Required Files **
```
-smapFile:			Annotated_SV.smap 
-cnvFile:			CNV.txt 
-aneuploidyFile:	Aneuploidy.txt 
-outFile: 			output_file_prefix.txt 
-centromereFile: 	hg38_centromeres.txt (or other appropriate build)
-chromLengthFile: 	hg38_chrom_length.txt (or other appropriate build)
-cnvMaskFile: 		hg38_cnv_masks.bed (or other appropriate build)
```

You can find the files for centromere, chromosome length, and cnv mask for hg38 and hg19 in this repo (*hg38_input_files* and *hg19_input_files* folders). The other files will need to be obtained from a **Rare Variant Analysis**.

If running the RVA on the command line, the files can be found in the output here:
```
smapFile: output/data/annotation/variants_combine_filters_inMoleRefine1.smap
cnvFile: output/data/alignmolvref/copynumber/cnv_calls_exp.txt
aneuploidyFile: output/data/alignmolvref/copynumber/cnv_chrAneuploidy_exp.txt
```

If downloading the files from a RVA run on Access, you will want to SHOW ALL variants in Access before downloading the smap, cnv and aneuploidy files.


You should also make sure to have a run log to check for errors.
```
&> run_log.txt
```


** Execution Command Example **
```
perl HRD_project/calculate_HRD_LST_LOH_TAI_score_smap_cnv_3.7.pl -smapFile HRD_project/bnx_files/CLL_7_-_Rare_Variant_Analysis_9_23_2022_17_57_16_Annotated_SV.smap -cnvFile HRD_project/bnx_files/CLL_7_-_Rare_Variant_Analysis_9_23_2022_17_57_16_CNV.txt -aneuploidyFile HRD_project/bnx_files/CLL_7_-_Rare_Variant_Analysis_9_23_2022_17_57_16_Aneuploidy.txt -outFile HRD_project/HRD_scores/CLL_7_HRD_scores.txt -centromereFile HRD_project/hg38_centromeres.txt -chromLengthFile HRD_project/chrom_length_only.txt -cnvMaskFile HRD_project/hg38_cnv_masks.bed &> HRD_project/logs/CLL_7_run_log.txt
```


### Output Example ###
---
Below is what an output file will look like. The first column is the chromosome, the next three columns are the measurments of HRD outlined above, and the final column is the number of intrachromosomal fusions.
```
chrom	LST	LOH	TAI	numIntrachromFusion_chromothripsis
1	5	1	0	4
2	1	0	1	0
3	3	0	1	0
4	0	0	0	0
5	0	0	0	0
6	11	2	1	3
7	1	0	0	0
8	4	1	1	0
9	1	0	0	1
10	0	0	0	0
11	1	0	1	0
12	2	0	1	1
13	5	1	0	2
14	1	0	0	0
15	0	0	0	0
16	1	1	0	0
17	5	1	1	0
18	0	0	0	0
19	0	0	0	0
20	0	0	0	0
21	0	0	0	0
22	0	0	0	0
23	4	0	1	0
24	0	0	0	0
```


### Publication ###
---

The paper that the script was written for is currently being peer-reviewed

### Follow Ups & Questions ###
---
The HRD script was written by Andy Pang, and this README and files have been uploaded by Kelsea Chang. For any questions, please reach out to either party.