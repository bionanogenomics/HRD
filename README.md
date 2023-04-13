![Bionano logo](images/Bionano-Logo.png?raw=true)

# Homologous Recombination Deficiency Script #

This script implements an automated algorithm to calculate a homologous recombination deficiency (HRD) score. The HRD score which is the summation of the three HRD signatures:
- Loss of heterozygosity (HRD-LOH): the number of regions with a loss  >15 Mbp but shorter than the whole chromosome
- Telomeric allelic imbalance (TAI): the number of regions of gain and loss >10Mbp that extend to a subtelomere but do not cross the centromere
- Large-scale state transition LST: the number of chromosomal breakpoints whose SV size >10Mb but not the whole chromosome


### Setup ###
---
Current Implementation
```
perl calculate_HRD_LST_LOH_TAI_score_smap_cnv_3.7.pl --smapFile --cnvFile --aneuploidyFile --outFile --centromereFile --chromLengthFile --cnvMaskFile
```

** Required Files **
```
-cnvFile: CNV.txt 
-aneuploidyFile: Aneuploidy.txt 
-outFile: output_file_prefix.txt 
-centromereFile: hg38_centromeres.txt (or other appropriate build)
-chromLengthFile: chrom_length_only.txt 
-cnvMaskFile: hg38_cnv_masks.bed (or other appropriate build)
```
