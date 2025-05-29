# bromecast
Data processing and visualization for landscape genomics project of cheatgrass (*Bromus tectorum*). 

# Directory Stucture

```txt
bromecast
├── 01_data/                              
│   ├── bioclim_and_other/                # Bioclimatic and auxiliary data
│   │   ├── chelsav2/                     # Contains CHELSA V2 climate data. Can be found here: https://chelsa-climate.org/bioclim/
│   │   ├── 307tips.csv
│   │   ├── BRTE_chr_bp.txt               # Chromosome lengths
│   │   ├── BRTEcg_genotypesCode.csv
│   │   ├── BRTEclim.csv
│   │   ├── BioclimateOfOrigin_AllGenotypes.csv
│   │   ├── BromecastSites.csv
│   │   └── assigned_genotypes.csv
│   ├── snps/                           
│   │   ├── BRTE307_LDfilteredSNPS.bed  # Linkage disequilibrium-filtered SNP data of samples in both native and invaded ranges. 
│   │   ├── BRTEcg_SNPs.bed             # Linkage disequilibrium-filtered SNP data of common garden samples
│   │   ├── aug_invaded_Z_and_D.RData
│   │   └── invaded_Z_and_D.RData
│   └── processing_scripts/              
│       ├── invaded_range_bioclim_data_equal_dist.R
│       ├── make_invaded_Z_and_D.R
│       ├── make_invaded_augmented_Z_and_D.R
│       ├── make_native_Z_and_D.R
│       └── native_range_bioclim_data_equal_distance.R
├── 02_visualize/                         # Scripts for visualizing data
│   ├── invaded_range_gen_distance_plot.R
│   ├── native_range_gen_distance_plot.R
│   ├── plot_location_types.R
├── README.md                             
```
