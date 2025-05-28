library(readr)
library(tidyverse)

###
### Construct native range Z and D matrices
###

# All 307 samples, both native and invaded 
native <- read_csv(here('01_data', 'data', 'bioclim_and_other', 'BRTEclim.csv')) %>% filter(range == 'native')

# IBS.id column contains snp matrix columns (starting from the 3rd col)
cols <- read_csv(here('01_data', 'data', 'bioclim_and_other', '307tips.csv'))
native <- left_join(native, cols, by = 'PopNum') %>% filter(!is.na(Longitude)) 
native <- native %>% arrange(IBS.id) %>% mutate(snp_map_col = IBS.id + 3) # + 3 to account for the first 3 non-genotype columns

# Extract native genotypes
all_SNPs <- as.data.frame(read.table(here('01_data', 'data', 'snps', 'BRTE307_LDfilteredSNPs.bed'), 
                                     header = FALSE, sep =",", stringsAsFactors=FALSE))

# Extract samples for native genotypes
native_snps <- all_SNPs[,c(1:3,native$snp_map_col)] 

# Sample 300 snps/chrom

chroms <- c(1:7)
bin_num <- 1
num_snps_each_bin <- 300

# Extract the SNP name and create a new df with the chromosome number and base pair number
snp_labels <- native_snps[,1]
snp_data <- data.frame(snp_name = snp_labels) %>%
  mutate(chrom = as.numeric(str_split_i(snp_name, "_", i = 2)),
         bp = as.numeric(str_split_i(snp_name, "_", i = 3)))

chrom_lengths <- read.table(here('01_data', 'data', 'bioclim_and_other', 'BRTE_chr_bp.txt')) %>%
  rename(chr = V1, length = V2) %>%
  mutate(chunk_length = length %/% bin_num)

# Randomly sample 2100 total loci
# Use seed for reproducibility
set.seed(2) 
snp_samps <- c()

for(i in chroms){
  # Get unique chunk_size and initialize boundaries 
  chunk_size <- chrom_lengths$chunk_length[chrom_lengths$chr == paste0("chr", i)]
  beg <- 0 
  end <- chunk_size 
  chrom_i <- snp_data %>% filter(chrom == i)
  
  # Chunk snps by chrom, and sample one locus from each chunk 
  for(j in 1:bin_num){
    
    chrom_i_samp <- chrom_i %>%
      filter(bp >= beg & bp <= end) %>%
      dplyr::select(snp_name) %>% 
      sample_n(size = num_snps_each_bin, replace = FALSE) %>%
      pull(snp_name)
    
    snp_samps <- c(snp_samps, chrom_i_samp)
    
    # Update boundaries 
    beg <- beg + chunk_size
    end <- end + chunk_size
  }
}

# Extract indices of sampled loci
keep_idx <- which(native_snps$V1 %in% snp_samps)

# Find genetic distances 
native_Z <- native_snps[keep_idx, -c(1:3)]
locus <- t(native_Z) # select 1:n_u rows, and unselect first 3 columns
Mydata <- df2genind(locus, ploidy = 2, loc.names = native_snps[keep_idx,1], ind.names = native$NewSiteCode, sep=" ")
Mydata2 <- genind2loci(Mydata)
distgenDIFF <- dist.gene(Mydata2, method="pairwise")
native_D <- as.matrix(distgenDIFF)
max(distgenDIFF) # 1017

# Visualize
hist(distgenDIFF)

# Save 
save_path <- here('01_data', 'data', 'snps', 'native_Z_and_D.RData')
save(file = save_path, native_D, native_Z)
