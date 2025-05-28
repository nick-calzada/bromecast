library(dplyr)
library(sf)
library(stringr)
library(geosphere)

###
### Construct augmented genetic distance matrices 
###

# Load in locations and codes
genotype_locations <- 
  read_csv(here('01_data', 'data', 'bioclim_and_other', 'BioclimateOfOrigin_AllGenotypes.csv')) %>%
  dplyr::select("source" = site_code, genotype, "longitude" = lon, "latitude" = lat)

# Join with key matrix for connecting with genotype matrix 
genotype_codes <- 
  read_csv(here('01_data', 'data', 'bioclim_and_other', 'BRTEcg_genotypesCode.csv')) %>%
  left_join(genotype_locations)  %>%
  filter(!is.na(longitude)) %>%
  arrange(SNPmatrix_column) %>%
  st_as_sf(.,
           coords = c("longitude", "latitude"), 
           crs = 4326) %>%
  ### Keep lat, long as columns too 
  dplyr::mutate(longitude = sf::st_coordinates(.)[,1],
                latitude = sf::st_coordinates(.)[,2])

# Gather assigned satellite sites
assigned <- read_csv(here('01_data', 'data', 'bioclim_and_other', 'assigned_genotypes.csv')) %>% 
  mutate(dist = sqrt(Lat_Diff ** 2 + Lon_Diff ** 2))

# Gather indices of smallest distances
idx <- which(!is.na(assigned$genotype) & assigned$dist < 0.05) # make sure same order as in assigned_genotypes.csv
filt_assigned <- assigned[idx, ]

# Join with SNP mapping column
joined <- left_join(filt_assigned, genotype_codes, by = 'genotype') %>% rename(source=source.x) 

# Read in SNP data
SNPs <- as.data.frame(read.table(here('01_data', 'data', 'snps', 'BRTEcg_SNPs.bed'), header = FALSE, sep=",",stringsAsFactors=FALSE))

# For aug_mat, concatenate cg and assigned sat site cols
snp_cols_aug <- c(genotype_codes$SNPmatrix_column, joined$SNPmatrix_column) # for the augmented matrix

# Extract snp_cols 
SNPs_aug <- SNPs[,c(1:3,snp_cols_aug)]

# Randomly sample by chromosome (2100 loci total)
chroms <- c(1:7)
bin_num <- 1
num_snps_each_bin <- 300
total_snps <- bin_num*num_snps_each_bin*7

snp_labels <- SNPs_aug[,1]

snp_data <- data.frame(snp_name = snp_labels) %>%
  mutate(chrom = as.numeric(str_split_i(snp_name, "_", i = 2)),
         bp = as.numeric(str_split_i(snp_name, "_", i = 3)))

chrom_lengths <- read.table(here('01_data', 'data', 'bioclim_and_other','BRTE_chr_bp.txt')) %>%
  rename(chr = V1, length = V2) %>%
  mutate(chunk_length = length %/% bin_num)


set.seed(2) # Seed used for april_25_aug_mat...RData
snp_samps <- c()
for(i in chroms){
  
  print(glue('Sampling from chromosome {i}...'))
  chunk_size <- chrom_lengths$chunk_length[chrom_lengths$chr == paste0("chr", i)]
  
  beg <- 0 
  end <- chunk_size 
  chrom_i <- snp_data %>% filter(chrom == i)
  
  # Chunk loci by chromosome, and sample one locus from each chunk 
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

# Extract the augmented Z matrix
keep_idx <- which(SNPs_aug$V1 %in% snp_samps)
Z_aug <- SNPs_aug[keep_idx, -c(1:3)]
locus_aug <- t(Z_aug)
Mydata_aug <- df2genind(locus_aug, ploidy = 2, loc.names = SNPs_aug[keep_idx,1], ind.names = c(genotype_codes$source,joined$source), sep="")

# Calculate allele sharing distance for augmented matrix
Mydata2_aug <- genind2loci(Mydata_aug)
distgenDIFF_aug <- dist.gene(Mydata2_aug, method="pairwise")
D_aug <- as.matrix(distgenDIFF_aug)

# Visualize genetic distances 
hist(distgenDIFF_aug)

# For the augmented matrix, rename row and column names: 
row_col_names <- c(genotype_codes$source, joined$source)
rownames(D_aug) <- row_col_names
colnames(D_aug) <- row_col_names

# Save
save_path <- here('01_data', 'data', 'snps', 'aug_invaded_Z_and_D.RData')
save(file = save_path, Z_aug, D_aug)
