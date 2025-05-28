library(dplyr)
library(sf)
library(stringr)
library(here)

###
### Construct genetic distance matrices (Z(C) and D(C))
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

# Read in SNP data
SNPs <- as.data.frame(read.table(here('01_data', 'data', 'snps', 'BRTEcg_SNPs.bed'), header = FALSE, sep=",",stringsAsFactors=FALSE))

# For D, take just the cg observations (for D(C))
snp_cols <- genotype_codes$SNPmatrix_column 

# Extract snp_cols 
SNPs_1 <- SNPs[,c(1:3,snp_cols)] ## Columns 1:3 provide line name, reference and alternate allele.

# randomly sample by chromosome (2100 loci total)
chroms <- c(1:7)
bin_num <- 1
num_snps_each_bin <- 300
total_snps <- bin_num*num_snps_each_bin*7

snp_labels <- SNPs_1[,1]

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

keep_idx <- which(SNPs_1$V1 %in% snp_samps)

# Extract the Z matrix for the CG sites (Z(CG))
invaded_Z <- SNPs_1[keep_idx, -c(1:3)]
locus <- t(invaded_Z)
Mydata <- df2genind(locus, ploidy = 2, loc.names = SNPs_1[keep_idx,1], ind.names = genotype_codes$source, sep="")

# Calculate allele sharing distance for CG sites
Mydata2 <- genind2loci(Mydata)
distgenDIFF <- dist.gene(Mydata2, method="pairwise")
invaded_D <- as.matrix(distgenDIFF)

# Visualize genetic distances 
# hist(distgenDIFF)

save_path <- here('01_data', 'data', 'snps','invaded_Z_and_D.RData')
save(file = save_path, invaded_Z, invaded_D)
