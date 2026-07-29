############################################################
### MITO-NUCLEAR ASSOCIATION ANALYSIS
############################################################

############################################################
### 1. LOAD PACKAGES
############################################################

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

pkgs <- c(
  "SNPRelate",
  "SeqArray",
  "GenomicRanges",
  "rtracklayer",
  "tidyverse",
  "lme4"
)

for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE))
    BiocManager::install(p)
}

library(SNPRelate)
library(SeqArray)
library(GenomicRanges)
library(rtracklayer)
library(tidyverse)
library(lme4)

############################################################
### 2. CONVERT VCF TO GDS (RUN ONCE)
############################################################

if (!file.exists("nuclear.gds")) {
  
  snpgdsVCF2GDS(
    vcf.fn = "mallard.CDS_biallelic_snps.vcf.gz",
    out.fn = "nuclear.gds",
    method = "biallelic.only"
  )
  
}

gds <- snpgdsOpen("nuclear.gds")

############################################################
### 3. LOAD METADATA
############################################################

mt <- read_csv(
  "mito_haplotypes.csv",
  show_col_types = FALSE
)

Q <- read_csv(
  "admixture_Qmatrix.csv",
  show_col_types = FALSE
)

meta <- read_csv(
  "sample_metadata.csv",
  show_col_types = FALSE
)

############################################################
### 4. LOAD MITOCHONDRIAL-FUNCTION GENE SET
############################################################
mito_gff <- read.table(
  "mallard_omics_numt_CDS.gff",
  sep = "\t",
  header = FALSE,
  comment.char = "#",
  stringsAsFactors = FALSE,
  fill = TRUE
)

if (ncol(mito_gff) == 9) {
  
  colnames(mito_gff) <- c(
    "chr",
    "feature",
    "start",
    "end",
    "score",
    "strand",
    "phase",
    "gene_id",
    "transcript_id"
  )
  
} else if (ncol(mito_gff) == 10) {
  
  colnames(mito_gff) <- c(
    "chr",
    "source",
    "feature",
    "start",
    "end",
    "score",
    "strand",
    "phase",
    "gene_id",
    "transcript_id"
  )
  
} else {
  
  stop("Unexpected GFF format")
  
}

mito_gff <- mito_gff %>%
  filter(feature == "CDS")

gene_ranges <- mito_gff %>%
  group_by(chr, gene_id) %>%
  summarise(
    start = min(start),
    end = max(end),
    .groups = "drop"
  )

gene_ranges$chr <- sub(
  "^chr",
  "",
  gene_ranges$chr
)

mito_gr <- GRanges(
  
  seqnames = gene_ranges$chr,
  
  ranges = IRanges(
    start = gene_ranges$start,
    end = gene_ranges$end
  ),
  
  gene_id = gene_ranges$gene_id
  
)

mito_ranges <- mito_gr

cat(
  "Number of mito genes:",
  length(mito_ranges),
  "\n"
)
############################################################
### 5. SNP INFORMATION
############################################################

snp_info <- tibble(
  
  snp_id = read.gdsn(
    index.gdsn(gds, "snp.id")
  ),
  
  chr = as.character(
    read.gdsn(
      index.gdsn(gds, "snp.chromosome")
    )
  ),
  
  pos = read.gdsn(
    index.gdsn(gds, "snp.position")
  )
  
)

sample_ids <- read.gdsn(
  index.gdsn(gds, "sample.id")
)

############################################################
### 6. SNP GRANGES
############################################################

snp_gr <- GRanges(
  
  seqnames = snp_info$chr,
  
  ranges = IRanges(
    start = snp_info$pos,
    width = 1
  ),
  
  snp_id = snp_info$snp_id
  
)

############################################################
### 7. ASSIGN MITO SNPs
############################################################

hits <- findOverlaps(
  snp_gr,
  mito_ranges
)

mito_snps <- unique(
  mcols(snp_gr)$snp_id[
    queryHits(hits)
  ]
)

bg_snps <- setdiff(
  mcols(snp_gr)$snp_id,
  mito_snps
)

############################################################
### 8. Assign mito-CDS SNPs
############################################################

hits <- findOverlaps(
  snp_gr,
  mito_ranges
)

mito_snps <- unique(
  mcols(snp_gr)$snp_id[
    queryHits(hits)
  ]
)

bg_snps <- setdiff(
  mcols(snp_gr)$snp_id,
  mito_snps
)

cat(
  "Mito SNPs:",
  length(mito_snps),
  "\n"
)

cat(
  "Background SNPs:",
  length(bg_snps),
  "\n"
)
############################################################
### 9. GENOTYPES
############################################################

mito_pruned <- unlist(
  snpgdsLDpruning(
    gds,
    snp.id = mito_snps,
    ld.threshold = 0.2
  ),
  use.names = FALSE
)

bg_pruned <- unlist(
  snpgdsLDpruning(
    gds,
    snp.id = bg_snps,
    ld.threshold = 0.2
  ),
  use.names = FALSE
)

cat(
  "Pruned mito SNPs:",
  length(mito_pruned),
  "\n"
)

cat(
  "Pruned background SNPs:",
  length(bg_pruned),
  "\n"
)
all_pruned <- unique(
  c(
    mito_pruned,
    bg_pruned
  )
)

geno_all <- snpgdsGetGeno(
  gds,
  snp.id = all_pruned
)

geno_all <- t(geno_all)

############################################################
### 10. PHENOTYPES
############################################################

pheno <- mt %>%
  inner_join(Q, by = "sample_id") %>%
  inner_join(meta, by = "sample_id") %>%
  mutate(
    mt_haplotype = factor(
      mt_haplotype,
      levels = c(
        "NWB",
        "OWA"
      )
    )
  )

sample_match <- match(
  pheno$sample_id,
  sample_ids
)

geno_all <- geno_all[, sample_match]

all_snp_ids <- all_pruned

############################################################
### 11. SNP association function
############################################################

test_snp <- function(g) {
  
  df <- pheno
  
  df$geno <- as.numeric(g)
  
  if (var(df$geno, na.rm = TRUE) == 0) {
    
    return(
      c(beta = NA, p = NA)
    )
    
  }
  
  fit <- tryCatch(
    
    lm(
      geno ~ mt_haplotype + Q2,
      data = df
    ),
    
    error = function(e) NULL
    
  )
  
  if (is.null(fit)) {
    
    return(
      c(beta = NA, p = NA)
    )
    
  }
  
  coefs <- coef(summary(fit))
  
  if (!"mt_haplotypeOWA" %in% rownames(coefs)) {
    
    return(
      c(beta = NA, p = NA)
    )
    
  }
  
  c(
    
    beta =
      coefs[
        "mt_haplotypeOWA",
        "Estimate"
      ],
    
    p =
      coefs[
        "mt_haplotypeOWA",
        "Pr(>|t|)"
      ]
    
  )
  
}
############################################################
### 12. Test mito SNPs
############################################################

mito_idx <- match(
  mito_pruned,
  all_snp_ids
)

mito_res <- t(
  
  sapply(
    
    mito_idx,
    
    function(i)
      test_snp(
        geno_all[i, ]
      )
    
  )
  
)

mito_results <- tibble(
  
  snp_id = mito_pruned,
  
  beta = mito_res[, "beta"],
  
  p = mito_res[, "p"]
  
) %>%
  
  mutate(
    FDR = p.adjust(
      p,
      method = "BH"
    )
  )

############################################################
### 13. REPEATED BACKGROUND SAMPLING
############################################################

obs_mito_sig <- sum(
  mito_results$p < 0.05,
  na.rm = TRUE
)

set.seed(123)

bg_sig_counts <- replicate(500, {
  
  bg_subset <- sample(
    bg_pruned,
    size = length(mito_pruned)
  )
  
  bg_idx <- match(
    bg_subset,
    all_snp_ids
  )
  
  bg_res <- t(
    
    sapply(
      bg_idx,
      function(i)
        test_snp(
          geno_all[i, ]
        )
    )
    
  )
  
  sum(
    bg_res[, "p"] < 0.05,
    na.rm = TRUE
  )
  
})

summary(bg_sig_counts)

mean(bg_sig_counts)

sd(bg_sig_counts)

empirical_p <- mean(
  bg_sig_counts >= obs_mito_sig
)

cat(
  "\nObserved mito significant SNPs:",
  obs_mito_sig,
  "\n"
)

cat(
  "Background mean:",
  mean(bg_sig_counts),
  "\n"
)

cat(
  "Background SD:",
  sd(bg_sig_counts),
  "\n"
)

cat(
  "Empirical P-value:",
  empirical_p,
  "\n"
)

hist(
  bg_sig_counts,
  breaks = 20,
  main = "Background Null Distribution",
  xlab = "Number of Significant SNPs"
)

abline(
  v = obs_mito_sig,
  col = "red",
  lwd = 3
)
############################################################
### 14. SINGLE BACKGROUND DRAW FOR FISHER TEST
############################################################

bg_subset <- sample(
  bg_pruned,
  size = length(mito_pruned)
)

bg_idx <- match(
  bg_subset,
  all_snp_ids
)

bg_res <- t(
  
  sapply(
    bg_idx,
    function(i)
      test_snp(
        geno_all[i, ]
      )
  )
  
)

bg_results <- tibble(
  snp_id = bg_subset,
  beta = bg_res[, "beta"],
  p = bg_res[, "p"]
)

############################################################
### 15. Fisher enrichment test
############################################################

n_mito_sig <- sum(
  mito_results$p < 0.05,
  na.rm = TRUE
)

n_bg_sig <- sum(
  bg_results$p < 0.05,
  na.rm = TRUE
)

enrichment_table <- matrix(
  
  c(
    n_mito_sig,
    nrow(mito_results) - n_mito_sig,
    
    n_bg_sig,
    nrow(bg_results) - n_bg_sig
  ),
  
  nrow = 2,
  byrow = TRUE
  
)

rownames(enrichment_table) <- c(
  "MitoGO",
  "Background"
)

colnames(enrichment_table) <- c(
  "Significant",
  "NotSignificant"
)

fisher_res <- fisher.test(
  enrichment_table
)

print(fisher_res)

############################################################
### 16. Candidate gene identification
############################################################

sig_snps <- mito_results %>%
  filter(FDR < 0.05)

if(nrow(sig_snps) > 0){
  
  sig_positions <- snp_info %>%
    filter(
      snp_id %in% sig_snps$snp_id
    )
  
  sig_gr <- GRanges(
    
    seqnames = sig_positions$chr,
    
    ranges = IRanges(
      start = sig_positions$pos,
      width = 1
    ),
    
    snp_id = sig_positions$snp_id
    
  )
  
  gene_hits <- findOverlaps(
    sig_gr,
    mito_gr
  )
  
  candidate_genes <- tibble(
    
    snp_id = mcols(sig_gr)$snp_id[
      queryHits(gene_hits)
    ],
    
    gene_id = mcols(mito_gr)$gene_id[
      subjectHits(gene_hits)
    ]
    
  )
  
  print(candidate_genes)
  
  write_csv(
    candidate_genes,
    "candidate_genes.csv"
  )
  
}

############################################################
### 16. Top nominal candidate genes
############################################################

top_hits <- mito_results %>%
  filter(!is.na(p)) %>%
  arrange(p) %>%
  head(20)

print(top_hits)

top_positions <- top_hits %>%
  left_join(
    snp_info,
    by = "snp_id"
  )

top_gr <- GRanges(
  
  seqnames = top_positions$chr,
  
  ranges = IRanges(
    start = top_positions$pos,
    width = 1
  ),
  
  snp_id = top_positions$snp_id
  
)

gene_hits <- findOverlaps(
  top_gr,
  mito_gr
)

candidate_genes <- tibble(
  
  snp_id =
    mcols(top_gr)$snp_id[
      queryHits(gene_hits)
    ],
  
  gene_id =
    mcols(mito_gr)$gene_id[
      subjectHits(gene_hits)
    ]
  
)

candidate_genes <- candidate_genes %>%
  left_join(
    top_hits,
    by = "snp_id"
  )

print(candidate_genes)

write_csv(
  candidate_genes,
  "top_candidate_genes.csv"
)

############################################################
### 18. SUMMARY
############################################################

cat(
  "\nObserved mito hits:",
  obs_mito_sig,
  "\n"
)

cat(
  "Background mean:",
  mean(bg_sig_counts),
  "\n"
)

cat(
  "Empirical P:",
  empirical_p,
  "\n"
)

print(fisher_res)

############################################################
### 19. CLOSE GDS
############################################################

snpgdsClose(gds)
