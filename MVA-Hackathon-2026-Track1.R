# =========================================================================
# RARE DISEASE, REAL KID: MVA HACKATHON 2026 - TRACK 1 PIPELINE
# Description: Targeted WGS VCF filtering and genotype extraction
# Target Genes: BUB1B, CEP57, TRIP13 (MVA Syndrome Panel)
# =========================================================================

# 1. ENVIRONMENT SETUP & PACKAGE INSTALLATION -----------------------------
# Install BiocManager if missing
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install necessary Bioconductor packages silently
BiocManager::install(c("VariantAnnotation", "TxDb.Hsapiens.UCSC.hg38.knownGene", 
                       "org.Hs.eg.db", "GenomicRanges"), ask = FALSE)

# Install reticulate for Python integration
if (!requireNamespace("reticulate", quietly = TRUE)) {
  install.packages("reticulate")
}
library(reticulate)

# 2. HUGGING FACE DATA DOWNLOAD -------------------------------------------
# Set up isolated Miniconda environment to avoid OS path conflicts
# install_miniconda() # Run once if miniconda is not installed
use_condaenv("r-reticulate", required = TRUE)
# py_install("huggingface_hub") # Run once to install the module

hf <- import("huggingface_hub")
token <- "YOUR_HF_TOKEN" # NOTE: Replace with your actual Hugging Face token

print("Downloading dataset files...")
# Uncomment below lines to download files if not already downloaded
# hf$hf_hub_download(repo_id="SageBio/mva-hackathon-2026-data", filename="WGS_EX2312012_HGWCNDSX7.vcf.gz", repo_type="dataset", local_dir="./data", token=token)
# hf$hf_hub_download(repo_id="SageBio/mva-hackathon-2026-data", filename="WGS_EX2312012_HGWCNDSX7.vcf.gz.tbi", repo_type="dataset", local_dir="./data", token=token)

vcf_file <- "./data/WGS_EX2312012_HGWCNDSX7.vcf.gz"

# 3. GENOMIC SUBSETTING ---------------------------------------------------
library(VariantAnnotation)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(GenomicRanges)

# Define target genes associated with MVA
target_genes <- c("BUB1B", "CEP57", "TRIP13")

# Get Entrez IDs
gene_ids <- select(org.Hs.eg.db, keys=target_genes, keytype="SYMBOL", columns="ENTREZID")

# Extract coordinates from hg38 reference
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
gene_coords <- genes(txdb, filter=list(gene_id=gene_ids$ENTREZID))

# Harmonize chromosome naming (Translate UCSC 'chr1' to NCBI '1' format)
seqlevelsStyle(gene_coords) <- "NCBI"

# 4. TARGETED VCF READING -------------------------------------------------
# Create read parameters based on gene coordinates
param <- ScanVcfParam(which = gene_coords)

print("Filtering target gene regions from the VCF file...")
targeted_vcf <- readVcf(vcf_file, "hg38", param=param)

print(paste("Filtered Variant Count:", length(rowRanges(targeted_vcf))))

# 5. GENOTYPE (GT) EXTRACTION & PRIORITIZATION ----------------------------
# Extract proband genotype
genotypes <- geno(targeted_vcf)$GT

# Compile final variant table
detailed_table <- data.frame(
  Chromosome = seqnames(rowRanges(targeted_vcf)),
  Position = start(rowRanges(targeted_vcf)),
  ID = rownames(genotypes),
  REF = as.character(ref(targeted_vcf)),
  ALT = as.character(unlist(alt(targeted_vcf))),
  Quality = qual(targeted_vcf),
  Genotype = genotypes[,1] 
)

# Clean up row names
rownames(detailed_table) <- NULL

# Filter for alternate alleles (0/1 heterozygous or 1/1 homozygous)
suspicious_variants <- subset(detailed_table, Genotype != "0/0")

print("--- FULL LIST OF SUSPICIOUS VARIANTS IN THE PROBAND ---")

print(suspicious_variants)