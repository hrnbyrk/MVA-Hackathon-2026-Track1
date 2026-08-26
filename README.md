# MVA-Hackathon-2026-Track1
**Variant Prediction pipeline in R for the Rare Disease MVA Hackathon 2026**

This repository contains the R and Python hybrid code pipeline developed for Track 1 of the "Rare Disease, Real Kid: MVA Hackathon 2026". The study integrates clinical phenotype mapping with targeted Whole Genome Sequencing (WGS) variant filtering to identify the causal driver mutations for an ultra-rare pediatric case.

---

### 🏆 Track 1 Submission

#### Methods Writeup
The proband's clinical features, including Rhabdomyosarcoma [HP:0002859], Nephrocalcinosis [HP:0000121], and Short stature [HP:0004322], strongly indicated Mosaic Variegated Aneuploidy (MVA) syndrome, a condition driven by chromosomal instability. 

To pinpoint the causal variant, we targeted known MVA-associated spindle assembly checkpoint genes (*BUB1B*, *CEP57*, *TRIP13*). We utilized the `VariantAnnotation` and `GenomicRanges` packages in R/Bioconductor to perform targeted subsetting of the 85GB WGS VCF file against hg38 reference coordinates. This approach efficiently reduced the dataset to 21 region-specific variants. 

Genotype (GT) isolation successfully highlighted several homozygous (1/1) INDELs in *BUB1B* and *TRIP13*. We prioritized these homozygous INDELs (e.g., rs57676380, rs200074915) due to their high likelihood of causing frameshifts, triggering nonsense-mediated decay, or producing truncated, non-functional spindle-checkpoint proteins.

#### Ranked Variant Predictions
1. **rs57676380** (Chromosome 15, *BUB1B* gene) - Homozygous (1/1) Deletion (TA -> T)
2. **rs200074915** (Chromosome 15, *BUB1B* gene) - Homozygous (1/1) Deletion (TA -> T)
3. **rs57347691** (Chromosome 15, *BUB1B* gene) - Homozygous (1/1) Insertion (G -> GT)
4. **rs33953938** (Chromosome 11, *TRIP13* gene) - Homozygous (1/1) Insertion (T -> TTA)
5. **rs5793747** (Chromosome 11, *TRIP13* gene) - Homozygous (1/1) Insertion (C -> CT)
6. **rs4924431** (Chromosome 15, *BUB1B* gene) - Homozygous (1/1) SNP (G -> A)

---

### 📂 Datasets Used
*   **Source:** Hugging Face gated repository (`SageBio/mva-hackathon-2026-data`).
*   **Size:** Approximately 85 GB total dataset size.
*   **Key Files Processed:** 
    *   `WGS_EX2312012_HGWCNDSX7.vcf.gz` (WGS Variant Call Format, ~315 MB)
    *   `WGS_EX2312012_HGWCNDSX7.vcf.gz.tbi` (VCF Index)
    *   `Challenge_Clinical_Phenotype_1.docx` (Clinical references)

---

### 🔬 Methodology
1. **Data Retrieval Pipeline:** Utilized the `huggingface_hub` Python module via R (`reticulate`) to pass user authentication tokens and securely download gated clinical and genomic files.
2. **Genomic Subsetting:** Used `TxDb.Hsapiens.UCSC.hg38.knownGene` to extract precise start/end coordinates for the target genes (*BUB1B*, *CEP57*, *TRIP13*).
3. **Targeted Filtering:** Applied `ScanVcfParam` to strictly read the genetic loci of target genes, dynamically converting chromosome nomenclature between UCSC (`chr1`) and NCBI/Ensembl (`1`) formats.
4. **Prioritization:** Filtered the subset for homozygous (1/1) INDELs (insertions and deletions) likely to cause frameshifts and trigger nonsense-mediated decay.

---

### 💻 System Prerequisites
*   **R & RStudio:** R version 4.4.x or newer.
*   **RTools:** Required for building Bioconductor packages from source.
*   **Python/Miniconda:** Handled automatically by the script via the `reticulate` package.
*   **Hugging Face Account:** A valid Access Token with permissions to the `SageBio/mva-hackathon-2026-data` gated repository.

#### R Package Installation
```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
BiocManager::install(c("VariantAnnotation", "TxDb.Hsapiens.UCSC.hg38.knownGene", "org.Hs.eg.db", "GenomicRanges"), ask = FALSE)
install.packages("reticulate")
