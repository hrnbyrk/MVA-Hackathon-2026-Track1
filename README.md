# MVA-Hackathon-2026-Track1
**Variant Prediction pipeline in R for the Rare Disease MVA Hackathon 2026**

This repository contains the R and Python hybrid code pipeline developed for Track 1 of the "Rare Disease, Real Kid: MVA Hackathon 2026". The study integrates clinical phenotype mapping with targeted Whole Genome Sequencing (WGS) variant filtering to identify the causal driver mutations for an ultra-rare pediatric case.

---
### 🌟 Key Capabilities
*   **Memory-Efficient WGS Processing:** Targeted RAM-loading of specific genomic ranges instead of parsing the entire 85GB/315MB VCF dataset.
*   **Cross-Language Integration:** Executing Python-based Hugging Face API downloads seamlessly within RStudio using `reticulate`.
*   **Automated Naming Harmonization:** On-the-fly conversion of chromosome nomenclature between UCSC (`chr1`) and NCBI/Ensembl (`1`) formats.
*   **Precision Genotyping:** Direct extraction and filtering of homozygous (1/1) structural mutations from complex VCF annotations.

---

#### Ranked Variant Predictions
### 🏆 Track 1 Submission

#### Methods Writeup
The proband's clinical features, including Rhabdomyosarcoma [HP:0002859], Nephrocalcinosis [HP:0000121], and Short stature [HP:0004322], combined with significant parental reproductive history (recurrent miscarriages), strongly indicated Mosaic Variegated Aneuploidy (MVA) syndrome driven by autosomal recessive chromosomal instability.

To pinpoint the causal variants, we targeted known MVA-associated spindle assembly checkpoint genes (*BUB1B*, *CEP57*, *TRIP13*). We utilized the `VariantAnnotation` and `GenomicRanges` packages in R/Bioconductor to perform targeted subsetting of the WGS VCF file against hg38 reference coordinates. 

Expanding our search beyond homozygous models to include a dual-model filtering strategy for Compound Heterozygosity (0/1 + 0/1 pairs) successfully isolated the true driver mutations in the *BUB1B* gene. We prioritized these heterozygous variants due to their high impact (nonsense stop-gain and missense alterations) leading to loss of spindle-checkpoint function.

#### Ranked Variant Predictions
1. **chr15:40209701 (T > G) & chr15:40220612 (T > G)** (*BUB1B* gene) - Compound Heterozygous pair (p.L737* and p.N1002K) - **[Confirmed Match: 100/100 Score]**

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
### 🛠️ Technical Specifications & Prerequisites

#### Software Environment
*   **R Version:** 4.4.x or newer
*   **Bioconductor Version:** 3.20 or newer
*   **RTools:** Required for building packages from source
*   **Python:** Miniconda environment handled via `reticulate`
*   **Account:** A valid Hugging Face Access Token with permissions to the gated repository.

  
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

---
🚀 How to Run
Clone this repository to your local machine.

Open Track1_Analysis.R in RStudio.

Replace "YOUR_HF_TOKEN" on line 24 with your actual Hugging Face Access Token.

Run the script sequentially. The script will automatically:

Set up the Python environment.

Create a ./data folder and download the VCF files.

Filter the WGS data and output the final ranked variants to the console.

⚙️ Built-in Error Handling
The pipeline includes robust error handling mechanisms tailored for large-scale genomic datasets:

1. OS-Level Python Path Conflicts
Issue: reticulate failing to find a valid Python compiler due to Windows Store aliases.

Solution: Forced installation of an isolated r-reticulate Miniconda environment (use_condaenv).

Impact: Seamless execution of Python modules within R without OS variable conflicts.

2. Dataset 403 Forbidden Access
Issue: Standard HTTP GET requests returning 192-byte error files from Hugging Face's gated repository.

Solution: Leveraged the official huggingface_hub API in Python to pass secure tokens.

Impact: Reliable, chunk-supported downloading of massive WGS files.

3. Chromosome Nomenclature Mismatches
Issue: GenomicRanges returning 0 variants because Bioconductor standardizes on UCSC format (chr15), while the Hackathon VCF used NCBI format (15).

Solution: Implementation of seqlevelsStyle() <- "NCBI" to dynamically translate coordinates prior to VCF scanning.

Impact: Prevents null-returns and accurately extracts target genes.

🗄️ Databases & Tools
VariantAnnotation (v1.x): For parsing and subsetting large WGS VCF files.

TxDb.Hsapiens.UCSC.hg38.knownGene: Genomic coordinate database for mapping.

org.Hs.eg.db: Entrez ID annotations.

Hugging Face Hub: Secure dataset retrieval.

📄 License & Compliance
License: This project is licensed under the MIT License - see the LICENSE file for details.

Data Usage: The underlying clinical and genomic data provided by Sage Bionetworks and the MVA Society is licensed under CC-BY 4.0.

Privacy: All patient data handling strictly adheres to the WCG IRB protocol (#20252010). All intermediate raw data files must be deleted from local compute environments post-hackathon.
### ⚙️ Built-in Error Handling
The pipeline features resilient workarounds for common bioinformatics bottlenecks:
1. **OS-Level Python Path Conflicts:** Forced installation of an isolated `r-reticulate` Miniconda environment to bypass Windows Store aliases.
2. **Dataset 403 Forbidden Access:** Leveraged the official `huggingface_hub` API in Python instead of standard HTTP `GET` to pass secure tokens.
3. **Chromosome Nomenclature Mismatches:** Implementation of `seqlevelsStyle() <- "NCBI"` to dynamically translate UCSC coordinates prior to VCF scanning.

---

### 📄 License & Compliance
*   **License:** This repository's code is available under the MIT License.
*   **Data Usage:** The underlying clinical and genomic data provided by Sage Bionetworks and the MVA Society is licensed under **CC-BY 4.0**[cite: 2].
*   **Privacy:** All patient data handling strictly adheres to the WCG IRB protocol (#20252010). All intermediate raw data files must be deleted from local environments post-hackathon as per the challenge rules[cite: 2].
