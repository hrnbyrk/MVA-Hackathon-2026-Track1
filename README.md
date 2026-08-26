# MVA-Hackathon-2026-Track1
**Variant Prediction pipeline in R for the Rare Disease MVA Hackathon 2026**

### Methods Writeup
The proband's clinical features, including Rhabdomyosarcoma [HP:0002859], Nephrocalcinosis [HP:0000121], and Short stature [HP:0004322], strongly indicated Mosaic Variegated Aneuploidy (MVA) syndrome, a condition driven by chromosomal instability. 

To pinpoint the causal variant, we targeted known MVA-associated spindle assembly checkpoint genes (*BUB1B*, *CEP57*, *TRIP13*). We utilized the `VariantAnnotation` and `GenomicRanges` packages in R/Bioconductor to perform targeted subsetting of the WGS VCF file against hg38 reference coordinates. This approach efficiently reduced the dataset to 21 region-specific variants. 

Genotype (GT) isolation successfully highlighted several homozygous (1/1) INDELs in *BUB1B* and *TRIP13*. We prioritized these homozygous INDELs (e.g., rs57676380, rs200074915) due to their high likelihood of causing frameshifts, triggering nonsense-mediated decay, or producing truncated, non-functional spindle-checkpoint proteins.

### Ranked Variant Predictions
1. **rs57676380** (Chromosome 15, *BUB1B* gene) - Homozygous (1/1) Deletion (TA -> T)
2. **rs200074915** (Chromosome 15, *BUB1B* gene) - Homozygous (1/1) Deletion (TA -> T)
3. **rs57347691** (Chromosome 15, *BUB1B* gene) - Homozygous (1/1) Insertion (G -> GT)
4. **rs33953938** (Chromosome 11, *TRIP13* gene) - Homozygous (1/1) Insertion (T -> TTA)
5. **rs5793747** (Chromosome 11, *TRIP13* gene) - Homozygous (1/1) Insertion (C -> CT)
6. **rs4924431** (Chromosome 15, *BUB1B* gene) - Homozygous (1/1) SNP (G -> A)
