# Dissecting genetic variance structure and evaluating genomic prediction models for single-cross hybrids derived from Stiff Stalk and Non-Stiff Stalk maize populations

**Godoy dos Santos, J.C., Edwards, J., Lee, E., Mikel, M.A., Fernandes, S.B., Hirsch, C., Lipka, A.E., & Bohn, M.O.**  
(Submitted to *Theoretical and Applied Genetics*, 2025)

---

### Repository description

This repository is organized as follows:

- **`data/`** – contains the processed data sets used in the analyses:  
  - BLUEs from single- and multi-environment models (`single_env_analyses/` and `multi_env_analyses/`);  
  - Additive and dominance genomic relationship matrices considering all individuals combined, as well as separated by maturity group (`genomic_matrices/`).  

- **`scripts/`** – includes the R scripts used for:  
  - estimation of genetic variance components;
  - genomic prediction of single-cross hybrids under different cross-validation scenarios.
    
The **`scripts/`** folder also contains a helper script (`0_functions.R`) with user-defined functions that perform key analytical steps, including:  
    - (1) `fit_design2_models()` — fits GBLUP-based multi-kernel models for variance component estimation and prediction;  
    - (2) `get_covariance_matrices()` — computes genetic covariances between tested and untested hybrids;  
    - (3) `build_s_matrix()` — constructs the S matrix from the additive relationships of the parental groups. 

Note:
In all data files,
Pedigree1 refers to the seed parents (Stiff Stalk lines, SS),
Pedigree2 refers to the pollen parents (Non-Stiff Stalk lines, NSS), and
Pedigree identifies the single-cross hybrids derived from their combinations.

---

### Raw data

Phenotypic and genotypic data were obtained from the **Genomes to Fields (G2F)** initiative and are publicly accessible through the following links:

#### **Phenotypic data**
- **2016 season**  
  - File name: `g2f_2016_hybrid_data_clean.csv`  
  - Source: [G2F Planting Season 2016 v2](https://datacommons.cyverse.org/browse/iplant/home/shared/commons_repo/curated/GenomesToFields_2014_2017_v1/G2F_Planting_Season_2016_v2/a._2016_hybrid_phenotypic_data)

- **2017 season**  
  - File name: `g2f_2017_hybrid_data_clean.csv`  
  - Source: [G2F Planting Season 2017 v2](https://datacommons.cyverse.org/browse/iplant/home/shared/commons_repo/curated/GenomesToFields_2014_2017_v1/G2F_Planting_Season_2017_v2/a._2017_hybrid_phenotypic_data)

#### **Genotypic data**
- **Inbred genotypic panel (2014–2023)**  
  - File name: `inbreds_G2F_2014-2023_437k.vcf`  
  - Source: [G2F Inbred Genotypic Data 2014–2023 (437k SNPs)](https://datacommons.cyverse.org/browse/iplant/home/shared/commons_repo/curated/GenomesToFields_G2F_genotypic_data_2014_to_2023/inbreds_G2F_2014-2023_437k.vcf)

---
Any questions about the analyses, please, contact me!

Jenifer Camila

Email: jcg92@illinois.edu
