# A network medicine framework for positioning Myalgic Encephalomyelitis/Chronic Fatigue Syndrome (ME/CFS) within the human disease landscape

## Abstract

Understanding where ME/CFS sits among common diseases is a key step toward understanding its pathological mechanisms. CompareME is an R-based computational pipeline that constructs protein–protein interaction (PPI) networks for approximately 30 human diseases and systematically compares them using three orthogonal similarity metrics: gene-level overlap (Jaccard index), functional enrichment correlation (ORA correlation), and network topology (network separation, SAB). For each disease, candidate genes are retrieved from two evidence streams: GWAS credible sets (via the Open Targets platform, Locus-to-Gene scoring) and rare-variant studies (ClinGen curated genes and gene-burden tests). For ME/CFS specifically, a custom gene module is constructed by merging gene lists from two recent transcriptomic and proteomic studies (Zhang *et al.* 2025; Sardell *et al.* 2025). All pairwise similarity scores are benchmarked against a null distribution of 2,000 size-matched random disease modules, and statistical significance is evaluated empirically. The result is a quantitative map of pairwise disease similarity, with ME/CFS as a focal point.

## Methods

### Data source

I selected 28 common diseases, including neurological, psychiatric, metabolic, cardiovascular, inflammatory, and autoimmune conditions. The complete list, with full names, abbreviations, and identifiers (EFO or MONDO codes), is in Table 1. The list of diseases is passed to the script through [mydiseases.yml](main/mydiseases.yml).
ME/CFS is handled separately with the identifier `customCFS` (see below).

| Disease Full Name                        | Abbreviation | ID            |
| ---------------------------------------- | ------------ | ------------- |
| Alzheimer disease                        | AD           | MONDO_0004975 |
| Anxiety disorder                         | ANX          | EFO_0006788   |
| Arteriosclerosis disorder                | AS           | MONDO_0002277 |
| Asthma                                   | ASMA         | MONDO_0004979 |
| Attention deficit hyperactivity disorder | ADHD         | EFO_0003888   |
| Bipolar Disorder                         | BD           | MONDO_0004985 |
| Blood coagulation disease                | BCD          | EFO_0009314   |
| Chronic Fatigue Syndrome                 | CFS          | customCFS     |
| Chronic obstructive pulmonary disease    | COPD         | EFO_0000341   |
| Crohn disease                            | CD           | EFO_0000384   |
| Depressive Disorder                      | DD           | MONDO_0002050 |
| Diabetes Mellitus                        | DM           | EFO_0000400   |
| Epilepsy                                 | EPI          | EFO_0000474   |
| Heart failure                            | HF           | EFO_0003144   |
| Hypercholesterolemia                     | HC           | HP_0003124    |
| Hypertension                             | HTN          | EFO_0000537   |
| Lupus erythematosus                      | SLE          | MONDO_0004670 |
| Metabolic syndrome                       | MetS         | EFO_0000195   |
| Multiple Sclerosis                       | MS           | MONDO_0005301 |
| Obesity                                  | OB           | EFO_0001073   |
| Parkinson                                | PD           | MONDO_0005180 |
| Psoriasis                                | PSO          | EFO_0000676   |
| Post-traumatic stress disorder           | PTSD         | EFO_0001358   |
| Rheumatoid arthritis                     | RA           | EFO_0000685   |
| Schizophrenia                            | SCZ          | MONDO_0005090 |
| Sleep Disorder                           | SD           | EFO_0008568   |
| Ulcerative colitis                       | UlCo         | EFO_0000729   |
| Vasculitis                               | VAS          | EFO_0006803   |

<p align="left">
  <em>Table 1. Diseases included in the present study, in alphabetical order. </em>
</p>

For each disease except ME/CFS, the function `Targets4Disease()` queries the Open Targets GraphQL API (v4). Gene–disease associations are collected from multiple evidence sources including:

- genome-wide association studies (GWAS)
- ClinGen curated rare-variant evidence
- gene burden studies from sequencing data

Only genes meeting predefined evidence thresholds are retained. The default filtering parameters include:

| Parameter | Description | Default |
|---|---|---|
| L2G cutoff | minimum locus-to-gene score | 0.5 |
| ClinGen cutoff | minimum ClinGen evidence score | 0.5 |
| GeneBurden cutoff | minimum gene-burden score | 0.5 |
| Sample cutoff | minimum GWAS sample size | 0 |

<p align="left">
  <em>Table 2. Sources of the genes used to build the disease module of ME/CFS. </em>
</p>

For each disease, the pipeline retrieves associated genes using a programmatic query to the Open Targets platform. Each gene is assigned a list label (`GWAS`, `Rare`, or `GWAS/Rare`) and annotated with its STRING preferred name (via the STRING API) and its NCBI Entrez ID (via a local copy of `gene_info.gz`). For myalgic encephalomyelitis/chronic fatigue syndrome (ME/CFS), I built a custom disease module based on the results of the studies in Table 3. 

| Number of cases | Sequencing Method | Gene-Mapping Method    | Genes | Criteria | Reference |
|----------------:|:------------------|:-----------------------|---------|:---------|:----------|
|464              |WGS                | Deep Learning          | 115 | ICC-IOM  |([Zhang S 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12047926/))|
|14767            |Axiom UKB array    | Combinatorial analysis | 259 | CCC-IOM  |([Sardell JM 2025](https://www.medrxiv.org/content/10.64898/2025.12.01.25341362v1))|

<p align="left">
  <em>Table 3. Sources of the genes used to build the disease module of ME/CFS. </em>
</p>

### PPI network construction

The STRING v12.0 human PPI database (`9606.protein.links.v12.0`) is downloaded automatically on first run and filtered to interactions with a combined score ≥ 0.4 (configurable via `STRING.co`). For each disease, `GeneMatrix()` builds a weighted, symmetric adjacency matrix restricted to the disease gene set. The full union of all disease genes is also assembled into a background network used for inter-disease distance calculations. All matrices are stored as `.rds` files under `Modules/`.

### Random disease modules

For each real disease module of size *N*, 1,000 random modules are generated by sampling *N* genes uniformly at random from the pool of all disease genes (`myDiseaseGenes`). These random modules serve as the null distribution for separation (see below). Another 1,000 random diseases with size given by the average size of disease modules are generated using disease genes; they serve as null distribution for the similarity metric based on correlation between Z scores of over-representation analysis (see below). These modules are stored in the folder `Random/` and a zipped copy of it is available ([here](https://huggingface.co/datasets/PaoloMaccallini/CompareME)).

### Module characterisation

For every disease module, the following network properties are computed and compared against the corresponding random null:

| Property | Description |
|---|---|
| Module size | Number of nodes in the largest connected component |
| Mean shortest distance | Average weighted geodesic within the module |
| Mean degree | Average number of PPI edges per gene |
| Mean strength | Average weighted degree (sum of PPI scores per gene) |
| Relative strength | Mean strength / mean degree |

Empirical p-values are derived from the right tail (`P_upper`) or left tail (`P_lower`) of the random distribution. Results are saved to `Modules/Modules_analysis.csv`.

### Over-representation analysis (ORA)

For each disease module, `ORA.fun()` runs hypergeometric over-representation tests against KEGG, Reactome, GO Cellular Component (GO CC), and Disease Ontology gene sets using the `clusterProfiler` and `ReactomePA` packages. Separately, `Tissue.ORA()` computes tissue enrichment z-scores using `TissueEnrich`. Both analyses are also run on each of the 1,000 random modules of the same size to build a pathway-level null distribution.

### Pairwise disease similarity 

All pairwise comparisons are stored under `Comparisons/`.

**Jaccard Index** (gene overlap). It is a standard measure of genetic overlap between two diseases and it is calculated as:

$$J(A,B) = \frac{|A \cap B|}{|A \cup B|}$$

Statistical significance is assessed by a hypergeometric test against the universe of all disease genes. Results in `Comparisons/Jaccard/`.

**ORA correlation** (functional similarity). The Spearman correlation between the z-score vectors of two disease modules across all pathway and tissue terms is computed as a functional similarity score. For both disease modules, a comparison against all 1,000 random modules is performed to build a null correlation distribution of 2,000 correlation coefficients; significance is then the empirical p-value from the merged null of the two diseases being compared. An empirical upper-tail p-value is used to test for significance (custom function `P_upper`). Results in `Comparisons/Correlation/`.

**Network separation SAB** (topological similarity). I used the definition of separation between two gene networks proposed in ([Menche J et al. 2015](https://pmc.ncbi.nlm.nih.gov/articles/PMC4435741/)):

$$S_{AB} = \langle d_{AB} \rangle - \frac{\langle d_{AA} \rangle + \langle d_{BB} \rangle}{2}$$

where $\langle d_{AB} \rangle$ is the mean shortest path between genes of disease A and genes of disease B in the full disease interactome, computed with Dijkstra's algorithm, using the function `distances()` of the package `igraph`. Negative SAB indicates module overlap; positive $S_{AB}$ indicates topological separation. The null distribution is built in two steps: first, we calculate $S_{AB}$ between disease A and each one of the 1,000 random diseases of the same size of disease B; next, we perform the same calculations for disease B. This algorithm generates a distribution of 2,000 random separations. An empirical upper-tail p-value is used to test for significance (custom function `P_upper`). Results in `Comparisons/Separation/`.

### Cross-metric comparison

Pairwise regression (linear, quadratic, and cubic) is performed across all three similarity metrics to quantify their mutual consistency. Plots are saved to `Comparisons/`.

### Gene-level network properties

For each gene in each module, the within-module degree and total STRING interaction count are retrieved. A linear model of within-module degree ~ STRING degree is fitted. This analysis was performed to study the level of connectivity in the complete interactome of those genes that appear isolated in disease modules. Are they isolated because fewer interactions are known for them, overall?

## Results

### General properties of disease modules

Each disease network has been analysed in terms of its main network properties and compared with 1,000 random diseases of the same size, used to build the null distributions. Empirical one-tailed p-values have been used to test for significance, then a Benjamini-Hochberg correction was applied (Table )

| Disease | Vertices | Size | Size_% | p-val | Short_Dist | p-val | Degree | p-val | Strength | p-val | Rel_Strength | p-val |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Alzheimer disease | 134 | 79 | 0.59 | **0.0042** | 1.68 | 0.14 | 5.21 | **0.0014** | 3.29 | **0.0015** | 0.63 | 0.28 |
| Anxiety disorder | 95 | 63 | 0.66 | **0.0019** | 1.75 | 0.54 | 2.42 | **0.0014** | 1.28 | **0.0042** | 0.53 | 1 |
| Arteriosclerosis disorder | 170 | 116 | 0.68 | **0.0031** | 1.60 | **0.0056** | 5.15 | **0.0014** | 3.15 | **0.0015** | 0.61 | 0.57 |
| Asthma | 273 | 213 | 0.78 | **0.0019** | 1.72 | **0.0056** | 9.04 | **0.0014** | 5.59 | **0.0015** | 0.62 | 0.20 |
| Attention deficit hyperactivity disorder | 162 | 109 | 0.67 | **0.0031** | 1.95 | 0.22 | 2.89 | **0.0064** | 1.68 | **0.0089** | 0.58 | 1 |
| Bipolar Disorder | 107 | 50 | 0.47 | 0.056 | 1.65 | 0.31 | 2.00 | **0.0053** | 1.06 | **0.022** | 0.53 | 1 |
| Blood coagulation disease | 95 | 64 | 0.67 | **0.0019** | 1.92 | 0.66 | 9.47 | **0.0014** | 6.96 | **0.0015** | 0.73 | **0.014** |
| **Chronic Fatigue Syndrome** | 369 | 276 | 0.75 | 0.14 | 2.05 | 0.31 | 4.92 | **0.016** | 3.20 | **0.0053** | 0.65 | **0.014** |
| Chronic obstructive pulmonary disease | 120 | 82 | 0.68 | **0.0019** | 1.77 | 0.27 | 3.18 | **0.0014** | 1.89 | **0.0015** | 0.59 | 0.95 |
| Crohn disease | 177 | 141 | 0.80 | **0.0019** | 1.35 | **0.0056** | 10.53 | **0.0014** | 6.43 | **0.0015** | 0.61 | 0.54 |
| Depressive Disorder | 205 | 139 | 0.68 | **0.015** | 2.02 | 0.18 | 2.92 | **0.033** | 1.61 | 0.063 | 0.55 | 1 |
| Diabetes Mellitus | 789 | 712 | 0.90 | **0.0042** | 1.59 | **0.0056** | 11.92 | **0.0014** | 7.02 | **0.0015** | 0.59 | 1 |
| Epilepsy | 78 | 60 | 0.77 | **0.0019** | 1.46 | 0.52 | 7.18 | **0.0014** | 4.39 | **0.0015** | 0.61 | 0.59 |
| Heart failure | 118 | 88 | 0.75 | **0.0019** | 1.72 | 0.24 | 4.34 | **0.0014** | 2.66 | **0.0015** | 0.61 | 0.57 |
| Hypercholesterolemia | 136 | 89 | 0.65 | **0.0019** | 1.63 | 0.066 | 6.97 | **0.0014** | 4.58 | **0.0015** | 0.66 | **0.037** |
| Hypertension | 703 | 597 | 0.85 | 0.43 | 1.88 | 0.52 | 6.96 | 0.49 | 4.14 | 0.48 | 0.59 | 1 |
| Lupus erythematosus | 143 | 100 | 0.70 | **0.0019** | 1.83 | 0.16 | 5.93 | **0.0014** | 3.59 | **0.0015** | 0.61 | 0.57 |
| Metabolic syndrome | 121 | 82 | 0.68 | **0.0019** | 1.42 | 0.066 | 7.09 | **0.0014** | 4.62 | **0.0015** | 0.65 | 0.10 |
| Multiple Sclerosis | 93 | 63 | 0.68 | **0.0019** | 1.45 | 0.31 | 5.05 | **0.0014** | 3.03 | **0.0015** | 0.60 | 0.77 |
| Obesity | 295 | 235 | 0.80 | **0.0019** | 1.93 | 0.066 | 4.90 | **0.0014** | 2.82 | **0.0015** | 0.58 | 1 |
| Parkinson | 67 | 42 | 0.63 | **0.0019** | 1.41 | 0.65 | 5.22 | **0.0014** | 3.45 | **0.0015** | 0.66 | 0.28 |
| Post-traumatic stress disorder | 36 | 8 | 0.22 | 0.13 | 1.19 | 0.88 | 0.67 | 0.084 | 0.40 | 0.068 | 0.61 | 0.64 |
| Psoriasis | 142 | 112 | 0.79 | **0.0019** | 1.38 | **0.020** | 9.08 | **0.0014** | 5.61 | **0.0015** | 0.62 | 0.39 |
| Rheumatoid arthritis | 244 | 178 | 0.73 | **0.0031** | 1.68 | **0.0056** | 8.69 | **0.0014** | 5.27 | **0.0015** | 0.61 | 0.50 |
| Schizophrenia | 255 | 180 | 0.71 | **0.017** | 2.06 | 0.22 | 3.33 | **0.046** | 1.90 | 0.091 | 0.57 | 1 |
| Sleep Disorder | 257 | 176 | 0.68 | 0.068 | 1.91 | 0.056 | 3.35 | **0.038** | 1.88 | 0.093 | 0.56 | 1 |
| Ulcerative colitis | 146 | 89 | 0.61 | **0.0053** | 1.35 | **0.014** | 7.53 | **0.0014** | 4.63 | **0.0015** | 0.61 | 0.57 |
| Vasculitis | 47 | 27 | 0.57 | **0.0019** | 1.15 | 0.70 | 4.17 | **0.0014** | 2.85 | **0.0015** | 0.68 | 0.30 |


The pipeline produces the following outputs:

**Module files** (`Modules/`)

- `<Disease>.csv` — gene list with GWAS/rare variant evidence, STRING name, NCBI ID, PPI count, and within-module degree.
- `<Disease>.rds` — weighted adjacency matrix.
- `<Disease>.tiff/.jpeg` — network graph, nodes coloured by evidence type (GWAS vs. rare variant).
- `<Disease>_ORA.tsv` — pathway enrichment z-scores (KEGG, Reactome, GO, Disease Ontology).
- `<Disease>_Tissue_ORA.tsv` — tissue enrichment z-scores.
- `Modules_analysis.csv` — module-level network properties with empirical p-values for all diseases.
- `myDiseaseGenes.csv` — union catalogue of all disease genes.

**Comparison files** (`Comparisons/`)

- `Jaccard/Score.rds` — Jaccard index and –log₁₀(p) matrices; `Tree_Jaccard.tiff/jpeg` — hierarchical clustering dendrogram; per-disease scatter plots in `Jaccard/TIFF/` and `Jaccard/JPEG/`.
- `Correlation/Score.rds` — ORA Spearman correlation and –log₁₀(p) matrices; `Tree_ORA_cor.tiff/jpeg` — dendrogram; per-disease scatter plots.
- `Separation/Score.rds` — SAB score and –log₁₀(p) matrices; `Tree_SAB.tiff/jpeg` — dendrogram; per-disease scatter plots; null distribution histograms in `Separation/Distributions/`.
- Cross-metric regression plots (`Jaccard_Correlation_*.tiff`, `Jaccard_Separation_*.tiff`, `Correlation_Separation_*.tiff`).

**Random module files** (`Random/`)

- `<Disease>_1000_gene_matrix.rsd` — list of 1,000 random adjacency matrices per disease.
- `<i>_ORA.tsv` / `<i>_Tissue_ORA.tsv` — ORA results for each random module.

## About the pipeline

### Dependencies

The pipeline is written in R and requires the following packages: 

`dplyr`, `rentrez`, `httr`, `jsonlite`, `curl`, `biomaRt`, `stringr`, `DOSE`, `rstatix`, `ReactomePA`, `igraph`, `MASS`, `data.table`, `clusterProfiler`, `org.Hs.eg.db`, `pathview`, `enrichplot`, `GOplot`, `readxl`, `writexl`, `TissueEnrich`, `magick`, `yaml`, `calibrate`, `Matrix`, `parallel`

External databases are downloaded automatically on first run: STRING v12.0 (protein links and protein info for *Homo sapiens*) and NCBI `gene_info`.

### Usage

1. Edit `mydiseases.yml` to select the diseases of interest (uncomment lines to activate).
2. Source `Module_Func.R` and then run `Module_Main.R` sequentially from the repository root.

```r
source("Module_Func.R")
source("Module_Main.R")
```

Results accumulate incrementally: if a module or comparison file already exists it is skipped, so interrupted runs can be safely resumed.

## References

- Menche J. *et al.* Uncovering disease-disease relationships through the incomplete interactome. *Science* 2015.
- Zhang S. *et al.* 2025. PMC12047926.
- Sardell J.M. *et al.* 2025. medRxiv 2025.12.01.25341362.
- Open Targets Platform: https://platform.opentargets.org
- STRING v12.0: https://string-db.org

## License

MIT © Paolo Maccallini 2026
