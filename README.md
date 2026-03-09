# Methods

## Diseases

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

## Disease modules

### ME/CFS

For myalgic encephalomyelitis/chronic fatigue syndrome (ME/CFS), I built a custom disease module based on the results of the studies in Table 2. 

| Number of cases | Sequencing Method | Gene-Mapping Method    | Genes | Criteria | Reference |
|----------------:|:------------------|:-----------------------|---------|:---------|:----------|
|464              |WGS                | Deep Learning          | 115 | ICC-IOM  |([Zhang S 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12047926/))|
|14767            |Axiom UKB array    | Combinatorial analysis | 259 | CCC-IOM  |([Sardell JM 2025](https://www.medrxiv.org/content/10.64898/2025.12.01.25341362v1))|

<p align="left">
  <em>Table 2. Sources of the genes used to build the disease module of ME/CFS. </em>
</p>

### Other diseases

For all diseases except ME/CFS, the function `Targets4Disease()` queries the **Open Targets GraphQL API** (v4). Only protein-coding genes are retained; evidence is accepted from three sources:

- **GWAS credible sets**: genes are included when the Locus-to-Gene (L2G) score ≥ 0.5 and the GWAS sample size passes a configurable threshold (default: any size).
- **ClinGen**: curated gene–disease relationships with score ≥ 0.5.
- **Gene burden tests**: rare-variant collapsing analyses with score ≥ 0.5.

Each gene is assigned a list label (`GWAS`, `Rare`, or `GWAS/Rare`) and annotated with its STRING preferred name (via the STRING API) and its NCBI Entrez ID (via a local copy of `gene_info.gz`). This strategy was not adopted for ME/CFS because summary statistics from DecodeME (the biggest GWAS on ME/CFS) are not yet available in Open Targets.

### 3. PPI network construction

The **STRING v12.0** human PPI database (`9606.protein.links.v12.0`) is downloaded automatically on first run and filtered to interactions with a combined score ≥ 0.4 (configurable via `STRING.co`). For each disease, `GeneMatrix()` builds a weighted, symmetric adjacency matrix restricted to the disease gene set. The full union of all disease genes is also assembled into a background network used for inter-disease distance calculations. All matrices are stored as `.rds` files under `Modules/`.
