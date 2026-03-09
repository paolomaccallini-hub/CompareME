# Methods

## Diseases

I selected 28 common diseases, including neurological, psychiatric, metabolic, cardiovascular, inflammatory, and autoimmune conditions. The complete list, with full names, abbreviations, and identifiers (EFO or MONDO codes), is in Table 1. The list of diseases is passed to the script through [mydiseases.yml].
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

## Disease module for ME/CFS

For myalgic encephalomyelitis/chronic fatigue syndrome (ME/CFS), I built a custom disease module based on the results of the studies in Table 2.

| Number of cases | Sequencing Method | Gene-Mapping Method    | Genes | Criteria | Reference |
|----------------:|:------------------|:-----------------------|---------|:---------|:----------|
|464              |WGS                | Deep Learning          | 115 | ICC-IOM  |([Zhang S 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12047926/))|
|14767            |Axiom UKB array    | Combinatorial analysis | 259 | CCC-IOM  |([Sardell JM 2025](https://www.medrxiv.org/content/10.64898/2025.12.01.25341362v1))|

<p align="left">
  <em>Table 2. Sources of the genes used to build the disease module of ME/CFS. </em>
</p>

## Disease modules for the other diseases

