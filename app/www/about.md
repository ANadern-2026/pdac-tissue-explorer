<div class="about">

## What this is

The **PDAC Tissue Explorer** shows where a gene is expressed in human pancreatic ductal adenocarcinoma (PDAC) and how its expression differs between tumour and non-tumour pancreas, using three public datasets. It is meant for researchers and clinicians who want a quick, sample-level view without writing code. All figures and the underlying numbers can be downloaded.

## Data

| Source | What it contributes | Samples |
|---|---|---|
| **CPTAC PDAC** (Cao et al., *Cell* 2021), via LinkedOmics | Protein (TMT mass spectrometry, gene level) and RNA in tumours and matched adjacent normal tissue | 135 PDAC tumours; 64 protein and 19 RNA tumour/normal pairs |
| **TCGA-PAAD** and **GTEx** pancreas, uniformly processed by **recount3** (Wilks et al., *Genome Biol* 2021) | Bulk RNA of PDAC vs healthy pancreas | 147 PDAC (ductal histology), 4 adjacent normal, 330 healthy pancreas |
| **PDAC single-cell atlas** (Loveless et al., *Clin Cancer Res* 2025; 12 integrated studies) | Cell-type resolution: which cells express a gene, in which tissue state | ~620,000 cells from 234 samples: healthy donor, adjacent normal, primary tumour (treatment-naive and treated), metastases |

## Methods in brief

**CPTAC.** Protein values are log2 median-normalised intensities, RNA values are log2 RSEM upper-quartile, as provided by LinkedOmics. Adenosquamous carcinomas are excluded. Tumour and adjacent normal tissue are compared within the same patient (paired Wilcoxon signed-rank test; exact for small samples); the median change has a 95% bootstrap confidence interval. Proteins of low abundance are often quantified in only part of the patients; the app reports this.

**TCGA/GTEx.** recount3 gene counts were converted to TPM and log2(TPM + 1). Only TCGA-PAAD primary tumours of ductal histology are used. Genes expressed at TPM ≥ 1 in at least 10% of samples of any cohort are included. TCGA and GTEx are different studies; no cross-study batch correction is applied, so differences are descriptive. Composition scores (fibroblast, endothelial, epithelial, acinar, immune) are mean log2 TPM of marker genes.

**Single-cell atlas.** The published atlas object was filtered before summarisation: the cycling-cell cluster was removed and cells with high ambient-RNA contamination (decontX score above 0.30, or above the 98th percentile within each sample and cell type) were excluded; cells with fewer than 200 detected genes were dropped. Expression values are the atlas's log-normalised values (log1p of counts per 10,000). For every sample and cell type the app stores the mean normalised expression and the fraction of cells with any detected transcript. *Average expression* = log1p(mean CP10k). A sample contributes to a cell type only with at least 20 cells of that type. Cell types are the atlas's annotations; ductal cells are split into those assigned to the atlas's tumour subclusters and the remaining (mostly non-malignant) ductal cells. Fibroblast subtypes (myFb 1–5, iFb 1–3) and tumour-cell groups (classical-like, hybrid, basal-like) follow the atlas publication. "Share of transcripts" is the fraction of all transcripts of a gene, pooled over all cells of a context, that come from each cell type. Tissue states: healthy donor; adjacent normal; treatment-naive primary tumour; treatment-naive liver metastasis; primary tumour after chemotherapy or chemoradiation.

## Limitations

- Abundance of an RNA or protein does not show its activity, modification or function.
- Adjacent normal tissue from pancreatic-cancer surgery is inflamed and fibrotic; it is not healthy pancreas. Healthy pancreas (GTEx, donors) is mostly acinar tissue, so bulk tumour-vs-healthy differences partly reflect tissue composition.
- Single-cell RNA sequencing misses many low-abundance transcripts; "not detected" is not proof of absence.
- In the atlas, tissue state and treatment are partly confounded with the contributing study. Metastases are almost all liver metastases, and few metastatic samples contain enough fibroblasts.
- p-values are unadjusted and refer to the single comparison stated; with many genes, some will be small by chance.

## How to cite

If you use this resource, please cite it — Nadernezhad A., *PDAC Tissue Explorer*, Zenodo, [doi:10.5281/zenodo.23264879](https://doi.org/10.5281/zenodo.23264879); processed data [doi:10.5281/zenodo.23264998](https://doi.org/10.5281/zenodo.23264998) — **and the original data sources**: Cao et al., *Cell* 2021; Wilks et al., *Genome Biol* 2021 (recount3); The Cancer Genome Atlas Research Network, *Cancer Cell* 2017; the GTEx Consortium, *Science* 2020; Loveless et al., *Clin Cancer Res* 2025, and the twelve single-cell studies it integrates.

## Author

Ali Nadernezhad, Leibniz-Institut für Polymerforschung Dresden e.V. ([ORCID 0000-0002-3761-1601](https://orcid.org/0000-0002-3761-1601)). Source code: <https://github.com/ANadern-2026/pdac-tissue-explorer>. Live app: <https://anadern-2026-pdac-tissue-explorer.share.connect.posit.cloud/>.

## Development

Data-preparation code, app code and provenance logs are public in the source repository. Code was drafted with the assistance of an AI model (Anthropic Claude) and executed, checked and curated by the author. Code is released under the MIT license; the derived data tables under CC BY 4.0, in line with the licences of the source data.

</div>
