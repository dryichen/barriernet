[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)

# barriernet

barriernet compresses spatial tumor immune microenvironment (TIME) architecture into two portable, marker-level axes — **Barrier** (CAF/ECM + suppressive TAM) and **NetImmune** (APM + IFN + TLS + T-cell minus Barrier) — enabling prognostic stratification, ICI response prediction, and drug prioritization across bulk RNA, proteomics, IHC/IF, and spatial transcriptomics platforms.

## Install

You can run without installing:

```bash
python -m barriernet --help
```

Or install editable:

```bash
pip install -r requirements.txt
pip install -e .
```

## Quick start

### 1) Bulk RNA / proteomics (samples×genes or genes×samples)
```bash
barriernet bulk --expr expr.tsv --out out_bulk
```

Outputs:
- `barriernet_bulk_scores.tsv` (modules + Barrier/Immune/NetImmune)
- `barriernet_bulk_used_genes.tsv`
- `barriernet_bulk_panel_matrix.tsv` (only genes actually used)

### 2) IHC / IF panel table (ROI/cell-level allowed; aggregated to sample)
```bash
barriernet ihc --table ihc_table.tsv --sample-col sample_id --out out_ihc
```

Outputs:
- `barriernet_ihc_scores.tsv`
- `barriernet_ihc_used_markers.tsv`

### 3) Spatial Visium (uses `obs['array_row']`, `obs['array_col']`, `obsm['spatial']`)
```bash
barriernet spatial --h5ad your.with_state_and_ring.h5ad --out out_spatial
```

Outputs per input:
- `<stem>_spot_rings.tsv`
- `<stem>_ring_enrichment.tsv`
- `<stem>.with_barriernet.h5ad` (ring metrics written to obs)

## Geneset override

Default panels live in:
- RNA/spatial: `barriernet/data/geneset_rna.yaml`
- IHC/protein: `barriernet/data/geneset_ihc.yaml`

Override with your own:
```bash
barriernet bulk --expr expr.tsv --geneset my_geneset.yaml --out out_bulk
```

## Data availability

| Dataset | Accession | Platform |
|---------|-----------|----------|
| scRNA-seq (29 patients) | GSE183904 | 10x Chromium |
| Spatial transcriptomics (10 sections) | GSE251950 | 10x Visium |
| TCGA-STAD (n=395) | GDC Portal | RNA-seq |
| GSE84437 (n=433) | GSE84437 | Microarray |
| GSE62254/ACRG (n=300) | GSE62254 | Microarray |
| PDC-STAD proteomics (n=130) | PDC | Mass spec |
| ICI cohort (n=44) | PRJEB25780 | RNA-seq |
| GDSC2 (28 GC lines) | cancerrxgene.org | IC50 + expression |

## License

This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.
