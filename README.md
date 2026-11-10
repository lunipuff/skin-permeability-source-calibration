# Human skin permeability calibration using shared compounds

Analysis code and outputs for **From study variation to calibration: human skin permeability QSPR using shared compounds**.

The primary analysis pools 921 measurements of 349 compounds from 164 publications in the Cheruvu, HuSkinDB, and SkinPiX datasets. We hold out one publication at a time and fit a molecular quantitative structure–property relationship (QSPR) model to the remaining publications. Compounds measured in both the training data and the held-out publication serve as calibration anchors. We estimate a publication-specific offset from those anchors and evaluate it on compounds absent from training. We enumerate every possible panel of one to five anchors, where available, and compare one and five anchors on a fixed subset of publications.

## Reproduce the analysis

Use R with `Rscript` and Cairo PDF support, plus a Python environment with RDKit. The manuscript reports R version 4.5.2 and lme4 version 2.0-1. Install the required packages in an R session.

```r
install.packages(c("readxl", "lme4", "clue", "ggplot2", "patchwork"))
```

Activate the Python environment so the `python` command points to the interpreter with RDKit installed. Check that Python can load RDKit and that R can load Cairo before running the pipeline.

```bash
python -c "from rdkit import Chem, rdBase; print(rdBase.rdkitVersion)"
Rscript -e 'stopifnot(isTRUE(grDevices::grSoftVersion()[["cairo"]] != ""))'
```

**Running the pipeline clears existing contents of `data/processed/`, `R/results/`, `R/tables/`, `R/figures/`, and `R/additional_files/` before rebuilding them.** The runner preserves `R/figures/00_GA.pdf`, the manually prepared graphical abstract. Raw inputs and curation crosswalks are read in place.

```bash
Rscript run_all.R
```

The runner checks inputs and dependencies before cleanup, creates the required directories, and executes each stage in a separate R session. If a stage fails, inspect its log under `R/results/logs/`. The file `R/results/logs/run_all.log` records the run status, R session information, Python executable, and RDKit version.

## Data and curation

| Input | Source |
| --- | --- |
| [`cheruvu.xlsx`](data/raw/cheruvu.xlsx) | Cheruvu et al. (2022), [*Data in Brief*](https://doi.org/10.1016/j.dib.2022.108242) |
| [`huskinDB.csv`](data/raw/huskinDB.csv) | Stepanov et al. (2020), [*Scientific Data*](https://doi.org/10.1038/s41597-020-00764-z) |
| [`skinpix.xlsx`](data/raw/skinpix.xlsx) | Chedik et al. (2024), [*Scientific Data*](https://doi.org/10.1038/s41597-024-03026-4) |
| [`rdkit_descriptors.csv`](data/raw/rdkit_descriptors.csv) | Reference mapping used to propose structures for compounds with missing SMILES strings |

The [publication crosswalk](data/curation/80_source_crosswalk.csv) records publication identifiers and review decisions. The [compound crosswalk](data/curation/80_compound_crosswalk.csv) records compound review decisions and curated SMILES strings. The full runner uses the crosswalks in read-only mode.

Preparation harmonizes permeability to base-10 logarithms with permeability expressed in cm/h. RDKit generates canonical isomeric SMILES and calculates the descriptors used in the model, molecular weight and Crippen logP. Preparation also applies nine permeability corrections verified against the original publications, with a [record-level audit](R/results/80_verified_unit_corrections.csv).

When a publication appears in multiple datasets, the primary analysis retains its records from the highest-priority dataset with eligible records. The priority order is SkinPiX, HuSkinDB, then Cheruvu. Sensitivity analyses vary the database priority or apply agreement thresholds to overlapping publication–compound records.

## Pipeline

| Stage | Purpose |
| --- | --- |
| [80](R/80_prepare_pooled_priority_only.R) | Harmonize units and identities, apply curation, calculate descriptors, and build the pooled dataset. |
| [81](R/81_validate_pooled_calibration.R) | Fit the model by holding out publications, evaluate all available anchor panels, and summarize prediction error. |
| [82](R/82_pooled_source_influence.R) | Check how omitting one publication at a time from the evaluation summaries affects the estimated gains. |
| [83](R/83_run_sensitivity.R) | Rebuild the pool and repeat the evaluation under three overlap-agreement thresholds. |
| [84](R/84_diagnose_pooled_calibration.R) | Compare anchor-based offsets with the best retrospective correction and examine calibration failures. |
| [85](R/85_describe_pooled_variation.R) | Summarize dataset composition and metadata availability, and estimate variance components. |
| [86](R/86_test_compound_sensitivity.R) | Recalculate performance summaries for publications with at least two test compounds. |
| [87](R/87_database_priority_sensitivity.R) | Evaluate all six database priority orders, reusing results when orders produce identical datasets. |
| [99a](R/99a_make_tables.R) | Generate CSV tables and additional files for the manuscript. |
| [99b](R/99b_make_figures.R) | Generate Figures 1–4 as PDFs. |

Bootstrap settings and seeds are recorded with the analysis outputs. The workflow uses 2,000 bootstrap replicates for calibration and variance summaries and 20,000 for the publication omission analyses.

## Outputs

| Location | Contents |
| --- | --- |
| [`data/processed/80_pooled.csv`](data/processed/80_pooled.csv) | Primary analysis dataset with record, compound, and publication identifiers. |
| [`R/results/`](R/results/) | Uncalibrated test predictions, anchor-panel results, fit diagnostics, bootstrap summaries, and sensitivity analyses. Start with the [one-anchor summary](R/results/81_pooled_calibration_summary_all.csv) or the [paired one-versus-five comparison](R/results/81_pooled_one_vs_five_anchors.csv). |
| [`R/tables/`](R/tables/) | CSV inputs for the manuscript tables and figures. |
| [`R/additional_files/`](R/additional_files/) | Numbered CSV files for manuscript submission. |
| [`R/figures/`](R/figures/) | The graphical abstract and four analysis figures. |

### Additional files for submission

| Manuscript reference | CSV |
| --- | --- |
| Additional file 1 | [Publication identities and calibration eligibility](R/additional_files/Additional_file_1_publication_index.csv) |
| Additional file 2 | [Calibration performance and panel summaries by publication and anchor count](R/additional_files/Additional_file_2_calibration_by_publication.csv) |
| Additional file 3 | [Decomposition of calibration error](R/additional_files/Additional_file_3_failure_decomposition.csv) |

Stage 99a exports these files from the corresponding CSVs in `R/tables/`. To refresh the tables and additional files from existing analysis results, run this command from the repository root.

```bash
Rscript R/99a_make_tables.R
```

After a successful run, you can regenerate an output without clearing the full set of results. For example, run `Rscript R/99b_make_figures.R` from the repository root to regenerate the figures from the existing analysis outputs.

## Citation and license

When using this work, cite the associated manuscript and the original dataset publications listed above.

The original analysis code is distributed under the [MIT License](LICENSE). Third-party data and software are subject to their respective terms.

ChatGPT assisted with language editing, code refinement, and debugging. The authors reviewed the resulting text and code.
