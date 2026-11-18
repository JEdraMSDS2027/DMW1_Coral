# Coral bleaching: final report, preprocessing and PCA

The final study entry point is **`consolidated_notebook.ipynb`**. It executes the complete raw-data-to-PCA workflow and presents the rubric-aligned report: abstract, introduction, dataset description, methodology, results/discussion, conclusion, limitations/recommendations, APA references and supporting appendices. It does not require either supporting notebook to have been run first.

The detailed supporting notebooks remain **`coral_cleaning_preprocessing.ipynb`** for preparation and **`coral_pca.ipynb`** for PCA. Their existing content is preserved. PCA methodological corrections are documented in **`PCA_CHANGELOG.md`**.

An alternative reading version is **`consolidated_report_integrated.ipynb`**. It has no appendices: 18 compact tables and 10 selected figures are placed within the relevant report sections. Repeated raw-field listings, duplicate plots and visible legacy-test/fold logs are omitted while the calculations and numerical audits remain reproducible. Run `python scripts/execute_notebook.py consolidated_report_integrated.ipynb`; its HTML, tables, figures and model are saved separately under `outputs/integrated_report/`. The original final notebook and its HTML reading copy remain unchanged.

The combined export preserves **23,069 sampling events and 62 columns**, including **35 processed environmental inputs** for PCA. Categorical labels, comments, and environmental fields excluded from PCA remain available as metadata. The supporting preprocessing notebook stops before PCA; the final study notebook executes both stages and retains all meaningful supporting analyses in appendices.

## Project structure

```text
consolidated_notebook.ipynb         Final end-to-end technical report
consolidated_report_integrated.ipynb  Alternative integrated report without appendices
coral_cleaning_preprocessing.ipynb   Main report, with saved tables and figures
coral_pca.ipynb                     Main consolidated PCA report
coral_pipeline/
    data.py                        Schema, input loading, dictionary, roles, exports
    cleaning.py                    Types, validity, sampling events, cohort, feature selection
    preprocessing.py               Fitted medians, power transforms, scaling, reuse
    reporting.py                   Quality tables, plots, readiness checks
    pca.py                         Extraction, retention, sign orientation, source/stability refits
    pca_statistics.py              Existing effects/tests, temporal/CV checks, overlap correction
    pca_plots.py                   PCA figures returned without display/export side effects
    pca_reporting.py               Summary tables, numerical comparison, narrative, exports
    final_report.py                Complete raw profiles, stage percentages, captions, report prose
scripts/execute_notebook.py         Execute in a fresh kernel and export HTML
tests/test_cleaning.py              Tests of consequential cleaning rules
requirements.txt                   Versions used for execution
source_manifest.json               Fingerprints of original project inputs
PCA_CHANGELOG.md                    Structural changes and methodological corrections
outputs/pca/                       Analytical tables, 16 PCA + 4 preparation/workflow figures, model, metadata
outputs/final_report/              Final HTML report, computed abstract and conclusion
outputs/integrated_report/         Separate integrated reading copy and numerical exports
processed/
    coral_pca_features.csv          Numeric-only PCA inputs, 23,069 × 35
    coral_pca_ready.csv             All sample-level columns, 23,069 × 62
    coral_cleaned.csv               Same 62 columns; selected inputs in original units after imputation
    coral_sample_metadata.csv       All 27 non-PCA fields, in the same row order
    coral_substrate_records.csv     Original substrate rows for retained events, linked by Sample_ID
    coral_preprocessing_metadata.json
    coral_transformers.joblib
    coral_cleaning_preprocessing.html
    figures/                       Three preprocessing figures
    ...                            Quality and preprocessing audit CSVs
```

The original files remain in the project root:

| Source | Use in consolidation |
|---|---|
| `coral.csv` | Read-only raw records; matches the original notebook's SHA-256 |
| `data_dictionary.md` | Variable meanings, units, and analytical-role notes |
| `coral_eda_cleaning_preprocessing.ipynb` | Original cleaning decisions, code, and saved quality findings |
| `coral_data_cleaning_summary.md` | Reconciliation of the original decisions and expected cohort |
| `coral_pca_ready.csv` | Previous standardised export, retained for comparison |
| `coral_pca_analysis.ipynb` | Historical downstream notebook, inspected for input expectations only |

The two older notebooks are historical references. `coral_pca.ipynb` consolidates the newer methodology notebook and the distinct analyses in `coral_pca_analysis.ipynb`. The older notebooks themselves remain unchanged. No Python modules were present in the original project folder.

## Reproduce the report

Use Python 3.13 and run from the project directory:

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\.venv\Scripts\python.exe scripts/execute_notebook.py consolidated_notebook.ipynb
.\.venv\Scripts\python.exe scripts/execute_notebook.py
.\.venv\Scripts\python.exe scripts/execute_notebook.py coral_pca.ipynb
.\.venv\Scripts\python.exe -m pytest -q
```

Alternatively, open the report in Jupyter or VS Code with that environment and run all cells from the project directory. The script uses its own interpreter for the kernel and sets the working directory explicitly. Temporary Jupyter state stays in the ignored `.jupyter/` directory. A headless run needs local loopback ports for Jupyter messaging. Derived files in `processed/` are regenerated; the raw dataset and previous root-level export are never overwritten.

The first execution command above produces the complete final report. The next two commands are optional executions of the preserved supporting notebooks. The final-report script refreshes the opening abstract from the completed numerical results. For an interactive Jupyter run, a current computed abstract is also produced at the end and written to `outputs/final_report/abstract.md`.

The source is checksum-verified BCO-DMO Version 2, DOI **10.26008/1912/bco-dmo.773466.2**. Download its `global_bleaching_environmental.csv` as `coral.csv` to reconstruct the study locally; the expected MD5 is `42a6f10a6efa9e6b926d04b86b45a19f`. The final report cites original collection methods and explicitly discloses limited documentation for two small source labels and uncertain unit labels. It preserves stored measurements rather than applying speculative conversions.

The final report uses numbered captions and a specific workflow diagram. Main sections show core evidence; Appendices A–F retain full profiles, transformation audits, legacy statistical/CV comparisons, all supporting figures, sensitivity/time analyses and reproducibility checks. The completed project has **33 passing tests**. The final notebook contains no study-author name or personal absolute paths. Follow the evaluation rubric's instruction to omit dataset files from the assessed bundle while providing the notebook, modules, dictionary, requirements and data-access instructions.

## Cleaning method

1. Replace `nd` markers and explicitly parse documented numeric fields; reject unexpected tokens and invalid ranges.
2. Verify within-event consistency, remove exact duplicate records if any, collapse substrate records to one row per `Sample_ID`, and average differing bleaching outcomes within the same event. Do not merge distinct event IDs with identical measurements. Preserve substrate-specific labels and cover values in `coral_substrate_records.csv`, because they do not have one unambiguous value per event.
3. Preserve the original labelled cohort: exclude events without observed bleaching and retain fixed groups, Low 0–10%, Moderate >10–50%, High >50–100%. Never impute the outcome.
4. Validate dates and retain only sampling years **1998 onward**, excluding 1997 and earlier before feature selection or fitted preprocessing. Keep identifiers, dates, coordinates, outcomes, and missingness flags outside the feature matrix. Retain nominal labels and free text in the combined CSV as metadata; exclude them only from PCA inputs. No categorical encoding is introduced.
5. Exclude the constant `SSTA_Mean` and near-duplicate `TSA_Standard_Deviation` (absolute Spearman rho ≥ 0.999) from PCA inputs, preserving their original values and missingness in the combined CSV. Retain the other correlated environmental families as inputs.
6. Exclude events missing more than half of the selected inputs. Fit column medians on the remaining events and preserve `Depth_m_was_missing` as metadata.
7. Apply Yeo–Johnson where absolute sample skewness exceeds 0.75, with recorded numerical rescaling where needed, then standardise all inputs to mean zero and population SD one. Preserve plausible extremes without clipping or deletion.

The cohort reconciles as **41,361 substrate records → 27,005 events → 23,271 labelled events → 23,126 events from 1998 onward → 23,069 retained events**. There are eight events with differing observed bleaching values, 145 labelled events excluded by year, 57 sparse events excluded, and 1,608 depths imputed. Thirty-one features receive a power transform; four are standardised directly. The final recorded date range is 1998–2019. The excluded labelled event IDs and dates are saved in `processed/year_exclusions.csv`.

## Using the derived data

```python
import json
from pathlib import Path
import pandas as pd

output = Path("processed")
metadata = json.loads((output / "coral_preprocessing_metadata.json").read_text())
data = pd.read_csv(output / "coral_pca_ready.csv")
X = data[metadata["pca_input_columns"]]
```

`coral_pca_features.csv` contains the same 35-column `X` without metadata. It has no index column. Its row order matches `coral_pca_ready.csv`, `coral_cleaned.csv`, and `coral_sample_metadata.csv`; use `Sample_ID` from the combined or metadata table to track events. The 35 selected inputs in the cleaned table are median-imputed, and those inputs in the ready table are already transformed and scaled. Do not transform or scale them again by accident.

The combined table retains **60 of the 62 source columns**, plus `Bleaching_Group` and `Depth_m_was_missing`. The only source columns absent from this one-row-per-event table are `Substrate_Name` and `Percent_Cover`; their original one-to-many records are preserved separately in `coral_substrate_records.csv` for the same retained events. This detail table also includes the recorded substrate-row bleaching percentage and both event/site IDs. No arbitrary first-substrate selection, cover averaging, or flattening is introduced.

Only names listed in `pca_input_columns` are PCA inputs. Do not pass the full combined CSV, or every numeric column, into PCA: it also contains identifiers, outcomes, date parts, coordinates, flags, a constant, and a redundant measurement. Categories and comments remain readable text. The two excluded environmental fields remain in original units, and non-PCA missing values are left unchanged. The metadata JSON distinguishes columns physically absent from the combined table (`removed_columns`) from columns excluded only from PCA (`excluded_from_pca_columns`).

`coral_transformers.joblib` stores a dictionary with `medians`, `transformers`, and `pca_input_columns`. The reusable `coral_pipeline.preprocessing.transform_features()` applies this state to cleaned sample-level original-unit measurements without refitting. Apply the same quality, sampling-year, and sparse-profile rules first. Use the recorded package versions when reloading fitted state. JSON metadata contains ordered feature names, fitted transformation details, exclusions, population scope, software versions, and output hashes; audit CSVs retain full numerical detail.

## Assumptions and changes

The underlying cleaning and transformation rules are preserved, with the requested sampling-year cutoff of 1998 inclusive. The combined table now retains every source column that can be represented at the event level. Exclusion from the PCA feature list no longer causes a column to disappear from the combined export. The implementation includes a strict schema, failures on conflicting retained sample fields or invalid dates, an explicit numeric-only matrix, an outcome-aggregation audit, saved provenance, and validation that does not fit PCA. Near-perfect rank correlation is described as near redundancy rather than exact numerical equality.

The new combined export retains the 1998-onward subset of event IDs from the supplied historical export. Interpretation values for those events are unchanged. Feature selection and fitted medians, power transforms, and scaling are recalculated on the restricted cohort; standardised values therefore change. Under the unchanged skewness threshold, `Cyclone_Frequency` now receives standardisation only.

The report corrects two overstatements in the original explanation: median imputation can affect covariance, and selective transformation does not establish normality or guarantee freedom from influential observations. The raw data pass all basic measurement checks; the final table passes all readiness checks. Tests cover record weighting, within-event conflicts, missingness and group boundaries, date disagreement, inclusive sampling-year boundaries, training-only transform reuse, metadata preservation, and substrate-detail linkage.

Validation completed: all 11 notebook code cells executed sequentially without errors or warning outputs, 18 tests passed, and all 20 readiness checks passed. The three figures were visually inspected. Saved CSV row alignment, reuse of the serialised transformations, output fingerprints, unchanged source fingerprints, preservation of all sample-level columns, and unchanged PCA inputs were checked independently.

The global fit is appropriate for the intended descriptive PCA of this cohort. For later predictive evaluation, split sample-level data by site, and by time when needed, before fitting feature selection, imputation, power transforms, scaling, or PCA. Fit each step inside the training folds. The exported globally fitted table must not be used to claim unbiased predictive performance. Historical thermal summaries also require temporal-availability review before forecasting.

Carry forward the labelled-cohort and sampling-period selection, repeated-site dependence, source-related depth missingness, imputation effects, approximate day precision, residual zero-heavy skew, and unequal weighting of correlated thermal families. Source and exposure labels in the sidecar support later sensitivity checks. Windspeed units remain as supplied pending verification; no unsupported conversion is applied.

## PCA report and outputs

The PCA notebook uses only the ordered 35-column preprocessing metadata list. It fails explicitly on missing features, invalid years/groups, non-finite or constant inputs, or nonstandardised measurements; it does not infer numeric features or silently repeat preprocessing. The loaded 62-column table retains all categorical and text fields for interpretation. The original 1998–2019 cohort remains unchanged.

The report connects component retention and interpretation to effect sizes, original Welch/FDR and Hotelling comparisons, source interpretation, within-source contrasts, source-specific fits, variable-level effects, grouped prediction, robustness checks, imputation sensitivity, and time/hemisphere interpretation. The old equal-weight-SD variable effect and newer sample-size-weighted Cohen's d are labelled separately.

The output folder contains full-precision tables, named PNG figures, a sign-oriented `pca_model.joblib`, `run_metadata.json`, generated `findings.md`, and an HTML reading copy. Exported PCA scores retain the full interpretation table alongside the retained component scores. Load the model with the recorded environment and use `model.transform(X.to_numpy())`; it returns all 35 component scores, of which the manifest's selected count is retained for the report.

The parallel rule and original `legacy_min` option both retain seven components here (75.9% variance). Independent-site Welch results remain labelled legacy because 196 sites occur in both groups. The additional site-cluster intervals preserve that overlap. Legacy globally fitted CV and fold-fitted CV are shown separately. See the change log for estimands, assumptions and numerical effects. The corrected CV remains conditional on the fixed study cohort, feature names and chosen dimension.

The consolidated PCA notebook has been executed in a fresh kernel; the expanded project suite has 33 passing tests. `tests/pca_reference.json` records original numerical reference values for this input fingerprint, and `reference_validation.csv` reports the comparisons. If the input changes, the notebook marks the reference as incomparable rather than claiming agreement.
