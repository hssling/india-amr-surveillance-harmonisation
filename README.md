# Do India's human and animal antimicrobial resistance surveillance networks speak the same language?

Code for a harmonisation audit and concordance analysis of the ICMR Antimicrobial Resistance Surveillance and Research
Network (AMRSN) and the ICAR-FAO Indian Network for Fishery and Animal Antimicrobial Resistance (INFAAR), 2019-2022.
Study by Siddalingaiah H S.

**Status:** v1.0.0, the code release accompanying submission of the manuscript (November 2026). The analysis plan was
frozen and archived on 28 November 2026 (v0.1.0) before any human-animal comparison was computed; every departure from
it, and every analysis added later, is listed in `docs/DEVIATIONS.md`. Extracted datasets and analysis outputs will be
added to this archive on acceptance of the manuscript.

**Archive:** Zenodo concept DOI [10.5281/zenodo.24601806](https://doi.org/10.5281/zenodo.24601806) (always resolves to
the latest version; analysis plan v0.1.0: [10.5281/zenodo.24601808](https://doi.org/10.5281/zenodo.24601808)).

## Source documents (public; not redistributed here)

Place these in `data_raw/` (file names as below). SHA-256 prefixes identify the exact versions analysed.

| File | Document | SHA-256 (first 16) |
|---|---|---|
| `amrsn_2019.pdf` | ICMR AMRSN annual report 2019 (icmr.gov.in, Final_AMRSN_Annual_Report_2019_29072020.pdf) | 58d2c5fa3cde7a97 |
| `amrsn_2020.pdf` | ICMR AMRSN annual report 2020 (icmr.gov.in, AMRSN_annual_report_2020_1.pdf) | 345a560cf8bcc361 |
| `amrsn_2021.pdf` | ICMR AMRSN annual report 2021 (icmr.gov.in, AMR_Annual_Report_2021.pdf) | 267eb96a8814b19f |
| `amrsn_2022.pdf` | ICMR AMRSN annual report 2022 (icmr.gov.in, AMRSN_Annual_Report_2022.pdf) | cf669d2a8193cdd2 |
| `infaar_2019_22.pdf` | FAO and ICAR (2024), INFAAR analytical report 2019-2022, doi:10.4060/cc9303en | a81afc8f073bbb7e |

The 2019 ICMR standard operating procedures (human and veterinary bacteriology), the EU Decision 2020/1729 and the WHO
(2017) guidance used in the audit are cited in the manuscript; they are not inputs to the code.

## Pipeline

Run from the repository root with Python 3 (PyMuPDF, pandas, numpy, scipy, networkx, matplotlib):

| Step | Script | Does |
|---|---|---|
| 1 | `src/01_extract_infaar.py` | Reads INFAAR bar heights from the vector drawing instructions of the report |
| 2 | `src/02_label_infaar.py` | Assigns bars to figure, source and antimicrobial |
| 3 | `src/03_validate_infaar.py` | Pre-specified validation against every percentage stated in the report text |
| 4 | `src/04_pixel_extract_infaar.py` | Second reading of every bar from the rendered page |
| 5 | `src/05_extract_amrsn.py` | Reads AMRSN tables by cell position; checks every count against its printed percentage |
| 6 | `src/06_amrsn_primary.py` | Primary human series; checks against the reports' text and later revisions |
| 7 | `src/07_concordance.py` | Human-animal comparisons for shared organism-antimicrobial pairs |
| 8 | `src/08_audit_panels.py` | Panel overlap (Jaccard) and antimicrobial-class coverage |
| 9 | `src/09_min_panel.py` | Minimum common panel; EU 2020/1729 and WHO 2017 benchmarks |
| 10 | `src/10_network_analyses.py` | Measurement DAG, panel-network modularity, resistance-profile network, Newcombe intervals |
| 11 | `src/11_figures.py` | Figures |

Curated inputs: `docs/infaar_validation_set.csv` and `docs/infaar_pooled_statements.csv` (text-stated INFAAR values,
quoted), `docs/amrsn_text_checks.csv` (text-stated AMRSN values), `docs/amrsn_urine_image_transcription.csv` (two image
tables transcribed with their centre rows).

## Licence

Code: MIT (`LICENSE`). Derived aggregate outputs, when released: CC BY 4.0 (`LICENSE-DATA`).

## Citation

See `CITATION.cff`.
