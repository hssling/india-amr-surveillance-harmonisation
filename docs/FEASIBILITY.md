# Feasibility check, 28 November 2026

**Verdict: GO.** The study can be done with the data on hand, and the feasibility work already shows its central
finding: India's human (ICMR-AMRSN) and animal (ICAR-FAO INFAAR) surveillance networks test largely different drug
panels.

## INFAAR analytical report 2019-2022 (FAO-ICAR, 2024; doi:10.4060/cc9303en)

- 32 pages. The numbered tables describe the surveillance (organism-antimicrobial combinations, isolates by sector,
  samples); resistance results are in bar charts without data labels, with a few exact values in the text.
- Charts are vector graphics: page 16 (Figure 11, E. coli in cattle, buffalo, goats, sheep, poultry and pigs,
  12 antibiotics) holds exactly 72 filled bars (6 x 12), so bar heights can be extracted from the geometry.
- A prototype extraction reproduced 6 of 7 values stated in the text within 0.6 percentage points (goat cefotaxime
  41.6 vs 41%; goat chloramphenicol 5.7 vs 6%; goat ampicillin 26.2 vs 26%; sheep cefotaxime 37.0 vs ~37%; sheep
  chloramphenicol 5.1 vs ~5%; sheep imipenem 10.6 vs ~10%). The seventh (goat amikacin 21.3 vs 35%) and missing labels
  for the first bars of the right-hand panels show a bar-to-label matching error in the prototype, to be fixed.
- Some fisheries figures are raster images and need digitising.
- Method for the study: vector extraction where possible; every extracted value checked against every value stated in
  the text; an independent second extraction (manual digitising) on a random subset; agreement reported.

## Overlap with ICMR-AMRSN (extract in D:\onehealth amr meta nma\data_processed\amrsn_long.csv)

- AMRSN E. coli 2019-2022 drugs: amikacin, cefepime, cefotaxime, ceftazidime, ceftriaxone, ciprofloxacin, colistin,
  ertapenem, gentamicin, imipenem, levofloxacin, meropenem, piperacillin-tazobactam, tigecycline.
- INFAAR E. coli panel: AMC, AMK, AMP, CAX (code to confirm), CHL, COT, CTR, CTX, ENR, IPM, NAL, TET.
- Shared drugs: 4 of 12 (amikacin, cefotaxime, ceftriaxone, imipenem). Fluoroquinolones comparable only as a class
  (enrofloxacin in animals, ciprofloxacin in humans). S. aureus/MRSA is reported by both.

## Next

Write the protocol: harmonisation audit (organisms, drug panels, interpretive standards, sampling frames, reporting
units, years) and concordance for the shared organism-drug pairs, 2019-2022.
