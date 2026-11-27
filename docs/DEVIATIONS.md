# Deviations from and additions to the analysis plan (v1.0, frozen 28 November 2026)

Each entry gives what changed, when, and why. Results are not reported here; they are released with the submission.

## A1. Investigation of validation disagreements (added 28 November 2026, after the validation was run)

The pre-specified criterion (section 4.2: at least 90% of text-stated percentages within 2 percentage points of the
extracted bar) is applied unchanged, with each value compared with the panel named in the report's figure captions.
To investigate the disagreements, as section 4.2 requires, two checks were added:

1. **Exchanged panels.** For Figure 10, the comparison is repeated with panels C and D exchanged, to test whether the
   text and the captions describe the panels in opposite order.
2. **Pooled statements.** Where the text gives a percentage for all sources combined and the figure legend gives the
   number of isolates from each source (Figures 2-6), the isolate-weighted mean of the extracted per-source values is
   compared with the pooled percentage (docs/infaar_pooled_statements.csv). Where a per-source text value disagrees
   with its bar, the per-source value implied by the pooled statement is back-calculated. The weighting assumes every
   isolate was tested against every antimicrobial, so this check is approximate.

Neither check changes the acceptance decision, which rests on the pre-specified comparison alone.

## A2. Labelling of bars (28 November 2026)

Series in grouped charts are assigned by fill colour, ordered by position within the category slot, after assignment
by position alone mislabelled bars whose category label is not centred on its group. A series with no bar in a
category is recorded as 0% and flagged as inferred, because bars of zero height are not drawn.

## A3. Second extraction (28 November 2026)

The pixel-based extraction reads bar heights from the rendered page but shares the value-scale calibration and the bar
positions and colours with the vector extraction; it is therefore a check on the reading of bar heights, not a fully
independent extraction. Two pixel columns either side of each bar's centre are read, because most charts draw an error bar
(undefined in the report) on each bar's centre line.

## C. Clarifications made during AMRSN extraction (28 November 2026), before any human-animal comparison was computed

C1. **AMRSN specimen category for E. coli.** No AMRSN table combines all specimens for E. coli. The rule in section
4.1 (the specimen category with the most isolates) is applied among categories for which susceptibility is reported:
"all specimens except urine and faeces" in every year (2019: the report has no urine susceptibility table for
Enterobacterales; 2020-2022: non-urine isolates outnumber urine isolates in each report's specimen distribution
table). Urine isolates are analysed as a secondary series.

C2. **Report version.** Each year's value is taken from that year's own annual report. Later reports' trend tables
revise earlier years; the revisions are extracted and reported as a check.

C3. **Measure.** AMRSN reports percentage susceptible only. The human measure is non-susceptible (100 minus %S,
i.e. resistant plus intermediate); INFAAR reports percentage resistant without stating how intermediate results were
handled. The comparison is labelled accordingly.

C4. **Ceftriaxone and co-trimoxazole (E. coli).** AMRSN reports these only for urine isolates: ceftriaxone in 2022
only (Table 2.13, OPD, ward and ICU combined); co-trimoxazole in 2020-2022. Urine values for 2020 and 2021 are
printed as images and were transcribed by eye from the "Overall" row (flagged as image-read; percentages are printed
as whole numbers).

C5. **Methicillin resistance in S. aureus.** AMRSN: cefoxitin non-susceptibility (phenotypic). INFAAR aquaculture:
cefoxitin resistance (phenotypic, Figure 2). INFAAR livestock: MRSA confirmed by mecA or mecC detection, reported as a
count of the S. aureus tested. Compared with the definition stated for each.

C6. **Further pairs (section 4.3, "any further pair reported by both networks").** E. coli: ceftazidime
(aquaculture), co-trimoxazole (urine series). S. aureus: penicillin, erythromycin, linezolid, tetracycline,
co-trimoxazole, ciprofloxacin (aquaculture, same drug) and fluoroquinolones as a class (livestock enrofloxacin
against human ciprofloxacin).

C7. **INFAAR denominators.** Per-drug numbers tested are not reported. Wilson intervals for INFAAR use the number of
isolates per source given in the figure legends or text (fisheries legends; livestock: cattle, buffalo, goat, poultry,
S. aureus) and are labelled approximate; no interval is given where no count is reported (sheep, pigs).

C8. **Further shared organisms (section 4.3), identified while auditing dimension 1 and before any comparison for
them was computed.** Coagulase-negative staphylococci (AMRSN: CoNS column of the staphylococcal tables; INFAAR:
Figure 3, fisheries, and Figure 10, livestock, where the panel assignment follows the report's text because the
captions and text disagree; see the extraction validation), Aeromonas species (AMRSN: all-specimens column of the
Aeromonas tables; INFAAR: Figure 7) and Vibrio species other than V. cholerae (AMRSN: "Vibrio spp." column,
flagged for small numbers; INFAAR: Figure 6). V. cholerae (AMRSN) and V. parahaemolyticus (INFAAR) are different
species and are not paired.
