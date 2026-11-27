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
