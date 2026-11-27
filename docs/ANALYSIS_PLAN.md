# Analysis plan: do India's human and animal antimicrobial resistance surveillance networks speak the same language?

A harmonisation audit and concordance analysis of the ICMR Antimicrobial Resistance Research and Surveillance
Network (AMRSN) and the ICAR-FAO Indian Network for Fishery and Animal Antimicrobial Resistance (INFAAR), 2019-2022.

Version 1.0, frozen 28 November 2026. Author: Siddalingaiah H S (Shridevi Institute of Medical Sciences and Research
Hospital, Tumkur, India). Any departure from this plan will be reported as a deviation with its reason.

**What had been done at freezing.** A feasibility check (docs/FEASIBILITY.md) extracted, as a test of method, the
E. coli resistance values for six food-animal species from INFAAR Figure 11, and listed the antimicrobials each
network reports for E. coli. No human-animal comparison of resistance levels had been computed. The author has
previously used an automated extract of AMRSN reports in another paper; that extract is not used here (see 4.1).

## 1. Background and aim

India's National Action Plan on Antimicrobial Resistance 2.0 (2025-2029) and the National One Health Mission call
for integrated surveillance across humans, animals and food. Integration requires that the networks measure the same
things in comparable ways. India has two national networks: AMRSN for human clinical isolates and INFAAR for
fisheries, livestock and poultry. No peer-reviewed study has examined INFAAR's data or compared the two networks.

Aims:
1. Audit how far AMRSN and INFAAR are harmonised in organisms, antimicrobial panels, testing and interpretation,
   sampling, isolate handling, reporting and quality assurance.
2. For the organism-antimicrobial pairs both networks report, describe resistance in human clinical isolates beside
   resistance in each animal and fishery source, 2019-2022.
3. Derive a minimum common panel that would allow routine integrated reporting.

## 2. Sources

| Source | Content | Years |
|---|---|---|
| ICMR AMRSN annual reports, 2019, 2020, 2021, 2022 | Human clinical isolates from network centres | 2019-2022 |
| ICAR-FAO INFAAR analytical report (FAO and ICAR, 2024; doi:10.4060/cc9303en) | Fisheries, livestock and poultry isolates | 2019-2022, pooled |

All are public documents. Each value used is recorded with its file, page and table or figure.

## 3. Part A: harmonisation audit

Seven dimensions, each rated for the two networks as **aligned**, **partly aligned**, **not aligned** or
**not reported**, with the supporting text quoted and its page recorded:

1. Organisms under surveillance.
2. Antimicrobial panel for each shared organism; overlap reported as counts and the Jaccard index.
3. Susceptibility testing method (disc diffusion, MIC, automated system) and interpretive standard with version
   (for example CLSI M100 for human isolates, CLSI VET01 for animal isolates, EUCAST), including the handling of
   intermediate results.
4. Sampling frame: clinical specimens from patients against samples from healthy animals, food or the environment;
   sites and geographical coverage.
5. Isolate selection and de-duplication (for example first isolate per patient or per sample).
6. Reporting: resistance as a percentage with or without denominators, stratification (specimen, site, species),
   time period, and public availability of isolate-level data.
7. Quality assurance: external quality assessment, reference laboratories, reference strains.

## 4. Part B: concordance for shared organism-antimicrobial pairs

### 4.1 Data extraction

- **AMRSN:** values are extracted afresh from the annual report PDFs (not from the earlier automated extract) for
  each shared pair and year: percentage resistant and the number of isolates tested, for all specimens combined where
  reported; otherwise the specimen category with the most isolates, stated for each value. Extraction uses the PDF
  text layer; any value read from a figure is flagged.
- **INFAAR:** bar heights are extracted from the vector geometry of each chart (bar top relative to the axis zero and
  100% marks); raster charts are digitised from the rendered page image. Denominators are taken from the figure
  titles or tables where given.

### 4.2 Validation of the extraction (criteria fixed in advance)

- Every percentage stated in the INFAAR report text is compared with the value extracted from the corresponding bar.
  The extraction is accepted if at least 90% of comparisons agree within 2 percentage points; every disagreement is
  investigated and reported.
- An independent second extraction by a different method (pixel-based reading of the rendered page) is made for all
  INFAAR bars; agreement is reported as the mean and maximum absolute difference.
- AMRSN values are checked against the report's own text where the text repeats them.

### 4.3 Pairs

- E. coli: amikacin, cefotaxime, ceftriaxone, imipenem (both networks); fluoroquinolones as a class only
  (ciprofloxacin in humans, enrofloxacin in animals), labelled as a class-level comparison.
- S. aureus: methicillin resistance, and any antimicrobial reported by both networks.
- Any further pair reported by both networks for 2019-2022 is added and listed.

### 4.4 Analysis

- For each pair: AMRSN percentage resistant by year and pooled over 2019-2022 (weighted by isolates where
  denominators are reported); INFAAR percentage resistant by source (each food-animal species; freshwater fish,
  shrimp and marine fish), with Wilson 95% confidence intervals where denominators are available.
- Difference (human minus animal) in percentage points for each pair and source.
- Ordering of antimicrobials by resistance within each network, compared descriptively.
- Interpretation is descriptive. Clinical isolates from patients and isolates from healthy animals come from different
  sampling frames and, possibly, different interpretive standards, so differences in level are not evidence of
  transmission or of its direction.

## 5. Outputs

- Table: harmonisation matrix (7 dimensions x 2 networks, rating and quoted evidence).
- Figure: organism-antimicrobial overlap between the networks.
- Table: concordance for shared pairs.
- Table: extraction validation.
- A proposed minimum common panel for integrated reporting.

## 6. Target journal

Indian Journal of Medical Research (no article processing charge).
