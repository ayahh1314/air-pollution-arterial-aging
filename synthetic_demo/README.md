# Fully synthetic API, arterial age, and MSS demonstration

All participant rows are entirely simulated. No UK Biobank participant records,
identifiers, fitted coefficients, or empirical distributions were used to
generate these inputs. Distribution parameters and relationships are arbitrary
illustrative choices. These files demonstrate the workflow and cannot reproduce
the numerical findings of the manuscript.

## Included data

All six participant tables contain the same 5,000 unique synthetic IDs and can
be joined one-to-one by participant_id. The combined table is a convenience
copy of the split inputs: do not append their rows together.

| File in data/demo | Contents |
| --- | --- |
| synthetic_participants.csv | Combined input to all three modules |
| synthetic_covariates.csv | Age, sex, race, center, education, income, smoking, drinking |
| synthetic_exposures.csv | Four correlated pollutant concentrations in original units |
| synthetic_followup.csv | Follow-up time, death status, recruitment and endpoint dates |
| synthetic_artery_proteins.csv | Five hypothetical proteins and disease_free indicator |
| synthetic_metabolites.csv | Twelve hypothetical metabolite measurements |
| data_dictionary.csv | Every input variable, type, definition, and coding |
| generation_summary.json | Seed, sample counts, versions, and generation provenance |

No missing measurements are included in this first demonstration. This keeps
attention on the analysis rather than introducing a new imputation method.
Five protein placeholders replace the much larger proteomics input. Twelve
metabolites include illustrative signal and noise variables. They have no real
biological identities. disease_free is simulated directly.

## Run the R analysis

The current versions of the three previously prepared analysis modules are
included under scripts. Extract this archive, open R in the project root, and
install dependencies if necessary:

```r
install.packages(c("survival", "caret", "glmnet", "mediation"))
source("run_demo.R")
```

Python is not needed for this step because the CSV data are already included.
The runner calculates API and arterial age gap from the inputs; those derived
quantities are deliberately NOT fabricated as input columns. It then constructs
MSS and fits interactions. R outputs, including sessionInfo(), are written to
results. Mediation uses the retained 200 bootstrap draws and can take time.

Screening thresholds are not relaxed to force a score. If no metabolite passes,
the available results and a status message are saved. Statistical significance
in constructed data is illustrative and not evidence for a scientific claim.

## Regenerate the supplied CSV data

The exact generator is supplied in Python because Python was available in the
preparation environment. It uses Python 3 and NumPy (generated with NumPy 2.3.5):

```sh
python -m pip install numpy==2.3.5
python scripts/00_generate_demo_data.py
```

The output location is based on the script's directory, not the current working
directory. The fixed seed is 20261115. The generator overwrites demonstration
inputs only. It does not inspect or load external data. No repeated seed search
or significance-based regeneration was performed.
