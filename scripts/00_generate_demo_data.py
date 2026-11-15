"""Generate fully synthetic inputs, with no real participant records.

Run from any directory: python scripts/00_generate_demo_data.py
Requires Python 3 and NumPy. All coefficients below are arbitrary illustrative
choices, not estimates obtained from UK Biobank or the manuscript.
"""
from pathlib import Path
import csv
import json
from datetime import date, timedelta
import numpy as np

SEED = 20261115
N = 5000
END = date(2023, 2, 1)
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "demo"


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    pid = np.array([f"SYN{i:06d}" for i in range(1, N + 1)])
    age = rng.uniform(40, 70, N)
    sex = rng.integers(0, 2, N)
    race = rng.choice([1, 2, 3, 4, 5], N, p=[.72, .07, .10, .07, .04])
    center = rng.choice([10001, 10002, 11001, 11002, 11003, 11004,
                         11022, 11023, 10003], N,
                        p=[.15, .15, .15, .15, .15, .16, .03, .03, .03])
    income = rng.choice([1, 2, 3, 4, 5, 6], N, p=[.15, .22, .24, .18, .13, .08])
    education = rng.choice([1, 2, 3], N, p=[.55, .40, .05])
    # Demo-only categories. They do not assert the coding of the real dataset.
    smoking = rng.choice([0, 1, 2], N, p=[.58, .30, .12])
    drinking = rng.choice([0, 1, 2], N, p=[.15, .30, .55])
    air = rng.normal(0, 1, N)
    vulnerability = rng.normal(0, 1, N)
    pm25 = np.maximum(4, 10 + 1.5 * air + rng.normal(0, .35, N))
    coarse = np.maximum(1, 7 + .7 * air + rng.normal(0, .5, N))
    pm10 = pm25 + coarse
    no2 = np.maximum(3, 25 + 6 * air + rng.normal(0, 1.5, N))
    nox = no2 + np.maximum(1, 14 + 3 * air + rng.normal(0, 1.5, N))
    # Simulate a prepared disease-free indicator, not clinical diagnoses/dates.
    disease_free = rng.binomial(1, .65, N)
    recruitment = [date(2006, 1, 1) + timedelta(days=int(x))
                   for x in rng.integers(0, 1826, N)]
    max_followup = np.array([(END - d).days for d in recruitment])
    # Illustrative exponential survival with positive pollution association.
    log_risk = (.055 * (age - 55) + .30 * sex + .35 * air
                + .20 * (smoking == 1) + .45 * (smoking == 2))
    event_wait = np.maximum(1, np.ceil(
        rng.exponential(1 / (.000035 * np.exp(log_risk))))).astype(int)
    # Random loss dates plus the shared administrative end date.
    lost = rng.binomial(1, .04, N).astype(bool)
    loss_wait = rng.integers(365, 4000, N)
    censor_wait = np.where(lost, np.minimum(max_followup, loss_wait), max_followup)
    death_event = (event_wait <= censor_wait).astype(int)
    followup = np.minimum(event_wait, censor_wait)
    endpoint = [d + timedelta(days=int(t)) for d, t in zip(recruitment, followup)]

    # Shared latent aging perturbation makes demonstration proteins partly
    # age-related and partly associated with pollution/metabolic variation.
    # These constructed associations are not scientific findings.
    perturbation = (1.7 * air + 1.8 * vulnerability
                    + .55 * air * np.tanh(vulnerability))
    cols = {
        "participant_id": pid, "age": age, "sex": sex, "race": race,
        "center": center, "income": income, "education": education,
        "smoking": smoking, "drinking": drinking, "pm25": pm25,
        "pm10": pm10, "no2": no2, "nox": nox,
        "followup_days": followup, "death_event": death_event,
        "disease_free": disease_free,
        "recruitment_date": np.array([d.isoformat() for d in recruitment]),
        "followup_end_date": np.array([d.isoformat() for d in endpoint]),
        "administrative_end_date": np.repeat(END.isoformat(), N),
    }
    for j, loading in enumerate([.10, .13, .16, .12, .14], 1):
        cols[f"artery_protein_{j}"] = (
            loading * (age - 55 + perturbation) + .08 * sex
            + rng.normal(0, .45, N))
    # Eight signal metabolites, four noise metabolites; no threshold chasing.
    # The same fixed seed is used once, not resampled until tests are significant.
    for j in range(1, 13):
        if j <= 8:
            sign = 1 if j <= 6 else -1
            value = (sign * (.60 * air + .75 * vulnerability)
                     + .01 * (age - 55) + rng.normal(0, .50, N))
        else:
            value = rng.normal(0, 1, N)
        cols[f"metabolite_{j}"] = value

    def write_table(filename, names):
        with (OUT / filename).open("w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)
            writer.writerow(names)
            for i in range(N):
                row = []
                for name in names:
                    v = cols[name][i]
                    row.append(format(float(v), ".12g") if isinstance(v, np.floating)
                               else v.item() if isinstance(v, np.generic) else v)
                writer.writerow(row)

    covs = ["participant_id", "age", "sex", "race", "center", "income",
            "education", "smoking", "drinking"]
    write_table("synthetic_participants.csv", list(cols))
    write_table("synthetic_covariates.csv", covs)
    write_table("synthetic_exposures.csv", ["participant_id", "pm25", "pm10", "no2", "nox"])
    write_table("synthetic_followup.csv", ["participant_id", "followup_days", "death_event",
        "recruitment_date", "followup_end_date", "administrative_end_date"])
    write_table("synthetic_artery_proteins.csv", ["participant_id", "disease_free"] +
                [f"artery_protein_{j}" for j in range(1, 6)])
    write_table("synthetic_metabolites.csv", ["participant_id"] +
                [f"metabolite_{j}" for j in range(1, 13)])

    definitions = {
        "participant_id": ("string", "Fully synthetic unique identifier", "SYN prefix; no UKB identifiers"),
        "age": ("numeric", "Age at recruitment", "years; 40-70"),
        "sex": ("integer", "Sex", "0=female; 1=male"),
        "race": ("integer", "Race category", "1=White; 2=Mixed; 3=Asian; 4=Black; 5=Other"),
        "center": ("integer", "Illustrative assessment center", "categorical; excluded centers retained for demonstration"),
        "income": ("integer", "Income category", "1-5=ordered categories; 6=missing/unknown category"),
        "education": ("integer", "Prepared education category", "1=lower; 2=upper per original mapping; 3=missing/other"),
        "smoking": ("integer", "Demo smoking category", "0=never; 1=former; 2=current; demo mapping only"),
        "drinking": ("integer", "Demo drinking category", "0=never; 1=less than weekly; 2=at least weekly; demo mapping only"),
        "followup_days": ("integer", "Time from recruitment to death/censoring", "days; positive"),
        "death_event": ("integer", "All-cause death status", "0=censored; 1=death"),
        "disease_free": ("integer", "Simulated prepared training eligibility", "0=ineligible; 1=eligible; not computed from disease dates"),
    }
    for name in ["pm25", "pm10", "no2", "nox"]:
        definitions[name] = ("numeric", f"Synthetic annual {name} concentration", "micrograms/m^3")
    for name in ["recruitment_date", "followup_end_date", "administrative_end_date"]:
        definitions[name] = ("date", name.replace("_", " "), "YYYY-MM-DD; follow-up ends by 2023-02-01")
    for j in range(1, 6):
        definitions[f"artery_protein_{j}"] = ("numeric", "Hypothetical artery-related protein", "arbitrary assay units; no real protein identity")
    for j in range(1, 13):
        definitions[f"metabolite_{j}"] = ("numeric", "Hypothetical preprocessed metabolite", "arbitrary centered units; not forced to sample SD=1")
    with (OUT / "data_dictionary.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["variable", "type", "description", "units_or_coding"])
        for name in cols:
            writer.writerow([name, *definitions[name]])
    summary = {
        "fully_synthetic": True, "seed": SEED, "n": N,
        "numpy_version": np.__version__, "protein_count": 5, "metabolite_count": 12,
        "deaths": int(death_event.sum()), "disease_free": int(disease_free.sum()),
        "administrative_end_date": END.isoformat(),
        "excluded_center_rows": int(np.isin(center, [11022, 11023, 10003]).sum()),
        "missing_values": 0,
        "source": "Fixed distributions and arbitrary coefficients; no UK Biobank data",
    }
    (OUT / "generation_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
