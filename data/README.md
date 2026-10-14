# Input data: not redistributed

The original grouped-binomial dataset was supplied for undergraduate MATH3091 coursework. Public redistribution rights have **not** been confirmed, so `immune.csv` is intentionally **not included** in this GitHub-ready portfolio.

To run the full project, place your own authorised copy at:

`data/immune.csv`

Expected columns (case-sensitive):

| Column | Meaning |
| --- | --- |
| `immune` | Positive biomarker level |
| `disease` | Number of people with the disease at that biomarker level |
| `nodisease` | Number of people without the disease at that level |

The coursework dataset contains **12 aggregated biomarker-level rows**, representing **418 individuals** in total. It is grouped count data, not an individual-level patient file.

The scripts use `cbind(disease, nodisease)` as the response to a binomial GLM, with each row representing a group of observations.

If you do not have permission to distribute the dataset, keep `immune.csv` only on your own computer and do not add it to GitHub. The provided `.gitignore` excludes `/data/*.csv`.
