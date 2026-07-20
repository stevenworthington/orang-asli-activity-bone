# dags/

Code for the manuscript's four causal-DAG figures, paired by base name: for each DAG, the dagitty structural definition (`.dag`) and the TikZ source (`.tex`). Rendered figures (`.pdf`/`.png`) live in `outputs/figures/final/`, not here.

| Manuscript figure | TikZ source | Structural definition |
|---|---|---|
| Supplementary Figure 1 (mediator DAG) | `mediator.tex` | `mediator.dag` |
| Figure 2 (simplified mediator DAG) | `mediator-simplified.tex` | `mediator.dag` |
| Supplementary Figure 2 (confounder DAG) | `confounder.tex` | `confounder.dag` |
| Supplementary Figure 3 (industrialization DAG) | `urb-bone.tex` | `urb-bone.dag` |

The mediator structure appears in two figures (full Supp Fig 1 + simplified main-text Fig 2), so `mediator.tex` and `mediator-simplified.tex` share `mediator.dag`. The `.tex` are standalone hand-authored TikZ (XeLaTeX, run twice); they do not read the `.dag`. The `.dag` are the machine-readable dagitty encodings of the same structures.

## Adjustment sets

- **Mediator DAG** (primary, main text): MSAS = {age, sex, pregnancy/lactation, smoking, alcohol, functional status} + village identity (fixed effect).
- **Confounder DAG** (sensitivity, SM8): the mediator MSAS plus fat mass + lean body mass.
- **Industrialization DAG**: MSAS = {age, sex}.

## Archived

Superseded/auxiliary DAG material in `.graveyard/2026-06-23-r-dag-generators/`: the R/ggdag generators (`manuscript-dag-figure.R`, `urb-bone-figure.R`), the measurement-model sidecar (`pa-bone-metadata.yaml`), the retired dual-DAG render, and the abandoned Typst experiment. Prior-framing material is in `.graveyard/2026-05-26-scope-narrowing/dags/`.

## Convention

See `~/.config/agents/CONVENTIONS.md` § "DAGs convention".
