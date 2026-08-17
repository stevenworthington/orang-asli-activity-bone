# dags/

Sources for the manuscript's causal-DAG figures, paired by base name: for each DAG, the dagitty structural definition (`.dag`) and the TikZ source (`.tex`).

| Manuscript figure | TikZ source | Structural definition |
|---|---|---|
| Figure 1 (all three hypotheses, 3 panels) | `three-hypotheses-simplified.tex` | all three `.dag` |
| Supplementary Figure 1 (mediator DAG) | `mediator.tex` | `mediator.dag` |
| Supplementary Figure 2 (confounder DAG) | `confounder.tex` | `confounder.dag` |
| Supplementary Figure 3 (industrialization DAG) | `urb-bone.tex` | `urb-bone.dag` |

The mediator structure appears in more than one figure, so `mediator.tex` and panels B/C of the composite both express `mediator.dag`.

**The `.tex` and the `.dag` are independent expressions of the same structures.** The `.tex` files are hand-authored standalone TikZ with tuned node positions; they do not read the `.dag`. The `.dag` files are the machine-readable dagitty encodings, for checking adjustment sets programmatically. Regenerating a figure's layout from its `.dag` would discard the hand-tuning, so the two are kept in step by hand.

## Adjustment sets

- **Mediator DAG** (primary, main text): MSAS = {age, sex, pregnancy/lactation, smoking, alcohol, functional status} + community identity (fixed effect).
- **Confounder DAG** (sensitivity): the mediator MSAS plus fat mass + lean body mass.
- **Industrialization DAG**: MSAS = {age, sex}.

## The composite figure

`three-hypotheses-simplified.tex` states the identification strategy for all three hypotheses in one figure.

**Panel order runs 1, 2, 3.** Panel A (top) carries hypotheses 1 and 2, which share the industrialization exposure. Panels B and C share the lower row because they are the two **competing** structures for the same hypothesis (3) and should be directly comparable — they differ only in the colour of the body-composition node and the direction of its arrows.

A 1×3 side-by-side arrangement is not usable: it comes out ~24 cm wide, which cannot reach a journal page without shrinking the type below legibility.

Node-name suffixes in the source (`ageA`, `plateB`, …) match the displayed panel letters, so a reorder has to touch both or neither — the two cannot silently drift apart. Panel A's bone outcome sits well right of the activity outcome to leave the hypothesis-3 edge label a clear run.

**Panel A's plate and arrows.** Every plate in the figure fits only its label and its adjustment-set nodes, with the same 11pt inner xsep, so all three relate to their contents identically. Panel A's exposure and outcome are therefore positioned by two explicit coordinates (`edgeLA` / `edgeRA`) rather than by anchoring to the plate, since the plate does not span the panel. Its two outer confounding arrows angle out of the plate's bottom corners; the arrow into the activity outcome stays vertical, that node being directly beneath the plate.

Those two angled arrows terminate at `indA.north` and `boneA.north`, **not** at `(indA)` / `(boneA)`. Ending at the node aims the line at the node's centre but clips it at the border, which for a line arriving from inboard lands the arrowhead about 4 mm off centre along the top edge. Naming `.north` is what makes them read as centred.

**Panel identifiers** (`\panelid`) are set at 13.5pt against the 9pt panel title — half again its size — in the same sans face but *not* bold. At that size the letter already dominates the line, and bolding it too makes it compete with the title instead of labelling it. The macro pins `\baselineskip` to the title's own 11pt so the taller glyph does not open the gap to the italic subtitle beneath.

## Legends

Every `.tex` draws its own legend, via the swatch helpers (`\lnode`, `\ledge`, `\ledgedash`, and the `\lentry` / `\eentry` wrappers) defined in its preamble. Three things to preserve if these are edited:

- **The swatches take a style name and resolve against each file's own `\tikzset`**, so a legend cannot drift from what the figure actually draws. Each legend lists only the roles and edge classes that file uses — `urb-bone.tex` has no latent nodes and calls its white nodes mediators, because that is what they are there.
- **`\lentry` / `\eentry` wrap each entry in an `\mbox`** so a line can only break *between* entries. Without it the legend hyphenates mid-phrase ("confound-/ing path").
- **Swatch-vs-text vertical centring is set per swatch type**, by `\swnode`, `\swdash` and `\swrule`. They do **not** share one value, and that is deliberate: the three macros build different objects (a tikz node, a tikz path, a plain rule) and they do not land alike for a given offset in these files' context.

  If they need retuning, **sweep them in the figure itself, never in a reduced test case.** A reduced case lacks these files' global `font=\sffamily\footnotesize` in `\tikzset`, which the nested `\tikz` swatches inherit; it will happily centre all three there while leaving the node swatches about 3 pt low here. Measure by rendering at high resolution and comparing each swatch's centre against its own adjacent label, taking the text band nearest the swatch (a fixed-height window catches the neighbouring legend line and hides the error). Allow for each label's own bias: a word with a descender sits lower than its x-height centre, one with ascenders and none sits higher.

`\setsansfont` names `BoldFont={Inter 18pt SemiBold}` explicitly, because plain Inter ships no bold shape and `\textbf` otherwise falls back silently to regular. Any style a legend swatch uses must be declared in a preamble `\tikzset` rather than in a `tikzpicture` option list — a nested `\tikz` does not inherit styles declared on the outer picture.

## Building

XeLaTeX, run twice so the `fit` plates settle:

```sh
xelatex <stem>.tex && xelatex <stem>.tex
```

Requires the Inter font family, including *Inter 18pt SemiBold*. Substitute another sans via `\setsansfont` if Inter is unavailable.
