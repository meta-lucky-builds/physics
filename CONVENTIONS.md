# Conventions

These rules apply to every section of this repo.

## Layout

Each section is a folder holding:

- `README.md`, which answers every question assigned to the section, in ID order.
- `code/verify_*.py`, scripts that compute or check results. Each must finish in
  under 2 minutes, print what it checked and exit non-zero if a check fails.
- `figures/` (optional), PNGs under 300 KB that a script generates and the
  README embeds.

Scripts run with `python <section>/code/verify_x.py` from the repo root, using only
the packages in `requirements.txt`. Shared constants live in
[`common/constants.py`](common/constants.py); a script imports them after adding
the repo root to `sys.path`.

## Answer format

Each answer starts with a heading containing its exact ID, followed by these parts:

```markdown
### F05: Centre of mass of a donut

**Status:** Derived + Code

**Short answer:** one to three plain sentences that answer what was asked.

**Starting from:** the assumptions and laws the argument takes as given.

**Derivation:** the step-by-step argument, with every step justified.

**Check:** what the code confirms, with a link to the script and the actual output.

**Limits:** where the result stops applying and what was left out.
```

Omit **Check** when there is no code. Questions that share one derivation may each
point to it, but each still needs its own heading, status and direct answer.

## Status labels

Use one or more of these, and only when true:

| Label | Meaning |
|---|---|
| Derived | Worked through here from stated principles. |
| Code | Computed or simulated by a script in this repo that passes when run. |
| Established | Standard or measured result that is cited, not derived here. |
| Estimate | Order-of-magnitude reasoning with stated inputs. |
| Open | Physicists or scientists don't have an agreed answer. Say what is known and what isn't. |
| Conceptual | Not provable by derivation: a question of definition, interpretation or philosophy. Give the strongest argument and say it isn't a proof. |
| Speculative | Goes beyond tested physics. |

## Rigor rules

1. Never present something as proved when it isn't. Many questions here have no
   proof, and "Open" or "Conceptual" is the right answer for those.
2. "First principles" means naming what you start from, because every derivation has
   axioms. State them.
3. Correct wrong premises directly, then answer the best version of the question.
   Also credit what the questioner got right.
4. Code must test the claim independently, for example by simulating the system
   rather than re-evaluating the formula being checked. Say what the test can't show.
5. Never hard-code a result into a check or weaken a test to make it pass. If a script
   disagrees with a derivation, find out which is wrong and report it.
6. Cite a source for every measured value and result not derived here, using only
   sources you have confirmed exist. Without a confirmed source, say the claim is from
   general knowledge and unverified.
7. Keep every number traceable: give units, inputs and precision that fits the inputs.
8. Math uses GitHub syntax: `$...$` inline and `$$...$$` on separate lines for display.
