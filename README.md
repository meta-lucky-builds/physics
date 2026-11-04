# Physics from first principles

171 questions about how the world works, each answered with a proof, a derivation or
code that tests the claim. Where no proof is possible, the answer says so.

The original wording of every question is in [QUESTIONS.md](QUESTIONS.md).

## Sections

| Section | Questions | Scripts | Covers |
|---|---|---|---|
| [01 Forces and rigid bodies](01-forces-and-rigid-bodies/README.md) | F01–F10 | 9 | Doors, friction, moving roads, centre of mass, foundations, tipping, spinning Earth |
| [02 Elasticity and collisions](02-elasticity-and-collisions/README.md) | E01–E09 | 6 | Springs, why things bounce, trampolines, cricket, endless collisions |
| [03 Chemistry, heat and radiation](03-chemistry-thermo-radiation/README.md) | T01–T21 | 8 | Bond energy, inductive effect, Graham's law, fridges, bulb colour, entropy, body heat |
| [04 Light, relativity and quantum](04-light-relativity-quantum/README.md) | L01–L16 | 11 | Mass of light, measuring c, gravity fields, uncertainty, faster-than-light claims, clocks |
| [05 The antimatter engine](05-antimatter-engine/README.md) | A01–A21 | 8 | Pair production, a 38 pK cold sink, the closed annihilation cycle and its energy account |
| [06 Black holes and spacetime](06-black-holes-and-spacetime/README.md) | B01–B48 | 6 | Curvature, accretion, wormholes, entropy, evaporation, the notebook paradox |
| [07 Cosmology, time and mind](07-cosmology-time-mind-foundations/README.md) | C01–C46 | 11 | Expansion, time travel, block universe, Boltzmann brains, thoughts, life elsewhere |

## How to read an answer

Each answer has a status saying what kind of support it has. The labels are defined in
[CONVENTIONS.md](CONVENTIONS.md):

- **Derived** and **Code** mean the result is worked out or tested in this repo.
- **Established** and **Estimate** mean it rests on cited results or stated assumptions.
- **Open**, **Conceptual** and **Speculative** mean it can't be settled by proof today.

Sources that couldn't be fully confirmed are flagged in each section. None of this
has been peer reviewed.

## Run the code

```bash
pip install -r requirements.txt
python tools/check.py
```

`tools/check.py` confirms every question has an answer and runs every verification
script. Add `--no-run` to skip the scripts, which take several minutes in total. You can
also run one script directly, for example:

```bash
python 06-black-holes-and-spacetime/code/verify_notebook.py
```
