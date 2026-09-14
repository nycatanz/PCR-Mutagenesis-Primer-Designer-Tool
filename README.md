# PCR Mutagenesis Primer Design

Back-to-back primers for site-directed mutagenesis by whole-plasmid PCR —
KOD One, then KLD.

**Live tool:** https://nycatanz.github.io/PCR-Mutagenesis-Primer-Designer-Tool/

Melting temperatures are on the KOD One scale: SantaLucia (1998) nearest-neighbour
with the Owczarzy (2008) divalent correction, at 500 nM oligo / 50 mM Na⁺ /
1.5 mM Mg²⁺. The same primer read *without* Mg²⁺ is about 5.7 °C lower. Set your annealing temperature from these
numbers, not from another calculator's.

> **Not yet validated at the bench.** No primer designed by this tool has been
> ordered and sequenced. Treat the output as a careful first draft.

## How it works

`design_primers.py` is ordinary Python, standard library only. The web page runs
it unmodified in your browser via [Pyodide](https://pyodide.org) — CPython
compiled to WebAssembly. There is no server and nothing is uploaded anywhere;
your sequence never leaves your machine.

That also means there is exactly one implementation of the melting-temperature
model. The page is a form; all the science is in the Python.

## Files

| | |
|---|---|
| `index.html` | the page |
| `design_primers.py` | the design script, run as-is |
| `bglb.gb` | BglB in pET29b(+) — residue 1 at base 5243 |
| `ldh.gb` | LDH in pET29b(+) |

## Command line

The script works on its own, with more control than the page exposes:

```bash
python3 design_primers.py bglb.gb --cds-start 5243 L171N C167A E164Q
python3 design_primers.py bglb.gb --info          # what's in the file
python3 design_primers.py --help                  # every constraint
```

See the full documentation in the project README for the constraint reasoning,
coverage figures, and known limits.
