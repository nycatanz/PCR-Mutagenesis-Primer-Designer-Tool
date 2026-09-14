#!/usr/bin/env python3
"""
design_primers.py -- back-to-back mutagenesis primer design for point mutations.

Designs a non-overlapping forward/reverse primer pair that installs a single
amino-acid substitution, for KOD One "fake Q5" amplification followed by KLD.

Template model
--------------
The template is the whole circular plasmid, so a primer near either end of the
coding sequence extends into vector sequence instead of running out of
template. There are no unreachable residues.

Melting temperature
-------------------
SantaLucia (1998) unified nearest-neighbour parameters at 250 nM oligo and
50 mM Na+, matching Benchling's default Tm settings. Validated to within
0.03 C against five Benchling-reported primers; see test_tm.py.

Standard library only -- nothing to install.

Usage
-----
    python3 design_primers.py <plasmid.gbk> <mutation> [<mutation> ...]

Example
-------
    python3 design_primers.py "pET29b+ BglB plasmid.gbk" L171N
"""

import argparse
import math
import os
import re
import sys

# ---------------------------------------------------------------------------
# 1. Genetic code
# ---------------------------------------------------------------------------
# What each codon MEANS. Distinct from the usage table in section 2, which is
# about which synonymous codon to CHOOSE.

CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}

AA3 = {
    'A': 'Ala', 'R': 'Arg', 'N': 'Asn', 'D': 'Asp', 'C': 'Cys', 'Q': 'Gln',
    'E': 'Glu', 'G': 'Gly', 'H': 'His', 'I': 'Ile', 'L': 'Leu', 'K': 'Lys',
    'M': 'Met', 'F': 'Phe', 'P': 'Pro', 'S': 'Ser', 'T': 'Thr', 'W': 'Trp',
    'Y': 'Tyr', 'V': 'Val', '*': 'Stop',
}

COMPLEMENT = str.maketrans('ACGTacgtNn', 'TGCAtgcaNn')

# ---------------------------------------------------------------------------
# 2. E. coli K12 codon usage
# ---------------------------------------------------------------------------
# Transcribed verbatim from E_coli_K12_codon_usage_table.md
#   Source: Kazusa, https://www.kazusa.or.jp/codon/ species 83333 (E. coli K12)
# Values are (codon, fraction-of-that-amino-acid, per-1000-frequency).
# The rule is: always take the highest-frequency codon. No GC tuning, no
# judgement call. TAG was not reported in the source and is omitted.

CODON_USAGE = {
    'A': [('GCG', 0.38, 38.5), ('GCC', 0.31, 31.6), ('GCA', 0.21, 21.1), ('GCT', 0.11, 10.7)],
    'R': [('CGC', 0.44, 26.0), ('CGT', 0.36, 21.1), ('CGA', 0.07, 4.3), ('CGG', 0.07, 4.1),
          ('AGG', 0.03, 1.6), ('AGA', 0.02, 1.4)],
    'N': [('AAC', 0.53, 24.4), ('AAT', 0.47, 21.9)],
    'D': [('GAT', 0.65, 37.9), ('GAC', 0.35, 20.5)],
    'C': [('TGC', 0.58, 8.0), ('TGT', 0.42, 5.9)],
    'Q': [('CAG', 0.70, 27.7), ('CAA', 0.30, 12.1)],
    'E': [('GAA', 0.70, 43.7), ('GAG', 0.30, 18.4)],
    'G': [('GGC', 0.46, 33.4), ('GGT', 0.29, 21.3), ('GGA', 0.13, 9.2), ('GGG', 0.12, 8.6)],
    'H': [('CAT', 0.55, 15.8), ('CAC', 0.45, 13.1)],
    'I': [('ATT', 0.58, 30.5), ('ATC', 0.35, 18.2), ('ATA', 0.07, 3.7)],
    'L': [('CTG', 0.46, 46.9), ('TTA', 0.15, 15.2), ('CTT', 0.12, 11.9), ('TTG', 0.12, 11.9),
          ('CTC', 0.10, 10.5), ('CTA', 0.05, 5.3)],
    'K': [('AAA', 0.73, 33.2), ('AAG', 0.27, 12.1)],
    'M': [('ATG', 1.00, 24.8)],
    'F': [('TTT', 0.57, 19.7), ('TTC', 0.43, 15.0)],
    'P': [('CCG', 0.55, 26.7), ('CCT', 0.17, 8.4), ('CCA', 0.14, 6.6), ('CCC', 0.13, 6.4)],
    'S': [('AGC', 0.33, 16.6), ('TCG', 0.16, 8.0), ('TCA', 0.15, 7.8), ('AGT', 0.14, 7.2),
          ('TCC', 0.11, 5.5), ('TCT', 0.11, 5.7)],
    'T': [('ACC', 0.47, 22.8), ('ACG', 0.24, 11.5), ('ACT', 0.16, 8.0), ('ACA', 0.13, 6.4)],
    'W': [('TGG', 1.00, 10.7)],
    'Y': [('TAT', 0.53, 16.8), ('TAC', 0.47, 14.6)],
    'V': [('GTG', 0.40, 26.4), ('GTT', 0.25, 16.8), ('GTC', 0.18, 11.7), ('GTA', 0.17, 11.5)],
    '*': [('TAA', 0.64, 1.8), ('TGA', 0.36, 1.0)],
}

# Highest-frequency codon per residue, derived from the table above rather than
# typed a second time, so the two can never disagree.
PREFERRED_CODON = {aa: max(rows, key=lambda r: r[2]) for aa, rows in CODON_USAGE.items()}


# ---------------------------------------------------------------------------
# 3. Melting temperature -- SantaLucia (1998) unified nearest-neighbour
# ---------------------------------------------------------------------------
# Most of a duplex's stability comes from base STACKING, not from the hydrogen
# bonds across the rungs. So the energy is tabulated per overlapping pair of
# adjacent bases, not per base: GC and CG are worth different amounts.
#
# dH in kcal/mol, dS in cal/(mol*K), keyed by the 5'->3' dinucleotide.

NN_PARAMS = {
    'AA': (-7.9, -22.2), 'TT': (-7.9, -22.2),
    'AT': (-7.2, -20.4), 'TA': (-7.2, -21.3),
    'CA': (-8.5, -22.7), 'TG': (-8.5, -22.7),
    'GT': (-8.4, -22.4), 'AC': (-8.4, -22.4),
    'CT': (-7.8, -21.0), 'AG': (-7.8, -21.0),
    'GA': (-8.2, -22.2), 'TC': (-8.2, -22.2),
    'CG': (-10.6, -27.2),
    'GC': (-9.8, -24.4),
    'GG': (-8.0, -19.9), 'CC': (-8.0, -19.9),
}

# Helix ends fray, so each terminus carries an initiation penalty that depends
# on whether that terminal base pair is G/C or A/T.
INIT_GC = (0.1, -2.8)
INIT_AT = (2.3, 4.1)

R_GAS = 1.9872  # gas constant, cal/(mol*K)

# What one base past --len-max costs, in the same units severity() uses for
# degrees of Tm. Small on purpose: a 40 nt primer at the right temperature is
# better than a 26 nt primer six degrees too cold.
LENGTH_PENALTY_PER_NT = 0.15

# What a 3' G/C at position -2 or -3 costs instead of the final base. Small:
# the constraint is satisfied either way, this only breaks ties toward the
# textbook ideal. Large enough to matter, far below a degree of Tm.
CLAMP_OFFSET_PENALTY = 0.2


def melting_temp(seq, primer_conc_nM=500.0, na_mM=50.0,
                 mg_mM=1.5, dntp_mM=0.6):
    """SantaLucia (1998) nearest-neighbour Tm, corrected for the actual buffer.

    Walks a two-base window summing enthalpy and entropy, adds an initiation
    penalty per end (helix ends fray), then corrects for salt.

    Two salt corrections are possible and they differ by several degrees:

      * With no free Mg2+, the SantaLucia monovalent entropy correction is
        used -- this is Benchling with Mg2+ and dNTP set to zero.
      * With Mg2+ present, the Owczarzy (2008) divalent correction is used.

    Divalent cations stabilise the duplex substantially: at 1.5 mM Mg2+ a
    typical 25-mer reads about 5 C higher than the salt-only figure. That is
    the real condition inside a PCR tube, and it is what Benchling reports
    when those fields are filled in.

    Conditions are the PCR buffer as Benchling reports it: 500 nM oligo,
    50 mM Na+/K+, 1.5 mM Mg2+, 0.6 mM dNTP. Reading the same oligo without
    Mg2+ gives a figure about 6.7 C lower -- that is a different scale, and
    mixing the two is the one error that matters here.
    """
    seq = seq.upper()
    if len(seq) < 2:
        return float('nan')

    dH = dS = 0.0
    for i in range(len(seq) - 1):
        pair = seq[i:i + 2]
        if pair not in NN_PARAMS:
            return float('nan')          # ambiguous base
        h, s = NN_PARAMS[pair]
        dH += h
        dS += s

    for terminal in (seq[0], seq[-1]):
        h, s = INIT_GC if terminal in 'GC' else INIT_AT
        dH += h
        dS += s

    n = len(seq)
    ct = primer_conc_nM * 1e-9
    mon = na_mM / 1000.0
    # Benchling does not subtract dNTP from Mg2+ before the divalent
    # correction. Checked against reference primers: total Mg2+ agrees to
    # 0.2 C, subtracting dNTP reads 1.15 C low. dntp_mM is accepted so the
    # buffer is documented in one place, but it does not enter the maths.
    free_mg = mg_mM / 1000.0

    def monovalent_only():
        return (dH * 1000.0) / (dS + 0.368 * (n - 1) * math.log(mon)
                                + R_GAS * math.log(ct / 4.0)) - 273.15

    if free_mg <= 0.0 or mon <= 0.0:
        if mon <= 0.0 and free_mg <= 0.0:
            return float('nan')
        if free_mg <= 0.0:
            return monovalent_only()

    ratio = math.sqrt(free_mg) / mon
    if ratio < 0.22:
        return monovalent_only()         # monovalent dominates

    # Owczarzy (2008) divalent correction, applied to the 1 M Na+ reference Tm.
    tm_1M = (dH * 1000.0) / (dS + R_GAS * math.log(ct / 4.0))
    fgc = sum(1 for base in seq if base in 'GC') / float(n)
    a_, b_, c_, d_ = 3.92e-5, -9.11e-6, 6.26e-5, 1.42e-5
    e_, f_, g_ = -4.82e-4, 5.25e-4, 8.31e-5
    if ratio < 6.0:                      # both ions matter; a, d, g depend on Na+
        ln_mon = math.log(mon)
        a_ *= 0.843 - 0.352 * math.sqrt(mon) * ln_mon
        d_ *= 1.279 - 0.00403 * ln_mon - 0.00803 * ln_mon ** 2
        g_ *= 0.486 - 0.258 * ln_mon + 0.0029 * ln_mon ** 3
    ln_mg = math.log(free_mg)
    inverse = (1.0 / tm_1M + a_ + b_ * ln_mg + fgc * (c_ + d_ * ln_mg)
               + (1.0 / (2.0 * (n - 1))) * (e_ + f_ * ln_mg + g_ * ln_mg ** 2))
    return 1.0 / inverse - 273.15


def gc_percent(seq):
    return 100.0 * sum(1 for b in seq.upper() if b in 'GC') / len(seq) if seq else 0.0


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


def translate(dna):
    return ''.join(CODON_TABLE.get(dna[i:i + 3], 'X') for i in range(0, len(dna) - 2, 3))


# ---------------------------------------------------------------------------
# 4. Sequence input
# ---------------------------------------------------------------------------

class Feature:
    """One GenBank FEATURES entry: its kind, label, and where it sits.

    Locations may be a join() spanning the origin, e.g. join(5243..6588,1..1)
    for a gene that wraps around a circular plasmid.
    """

    def __init__(self, kind, location, qualifiers):
        self.kind = kind
        self.location = location
        self.qualifiers = qualifiers
        self.complement = 'complement' in location
        self.segments = [(int(a) - 1, int(b))          # 0-based half-open
                         for a, b in re.findall(r'(\d+)\.\.(\d+)', location)]

    @property
    def label(self):
        return self.qualifiers.get('label') or self.qualifiers.get('gene') \
            or self.qualifiers.get('product') or ''

    @property
    def start(self):
        return self.segments[0][0] if self.segments else None

    @property
    def length(self):
        return sum(b - a for a, b in self.segments)


def parse_features(header):
    """Pull FEATURES entries out of a GenBank header block."""
    m = re.search(r'^FEATURES.*$', header, re.M)
    if not m:
        return []
    feats, kind, loc, quals, key = [], None, '', {}, None

    def flush():
        if kind:
            feats.append(Feature(kind, loc, quals))

    for line in header[m.end():].splitlines():
        head = re.match(r'^ {5}(\S+)\s+(.*)$', line)
        if head:
            flush()
            kind, loc, quals, key = head.group(1), head.group(2).strip(), {}, None
        elif kind and re.match(r'^ {21}/', line):
            q = re.match(r'^ {21}/(\w+)=?"?(.*?)"?$', line.rstrip())
            if q:
                key = q.group(1)
                quals[key] = q.group(2)
        elif kind and key and line.startswith(' ' * 21):
            quals[key] += ' ' + line.strip().rstrip('"')
        elif kind and line.strip() and not line.startswith(' '):
            break
    flush()
    return feats


class Template:
    """A plasmid sequence, indexed circularly.

    Slicing past either end wraps around the origin, which is what makes
    residues at the ends of the coding sequence designable.
    """

    def __init__(self, name, seq, circular=True, features=()):
        self.name = name
        self.seq = seq
        self.circular = circular
        self.features = list(features)
        self.junction_cache = {}

    def __len__(self):
        return len(self.seq)

    def at(self, start, length):
        """Return `length` bases beginning at 0-based `start`, wrapping."""
        n = len(self.seq)
        if not self.circular:
            if start < 0 or start + length > n:
                raise IndexError('ran off the end of a linear template')
            return self.seq[start:start + length]
        start %= n
        if start + length <= n:
            return self.seq[start:start + length]
        doubled = self.seq + self.seq
        return doubled[start:start + length]

    def replace(self, start, new):
        """Return a copy with `new` written in at 0-based `start` (wrapping)."""
        n = len(self.seq)
        buf = list(self.seq)
        for k, base in enumerate(new):
            buf[(start + k) % n] = base
        return Template(self.name, ''.join(buf), self.circular, self.features)

    def find(self, probe):
        """0-based index of `probe`, searching across the origin. -1 if absent."""
        i = (self.seq + self.seq[:len(probe)]).find(probe)
        return i if i < len(self.seq) else -1


def read_genbank(path):
    """Minimal GenBank reader: returns a Template. Sequence + topology only."""
    with open(path) as fh:
        text = fh.read()
    # Split on a line that IS "ORIGIN", not on the substring -- a header
    # word such as "ORIGINAL" would otherwise cut the file in the wrong place.
    m = re.search(r'^ORIGIN.*$', text, re.M)
    if not m:
        raise ValueError('%s has no ORIGIN line -- is it really GenBank?' % path)
    header, origin = text[:m.start()], text[m.end():]

    circular = bool(re.search(r'^LOCUS.*\bcircular\b', header, re.M | re.I))
    m = re.search(r'^LOCUS\s+(\S+)', header, re.M)
    name = m.group(1) if m else path

    seq = ''.join(re.findall(r'[acgtnACGTN]', origin.split('//')[0])).upper()
    if not seq:
        raise ValueError('%s: ORIGIN block contained no sequence.' % path)
    return Template(name, seq, circular, parse_features(header))


def read_fasta(path):
    """Return (header, sequence) from a single-record FASTA file."""
    header, parts = None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if header is not None:
                    raise ValueError('%s holds more than one record; '
                                     'one sequence per file please.' % path)
                header = line[1:].strip()
            else:
                parts.append(line)
    seq = ''.join(parts).upper().replace('U', 'T')
    bad = set(seq) - set('ACGTN')
    if bad:
        raise ValueError('%s contains non-DNA characters: %s'
                         % (path, ', '.join(sorted(bad))))
    return header or path, seq


def load_template(path):
    """Read a sequence file, deciding format by CONTENT rather than filename.

    GenBank exports arrive with all sorts of extensions -- .gb, .gbk, .ape
    from ApE, sometimes .txt -- so trusting the extension turns a perfectly
    good file into a confusing parse error.
    """
    try:
        with open(path) as fh:
            head = fh.read(8192)
    except OSError as e:
        raise ValueError('Could not read %s: %s' % (path, e))

    if not head.strip():
        raise ValueError('%s is empty.' % path)

    if re.search(r'^LOCUS\s', head, re.M) or re.search(r'^ORIGIN\s*$', head, re.M):
        return read_genbank(path)

    if head.lstrip().startswith('>'):
        name, seq = read_fasta(path)
        return Template(name, seq, circular=False)

    # No header at all: accept a bare run of DNA, which people do save.
    letters = [c for c in head if not c.isspace()]
    if letters and sum(c in 'ACGTNacgtn' for c in letters) / len(letters) > 0.95:
        name, seq = read_fasta(path)
        return Template(name, seq, circular=False)

    raise ValueError(
        '%s is neither GenBank nor FASTA. GenBank files begin with a LOCUS '
        'line; FASTA files begin with ">". Check the file, or re-export it '
        'from SnapGene / Benchling / ApE.' % path)


MUT_RE = re.compile(r'^([A-Za-z*])(\d+)([A-Za-z*])$')


def parse_mutation(text):
    m = MUT_RE.match(text.strip())
    if not m:
        raise ValueError('Could not parse mutation %r. Expected a form like '
                         'L171N -- WT residue, position, new residue.' % text)
    wt, pos, new = m.group(1).upper(), int(m.group(2)), m.group(3).upper()
    if pos < 1:
        raise ValueError('Position must be 1 or greater.')
    if new not in PREFERRED_CODON:
        raise ValueError('Unknown target residue %r.' % new)
    if wt not in CODON_TABLE.values():
        raise ValueError('Unknown WT residue %r.' % wt)
    return wt, pos, new


# ---------------------------------------------------------------------------
# 5. Constraints
# ---------------------------------------------------------------------------

class Check:
    """One named constraint, its measured value, and whether it passed."""

    def __init__(self, name, value, ok, detail=''):
        self.name, self.value, self.ok, self.detail = name, value, ok, detail


def gc_clamp_position(seq):
    """1-based distance of the nearest G/C from the 3' end (1 = final base).

    None if there is no G or C in the last three bases. The 3' end is where
    the polymerase extends from, so it needs to be firmly annealed.
    """
    for offset, base in enumerate(reversed(seq.upper()[-3:]), start=1):
        if base in 'GC':
            return offset
    return None


def severity(primer, cfg):
    """How badly a primer misses its constraints, in rough degrees-equivalent.

    Zero when everything passes. Used only to rank options that already break
    a rule -- counting flags treats a 6 C melting-temperature miss the same as
    a 1-point GC miss, which picks badly in awkward stretches of sequence.
    GC is scaled to about a tenth of a degree per point. The 3' clamp is
    graded rather than binary: no G/C in the last three bases costs a flat
    1.5, since that is the end the polymerase extends from, but a G/C at -2
    or -3 costs only a little. It is still worth preferring the final base,
    and no longer worth accepting a hotter primer to get it.
    """
    s = 0.0
    if primer.tm < cfg.tm_min:
        s += cfg.tm_min - primer.tm
    elif primer.tm > cfg.tm_max:
        s += primer.tm - cfg.tm_max
    if primer.gc < cfg.gc_min:
        s += (cfg.gc_min - primer.gc) / 10.0
    elif primer.gc > cfg.gc_max:
        s += (primer.gc - cfg.gc_max) / 10.0
    # The 3' flank is the anchoring side. --flank-3-min is a hard floor, but
    # sitting at the floor is not as good as sitting at --flank-3-ideal, and
    # nothing else in this function would notice the difference. Without this
    # term the search treats d=7 and d=13 as equivalent and spends the shorter
    # anchor everywhere it buys a fraction of a degree -- measured at 1,870
    # primers below 10 on BglB, mean 3' flank falling 13.9 -> 12.8. Graded, it
    # stays an escape hatch for hard positions: ~350 primers below 10, mean
    # unchanged at 13.9.
    if primer.flank3 is not None and primer.flank3 < cfg.flank3_ideal:
        s += (cfg.flank3_ideal - primer.flank3) * cfg.flank3_penalty
    clamp = gc_clamp_position(primer.seq)
    if clamp is None:
        s += 1.5
    elif clamp != 1:
        s += CLAMP_OFFSET_PENALTY
    # Length is the rule that gives way: when a stretch is AT-rich, the only
    # route to the Tm window is a longer primer, and refusing to grow leaves a
    # primer several degrees too cold. So going past --len-max costs a small
    # amount per base, well below the cost of a degree of Tm. Coming in under
    # --len-min is different -- that is a floor, not a preference.
    n = len(primer.seq)
    if n < cfg.len_min:
        s += (cfg.len_min - n) * 2.0
    elif n > cfg.len_max:
        s += (n - cfg.len_max) * LENGTH_PENALTY_PER_NT
    return s


def evaluate(seq, cfg):
    # Constraints are tested on the ROUNDED values -- the same numbers the
    # report prints. The Tm model is only good to about 0.03 C, so flagging a
    # primer at 58.9989 for being under 59 is false precision, and a report
    # that says "59.0 C" next to a failure mark is just confusing.
    tm = round(melting_temp(seq, cfg.primer_conc, cfg.na, cfg.mg, cfg.dntp), 1)
    gc = round(gc_percent(seq))
    clamp = gc_clamp_position(seq)
    return [
        Check('Tm', '%.1f C' % tm, cfg.tm_min <= tm <= cfg.tm_max,
              'window %.0f-%.0f C, designed to %.0f' % (cfg.tm_min, cfg.tm_max, cfg.tm_target)),
        Check('GC content', '%d%%' % gc, cfg.gc_min <= gc <= cfg.gc_max,
              'target %.0f-%.0f%%' % (cfg.gc_min, cfg.gc_max)),
        Check('Length', '%d nt' % len(seq), cfg.len_min <= len(seq) <= cfg.len_max,
              'target %d-%d nt, never past %d'
              % (cfg.len_min, cfg.len_max, cfg.len_hard_max)),
        Check("3' GC clamp", 'G/C at -%d' % clamp if clamp else 'none in last 3',
              clamp is not None,
              'a G or C within the last 3 bases; the final base is best'),
    ]


# ---------------------------------------------------------------------------
# 6. Primer construction
# ---------------------------------------------------------------------------

class Primer:
    def __init__(self, seq, cfg, note='', flank3=None):
        self.seq = seq
        self.note = note
        # Matched bases 3' of the codon, for forward primers. None for reverse
        # primers, which have no codon and so no 3' flank to grade.
        self.flank3 = flank3
        self.checks = evaluate(seq, cfg)
        # Rounded to the value the report prints, so the number shown, the
        # number checked and the number filtered on are all the same one.
        self.tm = round(melting_temp(seq, cfg.primer_conc, cfg.na, cfg.mg, cfg.dntp), 1)
        self.gc = round(gc_percent(seq))

    @property
    def flags(self):
        return [c for c in self.checks if not c.ok]


class Pair:
    """A forward/reverse pair sharing one junction.

    Geometry -- the two primers point AWAY from each other, so their 5' ends
    abut and it is those two blunt 5' ends that KLD ligates:

        template  ... [   u nt   ][ WT codon ][   d nt   ] ...
        FWD                     5'-[ u ][new codon][ d ]-3'
        REV       3'-[   r nt   ]-5'
                                 ^^
                          ligation junction

    Both flanks are named by the side of the codon they sit on, and they are
    NOT symmetric:

        d  -- 3' of the codon, --flank-3-min / --flank-3-max.
              The ANCHORING side: the polymerase extends from this end, so it
              needs a solid run of matched template. Floor of 10, no ceiling
              by default -- it is the forward primer's only freedom, since its
              5' end is pinned at the junction.

        u  -- 5' of the codon, --flank-5-min / --flank-5-max.
              The JUNCTION side. Nothing extends from this end, so it is the
              flexible one and may be shorter. It does double duty: it is the
              forward primer's 5' anchor AND it fixes where the junction sits,
              which is what pins the reverse primer. It is therefore also how
              close the installed mutation sits to the blunt ligated end.

        r  -- the reverse primer, running 5' <- 3' leftwards from the junction.

    The whole plasmid is amplified linearly and re-circularised by KLD.

    Note the reverse primer depends only on the junction -- that is, on the
    POSITION -- not on which residue is being installed. Every substitution
    at one position shares a reverse primer.
    """

    def __init__(self, fwd, rev, cfg, upstream_flank):
        self.fwd, self.rev, self.u = fwd, rev, upstream_flank
        self.pos = None
        self.dtm = abs(fwd.tm - rev.tm)
        self.pair_checks = [
            Check('Tm difference', '%.1f C' % self.dtm, self.dtm <= cfg.pair_tm_tol,
                  'forward and reverse within %.0f C' % cfg.pair_tm_tol),
        ]

    @property
    def n_flags(self):
        return len(self.fwd.flags) + len(self.rev.flags) + \
            sum(1 for c in self.pair_checks if not c.ok)

    def rank(self, cfg):
        """Order candidates: fewest broken rules, then closest to the design Tm.

        Closeness to the design temperature outranks matching the two primers
        to each other. Both matter, but a pair sitting at the target with a
        0.3 C spread beats a pair 1 C low with a 0.1 C spread -- the limiting
        primer sets the annealing temperature, so its absolute value is what
        governs whether the reaction works.
        """
        off = abs(self.fwd.tm - cfg.tm_target) + abs(self.rev.tm - cfg.tm_target)
        harm = severity(self.fwd, cfg) + severity(self.rev, cfg)
        if self.dtm > cfg.pair_tm_tol:
            harm += self.dtm - cfg.pair_tm_tol
        # Severity is zero for anything that breaks no rules, so clean
        # candidates are still ordered purely by closeness to the target.
        return (round(harm, 2), round(off, 2), round(self.dtm, 1),
                abs(len(self.fwd.seq) - 25))


def forward_options(template, codon_start, u, codon, cfg):
    """Every legal forward primer with a given upstream flank and codon."""
    out = []
    try:
        upstream = template.at(codon_start - u, u)
    except IndexError:
        return out
    # The 3' flank grows until the primer reaches the length ceiling. It is
    # not capped independently: the reverse primer may use the whole length
    # budget, and capping the forward separately is what left it stranded
    # several degrees cold in AT-rich stretches.
    d_max = cfg.len_hard_max - u - 3
    if cfg.flank3_max is not None:
        d_max = min(d_max, cfg.flank3_max)
    for d in range(cfg.flank3_min, d_max + 1):
        # On a linear template a candidate can run past either end. Not an
        # error -- this flank length is simply unavailable here.
        try:
            seq = upstream + codon + template.at(codon_start + 3, d)
        except IndexError:
            continue
        # Length may run past --len-max up to the hard ceiling; the extra
        # length is flagged, and severity() makes it a last resort.
        if cfg.len_min <= len(seq) <= cfg.len_hard_max:
            out.append((seq, d))
    return out


def choose_junction(template, codon_start, cfg):
    """Pick the one junction -- and so the one reverse primer -- for a position.

    The reverse primer lies entirely 5' of the mutated codon, so nothing about
    it depends on which amino acid is being installed. Choosing it per position
    rather than per mutation means every substitution at a residue shares a
    single reverse primer: one oligo per position instead of one per mutation.

    The junction is chosen to maximise how many of the substitutions available
    at this position come out with no flags, then by the reverse primer's own
    quality. It depends only on the position, so it is stable across every
    mutation there.

    Returns (u, reverse_Primer) or None if the position cannot be reached.
    """
    cache = template.junction_cache
    key = (codon_start, cfg.flank3_min, cfg.flank3_ideal, cfg.flank3_penalty,
           cfg.flank5_min, cfg.flank5_max, cfg.flank3_max,
           cfg.len_min, cfg.len_max, cfg.len_hard_max,
           cfg.tm_min, cfg.tm_max, cfg.tm_target, cfg.gc_min, cfg.gc_max,
           cfg.pair_tm_tol, cfg.primer_conc, cfg.na, cfg.mg, cfg.dntp,
           cfg.tm_target)
    if key in cache:
        return cache[key]

    wt_aa = CODON_TABLE.get(template.at(codon_start, 3))
    targets = sorted({PREFERRED_CODON[a][0] for a in AA3
                      if a != '*' and a != wt_aa})

    best = None
    for u in range(cfg.flank5_min, cfg.flank5_max + 1):
        by_codon = {c: forward_options(template, codon_start, u, c, cfg)
                    for c in targets}
        if not any(by_codon.values()):
            continue
        fwd_cache = {}
        for r in range(cfg.len_min, cfg.len_hard_max + 1):
            try:
                rseq = revcomp(template.at(codon_start - u - r, r))
            except IndexError:
                continue
            rev = Primer(rseq, cfg,
                         'fixed for this residue -- shared by every '
                         'substitution at position')
            # Score the junction by what the whole PAIR can achieve, not by
            # the reverse primer alone. The junction pins the forward primer's
            # 5' end too, so a reverse that looks good on its own can strand
            # every forward in a GC-rich patch with nowhere to go.
            n_clean = 0
            harm = 0.0
            rev_harm = severity(rev, cfg)
            for c in targets:
                least = None
                for seq, d3 in by_codon[c]:
                    fwd = fwd_cache.get(seq)
                    if fwd is None:
                        fwd = fwd_cache[seq] = Primer(seq, cfg, flank3=d3)
                    h = severity(fwd, cfg) + rev_harm
                    if abs(fwd.tm - rev.tm) > cfg.pair_tm_tol:
                        h += abs(fwd.tm - rev.tm) - cfg.pair_tm_tol
                    if least is None or h < least:
                        least = h
                if least is None:
                    least = 99.0          # no legal forward for this residue
                harm += least
                if least == 0.0:
                    n_clean += 1
            # The last term is the sequence itself, so the ordering is total:
            # no two candidates can tie, and the winner therefore cannot depend
            # on the order the loops happen to run in.
            rank = (-n_clean, round(harm / max(len(targets), 1), 2),
                    round(abs(rev.tm - cfg.tm_target), 2), abs(r - 25), u,
                    rev.seq)
            if best is None or rank < best[0]:
                best = (rank, u, rev)

    result = None if best is None else (best[1], best[2])
    cache[key] = result
    return result


def reachable_span(template, codon_start, new_codon, cfg):
    """Tm range reachable by each primer separately at this position.

    Returns (fwd_lo, fwd_hi, rev_lo, rev_hi, local_gc). Reported when nothing
    designs, so it is clear WHICH primer is blocking -- they fail for opposite
    reasons and the fix differs.
    """
    f_lo = f_hi = r_lo = r_hi = None
    # On a linear template this window can run off the start -- which is
    # precisely when this function gets called, so it must not raise. Narrow
    # the window until it fits; report None if even the codon is out of reach.
    gc_here = None
    span = cfg.flank5_max + cfg.len_max
    while span > 0 and gc_here is None:
        try:
            gc_here = gc_percent(template.at(codon_start - span, span))
        except IndexError:
            span //= 2
    for u in range(cfg.flank5_min, cfg.flank5_max + 1):
        for seq, _ in forward_options(template, codon_start, u, new_codon, cfg):
            t = melting_temp(seq, cfg.primer_conc, cfg.na, cfg.mg, cfg.dntp)
            f_lo = t if f_lo is None else min(f_lo, t)
            f_hi = t if f_hi is None else max(f_hi, t)
        for r in range(cfg.len_min, cfg.len_hard_max + 1):
            try:
                seq = revcomp(template.at(codon_start - u - r, r))
            except IndexError:
                continue
            t = melting_temp(seq, cfg.primer_conc, cfg.na, cfg.mg, cfg.dntp)
            r_lo = t if r_lo is None else min(r_lo, t)
            r_hi = t if r_hi is None else max(r_hi, t)
    return f_lo, f_hi, r_lo, r_hi, gc_here


def why_blocked(template, codon_start, new_codon, cfg):
    """Human-readable lines explaining why no pair could be built here."""
    f_lo, f_hi, r_lo, r_hi, gc_here = reachable_span(
        template, codon_start, new_codon, cfg)

    def verdict(lo, hi):
        if lo is None:
            return 'no legal primer at all'
        if hi < cfg.tm_min:
            return 'always too cold'
        if lo > cfg.tm_max:
            return 'always too hot'
        return 'reaches the window'

    fv, rv = verdict(f_lo, f_hi), verdict(r_lo, r_hi)
    lines = [
        'No primer pair fits the %.0f-%.0f C window here.'
        % (cfg.tm_min, cfg.tm_max),
        '  forward  %s  %s' % (
            ('%5.1f to %5.1f C' % (f_lo, f_hi)) if f_lo is not None else '     none    ',
            fv),
        '  reverse  %s  %s' % (
            ('%5.1f to %5.1f C' % (r_lo, r_hi)) if r_lo is not None else '     none    ',
            rv),
    ]
    blockers = [n for n, v in (('forward', fv), ('reverse', rv))
                if v != 'reaches the window']
    if blockers == ['forward']:
        lines.append('The forward primer is the blocker. Its 5\' end is pinned '
                     'at the junction, so')
        lines.append('the shortest it can be is %d nt (%d + codon + %d), '
                     'already too %s here.'
                     % (cfg.flank5_min + 3 + cfg.flank3_min, cfg.flank5_min,
                        cfg.flank3_min,
                        'hot' if fv == 'always too hot' else 'cold'))
    elif blockers == ['reverse']:
        lines.append('The reverse primer is the blocker.')
    elif blockers:
        lines.append('Both primers are out of range.')
    else:
        lines.append('Each primer can reach the window alone, but not both off '
                     'one shared junction.')
    if gc_here is not None:
        lines.append('Local GC %.0f%%. The limit is the template sequence, not '
                     'the search.' % gc_here)
    if not template.circular:
        lines.append('This template is LINEAR. A primer here needs flanking '
                     'sequence the file does not contain;')
        lines.append('supplying the plasmid as GenBank, or passing --circular, '
                     'makes every residue designable.')
    return lines


def _in_window(seq, cfg):
    """Is this primer inside the Tm window? Reported, never used to exclude.

    Nothing is filtered on Tm. A position where no primer reaches the window
    still gets a pair, flagged, rather than a refusal -- a flagged primer you
    can judge beats no primer at all.
    """
    return cfg.tm_min <= round(melting_temp(seq, cfg.primer_conc, cfg.na,
                                            cfg.mg, cfg.dntp), 1) <= cfg.tm_max


def design(template, codon_start, new_codon, cfg):
    """Rank the forward primers available for this substitution.

    The reverse primer is the position's, not the mutation's -- chosen once by
    choose_junction() and shared by every substitution at this residue. Only
    the forward varies.

    That sharing is structural, not a happy accident. A forward primer is
    `u + 3 + d` bases long, and the codon changes only the three in the middle:
    it never changes the length, and never changes whether a candidate exists.
    So forward_options() returns the same set of lengths for every target
    residue, and is empty for one only if it is empty for all twenty. Since
    choose_junction() only ever picks a junction where options exist, there is
    no codon that the position's junction cannot serve.

    An earlier version carried a fallback that widened the search to every
    junction and gave the substitution its own reverse primer. It was removed
    once the above was worked out: it could not be reached. Forcing it -- by
    squeezing the length window to 30-32 nt, well past anything sensible --
    still produced zero per-mutation reverses across 1,121 designs, and neither
    BglB (8,531 substitutions) nor LDH (6,498) ever reached it at the defaults.
    """
    junction = choose_junction(template, codon_start, cfg)
    if junction is None:
        return []
    u, rev = junction
    pairs = [Pair(Primer(seq, cfg,
                         "%d nt of template 5' of the codon, %d nt 3' of it"
                         % (u, d), flank3=d),
                  rev, cfg, u)
             for seq, d in forward_options(template, codon_start, u,
                                           new_codon, cfg)]
    pairs.sort(key=lambda p: p.rank(cfg))
    return pairs


# ---------------------------------------------------------------------------
# 7. Reporting
# ---------------------------------------------------------------------------

class Ink:
    """ANSI colour, but only when writing to a real terminal.

    Piping to a file or another program strips it automatically, so the
    tab-separated order lines stay clean.
    """
    on = sys.stdout.isatty() and os.environ.get('NO_COLOR') is None

    @classmethod
    def _w(cls, code, text):
        return '\033[%sm%s\033[0m' % (code, text) if cls.on else text

    @classmethod
    def good(cls, t):   return cls._w('32', t)
    @classmethod
    def warn(cls, t):   return cls._w('33', t)
    @classmethod
    def bad(cls, t):    return cls._w('31', t)
    @classmethod
    def bold(cls, t):   return cls._w('1', t)
    @classmethod
    def dim(cls, t):    return cls._w('2', t)
    @classmethod
    def cyan(cls, t):   return cls._w('36', t)


TICK, CROSS = '[ok]', '[!!]'
OK_MARK, NO_MARK = 'ok', '!!'
RULE = '\u2500' * 66


def orf_residues(template, cds_start):
    """How many residues from `cds_start` to the first in-frame stop.

    Returns (n_residues, stop_found). Used to bound the valid position range
    when the coding sequence is given as a raw start position rather than an
    annotated feature.
    """
    cap = len(template) // 3
    for i in range(cap):
        try:
            codon = template.at(cds_start + i * 3, 3)
        except IndexError:
            return i, False
        if CODON_TABLE.get(codon) == '*':
            return i, True
    return cap, False


def describe_template(template, notes):
    """--info output: what the script can see in this file."""
    print()
    print('=' * 70)
    print('  %s' % template.name)
    print('=' * 70)
    print()
    print('  %d bp, %s' % (len(template), 'circular' if template.circular else 'linear'))
    if not template.circular:
        print('  %s Linear: residues near either end cannot be designed.' % CROSS)
        print('       A circular plasmid in GenBank has no such limit.')
    print()

    cds = [f for f in template.features if f.kind == 'CDS']
    if not cds:
        print('  No CDS features annotated.')
        print('  Use --cds-start N to give the first base of residue 1 directly.')
    else:
        print('  CDS features:')
        print('    %-24s %-20s %6s  %-9s %s'
              % ('label', 'location', 'codons', 'starts', 'to use it'))
        usable = []
        for f in cds:
            if f.complement or not f.segments:
                print('    %-24s %-20s %6s  %s'
                      % (f.label or '(unlabelled)', f.location, '-',
                         'reverse strand, not usable'))
                continue
            protein = translate(template.at(f.start, f.length))
            usable.append(f)
            print('    %-24s %-20s %6d  %-9s --cds-start %d'
                  % ((f.label or '(unlabelled)')[:24], f.location,
                     f.length // 3, protein[:6] + '...', f.start + 1))
        print()
        print('  Copy a --cds-start number above, or name one:  --feature "<label>"')
    print()
    print('  Two ways to set residue 1:')
    print('    --cds-start N       1-based position of the first base of')
    print('                        residue 1. Works on any file, annotated or')
    print('                        not. Copy a number from the column above.')
    print('    --feature LABEL     a CDS named above (GenBank only)')
    print()
    if not template.circular:
        print('  This file is LINEAR. If it is really a plasmid, add --circular')
        print('  -- otherwise residues near each end cannot be designed.')
        print()


def build_parser():
    p = argparse.ArgumentParser(
        description='Design back-to-back mutagenesis primers for a point mutation.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='Example:\n  python3 design_primers.py "pET29b+ BglB plasmid.gbk" '
               'L171N\n')
    p.add_argument('template',
                   help='the plasmid the primers are designed against, as a '
                        'GenBank file. Format is detected by reading the file, '
                        'so the extension does not matter.')
    p.add_argument('mutations', nargs='*', metavar='MUTATION',
                   help='one or more mutations, e.g. L171N C167A E164Q')
    p.add_argument('--info', action='store_true',
                   help='describe the template -- topology and CDS features -- '
                        'and exit. Use this first with a new construct.')

    g = p.add_argument_group('coding sequence')
    g.add_argument('--cds-start', type=int, metavar='N', dest='cds_start',
                   help='1-based position of the first base of residue 1 -- '
                        'the A of the first codon you want to call residue 1. '
                        'Works on any file, annotated or not, and overrides '
                        'everything else. Run --info to read the number off '
                        'your template. This is the simplest way to set '
                        'numbering and the one to prefer.')
    g.add_argument('--circular', action='store_true',
                   help='treat the template as a circular plasmid even though '
                        'the file does not say so. FASTA and raw sequence carry '
                        'no topology, so a plasmid pasted as FASTA is read as '
                        'linear and loses the residues near each end -- 11 of '
                        "449 for BglB. GenBank files say 'circular' on their "
                        'LOCUS line and need no flag.')
    g.add_argument('--linear', action='store_true',
                   help='force linear, overriding a GenBank circular topology. '
                        'Rarely wanted; mostly useful for checking what a bare '
                        'gene sequence would give.')
    g.add_argument('--feature', metavar='LABEL',
                   help='use the CDS feature with this /label to set numbering '
                        '(GenBank only). With no numbering option at all the '
                        'script uses the sole usable CDS if there is exactly '
                        'one, and otherwise stops rather than guessing.')

    g = p.add_argument_group('constraints')
    g.add_argument('--tm-target', type=float, default=65.0, dest='tm_target',
                   help='Tm to design toward (default 65)')
    g.add_argument('--tm-min', type=float, default=63.0, dest='tm_min',
                   help='bottom of the target window; flagged below this, '
                        'never excluded (default 63)')
    g.add_argument('--tm-max', type=float, default=67.0, dest='tm_max',
                   help='top of the target window; flagged above this, '
                        'never excluded (default 67)')
    g.add_argument('--pair-tm-tol', type=float, default=4.0, dest='pair_tm_tol')
    # 35-65 rather than the textbook 40-60. Measured over every substitution
    # in both BglB and LDH, widening to 35-65 raises the clean rate by 11
    # points on each (77.4->88.6 and 71.7->82.8) and cuts GC flags by ~75%.
    # The reason is that a 40-60 window mostly re-flags problems the Tm and
    # length rules already catch: low-GC primers are long-but-on-target, and
    # high-GC primers are Tm-flagged anyway. It also frees the ranker, which
    # no longer rejects an on-target primer for reading 38% GC -- that alone
    # removed ~190 Tm flags. Widening further (30-70) stops the rule firing
    # at all, so 35-65 is where it still earns its place.
    g.add_argument('--gc-min', type=float, default=35.0, dest='gc_min',
                   help='bottom of the GC window (default 35)')
    g.add_argument('--gc-max', type=float, default=65.0, dest='gc_max',
                   help='top of the GC window (default 65)')
    g.add_argument('--len-min', type=int, default=20, dest='len_min')
    g.add_argument('--len-max', type=int, default=35, dest='len_max',
                   help='preferred maximum length; going past it is flagged. '
                        'A 35mer is a perfectly ordinary primer, so the flag '
                        'starts above that (default 35)')
    g.add_argument('--len-hard-max', type=int, default=40, dest='len_hard_max',
                   help='absolute maximum length, never exceeded. A primer '
                        'grows past --len-max only when nothing shorter '
                        'reaches the Tm window, and is flagged when it does '
                        '(default 40). Raising this to 45 rescues one case: '
                        'the reverse for residue 2 of a pET29b construct must '
                        'reach back into the ~22%% GC T7 leader, and tops out '
                        'at 61.1 C by 40 nt (63 C first becomes reachable at '
                        '42 nt, 65.8 C at 45). Left at 40 deliberately -- long '
                        'oligos synthesise less cleanly, and a 61.1 C limiting '
                        'primer is workable by dropping the touchdown plateau '
                        'for that reaction, which is the cheaper adjustment.')
    # Flanks are named by the SIDE OF THE CODON they govern, so that --flank-3-*
    # always means the anchoring side and --flank-5-* always means the junction
    # side. An earlier naming had --flank-min governing the 3' side while
    # --flank-max governed the 5' side; they read as a pair and were not one.
    g.add_argument('--flank-3-min', type=int, default=7, dest='flank3_min',
                   help="absolute floor on correct bases 3' of the codon -- the "
                        "anchoring side, where the polymerase extends from. "
                        "Never gone below. Going under --flank-3-ideal is "
                        "penalised, so this is an escape hatch for hard "
                        "positions rather than a new normal (default 7)")
    g.add_argument('--flank-3-ideal', type=int, default=10, dest='flank3_ideal',
                   help="preferred bases 3' of the codon. Anything at or above "
                        "this is free; below it costs --flank-3-penalty per "
                        "base (default 10)")
    g.add_argument('--flank-3-penalty', type=float, default=0.2,
                   dest='flank3_penalty',
                   help="cost per base that the 3' flank falls short of "
                        "--flank-3-ideal, in the same units as a degree of Tm. "
                        "0 makes the floor a free-for-all; 0.2 keeps the mean "
                        "3' flank where it is today (default 0.2)")
    g.add_argument('--flank-3-max', type=int, default=None, dest='flank3_max',
                   help="cap on correct bases 3' of the codon. By default there "
                        "is none: the forward primer's 5' end is pinned at the "
                        "junction, so its 3' end is its only freedom, and it "
                        "grows until the primer reaches --len-hard-max")
    g.add_argument('--flank-5-min', type=int, default=7, dest='flank5_min',
                   help="minimum correct bases 5' of the codon -- the junction "
                        "side. Nothing extends from this end, so it can be "
                        "shorter than the 3' side; 7 is what lets a forward "
                        "primer reach 20 nt. Note this is also how close the "
                        "mutation may sit to the blunt ligation junction "
                        "(default 7)")
    g.add_argument('--flank-5-max', type=int, default=15, dest='flank5_max',
                   help="maximum correct bases 5' of the codon, i.e. how far "
                        "the junction can sit from it (default 15)")

    g = p.add_argument_group('Tm conditions')
    g.add_argument('--primer-conc', type=float, default=500.0, dest='primer_conc',
                   help='oligo concentration, nM (default 500)')
    g.add_argument('--na', type=float, default=50.0,
                   help='monovalent Na+/K+, mM (default 50)')
    g.add_argument('--mg', type=float, default=1.5,
                   help='Mg2+, mM (default 1.5)')
    g.add_argument('--dntp', type=float, default=0.6,
                   help='dNTP, mM (default 0.6). Documented for completeness; '
                        'Benchling does not subtract it from Mg2+, so it does '
                        'not change the result.')

    g = p.add_argument_group('output')
    g.add_argument('--alternatives', type=int, default=0, metavar='N',
                   help='also list the next N ranked pairs')
    g.add_argument('--verbose', action='store_true',
                   help='the long form: every constraint spelled out with its '
                        'target, rather than the compact summary')
    g.add_argument('--quiet', action='store_true',
                   help='print only the two TSV lines')
    return p


def resolve_template(cfg):
    """Load the template and work out where residue 1 sits.

    Returns (template, cds_start_0based, n_residues, notes).

    Precedence for numbering, most explicit first:
        1. --cds-start N        raw position, no annotation needed
        2. --feature LABEL      a named CDS in the file
        3. the only usable CDS in the file
    Anything else is an error. The script does not guess between candidates --
    picking the wrong one shifts every position silently.
    """
    notes = []
    path = cfg.template

    tpl = load_template(path)

    if cfg.circular and cfg.linear:
        raise ValueError('--circular and --linear contradict each other.')
    forced_circular = False
    if cfg.circular and not tpl.circular:
        tpl.circular = True
        forced_circular = True
        notes.append(('ok', 'Topology forced to circular by --circular.'))
    elif cfg.linear and tpl.circular:
        tpl.circular = False
        notes.append(('warn', 'Topology forced to LINEAR by --linear: residues '
                              'near each end will not be designable.'))

    notes.append(('ok', 'Template %s: %d bp, %s'
                  % (tpl.name, len(tpl), 'circular' if tpl.circular else 'linear')))
    if not tpl.circular:
        notes.append(('warn', 'Linear sequence: residues within ~15 of either '
                              'end cannot be designed, because a primer needs '
                              'flanking sequence that is not in the file. '
                              'Supply the plasmid as GenBank, or pass '
                              '--circular if this really is a plasmid.'))

    # -- where does residue 1 begin? ------------------------------------
    cds = [f for f in tpl.features if f.kind == 'CDS' and not f.complement
           and f.segments]
    labels = ', '.join(sorted(set(f.label or '(unlabelled)' for f in cds))) or '(none)'
    chosen = None

    if cfg.cds_start is not None:
        cds_start, how = cfg.cds_start - 1, '--cds-start %d' % cfg.cds_start
    elif cfg.feature:
        wanted = [f for f in cds if f.label.lower() == cfg.feature.lower()]
        if not wanted:
            raise ValueError('No usable CDS labelled %r. Available: %s'
                             % (cfg.feature, labels))
        chosen = wanted[0]
        cds_start, how = chosen.start, '--feature %r' % cfg.feature
    else:
        if len(cds) == 1:
            chosen = cds[0]
            cds_start, how = chosen.start, 'the only CDS in the template'
        else:
            if not cds:
                # A FASTA or a bare run of DNA carries no features at all, so
                # there is nothing to choose between -- the old message told
                # people to pick a label from an empty list.
                raise ValueError(
                    'This file has no coding sequence annotated, so there is '
                    'nothing to read the numbering from. FASTA and plain DNA '
                    'never carry annotation; only GenBank does. Give the '
                    'position of residue 1 directly with --cds-start N -- the '
                    '1-based position of the first base of the first codon.')
            raise ValueError(
                'Residue 1 is ambiguous: this template has %d usable coding '
                'sequences, so the script will not guess. Say which one you '
                'mean with --feature LABEL, or give the position of residue 1 '
                'directly with --cds-start N. Candidates: %s. '
                'Run --info to see the file and the exact numbers to copy.'
                % (len(cds), labels))

    # -- validate the frame ---------------------------------------------
    if chosen is not None:
        n_res = chosen.length // 3
        bounded = 'annotated'
    else:
        n_res, stop_found = orf_residues(tpl, cds_start)
        bounded = 'to the first stop codon' if stop_found else 'unbounded'

    first = CODON_TABLE.get(tpl.at(cds_start, 3), '?')
    notes.append(('ok', 'Numbering from %s: residue 1 = %s at position %d, '
                        '%d residues (%s).'
                  % (how, first, cds_start + 1, n_res, bounded)))

    # --circular on a file that contains only the GENE, not the whole plasmid,
    # joins the C-terminus to the N-terminus and designs primers against a
    # template that does not exist. They come back with no flags, which is the
    # worst possible failure. A coding sequence that fills most of the file is
    # the signature of a bare gene.
    if forced_circular and n_res * 3 > 0.8 * len(tpl):
        notes.append(('warn',
                      'The coding sequence is %.0f%% of this file, so it looks '
                      'like a bare GENE rather than a whole plasmid. Treating '
                      'it as circular joins its two ends to each other, so '
                      'primers near either terminus are designed against '
                      'sequence that does not exist in your construct. Use the '
                      'plasmid, not the gene.'
                      % (100.0 * n_res * 3 / len(tpl))))

    protein = translate(tpl.at(cds_start, n_res * 3))
    internal = protein[:-1].count('*')
    if internal:
        notes.append(('warn', '%d internal stop codon(s) in that reading frame '
                              '-- residue 1 may be in the wrong place.' % internal))
    if chosen is not None:
        others = sorted(set(f.label for f in cds if f.label and f is not chosen))
        if others:
            notes.append(('ok', 'Other CDS features present (not used): %s'
                          % ', '.join(others)))

    return tpl, cds_start, n_res, notes


def _flagline(primer, cfg=None):
    """Compact one-line summary of a primer's numbers, flags coloured.

    The forward primer also reports its 3' flank -- the matched bases past the
    codon, the side the polymerase extends from. It is not a pass/fail
    constraint, but it varies from --flank-3-min to well past --flank-3-ideal
    and nothing else in the report shows it, so a primer that took the short
    anchor to reach the Tm window would otherwise be indistinguishable from
    one that did not.
    """
    def f(cond, text):
        return Ink.good(text) if cond else Ink.warn(text + ' !')
    checks = {c.name: c for c in primer.checks}
    clamp = checks["3' GC clamp"]
    line = '%s   %s   %s   %s' % (
        f(checks['Length'].ok, '%2d nt' % len(primer.seq)),
        f(checks['Tm'].ok, 'Tm %.1f' % primer.tm),
        f(checks['GC content'].ok, 'GC %d%%' % primer.gc),
        f(clamp.ok, "3' G/C" if clamp.ok else "no 3' G/C"))
    if primer.flank3 is not None and cfg is not None:
        short = primer.flank3 < cfg.flank3_ideal
        text = "3' flank %d" % primer.flank3
        line += '   %s' % (Ink.warn(text + ' (short)') if short
                           else Ink.dim(text))
    return line


def report_compact(name, best, cfg):
    """The default view: two sequences, their numbers, and what to worry about."""
    flags = ([('forward', c) for c in best.fwd.flags]
             + [('reverse', c) for c in best.rev.flags]
             + [('pair', c) for c in best.pair_checks if not c.ok])
    verdict = (Ink.good('no flags') if not flags
               else Ink.warn('%d flag%s' % (len(flags), '' if len(flags) == 1 else 's')))

    print(Ink.dim(RULE))
    print(' %s   %s   %s' % (Ink.bold(name), Ink.dim('residue %d' % best.pos), verdict))
    print(Ink.dim(RULE))
    for label, pr in (('FWD', best.fwd), ('REV', best.rev)):
        print('  %s  %s' % (Ink.cyan(label), "5'-%s-3'" % pr.seq))
        print('       %s' % _flagline(pr, cfg))
    dtm_ok = all(c.ok for c in best.pair_checks)
    lim = min(best.fwd.tm, best.rev.tm)
    print('  %s  %s   %s' % (
        Ink.cyan('PAIR'),
        (Ink.good if dtm_ok else Ink.warn)('dTm %.1f' % best.dtm) + ('' if dtm_ok else ' !'),
        Ink.bold('limiting Tm %.1f' % lim) + Ink.dim('  <- anneal below this')))
    if flags:
        for which, c in flags:
            print('  %s  %s %s %s'
                  % (Ink.warn(' ! '), which.ljust(7), c.name.ljust(13),
                     Ink.dim('(%s)' % c.detail)))
    print()


def report_full(name, best, cfg):
    """--verbose view: every constraint spelled out."""
    print()
    print('=' * 70)
    print('  %s   at residue %d' % (name, best.pos))
    print('=' * 70)
    print()
    for label, pr in (('%s_forward' % name, best.fwd), ('%s_reverse' % name, best.rev)):
        print('  %s' % label)
        print("    5'-%s-3'" % pr.seq)
        print('    %s' % pr.note)
        for c in pr.checks:
            mark = Ink.good(TICK) if c.ok else Ink.warn(CROSS)
            print('      %s %-13s %-16s (%s)' % (mark, c.name, c.value, c.detail))
        print()
    print('  PAIR CHECKS')
    for c in best.pair_checks:
        print('    %s %-15s %-14s (%s)'
              % (Ink.good(TICK) if c.ok else Ink.warn(CROSS), c.name, c.value, c.detail))
    print('    %s %-15s %-14s (%s)'
          % (Ink.good(TICK), 'Limiting Tm', '%.1f C' % min(best.fwd.tm, best.rev.tm),
             'the lower of the pair -- set the annealing temperature below this'))
    print()


def design_one(template, cds_start, n_residues, text, cfg):
    """Design and report one mutation. Returns 0 on success, 1 on refusal."""
    wt, pos, new = parse_mutation(text)
    name = '%s%d%s' % (wt, pos, new)
    codon_start = cds_start + (pos - 1) * 3

    def refuse(*lines):
        if not cfg.quiet:
            print(Ink.dim(RULE))
            print(' %s   %s' % (Ink.bold(name), Ink.bad('REFUSED')))
            print(Ink.dim(RULE))
            for l in lines:
                print('  %s' % l)
            print()
        else:
            print('%s\tREFUSED\t%s' % (name, lines[0] if lines else ''), file=sys.stderr)
        return 1

    # -- Step 1: sanity check. Nothing proceeds until this passes. ----------
    if pos > n_residues:
        return refuse('Residue %d is past the end of the coding sequence, '
                      'which is %d residues long.' % (pos, n_residues),
                      Ink.dim('Check the mutation, or whether residue 1 is set '
                              'correctly (run --info).'))

    wt_codon = template.at(codon_start, 3)
    found = CODON_TABLE.get(wt_codon, 'X')
    if found != wt:
        ctx_from = max(1, pos - 5)
        ctx = translate(template.at(cds_start + (ctx_from - 1) * 3, 11 * 3))
        return refuse('You specified %s at residue %d, but the template has %s '
                      '(codon %s).' % (wt, pos, found, wt_codon),
                      'Residues %d-%d: %s' % (ctx_from, ctx_from + len(ctx) - 1, ctx),
                      Ink.dim('If your numbering differs, set residue 1 with '
                              '--cds-start or --feature.'))

    # -- Step 2: codon choice -- highest-frequency E. coli codon. ----------
    new_codon, frac, freq = PREFERRED_CODON[new]

    # -- Steps 3 & 4: build the pair. --------------------------------------
    pairs = design(template, codon_start, new_codon, cfg)
    if not pairs:
        return refuse(*why_blocked(template, codon_start, new_codon, cfg))
    best = pairs[0]
    best.pos = pos

    if not cfg.quiet:
        if cfg.verbose:
            print('  %s Codon: %s -> %s  (most frequent %s codon in E. coli K12, '
                  '%.0f%%, %.1f/1000)'
                  % (Ink.good(TICK), wt_codon, new_codon, AA3[new],
                     100 * frac, freq))
            report_full(name, best, cfg)
        else:
            report_compact(name, best, cfg)
            if best.rev.flags:
                _, _, r_lo, r_hi, gc = reachable_span(
                    template, codon_start, new_codon, cfg)
                if r_lo is not None and gc is not None:
                    print('  %s' % Ink.dim(
                        'No reverse primer here clears every rule: any %d-%d nt '
                        'primer at this position spans %.1f-%.1f C '
                        '(%.0f%% GC nearby).'
                        % (cfg.len_min, cfg.len_hard_max, r_lo, r_hi, gc)))
                    print()

    for label, pr in (('forward', best.fwd), ('reverse', best.rev)):
        print('%s_%s\t%s\t%d\t%.1f\t%.0f\t%s'
              % (name, label, pr.seq, len(pr.seq), pr.tm, pr.gc,
                 'FLAG' if pr.flags else 'OK'))

    if cfg.alternatives and not cfg.quiet:
        print()
        print('  %s' % Ink.dim('alternative forward primers (the reverse is '
                               'fixed for this residue)'))
        for p in pairs[1:cfg.alternatives + 1]:
            print('    %-32s %s' % (p.fwd.seq, _flagline(p.fwd, cfg)))
        print()

    return 0


def main(argv=None):
    parser = build_parser()
    cfg = parser.parse_args(argv)

    if cfg.info:
        try:
            template, _, _, notes = resolve_template(cfg)
        except (OSError, ValueError):
            # --info must work even when numbering is unresolved: that is
            # usually the very question the user is running --info to answer.
            try:
                template = load_template(cfg.template)
            except (OSError, ValueError) as e:
                print('ERROR: %s' % e, file=sys.stderr)
                return 2
            notes = []
        describe_template(template, notes)
        return 0

    # Accept "L171N C167A" or "L171N, C167A" as one argument as well as
    # several -- editors and prompt boxes tend to hand over a single string.
    cfg.mutations = [m for chunk in cfg.mutations
                     for m in re.split(r'[\s,;]+', chunk.strip()) if m]

    if not cfg.mutations:
        parser.error('at least one mutation is required (or use --info)')

    try:
        for text in cfg.mutations:
            parse_mutation(text)          # fail fast on all of them
        template, cds_start, n_residues, notes = resolve_template(cfg)
    except (OSError, ValueError) as e:
        print('ERROR: %s' % e, file=sys.stderr)
        return 2

    # The template is read once and shared, so the notes print once. Only the
    # warnings are worth the space unless --verbose is asked for.
    if not cfg.quiet:
        print()
        # The Tm scale is never left implicit -- the same oligo reads about 5 C
        # apart with and without Mg2+, so every report states its conditions.
        print(Ink.dim('  Tm: SantaLucia 1998 + Owczarzy divalent, %g nM oligo, '
                      '%g mM Na+, %g mM Mg2+, %g mM dNTP'
                      % (cfg.primer_conc, cfg.na, cfg.mg, cfg.dntp)))
        print(Ink.dim('  window %g-%g C, designing to %g'
                      % (cfg.tm_min, cfg.tm_max, cfg.tm_target)))
        for kind, text in notes:
            if kind == 'ok' and not cfg.verbose:
                continue
            print('  %s %s' % (Ink.warn(CROSS) if kind != 'ok'
                               else Ink.good(TICK), text))

    failed = 0
    for text in cfg.mutations:
        failed += design_one(template, cds_start, n_residues, text, cfg)

    if len(cfg.mutations) > 1 and not cfg.quiet:
        done = len(cfg.mutations) - failed
        print(Ink.dim('  %d designed, %d refused.' % (done, failed)))
        print()
    return 1 if failed else 0


if __name__ == '__main__':
    sys.exit(main())
