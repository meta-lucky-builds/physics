"""C42-C46: what '2% different' does and doesn't say; the Bayes structure of an undetectable hypothesis."""
import numpy as np
from _harness import Checks

C = Checks("Sequence percent and undetectable manipulation")

C.section("Genome arithmetic (Chimpanzee Sequencing and Analysis Consortium 2005)")
aligned = 2.4e9
snv = 0.0123
print(f"  aligned sequence {aligned:.1e} bp, substitutions at 1.23%: {aligned*snv:.1e} (paper: about 35 million)")
C.check("1.23% of the aligned sequence reproduces about 30-35 million", 2.5e7 < aligned * snv < 3.6e7)
indel_bp = 90e6
print(f"  indel sequence total about 90 Mb in about 5 million events = {indel_bp/5e6:.0f} bp per event on average")
print(f"  paper's total: 1.23% by substitutions, about 3% more by indels. Events: {35/5:.0f} times more substitutions, "
      f"bases: {indel_bp/(aligned*snv):.1f} times more in indels")
C.check("event count and base count rank the two kinds of change oppositely: the percent depends on the unit",
        35e6 > 5e6 and indel_bp > aligned * snv)
print("  Yoo et al. 2025 complete genomes: 12.5-27.3% of an ape genome fails simple one-to-one alignment (all apes)")
print("  Same pair, three legitimate numbers: about 1.2%, about 4%, above 12%. Which is 'the' difference? None is.")

C.section("Same percent, opposite function: real genetic code")
bases = "TCAG"
aa = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
code = {a + b + c_: aa[16 * i + 4 * j + k] for i, a in enumerate(bases) for j, b in enumerate(bases)
        for k, c_ in enumerate(bases)}


def translate(s):
    out = []
    for i in range(0, len(s) - 2, 3):
        p = code[s[i:i + 3]]
        out.append(p)
        if p == "*":
            break
    return "".join(out)


rng = np.random.default_rng(2)
non_stop = [k for k, v in code.items() if v != "*"]
n_cod = 500
seq = "".join(rng.choice(non_stop, size=n_cod)) + "TAA"
L = len(seq)
target = int(round(0.02 * L))


def protein_identity(p, q):
    n = max(len(p), len(q))
    return sum(a == b for a, b in zip(p, q)) / n


prot = translate(seq)
# (a) synonymous-only edits at 2% of bases
s_syn = list(seq)
changed = 0
for pos in rng.permutation(n_cod):
    cod = seq[3 * pos:3 * pos + 3]
    alts = [k for k, v in code.items() if v == code[cod] and k != cod]
    alt1 = [k for k in alts if sum(x != y for x, y in zip(k, cod)) == 1]
    if alt1 and changed < target:
        s_syn[3 * pos:3 * pos + 3] = list(alt1[0])
        changed += 1
s_syn = "".join(s_syn)
# (b) same number of changed bases, spread at random
s_rand = list(seq)
idx = rng.choice(L - 3, size=target, replace=False)
for i in idx:
    s_rand[i] = rng.choice([b for b in bases if b != s_rand[i]])
s_rand = "".join(s_rand)
# (c) one single-base deletion near the start, nearly identical as a string
s_fs = seq[:30] + seq[31:]
# (d) one nonsense change at the start, only one letter different
for pos in range(3, 60, 3):
    cod = seq[pos:pos + 3]
    stops = [k for k, v in code.items() if v == "*" and sum(x != y for x, y in zip(k, cod)) == 1]
    if stops:
        s_stop = seq[:pos] + stops[0] + seq[pos + 3:]
        break
for name, s in [("synonymous only, 2% of bases", s_syn), ("random, 2% of bases", s_rand),
                ("1-base frameshift, 0.07% of bases", s_fs), ("1-base stop gain, 0.07% of bases", s_stop)]:
    diffs = sum(a != b for a, b in zip(seq, s)) if len(s) == L else 1
    print(f"  {name:36s}: base differences {diffs:3d}/{L}, protein identity {100*protein_identity(prot, translate(s)):6.2f}%")
C.check("2% synonymous change leaves the protein 100% identical", protein_identity(prot, translate(s_syn)) == 1.0)
C.check("one base (0.07%) can destroy over 90% of a protein", protein_identity(prot, translate(s_fs)) < 0.2 or
        protein_identity(prot, translate(s_stop)) < 0.1)

C.section("Brains are not 2% apart either")
print(f"  neurons: human 86 billion, chimpanzee predicted 28 billion (Herculano-Houzel 2012): ratio {86/28:.1f}")
print("  So 1.23% sequence divergence can coexist with a threefold neuron count difference.")
C.check("the neuron ratio is not 1.02", abs(86 / 28 - 1.02) > 1)
print("  Scale test: IQ is an interval scale built by norming, not a ratio. 2% above something depends on the zero:")
for shift in [0, 100, -50]:
    print(f"    add {shift:4d} to every score: 100 -> {100+shift}, 130 -> {130+shift}, ratio {(130+shift)/(100+shift):.3f}")
C.check("ratios on a shiftable scale are not meaningful", (130 / 100) != (230 / 200))

C.section("Bayes: hypotheses and their predictions")


def update(prior_h, like_h, like_alt):
    return prior_h * like_h / (prior_h * like_h + (1 - prior_h) * like_alt)


# H0: ordinary history, with data D drawn from some model.  H_mimic: secret manipulators who
# exactly reproduce every observation H0 predicts. Then P(D|H_mimic) = P(D|H0) for all D.
rng = np.random.default_rng(10)
K = 1000
p0 = rng.dirichlet(np.ones(K))
p_mimic = p0.copy()
post = 0.01
for d in rng.choice(K, size=10000, p=p0):
    post = update(post, p_mimic[d], p0[d])
print(f"  H_mimic: prior 0.0100; after 10000 observations posterior {post:.4f}")
C.check("the mimicking hypothesis can't be moved in either direction", abs(post - 0.01) < 1e-9)
# H_any: 'they could produce anything', spreads likelihood uniformly over the outcomes.
p_any = np.ones(K) / K
post_any = 0.01
logbf = 0.0
for d in rng.choice(K, size=10000, p=p0):
    logbf += np.log(p_any[d]) - np.log(p0[d])
print(f"  H_any: spreads over all outcomes; log10 Bayes factor after 10000 observations = {logbf/np.log(10):.1f}")
C.check("the vague 'could do anything' version loses to a sharper model", logbf < 0)
# H_trace: manipulation leaves a detectable anomaly with probability q per observation.
q = 0.001
n = np.array([100, 1000, 10000, 100000])
bf = (1 - q)**n
print("  H_trace: manipulators leave an anomaly with chance q = 0.001 per independent observation")
for ni, b in zip(n, bf):
    print(f"    {ni:7d} clean observations: Bayes factor {b:.2e}; posterior from 0.01 prior {update(0.01, b, 1):.2e}")
C.check("any version with detectable traces gets testable: clean data lowers it", bf[-1] < 1e-20)
C.finish()
