"""L08, L09: uncertainty relation.

- sympy: [x, p] = i hbar on an arbitrary function.
- Fourier: sigma_x sigma_p computed by FFT for many wave packets, including random ones.
- Robertson's bound for random Hermitian operators and states, in finite dimensions.
- Sampling: separate copies measured for x only or p only. Nothing is disturbed, and
  the spreads still obey the bound.
"""
import numpy as np
import sympy as sp

from checklib import check, finish, section

rng = np.random.default_rng(5)

section("Commutator")
x, hb = sp.symbols("x hbar", real=True, positive=True)
f = sp.Function("f")(x)
p = lambda g: -sp.I * hb * sp.diff(g, x)
comm = sp.simplify(x * p(f) - p(x * f))
check("[x, p] f = i hbar f", sp.simplify(comm - sp.I * hb * f) == 0, f"result {comm}")

section("Fourier spreads (hbar = 1)")
N = 2**13
L = 200.0
dx = L / N
xs = (np.arange(N) - N / 2) * dx
ks = 2 * np.pi * np.fft.fftfreq(N, dx)


def spreads(psi):
    """sigma_x and sigma_p from a sampled wavefunction, momentum via FFT."""
    psi = psi / np.sqrt(np.sum(np.abs(psi) ** 2) * dx)
    Px = np.abs(psi) ** 2 * dx
    mx = np.sum(xs * Px)
    sx = np.sqrt(np.sum((xs - mx) ** 2 * Px))
    phi = np.fft.fft(psi)
    Pk = np.abs(phi) ** 2
    Pk = Pk / Pk.sum()
    mk = np.sum(ks * Pk)
    sk = np.sqrt(np.sum((ks - mk) ** 2 * Pk))
    return sx, sk, psi, Pk


for s in [0.3, 1.0, 4.0]:
    sx, sk, *_ = spreads(np.exp(-(xs**2) / (4 * s**2)))
    check(f"Gaussian sigma = {s}: product = 1/2", abs(sx * sk - 0.5) < 1e-8, f"sx sk = {sx * sk:.10f}")
sx, sk, *_ = spreads(np.exp(-(xs**2) / 4 + 1j * 0.5 * xs**2 + 3j * xs))
print(f"  Chirped Gaussian: sx sk = {sx * sk:.4f} (same |psi|, but phase structure widens momentum)")
check("chirp raises the product above 1/2", sx * sk > 0.6)
a = 2.0
psi = np.where(np.abs(xs) < a, np.cos(np.pi * xs / (2 * a)), 0.0)
sx, sk, *_ = spreads(psi)
exact = 0.5 * np.sqrt(np.pi**2 / 3 - 2)
check("infinite well ground state matches sqrt(pi^2/3 - 2)/2", abs(sx * sk / exact - 1) < 2e-3, f"FFT {sx * sk:.5f}, exact {exact:.5f}")
sx, sk, *_ = spreads(np.exp(-((xs - 5) ** 2) / 2) + np.exp(-((xs + 5) ** 2) / 2))
print(f"  Two separated lumps: sx sk = {sx * sk:.3f}")
sx, sk, *_ = spreads(1 / np.cosh(xs))
print(f"  sech packet: sx sk = {sx * sk:.4f}")
check("sech packet respects the bound", sx * sk >= 0.5)

low = np.inf
for _ in range(1000):
    n = rng.integers(1, 6)
    psi = np.zeros(N, complex)
    for _j in range(n):
        c = rng.normal() + 1j * rng.normal()
        psi += c * np.exp(-((xs - rng.uniform(-6, 6)) ** 2) / (4 * rng.uniform(0.3, 2) ** 2) + 1j * rng.uniform(-3, 3) * xs)
    sx, sk, *_ = spreads(psi)
    low = min(low, sx * sk)
check("1000 random packets: smallest product is at least 1/2", low >= 0.5 - 1e-9, f"smallest {low:.8f}")

section("Robertson bound with random operators and states")


def rand_herm(d):
    m = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d))
    return (m + m.conj().T) / 2


worst = np.inf
tight = np.inf
n_trials = 5000
for _ in range(n_trials):
    d = rng.integers(2, 9)
    A, B = rand_herm(d), rand_herm(d)
    v = rng.normal(size=d) + 1j * rng.normal(size=d)
    v /= np.linalg.norm(v)
    a = (v.conj() @ A @ v).real
    b = (v.conj() @ B @ v).real
    sa = np.sqrt((v.conj() @ A @ A @ v).real - a**2)
    sb = np.sqrt((v.conj() @ B @ B @ v).real - b**2)
    rhs = 0.5 * abs(v.conj() @ (A @ B - B @ A) @ v)
    worst = min(worst, sa * sb - rhs)
    # Stronger Schroedinger form: add the anticommutator term
    anti = 0.5 * (v.conj() @ (A @ B + B @ A) @ v).real - a * b
    tight = min(tight, (sa * sb) ** 2 - rhs**2 - anti**2)
check(f"sigma_A sigma_B >= |<[A,B]>|/2 in {n_trials} random trials", worst >= -1e-12, f"smallest margin {worst:.2e}")
check("Schroedinger's stronger form also holds", tight >= -1e-10, f"smallest margin {tight:.2e}")
# Commuting operators: no bound at all.
D = np.diag([1.0, 2.0, 3.0])
D2 = np.diag([5.0, -1.0, 0.5])
v = np.array([1, 0, 0], complex)
sa = np.sqrt((v @ D @ D @ v).real - (v @ D @ v).real ** 2)
sb = np.sqrt((v @ D2 @ D2 @ v).real - (v @ D2 @ v).real ** 2)
check("commuting operators allow both spreads to be zero", sa == 0 and sb == 0)

section("State spread with no disturbance: separate copies")
sg = 0.7
psi = np.exp(-(xs**2) / (4 * sg**2)) * (1 + 0.6 * np.sin(1.3 * xs))
sx0, sk0, psi_n, Pk = spreads(psi)
Px = np.abs(psi_n) ** 2 * dx
Px /= Px.sum()
n = 200000
sample_x = rng.choice(xs, size=n, p=Px)  # copies measured for position only
sample_k = rng.choice(ks, size=n, p=Pk)  # different copies measured for momentum only
prod = sample_x.std() * sample_k.std()
check("independent ensembles reproduce the state's spread product", abs(prod / (sx0 * sk0) - 1) < 0.01, f"sampled {prod:.4f}, exact {sx0 * sk0:.4f}")
check("that product is above 1/2 though no copy was measured twice", prod > 0.5)

finish()
