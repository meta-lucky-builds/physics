"""Physical constants in SI units, from scipy's CODATA tables where available."""
from scipy import constants as _c

c = _c.c                  # speed of light, m/s
G = _c.G                  # gravitational constant, m^3 kg^-1 s^-2
h = _c.h                  # Planck constant, J s
hbar = _c.hbar            # reduced Planck constant, J s
k_B = _c.k                # Boltzmann constant, J/K
sigma_SB = _c.sigma       # Stefan-Boltzmann constant, W m^-2 K^-4
e = _c.e                  # elementary charge, C
m_e = _c.m_e              # electron mass, kg
m_p = _c.m_p              # proton mass, kg
N_A = _c.N_A              # Avogadro constant, 1/mol
R = _c.R                  # molar gas constant, J mol^-1 K^-1
g_0 = _c.g                # standard gravity, m/s^2
eV = _c.eV                # electronvolt, J

# Astronomical values (IAU 2015 nominal values and J2000 conventions).
M_sun = 1.98841e30        # solar mass, kg
M_earth = 5.9722e24       # Earth mass, kg
R_earth = 6.3781e6        # Earth equatorial radius, m
year = 365.25 * 86400     # Julian year, s
ly = c * year             # light-year, m
