"""
Generate a mass-concentration lookup table at z=0 using the Diemer (2019)
model from the colossus package.

Output: concentration_table.bin
  - A simple binary file with a 4-byte header (uint32 number of points N),
    followed by N float64 log10(M [Msun/h]) values, then N float64 c values.

Dependencies:
    pip install colossus numpy
"""

import struct
import numpy as np
from colossus.cosmology import cosmology
from colossus.halo import concentration

# ── Cosmology ────────────────────────────────────────────────────────────────
# Planck 2018 (Planck18); swap to any colossus-supported cosmology as needed.
COSMOLOGY = "planck18"
cosmology.setCosmology(COSMOLOGY)

# ── Sampling ─────────────────────────────────────────────────────────────────
N_POINTS   = 500
LOG_M_MIN  = 8.0    # log10(M / [Msun/h])
LOG_M_MAX  = 12.0
REDSHIFT   = 0.0
MODEL      = "diemer19"
MDEF       = "200c"   # mass definition expected by Diemer+19

log_masses = np.linspace(LOG_M_MIN, LOG_M_MAX, N_POINTS)
masses     = 10.0 ** log_masses   # Msun/h

# ── Evaluate ─────────────────────────────────────────────────────────────────
# concentration.concentration returns a scalar or array depending on input.
concentrations = concentration.concentration(masses, MDEF, REDSHIFT, model=MODEL)

print(f"Model   : {MODEL}")
print(f"Cosmology: {COSMOLOGY}")
print(f"Redshift: {REDSHIFT}")
print(f"Mass def: {MDEF}")
print(f"N points: {N_POINTS}")
print(f"log10(M) range : [{LOG_M_MIN}, {LOG_M_MAX}]")
print(f"c end poitns   : [{concentrations[0]:.3f}, {concentrations[-1]:.3f}]")
print(f"c range        : [{concentrations.min():.3f}, {concentrations.max():.3f}]")

# ── Write binary ─────────────────────────────────────────────────────────────
# Layout:
#   [0..3]   : uint32  N          (number of sample points)
#   [4..4+8N): float64 log10_mass (N values, ascending)
#   [4+8N..) : float64 conc       (N values, corresponding concentrations)
OUT_FILE = "concentration_table.bin"

with open(OUT_FILE, "wb") as f:
    f.write(struct.pack("<I", N_POINTS))          # little-endian uint32
    f.write(log_masses.astype("<f8").tobytes())   # little-endian float64
    f.write(concentrations.astype("<f8").tobytes())

print(f"\nWrote {OUT_FILE}  ({4 + 2 * N_POINTS * 8} bytes)")
