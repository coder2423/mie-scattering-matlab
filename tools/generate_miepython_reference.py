"""Regenerate pinned cross-language fixtures; MATLAB runtime needs no Python."""
import json
import os
import argparse
from pathlib import Path
import subprocess
import sys

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--miepython-source", type=Path, required=True,
                    help="External miepython Git checkout at the pinned commit.")
args = parser.parse_args()
REPO = args.miepython_source.resolve()
os.environ["MIEPYTHON_USE_JIT"] = "0"
sys.path.insert(0, str(REPO))
import miepython as mie  # noqa: E402

PINNED_COMMIT = "f4101f71aa6fbbf8a8e835b45778cee1b5cb637b"
commit = subprocess.check_output(["git", "-C", str(REPO), "rev-parse", "HEAD"], text=True).strip()
if commit != PINNED_COMMIT:
    raise RuntimeError(f"Reference checkout differs from pinned commit: {commit}")

theta = np.deg2rad([0, 0.001, 0.01, 0.1, 1, 5, 15, 30, 60, 90, 120, 150, 179, 180])
inputs = [(1.5 + 0j, v) for v in [1e-6, 1e-3, 0.0665, 0.1, 1, 10, 100, 1000, 10000]]
inputs += [(1.5 - 0.1j, v) for v in [1e-6, 0.0665, 1, 2, 10, 1000]]
inputs += [(0.75 + 0j, 10), (4/3 + 0j, 50), (1.5 - 1j, 2), (1.1 - 25j, 2),
           (1.33 - 1e-5j, 100), (1.5 - 1j, 10000), (1.0001 + 0j, 0.2),
           (0.1 - 0.01j, 5), (1 + 0j, 10)]
inputs += [(1.591 + 0j, np.pi * d / wavelength)
           for wavelength in [450e-9, 632e-9] for d in [10e-6, 100e-6, 200e-6]]
rng = np.random.default_rng(20260930)
for _ in range(24):
    inputs.append((complex(rng.uniform(0.5, 2.5), -rng.uniform(0, 0.8)),
                   float(10 ** rng.uniform(-2, 2.5))))

cases = []
for m, x in inputs:
    qext, qsca, qback, g = mie.efficiencies_mx(m, x)
    s1, s2 = mie.S1_S2(m, x, np.cos(theta), norm="wiscombe")
    # an_bn's public phase is conjugate to the coefficients used in S1_S2.
    a, b = mie.an_bn(m, x)
    cases.append(dict(mReal=m.real, mImag=m.imag, x=x,
                      efficiencies=[qext, qsca, qback, g],
                      S1Real=s1.real.tolist(), S1Imag=s1.imag.tolist(),
                      S2Real=s2.real.tolist(), S2Imag=s2.imag.tolist(),
                      aReal=a[:2].real.tolist(), aImag=a[:2].imag.tolist(),
                      bReal=b[:2].real.tolist(), bImag=b[:2].imag.tolist()))

# Independent ensemble reference, including a non-air host and volume density.
d = np.geomspace(0.05e-6, 12e-6, 65)
f_volume = np.exp(-0.5 * (np.log(d / 2e-6) / 0.45) ** 2) / d
spacing = np.diff(d)
w = np.concatenate(([spacing[0]], spacing[:-1] + spacing[1:], [spacing[-1]])) / 2
number_weights = w * f_volume / (np.pi * d**3 / 6)
number_weights /= number_weights.sum()
lambda0, n_medium, n_particle = 633e-9, 1.33, 1.59 - 0.01j
k = 2 * np.pi * n_medium / lambda0
ensemble = np.zeros_like(theta)
for dj, weight in zip(d, number_weights):
    intensity = mie.i_unpolarized(n_particle / n_medium, k * dj / 2,
                                 np.cos(theta), norm="wiscombe")
    ensemble += weight * intensity / k**2
data = dict(packageVersion=mie.__version__, commit=commit, theta=theta.tolist(), cases=cases,
            ensemble=dict(diameter=d.tolist(), volumeDensity=f_volume.tolist(),
                          lambda0=lambda0, nMedium=n_medium,
                          nParticleReal=n_particle.real, nParticleImag=n_particle.imag,
                          differentialCrossSection=ensemble.tolist()))
target = ROOT / "tests" / "fixtures" / "miepython_reference.json"
with target.open("w", encoding="utf-8", newline="\n") as stream:
    stream.write(json.dumps(data, indent=2, allow_nan=False))
print(f"Saved {len(cases)} cases and one ensemble: miepython {mie.__version__}, {commit}")
