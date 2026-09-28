"""Two-element Windkessel dynamics with an optional proximal resistance.

Units: seconds, mmHg, mL; R and Z in mmHg*s/mL, C in mL/mmHg.
The generated inflow is a prescribed sinusoidal pump, not a heart model.
"""
from dataclasses import dataclass
import numpy as np
from scipy.signal import lfilter
from scipy.optimize import least_squares

@dataclass(frozen=True)
class Patient:
    resistance: float
    compliance: float
    stroke_volume: float
    heart_rate: float
    impedance: float = 0.0

    def __post_init__(self):
        vals = [self.resistance, self.compliance, self.stroke_volume, self.heart_rate]
        if not all(np.isfinite(vals)) or min(vals) <= 0:
            raise ValueError("Physiological parameters must be positive and finite")
        if not np.isfinite(self.impedance) or self.impedance < 0:
            raise ValueError("Impedance must be finite and nonnegative")


def flow(t, heart_rate, stroke_volume, multiplier=1.0):
    """Nonnegative prescribed inflow with mean HR*SV/60 (mL/s)."""
    return multiplier * heart_rate * stroke_volume / 60 * (1 + 0.8*np.sin(2*np.pi*heart_rate*t/60))


def simulate(t, patient, multiplier=1.0, initial_pressure=85.0, venous_pressure=5.0):
    """Exact zero-order-hold capacitor update; input sampled at left endpoints.

    Returns arterial pressure Pc + Z*Q and inflow. Initial pressure is Pc.
    A nonzero Z creates a three-element generator for model-mismatch tests.
    """
    t = np.asarray(t, dtype=float)
    if len(t) < 2 or not np.all(np.isfinite(t)) or not np.allclose(np.diff(t), t[1]-t[0]) or t[1] <= t[0]:
        raise ValueError("Need a finite, increasing, uniform time grid")
    q = flow(t, patient.heart_rate, patient.stroke_volume, multiplier)
    if np.shape(q) != t.shape or not np.all(np.isfinite(q)) or np.any(q < 0):
        raise ValueError("Flow multiplier must produce finite nonnegative flow")
    a = np.exp(-(t[1]-t[0])/(patient.resistance*patient.compliance))
    drive = venous_pressure + patient.resistance*q[:-1]
    tail, _ = lfilter([1-a], [1, -a], drive, zi=[a*initial_pressure])
    pc = np.r_[initial_pressure, tail]
    return pc + patient.impedance*q, q


def fit_patient(t, observed, hr, sv):
    """Fit R, C and initial Pc using available observations and known HR/SV.

    Missing observations are NaN. No future observations may be supplied.
    """
    observed = np.asarray(observed)
    mask = np.isfinite(observed)
    if mask.sum() < 8 or len(observed) != len(t):
        raise ValueError("At least eight observations on the supplied grid required")
    def residual(x):
        p, _ = simulate(t, Patient(x[0], x[1], sv, hr), initial_pressure=x[2])
        return p[mask]-observed[mask]
    fit = least_squares(residual, [1.1, 1.5, 85], bounds=([.2,.2,30],[3.5,5,180]), max_nfev=150)
    if not fit.success:
        raise RuntimeError(f"Parameter fit failed: {fit.message}")
    return fit.x
