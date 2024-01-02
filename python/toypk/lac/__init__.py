from .db import db
from scipy.interpolate import Akima1DInterpolator
import numpy as np


def get_linear_attenuation_coefficient(
    energy_points: list[float], densities: dict[str, float]
) -> list[float]:
    """
    Parameters
    ----------
    energy_points : list[float]
        photon ray energy [MeV]
    densities: dit[str, float]
        the key is the element like "H", the value is the density [g/cm3]

    Returns
    -------
    lac: list[float]
        Akima-interpolated linear attenuation coefficients
    """
    lac = np.zeros(len(energy_points))
    for element, rho in densities.items():
        a = Akima1DInterpolator(db[element]["e"], db[element]["mu_rho"])
        lac += np.array(list(map(lambda e: rho * a(e), energy_points)))
    return lac.tolist()
