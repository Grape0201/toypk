from .db import db
from scipy.interpolate import Akima1DInterpolator


def get_linear_attenuation_coefficient(
    energy: float, densities: dict[str, float]
) -> float:
    """
    Parameters
    ----------
    energy : float
        photon ray energy [MeV]
    densities: dit[str, float]
        the key is the element like H, the value is the density [g/cm3]

    Returns
    -------
    lac: float
        linear attenuation coefficient
    """
    lac = 0.0
    for element, rho in densities.items():
        mu_rho = Akima1DInterpolator(db[element]["e"], db[element]["mu_rho"])(energy)
        lac += mu_rho * rho
    return lac
