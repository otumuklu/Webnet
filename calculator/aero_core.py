"""Core compressible-flow relations for the HAVA website calculator.

The browser calculator uses a JavaScript port of these formulas so the static
website remains interactive. Keep new Python-only modules, such as hypersonics
models, similarly small: inputs in, a dictionary of named outputs out.
"""

from __future__ import annotations

from math import asin, atan, isfinite, pi, pow, sin, sqrt
from typing import Callable


DEG = pi / 180.0
RAD = 180.0 / pi


def _check_gamma(gamma: float) -> None:
    if gamma <= 1.0:
        raise ValueError("gamma must be greater than 1")


def _check_number(value: float, name: str = "value") -> None:
    if not isfinite(value):
        raise ValueError(f"{name} must be finite")


def tt0(gamma: float, mach: float) -> float:
    return pow(1.0 + (gamma - 1.0) * mach * mach / 2.0, -1.0)


def pp0(gamma: float, mach: float) -> float:
    return pow(1.0 + (gamma - 1.0) * mach * mach / 2.0, -gamma / (gamma - 1.0))


def rr0(gamma: float, mach: float) -> float:
    return pow(1.0 + (gamma - 1.0) * mach * mach / 2.0, -1.0 / (gamma - 1.0))


def tts(gamma: float, mach: float) -> float:
    return tt0(gamma, mach) * (gamma + 1.0) / 2.0


def pps(gamma: float, mach: float) -> float:
    return pp0(gamma, mach) * pow((gamma + 1.0) / 2.0, gamma / (gamma - 1.0))


def rrs(gamma: float, mach: float) -> float:
    return rr0(gamma, mach) * pow((gamma + 1.0) / 2.0, 1.0 / (gamma - 1.0))


def area_ratio(gamma: float, mach: float) -> float:
    return (1.0 / rrs(gamma, mach)) * sqrt(1.0 / tts(gamma, mach)) / mach


def prandtl_meyer(gamma: float, mach: float) -> float:
    if mach < 1.0:
        raise ValueError("Prandtl-Meyer angle is defined for Mach >= 1")
    scale = sqrt((gamma + 1.0) / (gamma - 1.0))
    arg = sqrt((gamma - 1.0) / (gamma + 1.0) * (mach * mach - 1.0))
    return (scale * atan(arg) - atan(sqrt(mach * mach - 1.0))) * RAD


def normal_mach_after_shock(gamma: float, mach1: float) -> float:
    return sqrt(
        (1.0 + 0.5 * (gamma - 1.0) * mach1 * mach1)
        / (gamma * mach1 * mach1 - 0.5 * (gamma - 1.0))
    )


def solve_bisection(
    function: Callable[[float], float],
    low: float,
    high: float,
    tolerance: float = 1.0e-10,
    max_iterations: int = 200,
) -> float:
    f_low = function(low)
    f_high = function(high)
    if abs(f_low) < tolerance:
        return low
    if abs(f_high) < tolerance:
        return high
    if f_low * f_high > 0:
        raise ValueError("could not bracket a solution")

    for _ in range(max_iterations):
        mid = 0.5 * (low + high)
        f_mid = function(mid)
        if abs(f_mid) < tolerance or abs(high - low) < tolerance:
            return mid
        if f_low * f_mid <= 0:
            high = mid
            f_high = f_mid
        else:
            low = mid
            f_low = f_mid

    return 0.5 * (low + high)


def mach_from_area(gamma: float, area: float, branch: str) -> float:
    if area < 1.0:
        raise ValueError("A/A* must be greater than or equal to 1")
    if abs(area - 1.0) < 1.0e-12:
        return 1.0
    if branch == "sub":
        return solve_bisection(lambda mach: area_ratio(gamma, mach) - area, 1.0e-8, 0.999999)
    if branch == "sup":
        return solve_bisection(lambda mach: area_ratio(gamma, mach) - area, 1.000001, 100.0)
    raise ValueError("branch must be 'sub' or 'sup'")


def mach_from_prandtl_meyer(gamma: float, angle_degrees: float) -> float:
    max_angle = (sqrt((gamma + 1.0) / (gamma - 1.0)) - 1.0) * 90.0
    if angle_degrees <= 0.0 or angle_degrees >= max_angle:
        raise ValueError(f"Prandtl-Meyer angle must be between 0 and {max_angle}")
    return solve_bisection(
        lambda mach: prandtl_meyer(gamma, mach) - angle_degrees,
        1.000001,
        100.0,
    )


def isentropic_from_input(gamma: float, input_name: str, value: float) -> dict[str, float | None]:
    _check_gamma(gamma)
    _check_number(value)

    if input_name == "mach":
        if value <= 0.0:
            raise ValueError("Mach number must be greater than 0")
        mach = value
    elif input_name == "tt0":
        if value <= 0.0 or value >= 1.0:
            raise ValueError("T/T0 must be between 0 and 1")
        mach = sqrt(2.0 * (1.0 / value - 1.0) / (gamma - 1.0))
    elif input_name == "pp0":
        if value <= 0.0 or value >= 1.0:
            raise ValueError("p/p0 must be between 0 and 1")
        mach = sqrt(2.0 * (pow(1.0 / value, (gamma - 1.0) / gamma) - 1.0) / (gamma - 1.0))
    elif input_name == "rr0":
        if value <= 0.0 or value >= 1.0:
            raise ValueError("rho/rho0 must be between 0 and 1")
        mach = sqrt(2.0 * (pow(1.0 / value, gamma - 1.0) - 1.0) / (gamma - 1.0))
    elif input_name == "area-sub":
        mach = mach_from_area(gamma, value, "sub")
    elif input_name == "area-sup":
        mach = mach_from_area(gamma, value, "sup")
    elif input_name == "mach-angle":
        if value <= 0.0 or value >= 90.0:
            raise ValueError("Mach angle must be between 0 and 90 degrees")
        mach = 1.0 / sin(value * DEG)
    elif input_name == "pm-angle":
        mach = mach_from_prandtl_meyer(gamma, value)
    else:
        raise ValueError(f"unknown isentropic input: {input_name}")

    return {
        "Mach number": mach,
        "Mach angle": asin(1.0 / mach) * RAD if mach >= 1.0 else None,
        "P-M angle": prandtl_meyer(gamma, mach) if mach >= 1.0 else None,
        "T/T0": tt0(gamma, mach),
        "p/p0": pp0(gamma, mach),
        "rho/rho0": rr0(gamma, mach),
        "T/T*": tts(gamma, mach),
        "p/p*": pps(gamma, mach),
        "rho/rho*": rrs(gamma, mach),
        "A/A*": area_ratio(gamma, mach),
    }


def normal_shock_ratios(gamma: float, mach1: float) -> dict[str, float]:
    mach2 = normal_mach_after_shock(gamma, mach1)
    p2p1 = 1.0 + 2.0 * gamma / (gamma + 1.0) * (mach1 * mach1 - 1.0)
    p02p01 = pp0(gamma, mach1) / pp0(gamma, mach2) * p2p1

    return {
        "M1": mach1,
        "M2": mach2,
        "p2/p1": p2p1,
        "rho2/rho1": rr0(gamma, mach2) / rr0(gamma, mach1) * p02p01,
        "T2/T1": tt0(gamma, mach2) / tt0(gamma, mach1),
        "p02/p01": p02p01,
        "p1/p02": pp0(gamma, mach1) / p02p01,
    }


def normal_shock_from_input(gamma: float, input_name: str, value: float) -> dict[str, float]:
    _check_gamma(gamma)
    _check_number(value)

    if input_name == "m1":
        if value <= 1.0:
            raise ValueError("M1 must be greater than 1")
        mach1 = value
    elif input_name == "m2":
        min_m2 = sqrt((gamma - 1.0) / (2.0 * gamma))
        if value <= min_m2 or value >= 1.0:
            raise ValueError(f"M2 must be between {min_m2} and 1")
        mach1 = sqrt((1.0 + 0.5 * (gamma - 1.0) * value * value) / (gamma * value * value - 0.5 * (gamma - 1.0)))
    elif input_name == "p2p1":
        if value <= 1.0:
            raise ValueError("p2/p1 must be greater than 1")
        mach1 = sqrt((value - 1.0) * (gamma + 1.0) / (2.0 * gamma) + 1.0)
    elif input_name == "r2r1":
        max_density_ratio = (gamma + 1.0) / (gamma - 1.0)
        if value <= 1.0 or value >= max_density_ratio:
            raise ValueError(f"rho2/rho1 must be between 1 and {max_density_ratio}")
        mach1 = sqrt(2.0 * value / (gamma + 1.0 - value * (gamma - 1.0)))
    elif input_name == "t2t1":
        if value <= 1.0:
            raise ValueError("T2/T1 must be greater than 1")
        mach1 = solve_bisection(lambda mach: normal_shock_ratios(gamma, mach)["T2/T1"] - value, 1.000001, 100.0)
    elif input_name == "p02p01":
        if value <= 0.0 or value >= 1.0:
            raise ValueError("p02/p01 must be between 0 and 1")
        mach1 = solve_bisection(lambda mach: normal_shock_ratios(gamma, mach)["p02/p01"] - value, 1.000001, 100.0)
    elif input_name == "p1p02":
        max_value = pow((gamma + 1.0) / 2.0, -gamma / (gamma - 1.0))
        if value <= 0.0 or value >= max_value:
            raise ValueError(f"p1/p02 must be between 0 and {max_value}")
        mach1 = solve_bisection(lambda mach: normal_shock_ratios(gamma, mach)["p1/p02"] - value, 1.000001, 100.0)
    else:
        raise ValueError(f"unknown normal-shock input: {input_name}")

    return normal_shock_ratios(gamma, mach1)


UNIVERSAL_GAS_CONSTANT = 8314.462618


GAS_PROPERTIES = {
    "air": {"label": "Air", "gamma": 1.4, "molecular_weight": 28.965},
    "co2": {"label": "CO2", "gamma": 1.289, "molecular_weight": 44.01},
    "argon": {"label": "Argon", "gamma": 1.667, "molecular_weight": 39.948},
    "nitrogen": {"label": "N2", "gamma": 1.4, "molecular_weight": 28.0134},
    "oxygen": {"label": "O2", "gamma": 1.4, "molecular_weight": 31.9988},
    "helium": {"label": "Helium", "gamma": 1.667, "molecular_weight": 4.0026},
}


def gas_constant_from_molecular_weight(molecular_weight: float) -> float:
    """Return specific gas constant in J/(kg K) from MW in kg/kmol or g/mol."""
    _check_number(molecular_weight, "molecular_weight")
    if molecular_weight <= 0.0:
        raise ValueError("molecular_weight must be greater than 0")
    return UNIVERSAL_GAS_CONSTANT / molecular_weight


def speed_of_sound(gamma: float, gas_constant: float, temperature: float) -> float:
    """Return speed of sound in m/s for a perfect gas, a = sqrt(gamma R T)."""
    _check_gamma(gamma)
    _check_number(gas_constant, "gas_constant")
    _check_number(temperature, "temperature")
    if gas_constant <= 0.0:
        raise ValueError("gas_constant must be greater than 0")
    if temperature <= 0.0:
        raise ValueError("temperature must be greater than 0 K")
    return sqrt(gamma * gas_constant * temperature)


def speed_of_sound_from_gas(gas: str, temperature: float) -> dict[str, float | str]:
    properties = GAS_PROPERTIES[gas]
    gas_constant = gas_constant_from_molecular_weight(properties["molecular_weight"])
    speed = speed_of_sound(properties["gamma"], gas_constant, temperature)
    return {
        "Gas": properties["label"],
        "gamma": properties["gamma"],
        "MW, kg/kmol": properties["molecular_weight"],
        "R, J/(kg K)": gas_constant,
        "T, K": temperature,
        "a, m/s": speed,
        "a, ft/s": speed * 3.280839895,
    }
