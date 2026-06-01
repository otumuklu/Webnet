#!/usr/bin/env python3
"""Perfect-gas conical shock calculator.

This is a Python port of the conical-shock/Taylor-Maccoll routines in
``calculator.js``. Angles passed to the public API are in degrees.
"""

from __future__ import annotations

import argparse
import json
import math
from typing import Callable


DEG = math.pi / 180.0
RAD = 180.0 / math.pi


class ConicalShockError(ValueError):
    """Raised when the requested conical shock solution is not attached."""


def tt0(gamma: float, mach: float) -> float:
    return (1.0 + ((gamma - 1.0) / 2.0) * mach * mach) ** -1.0


def pp0(gamma: float, mach: float) -> float:
    return (1.0 + ((gamma - 1.0) / 2.0) * mach * mach) ** (
        -gamma / (gamma - 1.0)
    )


def rr0(gamma: float, mach: float) -> float:
    return (1.0 + ((gamma - 1.0) / 2.0) * mach * mach) ** (
        -1.0 / (gamma - 1.0)
    )


def normal_mach_after_shock(gamma: float, mach1: float) -> float:
    numerator = 1.0 + 0.5 * (gamma - 1.0) * mach1 * mach1
    denominator = gamma * mach1 * mach1 - 0.5 * (gamma - 1.0)
    return math.sqrt(numerator / denominator)


def normal_shock_ratios(gamma: float, mach1: float) -> dict[str, float]:
    mach2 = normal_mach_after_shock(gamma, mach1)
    p2p1 = 1.0 + (2.0 * gamma / (gamma + 1.0)) * (mach1 * mach1 - 1.0)
    p02p01 = (pp0(gamma, mach1) / pp0(gamma, mach2)) * p2p1

    return {
        "M1": mach1,
        "M2": mach2,
        "p2/p1": p2p1,
        "rho2/rho1": (rr0(gamma, mach2) / rr0(gamma, mach1)) * p02p01,
        "T2/T1": tt0(gamma, mach2) / tt0(gamma, mach1),
        "p02/p01": p02p01,
        "p1/p02": pp0(gamma, mach1) / p02p01,
    }


def solve_bisection(
    fn: Callable[[float], float],
    low: float,
    high: float,
    tolerance: float = 1e-10,
    max_iterations: int = 200,
) -> float:
    f_low = fn(low)
    f_high = fn(high)

    if abs(f_low) < tolerance:
        return low
    if abs(f_high) < tolerance:
        return high
    if f_low * f_high > 0.0:
        raise ConicalShockError("Could not bracket a solution for this input.")

    for _ in range(max_iterations):
        mid = 0.5 * (low + high)
        f_mid = fn(mid)
        if abs(f_mid) < tolerance or abs(high - low) < tolerance:
            return mid
        if f_low * f_mid <= 0.0:
            high = mid
            f_high = f_mid
        else:
            low = mid
            f_low = f_mid

    return 0.5 * (low + high)


def nondimensional_velocity(gamma: float, mach: float) -> float:
    return math.sqrt(((gamma - 1.0) * mach * mach) / (2.0 + (gamma - 1.0) * mach * mach))


def mach_from_nondimensional_velocity(gamma: float, velocity_magnitude: float) -> float:
    velocity_squared = velocity_magnitude * velocity_magnitude
    if velocity_squared <= 0.0 or velocity_squared >= 1.0:
        raise ConicalShockError("Taylor-Maccoll integration reached an invalid velocity state.")
    return math.sqrt((2.0 * velocity_squared) / ((gamma - 1.0) * (1.0 - velocity_squared)))


def taylor_maccoll_derivative(
    gamma: float,
    theta: float,
    radial_velocity: float,
    polar_velocity: float,
) -> tuple[float, float]:
    velocity_squared = radial_velocity * radial_velocity + polar_velocity * polar_velocity
    thermal_term = ((gamma - 1.0) / 2.0) * (1.0 - velocity_squared)
    denominator = thermal_term - polar_velocity * polar_velocity
    sin_theta = math.sin(theta)

    if (
        not math.isfinite(thermal_term)
        or not math.isfinite(denominator)
        or abs(denominator) < 1e-10
        or abs(sin_theta) < 1e-10
    ):
        raise ConicalShockError("Taylor-Maccoll integration failed for this cone angle.")

    cot_theta = math.cos(theta) / sin_theta
    numerator = (
        polar_velocity * polar_velocity * radial_velocity
        - thermal_term * (2.0 * radial_velocity + polar_velocity * cot_theta)
    )
    return polar_velocity, numerator / denominator


def integrate_taylor_maccoll(
    gamma: float,
    beta: float,
    cone_angle: float,
    radial_velocity: float,
    polar_velocity: float,
) -> tuple[float, float]:
    steps = max(24, math.ceil(abs(beta - cone_angle) / 0.00035))
    step = (cone_angle - beta) / steps
    theta = beta
    vr = radial_velocity
    vt = polar_velocity

    for _ in range(steps):
        k1r, k1t = taylor_maccoll_derivative(gamma, theta, vr, vt)
        k2r, k2t = taylor_maccoll_derivative(
            gamma, theta + step / 2.0, vr + (step * k1r) / 2.0, vt + (step * k1t) / 2.0
        )
        k3r, k3t = taylor_maccoll_derivative(
            gamma, theta + step / 2.0, vr + (step * k2r) / 2.0, vt + (step * k2t) / 2.0
        )
        k4r, k4t = taylor_maccoll_derivative(
            gamma, theta + step, vr + step * k3r, vt + step * k3t
        )

        vr += (step / 6.0) * (k1r + 2.0 * k2r + 2.0 * k3r + k4r)
        vt += (step / 6.0) * (k1t + 2.0 * k2t + 2.0 * k3t + k4t)
        theta += step

        if not math.isfinite(vr) or not math.isfinite(vt) or vr * vr + vt * vt >= 1.0:
            raise ConicalShockError("Taylor-Maccoll integration failed for this cone angle.")

    return vr, vt


def conical_shock_state_at_beta(
    gamma: float,
    mach1: float,
    beta: float,
    cone_angle: float,
) -> dict[str, object]:
    velocity1 = nondimensional_velocity(gamma, mach1)
    normal_mach1 = mach1 * math.sin(beta)
    if normal_mach1 <= 1.0:
        raise ConicalShockError("Shock angle must produce a supersonic normal Mach number.")

    shock_ratios = normal_shock_ratios(gamma, normal_mach1)
    density_ratio = shock_ratios["rho2/rho1"]
    radial_velocity = velocity1 * math.cos(beta)
    polar_velocity = -velocity1 * math.sin(beta) / density_ratio
    surface = integrate_taylor_maccoll(
        gamma, beta, cone_angle, radial_velocity, polar_velocity
    )

    return {
        "radial_velocity": radial_velocity,
        "polar_velocity": polar_velocity,
        "surface_radial_velocity": surface[0],
        "surface_polar_velocity": surface[1],
        "shock_ratios": shock_ratios,
    }


def cone_shock_beta(gamma: float, mach1: float, cone_angle: float) -> float | None:
    mach_angle = math.asin(1.0 / mach1)
    beta_low = max(mach_angle, cone_angle) + 1e-5
    beta_high = math.pi / 2.0 - 1e-5
    previous_beta = None
    previous_value = None

    for i in range(261):
        beta = beta_low + (i / 260.0) * (beta_high - beta_low)
        try:
            value = conical_shock_state_at_beta(
                gamma, mach1, beta, cone_angle
            )["surface_polar_velocity"]
        except ConicalShockError:
            continue

        if abs(value) < 1e-8:
            return beta
        if previous_value is not None and previous_value * value < 0.0:
            return solve_bisection(
                lambda candidate: conical_shock_state_at_beta(
                    gamma, mach1, candidate, cone_angle
                )["surface_polar_velocity"],
                previous_beta,
                beta,
                1e-9,
                120,
            )

        previous_beta = beta
        previous_value = value

    return None


def calculate_conical_shock(
    mach1: float,
    cone_angle_degrees: float,
    gamma: float = 1.4,
    pressure1: float | None = None,
    temperature1: float | None = None,
    gas_constant: float = 287.05,
) -> dict[str, float]:
    """Return attached conical shock angle and downstream flow parameters."""

    if gamma <= 1.0:
        raise ConicalShockError("Gamma must be greater than 1.")
    if mach1 <= 1.0 or not math.isfinite(mach1):
        raise ConicalShockError("M1 must be greater than 1.")
    if cone_angle_degrees <= 0.0 or not math.isfinite(cone_angle_degrees):
        raise ConicalShockError("Cone half-angle must be greater than 0 degrees.")
    if cone_angle_degrees >= 89.0:
        raise ConicalShockError("Cone half-angle must be less than 89 degrees.")
    if (pressure1 is None) != (temperature1 is None):
        raise ConicalShockError("Provide both pressure1 and temperature1, or neither.")
    if pressure1 is not None and pressure1 <= 0.0:
        raise ConicalShockError("pressure1 must be greater than 0.")
    if temperature1 is not None and temperature1 <= 0.0:
        raise ConicalShockError("temperature1 must be greater than 0.")

    cone_angle = cone_angle_degrees * DEG
    beta = cone_shock_beta(gamma, mach1, cone_angle)
    if beta is None:
        raise ConicalShockError(
            f"No attached conical shock found for M1 = {mach1:g} "
            f"and cone angle = {cone_angle_degrees:g} degrees."
        )
    if beta <= cone_angle:
        raise ConicalShockError(
            "Invalid conical shock solution: shock angle must be greater than the cone half-angle."
        )

    state = conical_shock_state_at_beta(gamma, mach1, beta, cone_angle)
    shock_velocity = math.hypot(state["radial_velocity"], state["polar_velocity"])
    surface_velocity = math.hypot(
        state["surface_radial_velocity"], state["surface_polar_velocity"]
    )
    shock_mach = mach_from_nondimensional_velocity(gamma, shock_velocity)
    surface_mach = mach_from_nondimensional_velocity(gamma, surface_velocity)

    shock_ratios = state["shock_ratios"]
    shock_temperature_ratio = shock_ratios["T2/T1"]
    shock_pressure_ratio = shock_ratios["p2/p1"]
    shock_density_ratio = shock_ratios["rho2/rho1"]
    surface_to_shock_temperature_ratio = (1.0 - surface_velocity * surface_velocity) / (
        1.0 - shock_velocity * shock_velocity
    )
    surface_to_shock_pressure_ratio = surface_to_shock_temperature_ratio ** (
        gamma / (gamma - 1.0)
    )
    surface_to_shock_density_ratio = surface_to_shock_temperature_ratio ** (
        1.0 / (gamma - 1.0)
    )

    result = {
        "M1": mach1,
        "gamma": gamma,
        "cone_half_angle_deg": cone_angle_degrees,
        "shock_angle_deg": beta * RAD,
        "Mn1": mach1 * math.sin(beta),
        "M2_behind_shock": shock_mach,
        "p2/p1": shock_pressure_ratio,
        "T2/T1": shock_temperature_ratio,
        "rho2/rho1": shock_density_ratio,
        "Mc_cone_surface": surface_mach,
        "pc/p1": shock_pressure_ratio * surface_to_shock_pressure_ratio,
        "Tc/T1": shock_temperature_ratio * surface_to_shock_temperature_ratio,
        "rhoc/rho1": shock_density_ratio * surface_to_shock_density_ratio,
    }

    if pressure1 is not None and temperature1 is not None:
        rho1 = pressure1 / (gas_constant * temperature1)
        result.update(
            {
                "p1": pressure1,
                "T1": temperature1,
                "rho1": rho1,
                "p2": pressure1 * result["p2/p1"],
                "T2": temperature1 * result["T2/T1"],
                "rho2": rho1 * result["rho2/rho1"],
                "pc": pressure1 * result["pc/p1"],
                "Tc": temperature1 * result["Tc/T1"],
                "rhoc": rho1 * result["rhoc/rho1"],
            }
        )

    return result


def print_table(values: dict[str, float]) -> None:
    width = max(len(key) for key in values)
    for key, value in values.items():
        print(f"{key:<{width}} : {value:.10g}")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Calculate attached perfect-gas conical shock angles and flow parameters."
    )
    parser.add_argument("--mach", "-M", type=float, required=True, help="Freestream Mach number M1.")
    parser.add_argument(
        "--cone-angle",
        "-c",
        type=float,
        required=True,
        help="Cone half-angle in degrees.",
    )
    parser.add_argument("--gamma", "-g", type=float, default=1.4, help="Specific heat ratio.")
    parser.add_argument("--p1", type=float, help="Optional freestream static pressure.")
    parser.add_argument("--T1", type=float, help="Optional freestream static temperature.")
    parser.add_argument(
        "--R",
        type=float,
        default=287.05,
        help="Gas constant for dimensional density, default air in J/(kg K).",
    )
    parser.add_argument("--json", action="store_true", help="Print results as JSON.")
    args = parser.parse_args()

    try:
        result = calculate_conical_shock(
            mach1=args.mach,
            cone_angle_degrees=args.cone_angle,
            gamma=args.gamma,
            pressure1=args.p1,
            temperature1=args.T1,
            gas_constant=args.R,
        )
    except ConicalShockError as exc:
        parser.error(str(exc))

    if args.json:
        print(json.dumps(result, indent=2))
    else:
        print_table(result)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
