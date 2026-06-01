#!/usr/bin/env python3
"""Combined Hansen equilibrium-air and conical shock calculator.

This script keeps the fast scalar Hansen chemistry/thermodynamics in the same
file as the Taylor-Maccoll conical-shock solver. It can run a perfect-gas cone
shock, evaluate Hansen equilibrium air at any state, or solve the reacting
Hansen cone-shock model used by ``calculator.js``.
"""

from __future__ import annotations

import argparse
import json
import math
from functools import lru_cache
from typing import Callable


DEG = math.pi / 180.0
RAD = 180.0 / math.pi
ATM_TO_PA = 101325.0
R_UNIVERSAL = 8.31446261815324
DEFAULT_XI = 0.25
SPECIES = ("N2", "O2", "N", "O", "O+", "N+", "e-")
MOLAR_MASS = (28.02, 32.00, 14.01, 16.00, 16.00, 14.01, 5.485799e-4)
SPECIFIC_R = tuple(8314.462618 / mw for mw in MOLAR_MASS)
FORMATION_OVER_R = {
    "N2": 0.0,
    "O2": 0.0,
    "O": 59000.0 / 2.0,
    "N": 113200.0 / 2.0,
    "O+": 59000.0 / 2.0 + 158000.0,
    "N+": 113200.0 / 2.0 + 168800.0,
    "e-": 0.0,
}


class SolverError(ValueError):
    """Raised when a requested gas-dynamics or Hansen solve is invalid."""


def safe_exp(value: float) -> float:
    return math.exp(max(-700.0, min(700.0, value)))


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
        raise SolverError("Could not bracket a solution for this input.")

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


def tt0(gamma: float, mach: float) -> float:
    return (1.0 + 0.5 * (gamma - 1.0) * mach * mach) ** -1.0


def pp0(gamma: float, mach: float) -> float:
    return (1.0 + 0.5 * (gamma - 1.0) * mach * mach) ** (-gamma / (gamma - 1.0))


def rr0(gamma: float, mach: float) -> float:
    return (1.0 + 0.5 * (gamma - 1.0) * mach * mach) ** (-1.0 / (gamma - 1.0))


def normal_mach_after_shock(gamma: float, mach1: float) -> float:
    return math.sqrt(
        (1.0 + 0.5 * (gamma - 1.0) * mach1 * mach1)
        / (gamma * mach1 * mach1 - 0.5 * (gamma - 1.0))
    )


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


def nondimensional_velocity(gamma: float, mach: float) -> float:
    return math.sqrt(((gamma - 1.0) * mach * mach) / (2.0 + (gamma - 1.0) * mach * mach))


def mach_from_nondimensional_velocity(gamma: float, velocity: float) -> float:
    v2 = velocity * velocity
    if v2 <= 0.0 or v2 >= 1.0:
        raise SolverError("Taylor-Maccoll integration reached an invalid velocity state.")
    return math.sqrt((2.0 * v2) / ((gamma - 1.0) * (1.0 - v2)))


def taylor_maccoll_derivative(gamma: float, theta: float, vr: float, vt: float) -> tuple[float, float]:
    v2 = vr * vr + vt * vt
    thermal = 0.5 * (gamma - 1.0) * (1.0 - v2)
    denominator = thermal - vt * vt
    sin_theta = math.sin(theta)
    if (
        not math.isfinite(thermal)
        or not math.isfinite(denominator)
        or abs(denominator) < 1e-10
        or abs(sin_theta) < 1e-10
    ):
        raise SolverError("Taylor-Maccoll integration failed for this cone angle.")
    cot_theta = math.cos(theta) / sin_theta
    numerator = vt * vt * vr - thermal * (2.0 * vr + vt * cot_theta)
    return vt, numerator / denominator


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
            gamma, theta + step / 2.0, vr + step * k1r / 2.0, vt + step * k1t / 2.0
        )
        k3r, k3t = taylor_maccoll_derivative(
            gamma, theta + step / 2.0, vr + step * k2r / 2.0, vt + step * k2t / 2.0
        )
        k4r, k4t = taylor_maccoll_derivative(gamma, theta + step, vr + step * k3r, vt + step * k3t)
        vr += step * (k1r + 2.0 * k2r + 2.0 * k3r + k4r) / 6.0
        vt += step * (k1t + 2.0 * k2t + 2.0 * k3t + k4t) / 6.0
        theta += step
        if not math.isfinite(vr) or not math.isfinite(vt) or vr * vr + vt * vt >= 1.0:
            raise SolverError("Taylor-Maccoll integration failed for this cone angle.")
    return vr, vt


def conical_shock_state_at_beta(gamma: float, mach1: float, beta: float, cone_angle: float) -> dict[str, object]:
    velocity1 = nondimensional_velocity(gamma, mach1)
    normal_mach1 = mach1 * math.sin(beta)
    if normal_mach1 <= 1.0:
        raise SolverError("Shock angle must produce a supersonic normal Mach number.")
    shock_ratios = normal_shock_ratios(gamma, normal_mach1)
    density_ratio = shock_ratios["rho2/rho1"]
    radial_velocity = velocity1 * math.cos(beta)
    polar_velocity = -velocity1 * math.sin(beta) / density_ratio
    surface_radial, surface_polar = integrate_taylor_maccoll(
        gamma, beta, cone_angle, radial_velocity, polar_velocity
    )
    return {
        "radial_velocity": radial_velocity,
        "polar_velocity": polar_velocity,
        "surface_radial_velocity": surface_radial,
        "surface_polar_velocity": surface_polar,
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
            value = conical_shock_state_at_beta(gamma, mach1, beta, cone_angle)["surface_polar_velocity"]
        except SolverError:
            continue
        if abs(value) < 1e-8:
            return beta
        if previous_value is not None and previous_value * value < 0.0:
            return solve_bisection(
                lambda b: conical_shock_state_at_beta(gamma, mach1, b, cone_angle)["surface_polar_velocity"],
                previous_beta,
                beta,
                1e-9,
                120,
            )
        previous_beta = beta
        previous_value = value
    return None


def calculate_perfect_cone(
    mach1: float,
    cone_angle_degrees: float,
    gamma: float = 1.4,
    pressure1: float | None = None,
    temperature1: float | None = None,
    gas_constant: float = 287.05,
) -> dict[str, float]:
    if gamma <= 1.0:
        raise SolverError("Gamma must be greater than 1.")
    if mach1 <= 1.0 or not math.isfinite(mach1):
        raise SolverError("M1 must be greater than 1.")
    if cone_angle_degrees <= 0.0 or cone_angle_degrees >= 89.0:
        raise SolverError("Cone half-angle must be between 0 and 89 degrees.")
    if (pressure1 is None) != (temperature1 is None):
        raise SolverError("Provide both p1 and T1, or neither.")

    cone_angle = cone_angle_degrees * DEG
    beta = cone_shock_beta(gamma, mach1, cone_angle)
    if beta is None:
        raise SolverError("No attached perfect-gas conical shock found.")
    state = conical_shock_state_at_beta(gamma, mach1, beta, cone_angle)
    shock_velocity = math.hypot(state["radial_velocity"], state["polar_velocity"])
    surface_velocity = math.hypot(state["surface_radial_velocity"], state["surface_polar_velocity"])
    shock_mach = mach_from_nondimensional_velocity(gamma, shock_velocity)
    surface_mach = mach_from_nondimensional_velocity(gamma, surface_velocity)
    ratios = state["shock_ratios"]
    tc_t2 = (1.0 - surface_velocity * surface_velocity) / (1.0 - shock_velocity * shock_velocity)
    pc_p2 = tc_t2 ** (gamma / (gamma - 1.0))
    rhoc_rho2 = tc_t2 ** (1.0 / (gamma - 1.0))
    result = {
        "M1": mach1,
        "gamma": gamma,
        "cone_half_angle_deg": cone_angle_degrees,
        "shock_angle_deg": beta * RAD,
        "Mn1": mach1 * math.sin(beta),
        "M2_behind_shock": shock_mach,
        "p2/p1": ratios["p2/p1"],
        "T2/T1": ratios["T2/T1"],
        "rho2/rho1": ratios["rho2/rho1"],
        "Mc_cone_surface": surface_mach,
        "pc/p1": ratios["p2/p1"] * pc_p2,
        "Tc/T1": ratios["T2/T1"] * tc_t2,
        "rhoc/rho1": ratios["rho2/rho1"] * rhoc_rho2,
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


def feed_fractions(xi: float) -> tuple[float, float]:
    if xi <= 0.0 or not math.isfinite(xi):
        raise SolverError("O/N atomic ratio xi must be positive.")
    return 1.0 / (1.0 + xi), xi / (1.0 + xi)


def mixture_gas_constant(xi: float) -> float:
    x_n2, x_o2 = feed_fractions(xi)
    molecular_mass_kg = (x_n2 * MOLAR_MASS[0] + x_o2 * MOLAR_MASS[1]) / 1000.0
    return R_UNIVERSAL / molecular_mass_kg


@lru_cache(maxsize=20000)
def ln_qp(species: str, temperature: float) -> float:
    t = temperature
    if t <= 0.0:
        raise SolverError("Temperature must be greater than 0 K.")
    if species == "N2":
        return 3.5 * math.log(t) - 0.42 - math.log(1.0 - safe_exp(-3390.0 / t))
    if species == "O2":
        return (
            3.5 * math.log(t)
            + 0.11
            - math.log(1.0 - safe_exp(-2270.0 / t))
            + math.log(3.0 + 2.0 * safe_exp(-11390.0 / t) + safe_exp(-18990.0 / t))
        )
    if species == "O":
        return (
            2.5 * math.log(t)
            + 0.50
            + math.log(
                5.0
                + 3.0 * safe_exp(-228.0 / t)
                + safe_exp(-326.0 / t)
                + 5.0 * safe_exp(-22800.0 / t)
                + safe_exp(-48600.0 / t)
            )
        )
    if species == "N":
        return (
            2.5 * math.log(t)
            + 0.30
            + math.log(4.0 + 10.0 * safe_exp(-27700.0 / t) + 6.0 * safe_exp(-41500.0 / t))
        )
    if species == "O+":
        return (
            2.5 * math.log(t)
            + 0.50
            + math.log(4.0 + 10.0 * safe_exp(-38600.0 / t) + 6.0 * safe_exp(-58200.0 / t))
        )
    if species == "N+":
        return (
            2.5 * math.log(t)
            + 0.30
            + math.log(
                1.0
                + 3.0 * safe_exp(-70.6 / t)
                + 5.0 * safe_exp(-188.9 / t)
                + 5.0 * safe_exp(-22000.0 / t)
                + safe_exp(-47000.0 / t)
                + 5.0 * safe_exp(-67900.0 / t)
            )
        )
    if species == "e-":
        return 2.5 * math.log(t) - 14.24
    raise SolverError(f"Unknown Hansen species: {species}.")


def equilibrium_constants(temperature: float, xi: float = DEFAULT_XI) -> tuple[float, float, float]:
    x_n2, x_o2 = feed_fractions(xi)
    t = temperature
    kp_o2 = safe_exp(-59000.0 / t + 2.0 * ln_qp("O", t) - ln_qp("O2", t))
    kp_n2 = safe_exp(-113200.0 / t + 2.0 * ln_qp("N", t) - ln_qp("N2", t))
    kp_o_ion = safe_exp(-158000.0 / t + ln_qp("O+", t) + ln_qp("e-", t) - ln_qp("O", t))
    kp_n_ion = safe_exp(-168800.0 / t + ln_qp("N+", t) + ln_qp("e-", t) - ln_qp("N", t))
    return kp_o2, kp_n2, x_o2 * kp_o_ion + x_n2 * kp_n_ion


def hansen_epsilons(temperature: float, pressure_atm: float, xi: float = DEFAULT_XI) -> tuple[float, float, float]:
    kp_o2, kp_n2, kp_ion = equilibrium_constants(temperature, xi)
    a1 = 1.0 + 4.0 * pressure_atm / kp_o2
    a2 = 1.0 + 4.0 * pressure_atm / kp_n2
    ion_denominator = 1.0 + pressure_atm / kp_ion
    eps1 = (-0.8 + math.sqrt(0.64 + 0.8 * a1)) / (2.0 * a1) if math.isfinite(a1) else 0.0
    eps2 = (-0.4 + math.sqrt(0.16 + 3.84 * a2)) / (2.0 * a2) if math.isfinite(a2) else 0.0
    eps3 = ion_denominator ** -0.5 if math.isfinite(ion_denominator) else 0.0
    return eps1, eps2, eps3


@lru_cache(maxsize=50000)
def hansen_air_cached(temperature: float, pressure_atm: float, xi: float = DEFAULT_XI) -> tuple[float, ...]:
    if pressure_atm <= 0.0 or not math.isfinite(pressure_atm):
        raise SolverError("Pressure must be greater than 0 atm.")
    x_n2_feed, x_o2_feed = feed_fractions(xi)
    if temperature <= 1200.0:
        return (1.0, 0.0, 0.0, 0.0, x_n2_feed, x_o2_feed, 0.0, 0.0, 0.0, 0.0, 0.0)
    eps1_air, eps2_air, eps3 = hansen_epsilons(temperature, pressure_atm, xi)
    eps1 = min((eps1_air * x_o2_feed) / 0.2, x_o2_feed)
    eps2 = min((eps2_air * x_n2_feed) / 0.8, x_n2_feed)
    n_share = x_n2_feed / (x_n2_feed + x_o2_feed)
    o_share = x_o2_feed / (x_n2_feed + x_o2_feed)
    z = 1.0 + eps1 + eps2 + 2.0 * eps3
    xe = 2.0 * eps3 / z
    return (
        z,
        eps1,
        eps2,
        eps3,
        max(0.0, (x_n2_feed - eps2) / z),
        max(0.0, (x_o2_feed - eps1) / z),
        max(0.0, (2.0 * eps2 - 2.0 * n_share * eps3) / z),
        max(0.0, (2.0 * eps1 - 2.0 * o_share * eps3) / z),
        max(0.0, o_share * xe),
        max(0.0, n_share * xe),
        max(0.0, xe),
    )


def hansen_air(temperature: float, pressure_atm: float = 1.0, xi: float = DEFAULT_XI) -> dict[str, float]:
    z, eps1, eps2, eps3, x_n2, x_o2, x_n, x_o, x_op, x_np, xe = hansen_air_cached(
        float(temperature), float(pressure_atm), float(xi)
    )
    return {
        "T": temperature,
        "p_atm": pressure_atm,
        "xi": xi,
        "Z": z,
        "eps1": eps1,
        "eps2": eps2,
        "eps3": eps3,
        "xN2": x_n2,
        "xO2": x_o2,
        "xN": x_n,
        "xO": x_o,
        "xO+": x_op,
        "xN+": x_np,
        "xe-": xe,
    }


def dln_qp_dt(species: str, temperature: float) -> float:
    step = 1e-3 * temperature
    return (ln_qp(species, temperature + step) - ln_qp(species, temperature - step)) / (2.0 * step)


@lru_cache(maxsize=20000)
def species_h_over_rt(species: str, temperature: float) -> float:
    return temperature * dln_qp_dt(species, temperature) + FORMATION_OVER_R[species] / temperature


def mixture_zh_over_rt(temperature: float, pressure_atm: float, xi: float) -> float:
    state = hansen_air(temperature, pressure_atm, xi)
    return state["Z"] * sum(
        state["xe-" if sp == "e-" else f"x{sp}"] * species_h_over_rt(sp, temperature)
        for sp in SPECIES
    )


def mixture_ze_over_rt(temperature: float, pressure_atm: float, xi: float) -> float:
    state = hansen_air(temperature, pressure_atm, xi)
    return state["Z"] * sum(
        state["xe-" if sp == "e-" else f"x{sp}"] * (species_h_over_rt(sp, temperature) - 1.0)
        for sp in SPECIES
    )


def pressure_for_constant_density(new_temperature: float, ref_temperature: float, ref_pressure_atm: float, xi: float) -> float:
    z_ref = hansen_air(ref_temperature, ref_pressure_atm, xi)["Z"]
    target = ref_pressure_atm / (z_ref * ref_temperature)
    pressure = ref_pressure_atm * new_temperature / ref_temperature
    for _ in range(30):
        z = hansen_air(new_temperature, pressure, xi)["Z"]
        f = pressure / (z * new_temperature) - target
        dp = 1e-5 * max(pressure, 1e-8)
        zp = hansen_air(new_temperature, pressure + dp, xi)["Z"]
        fp = (pressure + dp) / (zp * new_temperature) - target
        dfdp = (fp - f) / dp
        if dfdp == 0.0 or not math.isfinite(dfdp):
            break
        pressure = max(pressure - f / dfdp, 1e-12)
    return pressure


def hansen_cp_over_r(temperature: float, pressure_atm: float, xi: float) -> float:
    step = 1e-3 * temperature
    yp = mixture_zh_over_rt(temperature + step, pressure_atm, xi)
    ym = mixture_zh_over_rt(temperature - step, pressure_atm, xi)
    return ((temperature + step) * yp - (temperature - step) * ym) / (2.0 * step)


def hansen_cv_over_r(temperature: float, pressure_atm: float, xi: float) -> float:
    step = 1e-3 * temperature
    pp = pressure_for_constant_density(temperature + step, temperature, pressure_atm, xi)
    pm = pressure_for_constant_density(temperature - step, temperature, pressure_atm, xi)
    xp = mixture_ze_over_rt(temperature + step, pp, xi)
    xm = mixture_ze_over_rt(temperature - step, pm, xi)
    return ((temperature + step) * xp - (temperature - step) * xm) / (2.0 * step)


def hansen_phi_factor(temperature: float, pressure_atm: float, xi: float) -> float:
    z0 = hansen_air(temperature, pressure_atm, xi)["Z"]
    rho0_scaled = pressure_atm / (z0 * temperature)
    dp = 1e-4 * max(pressure_atm, 1e-8)
    p1 = pressure_atm + dp
    p2 = max(pressure_atm - dp, 1e-12)
    rho1_scaled = p1 / (hansen_air(temperature, p1, xi)["Z"] * temperature)
    rho2_scaled = p2 / (hansen_air(temperature, p2, xi)["Z"] * temperature)
    return (rho0_scaled / pressure_atm) * ((p1 - p2) / (rho1_scaled - rho2_scaled))


def mole_fractions(temperature: float, pressure_pa: float, xi: float) -> list[float]:
    x_n2_feed, x_o2_feed = feed_fractions(xi)
    if temperature <= 1200.0:
        return [x_n2_feed, x_o2_feed, 0.0, 0.0, 0.0, 0.0, 0.0]
    state = hansen_air(temperature, pressure_pa / ATM_TO_PA, xi)
    values = [state["xN2"], state["xO2"], state["xN"], state["xO"], state["xO+"], state["xN+"], state["xe-"]]
    total = sum(values)
    if total <= 0.0 or not math.isfinite(total):
        raise SolverError("Invalid reacting-air species fractions.")
    return [max(0.0, value) / total for value in values]


def mass_fractions_from_mole_fractions(x: list[float]) -> list[float]:
    mixture_mw = sum(value * MOLAR_MASS[index] for index, value in enumerate(x))
    if mixture_mw <= 0.0 or not math.isfinite(mixture_mw):
        raise SolverError("Invalid reacting-air mixture molecular weight.")
    return [value * MOLAR_MASS[index] / mixture_mw for index, value in enumerate(x)]


def mixture_gas_constant_from_mass_fractions(y: list[float]) -> float:
    return sum(value * SPECIFIC_R[index] for index, value in enumerate(y))


def partial_pressure_composition(pressure_pa: float, temperature: float, xi: float) -> tuple[list[float], list[float]]:
    x = mole_fractions(temperature, pressure_pa, xi)
    return x, mass_fractions_from_mole_fractions(x)


def frozen_mixture_enthalpy(temperature: float, mass_fractions: list[float]) -> float:
    vib_n2 = 3389.82
    vib_o2 = 2771.09
    return (
        mass_fractions[0] * (3.5 * SPECIFIC_R[0] * temperature + vib_n2 * SPECIFIC_R[0] / (safe_exp(vib_n2 / temperature) - 1.0))
        + mass_fractions[1] * (3.5 * SPECIFIC_R[1] * temperature + vib_o2 * SPECIFIC_R[1] / (safe_exp(vib_o2 / temperature) - 1.0))
        + sum(mass_fractions[i] * 2.5 * SPECIFIC_R[i] * temperature for i in range(2, 7))
    )


def equilibrium_enthalpy(temperature: float, pressure_atm: float, xi: float) -> float:
    _, y = partial_pressure_composition(pressure_atm * ATM_TO_PA, temperature, xi)
    if temperature <= 1200.0:
        return frozen_mixture_enthalpy(temperature, y)
    return mixture_gas_constant(xi) * temperature * mixture_zh_over_rt(temperature, pressure_atm, xi)


def equilibrium_state(temperature: float, pressure_pa: float, xi: float) -> dict[str, object]:
    x, y = partial_pressure_composition(pressure_pa, temperature, xi)
    gas_constant = mixture_gas_constant_from_mass_fractions(y)
    density = pressure_pa / (gas_constant * temperature)
    pressure_atm = pressure_pa / ATM_TO_PA
    phi_gamma = (
        1.4
        if temperature <= 1200.0
        else hansen_cp_over_r(temperature, pressure_atm, xi)
        / hansen_cv_over_r(temperature, pressure_atm, xi)
        * hansen_phi_factor(temperature, pressure_atm, xi)
    )
    sound_speed = math.sqrt(phi_gamma * pressure_pa / density)
    return {
        "mole_fractions": x,
        "mass_fractions": y,
        "gas_constant": gas_constant,
        "density": density,
        "phi_gamma": phi_gamma,
        "sound_speed": sound_speed,
    }


def normal_shock_full(mach1: float, temperature1: float, pressure1_pa: float, xi: float = DEFAULT_XI) -> dict[str, object]:
    if mach1 <= 1.0:
        raise SolverError("Normal-shock M1 must be greater than 1.")
    upstream = equilibrium_state(temperature1, pressure1_pa, xi)
    density1 = upstream["density"]
    velocity_normal = mach1 * upstream["sound_speed"]
    gamma_guess = max(upstream["phi_gamma"], 1.05)
    density_ratio_ideal = ((gamma_guess + 1.0) * mach1 * mach1) / ((gamma_guess - 1.0) * mach1 * mach1 + 2.0)
    pressure_ratio_ideal = 1.0 + (2.0 * gamma_guess / (gamma_guess + 1.0)) * (mach1 * mach1 - 1.0)
    temperature_guess = temperature1 * pressure_ratio_ideal / density_ratio_ideal
    h1 = equilibrium_enthalpy(temperature1, pressure1_pa / ATM_TO_PA, xi)
    p_ref = max(pressure1_pa, 1.0)
    h_ref = max(abs(h1) + 0.5 * velocity_normal * velocity_normal, 1.0)
    t_min = max(1.0, min(200.0, 0.75 * temperature1))
    t_max = max(60000.0, 8.0 * temperature1, 1.5 * temperature_guess)
    bounds = (
        math.log(t_min),
        math.log(t_max),
        math.log(density1 * 1.0001),
        math.log(density1 * 100.0),
    )
    x = [
        math.log(min(max(temperature_guess, t_min * 1.01), t_max * 0.98)),
        math.log(max(density1 * density_ratio_ideal, density1 * 1.01)),
    ]

    def clamp(state: list[float]) -> list[float]:
        return [
            min(max(state[0], bounds[0]), bounds[1]),
            min(max(state[1], bounds[2]), bounds[3]),
        ]

    def residuals(state: list[float]) -> tuple[float, float]:
        temperature2 = math.exp(state[0])
        density2 = math.exp(state[1])
        velocity2 = velocity_normal * density1 / density2
        pressure2 = pressure1_pa + density1 * velocity_normal * velocity_normal * (1.0 - density1 / density2)
        if pressure2 <= 0.0 or not math.isfinite(pressure2):
            return 1e6, 1e6
        _, y = partial_pressure_composition(pressure2, temperature2, xi)
        gas_constant = mixture_gas_constant_from_mass_fractions(y)
        eos_residual = (pressure2 - density2 * gas_constant * temperature2) / p_ref
        target_h2 = h1 + 0.5 * (velocity_normal * velocity_normal - velocity2 * velocity2)
        h2 = equilibrium_enthalpy(temperature2, pressure2 / ATM_TO_PA, xi)
        enthalpy_residual = (h2 - target_h2) / h_ref
        if not math.isfinite(eos_residual) or not math.isfinite(enthalpy_residual):
            return 1e6, 1e6
        return eos_residual, enthalpy_residual

    def norm(values: tuple[float, float]) -> float:
        return math.hypot(values[0], values[1])

    for _ in range(80):
        f = residuals(x)
        if norm(f) < 1e-9:
            break
        jac = [[0.0, 0.0], [0.0, 0.0]]
        for column in range(2):
            step = 1e-5
            xp = x[:]
            xp[column] += step
            fp = residuals(clamp(xp))
            jac[0][column] = (fp[0] - f[0]) / step
            jac[1][column] = (fp[1] - f[1]) / step
        det = jac[0][0] * jac[1][1] - jac[0][1] * jac[1][0]
        if abs(det) < 1e-18 or not math.isfinite(det):
            break
        dx = [
            (jac[1][1] * f[0] - jac[0][1] * f[1]) / det,
            (-jac[1][0] * f[0] + jac[0][0] * f[1]) / det,
        ]
        current_norm = norm(f)
        damping = 1.0
        for _trial in range(12):
            candidate = clamp([x[0] - damping * dx[0], x[1] - damping * dx[1]])
            if norm(residuals(candidate)) < current_norm:
                x = candidate
                break
            damping *= 0.5
        else:
            break

    final = residuals(x)
    if norm(final) > 1e-5:
        raise SolverError("Reacting normal shock solve did not converge for these inputs.")
    temperature2 = math.exp(x[0])
    density2 = math.exp(x[1])
    velocity2 = velocity_normal * density1 / density2
    pressure2 = pressure1_pa + density1 * velocity_normal * velocity_normal * (1.0 - density1 / density2)
    mole, mass = partial_pressure_composition(pressure2, temperature2, xi)
    downstream = equilibrium_state(temperature2, pressure2, xi)
    return {
        "mach1": mach1,
        "mach2": velocity2 / downstream["sound_speed"],
        "temperature1": temperature1,
        "pressure1": pressure1_pa,
        "density1": density1,
        "velocity1_normal": velocity_normal,
        "temperature2": temperature2,
        "pressure2": pressure2,
        "density2": density2,
        "velocity2_normal": velocity2,
        "mole_fractions": mole,
        "mass_fractions": mass,
    }


def reacting_cone_surface_at_beta(
    gamma: float,
    mach1: float,
    beta: float,
    cone_angle: float,
    normal_shock: dict[str, object],
) -> dict[str, float]:
    velocity1 = nondimensional_velocity(gamma, mach1)
    density_ratio = normal_shock["density2"] / normal_shock["density1"]
    radial_velocity = velocity1 * math.cos(beta)
    polar_velocity = -velocity1 * math.sin(beta) / density_ratio
    surface_radial, surface_polar = integrate_taylor_maccoll(
        gamma, beta, cone_angle, radial_velocity, polar_velocity
    )
    return {
        "radial_velocity": radial_velocity,
        "polar_velocity": polar_velocity,
        "surface_radial_velocity": surface_radial,
        "surface_polar_velocity": surface_polar,
    }


def calculate_reacting_cone(
    mach1: float,
    temperature1: float,
    pressure1_pa: float,
    cone_angle_degrees: float,
    xi: float = DEFAULT_XI,
) -> dict[str, object]:
    if mach1 <= 1.0:
        raise SolverError("M1 must be greater than 1.")
    if cone_angle_degrees <= 0.0 or cone_angle_degrees >= 89.0:
        raise SolverError("Cone half-angle must be between 0 and 89 degrees.")
    if temperature1 <= 0.0 or pressure1_pa <= 0.0:
        raise SolverError("T1 and p1 must be positive.")

    cone_angle = cone_angle_degrees * DEG
    upstream = equilibrium_state(temperature1, pressure1_pa, xi)
    gamma = min(max(upstream["phi_gamma"], 1.05), 1.67)
    velocity1 = mach1 * upstream["sound_speed"]
    beta_min = max(math.asin(1.0 / mach1), cone_angle) + 1e-5
    beta_max = math.pi / 2.0 - 1e-5
    shock_cache: dict[str, dict[str, object]] = {}

    def shock_for_beta(beta: float) -> dict[str, object]:
        key = f"{beta:.14g}"
        if key not in shock_cache:
            shock_cache[key] = normal_shock_full(mach1 * math.sin(beta), temperature1, pressure1_pa, xi)
        return shock_cache[key]

    def residual(beta: float) -> float:
        try:
            normal_mach = mach1 * math.sin(beta)
            if normal_mach <= 1.0:
                return math.nan
            return reacting_cone_surface_at_beta(
                gamma, mach1, beta, cone_angle, shock_for_beta(beta)
            )["surface_polar_velocity"]
        except SolverError:
            return math.nan

    previous_beta = None
    previous_residual = None
    bracket = None
    for i in range(141):
        beta = beta_min + (i / 140.0) * (beta_max - beta_min)
        value = residual(beta)
        if (
            previous_residual is not None
            and math.isfinite(previous_residual)
            and math.isfinite(value)
            and previous_residual * value <= 0.0
        ):
            bracket = (previous_beta, beta)
            break
        if math.isfinite(value):
            previous_beta = beta
            previous_residual = value

    if bracket is None:
        raise SolverError("No attached reacting cone shock found for this M1 and cone angle.")
    beta = solve_bisection(residual, bracket[0], bracket[1], 1e-9, 80)
    normal_mach = mach1 * math.sin(beta)
    normal_shock = shock_for_beta(beta)
    cone_velocity = reacting_cone_surface_at_beta(gamma, mach1, beta, cone_angle, normal_shock)
    tangential_velocity = velocity1 * math.cos(beta)
    velocity_behind_shock = math.hypot(tangential_velocity, normal_shock["velocity2_normal"])
    downstream = equilibrium_state(normal_shock["temperature2"], normal_shock["pressure2"], xi)
    mach_behind_shock = velocity_behind_shock / downstream["sound_speed"]
    shock_velocity = math.hypot(cone_velocity["radial_velocity"], cone_velocity["polar_velocity"])
    surface_velocity = math.hypot(cone_velocity["surface_radial_velocity"], cone_velocity["surface_polar_velocity"])
    surface_mach = mach_from_nondimensional_velocity(gamma, surface_velocity)
    tc_t2 = (1.0 - surface_velocity * surface_velocity) / (1.0 - shock_velocity * shock_velocity)
    pc_p2 = tc_t2 ** (gamma / (gamma - 1.0))
    surface_temperature = normal_shock["temperature2"] * tc_t2
    surface_pressure = normal_shock["pressure2"] * pc_p2
    surface = equilibrium_state(surface_temperature, surface_pressure, xi)

    shock_species_density = {
        sp: normal_shock["mass_fractions"][index] * normal_shock["density2"]
        for index, sp in enumerate(SPECIES)
    }
    surface_species_density = {
        sp: surface["mass_fractions"][index] * surface["density"]
        for index, sp in enumerate(SPECIES)
    }
    return {
        "M1": mach1,
        "xi": xi,
        "gamma_tm": gamma,
        "cone_half_angle_deg": cone_angle_degrees,
        "shock_angle_deg": beta * RAD,
        "Mn1": normal_mach,
        "M2_behind_shock": mach_behind_shock,
        "T2": normal_shock["temperature2"],
        "p2": normal_shock["pressure2"],
        "rho2": normal_shock["density2"],
        "V2": velocity_behind_shock,
        "Mc_cone_surface": surface_mach,
        "Tc": surface_temperature,
        "pc": surface_pressure,
        "rhoc": surface["density"],
        "mole_fractions_behind_shock": dict(zip(SPECIES, normal_shock["mole_fractions"])),
        "mass_fractions_behind_shock": dict(zip(SPECIES, normal_shock["mass_fractions"])),
        "species_density_behind_shock": shock_species_density,
        "mole_fractions_cone_surface": dict(zip(SPECIES, surface["mole_fractions"])),
        "mass_fractions_cone_surface": dict(zip(SPECIES, surface["mass_fractions"])),
        "species_density_cone_surface": surface_species_density,
    }


def add_hansen_states(perfect: dict[str, float], xi: float) -> dict[str, object]:
    if not {"p1", "T1", "p2", "T2", "pc", "Tc"}.issubset(perfect):
        raise SolverError("Hansen state output for a perfect cone requires p1 and T1.")
    return {
        **perfect,
        "hansen_freestream": hansen_air(perfect["T1"], perfect["p1"] / ATM_TO_PA, xi),
        "hansen_behind_shock": hansen_air(perfect["T2"], perfect["p2"] / ATM_TO_PA, xi),
        "hansen_cone_surface": hansen_air(perfect["Tc"], perfect["pc"] / ATM_TO_PA, xi),
    }


def print_table(values: dict[str, object]) -> None:
    flat = {key: value for key, value in values.items() if isinstance(value, (int, float))}
    width = max(len(key) for key in flat) if flat else 0
    for key, value in flat.items():
        print(f"{key:<{width}} : {value:.10g}")
    nested = {key: value for key, value in values.items() if isinstance(value, dict)}
    for name, mapping in nested.items():
        print(f"\n{name}:")
        inner_width = max(len(key) for key in mapping)
        for key, value in mapping.items():
            if isinstance(value, (int, float)):
                print(f"  {key:<{inner_width}} : {value:.10g}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    air = subparsers.add_parser("air", help="Evaluate Hansen equilibrium air.")
    air.add_argument("--T", type=float, required=True, help="Temperature in K.")
    air.add_argument("--p-atm", type=float, help="Pressure in atm.")
    air.add_argument("--p-pa", type=float, help="Pressure in Pa.")
    air.add_argument("--xi", type=float, default=DEFAULT_XI, help="O/N atomic ratio, default air = 0.25.")
    air.add_argument("--json", action="store_true", help="Print JSON instead of a table.")

    cone = subparsers.add_parser("cone", help="Perfect-gas conical shock, optionally with Hansen states.")
    cone.add_argument("--mach", "-M", type=float, required=True)
    cone.add_argument("--cone-angle", "-c", type=float, required=True)
    cone.add_argument("--gamma", "-g", type=float, default=1.4)
    cone.add_argument("--p1", type=float, help="Freestream pressure in Pa.")
    cone.add_argument("--T1", type=float, help="Freestream temperature in K.")
    cone.add_argument("--R", type=float, default=287.05)
    cone.add_argument("--hansen", action="store_true", help="Append Hansen states at p1/T1, p2/T2, and pc/Tc.")
    cone.add_argument("--xi", type=float, default=DEFAULT_XI)
    cone.add_argument("--json", action="store_true", help="Print JSON instead of a table.")

    reacting = subparsers.add_parser("reacting-cone", help="Reacting Hansen conical shock.")
    reacting.add_argument("--mach", "-M", type=float, required=True)
    reacting.add_argument("--cone-angle", "-c", type=float, required=True)
    reacting.add_argument("--T1", type=float, required=True)
    reacting.add_argument("--p1", type=float, required=True, help="Freestream pressure in Pa.")
    reacting.add_argument("--xi", type=float, default=DEFAULT_XI)
    reacting.add_argument("--json", action="store_true", help="Print JSON instead of a table.")
    args = parser.parse_args()

    try:
        if args.command == "air":
            pressure_atm = args.p_atm if args.p_atm is not None else (args.p_pa / ATM_TO_PA if args.p_pa else 1.0)
            result = hansen_air(args.T, pressure_atm, args.xi)
        elif args.command == "cone":
            result = calculate_perfect_cone(args.mach, args.cone_angle, args.gamma, args.p1, args.T1, args.R)
            if args.hansen:
                result = add_hansen_states(result, args.xi)
        else:
            result = calculate_reacting_cone(args.mach, args.T1, args.p1, args.cone_angle, args.xi)
    except SolverError as exc:
        parser.error(str(exc))

    if args.json:
        print(json.dumps(result, indent=2))
    else:
        print_table(result)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
