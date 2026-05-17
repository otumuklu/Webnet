import numpy as np
from scipy.optimize import root
from scipy.optimize import root_scalar, least_squares


def safe_exp(x):
    return np.exp(np.clip(x, -700.0, 700.0))


SPECIES_NAMES = ["N2", "O2", "N", "O", "O+", "N+", "e-"]
SPECIES_LABELS = [r"$N_2$", r"$O_2$", r"$N$", r"$O$", r"$O^+$", r"$N^+$", r"$e^-$"]
SPECIES_COLORS = ["blue", "red", "green", "purple", "orange", "black", "cyan"]
ATM_TO_PA = 101325.0
R_UNIVERSAL = 8.31446261815324

# Species order: [N2, O2, N, O, O+, N+, e-]
MOLAR_MASS = np.array([28.02, 32.00, 14.01, 16.00, 16.00, 14.01, 5.485799e-4], dtype=float)
R_VALS = np.array([296.72, 259.81, 593.43, 519.63, 519.63, 593.43, 15156338.18], dtype=float)

E0_OVER_R = {
    "N2": 0.0,
    "O2": 0.0,
    "O": 59000.0 / 2.0,
    "N": 113200.0 / 2.0,
    "O+": 59000.0 / 2.0 + 158000.0,
    "N+": 113200.0 / 2.0 + 168800.0,
    "e-": 0.0,
}


# ============================================================
# Hansen / Axioms high-temperature equilibrium air model
# Pressure is in atm, temperature is in K.
# Valid mainly for high-temperature air up to ~15000 K.
# ============================================================

def lnQp(species, T):
    """
    Pressure-standardized partition function ln(Qp).
    In Hansen notation, ln(Q) expressions include -ln(p).
    Therefore ln(Qp) = ln(Q) + ln(p), so the pressure term cancels.
    """
    T = np.asarray(T, dtype=float)

    if species == "N2":
        return 3.5*np.log(T) - 0.42 - np.log(1.0 - np.exp(-3390.0/T))

    if species == "O2":
        return (3.5*np.log(T) + 0.11
                - np.log(1.0 - np.exp(-2270.0/T))
                + np.log(3.0 + 2.0*np.exp(-11390.0/T)
                         + np.exp(-18990.0/T)))

    if species == "O":
        return (2.5*np.log(T) + 0.50
                + np.log(5.0 + 3.0*np.exp(-228.0/T)
                         + np.exp(-326.0/T)
                         + 5.0*np.exp(-22800.0/T)
                         + np.exp(-48600.0/T)))

    if species == "N":
        return (2.5*np.log(T) + 0.30
                + np.log(4.0 + 10.0*np.exp(-27700.0/T)
                         + 6.0*np.exp(-41500.0/T)))

    if species == "O+":
        return (2.5*np.log(T) + 0.50
                + np.log(4.0 + 10.0*np.exp(-38600.0/T)
                         + 6.0*np.exp(-58200.0/T)))

    if species == "N+":
        return (2.5*np.log(T) + 0.30
                + np.log(1.0 + 3.0*np.exp(-70.6/T)
                         + 5.0*np.exp(-188.9/T)
                         + 5.0*np.exp(-22000.0/T)
                         + np.exp(-47000.0/T)
                         + 5.0*np.exp(-67900.0/T)))

    if species == "e-":
        return 2.5*np.log(T) - 14.24

    raise ValueError(f"Unknown species: {species}")


def equilibrium_constants(T):
    """
    Hansen pressure equilibrium constants Kp.
    Units are consistent with pressure in atm.
    """
    T = np.asarray(T, dtype=float)

    lnKp1 = -59000.0/T + 2.0*lnQp("O", T) - lnQp("O2", T)
    lnKp2 = -113200.0/T + 2.0*lnQp("N", T) - lnQp("N2", T)

    lnKpOion = -158000.0/T + lnQp("O+", T) + lnQp("e-", T) - lnQp("O", T)
    lnKpNion = -168800.0/T + lnQp("N+", T) + lnQp("e-", T) - lnQp("N", T)

    Kp1 = np.exp(lnKp1)
    Kp2 = np.exp(lnKp2)
    KpOion = np.exp(lnKpOion)
    KpNion = np.exp(lnKpNion)

    # Hansen population-weighted ionization constant
    Kp3 = 0.2*KpOion + 0.8*KpNion

    return Kp1, Kp2, Kp3, KpOion, KpNion


def epsilons(T, p_atm):
    """
    Computes Hansen epsilon parameters:
    eps1: O2 dissociation parameter
    eps2: N2 dissociation parameter
    eps3: ionization parameter
    """
    Kp1, Kp2, Kp3, _, _ = equilibrium_constants(T)
    p = p_atm

    A1 = 1.0 + 4.0*p/Kp1
    eps1 = (-0.8 + np.sqrt(0.64 + 0.8*A1)) / (2.0*A1)

    A2 = 1.0 + 4.0*p/Kp2
    eps2 = (-0.4 + np.sqrt(0.16 + 3.84*A2)) / (2.0*A2)

    eps3 = (1.0 + p/Kp3)**(-0.5)

    return eps1, eps2, eps3


def feed_molecular_fractions_from_xi(xi):
    """
    xi is the O/N atomic ratio. Air is xi=0.25, giving N2/O2 = 0.8/0.2.
    """
    xi = float(xi)
    if xi <= 0.0 or not np.isfinite(xi):
        raise ValueError("xi must be a positive finite O/N atomic ratio.")

    xN2_feed = 1.0 / (1.0 + xi)
    xO2_feed = xi / (1.0 + xi)

    return xN2_feed, xO2_feed


def mixture_gas_constant_from_feed(xi):
    xN2_feed, xO2_feed = feed_molecular_fractions_from_xi(xi)
    M_feed_kg = (xN2_feed * MOLAR_MASS[0] + xO2_feed * MOLAR_MASS[1]) / 1000.0
    return R_UNIVERSAL / M_feed_kg


def hansen_air(T, p_atm=1.0, xi=0.25):
    """
    Returns compressibility factor and species mole fractions.
    The Hansen constants are exact for air (xi=0.25). For other xi values,
    the same equilibrium model is scaled by the requested N2/O2 feed ratio.
    """
    eps1_air, eps2_air, eps3 = epsilons(T, p_atm)
    xN2_feed, xO2_feed = feed_molecular_fractions_from_xi(xi)

    eps1 = np.minimum(eps1_air * xO2_feed / 0.2, xO2_feed)
    eps2 = np.minimum(eps2_air * xN2_feed / 0.8, xN2_feed)

    N_atom_share = xN2_feed / (xN2_feed + xO2_feed)
    O_atom_share = xO2_feed / (xN2_feed + xO2_feed)

    Z = 1.0 + eps1 + eps2 + 2.0*eps3

    xO2 = (xO2_feed - eps1) / Z
    xN2 = (xN2_feed - eps2) / Z
    xO  = (2.0*eps1 - 2.0*O_atom_share*eps3) / Z
    xN  = (2.0*eps2 - 2.0*N_atom_share*eps3) / Z

    # Axioms paper gives x(O+ + N+) = x(e-) = 2 eps3 / Z.
    # Split ions by the requested oxygen/nitrogen feed ratio.
    xe  = 2.0*eps3 / Z
    xOp = O_atom_share*xe
    xNp = N_atom_share*xe

    return {
        "T": T,
        "p_atm": p_atm,
        "xi": xi,
        "eps1": eps1,
        "eps2": eps2,
        "eps3": eps3,
        "Z": Z,
        "xO2": xO2,
        "xN2": xN2,
        "xO": xO,
        "xN": xN,
        "xO+": xOp,
        "xN+": xNp,
        "xe-": xe,
    }


def dlnQp_dT(species, T):
    hT = 1e-3 * T
    return (lnQp(species, T + hT) - lnQp(species, T - hT)) / (2.0*hT)


def species_H_RT(species, T):
    return T*dlnQp_dT(species, T) + E0_OVER_R[species]/T


def species_E_RT(species, T):
    return species_H_RT(species, T) - 1.0


def mixture_ZH_RT(T, p_atm, xi=0.25):
    sol = hansen_air(T, p_atm, xi)
    return sol["Z"] * (
        sol["xO2"] * species_H_RT("O2", T) +
        sol["xN2"] * species_H_RT("N2", T) +
        sol["xO"]  * species_H_RT("O", T) +
        sol["xN"]  * species_H_RT("N", T) +
        sol["xO+"] * species_H_RT("O+", T) +
        sol["xN+"] * species_H_RT("N+", T) +
        sol["xe-"] * species_H_RT("e-", T)
    )


def mixture_ZE_RT(T, p_atm, xi=0.25):
    sol = hansen_air(T, p_atm, xi)
    return sol["Z"] * (
        sol["xO2"] * species_E_RT("O2", T) +
        sol["xN2"] * species_E_RT("N2", T) +
        sol["xO"]  * species_E_RT("O", T) +
        sol["xN"]  * species_E_RT("N", T) +
        sol["xO+"] * species_E_RT("O+", T) +
        sol["xN+"] * species_E_RT("N+", T) +
        sol["xe-"] * species_E_RT("e-", T)
    )


def pressure_for_constant_density(T_new, T_ref, p_ref, xi=0.25):
    Z_ref = hansen_air(T_ref, p_ref, xi)["Z"]
    target = p_ref / (Z_ref*T_ref)
    p_iter = p_ref * T_new/T_ref

    for _ in range(30):
        Z = hansen_air(T_new, p_iter, xi)["Z"]
        f = p_iter/(Z*T_new) - target

        dp = 1e-5 * max(p_iter, 1e-8)
        Zp = hansen_air(T_new, p_iter + dp, xi)["Z"]
        fp = (p_iter + dp)/(Zp*T_new) - target

        dfdp = (fp - f)/dp
        if dfdp == 0.0 or not np.isfinite(dfdp):
            break
        p_iter -= f/dfdp
        p_iter = max(p_iter, 1e-12)

    return p_iter


def cp_over_R(T, p_atm, xi=0.25):
    hT = 1e-3*T
    Yp = mixture_ZH_RT(T + hT, p_atm, xi)
    Ym = mixture_ZH_RT(T - hT, p_atm, xi)
    return ((T + hT)*Yp - (T - hT)*Ym) / (2.0*hT)


def cv_over_R(T, p_atm, xi=0.25):
    hT = 1e-3*T
    pp = pressure_for_constant_density(T + hT, T, p_atm, xi)
    pm = pressure_for_constant_density(T - hT, T, p_atm, xi)
    Xp = mixture_ZE_RT(T + hT, pp, xi)
    Xm = mixture_ZE_RT(T - hT, pm, xi)
    return ((T + hT)*Xp - (T - hT)*Xm) / (2.0*hT)


def phi_factor(T, p_atm, xi=0.25):
    Z0 = hansen_air(T, p_atm, xi)["Z"]
    rho0_scaled = p_atm/(Z0*T)

    dp = 1e-4 * max(p_atm, 1e-8)
    p1 = p_atm + dp
    p2 = max(p_atm - dp, 1e-12)

    Z1 = hansen_air(T, p1, xi)["Z"]
    Z2 = hansen_air(T, p2, xi)["Z"]

    rho1_scaled = p1/(Z1*T)
    rho2_scaled = p2/(Z2*T)
    dp_drho = (p1 - p2)/(rho1_scaled - rho2_scaled)

    return rho0_scaled/p_atm * dp_drho


def speed_sound_parameter(T, p_atm, xi=0.25):
    if T <= 1200.0:
        return 1.4

    Cp = cp_over_R(T, p_atm, xi)
    Cv = cv_over_R(T, p_atm, xi)
    gamma_eq = Cp/Cv
    Phi = phi_factor(T, p_atm, xi)
    return gamma_eq*Phi


def hansen_mixture_enthalpy(T, p_atm, xi=0.25):
    return mixture_gas_constant_from_feed(xi) * T * mixture_ZH_RT(T, p_atm, xi)


def equilibrium_enthalpy(T, p_atm, xi=0.25):
    p_local = p_atm * ATM_TO_PA
    c_eq, _ = partial_pressure(p_local, T, xi)

    if T <= 1200.0:
        return mixture_enthalpy(T, c_eq)

    return hansen_mixture_enthalpy(T, p_atm, xi)


def partial_pressure(p, T, xi):
    """
    Existing solver interface backed by the Hansen air model.
    Returns c and X in the species order [N2, O2, N, O, O+, N+, e-].
    """
    p_atm = p / ATM_TO_PA
    xN2_feed, xO2_feed = feed_molecular_fractions_from_xi(xi)

    if T <= 1200:
        X = np.array([xN2_feed, xO2_feed, 0.0, 0.0, 0.0, 0.0, 0.0], dtype=float)
    else:
        sol = hansen_air(T, p_atm, xi)
        X = np.array([
            sol["xN2"],
            sol["xO2"],
            sol["xN"],
            sol["xO"],
            sol["xO+"],
            sol["xN+"],
            sol["xe-"]
        ], dtype=float)
        X = np.maximum(X, 0.0)
        X_sum = np.sum(X)
        if X_sum <= 0.0 or not np.isfinite(X_sum):
            raise ValueError(f"Invalid Hansen mole fractions: {X, p, T}")
        X /= X_sum

    M_total = float(np.dot(X, MOLAR_MASS))

    if M_total == 0:
        raise ValueError(f"Current Values = {X, p, T}")

    c = X * (MOLAR_MASS / M_total)

    return c, X


# ==========================================================
# mixture_enthalpy
# ==========================================================
# Computes sensible enthalpy of mixture (no formation).
# ==========================================================

def mixture_enthalpy(T, c):
    c = np.array(c, dtype=float)

    V_N2 = 3389.82
    V_O2 = 2771.09

    h = (
        c[0] * (7.0 / 2.0 * R_VALS[0] * T + V_N2 * R_VALS[0] / (np.exp(V_N2 / T) - 1.0)) +
        c[1] * (7.0 / 2.0 * R_VALS[1] * T + V_O2 * R_VALS[1] / (np.exp(V_O2 / T) - 1.0)) +
        c[2] * (5.0 / 2.0 * R_VALS[2] * T) +
        c[3] * (5.0 / 2.0 * R_VALS[3] * T) +
        c[4] * (5.0 / 2.0 * R_VALS[4] * T) +
        c[5] * (5.0 / 2.0 * R_VALS[5] * T) +
        c[6] * (5.0 / 2.0 * R_VALS[6] * T)
    )

    return h


"""-----------------------------------------------------------------------------------------------------------"""

# Freestream conditions
T = 227.13
p = 1090.16
rho = 0.01672
a = 302.122
r = 0.1
h = equilibrium_enthalpy(T, p / ATM_TO_PA, 0.25)
gamma = 1.4
V = 2416


# ==========================================================
# shockwave_calculator
# ==========================================================
# Solves normal shock with equilibrium chemistry.
# Outer loop: density ratio iteration
# Inner loop: temperature solve via enthalpy balance
# ==========================================================

def normal_shock_full(Vn, xi=0.25):
    """
    Returns:
        T2, rho2, p2, V2, X_vals, c_vals
    """
    M = Vn / a
    rho2_rho1_ideal = (gamma + 1.0) * M**2 / ((gamma - 1.0) * M**2 + 2.0)
    p2_p1_ideal = 1.0 + 2.0 * gamma / (gamma + 1.0) * (M**2 - 1.0)
    T2_T1_ideal = p2_p1_ideal / rho2_rho1_ideal

    T2_guess = float(np.clip(T * T2_T1_ideal, 300.0, 40000.0))
    rho2_guess = float(max(rho * rho2_rho1_ideal, rho * 1.01))

    x0 = np.array([np.log(T2_guess), np.log(rho2_guess)], dtype=float)

    T2_min, T2_max = 200.0, 60000.0
    rho2_min = rho * 1.0001
    rho2_max = rho * 100.0

    lb = np.array([np.log(T2_min), np.log(rho2_min)], dtype=float)
    ub = np.array([np.log(T2_max), np.log(rho2_max)], dtype=float)

    p_ref = max(p, 1.0)
    h_ref = max(abs(h) + 0.5 * Vn**2, 1.0)

    def residuals(x):
        T2 = float(np.exp(x[0]))
        rho2 = float(np.exp(x[1]))

        V2 = Vn * rho / rho2
        p2 = p + rho * Vn**2 * (1.0 - rho / rho2)

        if (not np.isfinite(p2)) or (p2 <= 0.0):
            return np.array([1e6, 1e6], dtype=float)

        c_vals, X_vals = partial_pressure(p2, T2, xi)
        Rmix = float(np.dot(R_VALS, c_vals))

        r1 = (p2 - rho2 * Rmix * T2) / p_ref

        h2_target = h + 0.5 * (Vn**2 - V2**2)

        h2 = equilibrium_enthalpy(T2, p2 / ATM_TO_PA, xi)
        r2 = (h2 - h2_target) / h_ref

        if (not np.isfinite(r1)) or (not np.isfinite(r2)):
            return np.array([1e6, 1e6], dtype=float)

        return np.array([r1, r2], dtype=float)

    sol = least_squares(
        residuals,
        x0=x0,
        bounds=(lb, ub),
        xtol=1e-12,
        ftol=1e-12,
        gtol=1e-12,
        max_nfev=500
    )

    if not sol.success:
        T2_try, rho2_try = float(np.exp(sol.x[0])), float(np.exp(sol.x[1]))
        res = residuals(sol.x)
        raise RuntimeError(
            f"Normal shock solve failed: {sol.message}\n"
            f"  guess: T2={T2_guess:.3f}, rho2={rho2_guess:.6e}\n"
            f"  last : T2={T2_try:.3f}, rho2={rho2_try:.6e}\n"
            f"  res  : [{res[0]:.6e}, {res[1]:.6e}]"
        )

    T2_sol = float(np.exp(sol.x[0]))
    rho2_sol = float(np.exp(sol.x[1]))

    V2_sol = Vn * rho / rho2_sol
    p2_sol = p + rho * Vn**2 * (1.0 - rho / rho2_sol)

    c_vals, X_vals = partial_pressure(p2_sol, T2_sol, xi)

    return T2_sol, rho2_sol, p2_sol, V2_sol, X_vals, c_vals


def solve_ideal_oblique_shock(V1, theta, weak=True):
    M1 = V1 / a
    if M1 <= 1.0:
        raise ValueError("Oblique shock requires supersonic upstream flow (M1>1).")

    def theta_beta_mach(beta):
        sb = np.sin(beta)
        cb = np.cos(beta)
        denom = M1**2 * (gamma + np.cos(2.0 * beta)) + 2.0
        return 2.0 * cb / sb * (M1**2 * sb**2 - 1.0) / denom

    beta_min = np.arcsin(min(0.999999, 1.0 / max(M1, 1.000001))) + 1e-6
    beta_max = 0.5 * np.pi - 1e-6

    def g(beta):
        return np.tan(theta) - theta_beta_mach(beta)

    mid = 0.5 * (beta_min + beta_max)
    if weak:
        b1, b2 = beta_min, mid
    else:
        b1, b2 = mid, beta_max

    g1, g2 = g(b1), g(b2)
    if not (np.isfinite(g1) and np.isfinite(g2)):
        b1, b2 = beta_min, beta_max
        g1, g2 = g(b1), g(b2)

    if g1 * g2 > 0:
        b1, b2 = beta_min, beta_max
        g1, g2 = g(b1), g(b2)
        if g1 * g2 > 0:
            raise RuntimeError(
                "Failed to bracket a solution for beta from theta–beta–Mach relation. "
                "Theta may exceed the maximum deflection angle for this M1."
            )

    sol = root_scalar(g, bracket=[b1, b2], method='brentq', maxiter=200)
    if not sol.converged:
        raise RuntimeError("Failed to solve theta–beta–Mach relation for beta.")

    beta = float(sol.root)
    return beta


def solve_oblique_shock(Vin, theta, xi=0.25):
    def residuals(x):
        Vn1 = Vin * np.sin(x[0])
        T2_sol, rho2_sol, p2_sol, V2_sol, X_vals, c_vals = normal_shock_full(Vn1, xi=xi)
        Vn2 = rho / rho2_sol * Vn1
        res = np.tan(x[0] - theta) - Vn2 / Vn1 * np.tan(x[0])
        return np.array([res])

    b = solve_ideal_oblique_shock(Vin, theta)
    x0 = np.array([0.99 * b])

    sol = least_squares(
        residuals,
        x0=x0,
        xtol=1e-12,
        ftol=1e-12,
        gtol=1e-12,
        max_nfev=500
    )

    return float(sol.x[0])


# ==========================================================
#  EQUILIBRIUM SPEED OF SOUND FROM SPECIES
# ==========================================================

def mass_to_mole_fractions(c):
    c = np.array(c, dtype=float)

    s = np.sum(c)
    if s <= 0.0:
        raise ValueError("Species fractions sum to zero.")
    c = c / s

    n = c / MOLAR_MASS
    nsum = np.sum(n)
    if nsum <= 0.0:
        raise ValueError("Invalid mole-fraction conversion.")
    X = n / nsum

    return X


def infer_xi_from_mole_fractions(X):
    X = np.array(X, dtype=float)

    N_atoms = 2.0 * X[0] + X[2] + X[5]
    O_atoms = 2.0 * X[1] + X[3] + X[4]

    if N_atoms <= 0.0:
        raise ValueError("Nitrogen atoms are zero; cannot compute xi.")

    return O_atoms / N_atoms


def mixture_cp_frozen(T, c, dT=5.0):
    T1 = max(1.0, T - dT)
    T2 = T + dT

    h1 = mixture_enthalpy(T1, c)
    h2 = mixture_enthalpy(T2, c)

    cp = (h2 - h1) / (T2 - T1)

    if not np.isfinite(cp):
        raise ValueError("Non-finite cp.")

    return cp


def equilibrium_state_from_species(p_local, T_local, c_input):
    c_input = np.array(c_input, dtype=float)
    s = np.sum(c_input)
    if s <= 0.0:
        raise ValueError("Species fractions sum to zero.")
    c_input = c_input / s

    X = mass_to_mole_fractions(c_input)
    xi = infer_xi_from_mole_fractions(X)

    Rmix = float(np.dot(R_VALS, c_input))
    if Rmix <= 0.0 or not np.isfinite(Rmix):
        raise ValueError("Invalid Rmix.")

    rho_local = p_local / (Rmix * T_local)
    v = 1.0 / rho_local

    h_local = equilibrium_enthalpy(T_local, p_local / ATM_TO_PA, xi)

    e_local = h_local - p_local / rho_local

    if T_local <= 1200.0:
        cp = mixture_cp_frozen(T_local, c_input, dT=5.0)
        cv = cp - Rmix
    else:
        cp = cp_over_R(T_local, p_local / ATM_TO_PA, xi) * mixture_gas_constant_from_feed(xi)
        cv = cv_over_R(T_local, p_local / ATM_TO_PA, xi) * mixture_gas_constant_from_feed(xi)

    if cv <= 0.0 or not np.isfinite(cv):
        raise ValueError("cv <= 0 or non-finite.")

    gamma_f = cp / cv

    return {
        "c": c_input,
        "X": X,
        "Rmix": Rmix,
        "rho": rho_local,
        "v": v,
        "h": h_local,
        "e": e_local,
        "gamma_f": gamma_f
    }


def equilibrium_state_from_xi(p_local, T_local, xi):
    c_eq, X_eq = partial_pressure(p_local, T_local, xi)

    Rmix = float(np.dot(R_VALS, c_eq))
    if Rmix <= 0.0 or not np.isfinite(Rmix):
        raise ValueError("Invalid Rmix from equilibrium state.")

    rho_local = p_local / (Rmix * T_local)
    v = 1.0 / rho_local

    h_local = equilibrium_enthalpy(T_local, p_local / ATM_TO_PA, xi)

    e_local = h_local - p_local / rho_local

    if T_local <= 1200.0:
        cp = mixture_cp_frozen(T_local, c_eq, dT=5.0)
        cv = cp - Rmix
    else:
        cp = cp_over_R(T_local, p_local / ATM_TO_PA, xi) * mixture_gas_constant_from_feed(xi)
        cv = cv_over_R(T_local, p_local / ATM_TO_PA, xi) * mixture_gas_constant_from_feed(xi)

    if cv <= 0.0 or not np.isfinite(cv):
        raise ValueError("cv <= 0 or non-finite in equilibrium state.")

    gamma_f = cp / cv

    return {
        "c": c_eq,
        "X": X_eq,
        "Rmix": Rmix,
        "rho": rho_local,
        "v": v,
        "h": h_local,
        "e": e_local,
        "gamma_f": gamma_f
    }


def equilibrium_speed_of_sound_from_species(
        p_local, T_local,
        N2, O2, N, O, Op=0.0, Np=0.0, e=0.0,
        rel_dp=5e-3):

    c = np.array([N2, O2, N, O, Op, Np, e], dtype=float)

    if np.any(c < 0.0):
        raise ValueError("Negative species fraction.")

    s = np.sum(c)
    if s <= 0.0:
        raise ValueError("Species fractions sum to zero.")
    c = c / s

    s0 = equilibrium_state_from_species(p_local, T_local, c)
    xi = infer_xi_from_mole_fractions(s0["X"])
    phi_gamma = speed_sound_parameter(T_local, p_local / ATM_TO_PA, xi)
    ae2 = phi_gamma * p_local / s0["rho"]

    if not np.isfinite(ae2) or ae2 <= 0.0:
        raise ValueError("Negative or non-finite sound speed squared.")

    return np.sqrt(ae2)
# ==========================================================
#  TESTS AND PLOTS STARTS
# ==========================================================
def run_demo_calculations():
    """Run the original one-off checks and plot generation for this script."""
    try:
        ae = equilibrium_speed_of_sound_from_species(
            p_local=5000.0,
            T_local=4000.0,
            N2=0.55,
            O2=0.10,
            N=0.15,
            O=0.10,
            Op=0.05,
            Np=0.05,
            e=0.0
        )
        print("Equilibrium speed of sound =", ae, "m/s")
    except Exception as err:
        print("Single-point test failed:", err)

    n_points = 50
    deflection_angles = np.linspace(0.5, 30.0, n_points)
    shock_angles = np.zeros(n_points)

    for i in range(n_points):
        shock_angles[i] = solve_oblique_shock(3048.0, deflection_angles[i] / 180.0 * np.pi)

    plot_equilibrium_sound_speed_ratio()


def plot_equilibrium_sound_speed_ratio():
    import matplotlib.pyplot as plt

    xi = 0.25

    p_atm_list = [1e-4, 1e-3, 1e-2, 1e-1, 1.0, 10.0, 100.0]
    label_positions = {
        1e-4: 2500.0,
        1e-3: 3600.0,
        1e-2: 4700.0,
        1e-1: 6500.0,
        1.0: 7600.0,
        10.0: 9100.0,
        100.0: 12400.0,
    }
    atm_to_pa = ATM_TO_PA
    p_list = [p_atm * atm_to_pa for p_atm in p_atm_list]

    T_list = np.linspace(1.0, 15000.0, 500)

    skip_count = 0
    curves = []

    for p_atm, p_local in zip(p_atm_list, p_list):
        y_vals = []

        for T_local in T_list:
            try:
                c_eq, X_eq = partial_pressure(p_local, T_local, xi)

                ae_local = equilibrium_speed_of_sound_from_species(
                    p_local, T_local,
                    c_eq[0], c_eq[1], c_eq[2], c_eq[3], c_eq[4], c_eq[5], c_eq[6],
                    rel_dp=5e-3
                )

                s0 = equilibrium_state_from_species(p_local, T_local, c_eq)
                y = ae_local**2 / (p_local / s0["rho"])

                if np.isfinite(y) and (1.0 <= y <= 1.75):
                    y_vals.append(y)
                else:
                    y_vals.append(np.nan)

            except Exception:
                skip_count += 1
                y_vals.append(np.nan)

        y_vals = np.array(y_vals, dtype=float)
        curves.append((p_atm, y_vals))

    fig, ax = plt.subplots(figsize=(9, 5.7))

    for p_atm, y_vals in curves:
        ax.plot(T_list, y_vals, color="black", linewidth=1.15)

        label_T = label_positions[p_atm]
        label_y = np.interp(label_T, T_list, y_vals)
        ax.annotate(
            f"{p_atm:g}",
            xy=(label_T, label_y),
            xytext=(label_T - 650.0, label_y + 0.055),
            fontsize=11,
            ha="center",
            va="bottom",
            arrowprops=dict(arrowstyle="-", color="black", lw=0.8, shrinkA=2, shrinkB=2)
        )

    ax.set_xlabel("Temperature (K)", fontsize=15)
    ax.set_ylabel(r"$a_e^2 \rho / p$", fontsize=15)
    ax.set_xlim(0.0, 15000.0)
    ax.set_ylim(1.0, 1.7)
    ax.set_xticks([0, 5000, 10000, 15000])
    ax.set_yticks(np.arange(1.0, 1.71, 0.1))
    ax.minorticks_on()
    ax.grid(True, which="major", linestyle="-", linewidth=0.35, alpha=0.35)
    ax.grid(True, which="minor", linestyle=":", linewidth=0.3, alpha=0.25)
    ax.tick_params(direction="in", top=True, right=True, labelsize=12)
    fig.tight_layout()

    filename = "equilibrium_sound_speed_ratio_hansen_style.png"
    fig.savefig(filename, dpi=400, bbox_inches="tight")
    plt.close(fig)

    print(f"Figure saved as {filename}")

    color_fig, color_ax = plt.subplots(figsize=(10, 6))
    for p_atm, y_vals in curves:
        color_ax.plot(T_list, y_vals, linewidth=1.8, label=f"{p_atm:g} atm")

    color_ax.axhline(1.4, linestyle="--", linewidth=2.0)
    color_ax.text(3000.0, 1.405, "Frozen", fontsize=18, va="bottom")
    color_ax.set_xlabel("T, K", fontsize=20)
    color_ax.set_ylabel(r"$a_e^2/(p/\rho)$", fontsize=20)
    color_ax.set_xlim(0.0, 15000.0)
    color_ax.set_ylim(1.0, 1.7)
    color_ax.set_xticks([0, 3000, 6000, 9000, 12000, 15000])
    color_ax.tick_params(labelsize=16)
    color_ax.legend(fontsize=16)
    color_ax.grid(True, alpha=0.3)
    color_fig.tight_layout()

    color_filename = "equilibrium_sound_speed_ratio.png"
    color_fig.savefig(color_filename, dpi=300, bbox_inches="tight")
    plt.close(color_fig)

    print(f"Figure saved as {color_filename}")
    print("Skipped points:", skip_count)


# ==========================================================
#  TESTS AND PLOTS ENDS
# ==========================================================

# ==========================================================
#  PRANDTL-MEYER EXPANSION WITH EQUILIBRIUM REAL-GAS EFFECTS
#  O. Tumuklu, 2026-03-13
# ==========================================================

def _equilibrium_isentrope_temperature_step(
        p_prev, T_prev, p_new, xi,
        T_bounds=(40.0, 20000.0),
        max_tries=6,
        verbose=False):
    """
    March one step along an equilibrium isentrope using

        dh = v dp     with    ds = 0

    and solve for T_new at p_new.

    The residual is built from a trapezoidal approximation:

        h_new - h_prev - 0.5*(v_prev + v_new)*(p_new - p_prev) = 0.
    """
    st_prev = equilibrium_state_from_xi(p_prev, T_prev, xi)
    h_prev = st_prev["h"]
    v_prev = st_prev["v"]
    dp = p_new - p_prev

    T_low, T_high = T_bounds
    T_low = max(1.0, T_low)
    T_high = max(T_low + 1.0, T_high)

    logT_low = np.log(T_low)
    logT_high = np.log(T_high)

    def residual_logT(x):
        T_new = float(np.exp(x[0]))

        try:
            st_new = equilibrium_state_from_xi(p_new, T_new, xi)
            r = st_new["h"] - h_prev - 0.5 * (v_prev + st_new["v"]) * dp

            if not np.isfinite(r):
                return np.array([1.0e20], dtype=float)

            scale = max(abs(h_prev), abs(st_new["h"]), 1.0)
            return np.array([r / scale], dtype=float)

        except Exception:
            return np.array([1.0e20], dtype=float)

    guesses = [
        T_prev,
        max(T_low, min(T_high, 0.95 * T_prev)),
        max(T_low, min(T_high, 1.05 * T_prev)),
        max(T_low, min(T_high, 0.80 * T_prev)),
        max(T_low, min(T_high, 1.20 * T_prev)),
        0.5 * (T_low + T_high)
    ]

    best_sol = None
    best_cost = np.inf
    last_err = None

    for Tg in guesses[:max_tries]:
        x0 = np.array([np.log(Tg)], dtype=float)

        try:
            sol = least_squares(
                residual_logT,
                x0=x0,
                bounds=(np.array([logT_low]), np.array([logT_high])),
                xtol=1e-12,
                ftol=1e-12,
                gtol=1e-12,
                max_nfev=400
            )

            if sol.success and np.isfinite(sol.cost):
                if sol.cost < best_cost:
                    best_cost = sol.cost
                    best_sol = sol

                if sol.cost < 1.0e-24:
                    break

        except Exception as err:
            last_err = err

    if best_sol is None:
        raise ValueError(f"Failed to march isentrope temperature: {last_err}")

    T_best = float(np.exp(best_sol.x[0]))

    st_best = equilibrium_state_from_xi(p_new, T_best, xi)
    r_best = st_best["h"] - h_prev - 0.5 * (v_prev + st_best["v"]) * dp
    h_scale = max(abs(h_prev), abs(st_best["h"]), 1.0)

    if verbose:
        print(f"[isentrope step] "
              f"p_prev={p_prev:.6e}, p_new={p_new:.6e}, "
              f"T_prev={T_prev:.6f}, T_new={T_best:.6f}, "
              f"residual={r_best:.6e}, scaled={abs(r_best)/h_scale:.6e}, "
              f"cost={best_sol.cost:.6e}")

    if abs(r_best) / h_scale > 1.0e-8:
        raise ValueError(
            f"Least-squares isentrope step did not converge tightly enough: "
            f"T_new={T_best:.8f}, residual={r_best:.6e}, "
            f"scaled={abs(r_best)/h_scale:.6e}"
        )

    return T_best


def equilibrium_sound_speed_from_c(p_local, T_local, c, rel_dp=5e-3):
    """
    Convenience wrapper: equilibrium sound speed from a supplied mass-fraction vector.
    Species order is [N2, O2, N, O, O+, N+, e-].
    """
    c = np.array(c, dtype=float)
    s = np.sum(c)
    if s <= 0.0:
        raise ValueError("Species fractions sum to zero.")
    c = c / s

    return equilibrium_speed_of_sound_from_species(
        p_local, T_local,
        c[0], c[1], c[2], c[3], c[4], c[5], c[6],
        rel_dp=rel_dp
    )


def reacting_prandtl_meyer_given_theta(
        M1, T1, p1, xi, theta_deg,
        nsteps=400,
        pmin_ratio=1.0e-6,
        rel_dp_sound=5e-3,
        T_bounds=(150.0, 20000.0),
        verbose=False):
    """
    Equilibrium reacting-gas Prandtl-Meyer expansion with direct turning-angle input.

    Inputs
    ------
    M1        : upstream Mach number
    T1        : upstream temperature [K]
    p1        : upstream pressure [Pa]
    xi        : elemental ratio parameter
    theta_deg : requested flow deflection angle [deg]

    Returns
    -------
    dict with downstream state at the requested turning angle.
    """
    if M1 <= 1.0:
        raise ValueError("Prandtl-Meyer expansion requires M1 > 1.")
    if theta_deg < 0.0:
        raise ValueError("Deflection angle must be nonnegative for expansion.")
    if nsteps < 2:
        raise ValueError("nsteps must be at least 2.")

    # Upstream state is also built from equilibrium chemistry
    st1 = equilibrium_state_from_xi(p1, T1, xi)
    c1 = np.array(st1["c"], dtype=float)
    X1 = np.array(st1["X"], dtype=float)

    a1 = equilibrium_sound_speed_from_c(p1, T1, c1, rel_dp=rel_dp_sound)
    V1 = M1 * a1
    h0 = st1["h"] + 0.5 * V1 * V1

    if verbose:
        print("\n--- UPSTREAM STATE ---")
        print(f"M1 = {M1:.8f}")
        print(f"p1 = {p1:.8e} Pa")
        print(f"T1 = {T1:.8f} K")
        print(f"a1 = {a1:.8f} m/s")
        print(f"V1 = {V1:.8f} m/s")
        print(f"h1 = {st1['h']:.8e} J/kg")
        print(f"h0 = {h0:.8e} J/kg")
        print(f"c1_eq = {c1}")
        print(f"X1_eq = {X1}")
        print(f"xi = {xi:.8f}")

    p_end = max(p1 * pmin_ratio, 1.0e-2)
    p_grid = np.geomspace(p1, p_end, nsteps)

    states = [{
        "p": p1,
        "T": T1,
        "h": st1["h"],
        "rho": st1["rho"],
        "v": st1["v"],
        "a": a1,
        "V": V1,
        "M": M1,
        "c": c1.copy(),
        "X": X1.copy()
    }]

    theta_vals = [0.0]
    theta_target = np.deg2rad(theta_deg)
    theta_now = 0.0

    for i in range(1, len(p_grid)):
        prev = states[-1]
        p_new = float(p_grid[i])

        try:
            T_new = _equilibrium_isentrope_temperature_step(
                prev["p"], prev["T"], p_new, xi,
                T_bounds=T_bounds,
                verbose=False
            )
        except Exception:
            T_new = _equilibrium_isentrope_temperature_step(
                prev["p"], prev["T"], p_new, xi,
                T_bounds=(50.0, 40000.0),
                verbose=False
            )

        st_eq = equilibrium_state_from_xi(p_new, T_new, xi)
        c_new = np.array(st_eq["c"], dtype=float)
        X_new = np.array(st_eq["X"], dtype=float)
        a_new = equilibrium_sound_speed_from_c(p_new, T_new, c_new, rel_dp=rel_dp_sound)

        V2_sq = 2.0 * (h0 - st_eq["h"])
        if V2_sq <= 0.0:
            raise ValueError(
                f"Non-positive velocity squared encountered at step {i}: "
                f"2*(h0-h) = {V2_sq}"
            )

        V_new = np.sqrt(V2_sq)
        M_new = V_new / a_new

        if prev["M"] > 1.0 and M_new > 1.0 and prev["V"] > 0.0 and V_new > 0.0:
            f_prev = np.sqrt(prev["M"]**2 - 1.0) / prev["V"]
            f_new  = np.sqrt(M_new**2 - 1.0) / V_new
            dtheta = 0.5 * (f_prev + f_new) * (V_new - prev["V"])
        else:
            dtheta = 0.0

        theta_next = theta_now + dtheta

        current = {
            "p": p_new,
            "T": T_new,
            "h": st_eq["h"],
            "rho": st_eq["rho"],
            "v": st_eq["v"],
            "a": a_new,
            "V": V_new,
            "M": M_new,
            "c": c_new,
            "X": X_new
        }

        states.append(current)
        theta_vals.append(theta_next)

        if verbose:
            print(
                f"i={i:4d}, theta={np.degrees(theta_next):10.6f} deg, "
                f"p={p_new:12.6e} Pa, T={T_new:10.4f} K, M={M_new:10.6f}"
            )

        if theta_next >= theta_target:
            th0 = theta_now
            th1 = theta_next
            w = 1.0 if abs(th1 - th0) < 1e-14 else (theta_target - th0) / (th1 - th0)

            def interp_scalar(key):
                return prev[key] + w * (current[key] - prev[key])

            c2 = prev["c"] + w * (current["c"] - prev["c"])
            X2 = prev["X"] + w * (current["X"] - prev["X"])

            c2 = np.maximum(c2, 0.0)
            X2 = np.maximum(X2, 0.0)

            if np.sum(c2) > 0.0:
                c2 /= np.sum(c2)
            if np.sum(X2) > 0.0:
                X2 /= np.sum(X2)

            result = {
                "M2": interp_scalar("M"),
                "p2": interp_scalar("p"),
                "T2": interp_scalar("T"),
                "a2": interp_scalar("a"),
                "V2": interp_scalar("V"),
                "rho2": interp_scalar("rho"),
                "c2": c2,
                "X2": X2,
                "theta_deg": theta_deg,
                "xi": xi,
                "states": states,
                "theta_vals_deg": np.degrees(np.array(theta_vals))
            }
            return result

        theta_now = theta_next

    raise ValueError(
        "Target deflection angle was not reached. Increase nsteps or decrease pmin_ratio."
    )





def extract_pm_arrays(result):
    states = result["states"]
    theta_deg = np.array(result["theta_vals_deg"], dtype=float)

    p   = np.array([s["p"]   for s in states], dtype=float)
    T   = np.array([s["T"]   for s in states], dtype=float)
    rho = np.array([s["rho"] for s in states], dtype=float)
    v   = np.array([s["v"]   for s in states], dtype=float)
    h   = np.array([s["h"]   for s in states], dtype=float)
    a   = np.array([s["a"]   for s in states], dtype=float)
    V   = np.array([s["V"]   for s in states], dtype=float)
    M   = np.array([s["M"]   for s in states], dtype=float)

    c = np.array([s["c"] for s in states], dtype=float)
    X = np.array([s["X"] for s in states], dtype=float)

    return {
        "theta_deg": theta_deg,
        "p": p,
        "T": T,
        "rho": rho,
        "v": v,
        "h": h,
        "a": a,
        "V": V,
        "M": M,
        "c": c,
        "X": X
    }


def prandtl_meyer_nu(M, gamma=1.4):
    """
    Classical Prandtl-Meyer function for a calorically perfect gas.
    Returns nu in radians.
    """
    M = np.asarray(M, dtype=float)
    if np.any(M < 1.0):
        raise ValueError("Prandtl-Meyer function requires M >= 1.")

    term1 = np.sqrt((gamma + 1.0) / (gamma - 1.0))
    term2 = np.sqrt((gamma - 1.0) / (gamma + 1.0) * (M**2 - 1.0))
    term3 = np.sqrt(M**2 - 1.0)

    return term1 * np.arctan(term2) - np.arctan(term3)


def prandtl_meyer_inverse(nu_target, gamma=1.4, M_low=1.000001, M_high=100.0):
    """
    Solve nu(M)=nu_target for M using Brent's method.
    """
    def f(M):
        return prandtl_meyer_nu(M, gamma=gamma) - nu_target

    while f(M_high) < 0.0:
        M_high *= 2.0
        if M_high > 1.0e6:
            raise ValueError("Could not bracket M in prandtl_meyer_inverse.")

    sol = root_scalar(f, bracket=[M_low, M_high], method='brentq', maxiter=300)
    if not sol.converged:
        raise ValueError("Failed to invert Prandtl-Meyer function.")
    return float(sol.root)


def perfect_gas_pm_comparison(M1, theta_deg_array, gamma=1.4):
    """
    For a calorically perfect gas, compute M2, T2/T1, and p2/p1
    as functions of theta.
    """
    theta_deg_array = np.asarray(theta_deg_array, dtype=float)
    theta_rad_array = np.deg2rad(theta_deg_array)

    nu1 = prandtl_meyer_nu(M1, gamma=gamma)

    M2 = np.zeros_like(theta_rad_array)
    T2_T1 = np.zeros_like(theta_rad_array)
    p2_p1 = np.zeros_like(theta_rad_array)

    T0_over_T1 = 1.0 + 0.5 * (gamma - 1.0) * M1**2

    for i, th in enumerate(theta_rad_array):
        nu2 = nu1 + th
        M2[i] = prandtl_meyer_inverse(nu2, gamma=gamma)

        T0_over_T2 = 1.0 + 0.5 * (gamma - 1.0) * M2[i]**2
        T2_T1[i] = T0_over_T1 / T0_over_T2
        p2_p1[i] = (T2_T1[i])**(gamma / (gamma - 1.0))

    return {
        "theta_deg": theta_deg_array,
        "M2": M2,
        "T2_T1": T2_T1,
        "p2_p1": p2_p1
    }

def plot_pm_combined_one_graph(result, M1, gamma_pg=1.4, use_log_y=True):
    import matplotlib.pyplot as plt

    data = extract_pm_arrays(result)

    theta_deg = data["theta_deg"]
    p1 = data["p"][0]
    T1 = data["T"][0]

    # Reacting-air curves
    M2_reacting = data["M"]
    T1T2_reacting = T1 / data["T"]
    p1p2_reacting = p1 / data["p"]

    # Perfect-gas curves
    pg = perfect_gas_pm_comparison(M1, theta_deg, gamma=gamma_pg)

    T1T2_pg = 1.0 / pg["T2_T1"]
    p1p2_pg = 1.0 / pg["p2_p1"]

    plt.figure(figsize=(9,6))

    # ----------------------------
    # M2
    # ----------------------------
    plt.plot(theta_deg, M2_reacting,
             color="blue", linewidth=3,
             label=r"Reacting air $M_2$")

    plt.plot(pg["theta_deg"], pg["M2"],
             color="blue", linestyle="--", linewidth=3,
             marker='o', markerfacecolor='none', markevery=8,
             label=r"Perfect gas $M_2$")

    # ----------------------------
    # T1/T2
    # ----------------------------
    plt.plot(theta_deg, T1T2_reacting,
             color="red", linewidth=3,
             label=r"Reacting air $T_1/T_2$")

    plt.plot(pg["theta_deg"], T1T2_pg,
             color="red", linestyle="--", linewidth=3,
             marker='o', markerfacecolor='none', markevery=8,
             label=r"Perfect gas $T_1/T_2$")

    # ----------------------------
    # p1/p2
    # ----------------------------
    plt.plot(theta_deg, p1p2_reacting,
             color="green", linewidth=3,
             label=r"Reacting air $p_1/p_2$")

    plt.plot(pg["theta_deg"], p1p2_pg,
             color="green", linestyle="--", linewidth=3,
             marker='o', markerfacecolor='none', markevery=8,
             label=r"Perfect gas $p_1/p_2$")

    # Axis labels
    plt.xlabel(r"Turning Angle $\theta$ (deg)", fontsize=20)
    plt.ylabel(r"$M_2$ , $T_1/T_2$ , $p_1/p_2$", fontsize=20)

    if use_log_y:
        plt.yscale("log")

    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)

    from matplotlib.ticker import MultipleLocator
    plt.gca().xaxis.set_major_locator(MultipleLocator(5))

    plt.grid(True, which="both", linestyle="--", alpha=0.35)

    plt.legend(fontsize=14, ncol=2)

    plt.tight_layout()

    plt.savefig("pm_combined_one_graph.png",
                dpi=400,
                bbox_inches="tight")

    plt.show()

    print("Saved: pm_combined_one_graph.png")



# ==========================================================
#  ADDITIONAL PM FIGURES 
#  
# ==========================================================

def plot_pm_multifigure_summary(result, M1, gamma_pg=1.4,
                                use_mole_fractions=True,
                                species_filename="PM_species_summary.png",
                                thermo_filename="PM_thermo_summary.png"):
    import matplotlib.pyplot as plt

    """
    Appended plotting utility:
      Figure 1: M2, T2/T1, p2/p1 vs theta with reacting air vs perfect gas
      Figure 2: species fractions + V/a_e + rho/rho1 vs theta

    This does not touch the existing one-graph plot.
    """
    from matplotlib.ticker import MultipleLocator

    data = extract_pm_arrays(result)

    theta_deg = data["theta_deg"]
    p1 = data["p"][0]
    T1 = data["T"][0]
    rho1 = data["rho"][0]

    # Reacting-air curves
    M2_reacting = data["M"]
    T2T1_reacting = data["T"] / T1
    p2p1_reacting = data["p"] / p1

    # Perfect-gas comparison
    pg = perfect_gas_pm_comparison(M1, theta_deg, gamma=gamma_pg)

    # Species selection
    Y = data["X"] if use_mole_fractions else data["c"]
    ylab = "Mole fraction" if use_mole_fractions else "Mass fraction"

    # ------------------------------------------------------
    # FIGURE 1: Thermodynamic comparison
    # ------------------------------------------------------
    fig1, axes1 = plt.subplots(3, 1, figsize=(8, 12), sharex=True)

    # M2
    axes1[0].plot(theta_deg, M2_reacting,
                  color="blue", linewidth=3, label=r"Reacting air")
    axes1[0].plot(pg["theta_deg"], pg["M2"],
                  color="blue", linestyle="--", linewidth=3,
                  marker='o', markerfacecolor='none', markevery=8,
                  label=r"Perfect gas")
    axes1[0].set_ylabel(r"$M_2$", fontsize=18)
    axes1[0].grid(True, which="both", linestyle="--", alpha=0.35)
    axes1[0].legend(fontsize=13)

    # T2/T1
    axes1[1].plot(theta_deg, T2T1_reacting,
                  color="red", linewidth=3, label=r"Reacting air")
    axes1[1].plot(pg["theta_deg"], pg["T2_T1"],
                  color="red", linestyle="--", linewidth=3,
                  marker='o', markerfacecolor='none', markevery=8,
                  label=r"Perfect gas")
    axes1[1].set_ylabel(r"$T_2/T_1$", fontsize=18)
    axes1[1].grid(True, which="both", linestyle="--", alpha=0.35)
    axes1[1].legend(fontsize=13)

    # p2/p1
    axes1[2].plot(theta_deg, p2p1_reacting,
                  color="green", linewidth=3, label=r"Reacting air")
    axes1[2].plot(pg["theta_deg"], pg["p2_p1"],
                  color="green", linestyle="--", linewidth=3,
                  marker='o', markerfacecolor='none', markevery=8,
                  label=r"Perfect gas")
    axes1[2].set_ylabel(r"$p_2/p_1$", fontsize=18)
    axes1[2].set_xlabel(r"Turning Angle $\theta$ (deg)", fontsize=18)
    axes1[2].grid(True, which="both", linestyle="--", alpha=0.35)
    axes1[2].legend(fontsize=13)

    for ax in axes1:
        ax.tick_params(axis='both', labelsize=14)
        ax.xaxis.set_major_locator(MultipleLocator(5))

    fig1.tight_layout()
    fig1.savefig(thermo_filename, dpi=400, bbox_inches="tight")
    plt.close(fig1)

    # ------------------------------------------------------
    # FIGURE 2: Species + velocities + density
    # ------------------------------------------------------
    fig2, axes2 = plt.subplots(3, 1, figsize=(8, 13), sharex=True)

    # Species fractions
    for i in range(len(SPECIES_NAMES)):
        axes2[0].plot(theta_deg, Y[:, i],
                      color=SPECIES_COLORS[i],
                      linewidth=3,
                      label=SPECIES_LABELS[i])
    axes2[0].set_ylabel(ylab, fontsize=18)
    axes2[0].grid(True, which="both", linestyle="--", alpha=0.35)
    axes2[0].legend(fontsize=12, ncol=2)

    # Velocity and equilibrium sound speed
    axes2[1].plot(theta_deg, data["V"],
                  color="black", linewidth=3, label=r"$V$")
    axes2[1].plot(theta_deg, data["a"],
                  color="cyan", linewidth=3, label=r"$a_e$")
    axes2[1].set_ylabel(r"Velocity (m/s)", fontsize=18)
    axes2[1].grid(True, which="both", linestyle="--", alpha=0.35)
    axes2[1].legend(fontsize=13)

    # Density ratio
    axes2[2].plot(theta_deg, data["rho"] / rho1,
              color="brown",
              linestyle="--",
              marker="o",
              markevery=10,
              markersize=6,
              linewidth=3,
              label=r"$\rho/\rho_1$")

    axes2[2].set_ylabel(r"$\rho/\rho_1$", fontsize=18)
    axes2[2].set_xlabel(r"Turning Angle $\theta$ (deg)", fontsize=18)

    axes2[2].grid(True, which="both", linestyle="--", alpha=0.35)

    axes2[2].legend(fontsize=13)

    for ax in axes2:
        ax.tick_params(axis='both', labelsize=14)
        ax.xaxis.set_major_locator(MultipleLocator(5))

    fig2.tight_layout()
    fig2.savefig(species_filename, dpi=400, bbox_inches="tight")
    plt.close(fig2)

    print(f"Saved: {thermo_filename}")
    print(f"Saved: {species_filename}")


def plot_pm_species_only(result, use_mole_fractions=True,
                         filename="PM_species_only.png"):
    """
    Extra appended helper:
    species-only figure versus theta.
    """
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MultipleLocator

    data = extract_pm_arrays(result)
    theta_deg = data["theta_deg"]

    Y = data["X"] if use_mole_fractions else data["c"]
    ylab = "Mole fraction" if use_mole_fractions else "Mass fraction"

    plt.figure(figsize=(8, 5.5))

    for i in range(len(SPECIES_NAMES)):
        plt.plot(theta_deg, Y[:, i],
                 color=SPECIES_COLORS[i],
                 linewidth=3,
                 label=SPECIES_LABELS[i])

    plt.xlabel(r"Turning Angle $\theta$ (deg)", fontsize=18)
    plt.ylabel(ylab, fontsize=18)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.gca().xaxis.set_major_locator(MultipleLocator(5))
    plt.grid(True, which="both", linestyle="--", alpha=0.35)
    plt.legend(fontsize=13, ncol=2)
    plt.tight_layout()
    plt.savefig(filename, dpi=400, bbox_inches="tight")
    plt.close()

    print(f"Saved: {filename}")


# ==========================================================
#  APPENDED CALLS FOR FIGURES AND TESTS
# ==========================================================
if __name__ == "__main__":
    run_demo_calculations()

    # Example usage
    xi = 1

    res = reacting_prandtl_meyer_given_theta(
        M1=5.0,
        T1=3000.0,
        p1=1000.0,
        xi=xi,
        theta_deg=50.0,
        nsteps=500,
        verbose=True
    )

    plot_pm_multifigure_summary(
        res,
        M1=5.0,
        gamma_pg=1.4,
        use_mole_fractions=True,
        species_filename="PM_species_summary.png",
        thermo_filename="PM_thermo_summary.png"
    )

    plot_pm_species_only(
        res,
        use_mole_fractions=True,
        filename="PM_species_only.png"
    )

    print("M2 =", res["M2"])
    print("p2 =", res["p2"])
    print("T2 =", res["T2"])
    print("c2 =", res["c2"])
    print("X2 =", res["X2"])

    plot_pm_combined_one_graph(res, M1=5.0, gamma_pg=1.4, use_log_y=True)
