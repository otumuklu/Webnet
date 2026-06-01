import numpy as np
import matplotlib.pyplot as plt

# ============================================================
# Hansen / Axioms high-temperature equilibrium air model
# Pressure is in atm, temperature is in K.
# Valid mainly for high-temperature air up to ~15000 K.
# Developed by Ozgur Tumuklu at April, 28th 2026
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


def hansen_air(T, p_atm=1.0):
    """
    Returns compressibility factor and species mole fractions.
    """
    eps1, eps2, eps3 = epsilons(T, p_atm)

    Z = 1.0 + eps1 + eps2 + 2.0*eps3

    xO2 = (0.2 - eps1) / Z
    xN2 = (0.8 - eps2) / Z
    xO  = (2.0*eps1 - 0.4*eps3) / Z
    xN  = (2.0*eps2 - 1.6*eps3) / Z

    # Axioms paper gives x(O+ + N+) = x(e-) = 2 eps3 / Z.
    # Here split ions by 20/80 oxygen/nitrogen approximation.
    xe  = 2.0*eps3 / Z
    xOp = 0.2*xe
    xNp = 0.8*xe

    return {
        "T": T,
        "p_atm": p_atm,
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


# ============================================================
# Figure 5: speed-of-sound parameter a^2 rho / p = gamma Phi
# Numerical thermodynamic derivative version
# ============================================================

E0_over_R = {
    "N2": 0.0,
    "O2": 0.0,
    "O": 59000.0 / 2.0,
    "N": 113200.0 / 2.0,
    "O+": 59000.0 / 2.0 + 158000.0,
    "N+": 113200.0 / 2.0 + 168800.0,
    "e-": 0.0,
}

species_list = ["O2", "N2", "O", "N", "O+", "N+", "e-"]
x_keys = {
    "O2": "xO2",
    "N2": "xN2",
    "O": "xO",
    "N": "xN",
    "O+": "xO+",
    "N+": "xN+",
    "e-": "xe-",
}


def dlnQp_dT(species, T):
    h = 1e-3 * T
    return (lnQp(species, T + h) - lnQp(species, T - h)) / (2.0*h)


def species_H_RT(species, T):
    return T*dlnQp_dT(species, T) + E0_over_R[species]/T


def species_E_RT(species, T):
    # For ideal-gas species: H = E + RT
    return species_H_RT(species, T) - 1.0


def mixture_ZH_RT(T, p_atm):
    sol = hansen_air(T, p_atm)
    total = 0.0
    for sp in species_list:
        total += sol[x_keys[sp]] * species_H_RT(sp, T)
    return sol["Z"] * total


def mixture_ZE_RT(T, p_atm):
    sol = hansen_air(T, p_atm)
    total = 0.0
    for sp in species_list:
        total += sol[x_keys[sp]] * species_E_RT(sp, T)
    return sol["Z"] * total


def pressure_for_constant_density(T_new, T_ref, p_ref):
    """
    Constant density condition:
    p / (Z T) = constant.
    Solve p_new / (Z_new T_new) = p_ref / (Z_ref T_ref).
    """
    Z_ref = hansen_air(T_ref, p_ref)["Z"]
    target = p_ref / (Z_ref*T_ref)

    p = p_ref * T_new/T_ref

    for _ in range(30):
        Z = hansen_air(T_new, p)["Z"]
        f = p/(Z*T_new) - target

        dp = 1e-5 * max(p, 1e-8)
        Zp = hansen_air(T_new, p + dp)["Z"]
        fp = (p + dp)/(Zp*T_new) - target

        dfdp = (fp - f)/dp
        p -= f/dfdp
        p = max(p, 1e-12)

    return p


def cp_over_R(T, p_atm):
    h = 1e-3*T

    Yp = mixture_ZH_RT(T + h, p_atm)
    Ym = mixture_ZH_RT(T - h, p_atm)

    # Cp/R = d[(ZH/RT) T]/dT at constant pressure
    return ((T + h)*Yp - (T - h)*Ym) / (2.0*h)


def cv_over_R(T, p_atm):
    h = 1e-3*T

    pp = pressure_for_constant_density(T + h, T, p_atm)
    pm = pressure_for_constant_density(T - h, T, p_atm)

    Xp = mixture_ZE_RT(T + h, pp)
    Xm = mixture_ZE_RT(T - h, pm)

    # Cv/R = d[(ZE/RT) T]/dT at constant density
    return ((T + h)*Xp - (T - h)*Xm) / (2.0*h)


def phi_factor(T, p_atm):
    """
    Phi = (rho/p) (dp/drho)_T.
    Since rho proportional to p/(Z T), evaluate numerically.
    """
    Z0 = hansen_air(T, p_atm)["Z"]
    rho0_scaled = p_atm/(Z0*T)

    dp = 1e-4 * max(p_atm, 1e-8)

    p1 = p_atm + dp
    Z1 = hansen_air(T, p1)["Z"]
    rho1_scaled = p1/(Z1*T)

    p2 = p_atm - dp
    p2 = max(p2, 1e-12)
    Z2 = hansen_air(T, p2)["Z"]
    rho2_scaled = p2/(Z2*T)

    dp_drho = (p1 - p2)/(rho1_scaled - rho2_scaled)

    return rho0_scaled/p_atm * dp_drho


def speed_sound_parameter(T, p_atm):
    Cp = cp_over_R(T, p_atm)
    Cv = cv_over_R(T, p_atm)
    gamma = Cp/Cv
    Phi = phi_factor(T, p_atm)
    return gamma*Phi


# ============================================================
# Example: reproduce mole fractions and Z at 1 atm
# ============================================================

if __name__ == "__main__":

    T = np.linspace(2000.0, 15000.0, 500)
    p_atm = 1.0

    sol = hansen_air(T, p_atm)

    # Plot mole fractions
    plt.figure()
    for key in ["xO2", "xN2", "xO", "xN", "xO+", "xN+", "xe-"]:
        plt.plot(T, sol[key], label=key)

    plt.xlabel("Temperature [K]")
    plt.ylabel("Mole fraction")
    plt.title("Hansen equilibrium air mole fractions at 1 atm")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

    # Plot compressibility factor
    plt.figure()
    plt.plot(T, sol["Z"])
    plt.xlabel("Temperature [K]")
    plt.ylabel("Z")
    plt.title("Hansen compressibility factor at 1 atm")
    plt.grid(True)
    plt.tight_layout()
    plt.show()

    # Print one condition
    T0 = 9000.0
    out = hansen_air(T0, p_atm)

    print(f"T = {T0:.1f} K, p = {p_atm:.3g} atm")
    print(f"Z = {out['Z']:.6f}")
    for key in ["xO2", "xN2", "xO", "xN", "xO+", "xN+", "xe-"]:
        print(f"{key:4s} = {out[key]:.6e}")

# ============================================================
# Figure 5: speed-of-sound parameter a^2 rho / p = gamma Phi
# Numerical thermodynamic derivative version
# ============================================================

E0_over_R = {
    "N2": 0.0,
    "O2": 0.0,
    "O": 59000.0 / 2.0,
    "N": 113200.0 / 2.0,
    "O+": 59000.0 / 2.0 + 158000.0,
    "N+": 113200.0 / 2.0 + 168800.0,
    "e-": 0.0,
}

species_list = ["O2", "N2", "O", "N", "O+", "N+", "e-"]
x_keys = {
    "O2": "xO2",
    "N2": "xN2",
    "O": "xO",
    "N": "xN",
    "O+": "xO+",
    "N+": "xN+",
    "e-": "xe-",
}


def dlnQp_dT(species, T):
    h = 1e-3 * T
    return (lnQp(species, T + h) - lnQp(species, T - h)) / (2.0*h)


def species_H_RT(species, T):
    return T*dlnQp_dT(species, T) + E0_over_R[species]/T


def species_E_RT(species, T):
    # For ideal-gas species: H = E + RT
    return species_H_RT(species, T) - 1.0


def mixture_ZH_RT(T, p_atm):
    sol = hansen_air(T, p_atm)
    total = 0.0
    for sp in species_list:
        total += sol[x_keys[sp]] * species_H_RT(sp, T)
    return sol["Z"] * total


def mixture_ZE_RT(T, p_atm):
    sol = hansen_air(T, p_atm)
    total = 0.0
    for sp in species_list:
        total += sol[x_keys[sp]] * species_E_RT(sp, T)
    return sol["Z"] * total


def pressure_for_constant_density(T_new, T_ref, p_ref):
    """
    Constant density condition:
    p / (Z T) = constant.
    Solve p_new / (Z_new T_new) = p_ref / (Z_ref T_ref).
    """
    Z_ref = hansen_air(T_ref, p_ref)["Z"]
    target = p_ref / (Z_ref*T_ref)

    p = p_ref * T_new/T_ref

    for _ in range(30):
        Z = hansen_air(T_new, p)["Z"]
        f = p/(Z*T_new) - target

        dp = 1e-5 * max(p, 1e-8)
        Zp = hansen_air(T_new, p + dp)["Z"]
        fp = (p + dp)/(Zp*T_new) - target

        dfdp = (fp - f)/dp
        p -= f/dfdp
        p = max(p, 1e-12)

    return p


def cp_over_R(T, p_atm):
    h = 1e-3*T

    Yp = mixture_ZH_RT(T + h, p_atm)
    Ym = mixture_ZH_RT(T - h, p_atm)

    # Cp/R = d[(ZH/RT) T]/dT at constant pressure
    return ((T + h)*Yp - (T - h)*Ym) / (2.0*h)


def cv_over_R(T, p_atm):
    h = 1e-3*T

    pp = pressure_for_constant_density(T + h, T, p_atm)
    pm = pressure_for_constant_density(T - h, T, p_atm)

    Xp = mixture_ZE_RT(T + h, pp)
    Xm = mixture_ZE_RT(T - h, pm)

    # Cv/R = d[(ZE/RT) T]/dT at constant density
    return ((T + h)*Xp - (T - h)*Xm) / (2.0*h)


def phi_factor(T, p_atm):
    """
    Phi = (rho/p) (dp/drho)_T.
    Since rho proportional to p/(Z T), evaluate numerically.
    """
    Z0 = hansen_air(T, p_atm)["Z"]
    rho0_scaled = p_atm/(Z0*T)

    dp = 1e-4 * max(p_atm, 1e-8)

    p1 = p_atm + dp
    Z1 = hansen_air(T, p1)["Z"]
    rho1_scaled = p1/(Z1*T)

    p2 = p_atm - dp
    p2 = max(p2, 1e-12)
    Z2 = hansen_air(T, p2)["Z"]
    rho2_scaled = p2/(Z2*T)

    dp_drho = (p1 - p2)/(rho1_scaled - rho2_scaled)

    return rho0_scaled/p_atm * dp_drho


def speed_sound_parameter(T, p_atm):
    Cp = cp_over_R(T, p_atm)
    Cv = cv_over_R(T, p_atm)
    gamma = Cp/Cv
    Phi = phi_factor(T, p_atm)
    return gamma*Phi


# ============================================================
# Draw Figure 5
# ============================================================

T_range = np.linspace(1000.0, 15000.0, 250)
pressures = [1e-4, 1e-3, 1e-2, 1e-1, 1, 10, 100]

plt.figure()

for p in pressures:
    y = np.array([speed_sound_parameter(T, p) for T in T_range])
    plt.plot(T_range, y, label=f"{p:g} atm")

plt.xlabel("Temperature [K]")
plt.ylabel(r"$a^2 \rho / p = \gamma \Phi$")
plt.title("Figure 5: Hansen speed-of-sound parameter")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()
