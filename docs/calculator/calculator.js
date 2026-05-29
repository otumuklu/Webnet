const DEG = Math.PI / 180;
const RAD = 180 / Math.PI;
const UNIVERSAL_GAS_CONSTANT = 8314.462618;
const ATM_TO_PA = 101325.0;
const HANSEN_R_UNIVERSAL = 8.31446261815324;

const GAS_PROPERTIES = {
  air: { label: "Air", gamma: 1.4, molecularWeight: 28.965 },
  co2: { label: "CO₂", gamma: 1.289, molecularWeight: 44.01 },
  argon: { label: "Argon", gamma: 1.667, molecularWeight: 39.948 },
  nitrogen: { label: "N₂", gamma: 1.4, molecularWeight: 28.0134 },
  oxygen: { label: "O₂", gamma: 1.4, molecularWeight: 31.9988 },
  helium: { label: "Helium", gamma: 1.667, molecularWeight: 4.0026 },
};

const HANSEN_SPECIES = ["N2", "O2", "N", "O", "O+", "N+", "e-"];
const HANSEN_SPECIES_LABELS = {
  N2: "N<sub>2</sub>",
  O2: "O<sub>2</sub>",
  N: "N",
  O: "O",
  "O+": "O<sup>+</sup>",
  "N+": "N<sup>+</sup>",
  "e-": "e<sup>-</sup>",
};
const HANSEN_MOLAR_MASS = [28.02, 32.00, 14.01, 16.00, 16.00, 14.01, 5.485799e-4];
const HANSEN_SPECIFIC_GAS_CONSTANTS = HANSEN_MOLAR_MASS.map((molecularWeight) => (
  UNIVERSAL_GAS_CONSTANT / molecularWeight
));
const HANSEN_FORMATION_OVER_R = {
  N2: 0,
  O2: 0,
  O: 59000 / 2,
  N: 113200 / 2,
  "O+": 59000 / 2 + 158000,
  "N+": 113200 / 2 + 168800,
  "e-": 0,
};
const HANSEN_DEFAULT_XI = 0.25;
const HANSEN_UNDISSOCIATED_AIR_MW = 0.8 * HANSEN_MOLAR_MASS[0] + 0.2 * HANSEN_MOLAR_MASS[1];
const HANSEN_ELECTRON_MW = HANSEN_MOLAR_MASS[6];
const HANSEN_COLLISION_TABLE = [
  [500, 38.4, 0.946, 0.894, null, null, 0.877, 0.761],
  [1000, 34.9, 0.920, 0.838, null, null, 0.843, 0.703],
  [1500, 33.7, 0.889, 0.785, null, null, 0.817, 0.652],
  [2000, 33.2, 0.886, 0.742, null, null, 0.794, 0.611],
  [2500, 32.8, 0.846, 0.705, null, null, 0.775, 0.578],
  [3000, 32.6, 0.830, 0.675, null, null, 0.759, 0.551],
  [3500, 32.4, 0.815, 0.650, null, null, 0.745, 0.527],
  [4000, 32.3, 0.803, 0.628, null, null, 0.733, 0.507],
  [4500, 32.2, 0.792, 0.608, null, null, 0.722, 0.489],
  [5000, 32.1, 0.782, 0.591, null, null, 0.712, 0.473],
  [5500, 32.0, 0.773, 0.575, 0.397, 89.9, 0.703, 0.458],
  [6000, 32.0, 0.764, 0.561, 0.380, 75.6, 0.695, 0.445],
  [6500, 31.9, 0.757, 0.548, 0.366, 64.5, 0.688, 0.433],
  [7000, 31.9, 0.750, 0.536, 0.353, 55.7, 0.681, 0.422],
  [7500, 31.9, 0.743, 0.524, 0.342, 48.6, 0.674, 0.412],
  [8000, 31.8, 0.737, 0.514, 0.331, 42.8, 0.668, 0.402],
  [8500, 31.8, 0.731, 0.504, 0.321, 37.9, 0.662, 0.393],
  [9000, 31.8, 0.725, 0.495, 0.313, 33.8, 0.657, 0.385],
  [9500, 31.8, 0.720, 0.485, 0.304, 30.4, 0.652, 0.377],
  [10000, 31.8, 0.715, 0.478, 0.297, 27.4, 0.647, 0.370],
  [10500, 31.7, 0.710, 0.470, 0.290, 24.9, 0.642, 0.363],
  [11000, 31.7, 0.706, 0.463, 0.283, 22.7, 0.637, 0.356],
  [11500, 31.7, 0.701, 0.456, 0.281, 20.8, 0.633, 0.350],
  [12000, 31.7, 0.697, 0.448, 0.270, 19.0, 0.629, 0.342],
  [12500, 31.7, 0.693, 0.443, 0.266, 17.6, 0.625, 0.338],
  [13000, 31.7, 0.689, 0.437, 0.261, 16.27, 0.621, 0.332],
  [13500, 31.7, 0.684, 0.431, 0.256, 15.10, 0.618, 0.327],
  [14000, 31.7, 0.681, 0.426, 0.252, 14.00, 0.616, 0.322],
  [14500, 31.6, 0.420, 0.247, 0.247, 13.09, 0.613, 0.316],
  [15000, 31.6, 0.415, 0.243, 0.243, 12.24, 0.610, 0.312],
];

function gasConstantFromMolecularWeight(molecularWeight) {
  if (molecularWeight <= 0 || !Number.isFinite(molecularWeight)) {
    throw new Error("Molecular weight must be greater than 0.");
  }
  return UNIVERSAL_GAS_CONSTANT / molecularWeight;
}

function tt0(gamma, mach) {
  return Math.pow(1 + ((gamma - 1) / 2) * mach * mach, -1);
}

function pp0(gamma, mach) {
  return Math.pow(1 + ((gamma - 1) / 2) * mach * mach, -gamma / (gamma - 1));
}

function rr0(gamma, mach) {
  return Math.pow(1 + ((gamma - 1) / 2) * mach * mach, -1 / (gamma - 1));
}

function tts(gamma, mach) {
  return tt0(gamma, mach) * ((gamma + 1) / 2);
}

function pps(gamma, mach) {
  return pp0(gamma, mach) * Math.pow((gamma + 1) / 2, gamma / (gamma - 1));
}

function rrs(gamma, mach) {
  return rr0(gamma, mach) * Math.pow((gamma + 1) / 2, 1 / (gamma - 1));
}

function areaRatio(gamma, mach) {
  return (1 / rrs(gamma, mach)) * Math.sqrt(1 / tts(gamma, mach)) / mach;
}

function prandtlMeyer(gamma, mach) {
  const a = Math.sqrt((gamma + 1) / (gamma - 1));
  const b = Math.sqrt(((gamma - 1) / (gamma + 1)) * (mach * mach - 1));
  return (a * Math.atan(b) - Math.atan(Math.sqrt(mach * mach - 1))) * RAD;
}

function normalMachAfterShock(gamma, mach1) {
  return Math.sqrt((1 + 0.5 * (gamma - 1) * mach1 * mach1) / (gamma * mach1 * mach1 - 0.5 * (gamma - 1)));
}

function solveBisection(fn, low, high, tolerance = 1e-10, maxIterations = 200) {
  let fLow = fn(low);
  let fHigh = fn(high);

  if (Math.abs(fLow) < tolerance) return low;
  if (Math.abs(fHigh) < tolerance) return high;
  if (fLow * fHigh > 0) {
    throw new Error("Could not bracket a solution for this input.");
  }

  for (let i = 0; i < maxIterations; i += 1) {
    const mid = 0.5 * (low + high);
    const fMid = fn(mid);
    if (Math.abs(fMid) < tolerance || Math.abs(high - low) < tolerance) return mid;
    if (fLow * fMid <= 0) {
      high = mid;
      fHigh = fMid;
    } else {
      low = mid;
      fLow = fMid;
    }
  }

  return 0.5 * (low + high);
}

function machFromArea(gamma, area, branch) {
  if (area < 1) throw new Error("A/A* must be greater than or equal to 1.");
  if (Math.abs(area - 1) < 1e-12) return 1;

  if (branch === "sub") {
    return solveBisection((m) => areaRatio(gamma, m) - area, 1e-8, 0.999999);
  }
  return solveBisection((m) => areaRatio(gamma, m) - area, 1.000001, 100);
}

function machFromPrandtlMeyer(gamma, angleDegrees) {
  const maxAngle = (Math.sqrt((gamma + 1) / (gamma - 1)) - 1) * 90;
  if (angleDegrees <= 0 || angleDegrees >= maxAngle) {
    throw new Error(`Prandtl-Meyer angle must be between 0 and ${format(maxAngle)} degrees.`);
  }

  return solveBisection((m) => prandtlMeyer(gamma, m) - angleDegrees, 1.000001, 100);
}

function isentropicFromInput(gamma, input, value) {
  if (gamma <= 1) throw new Error("Gamma must be greater than 1.");
  if (!Number.isFinite(value)) throw new Error("Input value must be numeric.");

  let mach;
  if (input === "mach") {
    if (value <= 0) throw new Error("Mach number must be greater than 0.");
    mach = value;
  } else if (input === "tt0") {
    if (value <= 0 || value >= 1) throw new Error("T/T0 must be between 0 and 1.");
    mach = Math.sqrt((2 * (1 / value - 1)) / (gamma - 1));
  } else if (input === "pp0") {
    if (value <= 0 || value >= 1) throw new Error("p/p0 must be between 0 and 1.");
    mach = Math.sqrt((2 * (Math.pow(1 / value, (gamma - 1) / gamma) - 1)) / (gamma - 1));
  } else if (input === "rr0") {
    if (value <= 0 || value >= 1) throw new Error("rho/rho0 must be between 0 and 1.");
    mach = Math.sqrt((2 * (Math.pow(1 / value, gamma - 1) - 1)) / (gamma - 1));
  } else if (input === "area-sub") {
    mach = machFromArea(gamma, value, "sub");
  } else if (input === "area-sup") {
    mach = machFromArea(gamma, value, "sup");
  } else if (input === "mach-angle") {
    if (value <= 0 || value >= 90) throw new Error("Mach angle must be between 0 and 90 degrees.");
    mach = 1 / Math.sin(value * DEG);
  } else if (input === "pm-angle") {
    mach = machFromPrandtlMeyer(gamma, value);
  }

  return {
    "Mach number": mach,
    "Mach angle": mach >= 1 ? Math.asin(1 / mach) * RAD : null,
    "P-M angle": mach >= 1 ? prandtlMeyer(gamma, mach) : null,
    "T/T<sub>0</sub>": tt0(gamma, mach),
    "p/p<sub>0</sub>": pp0(gamma, mach),
    "ρ/ρ<sub>0</sub>": rr0(gamma, mach),
    "T/T<sup>*</sup>": tts(gamma, mach),
    "p/p<sup>*</sup>": pps(gamma, mach),
    "ρ/ρ<sup>*</sup>": rrs(gamma, mach),
    "A/A<sup>*</sup>": areaRatio(gamma, mach),
  };
}

function normalShockFromInput(gamma, input, value) {
  if (gamma <= 1) throw new Error("Gamma must be greater than 1.");
  if (!Number.isFinite(value)) throw new Error("Input value must be numeric.");

  let mach1;
  if (input === "m1") {
    if (value <= 1) throw new Error("M1 must be greater than 1.");
    mach1 = value;
  } else if (input === "m2") {
    const minM2 = Math.sqrt((gamma - 1) / (2 * gamma));
    if (value <= minM2 || value >= 1) throw new Error(`M2 must be between ${format(minM2)} and 1.`);
    mach1 = Math.sqrt((1 + 0.5 * (gamma - 1) * value * value) / (gamma * value * value - 0.5 * (gamma - 1)));
  } else if (input === "p2p1") {
    if (value <= 1) throw new Error("p2/p1 must be greater than 1.");
    mach1 = Math.sqrt(((value - 1) * (gamma + 1)) / (2 * gamma) + 1);
  } else if (input === "r2r1") {
    const maxDensityRatio = (gamma + 1) / (gamma - 1);
    if (value <= 1 || value >= maxDensityRatio) {
      throw new Error(`rho2/rho1 must be between 1 and ${format(maxDensityRatio)}.`);
    }
    mach1 = Math.sqrt((2 * value) / (gamma + 1 - value * (gamma - 1)));
  } else if (input === "t2t1") {
    if (value <= 1) throw new Error("T2/T1 must be greater than 1.");
    mach1 = solveBisection((m) => normalShockRatios(gamma, m).t2t1 - value, 1.000001, 100);
  } else if (input === "p02p01") {
    if (value <= 0 || value >= 1) throw new Error("p02/p01 must be between 0 and 1.");
    mach1 = solveBisection((m) => normalShockRatios(gamma, m).p02p01 - value, 1.000001, 100);
  } else if (input === "p1p02") {
    const maxValue = Math.pow((gamma + 1) / 2, -gamma / (gamma - 1));
    if (value <= 0 || value >= maxValue) throw new Error(`p1/p02 must be between 0 and ${format(maxValue)}.`);
    mach1 = solveBisection((m) => normalShockRatios(gamma, m).p1p02 - value, 1.000001, 100);
  }

  return normalShockRatios(gamma, mach1);
}

function normalShockRatios(gamma, mach1) {
  const mach2 = normalMachAfterShock(gamma, mach1);
  const p2p1 = 1 + (2 * gamma / (gamma + 1)) * (mach1 * mach1 - 1);
  const p02p01 = (pp0(gamma, mach1) / pp0(gamma, mach2)) * p2p1;

  return {
    "M<sub>1</sub>": mach1,
    "M<sub>2</sub>": mach2,
    "p<sub>2</sub>/p<sub>1</sub>": p2p1,
    "ρ<sub>2</sub>/ρ<sub>1</sub>": (rr0(gamma, mach2) / rr0(gamma, mach1)) * p02p01,
    "T<sub>2</sub>/T<sub>1</sub>": tt0(gamma, mach2) / tt0(gamma, mach1),
    "p<sub>02</sub>/p<sub>01</sub>": p02p01,
    "p<sub>1</sub>/p<sub>02</sub>": pp0(gamma, mach1) / p02p01,
    t2t1: tt0(gamma, mach2) / tt0(gamma, mach1),
    p02p01,
    p1p02: pp0(gamma, mach1) / p02p01,
  };
}

function normalShockDimensional(gamma, mach1, pressure1Pa, temperature1) {
  if (gamma <= 1) throw new Error("Gamma must be greater than 1.");
  if (mach1 <= 1 || !Number.isFinite(mach1)) throw new Error("M1 must be greater than 1.");
  if (pressure1Pa <= 0 || !Number.isFinite(pressure1Pa)) throw new Error("p1 must be greater than 0.");
  if (temperature1 <= 0 || !Number.isFinite(temperature1)) throw new Error("T1 must be greater than 0.");

  const ratios = normalShockRatios(gamma, mach1);
  const pressure2Pa = pressure1Pa * ratios["p<sub>2</sub>/p<sub>1</sub>"];
  const temperature2 = temperature1 * ratios.t2t1;
  const stagnationTemperature1 = temperature1 / tt0(gamma, mach1);
  const stagnationPressure1Pa = pressure1Pa / pp0(gamma, mach1);
  const stagnationPressure2Pa = pressure1Pa / ratios.p1p02;
  const stagnationTemperature2 = temperature2 / tt0(gamma, ratios["M<sub>2</sub>"]);
  const airGasConstant = gasConstantFromMolecularWeight(GAS_PROPERTIES.air.molecularWeight);
  const density1 = pressure1Pa / (airGasConstant * temperature1);

  return {
    __rowBreakBefore: ["M<sub>2</sub>"],
    "M<sub>1</sub>": mach1,
    "ρ<sub>1</sub>, kg/m<sup>3</sup>": density1,
    "p<sub>01</sub>, Pa": stagnationPressure1Pa,
    "T<sub>01</sub>, K": stagnationTemperature1,
    "M<sub>2</sub>": ratios["M<sub>2</sub>"],
    "p<sub>2</sub>, Pa": pressure2Pa,
    "T<sub>2</sub>, K": temperature2,
    "p<sub>02</sub>, Pa": stagnationPressure2Pa,
    "T<sub>02</sub>, K": stagnationTemperature2,
    "ρ<sub>2</sub>, kg/m<sup>3</sup>": pressure2Pa / (airGasConstant * temperature2),
  };
}

function obliqueThetaFromBeta(gamma, mach1, beta) {
  const sinBeta = Math.sin(beta);
  const numerator = 2 * (1 / Math.tan(beta)) * (mach1 * mach1 * sinBeta * sinBeta - 1);
  const denominator = mach1 * mach1 * (gamma + Math.cos(2 * beta)) + 2;
  return Math.atan(numerator / denominator);
}

function obliqueMaxTheta(gamma, mach1) {
  const betaMin = Math.asin(1 / mach1) + 1e-6;
  const betaMax = Math.PI / 2 - 1e-6;
  let bestBeta = betaMin;
  let bestTheta = 0;

  for (let i = 0; i <= 900; i += 1) {
    const beta = betaMin + (i / 900) * (betaMax - betaMin);
    const theta = obliqueThetaFromBeta(gamma, mach1, beta);
    if (theta > bestTheta) {
      bestTheta = theta;
      bestBeta = beta;
    }
  }

  return { theta: bestTheta, beta: bestBeta };
}

function solveObliqueBeta(gamma, mach1, theta, branch) {
  const maxPoint = obliqueMaxTheta(gamma, mach1);
  if (theta <= 0 || theta >= maxPoint.theta) return null;

  const betaMin = Math.asin(1 / mach1) + 1e-7;
  const betaMax = Math.PI / 2 - 1e-7;
  const low = branch === "weak" ? betaMin : maxPoint.beta;
  const high = branch === "weak" ? maxPoint.beta : betaMax;
  return solveBisection((beta) => obliqueThetaFromBeta(gamma, mach1, beta) - theta, low, high);
}

function obliqueShockSolution(gamma, mach1, thetaDegrees, branch) {
  const theta = thetaDegrees * DEG;
  const beta = solveObliqueBeta(gamma, mach1, theta, branch);
  if (!beta) return null;

  const m1n = mach1 * Math.sin(beta);
  const ratios = normalShockRatios(gamma, m1n);
  const m2n = ratios["M<sub>2</sub>"];
  const m2 = m2n / Math.sin(beta - theta);

  return {
    beta,
    theta,
    m1n,
    m2n,
    m2,
    ratios,
  };
}

function obliqueShockFromInput(gamma, mach1, thetaDegrees) {
  if (gamma <= 1) throw new Error("Gamma must be greater than 1.");
  if (mach1 <= 1 || !Number.isFinite(mach1)) throw new Error("M1 must be greater than 1.");
  if (thetaDegrees <= 0 || !Number.isFinite(thetaDegrees)) throw new Error("Turning angle must be greater than 0 degrees.");

  const maxPoint = obliqueMaxTheta(gamma, mach1);
  const maxThetaDegrees = maxPoint.theta * RAD;
  if (thetaDegrees >= maxThetaDegrees) {
    throw new Error(`No attached oblique shock. θ must be less than ${format(maxThetaDegrees)} degrees for M1 = ${format(mach1)}.`);
  }

  const weak = obliqueShockSolution(gamma, mach1, thetaDegrees, "weak");
  const strong = obliqueShockSolution(gamma, mach1, thetaDegrees, "strong");

  return {
    "M<sub>1</sub>": mach1,
    "θ (turn angle), deg.": thetaDegrees,
    "θ<sub>max</sub>, deg.": maxThetaDegrees,
    "β weak, deg.": weak.beta * RAD,
    "M<sub>2</sub> weak": weak.m2,
    "p<sub>2</sub>/p<sub>1</sub> weak": weak.ratios["p<sub>2</sub>/p<sub>1</sub>"],
    "T<sub>2</sub>/T<sub>1</sub> weak": weak.ratios["T<sub>2</sub>/T<sub>1</sub>"],
    "β strong, deg.": strong.beta * RAD,
    "M<sub>2</sub> strong": strong.m2,
    "p<sub>2</sub>/p<sub>1</sub> strong": strong.ratios["p<sub>2</sub>/p<sub>1</sub>"],
    "T<sub>2</sub>/T<sub>1</sub> strong": strong.ratios["T<sub>2</sub>/T<sub>1</sub>"],
  };
}

function obliqueShockDimensional(gamma, mach1, pressure1Pa, temperature1, thetaDegrees) {
  if (pressure1Pa <= 0 || !Number.isFinite(pressure1Pa)) throw new Error("p1 must be greater than 0.");
  if (temperature1 <= 0 || !Number.isFinite(temperature1)) throw new Error("T1 must be greater than 0.");

  const base = obliqueShockFromInput(gamma, mach1, thetaDegrees);
  const weak = obliqueShockSolution(gamma, mach1, thetaDegrees, "weak");
  const strong = obliqueShockSolution(gamma, mach1, thetaDegrees, "strong");
  const airGasConstant = gasConstantFromMolecularWeight(GAS_PROPERTIES.air.molecularWeight);

  const dimensionalBranch = (solution) => {
    const pressure2Pa = pressure1Pa * solution.ratios["p<sub>2</sub>/p<sub>1</sub>"];
    const temperature2 = temperature1 * solution.ratios["T<sub>2</sub>/T<sub>1</sub>"];
    return {
      pressure2Pa,
      temperature2,
      stagnationPressure2Pa: pressure2Pa / pp0(gamma, solution.m2),
      density2: pressure2Pa / (airGasConstant * temperature2),
    };
  };

  const weakDimensional = dimensionalBranch(weak);
  const strongDimensional = dimensionalBranch(strong);

  return {
    "M<sub>1</sub>": mach1,
    "θ (turn angle), deg.": thetaDegrees,
    "θ<sub>max</sub>, deg.": base["θ<sub>max</sub>, deg."],
    "β weak, deg.": weak.beta * RAD,
    "M<sub>2</sub> weak": weak.m2,
    "p<sub>2</sub> weak, Pa": weakDimensional.pressure2Pa,
    "T<sub>2</sub> weak, K": weakDimensional.temperature2,
    "p<sub>02</sub> weak, Pa": weakDimensional.stagnationPressure2Pa,
    "ρ<sub>2</sub> weak, kg/m<sup>3</sup>": weakDimensional.density2,
    "β strong, deg.": strong.beta * RAD,
    "M<sub>2</sub> strong": strong.m2,
    "p<sub>2</sub> strong, Pa": strongDimensional.pressure2Pa,
    "T<sub>2</sub> strong, K": strongDimensional.temperature2,
    "p<sub>02</sub> strong, Pa": strongDimensional.stagnationPressure2Pa,
    "ρ<sub>2</sub> strong, kg/m<sup>3</sup>": strongDimensional.density2,
  };
}

function nondimensionalVelocity(gamma, mach) {
  return Math.sqrt(((gamma - 1) * mach * mach) / (2 + (gamma - 1) * mach * mach));
}

function machFromNondimensionalVelocity(gamma, velocityMagnitude) {
  const velocitySquared = velocityMagnitude * velocityMagnitude;
  if (velocitySquared <= 0 || velocitySquared >= 1) {
    throw new Error("Taylor-Maccoll integration reached an invalid velocity state.");
  }
  return Math.sqrt((2 * velocitySquared) / ((gamma - 1) * (1 - velocitySquared)));
}

function taylorMaccollDerivative(gamma, theta, radialVelocity, polarVelocity) {
  const velocitySquared = radialVelocity * radialVelocity + polarVelocity * polarVelocity;
  const thermalTerm = ((gamma - 1) / 2) * (1 - velocitySquared);
  const denominator = thermalTerm - polarVelocity * polarVelocity;
  const sinTheta = Math.sin(theta);

  if (
    !Number.isFinite(thermalTerm)
    || !Number.isFinite(denominator)
    || Math.abs(denominator) < 1e-10
    || Math.abs(sinTheta) < 1e-10
  ) {
    throw new Error("Taylor-Maccoll integration failed for this cone angle.");
  }

  const cotTheta = Math.cos(theta) / sinTheta;
  const numerator = polarVelocity * polarVelocity * radialVelocity
    - thermalTerm * (2 * radialVelocity + polarVelocity * cotTheta);

  return {
    radial: polarVelocity,
    polar: numerator / denominator,
  };
}

function integrateTaylorMaccoll(gamma, beta, coneAngle, radialVelocity, polarVelocity) {
  const steps = Math.max(24, Math.ceil(Math.abs(beta - coneAngle) / 0.00035));
  const step = (coneAngle - beta) / steps;
  let theta = beta;
  let vr = radialVelocity;
  let vt = polarVelocity;

  for (let i = 0; i < steps; i += 1) {
    const k1 = taylorMaccollDerivative(gamma, theta, vr, vt);
    const k2 = taylorMaccollDerivative(
      gamma,
      theta + step / 2,
      vr + (step * k1.radial) / 2,
      vt + (step * k1.polar) / 2,
    );
    const k3 = taylorMaccollDerivative(
      gamma,
      theta + step / 2,
      vr + (step * k2.radial) / 2,
      vt + (step * k2.polar) / 2,
    );
    const k4 = taylorMaccollDerivative(
      gamma,
      theta + step,
      vr + step * k3.radial,
      vt + step * k3.polar,
    );

    vr += (step / 6) * (k1.radial + 2 * k2.radial + 2 * k3.radial + k4.radial);
    vt += (step / 6) * (k1.polar + 2 * k2.polar + 2 * k3.polar + k4.polar);
    theta += step;

    if (!Number.isFinite(vr) || !Number.isFinite(vt) || vr * vr + vt * vt >= 1) {
      throw new Error("Taylor-Maccoll integration failed for this cone angle.");
    }
  }

  return { radialVelocity: vr, polarVelocity: vt };
}

function conicalShockStateAtBeta(gamma, mach1, beta, coneAngle) {
  const velocity1 = nondimensionalVelocity(gamma, mach1);
  const normalMach1 = mach1 * Math.sin(beta);
  if (normalMach1 <= 1) throw new Error("Shock angle must produce a supersonic normal Mach number.");

  const shockRatios = normalShockRatios(gamma, normalMach1);
  const densityRatio = shockRatios["ρ<sub>2</sub>/ρ<sub>1</sub>"];
  const radialVelocity = velocity1 * Math.cos(beta);
  const polarVelocity = -velocity1 * Math.sin(beta) / densityRatio;
  const surface = integrateTaylorMaccoll(gamma, beta, coneAngle, radialVelocity, polarVelocity);

  return {
    radialVelocity,
    polarVelocity,
    surface,
    shockRatios,
  };
}

function coneShockBeta(gamma, mach1, coneAngle) {
  const machAngle = Math.asin(1 / mach1);
  const betaLow = Math.max(machAngle, coneAngle) + 1e-5;
  const betaHigh = Math.PI / 2 - 1e-5;

  let previousBeta = null;
  let previousValue = null;

  for (let i = 0; i <= 260; i += 1) {
    const beta = betaLow + (i / 260) * (betaHigh - betaLow);
    let value;

    try {
      value = conicalShockStateAtBeta(gamma, mach1, beta, coneAngle).surface.polarVelocity;
    } catch (exception) {
      continue;
    }

    if (Math.abs(value) < 1e-8) return beta;
    if (previousValue !== null && previousValue * value < 0) {
      return solveBisection(
        (candidate) => conicalShockStateAtBeta(gamma, mach1, candidate, coneAngle).surface.polarVelocity,
        previousBeta,
        beta,
        1e-9,
        120,
      );
    }

    previousBeta = beta;
    previousValue = value;
  }

  return null;
}

function coneShockFromInput(gamma, mach1, coneAngleDegrees) {
  if (gamma <= 1) throw new Error("Gamma must be greater than 1.");
  if (mach1 <= 1 || !Number.isFinite(mach1)) throw new Error("M1 must be greater than 1.");
  if (coneAngleDegrees <= 0 || !Number.isFinite(coneAngleDegrees)) {
    throw new Error("Cone half-angle must be greater than 0 degrees.");
  }
  if (coneAngleDegrees >= 89) throw new Error("Cone half-angle must be less than 89 degrees.");

  const coneAngle = coneAngleDegrees * DEG;
  const beta = coneShockBeta(gamma, mach1, coneAngle);
  if (!beta) {
    throw new Error(`No attached conical shock found for M1 = ${format(mach1)} and cone angle = ${format(coneAngleDegrees)} degrees.`);
  }
  if (beta <= coneAngle) {
    throw new Error("Invalid conical shock solution: shock angle must be greater than the cone half-angle.");
  }

  const state = conicalShockStateAtBeta(gamma, mach1, beta, coneAngle);
  const shockVelocity = Math.hypot(state.radialVelocity, state.polarVelocity);
  const surfaceVelocity = Math.hypot(state.surface.radialVelocity, state.surface.polarVelocity);
  const shockMach = machFromNondimensionalVelocity(gamma, shockVelocity);
  const surfaceMach = machFromNondimensionalVelocity(gamma, surfaceVelocity);
  const shockTemperatureRatio = state.shockRatios["T<sub>2</sub>/T<sub>1</sub>"];
  const shockPressureRatio = state.shockRatios["p<sub>2</sub>/p<sub>1</sub>"];
  const shockDensityRatio = state.shockRatios["ρ<sub>2</sub>/ρ<sub>1</sub>"];
  const surfaceToShockTemperatureRatio = (1 - surfaceVelocity * surfaceVelocity)
    / (1 - shockVelocity * shockVelocity);
  const surfaceToShockPressureRatio = Math.pow(
    surfaceToShockTemperatureRatio,
    gamma / (gamma - 1),
  );
  const surfaceToShockDensityRatio = Math.pow(
    surfaceToShockTemperatureRatio,
    1 / (gamma - 1),
  );

  return {
    __sectionBefore: {
      "M<sub>n,1</sub>": "Behind shock",
      "M<sub>c</sub>": "Cone surface",
    },
    "M<sub>1</sub>": mach1,
    "θ<sub>c</sub> cone half-angle, deg.": coneAngleDegrees,
    "β shock angle from axis, deg.": beta * RAD,
    "M<sub>n,1</sub>": mach1 * Math.sin(beta),
    "M<sub>2</sub>": shockMach,
    "p<sub>2</sub>/p<sub>1</sub>": shockPressureRatio,
    "T<sub>2</sub>/T<sub>1</sub>": shockTemperatureRatio,
    "ρ<sub>2</sub>/ρ<sub>1</sub>": shockDensityRatio,
    "M<sub>c</sub>": surfaceMach,
    "p<sub>c</sub>/p<sub>1</sub>": shockPressureRatio * surfaceToShockPressureRatio,
    "T<sub>c</sub>/T<sub>1</sub>": shockTemperatureRatio * surfaceToShockTemperatureRatio,
    "ρ<sub>c</sub>/ρ<sub>1</sub>": shockDensityRatio * surfaceToShockDensityRatio,
  };
}

function speedOfSoundFromInput(gamma, gasConstant, temperature, gas = "custom", molecularWeight = null) {
  if (gamma <= 1) throw new Error("Gamma must be greater than 1.");
  if (gasConstant <= 0 || !Number.isFinite(gasConstant)) throw new Error("R must be greater than 0.");
  if (temperature <= 0 || !Number.isFinite(temperature)) throw new Error("Temperature must be greater than 0 K.");

  const speed = Math.sqrt(gamma * gasConstant * temperature);

  return {
    Gas: GAS_PROPERTIES[gas]?.label ?? "Custom",
    "γ": gamma,
    ...(molecularWeight ? { "MW, kg/kmol": molecularWeight } : {}),
    "R, J/(kg K)": gasConstant,
    "T, K": temperature,
    "a, m/s": speed,
    "a, ft/s": speed * 3.280839895,
  };
}

function safeExp(value) {
  return Math.exp(Math.max(-700, Math.min(700, value)));
}

function interpolateNullable(low, high, fraction, index) {
  const lowValue = low[index];
  const highValue = high[index];
  if (lowValue === null && highValue === null) return null;
  if (lowValue === null) return highValue;
  if (highValue === null) return lowValue;
  return lowValue + fraction * (highValue - lowValue);
}

function hansenCollisionData(temperature) {
  const table = HANSEN_COLLISION_TABLE;
  if (temperature <= table[0][0]) return table[0];
  if (temperature >= table[table.length - 1][0]) return table[table.length - 1];

  for (let i = 0; i < table.length - 1; i += 1) {
    const low = table[i];
    const high = table[i + 1];
    if (temperature >= low[0] && temperature <= high[0]) {
      const fraction = (temperature - low[0]) / (high[0] - low[0]);
      return low.map((value, index) => (
        index === 0 ? temperature : interpolateNullable(low, high, fraction, index)
      ));
    }
  }

  return table[table.length - 1];
}

function hansenSpeciesClass(species) {
  if (species === "N2" || species === "O2") return "molecule";
  if (species === "e-") return "electron";
  if (species === "N+" || species === "O+") return "ion";
  return "atom";
}

function hansenCoulombLog(temperature, pressureAtm, electronMoleFraction, z) {
  const electronDensityM3 = electronMoleFraction * (pressureAtm * ATM_TO_PA)
    / (z * HANSEN_R_UNIVERSAL * temperature);
  if (electronDensityM3 <= 0 || !Number.isFinite(electronDensityM3)) return 10;

  const electronDensityCm3 = electronDensityM3 / 1e6;
  const temperatureEv = Math.max(temperature / 11604.51812, 1e-6);
  const raw = temperatureEv < 10
    ? 23.5 - Math.log(Math.sqrt(electronDensityCm3) * Math.pow(temperatureEv, -1.25))
      - Math.sqrt(1e-5 + Math.pow(Math.log(temperatureEv) - 2, 2) / 16)
    : 24 - Math.log(Math.sqrt(electronDensityCm3) / temperatureEv);
  return Math.max(2, Math.min(30, raw));
}

function hansenCrossSectionRatio(speciesI, speciesJ, collisionData, coulombLog, mode = "viscosity") {
  const classI = hansenSpeciesClass(speciesI);
  const classJ = hansenSpeciesClass(speciesJ);
  const classes = [classI, classJ];

  if (mode === "conductivity") {
    if (classI === "electron" || classJ === "electron") return 0;
    if (classI === "molecule" && classJ === "molecule") return 1;
    if (classes.includes("molecule")) return collisionData[6];
    return collisionData[7];
  }
  if (mode === "reaction") {
    if (classI === "electron" && classJ === "electron") return (collisionData[5] ?? 12.24) * coulombLog;
    if (classes.includes("electron")) return collisionData[4] ?? 0.243;
    if (classI === "molecule" && classJ === "molecule") return 1;
    if (classes.includes("molecule")) return collisionData[6];
    return collisionData[7];
  }

  if (classI === "electron" && classJ === "electron") {
    return (collisionData[5] ?? 12.24) * coulombLog;
  }
  if (classes.includes("electron")) {
    const otherClass = classI === "electron" ? classJ : classI;
    if (otherClass === "ion") return collisionData[4] ?? 0.243;
    return 0;
  }
  if (classI === "molecule" && classJ === "molecule") return 1;
  if (classes.includes("molecule")) return collisionData[2];
  return collisionData[3];
}

function lnQp(species, temperature) {
  const T = temperature;
  if (T <= 0 || !Number.isFinite(T)) throw new Error("Temperature must be greater than 0 K.");

  if (species === "N2") {
    return 3.5 * Math.log(T) - 0.42 - Math.log(1 - safeExp(-3390 / T));
  }
  if (species === "O2") {
    return 3.5 * Math.log(T) + 0.11
      - Math.log(1 - safeExp(-2270 / T))
      + Math.log(3 + 2 * safeExp(-11390 / T) + safeExp(-18990 / T));
  }
  if (species === "O") {
    return 2.5 * Math.log(T) + 0.50
      + Math.log(5 + 3 * safeExp(-228 / T) + safeExp(-326 / T)
        + 5 * safeExp(-22800 / T) + safeExp(-48600 / T));
  }
  if (species === "N") {
    return 2.5 * Math.log(T) + 0.30
      + Math.log(4 + 10 * safeExp(-27700 / T) + 6 * safeExp(-41500 / T));
  }
  if (species === "O+") {
    return 2.5 * Math.log(T) + 0.50
      + Math.log(4 + 10 * safeExp(-38600 / T) + 6 * safeExp(-58200 / T));
  }
  if (species === "N+") {
    return 2.5 * Math.log(T) + 0.30
      + Math.log(1 + 3 * safeExp(-70.6 / T) + 5 * safeExp(-188.9 / T)
        + 5 * safeExp(-22000 / T) + safeExp(-47000 / T) + 5 * safeExp(-67900 / T));
  }
  if (species === "e-") return 2.5 * Math.log(T) - 14.24;

  throw new Error(`Unknown Hansen species: ${species}.`);
}

function hansenEquilibriumConstants(temperature, xi = HANSEN_DEFAULT_XI) {
  const T = temperature;
  const { xN2Feed, xO2Feed } = hansenFeedFractions(xi);
  const kpO2 = safeExp(-59000 / T + 2 * lnQp("O", T) - lnQp("O2", T));
  const kpN2 = safeExp(-113200 / T + 2 * lnQp("N", T) - lnQp("N2", T));
  const kpOIon = safeExp(-158000 / T + lnQp("O+", T) + lnQp("e-", T) - lnQp("O", T));
  const kpNIon = safeExp(-168800 / T + lnQp("N+", T) + lnQp("e-", T) - lnQp("N", T));
  return { kpO2, kpN2, kpIon: xO2Feed * kpOIon + xN2Feed * kpNIon };
}

function hansenEpsilons(temperature, pressureAtm, xi = HANSEN_DEFAULT_XI) {
  const { kpO2, kpN2, kpIon } = hansenEquilibriumConstants(temperature, xi);
  const p = pressureAtm;
  const a1 = 1 + (4 * p) / kpO2;
  const eps1 = Number.isFinite(a1) ? (-0.8 + Math.sqrt(0.64 + 0.8 * a1)) / (2 * a1) : 0;
  const a2 = 1 + (4 * p) / kpN2;
  const eps2 = Number.isFinite(a2) ? (-0.4 + Math.sqrt(0.16 + 3.84 * a2)) / (2 * a2) : 0;
  const ionDenominator = 1 + p / kpIon;
  const eps3 = Number.isFinite(ionDenominator) ? Math.pow(ionDenominator, -0.5) : 0;
  return { eps1, eps2, eps3 };
}

function hansenFeedFractions(xi) {
  if (xi <= 0 || !Number.isFinite(xi)) throw new Error("O/N atomic ratio must be positive.");
  return { xN2Feed: 1 / (1 + xi), xO2Feed: xi / (1 + xi) };
}

function hansenMixtureGasConstant(xi) {
  const { xN2Feed, xO2Feed } = hansenFeedFractions(xi);
  const molecularMassKg = (xN2Feed * HANSEN_MOLAR_MASS[0] + xO2Feed * HANSEN_MOLAR_MASS[1]) / 1000;
  return HANSEN_R_UNIVERSAL / molecularMassKg;
}

function hansenFeedMolecularWeight(xi) {
  const { xN2Feed, xO2Feed } = hansenFeedFractions(xi);
  return xN2Feed * HANSEN_MOLAR_MASS[0] + xO2Feed * HANSEN_MOLAR_MASS[1];
}

function hansenXiFromSpeciesInput(data) {
  const speciesFeed = {
    n2: Number(data.get("n2Feed") ?? 0),
    o2: Number(data.get("o2Feed") ?? 0),
    n: Number(data.get("nFeed") ?? 0),
    o: Number(data.get("oFeed") ?? 0),
    nIon: Number(data.get("nIonFeed") ?? 0),
    oIon: Number(data.get("oIonFeed") ?? 0),
    electron: Number(data.get("electronFeed") ?? 0),
  };
  if (Object.values(speciesFeed).some((value) => !Number.isFinite(value) || value < 0)) {
    throw new Error("Starting species concentrations must be nonnegative numeric values.");
  }

  const nitrogenAtoms = 2 * speciesFeed.n2 + speciesFeed.n + speciesFeed.nIon;
  const oxygenAtoms = 2 * speciesFeed.o2 + speciesFeed.o + speciesFeed.oIon;
  if (nitrogenAtoms <= 0 || oxygenAtoms <= 0) {
    throw new Error("Starting species must include at least one nitrogen-bearing and one oxygen-bearing species.");
  }
  return oxygenAtoms / nitrogenAtoms;
}

function hansenAir(temperature, pressureAtm, xi = HANSEN_DEFAULT_XI) {
  if (pressureAtm <= 0 || !Number.isFinite(pressureAtm)) throw new Error("Pressure must be greater than 0 atm.");

  const { xN2Feed, xO2Feed } = hansenFeedFractions(xi);
  if (temperature <= 1200) {
    return {
      Z: 1,
      eps1: 0,
      eps2: 0,
      eps3: 0,
      xN2: xN2Feed,
      xO2: xO2Feed,
      xN: 0,
      xO: 0,
      "xO+": 0,
      "xN+": 0,
      "xe-": 0,
    };
  }

  const { eps1: eps1Air, eps2: eps2Air, eps3 } = hansenEpsilons(temperature, pressureAtm, xi);
  const eps1 = Math.min((eps1Air * xO2Feed) / 0.2, xO2Feed);
  const eps2 = Math.min((eps2Air * xN2Feed) / 0.8, xN2Feed);
  const nShare = xN2Feed / (xN2Feed + xO2Feed);
  const oShare = xO2Feed / (xN2Feed + xO2Feed);
  const z = 1 + eps1 + eps2 + 2 * eps3;
  const xe = (2 * eps3) / z;

  return {
    Z: z,
    eps1,
    eps2,
    eps3,
    xN2: Math.max(0, (xN2Feed - eps2) / z),
    xO2: Math.max(0, (xO2Feed - eps1) / z),
    xN: Math.max(0, (2 * eps2 - 2 * nShare * eps3) / z),
    xO: Math.max(0, (2 * eps1 - 2 * oShare * eps3) / z),
    "xO+": Math.max(0, oShare * xe),
    "xN+": Math.max(0, nShare * xe),
    "xe-": Math.max(0, xe),
  };
}

function dlnQpDT(species, temperature) {
  const step = 1e-3 * temperature;
  return (lnQp(species, temperature + step) - lnQp(species, temperature - step)) / (2 * step);
}

function speciesHOverRT(species, temperature) {
  return temperature * dlnQpDT(species, temperature) + HANSEN_FORMATION_OVER_R[species] / temperature;
}

function speciesCvOverR(species, temperature) {
  const step = 1e-3 * temperature;
  const energy = (T) => T * (speciesHOverRT(species, T) - 1);
  return (energy(temperature + step) - energy(temperature - step)) / (2 * step);
}

function mixtureZHOverRT(temperature, pressureAtm, xi) {
  const state = hansenAir(temperature, pressureAtm, xi);
  return state.Z * (
    state.xO2 * speciesHOverRT("O2", temperature)
    + state.xN2 * speciesHOverRT("N2", temperature)
    + state.xO * speciesHOverRT("O", temperature)
    + state.xN * speciesHOverRT("N", temperature)
    + state["xO+"] * speciesHOverRT("O+", temperature)
    + state["xN+"] * speciesHOverRT("N+", temperature)
    + state["xe-"] * speciesHOverRT("e-", temperature)
  );
}

function mixtureZEOverRT(temperature, pressureAtm, xi) {
  const state = hansenAir(temperature, pressureAtm, xi);
  return state.Z * HANSEN_SPECIES.reduce((sum, species) => {
    const key = species === "e-" ? "xe-" : `x${species}`;
    return sum + state[key] * (speciesHOverRT(species, temperature) - 1);
  }, 0);
}

function pressureForConstantDensity(newTemperature, referenceTemperature, referencePressureAtm, xi) {
  const zRef = hansenAir(referenceTemperature, referencePressureAtm, xi).Z;
  const target = referencePressureAtm / (zRef * referenceTemperature);
  let pressure = referencePressureAtm * (newTemperature / referenceTemperature);

  for (let i = 0; i < 30; i += 1) {
    const z = hansenAir(newTemperature, pressure, xi).Z;
    const f = pressure / (z * newTemperature) - target;
    const dp = 1e-5 * Math.max(pressure, 1e-8);
    const zp = hansenAir(newTemperature, pressure + dp, xi).Z;
    const fp = (pressure + dp) / (zp * newTemperature) - target;
    const dfdp = (fp - f) / dp;
    if (dfdp === 0 || !Number.isFinite(dfdp)) break;
    pressure = Math.max(pressure - f / dfdp, 1e-12);
  }

  return pressure;
}

function hansenCpOverR(temperature, pressureAtm, xi) {
  const step = 1e-3 * temperature;
  const yp = mixtureZHOverRT(temperature + step, pressureAtm, xi);
  const ym = mixtureZHOverRT(temperature - step, pressureAtm, xi);
  return ((temperature + step) * yp - (temperature - step) * ym) / (2 * step);
}

function hansenCvOverR(temperature, pressureAtm, xi) {
  const step = 1e-3 * temperature;
  const pp = pressureForConstantDensity(temperature + step, temperature, pressureAtm, xi);
  const pm = pressureForConstantDensity(temperature - step, temperature, pressureAtm, xi);
  const xp = mixtureZEOverRT(temperature + step, pp, xi);
  const xm = mixtureZEOverRT(temperature - step, pm, xi);
  return ((temperature + step) * xp - (temperature - step) * xm) / (2 * step);
}

function hansenPhiFactor(temperature, pressureAtm, xi) {
  const z0 = hansenAir(temperature, pressureAtm, xi).Z;
  const rho0Scaled = pressureAtm / (z0 * temperature);
  const dp = 1e-4 * Math.max(pressureAtm, 1e-8);
  const p1 = pressureAtm + dp;
  const p2 = Math.max(pressureAtm - dp, 1e-12);
  const rho1Scaled = p1 / (hansenAir(temperature, p1, xi).Z * temperature);
  const rho2Scaled = p2 / (hansenAir(temperature, p2, xi).Z * temperature);
  const dpdrho = (p1 - p2) / (rho1Scaled - rho2Scaled);
  return (rho0Scaled / pressureAtm) * dpdrho;
}

function hansenMoleFractionsArray(temperature, pressurePa, xi) {
  const { xN2Feed, xO2Feed } = hansenFeedFractions(xi);
  if (temperature <= 1200) return [xN2Feed, xO2Feed, 0, 0, 0, 0, 0];

  const state = hansenAir(temperature, pressurePa / ATM_TO_PA, xi);
  const values = [
    state.xN2,
    state.xO2,
    state.xN,
    state.xO,
    state["xO+"],
    state["xN+"],
    state["xe-"],
  ].map((value) => Math.max(0, value));
  const sum = values.reduce((total, value) => total + value, 0);
  if (sum <= 0 || !Number.isFinite(sum)) throw new Error("Invalid reacting-air species fractions.");
  return values.map((value) => value / sum);
}

function hansenMassFractionsFromMoleFractions(moleFractions) {
  const mixtureMolarMass = moleFractions.reduce((sum, value, index) => (
    sum + value * HANSEN_MOLAR_MASS[index]
  ), 0);
  if (mixtureMolarMass <= 0 || !Number.isFinite(mixtureMolarMass)) {
    throw new Error("Invalid reacting-air mixture molecular weight.");
  }
  return moleFractions.map((value, index) => value * HANSEN_MOLAR_MASS[index] / mixtureMolarMass);
}

function hansenPartialPressure(pressurePa, temperature, xi) {
  const moleFractions = hansenMoleFractionsArray(temperature, pressurePa, xi);
  return {
    moleFractions,
    massFractions: hansenMassFractionsFromMoleFractions(moleFractions),
  };
}

function hansenMixtureGasConstantFromMassFractions(massFractions) {
  return massFractions.reduce((sum, value, index) => (
    sum + value * HANSEN_SPECIFIC_GAS_CONSTANTS[index]
  ), 0);
}

function hansenFrozenMixtureEnthalpy(temperature, massFractions) {
  const vibrationalN2 = 3389.82;
  const vibrationalO2 = 2771.09;
  return (
    massFractions[0] * ((7 / 2) * HANSEN_SPECIFIC_GAS_CONSTANTS[0] * temperature
      + (vibrationalN2 * HANSEN_SPECIFIC_GAS_CONSTANTS[0]) / (safeExp(vibrationalN2 / temperature) - 1))
    + massFractions[1] * ((7 / 2) * HANSEN_SPECIFIC_GAS_CONSTANTS[1] * temperature
      + (vibrationalO2 * HANSEN_SPECIFIC_GAS_CONSTANTS[1]) / (safeExp(vibrationalO2 / temperature) - 1))
    + massFractions[2] * (5 / 2) * HANSEN_SPECIFIC_GAS_CONSTANTS[2] * temperature
    + massFractions[3] * (5 / 2) * HANSEN_SPECIFIC_GAS_CONSTANTS[3] * temperature
    + massFractions[4] * (5 / 2) * HANSEN_SPECIFIC_GAS_CONSTANTS[4] * temperature
    + massFractions[5] * (5 / 2) * HANSEN_SPECIFIC_GAS_CONSTANTS[5] * temperature
    + massFractions[6] * (5 / 2) * HANSEN_SPECIFIC_GAS_CONSTANTS[6] * temperature
  );
}

function hansenEquilibriumEnthalpy(temperature, pressureAtm, xi) {
  const composition = hansenPartialPressure(pressureAtm * ATM_TO_PA, temperature, xi);
  if (temperature <= 1200) return hansenFrozenMixtureEnthalpy(temperature, composition.massFractions);
  return hansenMixtureGasConstant(xi) * temperature * mixtureZHOverRT(temperature, pressureAtm, xi);
}

function hansenEquilibriumState(temperature, pressurePa, xi) {
  const composition = hansenPartialPressure(pressurePa, temperature, xi);
  const gasConstant = hansenMixtureGasConstantFromMassFractions(composition.massFractions);
  const density = pressurePa / (gasConstant * temperature);
  const pressureAtm = pressurePa / ATM_TO_PA;
  const phiGamma = temperature <= 1200
    ? 1.4
    : hansenCpOverR(temperature, pressureAtm, xi) / hansenCvOverR(temperature, pressureAtm, xi)
      * hansenPhiFactor(temperature, pressureAtm, xi);
  const soundSpeed = Math.sqrt(phiGamma * pressurePa / density);

  return {
    ...composition,
    gasConstant,
    density,
    phiGamma,
    soundSpeed,
  };
}

function normalShockFull(mach1, temperature1, pressure1Pa, xi = HANSEN_DEFAULT_XI) {
  if (mach1 <= 1 || !Number.isFinite(mach1)) throw new Error("M1 must be greater than 1.");
  if (temperature1 <= 0 || !Number.isFinite(temperature1)) throw new Error("T1 must be greater than 0 K.");
  if (pressure1Pa <= 0 || !Number.isFinite(pressure1Pa)) throw new Error("p1 must be greater than 0 Pa.");

  const upstream = hansenEquilibriumState(temperature1, pressure1Pa, xi);
  const density1 = upstream.density;
  const velocityNormal = mach1 * upstream.soundSpeed;
  const gammaGuess = Math.max(upstream.phiGamma, 1.05);

  const densityRatioIdeal = ((gammaGuess + 1) * mach1 * mach1) / ((gammaGuess - 1) * mach1 * mach1 + 2);
  const pressureRatioIdeal = 1 + (2 * gammaGuess / (gammaGuess + 1)) * (mach1 * mach1 - 1);
  const temperatureRatioIdeal = pressureRatioIdeal / densityRatioIdeal;
  const h1 = hansenEquilibriumEnthalpy(temperature1, pressure1Pa / ATM_TO_PA, xi);
  const pRef = Math.max(pressure1Pa, 1);
  const hRef = Math.max(Math.abs(h1) + 0.5 * velocityNormal * velocityNormal, 1);
  const downstreamTemperatureGuess = temperature1 * temperatureRatioIdeal;
  const temperatureLowerBound = Math.max(1, Math.min(200, 0.75 * temperature1));
  const temperatureUpperBound = Math.max(60000, 8 * temperature1, 1.5 * downstreamTemperatureGuess);

  const bounds = {
    logTMin: Math.log(temperatureLowerBound),
    logTMax: Math.log(temperatureUpperBound),
    logRhoMin: Math.log(density1 * 1.0001),
    logRhoMax: Math.log(density1 * 100),
  };
  let x = [
    Math.log(Math.min(Math.max(downstreamTemperatureGuess, temperatureLowerBound * 1.01), temperatureUpperBound * 0.98)),
    Math.log(Math.max(density1 * densityRatioIdeal, density1 * 1.01)),
  ];

  const clampState = (state) => [
    Math.min(Math.max(state[0], bounds.logTMin), bounds.logTMax),
    Math.min(Math.max(state[1], bounds.logRhoMin), bounds.logRhoMax),
  ];
  const residuals = (state) => {
    const temperature2 = Math.exp(state[0]);
    const density2 = Math.exp(state[1]);
    const velocity2 = velocityNormal * density1 / density2;
    const pressure2Pa = pressure1Pa + density1 * velocityNormal * velocityNormal * (1 - density1 / density2);
    if (pressure2Pa <= 0 || !Number.isFinite(pressure2Pa)) return [1e6, 1e6];

    const composition = hansenPartialPressure(pressure2Pa, temperature2, xi);
    const gasConstant = hansenMixtureGasConstantFromMassFractions(composition.massFractions);
    const eosResidual = (pressure2Pa - density2 * gasConstant * temperature2) / pRef;
    const targetH2 = h1 + 0.5 * (velocityNormal * velocityNormal - velocity2 * velocity2);
    const h2 = hansenEquilibriumEnthalpy(temperature2, pressure2Pa / ATM_TO_PA, xi);
    const enthalpyResidual = (h2 - targetH2) / hRef;

    if (!Number.isFinite(eosResidual) || !Number.isFinite(enthalpyResidual)) return [1e6, 1e6];
    return [eosResidual, enthalpyResidual];
  };
  const norm = (values) => Math.hypot(values[0], values[1]);

  for (let iteration = 0; iteration < 80; iteration += 1) {
    const f = residuals(x);
    if (norm(f) < 1e-9) break;

    const jacobian = [[], []];
    for (let column = 0; column < 2; column += 1) {
      const step = 1e-5;
      const xp = [...x];
      xp[column] += step;
      const fp = residuals(clampState(xp));
      jacobian[0][column] = (fp[0] - f[0]) / step;
      jacobian[1][column] = (fp[1] - f[1]) / step;
    }

    const det = jacobian[0][0] * jacobian[1][1] - jacobian[0][1] * jacobian[1][0];
    if (Math.abs(det) < 1e-18 || !Number.isFinite(det)) break;

    const dx = [
      (jacobian[1][1] * f[0] - jacobian[0][1] * f[1]) / det,
      (-jacobian[1][0] * f[0] + jacobian[0][0] * f[1]) / det,
    ];
    let damping = 1;
    let accepted = false;
    const currentNorm = norm(f);

    for (let trial = 0; trial < 12; trial += 1) {
      const candidate = clampState([x[0] - damping * dx[0], x[1] - damping * dx[1]]);
      if (norm(residuals(candidate)) < currentNorm) {
        x = candidate;
        accepted = true;
        break;
      }
      damping *= 0.5;
    }
    if (!accepted) break;
  }

  const finalResiduals = residuals(x);
  if (norm(finalResiduals) > 1e-5) {
    throw new Error("Reacting normal shock solve did not converge for these inputs.");
  }

  const temperature2 = Math.exp(x[0]);
  const density2 = Math.exp(x[1]);
  const velocity2 = velocityNormal * density1 / density2;
  const pressure2Pa = pressure1Pa + density1 * velocityNormal * velocityNormal * (1 - density1 / density2);
  const composition = hansenPartialPressure(pressure2Pa, temperature2, xi);
  const downstream = hansenEquilibriumState(temperature2, pressure2Pa, xi);
  const mach2 = velocity2 / downstream.soundSpeed;
  const speciesDensities = composition.massFractions.map((massFraction) => massFraction * density2);

  return {
    mach1,
    mach2,
    temperature1,
    pressure1Pa,
    density1,
    velocityNormal,
    temperature2,
    density2,
    pressure2Pa,
    velocity2,
    speciesDensities,
    moleFractions: composition.moleFractions,
    massFractions: composition.massFractions,
    summaryResults: {
      "M<sub>1</sub>": formatCompact(mach1),
      "M<sub>2</sub>": formatCompact(mach2),
      "ρ<sub>1</sub>, kg/m<sup>3</sup>": formatCompact(density1),
      "a<sub>1</sub>, m/s": formatCompact(upstream.soundSpeed),
      "V<sub>1</sub>, m/s": formatCompact(velocityNormal),
      "T<sub>2</sub>, K": formatCompact(temperature2),
      "p<sub>2</sub>, Pa": formatCompact(pressure2Pa),
      "ρ<sub>2</sub>, kg/m<sup>3</sup>": formatCompact(density2),
      "V<sub>2</sub>, m/s": formatCompact(velocity2),
    },
    speciesSections: [
      {
        label: "Species density downstream, kg/m<sup>3</sup>",
        rows: HANSEN_SPECIES.map((species, index) => [
          HANSEN_SPECIES_LABELS[species],
          formatCompact(speciesDensities[index]),
        ]),
      },
      {
        label: "Mole fraction (X) downstream",
        rows: HANSEN_SPECIES.map((species, index) => [
          HANSEN_SPECIES_LABELS[species],
          formatCompact(composition.moleFractions[index]),
        ]),
      },
      {
        label: "Mass fraction (Y) downstream",
        rows: HANSEN_SPECIES.map((species, index) => [
          HANSEN_SPECIES_LABELS[species],
          formatCompact(composition.massFractions[index]),
        ]),
      },
    ],
  };
}

function reactingObliqueShock(mach1, temperature1, pressure1Pa, thetaDegrees, xi = HANSEN_DEFAULT_XI) {
  if (mach1 <= 1 || !Number.isFinite(mach1)) throw new Error("M1 must be greater than 1.");
  if (thetaDegrees <= 0 || !Number.isFinite(thetaDegrees)) throw new Error("Turn angle must be greater than 0 degrees.");
  if (temperature1 <= 0 || !Number.isFinite(temperature1)) throw new Error("T1 must be greater than 0 K.");
  if (pressure1Pa <= 0 || !Number.isFinite(pressure1Pa)) throw new Error("p1 must be greater than 0 Pa.");

  const theta = thetaDegrees * DEG;
  const upstream = hansenEquilibriumState(temperature1, pressure1Pa, xi);
  const velocity1 = mach1 * upstream.soundSpeed;
  const betaMin = Math.asin(1 / mach1) + 1e-6;
  const betaMax = Math.PI / 2 - 1e-6;

  const residual = (beta) => {
    try {
      const normalMach1 = mach1 * Math.sin(beta);
      if (normalMach1 <= 1) return Number.NaN;
      const shock = normalShockFull(normalMach1, temperature1, pressure1Pa, xi);
      const normalVelocityRatio = shock.density1 / shock.density2;
      return Math.tan(beta - theta) - normalVelocityRatio * Math.tan(beta);
    } catch (exception) {
      return Number.NaN;
    }
  };

  let previousBeta = null;
  let previousResidual = null;
  let bracket = null;
  for (let i = 0; i <= 240; i += 1) {
    const beta = betaMin + (i / 240) * (betaMax - betaMin);
    const value = residual(beta);
    if (Number.isFinite(previousResidual) && Number.isFinite(value) && previousResidual * value <= 0) {
      bracket = [previousBeta, beta];
      break;
    }
    if (Number.isFinite(value)) {
      previousBeta = beta;
      previousResidual = value;
    }
  }

  if (!bracket) {
    throw new Error("No attached reacting oblique shock found for this M1 and turn angle.");
  }

  const beta = solveBisection(residual, bracket[0], bracket[1], 1e-10, 100);
  const normalMach1 = mach1 * Math.sin(beta);
  const normalShock = normalShockFull(normalMach1, temperature1, pressure1Pa, xi);
  const tangentialVelocity = velocity1 * Math.cos(beta);
  const normalVelocity2 = normalShock.velocity2;
  const velocity2 = Math.hypot(tangentialVelocity, normalVelocity2);
  const downstream = hansenEquilibriumState(normalShock.temperature2, normalShock.pressure2Pa, xi);
  const mach2 = velocity2 / downstream.soundSpeed;
  const speciesDensities = normalShock.massFractions.map((massFraction) => massFraction * normalShock.density2);

  return {
    summaryResults: {
      "M<sub>1</sub>": formatCompact(mach1),
      "θ, deg.": formatCompact(thetaDegrees),
      "β, deg.": formatCompact(beta * RAD),
      "M<sub>n,1</sub>": formatCompact(normalMach1),
      "M<sub>2</sub>": formatCompact(mach2),
      "T<sub>2</sub>, K": formatCompact(normalShock.temperature2),
      "p<sub>2</sub>, Pa": formatCompact(normalShock.pressure2Pa),
      "ρ<sub>2</sub>, kg/m<sup>3</sup>": formatCompact(normalShock.density2),
      "V<sub>2</sub>, m/s": formatCompact(velocity2),
    },
    speciesSections: [
      {
        label: "Species density downstream, kg/m<sup>3</sup>",
        rows: HANSEN_SPECIES.map((species, index) => [
          HANSEN_SPECIES_LABELS[species],
          formatCompact(speciesDensities[index]),
        ]),
      },
      {
        label: "Mole fraction (X) downstream",
        rows: HANSEN_SPECIES.map((species, index) => [
          HANSEN_SPECIES_LABELS[species],
          formatCompact(normalShock.moleFractions[index]),
        ]),
      },
      {
        label: "Mass fraction (Y) downstream",
        rows: HANSEN_SPECIES.map((species, index) => [
          HANSEN_SPECIES_LABELS[species],
          formatCompact(normalShock.massFractions[index]),
        ]),
      },
    ],
  };
}

function reactingConeShockSurfaceAtBeta(gamma, mach1, beta, coneAngle, normalShock) {
  const velocity1 = nondimensionalVelocity(gamma, mach1);
  const densityRatio = normalShock.density2 / normalShock.density1;
  const radialVelocity = velocity1 * Math.cos(beta);
  const polarVelocity = -velocity1 * Math.sin(beta) / densityRatio;
  const surface = integrateTaylorMaccoll(gamma, beta, coneAngle, radialVelocity, polarVelocity);

  return {
    radialVelocity,
    polarVelocity,
    surface,
  };
}

function reactingConeShock(mach1, temperature1, pressure1Pa, coneAngleDegrees, xi = HANSEN_DEFAULT_XI) {
  if (mach1 <= 1 || !Number.isFinite(mach1)) throw new Error("M1 must be greater than 1.");
  if (coneAngleDegrees <= 0 || !Number.isFinite(coneAngleDegrees)) {
    throw new Error("Cone half-angle must be greater than 0 degrees.");
  }
  if (coneAngleDegrees >= 89) throw new Error("Cone half-angle must be less than 89 degrees.");
  if (temperature1 <= 0 || !Number.isFinite(temperature1)) throw new Error("T1 must be greater than 0 K.");
  if (pressure1Pa <= 0 || !Number.isFinite(pressure1Pa)) throw new Error("p1 must be greater than 0 Pa.");

  const coneAngle = coneAngleDegrees * DEG;
  const upstream = hansenEquilibriumState(temperature1, pressure1Pa, xi);
  const gamma = Math.min(Math.max(upstream.phiGamma, 1.05), 1.67);
  const velocity1 = mach1 * upstream.soundSpeed;
  const betaMin = Math.max(Math.asin(1 / mach1), coneAngle) + 1e-5;
  const betaMax = Math.PI / 2 - 1e-5;
  const shockCache = new Map();

  const shockForBeta = (beta) => {
    const key = beta.toPrecision(14);
    if (!shockCache.has(key)) {
      shockCache.set(key, normalShockFull(mach1 * Math.sin(beta), temperature1, pressure1Pa, xi));
    }
    return shockCache.get(key);
  };

  const residual = (beta) => {
    try {
      const normalMach1 = mach1 * Math.sin(beta);
      if (normalMach1 <= 1) return Number.NaN;
      const normalShock = shockForBeta(beta);
      return reactingConeShockSurfaceAtBeta(gamma, mach1, beta, coneAngle, normalShock).surface.polarVelocity;
    } catch (exception) {
      return Number.NaN;
    }
  };

  let previousBeta = null;
  let previousResidual = null;
  let bracket = null;
  for (let i = 0; i <= 140; i += 1) {
    const beta = betaMin + (i / 140) * (betaMax - betaMin);
    const value = residual(beta);
    if (Number.isFinite(previousResidual) && Number.isFinite(value) && previousResidual * value <= 0) {
      bracket = [previousBeta, beta];
      break;
    }
    if (Number.isFinite(value)) {
      previousBeta = beta;
      previousResidual = value;
    }
  }

  if (!bracket) {
    throw new Error("No attached reacting cone shock found for this M1 and cone angle.");
  }

  const beta = solveBisection(residual, bracket[0], bracket[1], 1e-9, 80);
  if (beta <= coneAngle) {
    throw new Error("Invalid reacting cone shock solution: shock angle must be greater than the cone half-angle.");
  }
  const normalMach1 = mach1 * Math.sin(beta);
  const normalShock = shockForBeta(beta);
  const coneVelocityState = reactingConeShockSurfaceAtBeta(gamma, mach1, beta, coneAngle, normalShock);
  const tangentialVelocity = velocity1 * Math.cos(beta);
  const velocityBehindShock = Math.hypot(tangentialVelocity, normalShock.velocity2);
  const downstream = hansenEquilibriumState(normalShock.temperature2, normalShock.pressure2Pa, xi);
  const machBehindShock = velocityBehindShock / downstream.soundSpeed;
  const shockVelocity = Math.hypot(coneVelocityState.radialVelocity, coneVelocityState.polarVelocity);
  const surfaceVelocity = Math.hypot(
    coneVelocityState.surface.radialVelocity,
    coneVelocityState.surface.polarVelocity,
  );
  const surfaceMach = machFromNondimensionalVelocity(gamma, surfaceVelocity);
  const surfaceToShockTemperatureRatio = (1 - surfaceVelocity * surfaceVelocity)
    / (1 - shockVelocity * shockVelocity);
  const surfaceToShockPressureRatio = Math.pow(surfaceToShockTemperatureRatio, gamma / (gamma - 1));
  const surfaceTemperature = normalShock.temperature2 * surfaceToShockTemperatureRatio;
  const surfacePressurePa = normalShock.pressure2Pa * surfaceToShockPressureRatio;
  const surfaceState = hansenEquilibriumState(surfaceTemperature, surfacePressurePa, xi);
  const shockSpeciesDensities = normalShock.massFractions.map((massFraction) => (
    massFraction * normalShock.density2
  ));
  const surfaceSpeciesDensities = surfaceState.massFractions.map((massFraction) => (
    massFraction * surfaceState.density
  ));

  return {
    summaryResults: {
      __sectionBefore: {
        "M<sub>n,1</sub>": "Behind shock",
        "M<sub>c</sub>": "Cone surface",
      },
      "M<sub>1</sub>": formatCompact(mach1),
      "θ<sub>c</sub> cone half-angle, deg.": formatCompact(coneAngleDegrees),
      "β shock angle from axis, deg.": formatCompact(beta * RAD),
      "M<sub>n,1</sub>": formatCompact(normalMach1),
      "M<sub>2</sub>": formatCompact(machBehindShock),
      "T<sub>2</sub>, K": formatCompact(normalShock.temperature2),
      "p<sub>2</sub>, Pa": formatCompact(normalShock.pressure2Pa),
      "ρ<sub>2</sub>, kg/m<sup>3</sup>": formatCompact(normalShock.density2),
      "V<sub>2</sub>, m/s": formatCompact(velocityBehindShock),
      "M<sub>c</sub>": formatCompact(surfaceMach),
      "T<sub>c</sub>, K": formatCompact(surfaceTemperature),
      "p<sub>c</sub>, Pa": formatCompact(surfacePressurePa),
      "ρ<sub>c</sub>, kg/m<sup>3</sup>": formatCompact(surfaceState.density),
    },
    speciesSections: [
      {
        label: "Species density behind shock, kg/m<sup>3</sup>",
        rows: HANSEN_SPECIES.map((species, index) => [
          HANSEN_SPECIES_LABELS[species],
          formatCompact(shockSpeciesDensities[index]),
        ]),
      },
      {
        label: "Mole fraction (X) behind shock",
        rows: HANSEN_SPECIES.map((species, index) => [
          HANSEN_SPECIES_LABELS[species],
          formatCompact(normalShock.moleFractions[index]),
        ]),
      },
      {
        label: "Mass fraction (Y) behind shock",
        rows: HANSEN_SPECIES.map((species, index) => [
          HANSEN_SPECIES_LABELS[species],
          formatCompact(normalShock.massFractions[index]),
        ]),
      },
      {
        label: "Species density cone surface, kg/m<sup>3</sup>",
        rows: HANSEN_SPECIES.map((species, index) => [
          HANSEN_SPECIES_LABELS[species],
          formatCompact(surfaceSpeciesDensities[index]),
        ]),
      },
      {
        label: "Mole fraction (X) cone surface",
        rows: HANSEN_SPECIES.map((species, index) => [
          HANSEN_SPECIES_LABELS[species],
          formatCompact(surfaceState.moleFractions[index]),
        ]),
      },
      {
        label: "Mass fraction (Y) cone surface",
        rows: HANSEN_SPECIES.map((species, index) => [
          HANSEN_SPECIES_LABELS[species],
          formatCompact(surfaceState.massFractions[index]),
        ]),
      },
    ],
  };
}

function hansenMeanFreePathRatios(speciesMoleFractions, temperature, pressureAtm, state, mode = "viscosity") {
  const collisionData = hansenCollisionData(temperature);
  const coulombLog = hansenCoulombLog(temperature, pressureAtm, speciesMoleFractions["e-"], state.Z);
  const ratios = {};

  HANSEN_SPECIES.forEach((speciesI, i) => {
    const denominator = HANSEN_SPECIES.reduce((sum, speciesJ, j) => {
      const xj = speciesMoleFractions[speciesJ];
      if (xj <= 0) return sum;
      const sectionRatio = hansenCrossSectionRatio(speciesI, speciesJ, collisionData, coulombLog, mode);
      if (sectionRatio <= 0) return sum;
      return sum + xj * sectionRatio * Math.sqrt((1 + HANSEN_MOLAR_MASS[i] / HANSEN_MOLAR_MASS[j]) / 2);
    }, 0);
    ratios[speciesI] = denominator > 0 ? 1 / denominator : 0;
  });

  return ratios;
}

function hansenViscosityReference(temperature) {
  return 1.462e-6 * Math.sqrt(temperature) / (1 + 112 / temperature);
}

function hansenThermalConductivityReference(temperature, feedMolecularWeight = HANSEN_UNDISSOCIATED_AIR_MW) {
  return (19 / 4) * (UNIVERSAL_GAS_CONSTANT / feedMolecularWeight)
    * hansenViscosityReference(temperature);
}

function hansenViscosityRatio(speciesMoleFractions, meanFreePathRatios, feedMolecularWeight = HANSEN_UNDISSOCIATED_AIR_MW) {
  return HANSEN_SPECIES.reduce((sum, species, index) => {
    const x = speciesMoleFractions[species];
    if (x <= 0) return sum;
    return sum + x * Math.sqrt(HANSEN_MOLAR_MASS[index] / feedMolecularWeight)
      * meanFreePathRatios[species];
  }, 0);
}

function hansenMolecularConductivityRatio(speciesMoleFractions, meanFreePathRatios, temperature, feedMolecularWeight = HANSEN_UNDISSOCIATED_AIR_MW) {
  return HANSEN_SPECIES.reduce((sum, species, index) => {
    const x = speciesMoleFractions[species];
    if (x <= 0) return sum;
    const cvOverR = speciesCvOverR(species, temperature);
    const euckenFactor = (4 * cvOverR + 9) / 19;
    return sum + x * Math.sqrt(feedMolecularWeight / HANSEN_MOLAR_MASS[index])
      * meanFreePathRatios[species] * euckenFactor;
  }, 0);
}

function dlnKpDT(reaction, temperature) {
  const step = 1e-3 * temperature;
  const lnK = (T) => {
    const constants = hansenEquilibriumConstants(T);
    if (reaction === "O2") return Math.log(constants.kpO2);
    if (reaction === "N2") return Math.log(constants.kpN2);
    if (reaction === "OIon") {
      return -158000 / T + lnQp("O+", T) + lnQp("e-", T) - lnQp("O", T);
    }
    if (reaction === "NIon") {
      return -168800 / T + lnQp("N+", T) + lnQp("e-", T) - lnQp("N", T);
    }
    return Math.log(constants.kpIon);
  };
  return (lnK(temperature + step) - lnK(temperature - step)) / (2 * step);
}

function hansenReactiveConductivityRatio(
  speciesMoleFractions,
  temperature,
  reaction,
  coulombLog = 10,
  feedMolecularWeight = HANSEN_UNDISSOCIATED_AIR_MW,
) {
  const collisionData = hansenCollisionData(temperature);
  const stoichiometry = {
    O2: { O2: -1, O: 2 },
    N2: { N2: -1, N: 2 },
    OIon: { O: -1, "O+": 1, "e-": 1 },
    NIon: { N: -1, "N+": 1, "e-": 1 },
  }[reaction];
  const species = Object.keys(stoichiometry);
  const denominator = species.reduce((outerSum, speciesI) => {
    const i = HANSEN_SPECIES.indexOf(speciesI);
    const xi = Math.max(speciesMoleFractions[speciesI], 1e-30);
    const ai = stoichiometry[speciesI];
    return outerSum + species.reduce((innerSum, speciesJ) => {
      const j = HANSEN_SPECIES.indexOf(speciesJ);
      const aj = stoichiometry[speciesJ];
      const xj = Math.max(speciesMoleFractions[speciesJ], 1e-30);
      const sectionRatio = hansenCrossSectionRatio(speciesI, speciesJ, collisionData, coulombLog, "reaction");
      return innerSum + Math.sqrt(
        (HANSEN_MOLAR_MASS[i] * HANSEN_MOLAR_MASS[j])
        / (feedMolecularWeight * (HANSEN_MOLAR_MASS[i] + HANSEN_MOLAR_MASS[j]))
      ) * sectionRatio * ai * (ai * xj - aj * xi) / xi;
    }, 0);
  }, 0);

  if (Math.abs(denominator) < 1e-30) return 0;
  const slope = temperature * dlnKpDT(reaction, temperature);
  return (12 * Math.sqrt(2) / 95) * slope * slope / Math.abs(denominator);
}

function hansenTransportProperties(
  temperature,
  pressureAtm,
  state,
  speciesMoleFractions,
  cpOverR,
  xi = HANSEN_DEFAULT_XI,
) {
  const feedMolecularWeight = hansenFeedMolecularWeight(xi);
  const meanFreePathRatios = hansenMeanFreePathRatios(
    speciesMoleFractions,
    temperature,
    pressureAtm,
    state,
  );
  const viscosityRatio = hansenViscosityRatio(speciesMoleFractions, meanFreePathRatios, feedMolecularWeight);
  const molecularConductivityRatio = hansenMolecularConductivityRatio(
    speciesMoleFractions,
    meanFreePathRatios,
    temperature,
    feedMolecularWeight,
  );
  const coulombLog = hansenCoulombLog(temperature, pressureAtm, speciesMoleFractions["e-"], state.Z);
  const reactiveConductivityRatio = temperature <= 1200
    ? 0
    : hansenReactiveConductivityRatio(speciesMoleFractions, temperature, "O2", coulombLog, feedMolecularWeight)
      + hansenReactiveConductivityRatio(speciesMoleFractions, temperature, "N2", coulombLog, feedMolecularWeight)
      + hansenReactiveConductivityRatio(speciesMoleFractions, temperature, "OIon", coulombLog, feedMolecularWeight)
      + hansenReactiveConductivityRatio(speciesMoleFractions, temperature, "NIon", coulombLog, feedMolecularWeight);
  const conductivityRatio = molecularConductivityRatio + reactiveConductivityRatio;
  const viscosity = hansenViscosityReference(temperature) * viscosityRatio;
  const thermalConductivity = hansenThermalConductivityReference(temperature, feedMolecularWeight) * conductivityRatio;
  const prandtl = conductivityRatio > 0
    ? (4 / 19) * cpOverR * viscosityRatio / conductivityRatio
    : null;

  return {
    viscosityRatio,
    conductivityRatio,
    viscosity,
    thermalConductivity,
    prandtl,
  };
}

function hansenSpeciesConcentrationsFromInput(temperature, pressurePa, xi = HANSEN_DEFAULT_XI) {
  if (pressurePa <= 0 || !Number.isFinite(pressurePa)) throw new Error("Pressure must be greater than 0 Pa.");
  const pressureAtm = pressurePa / ATM_TO_PA;
  const state = hansenAir(temperature, pressureAtm, xi);
  const { xN2Feed, xO2Feed } = hansenFeedFractions(xi);
  const totalConcentration = (pressureAtm * ATM_TO_PA) / (state.Z * HANSEN_R_UNIVERSAL * temperature);
  const gasConstant = hansenMixtureGasConstant(xi);
  const cpOverR = temperature <= 1200 ? 3.5 : hansenCpOverR(temperature, pressureAtm, xi);
  const phiGamma = temperature <= 1200
    ? 1.4
    : cpOverR / hansenCvOverR(temperature, pressureAtm, xi)
      * hansenPhiFactor(temperature, pressureAtm, xi);
  const soundSpeed = Math.sqrt(phiGamma * gasConstant * temperature * state.Z);
  const speciesMoleFractions = {
    N2: state.xN2,
    O2: state.xO2,
    N: state.xN,
    O: state.xO,
    "O+": state["xO+"],
    "N+": state["xN+"],
    "e-": state["xe-"],
  };
  const mixtureMolarMass = HANSEN_SPECIES.reduce((sum, species, index) => (
    sum + speciesMoleFractions[species] * HANSEN_MOLAR_MASS[index]
  ), 0);
  const speciesRows = HANSEN_SPECIES.map((name, index) => {
    const moleFraction = speciesMoleFractions[name];
    return {
      name,
      label: HANSEN_SPECIES_LABELS[name],
      moleFraction,
      massFraction: mixtureMolarMass > 0 ? (moleFraction * HANSEN_MOLAR_MASS[index]) / mixtureMolarMass : 0,
      concentration: moleFraction * totalConcentration,
    };
  });
  const transport = hansenTransportProperties(
    temperature,
    pressureAtm,
    state,
    speciesMoleFractions,
    cpOverR,
    xi,
  );

  const summaryResults = {
    "T, K": formatCompact(temperature),
    "p, Pa": formatCompact(pressurePa),
    "N<sub>2</sub> (X)": formatCompact(xN2Feed),
    "O<sub>2</sub>(X)": formatCompact(xO2Feed),
    Z: formatCompact(state.Z),
    "C<sub>total</sub>, mol/m<sup>3</sup>": formatCompact(totalConcentration),
    "MW<sub>mix</sub>, g/mol": formatCompact(mixtureMolarMass),
    "a<sub>e</sub>, m/s": formatCompact(soundSpeed),
    "γ = a<sub>e</sub><sup>2</sup>ρ/p": formatCompact(phiGamma),
    "η, Pa s": formatCompact(transport.viscosity),
    "k, W/(m K)": formatCompact(transport.thermalConductivity),
    Pr: formatCompact(transport.prandtl),
  };
  const speciesSections = [
    {
      label: "Mole fraction (X)",
      rows: speciesRows.map((item) => [item.label, formatCompact(item.moleFraction)]),
    },
    {
      label: "Mass fraction (Y)",
      rows: speciesRows.map((item) => [item.label, formatCompact(item.massFraction)]),
    },
    {
      label: "Concentrations",
      rows: speciesRows.map((item) => [item.label, formatCompact(item.concentration)]),
    },
  ];

  return {
    summaryResults,
    speciesSections,
  };
}

function format(value) {
  if (value === null || value === undefined || Number.isNaN(value)) return "N/A";
  if (typeof value === "string") return value;
  if (Math.abs(value) >= 1e5 || (Math.abs(value) > 0 && Math.abs(value) < 1e-4)) {
    return value.toExponential(6)
      .replace("e", "E")
      .replace(/E([+-])0*(\d+)/, "E$1$2");
  }
  return value.toPrecision(7).replace(/\.?0+$/, "");
}

function formatCompact(value) {
  if (value === null || value === undefined || Number.isNaN(value)) return "N/A";
  if (typeof value === "string") return value;
  if (value === 0) return "0";

  const absValue = Math.abs(value);
  if (absValue >= 1000 || absValue < 1e-3) {
    return value.toExponential(3)
      .replace("e", "E")
      .replace(/E([+-])0*(\d+)/, "E$1$2");
  }

  return value.toPrecision(5).replace(/\.?0+$/, "");
}

function renderResults(target, results) {
  const rowBreakBefore = results.__rowBreakBefore ?? [];
  const sectionBefore = results.__sectionBefore ?? {};
  target.innerHTML = Object.entries(results)
    .filter(([label]) => !["t2t1", "p02p01", "p1p02", "__rowBreakBefore", "__sectionBefore"].includes(label))
    .map(([label, value]) => `
      ${sectionBefore[label] ? `<div class="result-section-label">${sectionBefore[label]}</div>` : ""}
      ${label === "β strong, deg." ? `<div class="result-branch-label">Strong Oblique Shock</div>` : ""}
      <div class="result-item${["β weak, deg.", "β strong, deg."].includes(label) || rowBreakBefore.includes(label) ? " starts-new-row" : ""}">
        <span class="result-label">${label}</span>
        <span class="result-value">${format(value)}</span>
      </div>
      ${label === "θ<sub>max</sub>, deg." ? `<div class="result-branch-label">Weak Oblique Shock</div>` : ""}
    `)
    .join("");
}

function renderResultItems(results) {
  const sectionBefore = results.__sectionBefore ?? {};
  return Object.entries(results)
    .filter(([label]) => label !== "__sectionBefore")
    .map(([label, value]) => `
      ${sectionBefore[label] ? `<div class="result-section-label">${sectionBefore[label]}</div>` : ""}
      <div class="result-item">
        <span class="result-label">${label}</span>
        <span class="result-value">${format(value)}</span>
      </div>
    `)
    .join("");
}

function renderSpeciesSections(sections) {
  return sections
    .map((section) => `
      <div class="result-section-label">${section.label}</div>
      ${section.rows.map(([label, value]) => `
        <div class="result-item">
          <span class="result-label">${label}</span>
          <span class="result-value">${format(value)}</span>
        </div>
      `).join("")}
    `)
    .join("");
}

function renderReactingAirResults(target, results) {
  target.className = "reacting-air-results";
  target.innerHTML = `
    <div class="result-section-label gas-properties-label">Gas Properties</div>
    <div class="result-grid reacting-summary-grid">
      ${renderResultItems(results.summaryResults)}
    </div>
    <div class="result-grid species-grid">
      ${renderSpeciesSections(results.speciesSections)}
    </div>
  `;
}

function bindCalculator(formId, errorId, resultsId, compute) {
  const form = document.getElementById(formId);
  const error = document.getElementById(errorId);
  const results = document.getElementById(resultsId);

  if (!form || !error || !results) return;

  const calculate = () => {
    error.textContent = "";

    try {
      const data = new FormData(form);
      const gamma = Number(data.get("gamma"));
      const input = data.get("input");
      const value = Number(data.get("value"));
      renderResults(results, compute(gamma, input, value));
    } catch (exception) {
      error.textContent = exception.message;
      results.innerHTML = "";
    }
  };

  form.addEventListener("submit", (event) => {
    event.preventDefault();
    calculate();
  });

  calculate();
}

document.querySelectorAll(".calc-tab").forEach((tab) => {
  tab.addEventListener("click", () => {
    document.querySelectorAll(".calc-tab").forEach((item) => item.classList.remove("is-active"));
    document.querySelectorAll(".calc-panel").forEach((panel) => panel.classList.remove("is-active"));
    tab.classList.add("is-active");
    document.getElementById(tab.dataset.panel).classList.add("is-active");
  });
});

bindCalculator("isentropic-form", "isentropic-error", "isentropic-results", isentropicFromInput);
bindCalculator("normal-shock-form", "normal-shock-error", "normal-shock-results", normalShockFromInput);

const normalShockDimensionalForm = document.getElementById("normal-shock-dimensional-form");
if (normalShockDimensionalForm) {
  const error = document.getElementById("normal-shock-dimensional-error");
  const results = document.getElementById("normal-shock-dimensional-results");
  const calculateNormalShockDimensional = () => {
    error.textContent = "";

    try {
      const data = new FormData(normalShockDimensionalForm);
      const mach = Number(data.get("mach"));
      const gamma = Number(data.get("gamma"));
      const pressure = Number(data.get("pressure"));
      const temperature = Number(data.get("temperature"));
      renderResults(results, normalShockDimensional(gamma, mach, pressure, temperature));
    } catch (exception) {
      error.textContent = exception.message;
      results.innerHTML = "";
    }
  };

  normalShockDimensionalForm.addEventListener("submit", (event) => {
    event.preventDefault();
    calculateNormalShockDimensional();
  });

  calculateNormalShockDimensional();
}

const reactingNormalShockForm = document.getElementById("reacting-normal-shock-form");
if (reactingNormalShockForm) {
  const error = document.getElementById("reacting-normal-shock-error");
  const results = document.getElementById("reacting-normal-shock-results");
  const calculateReactingNormalShock = () => {
    error.textContent = "";

    try {
      const data = new FormData(reactingNormalShockForm);
      const mach = Number(data.get("mach"));
      const temperature = Number(data.get("temperature"));
      const pressure = Number(data.get("pressure"));
      const xi = hansenXiFromSpeciesInput(data);
      renderReactingAirResults(results, normalShockFull(mach, temperature, pressure, xi));
    } catch (exception) {
      error.textContent = exception.message;
      results.innerHTML = "";
    }
  };

  reactingNormalShockForm.addEventListener("submit", (event) => {
    event.preventDefault();
    calculateReactingNormalShock();
  });

  calculateReactingNormalShock();
}

const reactingObliqueShockForm = document.getElementById("reacting-oblique-shock-form");
if (reactingObliqueShockForm) {
  const error = document.getElementById("reacting-oblique-shock-error");
  const results = document.getElementById("reacting-oblique-shock-results");
  const calculateReactingObliqueShock = () => {
    error.textContent = "";

    try {
      const data = new FormData(reactingObliqueShockForm);
      const mach = Number(data.get("mach"));
      const temperature = Number(data.get("temperature"));
      const pressure = Number(data.get("pressure"));
      const theta = Number(data.get("theta"));
      const xi = hansenXiFromSpeciesInput(data);
      renderReactingAirResults(results, reactingObliqueShock(mach, temperature, pressure, theta, xi));
    } catch (exception) {
      error.textContent = exception.message;
      results.innerHTML = "";
    }
  };

  reactingObliqueShockForm.addEventListener("submit", (event) => {
    event.preventDefault();
    calculateReactingObliqueShock();
  });

  calculateReactingObliqueShock();
}

const reactingConeShockForm = document.getElementById("reacting-cone-shock-form");
if (reactingConeShockForm) {
  const error = document.getElementById("reacting-cone-shock-error");
  const results = document.getElementById("reacting-cone-shock-results");
  const calculateReactingConeShock = () => {
    error.textContent = "";

    try {
      const data = new FormData(reactingConeShockForm);
      const mach = Number(data.get("mach"));
      const temperature = Number(data.get("temperature"));
      const pressure = Number(data.get("pressure"));
      const theta = Number(data.get("theta"));
      const xi = hansenXiFromSpeciesInput(data);
      renderReactingAirResults(results, reactingConeShock(mach, temperature, pressure, theta, xi));
    } catch (exception) {
      error.textContent = exception.message;
      results.innerHTML = "";
    }
  };

  reactingConeShockForm.addEventListener("submit", (event) => {
    event.preventDefault();
    calculateReactingConeShock();
  });

  calculateReactingConeShock();
}

const obliqueShockForm = document.getElementById("oblique-shock-form");
if (obliqueShockForm) {
  const error = document.getElementById("oblique-shock-error");
  const results = document.getElementById("oblique-shock-results");
  const calculateObliqueShock = () => {
    error.textContent = "";

    try {
      const data = new FormData(obliqueShockForm);
      const gamma = Number(data.get("gamma"));
      const mach = Number(data.get("mach"));
      const theta = Number(data.get("theta"));
      renderResults(results, obliqueShockFromInput(gamma, mach, theta));
    } catch (exception) {
      error.textContent = exception.message;
      results.innerHTML = "";
    }
  };

  obliqueShockForm.addEventListener("submit", (event) => {
    event.preventDefault();
    calculateObliqueShock();
  });

  calculateObliqueShock();
}

const obliqueShockDimensionalForm = document.getElementById("oblique-shock-dimensional-form");
if (obliqueShockDimensionalForm) {
  const error = document.getElementById("oblique-shock-dimensional-error");
  const results = document.getElementById("oblique-shock-dimensional-results");
  const calculateObliqueShockDimensional = () => {
    error.textContent = "";

    try {
      const data = new FormData(obliqueShockDimensionalForm);
      const gamma = Number(data.get("gamma"));
      const mach = Number(data.get("mach"));
      const pressure = Number(data.get("pressure"));
      const temperature = Number(data.get("temperature"));
      const theta = Number(data.get("theta"));
      renderResults(results, obliqueShockDimensional(gamma, mach, pressure, temperature, theta));
    } catch (exception) {
      error.textContent = exception.message;
      results.innerHTML = "";
    }
  };

  obliqueShockDimensionalForm.addEventListener("submit", (event) => {
    event.preventDefault();
    calculateObliqueShockDimensional();
  });

  calculateObliqueShockDimensional();
}

const coneShockForm = document.getElementById("cone-shock-form");
if (coneShockForm) {
  const error = document.getElementById("cone-shock-error");
  const results = document.getElementById("cone-shock-results");
  const calculateConeShock = () => {
    error.textContent = "";

    try {
      const data = new FormData(coneShockForm);
      const gamma = Number(data.get("gamma"));
      const mach = Number(data.get("mach"));
      const theta = Number(data.get("theta"));
      renderResults(results, coneShockFromInput(gamma, mach, theta));
    } catch (exception) {
      error.textContent = exception.message;
      results.innerHTML = "";
    }
  };

  coneShockForm.addEventListener("submit", (event) => {
    event.preventDefault();
    calculateConeShock();
  });

  calculateConeShock();
}

const speedForm = document.getElementById("speed-sound-form");
if (speedForm) {
  const gasSelect = speedForm.elements.gas;
  const gammaInput = speedForm.elements.gamma;
  const gasConstantInput = speedForm.elements.gasConstant;
  const molecularWeightInput = speedForm.elements.molecularWeight;
  const molecularWeightField = speedForm.querySelector(".mw-field");

  gasSelect.addEventListener("change", () => {
    const gas = GAS_PROPERTIES[gasSelect.value];
    const isCustom = gasSelect.value === "custom";
    gammaInput.readOnly = !isCustom;
    gasConstantInput.readOnly = true;
    molecularWeightField.classList.toggle("is-hidden", !isCustom);

    if (gas) {
      gammaInput.value = gas.gamma;
      molecularWeightInput.value = gas.molecularWeight;
      gasConstantInput.value = format(gasConstantFromMolecularWeight(gas.molecularWeight));
    } else {
      gasConstantInput.value = format(gasConstantFromMolecularWeight(Number(molecularWeightInput.value)));
    }
  });

  molecularWeightInput.addEventListener("input", () => {
    if (gasSelect.value !== "custom") return;
    try {
      gasConstantInput.value = format(gasConstantFromMolecularWeight(Number(molecularWeightInput.value)));
    } catch {
      gasConstantInput.value = "";
    }
  });

  gasSelect.dispatchEvent(new Event("change"));

  const error = document.getElementById("speed-sound-error");
  const results = document.getElementById("speed-sound-results");
  const calculateSpeed = () => {
    error.textContent = "";

    try {
      const data = new FormData(speedForm);
      const gas = data.get("gas");
      const gamma = Number(data.get("gamma"));
      const molecularWeight = Number(data.get("molecularWeight"));
      const gasConstant = gas === "custom"
        ? gasConstantFromMolecularWeight(molecularWeight)
        : Number(data.get("gasConstant"));
      const temperature = Number(data.get("temperature"));
      gasConstantInput.value = format(gasConstant);
      renderResults(results, speedOfSoundFromInput(gamma, gasConstant, temperature, gas, gas === "custom" ? molecularWeight : null));
    } catch (exception) {
      error.textContent = exception.message;
      results.innerHTML = "";
    }
  };

  speedForm.addEventListener("submit", (event) => {
    event.preventDefault();
    calculateSpeed();
  });

  gasSelect.addEventListener("change", calculateSpeed);
  calculateSpeed();
}

const reactingAirForm = document.getElementById("reacting-air-form");
if (reactingAirForm) {
  const error = document.getElementById("reacting-air-error");
  const results = document.getElementById("reacting-air-results");
  const calculateReactingAir = () => {
    error.textContent = "";

    try {
      const data = new FormData(reactingAirForm);
      const temperature = Number(data.get("temperature"));
      const pressure = Number(data.get("pressure"));
      const xi = hansenXiFromSpeciesInput(data);
      renderReactingAirResults(results, hansenSpeciesConcentrationsFromInput(temperature, pressure, xi));
    } catch (exception) {
      error.textContent = exception.message;
      results.innerHTML = "";
    }
  };

  reactingAirForm.addEventListener("submit", (event) => {
    event.preventDefault();
    calculateReactingAir();
  });

  calculateReactingAir();
}
