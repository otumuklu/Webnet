/* Browser-side finite-rate reacting-air model.
 * Park5 is the first transient model; Hansen remains the algebraic equilibrium model.
 * The same file runs in a Web Worker and exposes a small message API.
 */
(function () {
  if (typeof importScripts === "function") importScripts("reacting_air_data.js");
  const SPECIES = ["N2", "O2", "NO", "N", "O"];
  const W = { N2: 28.0134, O2: 31.9988, NO: 30.0061, N: 14.0067, O: 15.9994 };
  const R_UNIVERSAL = 8314.462618;
  const K_B = 1.380649e-23;
  const N_A = 6.02214076e26;
  const P_ATM = 101325;
  const VIB = ["N2", "O2", "NO"];
  const THETA = { N2: 3371, O2: 2256, NO: 2719 };
  const W_ION = { "N2+": 28.0129, "O2+": 31.9983, "NO+": 30.0056 };
  Object.assign(W, W_ION);
  Object.assign(THETA, { "N2+": 3371, "O2+": 2256, "NO+": 2719 });
  const HF = { N2: 0, O2: 0, NO: 10981, N: 56852, O: 29968 };
  const DISS = { N2: 3.36e7, O2: 1.54e7, NO: 2.09e7 };
  const CVTR = {
    N2: 2.5 * R_UNIVERSAL / W.N2, O2: 2.5 * R_UNIVERSAL / W.O2,
    NO: 2.5 * R_UNIVERSAL / W.NO, N: 1.5 * R_UNIVERSAL / W.N,
    O: 1.5 * R_UNIVERSAL / W.O,
  };

  const REACTIONS = [
    { lhs: { O2: 1 }, rhs: { O: 2 }, A: 2e18, beta: -1.5, Ta: 59500,
      eta: { N2: 1, O2: 1, NO: 1, N: 5, O: 5 },
      ni: [1e20, 1e21, 1e22, 1e23, 1e24, 1e25],
      A0: [1.8103, .91354, .64183, .55388, .52455, .50989],
      A1: [8.8685, 9.2238, 9.3331, 9.3678, 9.3793, 9.3851],
      A2: [3.5716, 2.2885, 1.9026, 1.7763, 1.7342, 1.7132],
      A3: [-7.3623, -6.7969, -6.6277, -6.572, -6.5534, -6.5441],
      A4: [.083861, .046338, .035151, .031445, .030209, .029591] },
    { lhs: { N2: 1 }, rhs: { N: 2 }, A: 7e18, beta: -1.6, Ta: 113220,
      eta: { N2: 1, O2: 1, NO: 1, N: 4.286, O: 4.286 },
      ni: [1e20, 1e21, 1e22, 1e23, 1e24, 1e25],
      A0: [3.4907, 2.0723, 1.606, 1.5351, 1.4766, 1.4766],
      A1: [7.73913, 8.2975, 8.481, 8.5139, 8.5369, 8.5369],
      A2: [4.0978, 2.0617, 1.3923, 1.2993, 1.2153, 1.2153],
      A3: [-12.728, -11.828, -11.533, -11.494, -11.457, -11.457],
      A4: [.07487, .015105, -.004543, -.00698, -.009444, -.009444] },
    { lhs: { NO: 1 }, rhs: { N: 1, O: 1 }, A: 5e12, beta: 0, Ta: 75500,
      eta: { N2: 1, O2: 1, NO: 1, N: 22, O: 22 },
      ni: [1e20, 1e21, 1e22, 1e23, 1e24, 1e25],
      A0: [2.1649, 1.0072, .63817, .55889, .515, .50765],
      A1: [6.986377, 7.44325, 7.58969, 7.62338, 7.64066, 7.64355],
      A2: [2.8508, 1.1911, .66336, .55396, .49096, .48042],
      A3: [-8.5422, -7.8098, -7.5773, -7.5304, -7.5025, -7.4979],
      A4: [.053043, .004394, -.011025, -.014089, -.015938, -.016247] },
  ];
  const EXCHANGE = [
    { lhs: { NO: 1, O: 1 }, rhs: { O2: 1, N: 1 }, A: 8.4e9, beta: 0, Ta: 19450,
      coeff: [[.35438,.093613,-.003732,.004815,-.009758,-.002428],[-1.8821,-1.7806,-1.7434,-1.7443,-1.7386,-1.7415],[-.72111,-1.0975,-1.2394,-1.2227,-1.2436,-1.2331],[-1.1797,-1.0128,-.94952,-.95824,-.949,-.95365],[-.030831,-.041949,-.046182,-.045545,-.046159,-.04585]] },
    { lhs: { N2: 1, O: 1 }, rhs: { NO: 1, N: 1 }, A: 6.4e14, beta: -1, Ta: 38400,
      coeff: [[1.3261,1.0653,.96794,.97646,.96188,.96921],[.75268,.85417,.89131,.89043,.89617,.89329],[1.2474,.87093,.7291,.74572,.72479,.73531],[-4.1857,-4.0188,-3.9555,-3.9642,-3.955,-3.9596],[.02184,.010721,.006488,.007123,.006509,.006818]] },
  ];

  function clamp(value, low, high) { return Math.max(low, Math.min(high, value)); }
  function arrh(A, beta, Ta, temperature) {
    const T = Math.max(temperature, 1);
    return A * Math.pow(T, beta) * Math.exp(-Ta / T);
  }
  function interpolate(values, density, grid) {
    if (density <= grid[0]) return values[0];
    if (density >= grid[grid.length - 1]) return values[values.length - 1];
    for (let i = 0; i < grid.length - 1; i += 1) {
      if (density <= grid[i + 1]) {
        const f = (density - grid[i]) / (grid[i + 1] - grid[i]);
        return values[i] + f * (values[i + 1] - values[i]);
      }
    }
    return values[values.length - 1];
  }
  function keq(temperature, pressure, reaction) {
    const T = Math.max(temperature, 1);
    const nD = pressure / (K_B * T);
    const z = 1e4 / T;
    const a = reaction.coeff || [reaction.A0, reaction.A1, reaction.A2, reaction.A3, reaction.A4];
    return Math.exp(interpolate(a[0], nD, reaction.ni) / z + interpolate(a[1], nD, reaction.ni)
      + interpolate(a[2], nD, reaction.ni) * Math.log(z) + interpolate(a[3], nD, reaction.ni) * z
      + interpolate(a[4], nD, reaction.ni) * z * z);
  }
  function vibEnergy(Tv, species) {
    const theta = THETA[species];
    const x = theta / Math.max(Tv, 1);
    return (R_UNIVERSAL / W[species]) * theta / Math.expm1(x);
  }
  function vibCv(Tv, species) {
    const theta = THETA[species];
    const x = theta / Math.max(Tv, 1);
    const ex = Math.exp(x);
    return (R_UNIVERSAL / W[species]) * x * x * ex / Math.pow(ex - 1, 2);
  }
  function temperatureFromEnergy(Ev, Y, rho, guess, vibSpecies = VIB) {
    let Tv = clamp(guess, 150, 50000);
    for (let i = 0; i < 40; i += 1) {
      let f = -Ev;
      let df = 0;
      vibSpecies.forEach((species) => { f += rho * Y[species] * vibEnergy(Tv, species); df += rho * Y[species] * vibCv(Tv, species); });
      const next = clamp(Tv - f / Math.max(df, 1e-8), 150, 50000);
      if (Math.abs(next - Tv) < 1e-7 * Math.abs(next) + 1e-8) return next;
      Tv = next;
    }
    return Tv;
  }
  function compositionFromMoles(xN2, xO2, species = SPECIES, supplied = {}) {
    const values = {}; species.forEach((s) => { values[s] = Math.max(Number(supplied[s] ?? 0), 0); });
    values.N2 = Math.max(Number(supplied.N2 ?? xN2), 0);
    values.O2 = Math.max(Number(supplied.O2 ?? xO2), 0);
    const sum = species.reduce((total, speciesName) => total + values[speciesName], 0);
    if (sum <= 0) throw new Error("N2 and O2 composition must contain positive total mole fraction.");
    const X = {}; species.forEach((s) => { X[s] = values[s] / sum; });
    const Wmix = species.reduce((sumW, s) => sumW + X[s] * W[s], 0);
    const Y = {}; species.forEach((s) => { Y[s] = X[s] * W[s] / Wmix; });
    return { X, Y, Wmix };
  }
  function products(concentrations, stoich) {
    return Object.entries(stoich).reduce((value, [s, power]) => value * Math.pow(Math.max(concentrations[s], 0), power), 1);
  }
  function addReaction(omega, reaction, net) {
    Object.entries(reaction.lhs).forEach(([s, nu]) => { omega[s] -= nu * net; });
    Object.entries(reaction.rhs).forEach(([s, nu]) => { omega[s] += nu * net; });
  }
  function step(state, rho, dt, twoTemperature) {
    const Y = {}; SPECIES.forEach((s, i) => { Y[s] = Math.max(state[i], 0); });
    const ySum = SPECIES.reduce((sum, s) => sum + Y[s], 0);
    SPECIES.forEach((s) => { Y[s] /= Math.max(ySum, 1e-30); });
    const T = clamp(state[5], 300, 50000);
    const Tv = twoTemperature ? temperatureFromEnergy(Math.max(state[6], 0), Y, rho, state[7]) : T;
    const c = {}; SPECIES.forEach((s) => { c[s] = rho * Y[s] / W[s]; });
    const cTotal = SPECIES.reduce((sum, s) => sum + c[s], 0);
    const pressure = cTotal * R_UNIVERSAL * T;
    const nD = {}; SPECIES.forEach((s) => { nD[s] = c[s] * N_A; });
    const Tpark = Math.pow(T, .7) * Math.pow(Math.max(Tv, 1), .3);
    const omega = {}; SPECIES.forEach((s) => { omega[s] = 0; });

    REACTIONS.forEach((reaction) => {
      const kf = arrh(reaction.A, reaction.beta, reaction.Ta, Tpark);
      const kr = arrh(reaction.A, reaction.beta, reaction.Ta, T) / Math.max(keq(T, pressure, reaction), 1e-300);
      const thirdBody = SPECIES.reduce((sum, s) => sum + (reaction.eta[s] || 0) * c[s], 0);
      const net = (kf * c[Object.keys(reaction.lhs)[0]] - kr * products(c, reaction.rhs)) * thirdBody;
      addReaction(omega, reaction, net);
    });
    EXCHANGE.forEach((reaction) => {
      const kf = arrh(reaction.A, reaction.beta, reaction.Ta, T);
      const kr = kf / Math.max(keq(T, pressure, { ni: REACTIONS[0].ni, coeff: reaction.coeff }), 1e-300);
      addReaction(omega, reaction, kf * products(c, reaction.lhs) - kr * products(c, reaction.rhs));
    });

    const RR = {}; SPECIES.forEach((s) => { RR[s] = omega[s] * W[s]; });
    let QVT = 0;
    VIB.forEach((v) => {
      const A = v === "N2" ? 221.53 : v === "O2" ? 135.91 : 49.5;
      const B = v === "N2" ? .029 : v === "O2" ? .03 : .042;
      const tauMW = (P_ATM / Math.max(pressure, 1)) * Math.exp(A * (Math.pow(T, -1 / 3) - B) - 18.42);
      const tauPark = 1 / Math.max(Math.sqrt(8 * R_UNIVERSAL * T / (Math.PI * W[v])) * 3e-21 * Math.pow(5e4 / T, 2) * nD[v], 1e-300);
      const tau = tauMW + tauPark;
      QVT += rho * Y[v] * (vibEnergy(T, v) - vibEnergy(Tv, v)) / tau;
    });
    const Scv = .3 * VIB.reduce((sum, s) => sum + DISS[s] * RR[s], 0);
    const Sh = -SPECIES.reduce((sum, s) => sum + (R_UNIVERSAL * HF[s] / W[s]) * RR[s], 0);
    const cvMix = SPECIES.reduce((sum, s) => sum + Y[s] * CVTR[s], 0);
    const sumRRE = SPECIES.reduce((sum, s) => sum + RR[s] * CVTR[s] * T, 0);
    const dT = twoTemperature ? (Sh - QVT - Scv - sumRRE) / Math.max(rho * cvMix, 1e-30) : (Sh - sumRRE) / Math.max(rho * cvMix, 1e-30);
    const dEv = twoTemperature ? QVT + Scv : 0;
    const next = state.slice(); SPECIES.forEach((s, i) => { next[i] = Math.max(0, Y[s] + dt * RR[s] / rho); });
    next[5] = clamp(T + dt * dT, 300, 50000);
    next[6] = Math.max(0, state[6] + dt * dEv);
    next[7] = Tv;
    return next;
  }
  function runPark5(input) {
    const T0 = Number(input.temperature); const Tv0 = Number(input.vibrationalTemperature);
    const pressure0 = Number(input.pressure); const tEnd = Number(input.endTime);
    if (!(T0 > 0 && Tv0 > 0 && pressure0 > 0 && tEnd > 0)) throw new Error("Temperatures, pressure, and end time must be positive.");
    const initial = compositionFromMoles(Number(input.n2Feed), Number(input.o2Feed), SPECIES, input.composition);
    const rho = pressure0 / ((R_UNIVERSAL / initial.Wmix) * T0);
    const state = SPECIES.map((s) => initial.Y[s]);
    state.push(T0, VIB.reduce((sum, s) => sum + rho * initial.Y[s] * vibEnergy(Tv0, s), 0), Tv0);
    const points = []; const count = 300; const targets = [];
    for (let i = 0; i < count; i += 1) targets.push(tEnd * Math.pow(1e-8, 1 - i / (count - 1)));
    let t = 0; let targetIndex = 0;
    const dtInit = Math.min(1e-12, tEnd / 1e6);
    const phase1 = Math.min(tEnd, dtInit * 1e3);
    const phase2 = Math.min(tEnd, dtInit * 1e6);
    while (t < tEnd - 1e-18) {
      while (targetIndex < targets.length && t >= targets[targetIndex]) {
        points.push({ t, state: state.slice() }); targetIndex += 1;
      }
      const dt = t < phase1 ? Math.min(dtInit, phase1 - t) : t < phase2 ? Math.min(dtInit * 1e3, phase2 - t) : Math.min(dtInit * 1e4, tEnd - t);
      const next = step(state, rho, dt, true); for (let i = 0; i < state.length; i += 1) state[i] = next[i]; t += dt;
    }
    points.push({ t: tEnd, state: state.slice() });
    const history = points.map((point) => {
      const Y = {}; SPECIES.forEach((s, i) => { Y[s] = point.state[i]; });
      const ySum = SPECIES.reduce((sum, s) => sum + Y[s], 0); SPECIES.forEach((s) => { Y[s] /= ySum; });
      const Tv = temperatureFromEnergy(point.state[6], Y, rho, point.state[7]);
      let invW = 0; SPECIES.forEach((s) => { invW += Y[s] / W[s]; }); const Wmix = 1 / invW;
      const X = {}; const concentration = {}; SPECIES.forEach((s) => { X[s] = Y[s] * Wmix / W[s]; concentration[s] = rho * Y[s] / W[s] * 1000; });
      return {
        t: point.t, T: point.state[5], Tv, X, Y, concentration,
        ...Object.fromEntries(SPECIES.map((s) => [`X_${s}`, X[s]])),
        pressure: rho * (R_UNIVERSAL / Wmix) * point.state[5],
      };
    });
    return { model: "Park 5", species: SPECIES, history, final: history[history.length - 1] };
  }

  const PARK11_SPECIES = ["N2", "O2", "NO", "N2+", "O2+", "NO+", "N", "O", "N+", "O+", "e-"];
  const PARK11_W = { N2: 28.0134, O2: 31.9988, NO: 30.0061, "N2+": 28.0129, "O2+": 31.9983, "NO+": 30.0056, N: 14.0067, O: 15.9994, "N+": 14.0062, "O+": 15.9989, "e-": 5.4858e-4 };
  const PARK11_HF = { N2: 0, O2: 0, NO: 10981, "N2+": 182953, "O2+": 140855, "NO+": 118370, N: 56852, O: 29968, "N+": 225730, "O+": 187996, "e-": 0 };
  const PARK11_THETA = { N2: 3371, O2: 2256, NO: 2719, "N2+": 3371, "O2+": 2256, "NO+": 2719 };
  const PARK11_DISS = { N2: 3.36e7, O2: 1.54e7, NO: 2.09e7, "N2+": 3e7, "O2+": 2.01e7, "NO+": 3.49e7 };
  const PARK11_VIB = Object.keys(PARK11_THETA);
  const PARK11_CV = {};
  PARK11_SPECIES.forEach((s) => { PARK11_CV[s] = (PARK11_THETA[s] ? 2.5 : 1.5) * R_UNIVERSAL / PARK11_W[s]; });
  PARK11_CV["e-"] = 1.5 * R_UNIVERSAL / PARK11_W["e-"];
  const PARK11_REACTIONS = [
    { lhs: { O2: 1 }, rhs: { O: 2 }, A: 2e18, beta: -1.5, Ta: 59500, kind: "diss" },
    { lhs: { N2: 1 }, rhs: { N: 2 }, A: 7e18, beta: -1.6, Ta: 113220, kind: "diss" },
    { lhs: { NO: 1 }, rhs: { N: 1, O: 1 }, A: 5e12, beta: 0, Ta: 75500, kind: "diss" },
    { lhs: { NO: 1, O: 1 }, rhs: { O2: 1, N: 1 }, A: 8.4e9, beta: 0, Ta: 19450, kind: "thermal" },
    { lhs: { N2: 1, O: 1 }, rhs: { NO: 1, N: 1 }, A: 6.4e14, beta: -1, Ta: 38400, kind: "thermal" },
    { lhs: { N: 1, "e-": 1 }, rhs: { "N+": 1, "e-": 2 }, A: 2.5e31, beta: -3.82, Ta: 168600, kind: "ion" },
    { lhs: { O: 1, "e-": 1 }, rhs: { "O+": 1, "e-": 2 }, A: 3.9e30, beta: -3.78, Ta: 158500, kind: "ion" },
    { lhs: { N: 1, O: 1 }, rhs: { "NO+": 1, "e-": 1 }, A: 5.3e9, beta: 0, Ta: 31900, kind: "assoc" },
    { lhs: { O: 1, "N+": 1 }, rhs: { NO: 1, N: 1 }, A: 2e10, beta: 0, Ta: 25000, kind: "thermal" },
  ];

  function genericKeq(reaction, species, weights, formation, temperature, pressure = 101325) {
    if (reaction.ni && reaction.A0 && reaction.A1 && reaction.A2 && reaction.A3 && reaction.A4) {
      return referenceKeq(reaction, temperature, pressure);
    }
    let deltaH = 0;
    Object.entries(reaction.rhs).forEach(([s, nu]) => { deltaH += nu * R_UNIVERSAL * formation[s] / weights[s]; });
    Object.entries(reaction.lhs).forEach(([s, nu]) => { deltaH -= nu * R_UNIVERSAL * formation[s] / weights[s]; });
    return Math.exp(clamp(-deltaH / (R_UNIVERSAL * Math.max(temperature, 1)), -500, 500));
  }
  function referenceKeq(reaction, temperature, pressure) {
    const T = Math.max(temperature, 1);
    const numberDensity = pressure / (K_B * T);
    const z = 1e4 / T;
    const coefficient = (values) => interpolate(values, numberDensity, reaction.ni);
    return Math.exp(coefficient(reaction.A0) / z + coefficient(reaction.A1)
      + coefficient(reaction.A2) * Math.log(z) + coefficient(reaction.A3) * z
      + coefficient(reaction.A4) * z * z);
  }
  function normalizeReferenceReaction(reaction) {
    const type = reaction.type;
    const kind = type === "dissociation" ? "diss"
      : type === "impactDissociation" ? "impactDiss"
        : type === "impactIonisation" ? "ion"
          : type === "associativeIonisation" ? "assoc"
              : type === "chargeExchange" ? "chargeExchange"
                : type === "dissociation_QK" ? "qkDiss"
              : type === "exchange_irr" ? "qkExchange" : "thermal";
    return { ...reaction, kind };
  }
  function referenceReactionSet(table, species) {
    return Object.values(table || {}).map(normalizeReferenceReaction).filter((reaction) => {
      const required = new Set([...Object.keys(reaction.lhs), ...Object.keys(reaction.rhs)]);
      return [...required].every((name) => species.includes(name));
    });
  }
  function genericProducts(concentrations, stoich) {
    return Object.entries(stoich).reduce((value, [s, power]) => value * Math.pow(Math.max(concentrations[s], 0), power), 1);
  }

  const KNAB_VV_DEFAULT = [0.01, 1.5e-18];
  const KNAB_VV = {
    "N2|O2": [0.01, 2.667e-19],
    "N2|NO": [0.01, 2.751e-19],
    "O2|NO": [0.01, 2.686e-19],
  };
  function knabCoeff(a, b) {
    return KNAB_VV[`${a}|${b}`] || KNAB_VV[`${b}|${a}`] || KNAB_VV_DEFAULT;
  }
  function computeKnabQVV(T, Tv, Y, nD, rho, vibSpecies, weights) {
    let source = 0;
    vibSpecies.forEach((i) => {
      const rhoI = rho * (Y[i] || 0);
      const evIT = vibEnergy(T, i);
      const evITv = vibEnergy(Tv, i);
      vibSpecies.forEach((j) => {
        if (i === j) return;
        const [p21, sigma12] = knabCoeff(i, j);
        const reducedWeight = weights[i] * weights[j] / (weights[i] + weights[j]);
        const relativeSpeed = Math.sqrt(8 * R_UNIVERSAL * Math.max(T, 1) / (Math.PI * reducedWeight));
        const evJT = vibEnergy(T, j);
        const evJTv = vibEnergy(Tv, j);
        const deltaEv = evIT / Math.max(evJT, 1e-300) * evJTv - evITv;
        source += rhoI * (nD[j] || 0) * sigma12 * p21 * relativeSpeed * deltaEv;
      });
    });
    return source;
  }

  function computeReferenceQVT(T, Tv, pressure, Y, c, nD, rho, species, vibSpecies, weights) {
    const tables = globalThis.REACTING_AIR_REFERENCE_TABLES;
    if (!tables || !tables.vt) return null;
    const constants = tables.vt_constants;
    const allSpecies = ["N2", "O2", "NO", "N2+", "O2+", "NO+", "N", "O", "N+", "O+", "e-"];
    const totalC = species.reduce((sum, name) => sum + c[name], 0);
    const electronFraction = (c["e-"] || 0) / Math.max(totalC, 1e-300);
    let source = 0;
    vibSpecies.forEach((vibrating) => {
      const evT = vibEnergy(T, vibrating);
      const evTv = vibEnergy(Tv, vibrating);
      let mixedRate = 0;
      species.forEach((collider) => {
        if (collider === "e-") return;
        const keySpecies = allSpecies.indexOf(collider) < allSpecies.indexOf(vibrating) ? collider : vibrating;
        const keyCollider = allSpecies.indexOf(collider) < allSpecies.indexOf(vibrating) ? vibrating : collider;
        const pair = tables.vt[`['${keySpecies}', '${keyCollider}']`];
        if (!pair) return;
        const nDcol = vibrating === collider ? nD[vibrating] : nD[vibrating] + nD[collider];
        const tauMW = (constants.p_atm / Math.max(pressure, 1)) * Math.exp(pair[0] * (Math.pow(T, -1 / 3) - pair[1]) - constants.offset);
        const vBar = Math.sqrt(8 * R_UNIVERSAL * T / (Math.PI * weights[keySpecies]));
        const parkCollision = vBar * constants.sigma1 * Math.pow(constants.sigma2 / T, 2) * nDcol;
        const tau = tauMW + 1 / Math.max(parkCollision, 1e-300);
        mixedRate += (c[collider] / Math.max(totalC, 1e-300)) / tau;
      });
      const tauMix = (1 - electronFraction) / Math.max(mixedRate, 1e-300);
      source += rho * Y[vibrating] * (evT - evTv) / tauMix;
    });
    return source;
  }
  function genericRate(reaction, T, Tv, pressure, species, weights, formation, mmt = false) {
    if (!mmt) {
      if (reaction.kind === "qkExchange") {
        return { kf: arrh(reaction.A, reaction.beta, reaction.Ta, T), kr: 0 };
      }
      if (reaction.kind === "qkDiss") return { kf: 0, kr: 0 };
      const parkTemperature = Math.pow(T, .7) * Math.pow(Math.max(Tv, 1), .3);
      if (reaction.kind === "diss") {
        const kf = arrh(reaction.A, reaction.beta, reaction.Ta, parkTemperature);
        const kr = arrh(reaction.A, reaction.beta, reaction.Ta, T) / Math.max(genericKeq(reaction, species, weights, formation, T, pressure), 1e-300);
        return { kf, kr };
      }
      if (reaction.kind === "impactDiss") {
        const kf = arrh(reaction.A, reaction.beta, reaction.Ta, parkTemperature);
        const kr = arrh(reaction.A, reaction.beta, reaction.Ta, Tv) / Math.max(genericKeq(reaction, species, weights, formation, Tv, pressure), 1e-300);
        return { kf, kr };
      }
      if (reaction.kind === "ion") {
        const kf = arrh(reaction.A, reaction.beta, reaction.Ta, Tv);
        return { kf, kr: kf / Math.max(genericKeq(reaction, species, weights, formation, Tv, pressure), 1e-300) };
      }
      if (reaction.kind === "assoc") {
        const kf = arrh(reaction.A, reaction.beta, reaction.Ta, T);
        const kr = arrh(reaction.A, reaction.beta, reaction.Ta, Tv) / Math.max(genericKeq(reaction, species, weights, formation, Tv, pressure), 1e-300);
        return { kf, kr };
      }
      const kf = arrh(reaction.A, reaction.beta, reaction.Ta, T);
      return { kf, kr: kf / Math.max(genericKeq(reaction, species, weights, formation, T, pressure), 1e-300) };
    }
    const diss = reaction.diss;
    const theta = { N2: 3414, O2: 2251, NO: 2744 }[diss];
    const td = { N2: 113200, O2: 59330, NO: 75360 }[diss];
    const qv = (value) => Math.expm1(-td / Math.max(value, 1)) / Math.expm1(-theta / Math.max(value, 1));
    const kArr = reaction.A * 1e-6 * Math.pow(T, reaction.beta) * Math.exp(-reaction.Ta / Math.max(T, 1));
    const U = 1 / (reaction.aU / Math.max(T, 1) + 1 / reaction.Ustar);
    const denominator = 1 / Math.max(Tv, 1) - 1 / Math.max(T, 1) - 1 / U;
    const Tf = Math.abs(denominator) > 1e-12 ? 1 / denominator : T;
    const correctionDenominator = qv(Tv) * qv(-U);
    const correction = Number.isFinite(correctionDenominator) && Math.abs(correctionDenominator) > 1e-300
      ? qv(T) * qv(Tf) / correctionDenominator : 1;
    const kf = Number.isFinite(correction) ? Math.max(0, kArr * correction * .5) : kArr * .5;
    const kr = kArr * .5 / Math.max(genericKeq(reaction, species, weights, formation, T, pressure), 1e-300);
    return { kf, kr };
  }
  function runGenericFiniteRate(input, config) {
    const species = config.species;
    const T0 = Number(input.temperature); const Tv0 = Number(input.vibrationalTemperature);
    const pressure0 = Number(input.pressure); const tEnd = Number(input.endTime);
    if (!(T0 > 0 && Tv0 > 0 && pressure0 > 0 && tEnd > 0)) throw new Error("Temperatures, pressure, and end time must be positive.");
    const xN2 = Math.max(Number(input.n2Feed), 0); const xO2 = Math.max(Number(input.o2Feed), 0);
    const supplied = input.composition || {};
    const X0 = {}; species.forEach((s) => { X0[s] = Math.max(Number(supplied[s] ?? 0), 0); });
    X0.N2 = Math.max(Number(supplied.N2 ?? xN2), 0); X0.O2 = Math.max(Number(supplied.O2 ?? xO2), 0);
    const compositionSum = species.reduce((sum, s) => sum + X0[s], 0);
    if (!(compositionSum > 0)) throw new Error("Initial species composition must contain a positive total.");
    species.forEach((s) => { X0[s] /= compositionSum; });
    const Wmix0 = species.reduce((sumW, s) => sumW + X0[s] * config.weights[s], 0);
    const rho = pressure0 / ((R_UNIVERSAL / Wmix0) * T0);
    const state = species.map((s) => X0[s] * config.weights[s] / Wmix0);
    const energy = config.vibSpecies.reduce((sum, s) => sum + rho * state[species.indexOf(s)] * vibEnergy(Tv0, s), 0);
    state.push(T0, energy, Tv0);
    const points = []; const targetCount = 260; const targets = [];
    for (let i = 0; i < targetCount; i += 1) targets.push(tEnd * Math.pow(1e-8, 1 - i / (targetCount - 1)));
    let t = 0; let targetIndex = 0;
    const dtInit = Math.min(1e-12, tEnd / 1e6);
    const phase1 = Math.min(tEnd, dtInit * 1e3);
    const phase2 = Math.min(tEnd, dtInit * 1e6);
    while (t < tEnd - 1e-18) {
      while (targetIndex < targets.length && t >= targets[targetIndex]) { points.push({ t, state: state.slice() }); targetIndex += 1; }
      const dt = t < phase1 ? Math.min(dtInit, phase1 - t)
        : t < phase2 ? Math.min(dtInit * 1e3, phase2 - t)
          : Math.min(dtInit * 1e4, tEnd - t);
      const Y = {}; species.forEach((s, i) => { Y[s] = Math.max(state[i], 0); }); const ySum = species.reduce((sum, s) => sum + Y[s], 0); species.forEach((s) => { Y[s] /= ySum; });
      const T = clamp(state[species.length], 300, 50000); const Tv = temperatureFromEnergy(Math.max(state[species.length + 1], 0), Y, rho, state[species.length + 2], config.vibSpecies);
      const c = {}; species.forEach((s) => { c[s] = rho * Y[s] / config.weights[s]; }); const pressure = species.reduce((sum, s) => sum + c[s], 0) * R_UNIVERSAL * T;
      const omega = {}; species.forEach((s) => { omega[s] = 0; });
      config.reactions.forEach((reaction) => {
        if (reaction.kind === "qkDiss") {
          Object.entries(reaction.colliders || {}).forEach(([collider, parameters]) => {
            if (!species.includes(collider)) return;
            const parkTemperature = Math.pow(T, .7) * Math.pow(Math.max(Tv, 1), .3);
            const net = arrh(parameters[0], parameters[1], parameters[2], parkTemperature)
              * genericProducts(c, reaction.lhs) * c[collider];
            Object.entries(reaction.lhs).forEach(([s, nu]) => { omega[s] -= nu * net; });
            Object.entries(reaction.rhs).forEach(([s, nu]) => { omega[s] += nu * net; });
          });
          return;
        }
        const rates = genericRate(reaction, T, Tv, pressure, species, config.weights, config.formation, config.mmt);
        let net = rates.kf * genericProducts(c, reaction.lhs) - rates.kr * genericProducts(c, reaction.rhs);
        if (reaction.kind === "diss" && reaction.eta) {
          net *= species.reduce((sum, name) => sum + (reaction.eta[name] || 0) * c[name], 0);
        }
        Object.entries(reaction.lhs).forEach(([s, nu]) => { omega[s] -= nu * net; });
        Object.entries(reaction.rhs).forEach(([s, nu]) => { omega[s] += nu * net; });
      });
      const rr = {}; species.forEach((s) => { rr[s] = omega[s] * config.weights[s]; });
      const nD = Object.fromEntries(species.map((s) => [s, c[s] * N_A]));
      const referenceQVT = computeReferenceQVT(T, Tv, pressure, Y, c, nD, rho, species, config.vibSpecies, config.weights);
      let qvt = referenceQVT ?? 0;
      if (referenceQVT === null) config.vibSpecies.forEach((s) => { const tau = 1e-7 * Math.max(1, P_ATM / Math.max(pressure, 1)) * Math.exp(80 * (Math.pow(Math.max(T, 300), -1 / 3) - .03)); qvt += rho * Y[s] * (vibEnergy(T, s) - vibEnergy(Tv, s)) / Math.max(tau, 1e-12); });
      const qvv = config.knabVV && input.knabVV ? computeKnabQVV(T, Tv, Y, nD, rho, config.vibSpecies, config.weights) : 0;
      const scv = config.vibSpecies.reduce((sum, s) => sum + .3 * (config.dissociation[s] || 0) * rr[s], 0);
      const sh = -species.reduce((sum, s) => sum + R_UNIVERSAL * config.formation[s] / config.weights[s] * rr[s], 0);
      const cvMix = species.reduce((sum, s) => sum + Y[s] * config.cv[s], 0); const carried = species.reduce((sum, s) => sum + rr[s] * config.cv[s] * T, 0);
      species.forEach((s, i) => { state[i] = Math.max(0, Y[s] + dt * rr[s] / rho); }); state[species.length] = clamp(T + dt * (sh - qvt - scv - carried) / Math.max(rho * cvMix, 1e-30), 300, 50000); state[species.length + 1] = Math.max(0, state[species.length + 1] + dt * (qvt + scv + qvv)); state[species.length + 2] = Tv; t += dt;
    }
    points.push({ t: tEnd, state: state.slice() });
    const history = points.map((point) => { const Y = {}; species.forEach((s, i) => { Y[s] = Math.max(point.state[i], 0); }); const sumY = species.reduce((sum, s) => sum + Y[s], 0); species.forEach((s) => { Y[s] /= sumY; }); const Tv = temperatureFromEnergy(point.state[species.length + 1], Y, rho, point.state[species.length + 2], config.vibSpecies); let inverseW = 0; species.forEach((s) => { inverseW += Y[s] / config.weights[s]; }); const Wmix = 1 / inverseW; const X = {}; const concentration = {}; species.forEach((s) => { X[s] = Y[s] * Wmix / config.weights[s]; concentration[s] = rho * Y[s] / config.weights[s] * 1000; }); return { t: point.t, T: point.state[species.length], Tv, X, Y, concentration, ...Object.fromEntries(species.map((s) => [`X_${s}`, X[s]])), pressure: rho * (R_UNIVERSAL / Wmix) * point.state[species.length] }; });
    return { model: config.label, species, history, final: history[history.length - 1] };
  }
  const PARK11_REFERENCE_REACTIONS = referenceReactionSet(globalThis.REACTING_AIR_REFERENCE_TABLES?.park, PARK11_SPECIES);
  const PARK11_CONFIG = { label: "Park 11", species: PARK11_SPECIES, weights: PARK11_W, formation: PARK11_HF, vibSpecies: PARK11_VIB, dissociation: PARK11_DISS, cv: PARK11_CV, reactions: PARK11_REFERENCE_REACTIONS.length ? PARK11_REFERENCE_REACTIONS : PARK11_REACTIONS, mmt: false, knabVV: true };
  const QK_REFERENCE_REACTIONS = referenceReactionSet(globalThis.REACTING_AIR_REFERENCE_TABLES?.qk, SPECIES);
  const QK5_CONFIG = {
    label: "QK 5",
    species: SPECIES,
    weights: W,
    formation: HF,
    vibSpecies: VIB,
    dissociation: DISS,
    cv: CVTR,
    reactions: QK_REFERENCE_REACTIONS.length ? QK_REFERENCE_REACTIONS : REACTIONS,
    mmt: false,
  };
  const MMT_SPECIES = ["N2", "O2", "NO", "N", "O"];
  const MMT_W = { N2: 28.0134, O2: 31.9988, NO: 30.0061, N: 14.0067, O: 15.9994 };
  const MMT_CV = { N2: 2.5 * R_UNIVERSAL / MMT_W.N2, O2: 2.5 * R_UNIVERSAL / MMT_W.O2, NO: 2.5 * R_UNIVERSAL / MMT_W.NO, N: 1.5 * R_UNIVERSAL / MMT_W.N, O: 1.5 * R_UNIVERSAL / MMT_W.O };
  const MMT_REACTIONS = [
    { lhs: { N2: 2 }, rhs: { N2: 1, N: 2 }, A: 1.5259e17, beta: -.40654, Ta: 113200, kind: "diss", diss: "N2", aU: .49725, Ustar: 265482 },
    { lhs: { N2: 1, O2: 1 }, rhs: { N: 2, O2: 1 }, A: 1.2415e18, beta: -.53233, Ta: 59330, kind: "diss", diss: "O2", aU: .71386, Ustar: -107879 },
    { lhs: { NO: 1, N: 1 }, rhs: { N: 2, O: 1 }, A: 2.5759e17, beta: -.59741, Ta: 75360, kind: "diss", diss: "NO", aU: .51679, Ustar: 3876179 },
    { lhs: { N2: 1, O: 1 }, rhs: { NO: 1, N: 1 }, A: 1.6339e11, beta: .76972, Ta: 37850, kind: "thermal" },
    { lhs: { NO: 1, O: 1 }, rhs: { O2: 1, N: 1 }, A: 1.99183e13, beta: .10007, Ta: 24427, kind: "thermal" },
    { lhs: { N2: 1, O2: 1 }, rhs: { NO: 2 }, A: 1.15292e18, beta: -.83459, Ta: 97485, kind: "thermal" },
  ];
  const MMT_CONFIG = { label: "MMT 5", species: MMT_SPECIES, weights: MMT_W, formation: { N2: 0, O2: 0, NO: 10981, N: 56852, O: 29968 }, vibSpecies: ["N2", "O2", "NO"], dissociation: { N2: 3.36e7, O2: 1.54e7, NO: 2.09e7 }, cv: MMT_CV, reactions: MMT_REACTIONS, mmt: true };


  const MMT_FULL = {"NASA9":{"N2":[[22103.71497,-381.846182,6.08273836,-0.00853091441,1.384646189e-05,-9.62579362e-09,2.519705809e-12,710.846086,-10.76003744],[587712.406,-2239.249073,6.06694922,-0.00061396855,1.491806679e-07,-1.923105485e-11,1.061954386e-15,12832.10415,-15.86640027],[831013916.0,-642073.354,202.0264635,-0.03065092046,2.486903333e-06,-9.70595411e-11,1.437538881e-15,4938707.04,-1672.09974]],"O2":[[-34255.6342,484.700097,1.119010961,0.00429388924,-6.83630052e-07,-2.0233727e-09,1.039040018e-12,-3391.45487,18.4969947],[-1037939.022,2344.830282,1.819732036,0.001267847582,-2.188067988e-07,2.053719572e-11,-8.19346705e-16,-16890.10929,17.38716506],[497529430.0,-286610.6874,66.9035225,-0.00616995902,3.016396027e-07,-7.4214166e-12,7.27817577e-17,2293554.027,-553.062161]],"NO":[[-11439.16503,153.6467592,3.43146873,-0.002668592368,8.48139912e-06,-7.68511105e-09,2.386797655e-12,9098.21441,6.72872549],[223901.8716,-1289.651623,5.43393603,-0.00036560349,9.88096645e-08,-1.416076856e-11,9.38018462e-16,17503.17656,-8.50166909],[-957530354.0,591243.448,-138.4566826,0.01694339403,-1.007351096e-06,2.912584076e-11,-3.29510935e-16,-4677501.24,1242.081216]],"N":[[0.0,0.0,2.5,0.0,0.0,0.0,0.0,56104.6378,4.193905036],[88765.0138,-107.12315,2.362188287,0.0002916720081,-1.7295151e-07,4.01265788e-11,-2.677227571e-15,56973.5133,4.865231506],[547518105.0,-310757.498,69.1678274,-0.00684798813,3.8275724e-07,-1.098367709e-11,1.277986024e-16,2550585.618,-584.8769753]],"O":[[-7953.6113,160.7177787,1.966226438,0.00101367031,-1.110415423e-06,6.5175075e-10,-1.584779251e-13,28403.62437,8.40424182],[261902.0262,-729.872203,3.31717727,-0.000428133436,1.036104594e-07,-9.43830433e-12,2.725038297e-16,33924.2806,-0.667958535],[177900426.4,-108232.8257,28.10778365,-0.002975232262,1.854997534e-07,-5.79623154e-12,7.191720164e-17,889094.263,-218.1728151]]},"reactions":[{"name":"2N2<->2N+N2","kind":"MMT","diss":"N2","A":1.5259e+17,"beta":-0.40654,"Ta":113200.0,"aU":0.49725,"Ustar":265482.0,"nu_f":[2.0,0.0,0.0,0.0,0.0],"nu_b":[1.0,0.0,0.0,2.0,0.0]},{"name":"N2+O2<->2N+O2","kind":"MMT","diss":"N2","A":1.6671e+20,"beta":-1.1484,"Ta":113200.0,"aU":0.369,"Ustar":376750.0,"nu_f":[1.0,1.0,0.0,0.0,0.0],"nu_b":[0.0,1.0,0.0,2.0,0.0]},{"name":"N2+NO<->2N+NO","kind":"MMT","diss":"N2","A":1.5259e+17,"beta":-0.40654,"Ta":113200.0,"aU":0.49725,"Ustar":265482.0,"nu_f":[1.0,0.0,1.0,0.0,0.0],"nu_b":[0.0,0.0,1.0,2.0,0.0]},{"name":"N2+N<->3N","kind":"MMT","diss":"N2","A":4.6825e+17,"beta":-0.5127,"Ta":113200.0,"aU":0.39547,"Ustar":805087.0,"nu_f":[1.0,0.0,0.0,1.0,0.0],"nu_b":[0.0,0.0,0.0,3.0,0.0]},{"name":"N2+O<->2N+O","kind":"MMT","diss":"N2","A":3.4678e+17,"beta":-0.50732,"Ta":113200.0,"aU":0.35075,"Ustar":-720429.0,"nu_f":[1.0,0.0,0.0,0.0,1.0],"nu_b":[0.0,0.0,0.0,2.0,1.0]},{"name":"O2+N2<->2O+N2","kind":"MMT","diss":"O2","A":1.2415e+18,"beta":-0.53233,"Ta":59330.0,"aU":0.71386,"Ustar":-107879.0,"nu_f":[1.0,1.0,0.0,0.0,0.0],"nu_b":[1.0,0.0,0.0,0.0,2.0]},{"name":"2O2<->2O+O2","kind":"MMT","diss":"O2","A":3.4529e+18,"beta":-0.59391,"Ta":59330.0,"aU":0.4611,"Ustar":54974.0,"nu_f":[0.0,2.0,0.0,0.0,0.0],"nu_b":[0.0,1.0,0.0,0.0,2.0]},{"name":"O2+NO<->2O+NO","kind":"MMT","diss":"O2","A":3.4529e+18,"beta":-0.59391,"Ta":59330.0,"aU":0.4611,"Ustar":54974.0,"nu_f":[0.0,1.0,1.0,0.0,0.0],"nu_b":[0.0,0.0,1.0,0.0,2.0]},{"name":"O2+N<->2O+N","kind":"MMT","diss":"O2","A":9.3045e+17,"beta":-0.44587,"Ta":59330.0,"aU":0.37965,"Ustar":-313988.0,"nu_f":[0.0,1.0,0.0,1.0,0.0],"nu_b":[0.0,0.0,0.0,1.0,2.0]},{"name":"O2+O<->3O","kind":"MMT","diss":"O2","A":8.2864e+17,"beta":-0.47445,"Ta":59330.0,"aU":0.36157,"Ustar":96143.0,"nu_f":[0.0,1.0,0.0,0.0,1.0],"nu_b":[0.0,0.0,0.0,0.0,3.0]},{"name":"NO+N2<->N+O+N2","kind":"thermal","A":358890000000000.0,"beta":-0.07945,"Ta":52644.0,"nu_f":[1.0,0.0,1.0,0.0,0.0],"nu_b":[1.0,0.0,0.0,1.0,1.0]},{"name":"NO+O2<->N+O+O2","kind":"thermal","A":358890000000000.0,"beta":-0.07945,"Ta":52644.0,"nu_f":[0.0,1.0,1.0,0.0,0.0],"nu_b":[0.0,1.0,0.0,1.0,1.0]},{"name":"2NO<->N+O+NO","kind":"thermal","A":358890000000000.0,"beta":-0.07945,"Ta":52644.0,"nu_f":[0.0,0.0,2.0,0.0,0.0],"nu_b":[0.0,0.0,1.0,1.0,1.0]},{"name":"NO+N<->2N+O","kind":"MMT","diss":"NO","A":2.5759e+17,"beta":-0.59741,"Ta":75360.0,"aU":0.51679,"Ustar":3876179.0,"nu_f":[0.0,0.0,1.0,1.0,0.0],"nu_b":[0.0,0.0,0.0,2.0,1.0]},{"name":"NO+O<->N+2O","kind":"MMT","diss":"NO","A":7.031e+16,"beta":-0.52778,"Ta":75360.0,"aU":0.301,"Ustar":-3007416.0,"nu_f":[0.0,0.0,1.0,0.0,1.0],"nu_b":[0.0,0.0,0.0,1.0,2.0]},{"name":"N2+O<->NO+N","kind":"thermal","A":163390000000.0,"beta":0.76972,"Ta":37850.0,"nu_f":[1.0,0.0,0.0,0.0,1.0],"nu_b":[0.0,0.0,1.0,1.0,0.0]},{"name":"NO+O<->O2+N","kind":"thermal","A":19918300000000.0,"beta":0.10007,"Ta":24427.0,"nu_f":[0.0,0.0,1.0,0.0,1.0],"nu_b":[0.0,1.0,0.0,1.0,0.0]},{"name":"N2+O2<->2NO","kind":"thermal","A":1.15292e+18,"beta":-0.83459,"Ta":97485.0,"nu_f":[1.0,1.0,0.0,0.0,0.0],"nu_b":[0.0,0.0,2.0,0.0,0.0]},{"name":"2NO<->2O+N2","kind":"thermal","A":1.1786e+16,"beta":-0.54728,"Ta":53783.0,"nu_f":[0.0,0.0,2.0,0.0,0.0],"nu_b":[1.0,0.0,0.0,0.0,2.0]},{"name":"2NO<->2N+O2","kind":"thermal","A":1.3688e+21,"beta":-1.6611,"Ta":106199.0,"nu_f":[0.0,0.0,2.0,0.0,0.0],"nu_b":[0.0,1.0,0.0,2.0,0.0]},{"name":"N2+O2<->NO+N+O","kind":"thermal","A":290530000000000.0,"beta":0.15022,"Ta":89282.0,"nu_f":[1.0,1.0,0.0,0.0,0.0],"nu_b":[0.0,0.0,1.0,1.0,1.0]}],"relax":{"N2|N2":[221.0,-24.83,33.3,-17.31],"N2|O2":[115.0,-20.25,-51.42,-15.26],"N2|NO":[221.0,-24.83,33.3,-17.31],"N2|N":[239.4,-28.52,32.34,-18.39],"N2|O":[149.5,-23.51,-1.201,-17.4],"O2|N2":[134.0,-22.37,-0.7123,-16.95],"O2|O2":[129.0,-22.29,-4.384,-16.7],"O2|NO":[129.0,-22.29,-4.384,-16.7],"O2|N":[81.22,-21.16,9.111,-17.61],"O2|O":[3.018,-17.74,-121.6,-15.82],"NO|N2":[127.2,-23.2,-27.61,-16.55],"NO|O2":[127.2,-23.2,-27.61,-16.55],"NO|NO":[127.2,-23.2,-27.61,-16.55],"NO|N":[64.31,-20.06,-62.33,-16.8],"NO|O":[-24.76,-17.16,-113.0,-20.7]},"ctrE":{"N2":20.786156545,"O2":20.786156545,"NO":20.786156545,"N":12.471693927,"O":12.471693927},"ctrH":{"N2":29.100619163,"O2":29.100619163,"NO":29.100619163,"N":20.786156545,"O":20.786156545},"href":{"N2":6.114502307264744e-06,"O2":-1.2807210958026514e-05,"NO":91270.78988350501,"N":472677.3062584215,"O":249173.58263372252},"cconst":{"N2":-8676.349597333947,"O2":-8676.34961625566,"NO":82594.44028005656,"N":466479.9136845298,"O":242976.19005983075}};
  const MMT_R = 8.314462618;
  const MMT_P_STD = 1.0e5;
  const MMT_T_REF = 298.15;
  const MMT_ATM = 101325.0;
  const MMT_FK_NB = 0.5;
  const MMT_FEPS_NB = 0.85;
  const MMT_THETA_V = { N2: 3414.0, O2: 2251.0, NO: 2744.0 };
  const MMT_T_DISS = { N2: 113200.0, O2: 59330.0, NO: 75360.0 };
  const MMT_DIATOMIC = ["N2", "O2", "NO"];
  const MMT_MOLAR_MASS = { N2: 28.0134e-3, O2: 31.9988e-3, NO: 30.0061e-3, N: 14.0067e-3, O: 15.9994e-3 };
  const MMT_IDX = { N2: 0, O2: 1, NO: 2, N: 3, O: 4 };

  function mmtCoeff(species, T) {
    const Tc = clamp(T, 200, 20000);
    const ranges = MMT_FULL.NASA9[species];
    return Tc <= 1000 ? ranges[0] : Tc <= 6000 ? ranges[1] : ranges[2];
  }
  function mmtCpFit(species, T) {
    const [a1, a2, a3, a4, a5, a6, a7] = mmtCoeff(species, T);
    return MMT_R * (a1 / (T * T) + a2 / T + a3 + a4 * T + a5 * T * T + a6 * T * T * T + a7 * T * T * T * T);
  }
  function mmtHFit(species, T) {
    const [a1, a2, a3, a4, a5, a6, a7, b1] = mmtCoeff(species, T);
    return MMT_R * (-a1 / T + a2 * Math.log(T) + a3 * T + a4 / 2 * T * T + a5 / 3 * T ** 3 + a6 / 4 * T ** 4 + a7 / 5 * T ** 5 + b1);
  }
  function mmtSFit(species, T) {
    const [a1, a2, a3, a4, a5, a6, a7, , b2] = mmtCoeff(species, T);
    return MMT_R * (-a1 / 2 * T ** -2 - a2 / T + a3 * Math.log(T) + a4 * T + a5 / 2 * T * T + a6 / 3 * T ** 3 + a7 / 4 * T ** 4 + b2);
  }
  function mmtSafeExpm1(x) { return Math.expm1(clamp(x, -700, 700)); }
  function mmtQv(Tp, TD, ThV) { return mmtSafeExpm1(-TD / Tp) / mmtSafeExpm1(-ThV / Tp); }
  function mmtEvel(species, T) { return mmtHFit(species, T) - MMT_FULL.ctrH[species] * (T - MMT_T_REF) - MMT_FULL.href[species]; }
  function mmtCpV(species, T) { return mmtCpFit(species, T) - MMT_FULL.ctrH[species]; }
  function mmtETwoTemp(species, T, Tv) { return mmtHFit(species, Tv) + MMT_FULL.ctrE[species] * T - MMT_FULL.ctrH[species] * Tv; }
  function mmtSolveTv(C, Ev, guess) {
    let Tv = clamp(guess || 1000, 150, 20000);
    for (let i = 0; i < 60; i += 1) {
      let f = -Ev;
      let df = 0;
      MMT_SPECIES.forEach((s, j) => { f += C[j] * mmtEvel(s, Tv); df += C[j] * mmtCpV(s, Tv); });
      const next = clamp(Tv - f / Math.max(df, 1e-8), 150, 20000);
      if (Math.abs(next - Tv) < 1e-6 * next + 1e-8) return next;
      Tv = next;
    }
    return Tv;
  }
  function mmtTemperature(C, Ev, E0) {
    let denom = 0;
    let constant = 0;
    MMT_SPECIES.forEach((s, i) => { denom += C[i] * MMT_FULL.ctrE[s]; constant += C[i] * MMT_FULL.cconst[s]; });
    return clamp((E0 - Ev - constant) / Math.max(denom, 1e-30), 100, 20000);
  }
  function mmtEquilibriumConstant(rxn, T) {
    let nuT = 0;
    let acc = 0;
    MMT_SPECIES.forEach((s, i) => {
      const nu = rxn.nu_b[i] - rxn.nu_f[i];
      nuT += nu;
      if (nu !== 0) acc += nu * (mmtHFit(s, T) - T * mmtSFit(s, T)) / (MMT_R * T);
    });
    return Math.pow(MMT_P_STD, nuT) * Math.exp(clamp(-acc, -700, 700)) * Math.pow(1 / (MMT_R * T), nuT);
  }
  function mmtRateConsts(rxn, T, Tv) {
    let kf;
    let kfEq;
    let TF = null;
    const kArr = rxn.A * 1e-6 * Math.pow(T, rxn.beta) * Math.exp(-rxn.Ta / Math.max(T, 1));
    if (rxn.kind === "MMT") {
      const theta = MMT_THETA_V[rxn.diss];
      const td = MMT_T_DISS[rxn.diss];
      let denomU = rxn.aU / Math.max(T, 1) + 1 / rxn.Ustar;
      if (Math.abs(denomU) < 1e-12) denomU = denomU < 0 ? -1e-12 : 1e-12;
      const U = 1 / denomU;
      let denomF = 1 / Math.max(Tv, 1) - 1 / Math.max(T, 1) - 1 / U;
      if (Math.abs(denomF) < 1e-12) denomF = denomF < 0 ? -1e-12 : 1e-12;
      TF = 1 / denomF;
      const zDen = mmtQv(Tv, td, theta) * mmtQv(-U, td, theta);
      const z = Math.abs(zDen) > 1e-300 ? (mmtQv(T, td, theta) * mmtQv(TF, td, theta)) / zDen : 1;
      kf = kArr * (Number.isFinite(z) ? z : 1) * MMT_FK_NB;
      kfEq = kArr * MMT_FK_NB;
    } else {
      kf = kArr;
      kfEq = kArr;
    }
    const kb = kfEq / Math.max(mmtEquilibriumConstant(rxn, T), 1e-300);
    return { kf: Number.isFinite(kf) ? Math.max(kf, 0) : 0, kb: Number.isFinite(kb) ? Math.max(kb, 0) : 0, TF };
  }
  function mmtTauPair(s, q, T, pAtm) {
    const [mlo, nlo, mhi, nhi] = MMT_FULL.relax[`${s}|${q}`];
    const tm13 = Math.pow(T, -1 / 3);
    return (Math.exp(clamp(mlo * tm13 + nlo, -700, 700)) + Math.exp(clamp(mhi * tm13 + nhi, -700, 700))) / Math.max(pAtm, 1e-8);
  }
  function mmtCombinedRelaxTime(C, T) {
    const ctot = C.reduce((a, b) => a + b, 0);
    const x = C.map((v) => v / Math.max(ctot, 1e-300));
    const pAtm = ctot * MMT_R * T / MMT_ATM;
    const tau = {};
    MMT_DIATOMIC.forEach((s) => {
      let inv = 0;
      MMT_SPECIES.forEach((q, i) => { inv += x[i] / Math.max(mmtTauPair(s, q, T, pAtm), 1e-300); });
      tau[s] = 1 / Math.max(inv, 1e-300);
    });
    let numer = 0;
    let denom = 0;
    MMT_DIATOMIC.forEach((s) => { numer += x[MMT_IDX[s]]; denom += x[MMT_IDX[s]] / Math.max(tau[s], 1e-300); });
    return numer / Math.max(denom, 1e-300);
  }
  function mmtProduct(C, nu) { return nu.reduce((v, power, i) => v * Math.pow(Math.max(C[i], 1e-30), power), 1); }
  function mmtRhs(Cin, Ev, E0, tvGuess) {
    const C = Cin.map((v) => Math.max(v, 1e-30));
    const Tv = mmtSolveTv(C, Ev, tvGuess);
    const T = mmtTemperature(C, Ev, E0);
    const dC = [0, 0, 0, 0, 0];
    let wchem = 0;
    MMT_FULL.reactions.forEach((rxn) => {
      const { kf, kb, TF } = mmtRateConsts(rxn, T, Tv);
      const Rr = kf * mmtProduct(C, rxn.nu_f) - kb * mmtProduct(C, rxn.nu_b);
      for (let i = 0; i < 5; i += 1) dC[i] += (rxn.nu_b[i] - rxn.nu_f[i]) * Rr;
      MMT_DIATOMIC.forEach((s) => {
        const idx = MMT_IDX[s];
        const net = rxn.nu_b[idx] - rxn.nu_f[idx];
        if (net === 0) return;
        let eps;
        if (rxn.kind === "MMT" && s === rxn.diss) {
          const theta = MMT_THETA_V[s];
          const td = MMT_T_DISS[s];
          eps = MMT_FEPS_NB * MMT_R * (theta / mmtSafeExpm1(theta / TF) - td / mmtSafeExpm1(td / TF));
        } else {
          eps = mmtEvel(s, Tv);
        }
        wchem += Rr * net * eps;
      });
    });
    const tauV = mmtCombinedRelaxTime(C, T);
    let wrelax = 0;
    MMT_SPECIES.forEach((s, i) => { wrelax += C[i] * (mmtEvel(s, T) - mmtEvel(s, Tv)); });
    return { dC, dEv: wrelax / Math.max(tauV, 1e-300) + wchem, T, Tv };
  }
  function runMMTFull(input) {
    const T0 = Number(input.temperature); const Tv0 = Number(input.vibrationalTemperature);
    const pressure0 = Number(input.pressure); const tEnd = Number(input.endTime);
    if (!(T0 > 0 && Tv0 > 0 && pressure0 > 0 && tEnd > 0)) throw new Error("Temperatures, pressure, and end time must be positive.");
    const supplied = input.composition || {};
    const X0 = {}; MMT_SPECIES.forEach((s) => { X0[s] = Math.max(Number(supplied[s] ?? 0), 0); });
    X0.N2 = Math.max(Number(supplied.N2 ?? input.n2Feed), 0); X0.O2 = Math.max(Number(supplied.O2 ?? input.o2Feed), 0);
    const xsum = MMT_SPECIES.reduce((sum, s) => sum + X0[s], 0);
    if (!(xsum > 0)) throw new Error("Initial MMT composition must contain a positive neutral mole fraction.");
    MMT_SPECIES.forEach((s) => { X0[s] /= xsum; });
    const ctot0 = pressure0 / (MMT_R * T0);
    const C = MMT_SPECIES.map((s) => X0[s] * ctot0);
    let Ev = MMT_SPECIES.reduce((sum, s, i) => sum + C[i] * mmtEvel(s, Tv0), 0);
    const E0 = MMT_SPECIES.reduce((sum, s, i) => sum + C[i] * mmtETwoTemp(s, T0, Tv0), 0);
    const points = []; const targets = [];
    for (let i = 0; i < 300; i += 1) targets.push(tEnd * Math.pow(1e-8, 1 - i / 299));
    let t = 0; let targetIndex = 0; let tvGuess = Tv0;
    const dtInit = Math.min(1e-12, tEnd / 1e6);
    const phase1 = Math.min(tEnd, dtInit * 1e3);
    const phase2 = Math.min(tEnd, dtInit * 1e6);
    while (t < tEnd - 1e-18) {
      while (targetIndex < targets.length && t >= targets[targetIndex]) { points.push({ t, C: C.slice(), Ev, tvGuess }); targetIndex += 1; }
      let dt = t < phase1 ? Math.min(dtInit, phase1 - t) : t < phase2 ? Math.min(dtInit * 1e3, phase2 - t) : Math.min(dtInit * 1e4, tEnd - t);
      const rhs = mmtRhs(C, Ev, E0, tvGuess);
      tvGuess = rhs.Tv;
      for (let i = 0; i < 5; i += 1) if (rhs.dC[i] < 0) dt = Math.min(dt, 0.2 * C[i] / -rhs.dC[i]);
      dt = Math.max(dt, Math.min(dtInit, tEnd - t));
      for (let i = 0; i < 5; i += 1) C[i] = Math.max(0, C[i] + dt * rhs.dC[i]);
      Ev += dt * rhs.dEv;
      t += dt;
    }
    points.push({ t: tEnd, C: C.slice(), Ev, tvGuess });
    const history = points.map((point) => {
      const Ck = point.C.map((v) => Math.max(v, 0));
      const ctot = Ck.reduce((a, b) => a + b, 0);
      const Tv = mmtSolveTv(Ck, point.Ev, point.tvGuess);
      const T = mmtTemperature(Ck, point.Ev, E0);
      const X = {}; const Y = {}; const concentration = {};
      MMT_SPECIES.forEach((s, i) => { X[s] = Ck[i] / Math.max(ctot, 1e-300); concentration[s] = Ck[i]; });
      const Wmix = MMT_SPECIES.reduce((sumW, s) => sumW + X[s] * MMT_MOLAR_MASS[s], 0);
      MMT_SPECIES.forEach((s) => { Y[s] = X[s] * MMT_MOLAR_MASS[s] / Math.max(Wmix, 1e-300); });
      return { t: point.t, T, Tv, X, Y, concentration, ...Object.fromEntries(MMT_SPECIES.map((s) => [`X_${s}`, X[s]])), pressure: ctot * MMT_R * T };
    });
    return { model: "MMT 5", species: MMT_SPECIES, history, final: history[history.length - 1] };
  }

  const api = { runPark5, runPark11: (input) => runGenericFiniteRate(input, PARK11_CONFIG), runQK5: (input) => runGenericFiniteRate(input, QK5_CONFIG), runMMT: runMMTFull };
  if (typeof globalThis !== "undefined") globalThis.ReactingAirModels = api;
  if (typeof document === "undefined" && typeof self !== "undefined") {
    self.onmessage = (event) => {
      try {
        const method = event.data.method || "runPark5";
        self.postMessage({ ok: true, result: api[method](event.data.input || event.data) });
      }
      catch (error) { self.postMessage({ ok: false, error: error.message }); }
    };
  }
}());
