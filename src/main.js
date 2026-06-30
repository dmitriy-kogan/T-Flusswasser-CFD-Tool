import * as THREE from 'three';
import { OrbitControls } from 'three/examples/jsm/controls/OrbitControls.js';

"use strict";

/* ---------- Sichtbare Fehleranzeige (statt stiller schwarzer Szene) ------ */
function showError(msg) {
  const e = document.getElementById('err');
  if (!e) return;
  e.className = 'show';
  e.textContent = '⚠ Fehler:\n' + msg;
}
window.addEventListener('error', ev => {
  showError((ev.message || '') + '\n' + (ev.filename || '') + ':' + (ev.lineno || '') +
    (ev.error && ev.error.stack ? '\n\n' + ev.error.stack : ''));
});
window.addEventListener('unhandledrejection', ev => {
  showError('Promise: ' + (ev.reason && ev.reason.message ? ev.reason.message : ev.reason));
});

/* =========================================================================
   FLUSSWASSERWÄRMEPUMPE – 3D-CFD-TOOL
   Aufbau:  (1) Parameter   (2) Rechengitter & Solver
            (3) Geometrie/Szene  (4) Visualisierung (Heatmaps, Partikel)
            (5) UI & abgeleitete Größen  (6) Diagramme  (7) Animationsschleife
   ========================================================================= */

/* ---------- (1) PARAMETER ------------------------------------------------ */
const P = {
  vFlow: 0.5,   // mittlere Strömungsgeschwindigkeit [m/s]
  depth: 5.0,   // maximale Wassertiefe (Flussmitte) [m]
  tIn: 8.0,     // Temperatur des gesamten Flusswassers [°C]
  power: 500,   // thermische Leistung [kW]  (Regler bis 20 000 kW = 20 MW)
  vSuct: 1.5,   // Ansauggeschwindigkeit im Rohr [m/s]
  dPipe: 0.8,   // Rohrdurchmesser [m]
  xCut: 120,    // Querschnittsposition [m]
  zCut: 1.5,    // Horizontalschnitt-Tiefe [m unter Oberfläche]
  yLong: 6      // Längsschnitt y-Position (quer zur Strömung) [m]
};

/* Geometrie des Flussabschnitts */
const L = 300, W = 50;                 // Länge (x), Breite (y) [m]
const X_IN = 20, X_OUT = 120;          // Position Einlass / Auslass [m]
const Y_BANK = 6;                      // Abstand der Rohrmündung vom Pumpenufer (y=0) [m]
const RHO = 1000, CP = 4186;           // Wasser: Dichte, spez. Wärmekapazität
const NU = 1.3e-6;                     // kinematische Viskosität [m²/s]
const ALPHA_MOL = 1.4e-7;              // molekulare Temperaturleitfähigkeit [m²/s]
const MINDEPTH = 0.4;                  // Resttiefe an den Ufern [m]
const WSURF = 10;                      // feste Höhe der Wasseroberfläche/Uferlinie [m]

/* Parabolisches Flussprofil: Oberfläche flach und FEST bei z = WSURF.
   Die Sohle ist parabelförmig (tiefste Stelle in der Mitte); beim Erhöhen der
   Tiefe wandert der tiefste Punkt nach unten (die Oberfläche bleibt fest). */
function localDepth(y) {
  const s = 2 * y / W - 1;
  return Math.max(P.depth * (1 - s * s), MINDEPTH);
}
function bedZ(y) { return WSURF - localDepth(y); }          // Sohle wandert mit P.depth nach unten
const ZSURF = () => WSURF;                                  // Wasseroberfläche (flach, fest)
/* Rohrmündungstiefe: ca. 0,5 m unter der Wasseroberfläche (in flachen Zonen begrenzt) */
function pipeZ() { return Math.max(WSURF - 0.5, bedZ(Y_BANK) + 0.2); }

/* Rechengitter (gröber für Echtzeit) */
const NX = 120, NY = 26, NZ = 8;
const dx = L / NX, dy = W / NY;
const NC = NX * NY * NZ;

/* Temperaturfelder (Doppelpuffer) */
let T = new Float32Array(NC);
let T2 = new Float32Array(NC);

const idx = (i, j, k) => i + NX * (j + NY * k);

/* ---------- (2) SOLVER --------------------------------------------------- */
/* abgeleitete Größen, bei jeder Parameteränderung aktualisiert */
const D = {};                          // derived values
function computeDerived() {
  D.area = Math.PI * Math.pow(P.dPipe / 2, 2);     // Rohrquerschnitt [m²]
  D.Qpump = D.area * P.vSuct;                      // Durchfluss [m³/s]
  D.dT = (P.power * 1000) / (RHO * CP * Math.max(D.Qpump, 1e-6)); // Temperaturabsenkung [K]
  D.alpha = ALPHA_MOL + 0.012 * P.vFlow * P.depth; // effektive (turbulente) Temp.-Leitf.
  D.Re = P.vFlow * P.depth / NU;
  D.Pe = P.vFlow * P.depth / D.alpha;
}

/* analytisches Geschwindigkeitsfeld an beliebiger Stelle (x,y,z) [m]
   = Grundströmung (Profil) + Ansaug-Senke + Auslass-Quelle (Strahl) */
function velocityAt(x, y, z, out) {
  const dl = localDepth(y), b = bedZ(y);
  // Grundprofil: parabolisch über Breite (langsamer am Ufer),
  // Potenzgesetz 1/7 über die lokale Tiefe (langsamer an der Sohle)
  const yy = (2 * y / W - 1);
  const fy = 1 - 0.5 * yy * yy;
  const sg = Math.min(Math.max((z - b) / Math.max(dl, 1e-3), 0.001), 1);  // 0=Sohle .. 1=Oberfläche
  const fz = Math.pow(sg, 1 / 7);
  let u = (z < b - 0.2) ? 0 : P.vFlow * fy * fz * 1.15;   // unter der Sohle keine Strömung
  let v = 0, w = 0;

  const zp = pipeZ();
  // Ansaug-Senke am Einlass – schwache, lokale Störung (erzeugt keine tote Spur)
  const ix = X_IN - x, iy = Y_BANK - y, iz = zp - z;
  const r2i = ix * ix + iy * iy + iz * iz + 8;
  const sink = Math.min(D.Qpump * 0.18 / r2i, 0.12);
  u += sink * ix; v += sink * iy * 0.4; w += sink * iz;

  // Auslass-Quelle – ebenfalls schwach (der sichtbare Einleitstrahl entsteht über
  // dedizierte Auslass-Partikel, nicht über ein starkes Geschwindigkeitsfeld)
  const ox = x - X_OUT, oy = y - Y_BANK, oz = z - zp;
  const r2o = ox * ox + oy * oy + oz * oz + 8;
  const src = Math.min(D.Qpump * 0.18 / r2o, 0.12);
  u += src * ox; v += src * oy * 0.4; w += src * oz;

  out.u = u; out.v = v; out.w = w;
}
const _vel = { u: 0, v: 0, w: 0 };

/* trilineare Abtastung des Temperaturfeldes in Sigma-Koordinaten */
function sampleT(x, y, z) {
  if (x < 0) return P.tIn;                       // Einström-Rand
  const dl = localDepth(y), b = bedZ(y);
  const sg = Math.min(Math.max((z - b) / Math.max(dl, 1e-3), 0), 1);
  let gx = x / dx - 0.5, gy = y / dy - 0.5, gz = sg * NZ - 0.5;
  gx = Math.min(Math.max(gx, 0), NX - 1.001);
  gy = Math.min(Math.max(gy, 0), NY - 1.001);
  gz = Math.min(Math.max(gz, 0), NZ - 1.001);
  const i0 = gx | 0, j0 = gy | 0, k0 = gz | 0;
  const fx = gx - i0, fy = gy - j0, fz = gz - k0;
  const i1 = i0 + 1, j1 = j0 + 1, k1 = k0 + 1;
  const c00 = T[idx(i0, j0, k0)] * (1 - fx) + T[idx(i1, j0, k0)] * fx;
  const c10 = T[idx(i0, j1, k0)] * (1 - fx) + T[idx(i1, j1, k0)] * fx;
  const c01 = T[idx(i0, j0, k1)] * (1 - fx) + T[idx(i1, j0, k1)] * fx;
  const c11 = T[idx(i0, j1, k1)] * (1 - fx) + T[idx(i1, j1, k1)] * fx;
  const c0 = c00 * (1 - fy) + c10 * fy, c1 = c01 * (1 - fy) + c11 * fy;
  return c0 * (1 - fz) + c1 * fz;
}
/* physische z-Höhe eines Gitterpunkts (Sigma -> Höhe) */
function cellZ(j, k) {
  const y = (j + 0.5) * dy;
  return bedZ(y) + (k + 0.5) / NZ * localDepth(y);
}

let tInletLocal = P.tIn, tOutlet = P.tIn;  // lokale Entnahme- / Auslasstemperatur

/* ein Solver-Schritt: semi-Lagrange-Advektion + Diffusion + Quellen/RB */
function solveStep(dt) {
  // lokale Entnahmetemperatur (= Flusstemperatur am Einlass) und Auslasstemperatur
  tInletLocal = sampleT(X_IN, Y_BANK, pipeZ());
  tOutlet = Math.max(tInletLocal - D.dT, 0);     // abgekühlt; kein Eis (>= 0 °C)

  // --- Advektion (semi-Lagrange) mit ortsabhängiger Tiefe ---
  for (let k = 0; k < NZ; k++) {
    for (let j = 0; j < NY; j++) {
      const y = (j + 0.5) * dy;
      for (let i = 0; i < NX; i++) {
        const x = (i + 0.5) * dx, z = cellZ(j, k);
        velocityAt(x, y, z, _vel);
        T2[idx(i, j, k)] = sampleT(x - _vel.u * dt, y - _vel.v * dt, z - _vel.w * dt);
      }
    }
  }
  // --- Diffusion (stabile konvexe Mittelung je Richtung) ---
  const wX = 0.05, wY = 0.30, wZ = 0.26;
  for (let k = 0; k < NZ; k++) {
    for (let j = 0; j < NY; j++) {
      for (let i = 0; i < NX; i++) {
        const c = idx(i, j, k);
        let val = T2[c];
        { const m = i > 0 ? T2[idx(i - 1, j, k)] : val, p = i < NX - 1 ? T2[idx(i + 1, j, k)] : val; val = val * (1 - wX) + (m + p) * 0.5 * wX; }
        { const m = j > 0 ? T2[idx(i, j - 1, k)] : val, p = j < NY - 1 ? T2[idx(i, j + 1, k)] : val; val = val * (1 - wY) + (m + p) * 0.5 * wY; }
        { const m = k > 0 ? T2[idx(i, j, k - 1)] : val, p = k < NZ - 1 ? T2[idx(i, j, k + 1)] : val; val = val * (1 - wZ) + (m + p) * 0.5 * wZ; }
        T2[c] = val;
      }
    }
  }

  // --- Quellterme & Randbedingungen ---
  // Einström-Rand x=0: gesamtes Flusswasser hat dieselbe Temperatur
  for (let k = 0; k < NZ; k++) {
    for (let j = 0; j < NY; j++) {
      T2[idx(0, j, k)] = P.tIn;
    }
  }
  // Auslass: Kaltwasser-Einleitung.
  const oi = Math.round(X_OUT / dx), oj = Math.round(Y_BANK / dy);
  const okSg = (pipeZ() - bedZ(Y_BANK)) / Math.max(localDepth(Y_BANK), 1e-6);   // Sigma der Mündung
  const ok = Math.min(NZ - 1, Math.max(0, Math.round(okSg * NZ - 0.5)));
  const zoneA = 5 * dy * 3 * (localDepth(Y_BANK) / NZ);            // Querschnitt der Einleitzone [m²]
  const dilute = D.Qpump / (D.Qpump + P.vFlow * zoneA * 0.55);    // lokaler Kaltwasseranteil (0..1)
  const core = Math.min(0.95, dilute * 1.7);                  // Kerngewicht der Einleitung
  for (let dk = -1; dk <= 1; dk++) {
    for (let dj = -2; dj <= 2; dj++) {
      for (let di = 0; di <= 3; di++) {
        const i = oi + di, j = oj + dj, k = ok + dk;
        if (i < 0 || i >= NX || j < 0 || j >= NY || k < 0 || k >= NZ) continue;
        const wgt = core * Math.exp(-(di * di * 0.12 + dj * dj * 0.5 + dk * dk * 0.9));
        const c = idx(i, j, k);
        T2[c] = T2[c] * (1 - wgt) + tOutlet * wgt;
      }
    }
  }

  const tmp = T; T = T2; T2 = tmp;
}

/* mehrere Schritte bis ~Gleichgewicht laufen lassen */
function relax(steps) {
  const dt = 1.2 * dx / Math.max(P.vFlow, 0.1);  // großzügiges dt (semi-Lagrange)
  for (let s = 0; s < steps; s++) solveStep(dt);
}

/* Farbskala: am physikalischen Bereich verankert (Auslasstemperatur … Flusstemperatur),
   damit die Abkühlung absolut ablesbar ist und das Verblassen der Fahne stromab sichtbar wird. */
let Tmin = 0, Tmax = 10;
function fieldStats() {
  Tmax = P.tIn;                       // Flusswasser = warmes Ende (rot)
  Tmin = Math.max(tOutlet, 0);        // eingeleitetes Wasser = kaltes Ende (blau)
  if (Tmax - Tmin < 0.4) { Tmin = Tmax - 0.4; } // Mindestspanne, falls kaum Abkühlung
}

/* ---------- Farbskala (blau→cyan→grün→gelb→orange→rot) ------------------- */
function tempColor(t) { // t in [0,1] → {r,g,b} 0..1
  const stops = [
    [0.00, [0.10, 0.25, 0.75]], [0.25, [0.10, 0.70, 0.85]], [0.45, [0.20, 0.80, 0.40]],
    [0.65, [0.95, 0.90, 0.25]], [0.82, [0.98, 0.55, 0.15]], [1.00, [0.90, 0.16, 0.16]]];
  t = Math.min(Math.max(t, 0), 1);
  for (let i = 1; i < stops.length; i++) {
    if (t <= stops[i][0]) {
      const a = stops[i - 1], b = stops[i];
      const f = (t - a[0]) / (b[0] - a[0]);
      return { r: a[1][0] + (b[1][0] - a[1][0]) * f, g: a[1][1] + (b[1][1] - a[1][1]) * f, b: a[1][2] + (b[1][2] - a[1][2]) * f };
    }
  }
  const e = stops[stops.length - 1][1]; return { r: e[0], g: e[1], b: e[2] };
}
const colNorm = t => tempColor((t - Tmin) / Math.max(Tmax - Tmin, 1e-6));

/* =========================================================================
   (3) GEOMETRIE / SZENE
   ========================================================================= */
const sceneEl = document.getElementById('scene');
let renderer;
try {
  renderer = new THREE.WebGLRenderer({ antialias: true });
  renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
  renderer.setClearColor(0x0b0f14, 1);
  renderer.domElement.style.cssText = 'display:block;width:100%;height:100%';
  sceneEl.appendChild(renderer.domElement);
} catch (e) {
  showError('WebGL konnte nicht initialisiert werden. Bitte einen Browser mit ' +
    'aktivierter WebGL-Unterstützung verwenden.\n\n' + e.message);
}

const scene = new THREE.Scene();
scene.background = new THREE.Color(0x0b0f14);
scene.fog = new THREE.Fog(0x0b0f14, 700, 1800);  // weit entfernt: blendet Geometrie nicht aus

const camera = new THREE.PerspectiveCamera(50, 1, 0.5, 5000);
camera.position.set(-90, 150, 220);
camera.lookAt(L / 2, 0, W / 2);

// OrbitControls
let controls;
if (typeof OrbitControls === 'function') {
  controls = new OrbitControls(camera, renderer.domElement);
  controls.enableDamping = true;
  controls.dampingFactor = 0.08;
  controls.target.set(L / 2, 0, W / 2);
} else {
  showError('OrbitControls nicht geladen – einfache Maussteuerung aktiv.');
  controls = makeFallbackControls(camera, renderer.domElement, new THREE.Vector3(L / 2, 0, W / 2));
}

/* einfache Orbit-Steuerung als Ersatz, falls OrbitControls fehlt */
function makeFallbackControls(cam, dom, target) {
  const sph = new THREE.Spherical().setFromVector3(cam.position.clone().sub(target));
  let dragging = false, px = 0, py = 0;
  dom.addEventListener('mousedown', e => { dragging = true; px = e.clientX; py = e.clientY; });
  window.addEventListener('mouseup', () => dragging = false);
  window.addEventListener('mousemove', e => {
    if (!dragging) return;
    sph.theta -= (e.clientX - px) * 0.005; sph.phi -= (e.clientY - py) * 0.005;
    sph.phi = Math.max(0.05, Math.min(Math.PI - 0.05, sph.phi)); px = e.clientX; py = e.clientY;
  });
  dom.addEventListener('wheel', e => {
    e.preventDefault();
    sph.radius *= (1 + Math.sign(e.deltaY) * 0.1);
  }, { passive: false });
  return { target, update() { cam.position.copy(target).add(new THREE.Vector3().setFromSpherical(sph)); cam.lookAt(target); } };
}

scene.add(new THREE.AmbientLight(0xffffff, 0.65));
const sun = new THREE.DirectionalLight(0xffffff, 0.8); sun.position.set(120, 260, 90); scene.add(sun);

/* Koordinaten-Mapping: Three.X = Fluss-x, Three.Y = Höhe (z), Three.Z = Fluss-y */
const V = (x, zh, y) => new THREE.Vector3(x, zh, y);

/* leichte Unregelmäßigkeit der Sohle */
function noiseBed(x, y) {
  return Math.sin(x * 0.06) * 0.16 + Math.cos(y * 0.18) * 0.12
    + Math.sin(x * 0.013 + y * 0.07) * 0.20;
}

/* generisches parametrisches Flächen-Mesh über (x in 0..L, y in 0..W);
   fz(x,y) liefert die Höhe (Three.Y) */
function gridMesh(nx, ny, fz, material) {
  const pos = [], ind = [], row = ny + 1;
  for (let a = 0; a <= nx; a++) {
    for (let b = 0; b <= ny; b++) {
      const x = a / nx * L, y = b / ny * W;
      pos.push(x, fz(x, y), y);
    }
  }
  for (let a = 0; a < nx; a++) {
    for (let b = 0; b < ny; b++) {
      const i0 = a * row + b, i1 = (a + 1) * row + b, i2 = (a + 1) * row + b + 1, i3 = a * row + b + 1;
      ind.push(i0, i1, i2, i0, i2, i3);
    }
  }
  const g = new THREE.BufferGeometry();
  g.setAttribute('position', new THREE.Float32BufferAttribute(pos, 3));
  g.setIndex(ind); g.computeVertexNormals();
  return new THREE.Mesh(g, material);
}

/* Text-Sprite-Helfer (Maßstab/Beschriftung) */
function textSprite(txt, color) {
  const c = document.createElement('canvas'); c.width = 256; c.height = 64;
  const ctx = c.getContext('2d'); ctx.fillStyle = color || '#8aa0b5';
  ctx.font = 'bold 40px Segoe UI, sans-serif'; ctx.textAlign = 'center'; ctx.textBaseline = 'middle';
  ctx.fillText(txt, 128, 32);
  const tex = new THREE.CanvasTexture(c);
  const sp = new THREE.Sprite(new THREE.SpriteMaterial({ map: tex, transparent: true, depthTest: false }));
  sp.scale.set(26, 6.5, 1); return sp;
}

/* Umgebung: Ufer/Boden auf Höhe der Wasserlinie + parabolisches Flussbett + Maßstab */
const envGroup = new THREE.Group(); scene.add(envGroup);
let bedMesh = null;
function buildEnv() {
  envGroup.clear();
  const land = new THREE.MeshLambertMaterial({ color: 0x2f3a2a, side: THREE.DoubleSide });
  function bank(zc) {
    const m = new THREE.Mesh(new THREE.PlaneGeometry(L, 44), land);
    m.rotation.x = -Math.PI / 2; m.position.set(L / 2, WSURF - 0.02, zc); return m;
  }
  envGroup.add(bank(-22)); envGroup.add(bank(W + 22));
  // parabolisches Flussbett (sandig)
  bedMesh = gridMesh(80, 28, (x, y) => bedZ(y) + noiseBed(x, y),
    new THREE.MeshLambertMaterial({ color: 0x8a7250, flatShading: true, side: THREE.DoubleSide }));
  envGroup.add(bedMesh);

  // Beschriftungen
  [0, 100, 200, 300].forEach(x => { const s = textSprite(x + ' m'); s.position.set(x, WSURF + 2, -6); envGroup.add(s); });
  { const s = textSprite('Breite 50 m'); s.position.set(-16, WSURF + 2, W / 2); envGroup.add(s); }
  { const s = textSprite('Strömung →', '#3fc1c9'); s.position.set(L / 2, WSURF + 2, W + 9); envGroup.add(s); }
}

/* Wasservolumen (semitransparent blau) mit parabolischem Boden, flacher Oberfläche */
let waterMesh = null;
function buildWater() {
  if (waterMesh) { scene.remove(waterMesh); waterMesh.traverse(o => { if (o.geometry) o.geometry.dispose(); }); }
  waterMesh = new THREE.Group();
  const mat = new THREE.MeshPhongMaterial({
    color: 0x1e6fb0, transparent: true, opacity: 0.24,
    depthWrite: false, shininess: 60, side: THREE.DoubleSide
  });
  // flache Oberfläche
  waterMesh.add(gridMesh(40, 20, () => WSURF - 0.02, mat));
  // Stirn- und Seitenflächen (Band zwischen Sohle und Oberfläche)
  const seg = 28;
  function strip(low, high) {
    const pos = [], ind = [];
    for (let b = 0; b <= seg; b++) { const a = low(b), o = high(b); pos.push(a.x, a.y, a.z, o.x, o.y, o.z); }
    for (let b = 0; b < seg; b++) { const i = b * 2; ind.push(i, i + 1, i + 2, i + 1, i + 3, i + 2); }
    const g = new THREE.BufferGeometry();
    g.setAttribute('position', new THREE.Float32BufferAttribute(pos, 3));
    g.setIndex(ind); g.computeVertexNormals();
    return new THREE.Mesh(g, mat);
  }
  for (const xp of [0, L])  // Stirnseiten (zeigen den parabolischen Querschnitt)
    waterMesh.add(strip(b => { const y = b / seg * W; return { x: xp, y: bedZ(y), z: y }; },
      b => { const y = b / seg * W; return { x: xp, y: WSURF, z: y }; }));
  for (const yp of [0, W])  // Seitenwände
    waterMesh.add(strip(b => { const x = b / seg * L; return { x: x, y: bedZ(yp), z: yp }; },
      b => { const x = b / seg * L; return { x: x, y: WSURF, z: yp }; }));
  scene.add(waterMesh);
}

/* Pumpengebäude + geschlossener Rohrkreislauf (Fluss → Pumpe → Fluss) */
const pumpGroup = new THREE.Group(); scene.add(pumpGroup);
let capIn = null, capOut = null;
let circuitPath = [], circuitPumpIdx = 0, circuitCum = [], circuitTotal = 1, pipeRadius = 0.4;
function buildCircuitLengths() {
  circuitCum = [0];
  for (let i = 0; i < circuitPath.length - 1; i++)
    circuitCum.push(circuitCum[i] + circuitPath[i].distanceTo(circuitPath[i + 1]));
  circuitTotal = Math.max(circuitCum[circuitCum.length - 1], 1e-6);
}
/* Position bei Bogenlängen-Anteil s∈[0,1] entlang des Pfades */
function samplePath(s, out) {
  const d = s * circuitTotal; let i = 0;
  while (i < circuitCum.length - 2 && circuitCum[i + 1] < d) i++;
  const segLen = Math.max(circuitCum[i + 1] - circuitCum[i], 1e-6);
  const t = (d - circuitCum[i]) / segLen;
  out.lerpVectors(circuitPath[i], circuitPath[i + 1], t);
  return circuitCum[i] + t * segLen;   // zurück: zurückgelegte Länge (für Pumpen-Index)
}
// Tiefe (Bogenlänge) ab der das Wasser abgekühlt ist
function circuitCoolFrac() { return (circuitCum[circuitPumpIdx] || 0) / circuitTotal; }
// Zylinder zwischen zwei Punkten
function tube(p1, p2, r, mat) {
  const dir = new THREE.Vector3().subVectors(p2, p1); const len = dir.length();
  const m = new THREE.Mesh(new THREE.CylinderGeometry(r, r, len, 14), mat);
  m.position.copy(p1).add(p2).multiplyScalar(0.5);
  m.quaternion.setFromUnitVectors(new THREE.Vector3(0, 1, 0), dir.clone().normalize());
  return m;
}
function tubePath(pts, r, mat) {
  const g = new THREE.Group();
  for (let i = 0; i < pts.length - 1; i++) g.add(tube(pts[i], pts[i + 1], r, mat));
  for (let i = 1; i < pts.length - 1; i++) {
    const s = new THREE.Mesh(new THREE.SphereGeometry(r, 10, 10), mat);
    s.position.copy(pts[i]); g.add(s);
  }   // glatte Knoten an den Ecken
  return g;
}
function buildPump() {
  pumpGroup.clear();
  const xMid = (X_IN + X_OUT) / 2, by = WSURF, bz = -14;
  // Gebäude auf dem Ufer
  const b = new THREE.Mesh(new THREE.BoxGeometry(22, 10, 16), new THREE.MeshLambertMaterial({ color: 0x9aa7b3 }));
  b.position.set(xMid, by + 5, bz); pumpGroup.add(b);
  const roof = new THREE.Mesh(new THREE.BoxGeometry(23, 1, 17), new THREE.MeshLambertMaterial({ color: 0x46525e }));
  roof.position.set(xMid, by + 10.2, bz); pumpGroup.add(roof);
  const lbl = textSprite('Wärmepumpe', '#eef4fa'); lbl.scale.set(36, 9, 1);
  lbl.position.set(xMid, by + 13, bz); pumpGroup.add(lbl);
  // Fernwärme-Abgang
  const heat = new THREE.MeshPhongMaterial({ color: 0xd9572e, shininess: 30 });
  pumpGroup.add(tubePath([V(xMid, by + 10, bz), V(xMid, by + 16, bz - 6)], 0.5, heat));
  const hl = textSprite('→ Fernwärme', '#ffae8a'); hl.scale.set(30, 7.5, 1);
  hl.position.set(xMid, by + 17.5, bz - 9); pumpGroup.add(hl);

  const r = Math.max(P.dPipe / 2, 0.15); pipeRadius = r;
  const zp = pipeZ();                  // Mündungstiefe ≈ 0,5 m unter der Oberfläche
  // durchsichtige Rohre
  const matIn = new THREE.MeshPhongMaterial({ color: 0x8fc0e6, shininess: 50, transparent: true, opacity: 0.22, depthWrite: false, side: THREE.DoubleSide });
  const matOut = new THREE.MeshPhongMaterial({ color: 0x6f9bff, shininess: 50, transparent: true, opacity: 0.22, depthWrite: false, side: THREE.DoubleSide });
  // Entnahme
  const inPath = [
    V(X_IN, zp, Y_BANK),
    V(X_IN, by + 5, Y_BANK),
    V(X_IN, by + 5, bz),
    V(xMid - 11, by + 5, bz)
  ];
  // Rückgabe
  const outPath = [
    V(xMid + 11, by + 5, bz),
    V(X_OUT, by + 5, bz),
    V(X_OUT, by + 5, Y_BANK),
    V(X_OUT, zp, Y_BANK)
  ];
  pumpGroup.add(tubePath(inPath, r, matIn));
  pumpGroup.add(tubePath(outPath, r, matOut));
  // vollständiger Kreislaufpfad
  circuitPath = inPath.concat([V(xMid, by + 5, bz)], outPath);
  circuitPumpIdx = inPath.length;     // ab hier ist das Medium abgekühlt
  buildCircuitLengths();
  // Mündungs-Marker
  capIn = new THREE.Mesh(new THREE.SphereGeometry(r * 1.4, 12, 12), new THREE.MeshBasicMaterial({ color: 0xffffff }));
  capIn.position.copy(inPath[0]); pumpGroup.add(capIn);
  capOut = new THREE.Mesh(new THREE.SphereGeometry(r * 1.4, 12, 12), new THREE.MeshBasicMaterial({ color: 0xffffff }));
  capOut.position.copy(outPath[outPath.length - 1]); pumpGroup.add(capOut);
}

/* =========================================================================
   (4) VISUALISIERUNG: Heatmap-Schnittebenen, Partikel-Stromlinien, Vektoren
   ========================================================================= */
function makeHeatPlane(wRes, hRes) {
  const cv = document.createElement('canvas'); cv.width = wRes; cv.height = hRes;
  const ctx = cv.getContext('2d');
  const img = ctx.createImageData(wRes, hRes);
  const tex = new THREE.CanvasTexture(cv);
  tex.magFilter = THREE.LinearFilter; tex.minFilter = THREE.LinearFilter;
  const mat = new THREE.MeshBasicMaterial({ map: tex, transparent: true, opacity: 0.95, side: THREE.DoubleSide });
  const mesh = new THREE.Mesh(new THREE.PlaneGeometry(1, 1), mat);
  return { cv, ctx, img, tex, mesh };
}
// kontinuierlich abgetastete Ebenen
const planeCross = makeHeatPlane(90, 60);    // Querschnitt (y × z)
const planeHoriz = makeHeatPlane(180, 60);   // Horizontalschnitt (x × y)
const planeLong = makeHeatPlane(200, 60);    // Längsschnitt  (x × z)
scene.add(planeCross.mesh, planeHoriz.mesh, planeLong.mesh);

/* Heatmap-Pixel füllen */
function paint(plane, getT) {
  const { img } = plane, w = img.width, h = img.height, d = img.data;
  for (let b = 0; b < h; b++) {
    for (let a = 0; a < w; a++) {
      const t = getT(a, b); const o = (b * w + a) * 4;
      if (t == null || t !== t) { d[o + 3] = 0; continue; }   // unter der Sohle: transparent
      const c = colNorm(t); d[o] = c.r * 255; d[o + 1] = c.g * 255; d[o + 2] = c.b * 255; d[o + 3] = 235;
    }
  }
  plane.ctx.putImageData(img, 0, 0); plane.tex.needsUpdate = true;
}
function updatePlanes() {
  if (planeCross.mesh.visible) paint(planeCross, (a, b) => {
    const y = (1 - a / (planeCross.cv.width - 1)) * W, z = WSURF - (b / (planeCross.cv.height - 1)) * P.depth;
    return (z < bedZ(y)) ? null : sampleT(P.xCut, y, z);
  });
  if (planeHoriz.mesh.visible) {
    const zh = WSURF - P.zCut; paint(planeHoriz, (a, b) => {
      const x = a / (planeHoriz.cv.width - 1) * L, y = b / (planeHoriz.cv.height - 1) * W;
      return (zh < bedZ(y)) ? null : sampleT(x, y, zh);
    });
  }
  if (planeLong.mesh.visible) paint(planeLong, (a, b) => {
    const x = a / (planeLong.cv.width - 1) * L, z = WSURF - (b / (planeLong.cv.height - 1)) * P.depth;
    return (z < bedZ(P.yLong)) ? null : sampleT(x, P.yLong, z);
  });
}
/* Schnittebenen im Raum positionieren/orientieren */
function placePlanes() {
  planeCross.mesh.geometry.dispose();
  planeCross.mesh.geometry = new THREE.PlaneGeometry(W, P.depth);
  planeCross.mesh.rotation.set(0, Math.PI / 2, 0);
  planeCross.mesh.position.set(P.xCut, WSURF - P.depth / 2, W / 2);
  planeHoriz.mesh.geometry.dispose();
  planeHoriz.mesh.geometry = new THREE.PlaneGeometry(L, W);
  planeHoriz.mesh.rotation.set(-Math.PI / 2, 0, 0);
  planeHoriz.mesh.position.set(L / 2, WSURF - P.zCut, W / 2);
  planeLong.mesh.geometry.dispose();
  planeLong.mesh.geometry = new THREE.PlaneGeometry(L, P.depth);
  planeLong.mesh.rotation.set(0, 0, 0);
  planeLong.mesh.position.set(L / 2, WSURF - P.depth / 2, P.yLong);
}

/* runde Partikel */
function discTexture() {
  const c = document.createElement('canvas'); c.width = c.height = 64;
  const g = c.getContext('2d');
  const grd = g.createRadialGradient(32, 32, 2, 32, 32, 30);
  grd.addColorStop(0, 'rgba(255,255,255,1)');
  grd.addColorStop(0.55, 'rgba(255,255,255,0.95)');
  grd.addColorStop(1, 'rgba(255,255,255,0)');
  g.fillStyle = grd; g.beginPath(); g.arc(32, 32, 30, 0, Math.PI * 2); g.fill();
  const t = new THREE.CanvasTexture(c); return t;
}
const DISC = discTexture();

/* Particle Tracer */
const NP = 2600;
const pPos = new Float32Array(NP * 3), pCol = new Float32Array(NP * 3);
const partGeo = new THREE.BufferGeometry();
partGeo.setAttribute('position', new THREE.BufferAttribute(pPos, 3));
partGeo.setAttribute('color', new THREE.BufferAttribute(pCol, 3));
const partMat = new THREE.PointsMaterial({ size: 1.7, map: DISC, alphaTest: 0.4, vertexColors: true, transparent: true, opacity: 0.95, sizeAttenuation: true, depthWrite: false });
const particles = new THREE.Points(partGeo, partMat);
particles.frustumCulled = false;
scene.add(particles);
const pX = new Float32Array(NP), pY = new Float32Array(NP), pZ = new Float32Array(NP);
function seedParticle(n, atInlet) {
  pX[n] = atInlet ? Math.random() * 3 : Math.random() * L;
  const y = 0.5 + Math.random() * (W - 1); pY[n] = y;
  pZ[n] = bedZ(y) + (0.08 + 0.86 * Math.random()) * localDepth(y);
}
for (let n = 0; n < NP; n++) seedParticle(n, false);
function updateParticles(dt) {
  const zp = pipeZ();
  for (let n = 0; n < NP; n++) {
    velocityAt(pX[n], pY[n], pZ[n], _vel);
    const sp = Math.hypot(_vel.u, _vel.v, _vel.w);
    pX[n] += _vel.u * dt; pY[n] += _vel.v * dt; pZ[n] += _vel.w * dt;
    const di = pX[n] - X_IN, dj = pY[n] - Y_BANK, dk = pZ[n] - zp;
    const sucked = (di * di + dj * dj + dk * dk) < 2.5;
    if (!(pX[n] < L) || pX[n] < 0 || pY[n] < 0 || pY[n] > W || pZ[n] < bedZ(pY[n]) || pZ[n] > WSURF || sp < 0.015 || sucked || !isFinite(pX[n]))
      seedParticle(n, true);
    const o = n * 3; pPos[o] = pX[n]; pPos[o + 1] = pZ[n]; pPos[o + 2] = pY[n];
    const c = colNorm(sampleT(pX[n], pY[n], pZ[n]));
    pCol[o] = c.r; pCol[o + 1] = c.g; pCol[o + 2] = c.b;
  }
  partGeo.attributes.position.needsUpdate = true;
  partGeo.attributes.color.needsUpdate = true;
}

/* Kreislauf-Partikel */
const NCP = 900;
const cpS = new Float32Array(NCP);
const cpOff = new Float32Array(NCP * 3);
const cpPos = new Float32Array(NCP * 3), cpCol = new Float32Array(NCP * 3);
for (let n = 0; n < NCP; n++) {
  cpS[n] = n / NCP;
  let ox, oy, oz, d2;
  do { ox = Math.random() * 2 - 1; oy = Math.random() * 2 - 1; oz = Math.random() * 2 - 1; d2 = ox * ox + oy * oy + oz * oz; } while (d2 > 1);
  cpOff[n * 3] = ox; cpOff[n * 3 + 1] = oy; cpOff[n * 3 + 2] = oz;
}
const circGeo = new THREE.BufferGeometry();
circGeo.setAttribute('position', new THREE.BufferAttribute(cpPos, 3));
circGeo.setAttribute('color', new THREE.BufferAttribute(cpCol, 3));
const circParticles = new THREE.Points(circGeo,
  new THREE.PointsMaterial({ size: 2.25, map: DISC, alphaTest: 0.35, vertexColors: true, transparent: true, opacity: 0.98, sizeAttenuation: true, depthWrite: false }));
circParticles.frustumCulled = false; scene.add(circParticles);
const _cp = new THREE.Vector3();
function updateCircuit(dt) {
  if (!circuitPath.length) return;
  const coolFrac = circuitCoolFrac();
  const cWarm = colNorm(tInletLocal), cCold = colNorm(tOutlet);
  const adv = dt * 0.05;
  const off = pipeRadius * 0.7;
  for (let n = 0; n < NCP; n++) {
    cpS[n] += adv; if (cpS[n] >= 1) cpS[n] -= 1;
    samplePath(cpS[n], _cp);
    const o = n * 3;
    cpPos[o] = _cp.x + cpOff[o] * off;
    cpPos[o + 1] = _cp.y + cpOff[o + 1] * off;
    cpPos[o + 2] = _cp.z + cpOff[o + 2] * off;
    const c = (cpS[n] < coolFrac) ? cWarm : cCold;
    cpCol[o] = c.r; cpCol[o + 1] = c.g; cpCol[o + 2] = c.b;
  }
  circGeo.attributes.position.needsUpdate = true;
  circGeo.attributes.color.needsUpdate = true;
}

/* Auslassstrahl */
const NJ = 1500;
const jX = new Float32Array(NJ), jY = new Float32Array(NJ), jZ = new Float32Array(NJ);
const jU = new Float32Array(NJ), jV = new Float32Array(NJ), jW = new Float32Array(NJ), jAge = new Float32Array(NJ);
const jPos = new Float32Array(NJ * 3), jCol = new Float32Array(NJ * 3);
const jetGeo = new THREE.BufferGeometry();
jetGeo.setAttribute('position', new THREE.BufferAttribute(jPos, 3));
jetGeo.setAttribute('color', new THREE.BufferAttribute(jCol, 3));
const jetParticles = new THREE.Points(jetGeo,
  new THREE.PointsMaterial({ size: 1.7, map: DISC, alphaTest: 0.4, vertexColors: true, transparent: true, opacity: 0.96, sizeAttenuation: true, depthWrite: false }));
jetParticles.frustumCulled = false; scene.add(jetParticles);
const JET_LIFE = 30;
function seedJet(n, fresh) {
  const zp = pipeZ();
  jX[n] = X_OUT + (Math.random() - 0.5) * 0.6;
  jY[n] = Y_BANK + (Math.random() - 0.5) * 0.6;
  jZ[n] = zp + (Math.random() - 0.5) * 0.6;
  const spd = P.vSuct * (0.45 + 0.55 * Math.random());
  const aXY = (Math.random() - 0.5) * 1.6;
  const aZ = (Math.random() - 0.5) * 1.1;
  const bx = 0.45, by = 0.90;
  const cs = Math.cos(aXY), sn = Math.sin(aXY);
  let dx = bx * cs - by * sn, dy = bx * sn + by * cs, dz = Math.tan(aZ) * 0.6;
  const inv = spd / Math.max(Math.hypot(dx, dy, dz), 1e-6);
  jU[n] = dx * inv; jV[n] = dy * inv; jW[n] = dz * inv;
  jAge[n] = fresh ? 0 : Math.random() * JET_LIFE;
}
for (let n = 0; n < NJ; n++) seedJet(n, false);
function updateJet(dt) {
  for (let n = 0; n < NJ; n++) {
    velocityAt(jX[n], jY[n], jZ[n], _vel);
    const decay = Math.exp(-jAge[n] / 8);
    jX[n] += (_vel.u + jU[n] * decay) * dt;
    jY[n] += (_vel.v + jV[n] * decay) * dt;
    jZ[n] += (_vel.w + jW[n] * decay) * dt;
    jAge[n] += dt * 0.5;
    if (jAge[n] > JET_LIFE || !(jX[n] < L) || jY[n] < 0 || jY[n] > W || jZ[n] < bedZ(jY[n]) || jZ[n] > WSURF || !isFinite(jX[n]))
      seedJet(n, true);
    const o = n * 3; jPos[o] = jX[n]; jPos[o + 1] = jZ[n]; jPos[o + 2] = jY[n];
    const fr = Math.min(jAge[n] / 14, 1);
    const tval = tOutlet * (1 - fr) + sampleT(jX[n], jY[n], jZ[n]) * fr;
    const c = colNorm(tval); jCol[o] = c.r; jCol[o + 1] = c.g; jCol[o + 2] = c.b;
  }
  jetGeo.attributes.position.needsUpdate = true;
  jetGeo.attributes.color.needsUpdate = true;
}

/* Geschwindigkeitsvektoren */
let vecGroup = new THREE.Group(); scene.add(vecGroup);
function buildVectors() {
  vecGroup.clear();
  const mat = new THREE.LineBasicMaterial({ color: 0x6fb7c9, transparent: true, opacity: 0.6 });
  const pts = [];
  for (let i = 2; i < NX; i += 8) {
    for (let j = 1; j < NY; j += 4) {
      const x = (i + 0.5) * dx, y = (j + 0.5) * dy, z = bedZ(y) + 0.6 * localDepth(y);
      velocityAt(x, y, z, _vel);
      const s = 6 / Math.max(1, Math.hypot(_vel.u, _vel.v, _vel.w));
      pts.push(x, z, y, x + _vel.u * s, z + _vel.w * s, y + _vel.v * s);
    }
  }
  const g = new THREE.BufferGeometry();
  g.setAttribute('position', new THREE.Float32BufferAttribute(pts, 3));
  vecGroup.add(new THREE.LineSegments(g, mat));
}

/* ---- Konturlinien (Marching Squares) ---- */
const MS_EDGES = {
  1: [[3, 0]], 2: [[0, 1]], 3: [[3, 1]], 4: [[1, 2]], 5: [[3, 0], [1, 2]], 6: [[0, 2]],
  7: [[3, 2]], 8: [[2, 3]], 9: [[2, 0]], 10: [[0, 1], [2, 3]], 11: [[2, 1]], 12: [[1, 3]], 13: [[1, 0]], 14: [[0, 3]]
};
function marchingSquares(val, gx, gy, level) {
  const seg = [];
  const ep = (edge, x, y, vtl, vtr, vbr, vbl) => {
    let t;
    if (edge === 0) { t = (level - vbl) / (vbr - vbl); return [x + (isFinite(t) ? t : 0.5), y + 1]; }
    if (edge === 1) { t = (level - vbr) / (vtr - vbr); return [x + 1, y + 1 - (isFinite(t) ? t : 0.5)]; }
    if (edge === 2) { t = (level - vtr) / (vtl - vtr); return [x + 1 - (isFinite(t) ? t : 0.5), y]; }
    t = (level - vtl) / (vbl - vtl); return [x, y + (isFinite(t) ? t : 0.5)];
  };
  for (let y = 0; y < gy - 1; y++) {
    for (let x = 0; x < gx - 1; x++) {
      const vtl = val(x, y), vtr = val(x + 1, y), vbr = val(x + 1, y + 1), vbl = val(x, y + 1);
      let ci = 0; if (vbl > level) ci |= 1; if (vbr > level) ci |= 2; if (vtr > level) ci |= 4; if (vtl > level) ci |= 8;
      const pairs = MS_EDGES[ci]; if (!pairs) continue;
      for (const pr of pairs) {
        const p = ep(pr[0], x, y, vtl, vtr, vbr, vbl), q = ep(pr[1], x, y, vtl, vtr, vbr, vbl);
        seg.push(p[0], p[1], q[0], q[1]);
      }
    }
  }
  return seg;
}
function contourLines(getScalar, levels, colorFn, zOff) {
  const GX = 90, GY = 34;
  const cache = []; for (let b = 0; b < GY; b++) {
    cache[b] = []; for (let a = 0; a < GX; a++)
      cache[b][a] = getScalar(a / (GX - 1) * L, b / (GY - 1) * W);
  }
  const val = (a, b) => cache[b][a];
  const group = new THREE.Group();
  for (const lv of levels) {
    const seg = marchingSquares(val, GX, GY, lv); if (!seg.length) continue;
    const pos = [];
    for (let s = 0; s < seg.length; s += 4) {
      pos.push(seg[s] / (GX - 1) * L, zOff, seg[s + 1] / (GY - 1) * W,
        seg[s + 2] / (GX - 1) * L, zOff, seg[s + 3] / (GY - 1) * W);
    }
    const g = new THREE.BufferGeometry();
    g.setAttribute('position', new THREE.Float32BufferAttribute(pos, 3));
    const col = colorFn(lv);
    group.add(new THREE.LineSegments(g, new THREE.LineBasicMaterial({ color: col, transparent: true, opacity: 0.95 })));
  }
  group.frustumCulled = false;
  return group;
}

/* Isothermen */
let isoGroup = null;
function buildIso() {
  if (isoGroup) { scene.remove(isoGroup); isoGroup.traverse(o => o.geometry && o.geometry.dispose()); }
  const lo = Tmin, hi = Tmax, levels = [];
  for (let n = 1; n <= 7; n++) levels.push(lo + (hi - lo) * n / 8);
  const cN = new THREE.Color();
  isoGroup = contourLines((x, y) => sampleT(x, y, WSURF - 0.05), levels,
    lv => { const c = colNorm(lv); cN.setRGB(c.r, c.g, c.b); return cN.getHex(); }, WSURF + 0.06);
  isoGroup.visible = document.getElementById('t_iso').checked;
  scene.add(isoGroup);
}

/* Strömungs-Isolinien */
let flowGroup = null;
function buildFlow() {
  if (flowGroup) { scene.remove(flowGroup); flowGroup.traverse(o => o.geometry && o.geometry.dispose()); }
  const speed = (x, y) => {
    const z = bedZ(y) + 0.85 * localDepth(y); velocityAt(x, y, z, _vel);
    return Math.hypot(_vel.u, _vel.v, _vel.w);
  };
  let vmax = 0.1; for (let a = 0; a < 12; a++) for (let b = 0; b < 8; b++) vmax = Math.max(vmax, speed(a / 11 * L, b / 7 * W));
  const levels = []; for (let n = 1; n <= 6; n++) levels.push(vmax * n / 7);
  flowGroup = contourLines(speed, levels, () => 0x8fe3ff, WSURF + 0.10);
  flowGroup.visible = document.getElementById('t_flow').checked;
  scene.add(flowGroup);
}

/* =========================================================================
   (5) UI – Regler, abgeleitete Werte, Schalter, KPIs
   ========================================================================= */
const fmt = (v, d = 1) => v.toLocaleString('de-DE', { minimumFractionDigits: d, maximumFractionDigits: d });
const UNIT = {
  vFlow: ' m/s', depth: ' m', tIn: ' °C', power: ' kW', vSuct: ' m/s',
  dPipe: ' m', xCut: ' m', zCut: ' m', yLong: ' m'
};
const DEC = { vFlow: 2, depth: 1, tIn: 1, power: 0, vSuct: 2, dPipe: 2, xCut: 0, zCut: 1, yLong: 1 };
const powMW = () => fmt(P.power / 1000, 2) + ' MW';

function refreshLabels() {
  document.querySelectorAll('.ctrl[data-p]').forEach(c => {
    const p = c.dataset.p;
    c.querySelector('.val').textContent = (p === 'power') ? powMW() : fmt(P[p], DEC[p]) + UNIT[p];
    c.querySelector('input').value = P[p];
  });
}
function refreshDerived() {
  document.getElementById('derivedFlow').innerHTML =
    `Reynolds-Zahl Re ≈ <b>${D.Re.toExponential(2)}</b> &nbsp;·&nbsp; turbulent<br>` +
    `Péclet-Zahl Pe ≈ <b>${D.Pe.toExponential(2)}</b><br>` +
    `eff. Temp.-Leitf. α ≈ <b>${D.alpha.toExponential(2)}</b> m²/s`;
  document.getElementById('derivedPump').innerHTML =
    `Rohrquerschnitt A = <b>${fmt(D.area, 3)}</b> m²<br>` +
    `Durchfluss Q = <b>${fmt(D.Qpump, 3)}</b> m³/s = <b>${fmt(D.Qpump * 1000, 0)}</b> l/s<br>` +
    `Temp.-Absenkung ΔT = <b>${fmt(D.dT, 2)}</b> K<br>` +
    `Auslasstemperatur ≈ <b>${fmt(tOutlet, 2)}</b> °C`;
}

/* Slider-Verdrahtung */
document.querySelectorAll('.ctrl[data-p] input').forEach(inp => {
  const p = inp.closest('.ctrl').dataset.p;
  inp.addEventListener('input', () => {
    P[p] = parseFloat(inp.value);
    onParamChange(p);
  });
});
function onParamChange(p) {
  computeDerived();
  if (p === 'depth') { buildEnv(); buildWater(); buildPump(); }
  if (p === 'dPipe') { buildPump(); }
  placePlanes();
  refreshLabels(); refreshDerived();
  relax(16);
  fieldStats(); buildVectors(); buildFlow();
}

/* Schalter (Toggles) */
const TG = {
  t_cross: planeCross.mesh, t_long: planeLong.mesh, t_horiz: planeHoriz.mesh,
  t_stream: particles, t_vec: vecGroup
};
Object.entries(TG).forEach(([id, obj]) => {
  const el = document.getElementById(id);
  obj.visible = el.checked;
  el.addEventListener('change', () => { obj.visible = el.checked; });
});
document.getElementById('t_stream').addEventListener('change', e => { circParticles.visible = e.target.checked; jetParticles.visible = e.target.checked; });
circParticles.visible = document.getElementById('t_stream').checked;
jetParticles.visible = document.getElementById('t_stream').checked;
document.getElementById('t_water').addEventListener('change', e => { if (waterMesh) waterMesh.visible = e.target.checked; });
document.getElementById('t_iso').addEventListener('change', e => { if (!isoGroup) buildIso(); isoGroup.visible = e.target.checked; });
document.getElementById('t_flow').addEventListener('change', e => { if (!flowGroup) buildFlow(); flowGroup.visible = e.target.checked; });

/* Ansichten */
function setView(v) {
  const cx = L / 2, cz = W / 2;
  const fov = camera.fov * Math.PI / 180;
  const fit = (L * 0.6) / Math.tan(fov / 2) / Math.max(camera.aspect, 1) + 120;
  const cy = WSURF * 0.5;
  if (v === 'top') camera.position.set(cx, fit * 0.9, cz + 0.1);
  else if (v === 'side') camera.position.set(cx, 30, cz + fit);
  else camera.position.set(cx - fit * 0.55, fit * 0.5, cz + fit * 0.7);
  controls.target.set(cx, cy, cz);
  camera.lookAt(cx, cy, cz);
  controls.update();
}
document.querySelectorAll('.views button').forEach(b =>
  b.addEventListener('click', () => setView(b.dataset.view)));

/* KPI-Panel */
function updateKPI() {
  const ref = P.tIn;
  let maxDrop = 0, recX = X_OUT;
  const jL = Math.round(Y_BANK / dy), ksurf = NZ - 1;
  for (let i = 0; i < NX; i++) {
    const t = T[idx(i, jL, ksurf)];
    const drop = ref - t; if (drop > maxDrop) maxDrop = drop;
  }
  const iO = Math.round(X_OUT / dx);
  recX = L;
  for (let i = iO; i < NX; i++) {
    const t = T[idx(i, jL, ksurf)];
    if (ref - t < 0.05) { recX = (i + 0.5) * dx; break; }
  }
  const recLen = Math.max(recX - X_OUT, 0);
  const Ariver = (2 / 3) * W * P.depth;
  const Qriver = Ariver * P.vFlow;
  const dTmix = (P.power * 1000) / (RHO * CP * Math.max(Qriver, 1e-6));
  const rows = [
    ['Entnommene Leistung', powMW()],
    ['Durchfluss Pumpe', fmt(D.Qpump, 3) + ' m³/s'],
    ['Durchfluss Fluss', fmt(Qriver, 1) + ' m³/s'],
    ['ΔT Rückgabe (Rohr)', fmt(D.dT, 2) + ' K'],
    ['Auslasstemperatur', fmt(tOutlet, 2) + ' °C'],
    ['Abkühlung stromab (durchmischt)', fmt(dTmix, 3) + ' K → ' + fmt(P.tIn - dTmix, 2) + ' °C'],
    ['Max. lokale Absenkung (Fahne)', fmt(maxDrop, 2) + ' K'],
    ['Reichweite Kaltwasserfahne', recLen >= L - X_OUT ? '> ' + fmt(L - X_OUT, 0) + ' m' : fmt(recLen, 0) + ' m'],
    ['Reynolds-Zahl', D.Re.toExponential(2)],
    ['Péclet-Zahl', D.Pe.toExponential(2)],
  ];
  document.getElementById('kpiBody').innerHTML = rows.map(r =>
    `<tr><td class="k">${r[0]}</td><td class="v">${r[1]}</td></tr>`).join('');
  const w = document.getElementById('warn'); let msg = '';
  if (D.dT > 10) msg = '⚠ ΔT sehr groß (' + fmt(D.dT, 1) + ' K): Durchfluss zu klein oder Leistung zu hoch.';
  if (tInletLocal - D.dT < 0) msg = '⚠ Rechnerische Auslasstemperatur < 0 °C – auf 0 °C begrenzt (Eisgefahr).';
  w.className = msg ? 'show' : ''; w.textContent = msg;
}

/* Colorbar zeichnen */
const cbar = document.getElementById('cbar'), cbx = cbar.getContext('2d');
function drawColorbar() {
  for (let y = 0; y < cbar.height; y++) {
    const t = 1 - y / cbar.height; const c = tempColor(t);
    cbx.fillStyle = `rgb(${c.r * 255 | 0},${c.g * 255 | 0},${c.b * 255 | 0})`;
    cbx.fillRect(0, y, cbar.width, 1);
  }
  const ticks = document.getElementById('cbarTicks'); ticks.innerHTML = '';
  for (let i = 0; i < 5; i++) {
    const t = Tmax - (Tmax - Tmin) * i / 4;
    const s = document.createElement('span'); s.textContent = fmt(t, 1); ticks.appendChild(s);
  }
}

/* =========================================================================
   (6) DIAGRAMME
   ========================================================================= */
function lineChart(canvas, series, xlabel, ylabel, xRange, yRange) {
  const ctx = canvas.getContext('2d');
  const W0 = canvas.clientWidth, H0 = canvas.clientHeight;
  if (canvas.width !== W0) canvas.width = W0; if (canvas.height !== H0) canvas.height = H0;
  ctx.clearRect(0, 0, W0, H0);
  const m = { l: 38, r: 8, t: 6, b: 20 };
  const pw = W0 - m.l - m.r, ph = H0 - m.t - m.b;
  ctx.strokeStyle = '#2a3a49'; ctx.lineWidth = 1;
  ctx.strokeRect(m.l, m.t, pw, ph);
  ctx.fillStyle = '#6b8096'; ctx.font = '9px Segoe UI'; ctx.textAlign = 'center';
  for (let g = 0; g <= 4; g++) {
    const yv = yRange[0] + (yRange[1] - yRange[0]) * g / 4;
    const py = m.t + ph - ph * g / 4;
    ctx.strokeStyle = '#1b2733'; ctx.beginPath(); ctx.moveTo(m.l, py); ctx.lineTo(m.l + pw, py); ctx.stroke();
    ctx.textAlign = 'right'; ctx.fillText(yv.toFixed(yRange[1] - yRange[0] < 3 ? 2 : 1), m.l - 4, py + 3);
  }
  ctx.textAlign = 'center';
  for (let g = 0; g <= 4; g++) {
    const xv = xRange[0] + (xRange[1] - xRange[0]) * g / 4;
    const px = m.l + pw * g / 4; ctx.fillText(xv.toFixed(0), px, H0 - 7);
  }
  ctx.fillStyle = '#52677a'; ctx.fillText(xlabel, m.l + pw / 2, H0 - 0.5);
  const sx = v => m.l + pw * (v - xRange[0]) / (xRange[1] - xRange[0]);
  const sy = v => m.t + ph - ph * (v - yRange[0]) / (yRange[1] - yRange[0]);
  series.forEach(s => {
    ctx.strokeStyle = s.color; ctx.lineWidth = 1.6; ctx.beginPath();
    s.data.forEach((d, n) => { const X = sx(d[0]), Y = sy(d[1]); n ? ctx.lineTo(X, Y) : ctx.moveTo(X, Y); });
    ctx.stroke();
  });
  if (series.length > 1) {
    ctx.textAlign = 'left'; let ly = m.t + 10;
    series.forEach(s => { ctx.fillStyle = s.color; ctx.fillText('— ' + s.name, m.l + 6, ly); ly += 11; });
  }
}
const chTx = document.getElementById('chTx'), chTy = document.getElementById('chTy'),
  chVz = document.getElementById('chVz'), chTz = document.getElementById('chTz');
function updateCharts() {
  const jMid = Math.round(NY / 2), jBank = Math.round(Y_BANK / dy), ksurf = NZ - 1;
  const iC = Math.min(Math.max(Math.round(P.xCut / dx), 0), NX - 1);
  const mid = [], bank = [];
  for (let i = 0; i < NX; i++) {
    const x = (i + 0.5) * dx;
    mid.push([x, T[idx(i, jMid, ksurf)]]); bank.push([x, T[idx(i, jBank, ksurf)]]);
  }
  const TLO = 0, THI = 30;
  lineChart(chTx, [{ name: 'Mittellinie', color: '#3fc1c9', data: mid },
  { name: 'Pumpenufer', color: '#ffb454', data: bank }],
    'x [m]', 'T', [0, L], [TLO, THI]);
  const ty = []; for (let j = 0; j < NY; j++) { const y = (j + 0.5) * dy; ty.push([y, T[idx(iC, j, ksurf)]]); }
  lineChart(chTy, [{ name: 'T(y)', color: '#54d98c', data: ty }], 'y [m]', 'T', [0, W], [TLO, THI]);
  const vz = []; let vmax = 0.1; const dlC = localDepth(W / 2), bC = bedZ(W / 2);
  for (let k = 0; k < NZ; k++) {
    const z = bC + (k + 0.5) / NZ * dlC; velocityAt(P.xCut, W / 2, z, _vel);
    const sp = Math.hypot(_vel.u, _vel.v, _vel.w); vmax = Math.max(vmax, sp); vz.push([z, sp]);
  }
  lineChart(chVz, [{ name: 'v(z)', color: '#9be7ff', data: vz }], 'Höhe z [m]', 'v', [bC, WSURF], [0, vmax * 1.1]);
  const tz = []; const dlB = localDepth((jBank + 0.5) * dy), bB = bedZ((jBank + 0.5) * dy);
  for (let k = 0; k < NZ; k++) { const z = bB + (k + 0.5) / NZ * dlB; tz.push([z, T[idx(iC, jBank, k)]]); }
  lineChart(chTz, [{ name: 'T(z)', color: '#ff8fab', data: tz }], 'Höhe z [m]', 'T', [bB, WSURF], [TLO, THI]);
}

/* =========================================================================
   (7) INITIALISIERUNG & ANIMATIONSSCHLEIFE
   ========================================================================= */
function resize() {
  let w = sceneEl.clientWidth, h = sceneEl.clientHeight;
  if (!w || !h) { w = window.innerWidth - 350; h = window.innerHeight - 280; }
  w = Math.max(w, 200); h = Math.max(h, 200);
  renderer.setSize(w, h, false);
  renderer.domElement.style.width = '100%';
  renderer.domElement.style.height = '100%';
  camera.aspect = w / h; camera.updateProjectionMatrix();
}
window.addEventListener('resize', resize);
if (window.ResizeObserver) { try { new ResizeObserver(() => resize()).observe(sceneEl); } catch (e) { } }

function init() {
  try {
    if (!renderer) return;
    computeDerived();
    T.fill(P.tIn);
    buildEnv(); buildWater(); buildPump(); placePlanes(); buildVectors();
    refreshLabels(); refreshDerived();
    relax(80);
    fieldStats(); buildIso(); buildFlow();
    resize(); setView('oblique'); drawColorbar();
    Object.entries(TG).forEach(([id, obj]) => obj.visible = document.getElementById(id).checked);
    circParticles.visible = document.getElementById('t_stream').checked;
    jetParticles.visible = document.getElementById('t_stream').checked;
    if (isoGroup) isoGroup.visible = document.getElementById('t_iso').checked;
    if (flowGroup) flowGroup.visible = document.getElementById('t_flow').checked;
    updatePlanes();
    controls.update();
    renderer.render(scene, camera);
    animate();
  } catch (e) { showError('init(): ' + (e && e.stack ? e.stack : e)); }
}

let frame = 0, animErr = false;
function animate() {
  requestAnimationFrame(animate);
  try {
    relax(2);
    fieldStats();
    updatePlanes();
    updateParticles(0.4);
    updateCircuit(0.4);
    updateJet(0.4);
    if (capIn) { const c = colNorm(tInletLocal); capIn.material.color.setRGB(c.r, c.g, c.b); }
    if (capOut) { const c = colNorm(tOutlet); capOut.material.color.setRGB(c.r, c.g, c.b); }
    if (frame % 6 === 0) {
      drawColorbar(); updateKPI(); updateCharts();
      const cv = renderer.domElement, hud = document.getElementById('hud');
      if (hud) hud.textContent = `Canvas ${cv.clientWidth}×${cv.clientHeight}px · Buffer ${cv.width}×${cv.height} · ` +
        `Objekte ${scene.children.length} · Kamera (${camera.position.x | 0},${camera.position.y | 0},${camera.position.z | 0})`;
    }
    if (frame % 24 === 0) {
      if (document.getElementById('t_iso').checked) buildIso();
      if (document.getElementById('t_flow').checked) buildFlow();
    }
    controls.update();
    renderer.render(scene, camera);
    frame++;
  } catch (e) { if (!animErr) { animErr = true; showError('animate(): ' + (e && e.stack ? e.stack : e)); } }
}

window.addEventListener('load', init);
// Auch sofort initialisieren falls das DOM schon fertig ist
if (document.readyState === 'complete' || document.readyState === 'interactive') {
  init();
}
