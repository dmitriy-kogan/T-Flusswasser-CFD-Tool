import * as THREE from 'three';
import { OrbitControls } from 'three/examples/jsm/controls/OrbitControls.js';

"use strict";
/* ---------- Sichtbare Fehleranzeige (statt stiller schwarzer Szene) ------ */
function showError(msg){
  const e=document.getElementById('err'); if(!e)return;
  e.className='show'; e.textContent='⚠ Fehler:\n'+msg;
}
window.addEventListener('error',ev=>{
  showError((ev.message||'')+'\n'+(ev.filename||'')+':'+(ev.lineno||'')+
    (ev.error&&ev.error.stack?'\n\n'+ev.error.stack:''));
});
window.addEventListener('unhandledrejection',ev=>{
  showError('Promise: '+(ev.reason&&ev.reason.message?ev.reason.message:ev.reason));
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
  tIn:   8.0,   // Temperatur des gesamten Flusswassers [°C]
  power: 500,   // thermische HEIZLEISTUNG P_heiz (Nutzwärme Fernwärme) [kW]
  cop:   4.0,   // Leistungszahl (COP) der Wärmepumpe [-]
  vSuct: 1.5,   // Ansauggeschwindigkeit im Rohr [m/s]
  dPipe: 0.8,   // Rohrdurchmesser [m]
  xCut:  220,   // Querschnittsposition [m]
  zCut:  1.5,   // Horizontalschnitt-Tiefe [m unter Oberfläche]
  yLong: 6      // Längsschnitt y-Position (quer zur Strömung) [m]
};

/* Geometrie des Flussabschnitts */
const L = 400, W = 50;                 // Länge (x), Breite (y) [m]  (100 m Vorlauf vor dem Einlass)
const X_IN = 120, X_OUT = 220;         // Position Einlass / Auslass [m]
const Y_BANK = 6;                      // Abstand der Rohrmündung vom Pumpenufer (y=0) [m]
const RHO = 1000, CP = 4186;           // Wasser: Dichte, spez. Wärmekapazität
const NU = 1.3e-6;                     // kinematische Viskosität [m²/s]
const ALPHA_MOL = 1.4e-7;              // molekulare Temperaturleitfähigkeit [m²/s]
const MINDEPTH = 0.4;                  // Resttiefe an den Ufern [m]
const WSURF = 10;                      // feste Höhe der Wasseroberfläche/Uferlinie [m]

/* Parabolisches Flussprofil: Oberfläche flach und FEST bei z = WSURF.
   Die Sohle ist parabelförmig (tiefste Stelle in der Mitte); beim Erhöhen der
   Tiefe wandert der tiefste Punkt nach unten (die Oberfläche bleibt fest). */
function localDepth(y){ const s=2*y/W-1; return Math.max(P.depth*(1-s*s), MINDEPTH); }
function bedZ(y){ return WSURF - localDepth(y); }          // Sohle wandert mit P.depth nach unten
const ZSURF = () => WSURF;                                  // Wasseroberfläche (flach, fest)
/* Rohrmündungstiefe: ca. 0,5 m unter der Wasseroberfläche (in flachen Zonen begrenzt) */
function pipeZ(){ return Math.max(WSURF - 0.5, bedZ(Y_BANK) + 0.2); }

/* Rechengitter (gröber für Echtzeit) */
/* Basisauflösung; über den Regler "Gitterauflösung" skalierbar (4 Stufen).
   dt ~ dx hält die Advektion stabil; die Mischung wird über PHYSIKALISCHE
   Diffusivitäten definiert (siehe solveStep) und ist damit gitterunabhängig. */
const NX0 = 160, NY0 = 26, NZ0 = 8;
const GRID_S = [0.75, 1, 1.5, 2];                       // Skalierungsfaktoren je Stufe
let NX = NX0, NY = NY0, NZ = NZ0;
let dx = L / NX, dy = W / NY;
let NC = NX * NY * NZ;
/* Zellvolumen je Quer-Spalte j (Sigma: alle NZ Schichten einer Säule gleich dick);
   wird bei Änderung von Wassertiefe oder Gitter neu berechnet */
let VCOL = new Float64Array(NY);
function updateVCOL(){
  if(VCOL.length!==NY) VCOL=new Float64Array(NY);
  for(let j=0;j<NY;j++){ const y=(j+0.5)*dy;
    VCOL[j] = dx*dy*Math.max(localDepth(y),MINDEPTH)/NZ; }
}
updateVCOL();

/* Temperaturfelder (Doppelpuffer) */
let T  = new Float32Array(NC);
let T2 = new Float32Array(NC);

const idx = (i,j,k) => i + NX*(j + NY*k);

/* ---------- (2) SOLVER --------------------------------------------------- */
/* abgeleitete Größen, bei jeder Parameteränderung aktualisiert */
const D = {};                          // derived values
function computeDerived(){
  D.area   = Math.PI * Math.pow(P.dPipe/2, 2);     // Rohrquerschnitt [m²]
  D.Qpump  = D.area * P.vSuct;                      // Durchfluss [m³/s]
  /* Thermodynamik der Wärmepumpe:
     P.power = thermische HEIZLEISTUNG P_heiz (Nutzwärme für das Fernwärmenetz).
     P_heiz = P_entzug + P_el, gekoppelt über die Leistungszahl COP:
       P_el      = P_heiz / COP            (elektrische Kompressorleistung)
       P_entzug  = P_heiz · (1 − 1/COP)    (dem Fluss entzogene Wärme)
     Nur P_entzug kühlt den Fluss ab. */
  D.Pel      = P.power / Math.max(P.cop, 1e-6);          // [kW]
  D.Pextract = P.power * (1 - 1/Math.max(P.cop, 1e-6));  // [kW]
  D.dT     = (D.Pextract*1000) / (RHO*CP*Math.max(D.Qpump,1e-6)); // Temperaturabsenkung [K]
  D.alpha  = ALPHA_MOL + 0.012 * P.vFlow * P.depth; // effektive (turbulente) Temp.-Leitf.
  D.Re     = P.vFlow * P.depth / NU;
  D.Pe     = P.vFlow * P.depth / D.alpha;
}

/* analytisches Geschwindigkeitsfeld an beliebiger Stelle (x,y,z) [m]
   = Grundströmung (Profil) + Ansaug-Senke + Auslass-Quelle (Strahl) */
function velocityAt(x,y,z, out){
  const dl=localDepth(y), b=bedZ(y);
  // Grundprofil: parabolisch über Breite (langsamer am Ufer),
  // Potenzgesetz 1/7 über die lokale Tiefe (langsamer an der Sohle)
  const yy=(2*y/W-1);
  const fy=1-0.5*yy*yy;
  const sg=Math.min(Math.max((z-b)/Math.max(dl,1e-3),0.001),1);  // 0=Sohle .. 1=Oberfläche
  const fz=Math.pow(sg,1/7);
  let u = (z<b-0.2) ? 0 : P.vFlow*fy*fz*1.15;   // unter der Sohle keine Strömung
  let v = 0, w = 0;

  const zp=pipeZ();
  // Ansaug-Senke am Einlass – PHYSIKALISCH skaliert: v_r = Q/(4π·r²)
  // (Punktsenke; die reale Einzugszone ist ein dünner Stromröhren-Schlauch A = Q/u)
  const SNK = D.Qpump/(4*Math.PI);
  const ix=X_IN-x, iy=Y_BANK-y, iz=zp-z;
  const r2i = ix*ix+iy*iy+iz*iz + 1.0;
  const sink = Math.min(SNK / (r2i*Math.sqrt(r2i)), 0.5);   // v_r/r -> Komponenten unten
  u += sink*ix; v += sink*iy; w += sink*iz;

  // Auslass-Quelle – Punktquelle gleicher Stärke (Nahfeld des Strahls
  // wird über die dedizierten Strahl-Partikel visualisiert)
  const ox=x-X_OUT, oy=y-Y_BANK, oz=z-zp;
  const r2o = ox*ox+oy*oy+oz*oz + 1.0;
  const src = Math.min(SNK / (r2o*Math.sqrt(r2o)), 0.5);
  u += src*ox; v += src*oy; w += src*oz;

  out.u=u; out.v=v; out.w=w;
}
const _vel={u:0,v:0,w:0};

/* trilineare Abtastung des Temperaturfeldes in Sigma-Koordinaten */
function sampleT(x,y,z){
  if (x < 0) return P.tIn;                       // Einström-Rand
  const dl=localDepth(y), b=bedZ(y);
  const sg=Math.min(Math.max((z-b)/Math.max(dl,1e-3),0),1);
  let gx = x/dx - 0.5, gy = y/dy - 0.5, gz = sg*NZ - 0.5;
  gx=Math.min(Math.max(gx,0),NX-1.001);
  gy=Math.min(Math.max(gy,0),NY-1.001);
  gz=Math.min(Math.max(gz,0),NZ-1.001);
  const i0=gx|0, j0=gy|0, k0=gz|0;
  const fx=gx-i0, fy=gy-j0, fz=gz-k0;
  const i1=i0+1, j1=j0+1, k1=k0+1;
  const c00=T[idx(i0,j0,k0)]*(1-fx)+T[idx(i1,j0,k0)]*fx;
  const c10=T[idx(i0,j1,k0)]*(1-fx)+T[idx(i1,j1,k0)]*fx;
  const c01=T[idx(i0,j0,k1)]*(1-fx)+T[idx(i1,j0,k1)]*fx;
  const c11=T[idx(i0,j1,k1)]*(1-fx)+T[idx(i1,j1,k1)]*fx;
  const c0=c00*(1-fy)+c10*fy, c1=c01*(1-fy)+c11*fy;
  return c0*(1-fz)+c1*fz;
}
/* physische z-Höhe eines Gitterpunkts (Sigma -> Höhe) */
function cellZ(j,k){ const y=(j+0.5)*dy; return bedZ(y) + (k+0.5)/NZ*localDepth(y); }

let tInletLocal = P.tIn, tOutlet = P.tIn;  // lokale Entnahme- / Auslasstemperatur

/* ein Solver-Schritt: semi-Lagrange-Advektion + Diffusion + Quellen/RB */
let _wzj=new Float64Array(0);   // Puffer für vertikale Diffusionsgewichte je Spalte
function solveStep(dt){
  // lokale Entnahmetemperatur (= Flusstemperatur am Einlass) und Auslasstemperatur
  tInletLocal = sampleT(X_IN, Y_BANK, pipeZ());
  tOutlet = Math.max(tInletLocal - D.dT, 0);     // abgekühlt; kein Eis (>= 0 °C)

  // --- Advektion (semi-Lagrange) mit ortsabhängiger Tiefe ---
  for (let k=0;k<NZ;k++){
    for (let j=0;j<NY;j++){
      const y=(j+0.5)*dy;
      for (let i=0;i<NX;i++){
        const x=(i+0.5)*dx, z=cellZ(j,k);
        velocityAt(x,y,z,_vel);
        T2[idx(i,j,k)] = sampleT(x-_vel.u*dt, y-_vel.v*dt, z-_vel.w*dt);
      }
    }
  }
  // --- Diffusion: KONSERVATIVE Fluss-Form (paarweiser Energieaustausch) ---
  // Gewichte aus PHYSIKALISCHEN Diffusivitäten (gitterunabhängig!):
  //   u* ≈ 0,1·v_m (Schubspannungsgeschw.), D_y = 0,6·h·u*, D_z = 0,067·h·u*,
  //   D_x klein (Längsvermischung dominiert die Advektion).
  //   w = 2·D·Δt/Δ² (Stabilitätskappe 0,45); vertikal je Spalte (Δz = h(y)/NZ).
  const uStar = 0.1*Math.max(P.vFlow,0.05);
  const Dy = 0.6*P.depth*uStar, Dz = 0.067*P.depth*uStar, Dx = 0.3*Dy;
  const wX = Math.min(2*Dx*dt/(dx*dx), 0.45);
  const wY = Math.min(2*Dy*dt/(dy*dy), 0.45);
  const WZJ = _wzj.length===NY ? _wzj : (_wzj=new Float64Array(NY));
  for(let j=0;j<NY;j++){
    const dz=Math.max(localDepth((j+0.5)*dy),MINDEPTH)/NZ;
    WZJ[j]=Math.min(2*Dz*dt/(dz*dz), 0.45);
  }
  // x-Richtung (Zellvolumina längs konstant)
  for (let k=0;k<NZ;k++) for (let j=0;j<NY;j++) for (let i=0;i<NX-1;i++){
    const a=idx(i,j,k), b=idx(i+1,j,k);
    const ex=wX*0.5*(T2[a]-T2[b]); T2[a]-=ex; T2[b]+=ex;
  }
  // y-Richtung (variable Säulenvolumina -> harmonisches Volumenmittel)
  for (let k=0;k<NZ;k++) for (let j=0;j<NY-1;j++){
    const Va=VCOL[j], Vb=VCOL[j+1], Vh=2*Va*Vb/(Va+Vb);
    for (let i=0;i<NX;i++){
      const a=idx(i,j,k), b=idx(i,j+1,k);
      const E=wY*0.5*Vh*(T2[a]-T2[b]); T2[a]-=E/Va; T2[b]+=E/Vb;
    }
  }
  // z-Richtung (Sigma-Schichten einer Säule gleich dick; Gewicht je Spalte)
  for (let j=0;j<NY;j++){ const wZ=WZJ[j];
    for (let k=0;k<NZ-1;k++) for (let i=0;i<NX;i++){
      const a=idx(i,j,k), b=idx(i,j,k+1);
      const ex=wZ*0.5*(T2[a]-T2[b]); T2[a]-=ex; T2[b]+=ex;
    }
  }

  // --- Quellterme & Randbedingungen ---
  // Einström-Rand x=0: gesamtes Flusswasser hat dieselbe Temperatur
  for (let k=0;k<NZ;k++) for (let j=0;j<NY;j++) T2[idx(0,j,k)] = P.tIn;
  // Auslass: ENERGIEERHALTENDE Kaltwasser-Einleitung (Quellterm q der Transportgleichung).
  // Pro Zeitschritt wird exakt das Wärmedefizit P_entzug·Δt in die Einleitzone eingebracht:
  //   T_Zelle -= P_entzug·Δt·ŵ / (ρ·c_p·V_Zelle),  Σŵ = 1,
  // geklippt bei T_aus (das Flusswasser kann lokal nicht kälter werden als das eingeleitete
  // Wasser). Der Strahl tritt senkrecht nach unten aus (Rohr von oben), prallt auf die Sohle
  // und breitet sich sohlnah aus -> Gewichte sind SOHLNAH zentriert (Prallstrahl).
  // Die Verdünnung mit der Strömung ergibt sich damit von selbst aus der Energiebilanz:
  // mehr Durchfluss am Auslass vorbei = gleiches Defizit auf mehr Wasser = wärmere Fahne.
  const oi = Math.round(X_OUT/dx), oj=Math.round(Y_BANK/dy);
  const tClip = Math.max(tOutlet, 0);
  {
    // kompakte Einmischzone des Prallstrahls: IN METERN definiert (≈5 m stromab,
    // ±2 m quer, volle Wassersäule sohlwärts gewichtet) -> gitterunabhängig
    const diN=Math.max(2,Math.round(5/dx)), djN=Math.max(1,Math.round(2/dy));
    let wsum=0; const cells=[];
    for (let dk=0;dk<NZ;dk++) for (let dj=-djN;dj<=djN;dj++) for (let di=0;di<=diN;di++){
      const i=oi+di, j=oj+dj, k=dk;
      if(i<0||i>=NX||j<0||j>=NY||k>=NZ) continue;
      const w = Math.exp(-( Math.pow(di*dx/2.0,2)*0.25
                          + Math.pow(dj*dy/1.0,2)*0.15
                          + (dk/NZ)*8*0.45 ));
      cells.push({c:idx(i,j,k), w, Vc:VCOL[j]}); wsum+=w;
    }
    // Defizit-Energie dieses Zeitschritts konservativ verteilen
    let Erest = D.Pextract*1000*dt;                    // [J] pro Solver-Schritt
    for(const cl of cells){
      const dTc = Erest>0 ? (D.Pextract*1000*dt)*(cl.w/wsum)/(RHO*CP*cl.Vc) : 0;
      const nv = T2[cl.c]-dTc;
      if(nv<tClip){ Erest -= (T2[cl.c]-tClip)*RHO*CP*cl.Vc; T2[cl.c]=tClip; }
      else        { Erest -= dTc*RHO*CP*cl.Vc;              T2[cl.c]=nv;   }
    }
    // Umlage: bei T_aus-Clipping verbleibende Energie auf Zellen mit Spielraum verteilen
    if(Erest > 1){
      let cap=0; for(const cl of cells) cap += Math.max(T2[cl.c]-tClip,0)*RHO*CP*cl.Vc;
      if(cap>1){ const f=Math.min(Erest/cap,1);
        for(const cl of cells) T2[cl.c] -= (T2[cl.c]-tClip)*f; }
    }
  }
  // (kein Luft-Wasser-Austausch: der Fluss ist überall gleich temperiert,
  //  damit ausschließlich der Effekt der Wärmepumpe sichtbar wird)

  const tmp=T; T=T2; T2=tmp;
}

/* mehrere Schritte bis ~Gleichgewicht laufen lassen */
function relax(steps){
  const dt = 1.2*dx/Math.max(P.vFlow,0.1);  // großzügiges dt (semi-Lagrange)
  for(let s=0;s<steps;s++) solveStep(dt);
}

/* Farbskala: FEST auf 0 … 25 °C (blau = 0 °C … rot = 25 °C, keine Autoskalierung). */
let Tmin=0, Tmax=25;
function fieldStats(){
  Tmin = 0; Tmax = 25;                  // feste Grenzen, unabhängig vom Feld
}

/* ---------- Farbskala (blau→cyan→grün→gelb→orange→rot) ------------------- */
function tempColor(t){ // t in [0,1] → {r,g,b} 0..1
  const stops=[
    [0.00,[0.10,0.25,0.75]],[0.25,[0.10,0.70,0.85]],[0.45,[0.20,0.80,0.40]],
    [0.65,[0.95,0.90,0.25]],[0.82,[0.98,0.55,0.15]],[1.00,[0.90,0.16,0.16]]];
  t=Math.min(Math.max(t,0),1);
  for(let i=1;i<stops.length;i++){
    if(t<=stops[i][0]){
      const a=stops[i-1],b=stops[i];const f=(t-a[0])/(b[0]-a[0]);
      return {r:a[1][0]+(b[1][0]-a[1][0])*f, g:a[1][1]+(b[1][1]-a[1][1])*f, b:a[1][2]+(b[1][2]-a[1][2])*f};
    }
  }
  const e=stops[stops.length-1][1]; return {r:e[0],g:e[1],b:e[2]};
}
const colNorm = t => tempColor((t-Tmin)/Math.max(Tmax-Tmin,1e-6));

/* =========================================================================
   (3) GEOMETRIE / SZENE
   ========================================================================= */
const sceneEl=document.getElementById('scene');
let renderer;
try{
  renderer=new THREE.WebGLRenderer({antialias:true});
  renderer.setPixelRatio(Math.min(window.devicePixelRatio,2));
  renderer.setClearColor(0x0b0f14,1);
  renderer.domElement.style.cssText='display:block;width:100%;height:100%';
  sceneEl.appendChild(renderer.domElement);
}catch(e){
  showError('WebGL konnte nicht initialisiert werden. Bitte einen Browser mit '+
            'aktivierter WebGL-Unterstützung verwenden.\n\n'+e.message);
}

const scene=new THREE.Scene();
scene.background=new THREE.Color(0x0b0f14);
scene.fog=new THREE.Fog(0x0b0f14, 700, 1800);  // weit entfernt: blendet Geometrie nicht aus

const camera=new THREE.PerspectiveCamera(50,1,0.5,5000);
camera.position.set(-90,150,220);
camera.lookAt(L/2,0,W/2);

// OrbitControls mit Fallback (falls CDN nicht verfügbar)
let controls;
if(typeof OrbitControls==='function'){
  controls=new OrbitControls(camera,renderer.domElement);
  controls.enableDamping=true; controls.dampingFactor=0.08;
  controls.target.set(L/2, 0, W/2);
}else{
  showError('OrbitControls nicht geladen – einfache Maussteuerung aktiv. '+
            'Bitte Internetverbindung/CDN prüfen.');
  controls=makeFallbackControls(camera,renderer.domElement,new THREE.Vector3(L/2,0,W/2));
}

/* einfache Orbit-Steuerung als Ersatz, falls OrbitControls fehlt */
function makeFallbackControls(cam,dom,target){
  const sph=new THREE.Spherical().setFromVector3(cam.position.clone().sub(target));
  let dragging=false,px=0,py=0;
  dom.addEventListener('mousedown',e=>{dragging=true;px=e.clientX;py=e.clientY;});
  window.addEventListener('mouseup',()=>dragging=false);
  window.addEventListener('mousemove',e=>{
    if(!dragging)return;
    sph.theta-=(e.clientX-px)*0.005; sph.phi-=(e.clientY-py)*0.005;
    sph.phi=Math.max(0.05,Math.min(Math.PI-0.05,sph.phi)); px=e.clientX;py=e.clientY;
  });
  dom.addEventListener('wheel',e=>{e.preventDefault();
    sph.radius*=(1+Math.sign(e.deltaY)*0.1);},{passive:false});
  return {target,update(){cam.position.copy(target).add(new THREE.Vector3().setFromSpherical(sph));cam.lookAt(target);}};
}

scene.add(new THREE.AmbientLight(0xffffff,0.65));
const sun=new THREE.DirectionalLight(0xffffff,0.8); sun.position.set(120,260,90); scene.add(sun);

/* Koordinaten-Mapping: Three.X = Fluss-x, Three.Y = Höhe (z), Three.Z = Fluss-y */
const V=(x,zh,y)=>new THREE.Vector3(x, zh, y);

/* leichte Unregelmäßigkeit der Sohle */
function noiseBed(x,y){ return Math.sin(x*0.06)*0.16+Math.cos(y*0.18)*0.12
   +Math.sin(x*0.013+y*0.07)*0.20; }

/* generisches parametrisches Flächen-Mesh über (x in 0..L, y in 0..W);
   fz(x,y) liefert die Höhe (Three.Y) */
function gridMesh(nx,ny,fz,material){
  const pos=[],ind=[],row=ny+1;
  for(let a=0;a<=nx;a++) for(let b=0;b<=ny;b++){ const x=a/nx*L, y=b/ny*W; pos.push(x,fz(x,y),y); }
  for(let a=0;a<nx;a++) for(let b=0;b<ny;b++){
    const i0=a*row+b,i1=(a+1)*row+b,i2=(a+1)*row+b+1,i3=a*row+b+1;
    ind.push(i0,i1,i2,i0,i2,i3);
  }
  const g=new THREE.BufferGeometry();
  g.setAttribute('position',new THREE.Float32BufferAttribute(pos,3));
  g.setIndex(ind); g.computeVertexNormals();
  return new THREE.Mesh(g,material);
}

/* Text-Sprite-Helfer (Maßstab/Beschriftung) */
function textSprite(txt,color){
  const c=document.createElement('canvas');c.width=256;c.height=64;
  const ctx=c.getContext('2d');ctx.fillStyle=color||'#8aa0b5';
  ctx.font='bold 40px Segoe UI, sans-serif';ctx.textAlign='center';ctx.textBaseline='middle';
  ctx.fillText(txt,128,32);
  const tex=new THREE.CanvasTexture(c);
  const sp=new THREE.Sprite(new THREE.SpriteMaterial({map:tex,transparent:true,depthTest:false}));
  sp.scale.set(26,6.5,1); return sp;
}

/* Umgebung: Ufer/Boden auf Höhe der Wasserlinie + parabolisches Flussbett + Maßstab */
const envGroup=new THREE.Group(); scene.add(envGroup);
let bedMesh=null;
function buildEnv(){
  envGroup.traverse(o=>{
    if(o.geometry) o.geometry.dispose();
    if(o.material){ if(o.material.map) o.material.map.dispose(); o.material.dispose(); }
  });
  envGroup.clear();
  const land=new THREE.MeshLambertMaterial({color:0x2f3a2a,side:THREE.DoubleSide});
  function bank(zc){ const m=new THREE.Mesh(new THREE.PlaneGeometry(L,44),land);
    m.rotation.x=-Math.PI/2; m.position.set(L/2,WSURF-0.02,zc); return m; }
  envGroup.add(bank(-22)); envGroup.add(bank(W+22));
  // parabolisches Flussbett (sandig)
  bedMesh=gridMesh(80,28,(x,y)=>bedZ(y)+noiseBed(x,y),
    new THREE.MeshLambertMaterial({color:0x8a7250,side:THREE.DoubleSide}));
  envGroup.add(bedMesh);
  // (kein Maßstab-Gitter auf der Wasseroberfläche)
  // Beschriftungen
  // Kilometrierung: kleine Schilder am Ufer (Pfosten + Tafel) statt schwebender Sprites
  {
    const poleMat=new THREE.MeshLambertMaterial({color:0x6b6f75});
    const mkSign=(txt,x,z)=>{
      // Pfosten
      const pole=new THREE.Mesh(new THREE.CylinderGeometry(0.42,0.42,7.2,10),poleMat);
      pole.position.set(x, WSURF+3.6, z); envGroup.add(pole);
      // Tafel mit Text (beidseitig sichtbar)
      const c=document.createElement('canvas'); c.width=256; c.height=96;
      const ctx=c.getContext('2d');
      ctx.fillStyle='#22303c'; ctx.fillRect(0,0,256,96);
      ctx.strokeStyle='#8aa0b5'; ctx.lineWidth=6; ctx.strokeRect(3,3,250,90);
      ctx.fillStyle='#dfe9f2';
      ctx.font='bold 48px Segoe UI, sans-serif'; ctx.textAlign='center'; ctx.textBaseline='middle';
      ctx.fillText(txt,128,50);
      const board=new THREE.Mesh(new THREE.PlaneGeometry(12.6,4.8),
        new THREE.MeshBasicMaterial({map:new THREE.CanvasTexture(c),side:THREE.DoubleSide}));
      board.position.set(x, WSURF+7.8, z);   // Tafel oben auf dem Pfosten
      envGroup.add(board);
    };
    [0,100,200,300,400].forEach(x=>mkSign(x+' m', x, -4));
  }
  // Wohnbebauung auf der Pumpenseite (für spätere Fernwärme-Anschlüsse):
  // Mehrfamilienhäuser (groß, Flachdach) und Einzelhäuser (klein, Spitzdach),
  // jeweils mit Fenstern und Tür auf der Flussseite (Fassaden-Textur)
  {
    const wallM=new THREE.MeshLambertMaterial({color:0xb9b2a6});
    const wallE=new THREE.MeshLambertMaterial({color:0xc9bfa8});
    const roofF=new THREE.MeshLambertMaterial({color:0x4a5560});   // Flachdach
    const roofS=new THREE.MeshLambertMaterial({color:0x8a4a3a});   // Spitzdach (Ziegel)
    // Fassade: Fenstergitter + Tür auf transparenter Canvas -> Ebene knapp vor der Wand
    const mkFacade=(wx,h,floors,cols,x,y,z)=>{
      const c=document.createElement('canvas'); c.width=256; c.height=256;
      const ctx=c.getContext('2d');
      const cw=256/cols, ch=256/(floors+0.4);
      for(let f=0;f<floors;f++)for(let k=0;k<cols;k++){
        // Erdgeschoss Mitte: Tür statt Fenster
        const ground=(f===floors-1), mid=(k===(cols>>1));
        if(ground&&mid){
          ctx.fillStyle='#4b3a28';
          ctx.fillRect(k*cw+cw*0.28, (f+0.18)*ch+ch*0.10, cw*0.44, ch*0.9);
        }else{
          ctx.fillStyle='#cfe4f5';
          ctx.fillRect(k*cw+cw*0.22, (f+0.18)*ch+ch*0.16, cw*0.56, ch*0.55);
          ctx.strokeStyle='#5c7186'; ctx.lineWidth=3;
          ctx.strokeRect(k*cw+cw*0.22, (f+0.18)*ch+ch*0.16, cw*0.56, ch*0.55);
        }
      }
      const m=new THREE.Mesh(new THREE.PlaneGeometry(wx-0.6,h-0.4),
        new THREE.MeshBasicMaterial({map:new THREE.CanvasTexture(c),transparent:true}));
      m.position.set(x, y, z);                     // knapp vor der flussseitigen Wand (+z)
      envGroup.add(m);
    };
    // Mehrfamilienhaus: großer Quader + Flachdach-Platte + Fassade
    const mkMFH=(x,z,wx,wz,h)=>{
      const b=new THREE.Mesh(new THREE.BoxGeometry(wx,h,wz),wallM);
      b.position.set(x,WSURF+h/2,z); envGroup.add(b);
      const r=new THREE.Mesh(new THREE.BoxGeometry(wx+0.8,0.7,wz+0.8),roofF);
      r.position.set(x,WSURF+h+0.35,z); envGroup.add(r);
      mkFacade(wx,h,Math.max(3,Math.round(h/3.4)),Math.max(4,Math.round(wx/4.5)),
               x, WSURF+h/2, z+wz/2+0.05);
    };
    // Einzelhaus: kleiner Quader + Pyramidendach + Fassade (1 Etage: 2 Fenster + Tür)
    const mkEFH=(x,z)=>{
      const b=new THREE.Mesh(new THREE.BoxGeometry(9,5,8),wallE);
      b.position.set(x,WSURF+2.5,z); envGroup.add(b);
      const r=new THREE.Mesh(new THREE.ConeGeometry(6.6,3.4,4),roofS);
      r.rotation.y=Math.PI/4;                       // Pyramide über den rechteckigen Grundriss
      r.position.set(x,WSURF+5+1.7,z); envGroup.add(r);
      mkFacade(9,5,1,3, x, WSURF+2.5, z+4+0.05);
    };
    // Mehrfamilienhäuser (stromauf der Pumpe)
    mkMFH( 40,-30, 26,12,14);
    mkMFH( 78,-32, 22,12,11);
    // Einzelhäuser (stromab der Pumpe)
    mkEFH(255,-27);
    mkEFH(285,-30);
    mkEFH(320,-26);
    mkEFH(355,-30);
  }
  // Strömungsrichtung: flacher 50-m-Pfeil auf dem Rasen (gegenüberliegendes Ufer),
  // daneben die Beschriftung "Strömungsrichtung" – beides liegend, entlang des Flusses
  {
    const arrMat=new THREE.MeshBasicMaterial({color:0x3fc1c9,side:THREE.DoubleSide});
    const sh=new THREE.Shape();                          // Pfeil: Schaft + Spitze, 50 m lang
    sh.moveTo(-25,-1.5); sh.lineTo(17,-1.5); sh.lineTo(17,-4.5); sh.lineTo(25,0);
    sh.lineTo(17,4.5);  sh.lineTo(17,1.5);  sh.lineTo(-25,1.5); sh.closePath();
    const arrow=new THREE.Mesh(new THREE.ShapeGeometry(sh),arrMat);
    arrow.rotation.x=-Math.PI/2;                         // flach auf den Rasen legen
    arrow.position.set(L/2, WSURF+0.06, W+8);
    envGroup.add(arrow);
    // Beschriftung daneben, ebenfalls flach, entlang des Pfeils lesbar
    const c=document.createElement('canvas'); c.width=1024; c.height=128;
    const ctx=c.getContext('2d');
    ctx.fillStyle='#3fc1c9';
    ctx.font='bold 88px Segoe UI, sans-serif'; ctx.textAlign='center'; ctx.textBaseline='middle';
    ctx.fillText('Strömungsrichtung',512,66);
    const lbl=new THREE.Mesh(new THREE.PlaneGeometry(46,5.75),
      new THREE.MeshBasicMaterial({map:new THREE.CanvasTexture(c),transparent:true,side:THREE.DoubleSide}));
    lbl.rotation.x=-Math.PI/2;
    lbl.position.set(L/2, WSURF+0.06, W+15);
    envGroup.add(lbl);
  }
}

/* Wasservolumen (semitransparent blau) mit parabolischem Boden, flacher Oberfläche */
let waterMesh=null;
function buildWater(){
  if(waterMesh){scene.remove(waterMesh);waterMesh.traverse(o=>{if(o.geometry)o.geometry.dispose();});}
  waterMesh=new THREE.Group();
  const mat=new THREE.MeshPhongMaterial({color:0x1e6fb0,transparent:true,opacity:0.24,
    depthWrite:false,shininess:60,side:THREE.DoubleSide});
  // flache Oberfläche
  waterMesh.add(gridMesh(40,20,()=>WSURF-0.02,mat));
  // Stirn- und Seitenflächen (Band zwischen Sohle und Oberfläche)
  const seg=28;
  function strip(low,high){
    const pos=[],ind=[];
    for(let b=0;b<=seg;b++){const a=low(b),o=high(b);pos.push(a.x,a.y,a.z,o.x,o.y,o.z);}
    for(let b=0;b<seg;b++){const i=b*2;ind.push(i,i+1,i+2,i+1,i+3,i+2);}
    const g=new THREE.BufferGeometry();
    g.setAttribute('position',new THREE.Float32BufferAttribute(pos,3));
    g.setIndex(ind); g.computeVertexNormals();
    return new THREE.Mesh(g,mat);
  }
  for(const xp of [0,L])  // Stirnseiten (zeigen den parabolischen Querschnitt)
    waterMesh.add(strip(b=>{const y=b/seg*W;return {x:xp,y:bedZ(y),z:y};},
                        b=>{const y=b/seg*W;return {x:xp,y:WSURF,z:y};}));
  for(const yp of [0,W])  // Seitenwände
    waterMesh.add(strip(b=>{const x=b/seg*L;return {x:x,y:bedZ(yp),z:yp};},
                        b=>{const x=b/seg*L;return {x:x,y:WSURF,z:yp};}));
  scene.add(waterMesh);
}

/* Pumpengebäude + geschlossener Rohrkreislauf (Fluss → Pumpe → Fluss) */
const pumpGroup=new THREE.Group(); scene.add(pumpGroup);

// Kreislaufpfad (Polylinie) für die Partikel im Rohr
let circuitPath=[], circuitPumpIdx=0, circuitCum=[], circuitTotal=1, pipeRadius=0.4;
function buildCircuitLengths(){
  circuitCum=[0];
  for(let i=0;i<circuitPath.length-1;i++)
    circuitCum.push(circuitCum[i]+circuitPath[i].distanceTo(circuitPath[i+1]));
  circuitTotal=Math.max(circuitCum[circuitCum.length-1],1e-6);
}
/* Position bei Bogenlängen-Anteil s∈[0,1] entlang des Pfades */
function samplePath(s,out){
  const d=s*circuitTotal; let i=0;
  while(i<circuitCum.length-2 && circuitCum[i+1]<d) i++;
  const segLen=Math.max(circuitCum[i+1]-circuitCum[i],1e-6);
  const t=(d-circuitCum[i])/segLen;
  out.lerpVectors(circuitPath[i], circuitPath[i+1], t);
  return circuitCum[i]+t*segLen;   // zurück: zurückgelegte Länge (für Pumpen-Index)
}
// Tiefe (Bogenlänge) ab der das Wasser abgekühlt ist
function circuitCoolFrac(){ return (circuitCum[circuitPumpIdx]||0)/circuitTotal; }
// Zylinder zwischen zwei Punkten
function tube(p1,p2,r,mat){
  const dir=new THREE.Vector3().subVectors(p2,p1); const len=dir.length();
  const m=new THREE.Mesh(new THREE.CylinderGeometry(r,r,len,14),mat);
  m.position.copy(p1).add(p2).multiplyScalar(0.5);
  m.quaternion.setFromUnitVectors(new THREE.Vector3(0,1,0),dir.clone().normalize());
  return m;
}
function tubePath(pts,r,mat){
  const g=new THREE.Group();
  for(let i=0;i<pts.length-1;i++) g.add(tube(pts[i],pts[i+1],r,mat));
  for(let i=1;i<pts.length-1;i++){ const s=new THREE.Mesh(new THREE.SphereGeometry(r,10,10),mat);
    s.position.copy(pts[i]); g.add(s); }   // glatte Knoten an den Ecken
  return g;
}
function buildPump(){
  // alte Geometrien/Texturen sauber freigeben, dann Gruppe leeren (kein Speicherleck)
  pumpGroup.traverse(o=>{
    if(o.geometry) o.geometry.dispose();
    if(o.material){ if(o.material.map) o.material.map.dispose(); o.material.dispose(); }
  });
  pumpGroup.clear();
  const xMid=(X_IN+X_OUT)/2, by=WSURF, bz=-14;
  // Gebäude auf dem Ufer
  const b=new THREE.Mesh(new THREE.BoxGeometry(22,10,16),new THREE.MeshLambertMaterial({color:0x9aa7b3}));
  b.position.set(xMid, by+5, bz); pumpGroup.add(b);
  const roof=new THREE.Mesh(new THREE.BoxGeometry(23,1,17),new THREE.MeshLambertMaterial({color:0x46525e}));
  roof.position.set(xMid, by+10.2, bz); pumpGroup.add(roof);
  // Beschriftung direkt AUF der flussseitigen Wand (Textur-Ebene, kein schwebendes Sprite)
  {
    const c=document.createElement('canvas'); c.width=512; c.height=128;
    const ctx=c.getContext('2d');
    ctx.fillStyle='#2c3640';
    ctx.font='bold 64px Segoe UI, sans-serif'; ctx.textAlign='center'; ctx.textBaseline='middle';
    ctx.fillText('Wärmepumpe',256,64);
    const tex=new THREE.CanvasTexture(c);
    const wall=new THREE.Mesh(new THREE.PlaneGeometry(18,4.5),
      new THREE.MeshBasicMaterial({map:tex,transparent:true}));
    wall.position.set(xMid, by+5.4, bz+8.06);   // knapp vor der Wand (z = bz+8)
    pumpGroup.add(wall);
  }
  /* Fernwärme-Kreislauf: zwei parallele Rohrstücke von der landseitigen Gebäudewand
     in den Vordergrund (zum Kunden). Wärmetauscher-Wirkungsgrad = 1,0 (rein visuell:
     feste Farben repräsentieren die Temperaturniveaus, keine Solver-Änderung nötig). */
  {
    const wallZ = bz-8;                                  // landseitige Wand (z = bz−8)
    // Wandbeschriftung (Textur-Ebene auf der landseitigen Wand, Blick nach −z)
    const wallText=(txt,color,x,y)=>{
      const c=document.createElement('canvas'); c.width=512; c.height=96;
      const ctx=c.getContext('2d');
      ctx.fillStyle=color;
      ctx.font='bold 64px Segoe UI, sans-serif'; ctx.textAlign='center'; ctx.textBaseline='middle';
      ctx.fillText(txt,256,48);
      const m=new THREE.Mesh(new THREE.PlaneGeometry(6,1.15),
        new THREE.MeshBasicMaterial({map:new THREE.CanvasTexture(c),transparent:true}));
      m.rotation.y=Math.PI;                              // zur Landseite (−z) ausgerichtet
      m.position.set(x, y, wallZ-0.06);
      pumpGroup.add(m);
    };
    // Vorlauf: heißes Medium (90 °C) verlässt das Gebäude – mattes kräftiges Rot
    const matVor  = new THREE.MeshLambertMaterial({color:0xff3b30});
    pumpGroup.add(tubePath([V(xMid+4, by+4, wallZ), V(xMid+4, by+4, wallZ-14)], 0.5, matVor));
    wallText('FW-Vorlauf','#a61712', xMid+4, by+5.6);
    // Rücklauf: abgekühltes Medium (60 °C) kommt vom Kunden zurück – mattes kräftiges Blau
    const matRueck= new THREE.MeshLambertMaterial({color:0x2b75ff});
    pumpGroup.add(tubePath([V(xMid-4, by+4, wallZ-14), V(xMid-4, by+4, wallZ)], 0.5, matRueck));
    wallText('FW-Rücklauf','#123f9e', xMid-4, by+5.6);
  }

  const r=Math.max(P.dPipe/2,0.15);
  const rVis=r*1.8;                    // sichtbarer Rohrradius: dicker, damit das Rohr Volumen hat
  pipeRadius=rVis;
  const zp=pipeZ();                  // Mündungstiefe ≈ 0,5 m unter der Oberfläche
  // durchsichtige Rohre, damit man die Partikel im Inneren über die ganze Länge sieht
  const matIn =new THREE.MeshPhongMaterial({color:0x8fc0e6,shininess:50,transparent:true,opacity:0.30,depthWrite:false,side:THREE.DoubleSide});
  const matOut=new THREE.MeshPhongMaterial({color:0x6f9bff,shininess:50,transparent:true,opacity:0.30,depthWrite:false,side:THREE.DoubleSide});
  // Entnahme (Mündung im Wasser → Pumpe): senkrechtes Eintauchen, dann erhöht zur Wandmitte
  const inPath =[
    V(X_IN,     zp,    Y_BANK),       // Mündung: Rohr endet ~0,5 m tief im Wasser
    V(X_IN,     by+5,  Y_BANK),       // senkrecht von oben ins Wasser (ohne Ufer zu berühren)
    V(X_IN,     by+5,  bz),           // erhöht quer übers Ufer hinweg bis zur Pumpenachse
    V(xMid-11,  by+5,  bz)            // mittig in die stromaufwärtige Pumpenwand
  ];
  // Rückgabe (Pumpe → Mündung): symmetrisch
  const outPath=[
    V(xMid+11,  by+5,  bz),           // mittig aus der stromabwärtigen Pumpenwand
    V(X_OUT,    by+5,  bz),
    V(X_OUT,    by+5,  Y_BANK),       // erhöht über das Wasser
    V(X_OUT,    zp,    Y_BANK)        // Mündung: Rohr endet ~0,5 m tief im Wasser
  ];
  pumpGroup.add(tubePath(inPath, rVis, matIn));
  pumpGroup.add(tubePath(outPath,rVis, matOut));
  // vollständiger Kreislaufpfad durch die Pumpe (für die Kreislauf-Partikel)
  circuitPath = inPath.concat([V(xMid, by+5, bz)], outPath);
  circuitPumpIdx = inPath.length;     // ab hier ist das Medium abgekühlt
  buildCircuitLengths();
}

/* =========================================================================
   (4) VISUALISIERUNG: Heatmap-Schnittebenen, Partikel-Stromlinien, Vektoren
   ========================================================================= */
/* generische Heatmap-Ebene über CanvasTexture */
function makeHeatPlane(wRes,hRes){
  const cv=document.createElement('canvas'); cv.width=wRes; cv.height=hRes;
  const ctx=cv.getContext('2d');
  const img=ctx.createImageData(wRes,hRes);
  const tex=new THREE.CanvasTexture(cv);
  tex.magFilter=THREE.LinearFilter; tex.minFilter=THREE.LinearFilter;
  const mat=new THREE.MeshBasicMaterial({map:tex,transparent:true,opacity:0.95,side:THREE.DoubleSide});
  const mesh=new THREE.Mesh(new THREE.PlaneGeometry(1,1),mat);
  return {cv,ctx,img,tex,mesh};
}
// kontinuierlich abgetastete Ebenen (glatt, mit Maskierung der Sohle)
const planeCross=makeHeatPlane(90,60);    // Querschnitt (y × z)
const planeHoriz=makeHeatPlane(180,60);   // Horizontalschnitt (x × y)
const planeLong =makeHeatPlane(200,60);   // Längsschnitt  (x × z) entlang der Strömung
scene.add(planeCross.mesh,planeHoriz.mesh,planeLong.mesh);

/* Heatmap-Pixel füllen; getT(a,b) liefert Temperatur oder null/NaN = transparent */
function paint(plane, getT){
  const {img}=plane, w=img.width, h=img.height, d=img.data;
  for(let b=0;b<h;b++) for(let a=0;a<w;a++){
    const t=getT(a,b); const o=(b*w+a)*4;
    if(t==null || t!==t){ d[o+3]=0; continue; }   // unter der Sohle: transparent
    const c=colNorm(t); d[o]=c.r*255; d[o+1]=c.g*255; d[o+2]=c.b*255; d[o+3]=235;
  }
  plane.ctx.putImageData(img,0,0); plane.tex.needsUpdate=true;
}
function updatePlanes(){
  // Querschnitt an x=xCut: y × z, parabolische Sohle maskiert.
  // Breite gespiegelt, damit die angezeigte Seite der tatsächlichen Flussseite entspricht
  if(planeCross.mesh.visible) paint(planeCross,(a,b)=>{
    const y=(1 - a/(planeCross.cv.width-1))*W, z=WSURF-(b/(planeCross.cv.height-1))*P.depth;
    return (z<bedZ(y)) ? null : sampleT(P.xCut,y,z);
  });
  // Horizontalschnitt in Tiefe zCut (z = Oberfläche - zCut)
  if(planeHoriz.mesh.visible){ const zh=WSURF-P.zCut; paint(planeHoriz,(a,b)=>{
    const x=a/(planeHoriz.cv.width-1)*L, y=b/(planeHoriz.cv.height-1)*W;
    return (zh<bedZ(y)) ? null : sampleT(x,y,zh);
  });}
  // Längsschnitt entlang der Strömung bei y=yLong: x × z, Sohle maskiert
  if(planeLong.mesh.visible) paint(planeLong,(a,b)=>{
    const x=a/(planeLong.cv.width-1)*L, z=WSURF-(b/(planeLong.cv.height-1))*P.depth;
    return (z<bedZ(P.yLong)) ? null : sampleT(x,P.yLong,z);
  });
}
/* Schnittebenen im Raum positionieren/orientieren */
function placePlanes(){
  planeCross.mesh.geometry.dispose();
  planeCross.mesh.geometry=new THREE.PlaneGeometry(W,P.depth);
  planeCross.mesh.rotation.set(0,Math.PI/2,0);
  planeCross.mesh.position.set(P.xCut,WSURF-P.depth/2,W/2);
  planeHoriz.mesh.geometry.dispose();
  planeHoriz.mesh.geometry=new THREE.PlaneGeometry(L,W);
  planeHoriz.mesh.rotation.set(-Math.PI/2,0,0);
  planeHoriz.mesh.position.set(L/2,WSURF-P.zCut,W/2);
  planeLong.mesh.geometry.dispose();
  planeLong.mesh.geometry=new THREE.PlaneGeometry(L,P.depth);
  planeLong.mesh.rotation.set(0,0,0);                       // x–z-Ebene (senkrecht, längs)
  planeLong.mesh.position.set(L/2,WSURF-P.depth/2,P.yLong);
}

/* runde (kugelartige) Partikel: weiche Kreis-Textur als Sprite */
function discTexture(){
  const c=document.createElement('canvas'); c.width=c.height=64;
  const g=c.getContext('2d');
  const grd=g.createRadialGradient(32,32,2,32,32,30);
  grd.addColorStop(0,'rgba(255,255,255,1)');
  grd.addColorStop(0.55,'rgba(255,255,255,0.95)');
  grd.addColorStop(1,'rgba(255,255,255,0)');
  g.fillStyle=grd; g.beginPath(); g.arc(32,32,30,0,Math.PI*2); g.fill();
  const t=new THREE.CanvasTexture(c); return t;
}
const DISC=discTexture();

/* ---------------------------------------------------------------------------
   Particle Tracer für die Strömung (masselose Partikel):
   - feste Anzahl, anfangs zufällig im Gebiet verteilt
   - je Frame Geschwindigkeit aus dem (interpolierten) FVM-Feld ablesen
   - Position per Euler-Integration: x += u·dt, y += v·dt, z += w·dt
   - verlässt ein Partikel das Gebiet / steckt fest -> Neustart am Einlass (x≈0)
   - Einfärbung nach lokaler Wassertemperatur (kalt = blau, warm = rot)
   --------------------------------------------------------------------------- */
const NP=2600;
const pPos=new Float32Array(NP*3), pCol=new Float32Array(NP*3);
const partGeo=new THREE.BufferGeometry();
partGeo.setAttribute('position',new THREE.BufferAttribute(pPos,3));
partGeo.setAttribute('color',new THREE.BufferAttribute(pCol,3));
const partMat=new THREE.PointsMaterial({size:1.7,map:DISC,alphaTest:0.4,vertexColors:true,transparent:true,opacity:0.95,sizeAttenuation:true,depthWrite:false});
const particles=new THREE.Points(partGeo,partMat);
particles.frustumCulled=false;        // nie wegblenden (auch beim Hineinzoomen sichtbar)
scene.add(particles);
const pX=new Float32Array(NP),pY=new Float32Array(NP),pZ=new Float32Array(NP);
function seedParticle(n,atInlet){
  pX[n] = atInlet ? Math.random()*3 : Math.random()*L;          // Einlass (x≈0) oder zufällig
  const y = 0.5 + Math.random()*(W-1); pY[n]=y;
  pZ[n] = bedZ(y) + (0.08+0.86*Math.random())*localDepth(y);    // in der Wassersäule
}
for(let n=0;n<NP;n++) seedParticle(n,false);
function updateParticles(dt){
  const zp=pipeZ();
  for(let n=0;n<NP;n++){
    velocityAt(pX[n],pY[n],pZ[n],_vel);                         // FVM-Geschwindigkeit interpolieren
    const sp=Math.hypot(_vel.u,_vel.v,_vel.w);
    pX[n]+=_vel.u*dt; pY[n]+=_vel.v*dt; pZ[n]+=_vel.w*dt;        // Euler-Integration
    const di=pX[n]-X_IN, dj=pY[n]-Y_BANK, dk=pZ[n]-zp;
    const sucked=(di*di+dj*dj+dk*dk)<2.5;                       // in die Ansaugung gezogen
    // Grenzbedingung: Gebiet verlassen / unter Sohle / steckt fest -> am Einlass neu
    if(!(pX[n]<L)||pX[n]<0||pY[n]<0||pY[n]>W||pZ[n]<bedZ(pY[n])||pZ[n]>WSURF||sp<0.015||sucked||!isFinite(pX[n]))
      seedParticle(n,true);
    const o=n*3; pPos[o]=pX[n]; pPos[o+1]=pZ[n]; pPos[o+2]=pY[n];
    const c=colNorm(sampleT(pX[n],pY[n],pZ[n]));                // Farbe = lokale Wassertemperatur
    pCol[o]=c.r; pCol[o+1]=c.g; pCol[o+2]=c.b;
  }
  partGeo.attributes.position.needsUpdate=true;
  partGeo.attributes.color.needsUpdate=true;
}

/* Kreislauf-Partikel: strömen durch die (durchsichtigen) Rohre und die Pumpe */
const NCP=900;
const cpS=new Float32Array(NCP);
const cpOff=new Float32Array(NCP*3);   // fester Querversatz (Einheitskugel) -> füllt den Rohrquerschnitt
const cpPos=new Float32Array(NCP*3), cpCol=new Float32Array(NCP*3);
for(let n=0;n<NCP;n++){
  cpS[n]=n/NCP;
  // zufälliger Punkt in einer Kugel (radial im Rohr verteilt)
  let ox,oy,oz,d2;
  do{ ox=Math.random()*2-1; oy=Math.random()*2-1; oz=Math.random()*2-1; d2=ox*ox+oy*oy+oz*oz; }while(d2>1);
  cpOff[n*3]=ox; cpOff[n*3+1]=oy; cpOff[n*3+2]=oz;
}
const circGeo=new THREE.BufferGeometry();
circGeo.setAttribute('position',new THREE.BufferAttribute(cpPos,3));
circGeo.setAttribute('color',new THREE.BufferAttribute(cpCol,3));
const circParticles=new THREE.Points(circGeo,
  new THREE.PointsMaterial({size:1.1,map:DISC,alphaTest:0.35,vertexColors:true,transparent:true,opacity:0.98,sizeAttenuation:true,depthWrite:false}));
circParticles.frustumCulled=false; scene.add(circParticles);
const _cp=new THREE.Vector3();
function updateCircuit(dt){
  if(!circuitPath.length) return;
  const coolFrac=circuitCoolFrac();
  const cWarm=colNorm(tInletLocal), cCold=colNorm(tOutlet);
  const adv=dt*0.05;                    // Vortriebsgeschwindigkeit entlang des Rohres
  const off=pipeRadius*0.35;   // deutlich innerhalb der (dickeren) Rohrwand
  for(let n=0;n<NCP;n++){
    cpS[n]+=adv; if(cpS[n]>=1) cpS[n]-=1;
    samplePath(cpS[n],_cp);
    const o=n*3;
    cpPos[o]  =_cp.x+cpOff[o]  *off;
    cpPos[o+1]=_cp.y+cpOff[o+1]*off;
    cpPos[o+2]=_cp.z+cpOff[o+2]*off;
    const c=(cpS[n]<coolFrac)?cWarm:cCold;   // nach der Pumpe: abgekühlt
    cpCol[o]=c.r; cpCol[o+1]=c.g; cpCol[o+2]=c.b;
  }
  circGeo.attributes.position.needsUpdate=true;
  circGeo.attributes.color.needsUpdate=true;
}

/* ---------------------------------------------------------------------------
   Auslassstrahl (physikalisch): Das Rohr mündet SENKRECHT von oben, der Strahl
   tritt daher NACH UNTEN aus – runder Freistrahl mit ~12° Halbwinkel und
   Geschwindigkeitsabfall ~ 1/s. Nach ~1–2 m trifft er auf die Sohle (Prallstrahl)
   und wird in einen radial auswärts gerichteten, abklingenden WANDSTRAHL
   umgelenkt, den die Flussströmung stromab trägt.
   --------------------------------------------------------------------------- */
const NJ=1500;
const jX=new Float32Array(NJ), jY=new Float32Array(NJ), jZ=new Float32Array(NJ);
const jU=new Float32Array(NJ), jV=new Float32Array(NJ), jW=new Float32Array(NJ);
const jAge=new Float32Array(NJ), jHit=new Uint8Array(NJ), jS=new Float32Array(NJ), jLife=new Float32Array(NJ);
const jPos=new Float32Array(NJ*3), jCol=new Float32Array(NJ*3);
const jetGeo=new THREE.BufferGeometry();
jetGeo.setAttribute('position',new THREE.BufferAttribute(jPos,3));
jetGeo.setAttribute('color',new THREE.BufferAttribute(jCol,3));
const jetParticles=new THREE.Points(jetGeo,
  new THREE.PointsMaterial({size:1.7,map:DISC,alphaTest:0.4,vertexColors:true,transparent:true,opacity:0.96,sizeAttenuation:true,depthWrite:false}));
jetParticles.frustumCulled=false; scene.add(jetParticles);
const JET_LIFE=34;
function seedJet(n,fresh){
  const zp=pipeZ(), rM=P.dPipe/2;
  // Start in der Rohrmündung (Kreisquerschnitt)
  const a=Math.random()*Math.PI*2, rr=rM*Math.sqrt(Math.random());
  jX[n]=X_OUT + Math.cos(a)*rr;
  jY[n]=Y_BANK + Math.sin(a)*rr;
  jZ[n]=zp;
  // Richtung: senkrecht NACH UNTEN, Freistrahl-Kegel ~±12°
  const th=Math.random()*Math.PI*2, tilt=Math.tan((Math.random()*12)*Math.PI/180);
  const spd=P.vSuct*(0.85+0.15*Math.random());
  const dxn=Math.cos(th)*tilt, dyn=Math.sin(th)*tilt, dzn=-1;
  const inv=spd/Math.hypot(dxn,dyn,dzn);
  jU[n]=dxn*inv; jV[n]=dyn*inv; jW[n]=dzn*inv;
  jHit[n]=0;
  jLife[n]=JET_LIFE*(0.45+1.1*Math.random());      // INDIVIDUELLE Lebensdauer -> Wolke
                                                    // blendet aus statt an einer Kante zu enden
  jAge[n]=fresh?0:Math.random()*jLife[n];           // beim Start zeitlich gestaffelt
  jS[n]=fresh?0:jAge[n]*P.vSuct*0.6;               // zurückgelegter Weg (für Entrainment)
}
for(let n=0;n<NJ;n++) seedJet(n,false);
function updateJet(dt){
  for(let n=0;n<NJ;n++){
    velocityAt(jX[n],jY[n],jZ[n],_vel);
    // Strahlgeschwindigkeit klingt ab: Freistrahl u_c ~ 1/s, Wandstrahl schneller
    const decay=Math.exp(-jAge[n]/(jHit[n]?5:9));
    // turbulente Dispersion: wächst, wenn der Strahlimpuls abklingt
    // (D_y ~ 0,1 m²/s -> Schrittweite sqrt(2·D·dt) ≈ 0,28 m, im Kern unterdrückt)
    const dj_=0.32*(1-decay);
    const rx=(Math.random()-0.5)*2*dj_, ry=(Math.random()-0.5)*2*dj_, rz=(Math.random()-0.5)*dj_;
    const mvx=(_vel.u + jU[n]*decay)*dt + rx,
          mvy=(_vel.v + jV[n]*decay)*dt + ry,
          mvz=(_vel.w + jW[n]*decay)*dt + rz;
    jX[n]+=mvx; jY[n]+=mvy; jZ[n]+=mvz;
    jS[n]+=Math.hypot(mvx,mvy,mvz);                 // zurückgelegte Weglänge s
    jAge[n]+=dt*0.5;
    // Prallstrahl: an der Sohle in radialen Wandstrahl umlenken
    const bz=bedZ(jY[n]);
    if(!jHit[n] && jZ[n]<=bz+0.25){
      jZ[n]=bz+0.25; jHit[n]=1;
      const rx=jX[n]-X_OUT, ry=jY[n]-Y_BANK;
      const rl=Math.max(Math.hypot(rx,ry),0.15);
      const spd=Math.hypot(jU[n],jV[n],jW[n])*0.85;   // Impulsverlust beim Aufprall
      jU[n]=spd*rx/rl; jV[n]=spd*ry/rl; jW[n]=0;      // radial auswärts entlang der Sohle
    }
    if(jAge[n]>jLife[n] || !(jX[n]<L) || jX[n]<0 || jY[n]<0 || jY[n]>W ||
       jZ[n]<bedZ(jY[n]) || jZ[n]>WSURF || !isFinite(jX[n]))
      seedJet(n,true);
    const o=n*3; jPos[o]=jX[n]; jPos[o+1]=jZ[n]; jPos[o+2]=jY[n];
    // Farbe = ENTRAINMENT-Mischtemperatur des Strahlwassers:
    // Im Potentialkern (s < ~6,2·d) ist das Strahlwasser noch unverdünnt bei T_aus;
    // danach sinkt der Quellwasseranteil wie c ≈ 6,2·d/s (runder Freistrahl).
    // Umgebung = lokales Feld -> stromab konvergiert die Farbe exakt zur Ebenenfarbe.
    const cFrac=Math.min(1, 6.2*P.dPipe/Math.max(jS[n],1e-3));
    const tmix=cFrac*tOutlet + (1-cFrac)*sampleT(jX[n],jY[n],jZ[n]);
    const c=colNorm(tmix);
    jCol[o]=c.r; jCol[o+1]=c.g; jCol[o+2]=c.b;
  }
  jetGeo.attributes.position.needsUpdate=true;
  jetGeo.attributes.color.needsUpdate=true;
}

/* ---------------------------------------------------------------------------
   Ansaugung am Einlass (physikalisch): Die Senkenströmung ist v_r = Q/(4π·r²)
   und gegen die Flussströmung schnell vernachlässigbar. Eingesaugt wird daher
   nur der dünne STROMRÖHREN-SCHLAUCH mit Querschnitt A = Q/u (Radius ~0,5–1 m)
   stromauf der Mündung. Die Partikel folgen dem ECHTEN Geschwindigkeitsfeld
   (inkl. physikalischer Senke, ohne künstliche Verstärkung): Partikel im
   Schlauch werden gefangen, die übrigen ziehen knapp vorbei.
   --------------------------------------------------------------------------- */
const NS=380;
const sX=new Float32Array(NS), sY=new Float32Array(NS), sZ=new Float32Array(NS);
const sEsc=new Float32Array(NS);        // individuelle (zufällige) Auslaufgrenze -> keine sichtbare Kante
const sPos=new Float32Array(NS*3), sCol=new Float32Array(NS*3);
const suctGeo=new THREE.BufferGeometry();
suctGeo.setAttribute('position',new THREE.BufferAttribute(sPos,3));
suctGeo.setAttribute('color',new THREE.BufferAttribute(sCol,3));
const suctParticles=new THREE.Points(suctGeo,
  new THREE.PointsMaterial({size:1.5,map:DISC,alphaTest:0.4,vertexColors:true,transparent:true,opacity:0.9,sizeAttenuation:true,depthWrite:false}));
suctParticles.frustumCulled=false; scene.add(suctParticles);
function captureRadius(){
  // Einzugsschlauch: A = Q/u  ->  r = sqrt(A/π); u = lokale Anströmung an der Mündung
  velocityAt(X_IN-6, Y_BANK, pipeZ(), _vel);
  const u=Math.max(Math.hypot(_vel.u,_vel.v,_vel.w), 0.05);
  return Math.sqrt(D.Qpump/(Math.PI*u));
}
function seedSuct(n){
  const zp=pipeZ(), rc=captureRadius();
  sEsc[n]=X_IN + 6 + Math.random()*20;               // zufällige Auslaufgrenze stromab
  for(let tries=0;tries<14;tries++){
    // stromauf säen: 65 % im Einzugsschlauch (werden gefangen),
    // 35 % im Ring darum (ziehen sichtbar vorbei -> zeigt die Schlauchgrenze)
    const inTube=Math.random()<0.65;
    const a=Math.random()*Math.PI*2;
    const rr=inTube ? rc*Math.sqrt(Math.random())
                    : rc*(1.15+1.6*Math.random());
    const x=X_IN - (2 + Math.random()*16);           // 2–18 m stromauf
    const y=Y_BANK + Math.cos(a)*rr;
    const z=zp     + Math.sin(a)*rr*0.7;
    if(y<0.3||y>W) continue;
    if(z<bedZ(y)+0.1||z>WSURF-0.05) continue;
    sX[n]=x; sY[n]=y; sZ[n]=z; return;
  }
  sX[n]=X_IN-6; sY[n]=Y_BANK; sZ[n]=pipeZ();
}
for(let n=0;n<NS;n++) seedSuct(n);
function updateSuction(dt){
  const zp=pipeZ();
  for(let n=0;n<NS;n++){
    velocityAt(sX[n],sY[n],sZ[n],_vel);              // ECHTES Feld (physikalische Senke)
    // leichtes turbulentes Zittern
    const jx=(Math.random()-0.5)*0.06, jy=(Math.random()-0.5)*0.06, jz=(Math.random()-0.5)*0.03;
    sX[n]+=(_vel.u+jx)*dt; sY[n]+=(_vel.v+jy)*dt; sZ[n]+=(_vel.w+jz)*dt;
    const r2=(sX[n]-X_IN)**2+(sY[n]-Y_BANK)**2+(sZ[n]-zp)**2;
    // an der Mündung eingesaugt oder vorbeigezogen -> neu säen
    if(r2<0.45 || sX[n]>sEsc[n] || sX[n]<X_IN-24 ||
       sZ[n]<bedZ(sY[n]) || sZ[n]>WSURF || !isFinite(sX[n]))
      seedSuct(n);
    const o=n*3; sPos[o]=sX[n]; sPos[o+1]=sZ[n]; sPos[o+2]=sY[n];
    const c=colNorm(sampleT(sX[n],sY[n],sZ[n]));     // Farbe = lokale Wassertemperatur
    sCol[o]=c.r; sCol[o+1]=c.g; sCol[o+2]=c.b;
  }
  suctGeo.attributes.position.needsUpdate=true;
  suctGeo.attributes.color.needsUpdate=true;
}

/* ---------------------------------------------------------------------------
   Analoge Thermometer auf Ein- und Auslassrohr: rundes Zifferblatt (0–25 °C,
   Skala wie die Farbskala) mit analogem Zeiger, per Canvas-Textur je Frame
   nachgeführt. Montiert auf kurzem Standrohr über der Rohrleitung.
   --------------------------------------------------------------------------- */
function makeGauge(x,z,label){
  const c=document.createElement('canvas'); c.width=c.height=256;
  const ctx=c.getContext('2d');
  const tex=new THREE.CanvasTexture(c);
  const grp=new THREE.Group();
  // Standrohr von der Rohrleitung (Höhe WSURF+5) nach oben
  const pole=new THREE.Mesh(new THREE.CylinderGeometry(0.16,0.16,3.0,10),
    new THREE.MeshLambertMaterial({color:0x6b6f75}));
  pole.position.set(x, WSURF+5+1.5, z); grp.add(pole);
  // Zifferblatt (beidseitig sichtbar)
  const dial=new THREE.Mesh(new THREE.PlaneGeometry(4.6,4.6),
    new THREE.MeshBasicMaterial({map:tex,transparent:true,side:THREE.DoubleSide}));
  dial.position.set(x, WSURF+5+4.6, z); grp.add(dial);
  scene.add(grp);
  function draw(temp){
    ctx.clearRect(0,0,256,256);
    // Gehäuse + Blatt
    ctx.beginPath(); ctx.arc(128,128,120,0,Math.PI*2);
    ctx.fillStyle='#1b2734'; ctx.fill();
    ctx.lineWidth=8; ctx.strokeStyle='#8aa0b5'; ctx.stroke();
    ctx.beginPath(); ctx.arc(128,128,104,0,Math.PI*2);
    ctx.fillStyle='#eef3f8'; ctx.fill();
    // Skala 0..25 °C über 240° (von -210° bis +30°, d.h. links unten -> rechts unten)
    const a0=-Math.PI*7/6, a1=Math.PI/6;
    for(let t=0;t<=25;t+=2.5){
      const a=a0+(a1-a0)*t/25, big=(t%5===0);
      const r1=big?78:88, r2=96;
      ctx.beginPath();
      ctx.moveTo(128+Math.cos(a)*r1,128+Math.sin(a)*r1);
      ctx.lineTo(128+Math.cos(a)*r2,128+Math.sin(a)*r2);
      ctx.lineWidth=big?5:2.5; ctx.strokeStyle='#33475c'; ctx.stroke();
      if(big){
        ctx.fillStyle='#33475c'; ctx.font='bold 20px Segoe UI'; 
        ctx.textAlign='center'; ctx.textBaseline='middle';
        ctx.fillText(t,128+Math.cos(a)*62,128+Math.sin(a)*62);
      }
    }
    // Beschriftung
    ctx.fillStyle='#5c7186'; ctx.font='bold 17px Segoe UI'; ctx.textAlign='center';
    ctx.fillText(label,128,158); ctx.fillText('°C',128,178);
    // Zeiger (analog)
    const tv=Math.max(0,Math.min(25,temp));
    const a=a0+(a1-a0)*tv/25;
    ctx.beginPath();
    ctx.moveTo(128-Math.cos(a)*16,128-Math.sin(a)*16);
    ctx.lineTo(128+Math.cos(a)*86,128+Math.sin(a)*86);
    ctx.lineWidth=7; ctx.lineCap='round'; ctx.strokeStyle='#d0342c'; ctx.stroke();
    ctx.beginPath(); ctx.arc(128,128,9,0,Math.PI*2); ctx.fillStyle='#22303c'; ctx.fill();
    tex.needsUpdate=true;
  }
  draw(P.tIn);
  return {draw, last:1e9};
}
const gaugeIn  = makeGauge(X_IN,  -6, 'Einlass');
const gaugeOut = makeGauge(X_OUT, -6, 'Auslass');
function updateGauges(){
  // nur bei sichtbarer Änderung neu zeichnen (Canvas-Arbeit sparen)
  if(Math.abs(gaugeIn.last -tInletLocal)>0.02){ gaugeIn.draw(tInletLocal); gaugeIn.last=tInletLocal; }
  if(Math.abs(gaugeOut.last-tOutlet)   >0.02){ gaugeOut.draw(tOutlet);    gaugeOut.last=tOutlet; }
}

/* ---------------------------------------------------------------------------
   Boot mit zwei Personen: gleitet in Flussmitte mit der lokalen
   Oberflächen-Strömungsgeschwindigkeit von Modellanfang bis -ende, verschwindet
   kurz und erscheint danach wieder am Anfang -> intuitives Gefühl für v_Fluss.
   --------------------------------------------------------------------------- */
const boat=new THREE.Group();
{
  /* Realistischer Maßstab (kleines Ruderboot):
     Rumpf 4,2 m lang, 1,5 m breit, Bordwand 0,55 m; liegt 0,15 m tief im Wasser.
     Personen SITZEND: Oberkörper ~0,7 m + Kopf -> Kopfhöhe ~1,3 m über Deck. */
  const hullMat=new THREE.MeshLambertMaterial({color:0x7a4f26,side:THREE.DoubleSide});
  // Rumpf als eine Form: flaches Heck, spitz zulaufender Bug (ExtrudeGeometry)
  const sh=new THREE.Shape();
  sh.moveTo(-2.1,-0.75); sh.lineTo(1.0,-0.75); sh.lineTo(2.1,0.0);
  sh.lineTo(1.0, 0.75);  sh.lineTo(-2.1,0.75); sh.closePath();
  const hullGeo=new THREE.ExtrudeGeometry(sh,{depth:0.55,bevelEnabled:false});
  const hull=new THREE.Mesh(hullGeo,hullMat);
  hull.rotation.x=-Math.PI/2; hull.position.y=0;   // Extrusion zeigt nach oben
  boat.add(hull);
  // Innenboden (dunkler)
  const floorGeo=new THREE.ExtrudeGeometry(sh,{depth:0.03,bevelEnabled:false});
  const floor=new THREE.Mesh(floorGeo,new THREE.MeshLambertMaterial({color:0x5b3a1c}));
  floor.rotation.x=-Math.PI/2; floor.scale.set(0.86,0.86,1);
  floor.position.y=0.08; boat.add(floor);
  // zwei Sitzbänke quer
  const benchMat=new THREE.MeshLambertMaterial({color:0x9a6b39});
  for(const bx of [0.8,-1.0]){
    const bench=new THREE.Mesh(new THREE.BoxGeometry(0.5,0.07,1.34),benchMat);
    bench.position.set(bx,0.34,0); boat.add(bench);
  }
  // zwei SITZENDE Personen (Oberkörper + Kopf, realistisch klein)
  const headMat=new THREE.MeshLambertMaterial({color:0xe8b88a});
  const mkPerson=(px,jacket)=>{
    const torso=new THREE.Mesh(new THREE.CylinderGeometry(0.20,0.26,0.68,10),
      new THREE.MeshLambertMaterial({color:jacket}));
    torso.position.set(px,0.34+0.34,0); boat.add(torso);
    const head=new THREE.Mesh(new THREE.SphereGeometry(0.21,10,10),headMat);
    head.position.set(px,0.34+0.68+0.20,0); boat.add(head);
  };
  mkPerson( 0.8,0xc9a227);   // Person vorn (gelbe Jacke)
  mkPerson(-1.0,0x2e7dd1);   // Person hinten (blaue Jacke)
  scene.add(boat);
}
let boatX=-6, boatWaitUntil=0;
function updateBoat(dt){
  const now=performance.now();
  if(boatX>L+6){                       // Modellende erreicht -> kurz warten, dann neu
    if(!boatWaitUntil){ boatWaitUntil=now+3000; boat.visible=false; }
    if(now>=boatWaitUntil){ boatWaitUntil=0; boatX=-6; boat.visible=true; }
  } else {
    velocityAt(Math.max(boatX,0), W/2, WSURF-0.3, _vel);   // Oberflächen-Strömung Flussmitte
    boatX += Math.max(_vel.u,0.05)*dt;
    boat.position.set(boatX, WSURF-0.15, W/2);             // liegt 0,15 m tief im Wasser
    boat.rotation.y = 0;                                    // Bug stets stromab
  }
}

/* Geschwindigkeitsvektoren: dichtes Pfeilfeld (Schaft + Pfeilspitze).
   Länge STRENG proportional zur Geschwindigkeit (Referenz = tatsächliches Feldmaximum),
   Helligkeit ebenfalls geschwindigkeitsabhängig -> korrekte quantitative Ablesung. */
let vecGroup=new THREE.Group(); scene.add(vecGroup);
function buildVectors(){
  vecGroup.traverse(o=>{ if(o.geometry) o.geometry.dispose(); });
  vecGroup.clear();
  // 1. Durchgang: Feld abtasten und Maximalgeschwindigkeit bestimmen
  const samples=[];
  let vmax=1e-6;
  for(let i=1;i<NX;i+=4) for(let j=0;j<NY;j+=2){
    const x=(i+0.5)*dx, y=(j+0.5)*dy, z=bedZ(y)+0.6*localDepth(y);
    velocityAt(x,y,z,_vel);
    const sp=Math.hypot(_vel.u,_vel.v,_vel.w);
    if(sp<1e-4) continue;
    samples.push({x,y,z,u:_vel.u,v:_vel.v,w:_vel.w,sp});
    if(sp>vmax) vmax=sp;
  }
  // 2. Durchgang: Pfeile bauen, Länge = 10 m * (v/vmax), Helligkeit ~ v/vmax
  const pts=[], cols=[];
  const LMAX=10;
  for(const s of samples){
    const f=s.sp/vmax;                       // 0..1
    const len=LMAX*f;
    if(len<0.5) continue;                    // praktisch stehendes Wasser: kein Pfeil
    const ux=s.u/s.sp, uy=s.v/s.sp, uz=s.w/s.sp;
    const tx=s.x+ux*len, ty=s.y+uy*len, tz=s.z+uz*len;
    const br=0.35+0.65*f;                    // Helligkeit
    const r=0.62*br, g=0.86*br, b=1.00*br;   // kühles Blau, schneller = heller
    // Schaft
    pts.push(s.x,s.z,s.y, tx,tz,ty); cols.push(r,g,b, r,g,b);
    // Pfeilspitze: zwei Flügel, horizontal senkrecht zur Richtung
    let px=-uy, py=ux;
    const pl=Math.hypot(px,py)||1; px/=pl; py/=pl;
    const hb=Math.max(len*0.28,0.5);
    const bx=tx-ux*hb, by=ty-uy*hb, bz=tz-uz*hb;
    pts.push(tx,tz,ty, bx+px*hb*0.45, bz, by+py*hb*0.45); cols.push(r,g,b, r,g,b);
    pts.push(tx,tz,ty, bx-px*hb*0.45, bz, by-py*hb*0.45); cols.push(r,g,b, r,g,b);
  }
  const g2=new THREE.BufferGeometry();
  g2.setAttribute('position',new THREE.Float32BufferAttribute(pts,3));
  g2.setAttribute('color',new THREE.Float32BufferAttribute(cols,3));
  vecGroup.add(new THREE.LineSegments(g2,
    new THREE.LineBasicMaterial({vertexColors:true,transparent:true,opacity:0.95})));
}

/* ---------------------------------------------------------------------------
   Numerisches Gitter (räumliche Diskretisierung des Fluidvolumens):
   - Draufsicht: NX×NY-Zellraster auf der Wasseroberfläche
   - Querschnitt an x=xCut: Sigma-Vernetzung der Wassersäule (NY×NZ) –
     terrainfolgende Schichten über der parabolischen Sohle
   - Längsschnitt an y=yLong: Sigma-Schichten in Strömungsrichtung (NX×NZ)
   --------------------------------------------------------------------------- */
let gridGroup=new THREE.Group(); gridGroup.visible=false; scene.add(gridGroup);
function buildGrid(){
  gridGroup.traverse(o=>{ if(o.geometry) o.geometry.dispose(); });
  gridGroup.clear();
  const pts=[];
  const zTop=WSURF+0.03;
  // (1) Oberflächenraster: Zellkanten i (quer) und j (längs)
  for(let i=0;i<=NX;i++){ const x=i*dx; pts.push(x,zTop,0, x,zTop,W); }
  for(let j=0;j<=NY;j++){ const y=j*dy; pts.push(0,zTop,y, L,zTop,y); }
  // (2) Querschnitt an x=xCut: Sigma-Schichten + vertikale Zellkanten
  const xc=Math.min(Math.max(P.xCut,0),L)+0.06;
  for(let k=0;k<=NZ;k++){                       // terrainfolgende Schichtlinien
    for(let j=0;j<NY;j++){
      const y1=j*dy, y2=(j+1)*dy;
      const z1=bedZ(y1)+k/NZ*localDepth(y1), z2=bedZ(y2)+k/NZ*localDepth(y2);
      pts.push(xc,z1,y1, xc,z2,y2);
    }
  }
  for(let j=0;j<=NY;j++){                       // vertikale Kanten Sohle -> Oberfläche
    const y=j*dy; pts.push(xc,bedZ(y),y, xc,WSURF,y);
  }
  // (3) Längsschnitt an y=yLong: Sigma-Schichten + vertikale Zellkanten
  const yl=Math.min(Math.max(P.yLong,0.01),W-0.01);
  const bl=bedZ(yl), dl=localDepth(yl);
  for(let k=0;k<=NZ;k++){ const z=bl+k/NZ*dl; pts.push(0,z,yl, L,z,yl); }
  for(let i=0;i<=NX;i++){ const x=i*dx; pts.push(x,bl,yl, x,WSURF,yl); }
  const g=new THREE.BufferGeometry();
  g.setAttribute('position',new THREE.Float32BufferAttribute(pts,3));
  gridGroup.add(new THREE.LineSegments(g,
    new THREE.LineBasicMaterial({color:0x7fd4dd,transparent:true,opacity:0.35,depthWrite:false})));
}

/* ---- Konturlinien (Marching Squares) an der Wasseroberfläche ---- */
// Kantentabelle: Bit bl=1, br=2, tr=4, tl=8 ; Kanten a=unten,b=rechts,c=oben,d=links
const MS_EDGES={1:[[3,0]],2:[[0,1]],3:[[3,1]],4:[[1,2]],5:[[3,0],[1,2]],6:[[0,2]],
  7:[[3,2]],8:[[2,3]],9:[[2,0]],10:[[0,1],[2,3]],11:[[2,1]],12:[[1,3]],13:[[1,0]],14:[[0,3]]};
// liefert Liniensegmente [x1,y1,x2,y2,...] in Gitterkoordinaten (0..gx-1, 0..gy-1)
function marchingSquares(val,gx,gy,level){
  const seg=[];
  const ep=(edge,x,y,vtl,vtr,vbr,vbl)=>{ // Kantenschnittpunkt
    let t;
    if(edge===0){ t=(level-vbl)/(vbr-vbl); return [x+ (isFinite(t)?t:0.5), y+1]; }      // unten
    if(edge===1){ t=(level-vbr)/(vtr-vbr); return [x+1, y+1-(isFinite(t)?t:0.5)]; }      // rechts
    if(edge===2){ t=(level-vtr)/(vtl-vtr); return [x+1-(isFinite(t)?t:0.5), y]; }        // oben
    t=(level-vtl)/(vbl-vtl); return [x, y+(isFinite(t)?t:0.5)];                           // links
  };
  for(let y=0;y<gy-1;y++) for(let x=0;x<gx-1;x++){
    const vtl=val(x,y), vtr=val(x+1,y), vbr=val(x+1,y+1), vbl=val(x,y+1);
    let ci=0; if(vbl>level)ci|=1; if(vbr>level)ci|=2; if(vtr>level)ci|=4; if(vtl>level)ci|=8;
    const pairs=MS_EDGES[ci]; if(!pairs) continue;
    for(const pr of pairs){
      const p=ep(pr[0],x,y,vtl,vtr,vbr,vbl), q=ep(pr[1],x,y,vtl,vtr,vbr,vbl);
      seg.push(p[0],p[1],q[0],q[1]);
    }
  }
  return seg;
}
// baut LineSegments aus mehreren Niveaus eines Oberflächen-Skalarfeldes
function contourLines(getScalar,levels,colorFn,zOff){
  const GX=90, GY=34;
  const cache=[]; for(let b=0;b<GY;b++){cache[b]=[];for(let a=0;a<GX;a++)
    cache[b][a]=getScalar(a/(GX-1)*L, b/(GY-1)*W);}
  const val=(a,b)=>cache[b][a];
  const group=new THREE.Group();
  for(const lv of levels){
    const seg=marchingSquares(val,GX,GY,lv); if(!seg.length) continue;
    const pos=[];
    for(let s=0;s<seg.length;s+=4){
      pos.push(seg[s]/(GX-1)*L, zOff, seg[s+1]/(GY-1)*W,
               seg[s+2]/(GX-1)*L, zOff, seg[s+3]/(GY-1)*W);
    }
    const g=new THREE.BufferGeometry();
    g.setAttribute('position',new THREE.Float32BufferAttribute(pos,3));
    const col=colorFn(lv);
    group.add(new THREE.LineSegments(g,new THREE.LineBasicMaterial({color:col,transparent:true,opacity:0.95})));
  }
  group.frustumCulled=false;
  return group;
}

/* Isothermen (Temperatur-Konturlinien) */
let isoGroup=null;
function buildIso(){
  if(isoGroup){scene.remove(isoGroup);isoGroup.traverse(o=>o.geometry&&o.geometry.dispose());}
  // Niveaus aus dem TATSÄCHLICHEN Wertebereich des Oberflächenfeldes bestimmen
  // (nicht aus der festen Farbskala – sonst schneidet kein Niveau das Feld)
  const fs=(x,y)=>sampleT(x,y,WSURF-0.05);
  let lo=1e9, hi=-1e9;
  for(let a=0;a<=40;a++)for(let b=0;b<=16;b++){
    const t=fs(a/40*L, b/16*W); if(t<lo)lo=t; if(t>hi)hi=t;
  }
  if(hi-lo<0.05){ const c=(lo+hi)/2; lo=c-0.025; hi=c+0.025; }  // Mindestspanne
  const levels=[];
  for(let n=1;n<=7;n++) levels.push(lo+(hi-lo)*n/8);   // 7 Niveaus im Feldbereich
  const cN=new THREE.Color();
  isoGroup=contourLines(fs, levels,
    lv=>{const c=colNorm(lv); cN.setRGB(c.r,c.g,c.b); return cN.getHex();}, WSURF+0.06);
  isoGroup.visible=document.getElementById('t_iso').checked;
  scene.add(isoGroup);
}

/* Strömungs-Isolinien (Iso-Linien gleicher Geschwindigkeit) */
let flowGroup=null;
function buildFlow(){
  if(flowGroup){scene.remove(flowGroup);flowGroup.traverse(o=>o.geometry&&o.geometry.dispose());}
  // Geschwindigkeitsbetrag nahe der Oberfläche
  const speed=(x,y)=>{ const z=bedZ(y)+0.85*localDepth(y); velocityAt(x,y,z,_vel);
    return Math.hypot(_vel.u,_vel.v,_vel.w); };
  // Niveaus aus dem aktuellen Maximum
  let vmax=0.1; for(let a=0;a<12;a++)for(let b=0;b<8;b++) vmax=Math.max(vmax,speed(a/11*L,b/7*W));
  const levels=[]; for(let n=1;n<=6;n++) levels.push(vmax*n/7);
  flowGroup=contourLines(speed, levels, ()=>0x8fe3ff, WSURF+0.10);
  flowGroup.visible=document.getElementById('t_flow').checked;
  scene.add(flowGroup);
}

/* =========================================================================
   (5) UI – Regler, abgeleitete Werte, Schalter, KPIs
   ========================================================================= */
const fmt=(v,d=1)=>v.toLocaleString('de-DE',{minimumFractionDigits:d,maximumFractionDigits:d});
const UNIT={vFlow:' m/s',depth:' m',tIn:' °C',power:' kW',vSuct:' m/s',
            dPipe:' m',xCut:' m',zCut:' m',yLong:' m',cop:''};
const DEC ={vFlow:2,depth:1,tIn:1,power:0,vSuct:2,dPipe:2,xCut:0,zCut:1,yLong:1,cop:1};
const powMW = ()=>fmt(P.power/1000,2)+' MW';   // Heizleistung in Megawatt
const fmtPow = kW => (kW>=1000) ? fmt(kW/1000,2)+' MW' : fmt(kW,0)+' kW';  // kW/MW automatisch

function refreshLabels(){
  document.querySelectorAll('.ctrl[data-p]').forEach(c=>{
    const p=c.dataset.p;
    c.querySelector('.val').textContent = (p==='power') ? powMW() : fmt(P[p],DEC[p])+UNIT[p];
    c.querySelector('input').value=P[p];
  });
}
function refreshDerived(){
  document.getElementById('derivedFlow').innerHTML=
    `Reynolds-Zahl Re ≈ <b>${D.Re.toExponential(2)}</b> &nbsp;·&nbsp; turbulent<br>`+
    `Péclet-Zahl Pe ≈ <b>${D.Pe.toExponential(2)}</b><br>`+
    `eff. Temp.-Leitf. α ≈ <b>${D.alpha.toExponential(2)}</b> m²/s`;
  document.getElementById('derivedPump').innerHTML=
    `Heizleistung P<sub>heiz</sub> = <b>${powMW()}</b> · COP = <b>${fmt(P.cop,1)}</b><br>`+
    `elektr. Leistung P<sub>el</sub> = P<sub>heiz</sub>/COP = <b>${fmtPow(D.Pel)}</b><br>`+
    `Entzug aus Fluss P<sub>entzug</sub> = P<sub>heiz</sub>·(1−1/COP) = <b>${fmtPow(D.Pextract)}</b><br>`+
    `Rohrquerschnitt A = <b>${fmt(D.area,3)}</b> m²<br>`+
    `Durchfluss Q = <b>${fmt(D.Qpump,3)}</b> m³/s = <b>${fmt(D.Qpump*1000,0)}</b> l/s<br>`+
    `Temp.-Absenkung ΔT = P<sub>entzug</sub>/(ρ·c<sub>p</sub>·Q) = <b>${fmt(D.dT,2)}</b> K<br>`+
    `Auslasstemperatur ≈ <b>${fmt(tOutlet,2)}</b> °C`;
}

/* Slider-Verdrahtung */
document.querySelectorAll('.ctrl[data-p] input').forEach(inp=>{
  const p=inp.closest('.ctrl').dataset.p;
  inp.addEventListener('input',()=>{
    P[p]=parseFloat(inp.value);
    onParamChange(p);
  });
});
/* ---------- Gitterauflösung umschalten (Regler "Numerik") ------------------ */
const GRID_NAMES=['Grob (0,75×)','Standard (1×)','Fein (1,5×)','Sehr fein (2×)'];
function updateGridInfo(lv){
  document.getElementById('gridVal').textContent=GRID_NAMES[lv-1];
  document.getElementById('gridInfo').innerHTML=
    `Zellen: <b>${NX} × ${NY} × ${NZ}</b> = <b>${NC.toLocaleString('de-DE')}</b> · `+
    `Δx = <b>${fmt(dx,2)}</b> m · Δy = <b>${fmt(dy,2)}</b> m · Δz = h(y)/${NZ}<br>`+
    `Mischung über physikalische Diffusivitäten (D<sub>y</sub> = 0,6·h·u*, `+
    `D<sub>z</sub> = 0,067·h·u*) – Ergebnisse damit gitterunabhängig.`;
}
function setGridLevel(lv){
  const s=GRID_S[lv-1];
  NX=Math.round(NX0*s); NY=Math.round(NY0*s); NZ=Math.round(NZ0*Math.min(s,2));
  dx=L/NX; dy=W/NY; NC=NX*NY*NZ;
  T=new Float32Array(NC).fill(P.tIn);
  T2=new Float32Array(NC);
  updateVCOL();
  relax(60);                              // neues Feld einschwingen
  fieldStats(); buildVectors(); buildFlow();
  if(isoGroup && document.getElementById('t_iso').checked) buildIso();
  if(gridGroup.visible || gridGroup.children.length>0) buildGrid();
  updateGridInfo(lv);
}
document.getElementById('gridRes').addEventListener('input',e=>setGridLevel(+e.target.value));

function onParamChange(p){
  computeDerived();
  if(p==='depth'){ updateVCOL(); buildEnv(); buildWater(); buildPump(); }
  if(p==='dPipe'){ buildPump(); }
  placePlanes();
  refreshLabels(); refreshDerived();
  // Felder Richtung neues Gleichgewicht iterieren (Rest erledigt die Schleife)
  relax(16);
  fieldStats(); buildVectors(); buildFlow();
  if(gridGroup.visible) buildGrid();   // Gitter folgt Tiefe/Schnittpositionen
}

/* Schalter (Toggles) */
const TG={t_cross:planeCross.mesh,t_long:planeLong.mesh,t_horiz:planeHoriz.mesh,
          t_stream:particles,t_vec:vecGroup};
Object.entries(TG).forEach(([id,obj])=>{
  const el=document.getElementById(id);
  obj.visible=el.checked;
  el.addEventListener('change',()=>{obj.visible=el.checked;});
});
// Kreislauf-Partikel folgen dem Stromlinien-Schalter
document.getElementById('t_stream').addEventListener('change',e=>{circParticles.visible=e.target.checked; jetParticles.visible=e.target.checked; suctParticles.visible=e.target.checked;});
circParticles.visible=document.getElementById('t_stream').checked;
jetParticles.visible=document.getElementById('t_stream').checked;
suctParticles.visible=document.getElementById('t_stream').checked;
document.getElementById('t_water').addEventListener('change',e=>{if(waterMesh)waterMesh.visible=e.target.checked;});
document.getElementById('t_iso').addEventListener('change',e=>{ if(!isoGroup)buildIso(); isoGroup.visible=e.target.checked; });
document.getElementById('t_grid').addEventListener('change',e=>{ if(e.target.checked && gridGroup.children.length===0) buildGrid(); gridGroup.visible=e.target.checked; });
document.getElementById('t_flow').addEventListener('change',e=>{ if(!flowGroup)buildFlow(); flowGroup.visible=e.target.checked; });

/* Ansichten */
function setView(v){
  const cx=L/2, cz=W/2;
  // Abstand so wählen, dass die gesamte Modellbreite/-länge sicher ins Bild passt
  const fov=camera.fov*Math.PI/180;
  const fit=(L*0.6)/Math.tan(fov/2)/Math.max(camera.aspect,1) + 120;
  const cy=WSURF*0.5;
  if(v==='top')        camera.position.set(cx, fit*0.9, cz+0.1);
  else if(v==='side')  camera.position.set(cx, 30, cz + fit);
  else /*oblique/reset*/camera.position.set(cx - fit*0.55, fit*0.5, cz + fit*0.7);
  controls.target.set(cx,cy,cz);
  camera.lookAt(cx,cy,cz);
  controls.update();
}
document.querySelectorAll('.views button').forEach(b=>
  b.addEventListener('click',()=>setView(b.dataset.view)));

/* KPI-Panel */
function updateKPI(){
  // maximale Temperaturabsenkung gegenüber Einström-/Umgebungsreferenz
  const ref=P.tIn;
  let maxDrop=0, recX=X_OUT;
  // Max. lokale Absenkung: über das GESAMTE 3D-Feld (das Defizit sitzt sohlnah!)
  for(let k=0;k<NZ;k++) for(let j=0;j<NY;j++) for(let i=0;i<NX;i++){
    const drop=ref-T[idx(i,j,k)]; if(drop>maxDrop)maxDrop=drop;
  }
  // Erholungslänge: ab Auslass stromab, bis die Fahne im GESAMTEN Querschnitt < 0.05 K
  const iO=Math.round(X_OUT/dx);
  recX=L;
  for(let i=iO;i<NX;i++){
    let mx=0;
    for(let k=0;k<NZ;k++) for(let j=0;j<NY;j++){
      const d=ref-T[idx(i,j,k)]; if(d>mx)mx=d;
    }
    if(mx<0.05){recX=(i+0.5)*dx;break;}
  }
  const recLen=Math.max(recX-X_OUT,0);
  // Energiebilanz: Flussdurchfluss (parabolischer Querschnitt) und durchmischte Abkühlung
  const Ariver = (2/3)*W*P.depth;                 // Fläche der Parabel
  const Qriver = Ariver * P.vFlow;                 // Volumenstrom Fluss [m³/s]
  const dTmix  = (D.Pextract*1000)/(RHO*CP*Math.max(Qriver,1e-6)); // P_entzug/(ρ·c_p·Q_Fluss)
  const rows=[
    ['Thermische Heizleistung', powMW()],
    ['Leistungszahl (COP)', fmt(P.cop,1)],
    ['Elektrische Leistung', fmtPow(D.Pel)],
    ['Entzugsleistung aus Fluss', fmtPow(D.Pextract)],
    ['Durchfluss Pumpe', fmt(D.Qpump,3)+' m³/s'],
    ['Durchfluss Fluss', fmt(Qriver,1)+' m³/s'],
    ['ΔT Rückgabe (Rohr)', fmt(D.dT,2)+' K'],
    ['Auslasstemperatur', fmt(tOutlet,2)+' °C'],
    ['Abkühlung stromab (durchmischt)', fmt(dTmix,3)+' K → '+fmt(P.tIn-dTmix,2)+' °C'],
    ['Max. lokale Absenkung (Fahne)', fmt(maxDrop,2)+' K'],
    ['Reichweite Kaltwasserfahne', recLen>=L-X_OUT?'> '+fmt(L-X_OUT,0)+' m':fmt(recLen,0)+' m'],
    ['Reynolds-Zahl', D.Re.toExponential(2)],
    ['Péclet-Zahl', D.Pe.toExponential(2)],
  ];
  document.getElementById('kpiBody').innerHTML=rows.map(r=>
    `<tr><td class="k">${r[0]}</td><td class="v">${r[1]}</td></tr>`).join('');
  // Warnungen
  const w=document.getElementById('warn'); let msg='';
  if(D.dT>10) msg='⚠ ΔT sehr groß ('+fmt(D.dT,1)+' K): Durchfluss zu klein oder Leistung zu hoch.';
  if(tInletLocal-D.dT<0) msg='⚠ Rechnerische Auslasstemperatur < 0 °C – auf 0 °C begrenzt (Eisgefahr).';
  w.className=msg?'show':''; w.textContent=msg;
}

/* Colorbar zeichnen */
const cbar=document.getElementById('cbar'), cbx=cbar.getContext('2d');
function drawColorbar(){
  for(let y=0;y<cbar.height;y++){
    const t=1-y/cbar.height; const c=tempColor(t);
    cbx.fillStyle=`rgb(${c.r*255|0},${c.g*255|0},${c.b*255|0})`;
    cbx.fillRect(0,y,cbar.width,1);
  }
  const ticks=document.getElementById('cbarTicks'); ticks.innerHTML='';
  for(let i=0;i<5;i++){
    const t=Tmax-(Tmax-Tmin)*i/4;
    const s=document.createElement('span'); s.textContent=fmt(t,1); ticks.appendChild(s);
  }
}

/* =========================================================================
   (6) DIAGRAMME (2D-Canvas, live)
   ========================================================================= */
function lineChart(canvas, series, xlabel, ylabel, xRange, yRange){
  const ctx=canvas.getContext('2d');
  const W0=canvas.clientWidth, H0=canvas.clientHeight;
  if(canvas.width!==W0)canvas.width=W0; if(canvas.height!==H0)canvas.height=H0;
  ctx.clearRect(0,0,W0,H0);
  const m={l:38,r:8,t:6,b:20};
  const pw=W0-m.l-m.r, ph=H0-m.t-m.b;
  // Achsen
  ctx.strokeStyle='#2a3a49'; ctx.lineWidth=1;
  ctx.strokeRect(m.l,m.t,pw,ph);
  ctx.fillStyle='#6b8096'; ctx.font='9px Segoe UI'; ctx.textAlign='center';
  // y-Gitter+Beschriftung
  for(let g=0;g<=4;g++){
    const yv=yRange[0]+(yRange[1]-yRange[0])*g/4;
    const py=m.t+ph-ph*g/4;
    ctx.strokeStyle='#1b2733'; ctx.beginPath();ctx.moveTo(m.l,py);ctx.lineTo(m.l+pw,py);ctx.stroke();
    ctx.textAlign='right'; ctx.fillText(yv.toFixed(yRange[1]-yRange[0]<3?2:1),m.l-4,py+3);
  }
  // x-Beschriftung
  ctx.textAlign='center';
  for(let g=0;g<=4;g++){
    const xv=xRange[0]+(xRange[1]-xRange[0])*g/4;
    const px=m.l+pw*g/4; ctx.fillText(xv.toFixed(0),px,H0-7);
  }
  ctx.fillStyle='#52677a'; ctx.fillText(xlabel,m.l+pw/2,H0-0.5);
  // Datenreihen
  const sx=v=>m.l+pw*(v-xRange[0])/(xRange[1]-xRange[0]);
  const sy=v=>m.t+ph-ph*(v-yRange[0])/(yRange[1]-yRange[0]);
  series.forEach(s=>{
    ctx.strokeStyle=s.color; ctx.lineWidth=1.6; ctx.beginPath();
    s.data.forEach((d,n)=>{const X=sx(d[0]),Y=sy(d[1]); n?ctx.lineTo(X,Y):ctx.moveTo(X,Y);});
    ctx.stroke();
  });
  // Legende
  if(series.length>1){
    ctx.textAlign='left'; let ly=m.t+10;
    series.forEach(s=>{ctx.fillStyle=s.color;ctx.fillText('— '+s.name,m.l+6,ly);ly+=11;});
  }
}
const chTx=document.getElementById('chTx'), chTy=document.getElementById('chTy'),
      chVz=document.getElementById('chVz'), chTz=document.getElementById('chTz');
function updateCharts(){
  const jMid=Math.round(NY/2), jBank=Math.round(Y_BANK/dy);
  const iC=Math.min(Math.max(Math.round(P.xCut/dx),0),NX-1);
  // Drei Vertikal-Auswertungen je Säule: Oberfläche (k=NZ-1), Tiefenmittel, Sohle (k=0).
  // Grund: Das Kaltwasser liegt nach dem Prallstrahl sohlnah – die drei Kurven zeigen
  // die vertikale Schichtung der Fahne direkt im 2D-Diagramm.
  const colAvg=(i,j)=>{let s=0; for(let k=0;k<NZ;k++) s+=T[idx(i,j,k)]; return s/NZ;};
  const CS={surf:'#7fd0ff', mean:'#ffb454', bed:'#ff5c5c'};
  // T(x) am Pumpenufer: Oberfläche / Tiefenmittel / Sohle
  const bxS = [];
  const bxM = [];
  const bxB = [];
  for(let i=0;i<NX;i++){const x=(i+0.5)*dx;
    bxS.push([x,T[idx(i,jBank,NZ-1)]]);
    bxM.push([x,colAvg(i,jBank)]);
    bxB.push([x,T[idx(i,jBank,0)]]);}
  let lo=Tmin-0.1, hi=Tmax+0.1;
  const TLO=0, THI=30;                 // feste Temperaturskala der Diagramme [°C]
  lineChart(chTx,[{name:'Oberfläche',color:CS.surf,data:bxS},
                  {name:'Tiefenmittel',color:CS.mean,data:bxM},
                  {name:'Sohle',color:CS.bed,data:bxB}],
            'x [m]','T',[0,L],[TLO,THI]);
  // T(y) am Querschnitt: Oberfläche / Tiefenmittel / Sohle
  const tyS = [];
  const tyM = [];
  const tyB = [];
  for(let j=0;j<NY;j++){const y=(j+0.5)*dy;
    tyS.push([y,T[idx(iC,j,NZ-1)]]);
    tyM.push([y,colAvg(iC,j)]);
    tyB.push([y,T[idx(iC,j,0)]]);}
  lineChart(chTy,[{name:'Oberfläche',color:CS.surf,data:tyS},
                  {name:'Tiefenmittel',color:CS.mean,data:tyM},
                  {name:'Sohle',color:CS.bed,data:tyB}],
            'y [m]','T',[0,W],[TLO,THI]);
  // v(z) über die Tiefe an x=xCut, Flussmitte (tiefste Stelle); Höhe z absolut (Oberfläche fest bei WSURF)
  const vz=[]; let vmax=0.1; const dlC=localDepth(W/2), bC=bedZ(W/2);
  for(let k=0;k<NZ;k++){ const z=bC+(k+0.5)/NZ*dlC; velocityAt(P.xCut,W/2,z,_vel);
    const sp=Math.hypot(_vel.u,_vel.v,_vel.w); vmax=Math.max(vmax,sp); vz.push([z,sp]); }
  lineChart(chVz,[{name:'v(z)',color:'#9be7ff',data:vz}],'Höhe z [m]','v',[bC,WSURF],[0,vmax*1.1]);
  // T(Tiefe) vertikales Profil am Querschnitt, Pumpenufer (lokale Tiefe)
  const tz=[]; const dlB=localDepth((jBank+0.5)*dy), bB=bedZ((jBank+0.5)*dy);
  for(let k=0;k<NZ;k++){ const z=bB+(k+0.5)/NZ*dlB; tz.push([z,T[idx(iC,jBank,k)]]); }
  lineChart(chTz,[{name:'T(z)',color:'#ff8fab',data:tz}],'Höhe z [m]','T',[bB,WSURF],[TLO,THI]);
  // aktuelle Daten für den CSV-Export registrieren
  EXPORTS.chTx={file:'T_x_Pumpenufer', head:['x [m]','T Oberfläche [°C]','T Tiefenmittel [°C]','T Sohle [°C]'],
                rows:bxM.map((p,i)=>[p[0],bxS[i][1],p[1],bxB[i][1]])};
  EXPORTS.chTy={file:'T_y_Querschnitt_x'+Math.round(P.xCut), head:['y [m]','T Oberfläche [°C]','T Tiefenmittel [°C]','T Sohle [°C]'],
                rows:tyM.map((p,j)=>[p[0],tyS[j][1],p[1],tyB[j][1]])};
  EXPORTS.chVz={file:'v_z_Profil_x'+Math.round(P.xCut), head:['Höhe z [m]','v [m/s]'], rows:vz};
  EXPORTS.chTz={file:'T_z_Profil_x'+Math.round(P.xCut), head:['Höhe z [m]','T [°C]'], rows:tz};
  // Spezifikationen für die SVG-Diagramme des Reports (Vektorgrafik)
  REPCHART.tx={series:[{name:'Oberfläche',color:'#2b7bb9',data:bxS},
                       {name:'Tiefenmittel',color:'#c77b1e',data:bxM},
                       {name:'Sohle',color:'#c0392b',data:bxB}],
               xl:'x [m]',yl:'T [°C]',xr:[0,L],yr:[TLO,THI]};
  REPCHART.ty={series:[{name:'Oberfläche',color:'#2b7bb9',data:tyS},
                       {name:'Tiefenmittel',color:'#c77b1e',data:tyM},
                       {name:'Sohle',color:'#c0392b',data:tyB}],
               xl:'y [m]',yl:'T [°C]',xr:[0,W],yr:[TLO,THI]};
  REPCHART.vz={series:[{name:'v(z)',color:'#2b7bb9',data:vz}],xl:'Höhe z [m]',yl:'v [m/s]',xr:[bC,WSURF],yr:[0,vmax*1.1]};
  REPCHART.tz={series:[{name:'T(z)',color:'#c2477d',data:tz}],xl:'Höhe z [m]',yl:'T [°C]',xr:[bB,WSURF],yr:[TLO,THI]};
}
const REPCHART={};

/* CSV-Export der Diagrammdaten (Semikolon-getrennt, Dezimalkomma – Excel/DE) */
const EXPORTS={};
function exportCSV(id){
  const e=EXPORTS[id]; if(!e) return;
  const num=v=>(Math.round(v*10000)/10000).toString().replace('.',',');
  const lines=[e.head.join(';')];
  for(const r of e.rows) lines.push(r.map(num).join(';'));
  const blob=new Blob(['\ufeff'+lines.join('\r\n')],{type:'text/csv;charset=utf-8'});
  const a=document.createElement('a');
  a.href=URL.createObjectURL(blob);
  a.download=e.file+'.csv';
  document.body.appendChild(a); a.click();
  setTimeout(()=>{URL.revokeObjectURL(a.href); a.remove();},0);
}
document.querySelectorAll('.chart .exp').forEach(btn=>{
  btn.addEventListener('click',()=>exportCSV(btn.dataset.c));
});

/* =========================================================================
   (6b) REPORT: vollständiger Bericht als zweiter Reiter, PDF via Druckdialog
   ========================================================================= */
const elReport=document.getElementById('report'), elScene3d=document.getElementById('scene');
const tab3d=document.getElementById('tab3d'), tabReport=document.getElementById('tabReport');
function showTab(which){
  const rep=(which==='report');
  elReport.hidden=!rep;
  tab3d.classList.toggle('on',!rep); tabReport.classList.toggle('on',rep);
}
tab3d.addEventListener('click',()=>showTab('3d'));
tabReport.addEventListener('click',()=>showTab('report'));

/* 3D-Schnappschuss: Kamera positionieren, rendern, in NATIVER Renderer-Auflösung erfassen */
function snap3d(camPos,target,w){
  const oldPos=camera.position.clone(), oldTgt=controls.target.clone();
  camera.position.copy(camPos); controls.target.copy(target);
  camera.lookAt(target); controls.update();
  renderer.render(scene,camera);
  const src=renderer.domElement;
  const cw=Math.min(w||src.width, src.width), chh=Math.round(cw*src.height/src.width);
  const c=document.createElement('canvas'); c.width=cw; c.height=chh;
  c.getContext('2d').drawImage(src,0,0,cw,chh);
  const url=c.toDataURL('image/jpeg',0.88);
  camera.position.copy(oldPos); controls.target.copy(oldTgt);
  camera.lookAt(oldTgt); controls.update();
  return url;
}

/* SVG-Liniendiagramm (Vektorgrafik für den Report – gestochen scharf in Ansicht und PDF) */
function svgChart(spec){
  const W0=860,H0=300,m={l:56,r:16,t:14,b:42};
  const pw=W0-m.l-m.r, ph=H0-m.t-m.b;
  const [x0,x1]=spec.xr,[y0,y1]=spec.yr;
  const X=v=>m.l+(v-x0)/(x1-x0)*pw, Y=v=>m.t+ph-(v-y0)/(y1-y0)*ph;
  const nf=v=>{const r=Math.round(v*100)/100;return r.toString().replace('.',',');};
  let s=`<svg class="rsvg" viewBox="0 0 ${W0} ${H0}" xmlns="http://www.w3.org/2000/svg" font-family="Segoe UI,Arial,sans-serif">`;
  s+=`<rect x="${m.l}" y="${m.t}" width="${pw}" height="${ph}" fill="#fdfefe" stroke="#8aa0b5" stroke-width="1.2"/>`;
  for(let g=0;g<=4;g++){
    const gy=m.t+ph*g/4, vy=y1-(y1-y0)*g/4;
    const gx=m.l+pw*g/4, vx=x0+(x1-x0)*g/4;
    if(g>0&&g<4){
      s+=`<line x1="${m.l}" y1="${gy}" x2="${m.l+pw}" y2="${gy}" stroke="#dde5ec"/>`;
      s+=`<line x1="${gx}" y1="${m.t}" x2="${gx}" y2="${m.t+ph}" stroke="#dde5ec"/>`;
    }
    s+=`<text x="${m.l-7}" y="${gy+4}" text-anchor="end" font-size="12" fill="#33475c">${nf(vy)}</text>`;
    s+=`<text x="${gx}" y="${m.t+ph+17}" text-anchor="middle" font-size="12" fill="#33475c">${nf(vx)}</text>`;
  }
  s+=`<text x="${m.l+pw/2}" y="${H0-8}" text-anchor="middle" font-size="12.5" fill="#33475c">${spec.xl}</text>`;
  s+=`<text x="15" y="${m.t+ph/2}" text-anchor="middle" font-size="12.5" fill="#33475c" transform="rotate(-90 15 ${m.t+ph/2})">${spec.yl}</text>`;
  spec.series.forEach((se,si)=>{
    const ptsS=se.data.map(p=>`${X(p[0]).toFixed(1)},${Y(Math.max(y0,Math.min(y1,p[1]))).toFixed(1)}`).join(' ');
    s+=`<polyline points="${ptsS}" fill="none" stroke="${se.color}" stroke-width="2.2" stroke-linejoin="round"/>`;
    if(spec.series.length>1){
      const lx=m.l+12, ly=m.t+16+si*18;
      s+=`<line x1="${lx}" y1="${ly-4}" x2="${lx+22}" y2="${ly-4}" stroke="${se.color}" stroke-width="2.5"/>`;
      s+=`<text x="${lx+28}" y="${ly}" font-size="12.5" fill="#33475c">${se.name}</text>`;
    }
  });
  return s+'</svg>';
}
/* Drehansicht: Frames rund um das Modell (2 Zoomstufen) für interaktives Drehen/Zoomen */
let TT={frames:[],rings:2,perRing:28,ring:0,idx:0};
function captureTurntable(){
  TT.frames=[]; 
  const cx=L/2, cy=WSURF*0.5, cz=W/2;
  const fov=camera.fov*Math.PI/180;
  const base=(L*0.6)/Math.tan(fov/2)/Math.max(camera.aspect,1)+120;
  const radii=[base*0.92, base*0.55], heights=[base*0.5, base*0.3];
  for(let r=0;r<TT.rings;r++){
    for(let i=0;i<TT.perRing;i++){
      const a=i/TT.perRing*Math.PI*2;
      const pos=new THREE.Vector3(cx+Math.cos(a)*radii[r], heights[r], cz+Math.sin(a)*radii[r]);
      TT.frames.push(snap3d(pos,new THREE.Vector3(cx,cy,cz),1400));
    }
  }
}
function ttSrc(){ return TT.frames[TT.ring*TT.perRing+TT.idx]; }

const RLBL={vFlow:'Strömungsgeschwindigkeit',depth:'Max. Wassertiefe (Flussmitte)',tIn:'Flusswassertemperatur',
  power:'Thermische Heizleistung',cop:'Leistungszahl (COP)',vSuct:'Ansauggeschwindigkeit (Rohr)',
  dPipe:'Rohrdurchmesser (Ein-/Auslass)',xCut:'Querschnitt x-Position',zCut:'Horizontalschnitt-Tiefe',yLong:'Längsschnitt y-Position'};
function rRow(k){ return `<tr><td>${RLBL[k]}</td><td>${fmt(P[k],DEC[k])}${UNIT[k]}</td></tr>`; }

function generateReport(){
  // Helfer für LaTeX: deutsches Dezimalkomma ohne Abstand ({,}), Leistung mit Einheit
  const mnum=s=>String(s).replace(/,/g,'{,}');
  const mpow=kW=>(kW>=1000)?mnum(fmt(kW/1000,2))+'\\ \\mathrm{MW}':mnum(fmt(kW,0))+'\\ \\mathrm{kW}';
  // aktuelle Anzeigen einsammeln
  const kpiHTML=document.getElementById('kpiBody').innerHTML;
  const warnHTML=document.getElementById('warn').innerHTML;
  // 3D-Ansichten (Kamera wird danach wiederhergestellt)
  const cx=L/2, cy=WSURF*0.5, cz=W/2;
  const fov=camera.fov*Math.PI/180;
  const fit=(L*0.6)/Math.tan(fov/2)/Math.max(camera.aspect,1)+120;
  const tgt=new THREE.Vector3(cx,cy,cz);
  const vOblique=snap3d(new THREE.Vector3(cx-fit*0.55,fit*0.5,cz+fit*0.7),tgt);
  const vTop    =snap3d(new THREE.Vector3(cx,fit*0.9,cz+0.1),tgt);
  const vSide   =snap3d(new THREE.Vector3(cx,30,cz+fit),tgt);
  const vOutlet =snap3d(new THREE.Vector3(X_OUT-45,38,Y_BANK+62),new THREE.Vector3(X_OUT+18,WSURF-2,Y_BANK+8));
  // Gitteransicht: Gitter temporär einblenden
  const gridWas=gridGroup.visible;
  if(gridGroup.children.length===0) buildGrid();
  gridGroup.visible=true;
  const vGrid=snap3d(new THREE.Vector3(P.xCut-60,34,W+58),new THREE.Vector3(P.xCut+10,WSURF-3,W/2));
  gridGroup.visible=gridWas;
  // Drehansicht erfassen
  captureTurntable();
  renderer.render(scene,camera);   // Live-Bild wiederherstellen

  const Qriver=(2/3)*W*P.depth*P.vFlow;
  const dTmix=(D.Pextract*1000)/(RHO*CP*Math.max(Qriver,1e-6));
  const now=new Date();
  const dstr=now.toLocaleDateString('de-DE',{day:'2-digit',month:'2-digit',year:'numeric'})+', '+
             now.toLocaleTimeString('de-DE',{hour:'2-digit',minute:'2-digit'})+' Uhr';

  const secs=[['s1','1. Zusammenfassung'],['s2','2. Eingangsparameter'],['s3','3. Thermodynamik der Wärmepumpe'],
    ['s4','4. Hydraulik & dimensionslose Kennzahlen'],['s5','5. Numerisches Modell & räumliche Diskretisierung'],
    ['s6','6. 3D-Ansichten der Berechnung'],['s7','7. Interaktive 3D-Ansicht'],
    ['s8','8. Ergebnisdiagramme'],['s9','9. Kennzahlen (vollständig)'],['s10','10. Hinweise & Modellgrenzen'],
    ['s11','11. Mathematische Modellbeschreibung']];

  document.getElementById('reportBody').innerHTML=`
    <h1>Flusswasserwärmepumpe – Simulationsbericht</h1>
    <div class="rsub">3D-CFD-Simulation (vereinfachtes hydraulisch-thermisches Advektions-Diffusions-Modell) · erstellt am ${dstr}</div>

    <h2>Inhaltsverzeichnis</h2>
    <div class="toc">${secs.map(s=>`<a data-t="${s[0]}" href="#${s[0]}">${s[1]}</a>`).join('')}</div>

    <div class="sec"><h2 id="s1">1. Zusammenfassung</h2>
    <p>Untersucht wird die Rückwirkung einer Flusswasserwärmepumpe (thermische Heizleistung
    <b>${powMW()}</b>, COP <b>${fmt(P.cop,1)}</b>) auf einen ${L} m langen, ${W} m breiten Flussabschnitt
    (max. Tiefe ${fmt(P.depth,1)} m, mittlere Strömungsgeschwindigkeit ${fmt(P.vFlow,2)} m/s,
    Wassertemperatur ${fmt(P.tIn,1)} °C). Dem Fluss werden <b>${fmtPow(D.Pextract)}</b> entzogen;
    das Rückgabewasser ist um <b>${fmt(D.dT,2)} K</b> kälter (Auslasstemperatur ${fmt(tOutlet,2)} °C).
    Vollständig durchmischt entspricht dies einer Abkühlung des Flusses um <b>${fmt(dTmix,3)} K</b>.</p></div>

    <div class="sec"><h2 id="s2">2. Eingangsparameter</h2>
    <h3>Strömung &amp; Fluss</h3>
    <table class="rt">${['vFlow','depth','tIn'].map(rRow).join('')}</table>
    <h3>Wärmepumpe</h3>
    <table class="rt">${['power','cop','vSuct','dPipe'].map(rRow).join('')}</table>
    <h3>Schnittebenen (Auswertepositionen)</h3>
    <table class="rt">${['xCut','yLong','zCut'].map(rRow).join('')}</table>
    <table class="rt"><tr><td>Position Einlass (Ansaugung)</td><td>x = ${X_IN} m</td></tr>
    <tr><td>Position Auslass (Rückgabe)</td><td>x = ${X_OUT} m</td></tr>
    <tr><td>Abstand Rohrmündungen vom Ufer</td><td>y = ${Y_BANK} m</td></tr></table></div>

    <div class="sec"><h2 id="s3">3. Thermodynamik der Wärmepumpe</h2>
    <p>Die thermische Heizleistung P<sub>heiz</sub> (Nutzwärme für das Fernwärmenetz) setzt sich aus der
    dem Fluss entzogenen Wärme und der elektrischen Kompressorleistung zusammen, gekoppelt über die
    Leistungszahl (COP):</p>
    <div class="frm">$$P_{\\mathrm{el}}=\\frac{P_{\\mathrm{heiz}}}{\\mathrm{COP}}
    =\\frac{${mnum(fmt(P.power/1000,2))}\\ \\mathrm{MW}}{${mnum(fmt(P.cop,1))}}=${mpow(D.Pel)}$$
    $$P_{\\mathrm{entzug}}=P_{\\mathrm{heiz}}\\cdot\\left(1-\\frac{1}{\\mathrm{COP}}\\right)=${mpow(D.Pextract)}$$
    $$\\text{Kontrolle:}\\quad P_{\\mathrm{el}}+P_{\\mathrm{entzug}}=P_{\\mathrm{heiz}}\\;\\checkmark$$</div>
    <p>Die Temperaturabsenkung des durch die Pumpe geförderten Wassers folgt aus der Energiebilanz:</p>
    <div class="frm">$$\\Delta T=\\frac{P_{\\mathrm{entzug}}}{\\rho\\; c_p\\; Q}
    =\\frac{${mnum(fmt(D.Pextract,0))}\\ \\mathrm{kW}}
    {${RHO}\\ \\mathrm{kg/m^{3}}\\;\\cdot\\;${CP}\\ \\mathrm{J/(kg\\cdot K)}\\;\\cdot\\;${mnum(fmt(D.Qpump,3))}\\ \\mathrm{m^{3}/s}}
    =${mnum(fmt(D.dT,2))}\\ \\mathrm{K}$$</div>
    <table class="rt">
    <tr><td>Rohrquerschnitt A = π·(d/2)²</td><td>${fmt(D.area,3)} m²</td></tr>
    <tr><td>Durchfluss Pumpe Q = A·v<sub>Ansaug</sub></td><td>${fmt(D.Qpump,3)} m³/s = ${fmt(D.Qpump*1000,0)} l/s</td></tr>
    <tr><td>Auslasstemperatur T<sub>aus</sub> = T<sub>ein,lokal</sub> − ΔT</td><td>${fmt(tOutlet,2)} °C</td></tr>
    <tr><td>Wärmetauscher-Übertragungswirkungsgrad (Annahme)</td><td>1,0</td></tr></table></div>

    <div class="sec"><h2 id="s4">4. Hydraulik &amp; dimensionslose Kennzahlen</h2>
    <table class="rt">
    <tr><td>Durchfluss Fluss Q<sub>Fluss</sub> = (2/3)·B·h<sub>max</sub>·v (Parabelprofil)</td><td>${fmt(Qriver,1)} m³/s</td></tr>
    <tr><td>Reynolds-Zahl Re = v·h/ν</td><td>${D.Re.toExponential(2)} (turbulent)</td></tr>
    <tr><td>Péclet-Zahl Pe = v·h/α</td><td>${D.Pe.toExponential(2)}</td></tr>
    <tr><td>eff. Temperaturleitfähigkeit α (turbulent, empirisch)</td><td>${D.alpha.toExponential(2)} m²/s</td></tr>
    <tr><td>Abkühlung stromab (vollständig durchmischt) = P<sub>entzug</sub>/(ρ·c<sub>p</sub>·Q<sub>Fluss</sub>)</td><td>${fmt(dTmix,3)} K</td></tr></table></div>

    <div class="sec"><h2 id="s5">5. Numerisches Modell &amp; räumliche Diskretisierung</h2>
    <p>Gelöst wird eine stationär iterierte Advektions-Diffusions-Gleichung für die Wassertemperatur
    auf einem strukturierten Finite-Volumen-Gitter mit terrainfolgenden (Sigma-)Koordinaten in der
    Vertikalen. Die Flusssohle ist parabelförmig (tiefste Stelle in der Flussmitte, Resttiefe
    ${MINDEPTH} m an den Ufern); die Wasseroberfläche liegt fest bei z = ${WSURF} m.</p>
    <table class="rt">
    <tr><th>Größe</th><th>Wert</th></tr>
    <tr><td>Gebietsgröße L × B × h<sub>max</sub></td><td>${L} m × ${W} m × ${fmt(P.depth,1)} m</td></tr>
    <tr><td>Zellenzahl NX × NY × NZ</td><td>${NX} × ${NY} × ${NZ} = ${(NX*NY*NZ).toLocaleString('de-DE')} Zellen</td></tr>
    <tr><td>Zellweite Δx (längs)</td><td>${fmt(dx,2)} m</td></tr>
    <tr><td>Zellweite Δy (quer)</td><td>${fmt(dy,2)} m</td></tr>
    <tr><td>Vertikal: Sigma-Schichten</td><td>${NZ} Schichten, lokale Dicke = h(y)/${NZ} (terrainfolgend)</td></tr>
    <tr><td>Zeitschritt (semi-Lagrange, CFL-frei)</td><td>Δt = 1,2·Δx / max(v, 0,1) = ${fmt(1.2*dx/Math.max(P.vFlow,0.1),1)} s</td></tr>
    <tr><td>Diffusion (konservative Fluss-Form)</td><td>aus physikalischen Diffusivitäten: D<sub>y</sub> = 0,6·h·u*, D<sub>z</sub> = 0,067·h·u*, D<sub>x</sub> = 0,3·D<sub>y</sub> (u* ≈ 0,1·v<sub>m</sub>) – gitterunabhängig</td></tr>
    <tr><td>Geschwindigkeitsprofil</td><td>quer: parabolisch (1 − 0,5·ŷ²) · vertikal: 1/7-Potenzgesetz</td></tr>
    <tr><td>Einleitung (Auslass)</td><td>verdünnte Quellzelle, Kernanteil strömungsabhängig (Verdünnungsmodell)</td></tr></table>
    <img class="rimg" src="${vGrid}">
    <div class="cap">Abb. 5.1 – Numerisches Gitter: Oberflächenraster (NX×NY) sowie Sigma-Vernetzung im Quer- und Längsschnitt.</div></div>

    <div class="sec"><h2 id="s6">6. 3D-Ansichten der Berechnung</h2>
    <img class="rimg" src="${vOblique}">
    <div class="cap">Abb. 6.1 – Schrägansicht: Gesamtmodell mit Wärmepumpe, Rohrleitungen, Bebauung und Temperaturfeld.</div>
    <img class="rimg rimg2" src="${vTop}"><img class="rimg rimg2" style="margin-left:2%" src="${vSide}">
    <div class="cap">Abb. 6.2 / 6.3 – Draufsicht und Seitenansicht.</div>
    <img class="rimg" src="${vOutlet}">
    <div class="cap">Abb. 6.4 – Detail am Auslass: Einleitstrahl und Kaltwasserfahne.</div></div>

    <div class="sec"><h2 id="s7">7. Interaktive 3D-Ansicht</h2>
    <div id="ttview"><img id="ttimg" src="${ttSrc()}" draggable="false">
      <div class="tthint">🖱 Ziehen = drehen · Mausrad = Zoomstufe</div></div>
    <div class="cap">Abb. 7.1 – Drehbare Modellansicht (${TT.perRing} Blickwinkel × ${TT.rings} Zoomstufen).
    Hinweis: interaktiv im Report-Reiter; im PDF erscheint der aktuell gewählte Blickwinkel als Standbild.</div></div>

    <div class="sec"><h2 id="s8">8. Ergebnisdiagramme</h2>
    ${svgChart(REPCHART.tx)}<div class="cap">Abb. 8.1 – T(x) am Pumpenufer: Oberflächen-, tiefengemittelte und Sohltemperatur (die Schichtung der sohlnahen Kaltwasserfahne ist direkt ablesbar).</div>
    ${svgChart(REPCHART.ty)}<div class="cap">Abb. 8.2 – T(y): Oberflächen-, tiefengemittelte und Sohltemperatur über die Flussbreite am Querschnitt x = ${fmt(P.xCut,0)} m.</div>
    ${svgChart(REPCHART.vz)}<div class="cap">Abb. 8.3 – v(z): Geschwindigkeitsprofil über die Tiefe (Flussmitte, 1/7-Potenzgesetz).</div>
    ${svgChart(REPCHART.tz)}<div class="cap">Abb. 8.4 – T(z): vertikales Temperaturprofil am Pumpenufer.</div>
    <p class="cap">Die zugrunde liegenden Daten können in der Anwendung je Diagramm als CSV exportiert werden.</p></div>

    <div class="sec"><h2 id="s9">9. Kennzahlen (vollständig)</h2>
    <table class="rt">${kpiHTML}</table>
    ${warnHTML?`<div class="note"><b>Warnhinweise:</b> ${warnHTML}</div>`:''}</div>

    <div class="sec"><h2 id="s10">10. Hinweise &amp; Modellgrenzen</h2>
    <p>Das Modell ist ein vereinfachtes Ingenieurmodell zur Veranschaulichung und Voreinschätzung:
    stationäre Strömung mit vorgegebenem Profil (keine Impulsgleichung), turbulente Vermischung über
    effektive, anisotrope Austauschkoeffizienten, Einleitung über ein Verdünnungsmodell in der Quellzelle,
    kein Wärmeaustausch mit Atmosphäre und Sohle, konstante Stoffwerte (ρ = ${RHO} kg/m³,
    c<sub>p</sub> = ${CP} J/(kg·K)). Für Genehmigungsplanungen sind gekoppelte 3D-CFD-Rechnungen und
    Naturmessungen erforderlich. Ergebnisse gelten für die oben dokumentierten Eingangsparameter.</p></div>

    <div class="sec"><h2 id="s11">11. Mathematische Modellbeschreibung</h2>
    <h3>11.1 Transportgleichung</h3>
    <p>Für die Wassertemperatur T(x, y, z) wird die Advektions-Diffusions-Gleichung gelöst
    (stationär iteriert, Quellterm q an der Einleitstelle):</p>
    <div class="frm">$$\\frac{\\partial T}{\\partial t}
    +u\\,\\frac{\\partial T}{\\partial x}+v\\,\\frac{\\partial T}{\\partial y}+w\\,\\frac{\\partial T}{\\partial z}
    =\\frac{\\partial}{\\partial x}\\!\\left(\\alpha_x\\frac{\\partial T}{\\partial x}\\right)
    +\\frac{\\partial}{\\partial y}\\!\\left(\\alpha_y\\frac{\\partial T}{\\partial y}\\right)
    +\\frac{\\partial}{\\partial z}\\!\\left(\\alpha_z\\frac{\\partial T}{\\partial z}\\right)+q$$</div>
    <p>Die Advektion wird semi-Lagrange'sch diskretisiert (Rückwärtsverfolgung entlang der Charakteristik
    mit linearer Interpolation, unbedingt stabil) mit dem Zeitschritt</p>
    <div class="frm">$$\\Delta t=\\frac{1{,}2\\;\\Delta x}{\\max\\left(v_m,\\;0{,}1\\ \\mathrm{m/s}\\right)}
    =${mnum(fmt(1.2*dx/Math.max(P.vFlow,0.1),1))}\\ \\mathrm{s}$$</div>
    <p>Die Diffusion ist als konservative Fluss-Form (paarweiser Energieaustausch) implementiert.
    Die Austauschgewichte werden je Auflösung aus PHYSIKALISCHEN Diffusivitäten berechnet und sind
    damit gitterunabhängig (Kappe 0,45 für Stabilität):</p>
    <div class="frm">$$w_{\\alpha}=\\min\\!\\left(\\frac{2\\,D_{\\alpha}\\,\\Delta t}{\\Delta_{\\alpha}^{2}},\\;0{,}45\\right),
    \\qquad D_y=0{,}6\\,h\\,u_{*},\\quad D_z=0{,}067\\,h\\,u_{*},\\quad D_x=0{,}3\\,D_y,
    \\quad u_{*}\\approx 0{,}1\\,v_m$$</div>
    <p>Die Quer- und Vertikalvermischung dominiert, wodurch sich die Kaltwasserfahne stromab
    verbreitert und auflöst; vertikal wird das Gewicht je Spalte aus der lokalen Schichtdicke
    Δz = h(y)/${NZ} bestimmt (dünne Uferschichten mischen schnell).</p>
    <h3>11.2 Geometrie und Sigma-Koordinaten</h3>
    <p>Die Sohle ist parabelförmig mit fester Wasseroberfläche bei z = ${WSURF} m:</p>
    <div class="frm">$$h(y)=\\max\\!\\Big(h_{\\max}\\,\\big(1-\\hat y^{\\,2}\\big),\\;${mnum(fmt(MINDEPTH,1))}\\ \\mathrm{m}\\Big),
    \\qquad \\hat y=\\frac{2y}{B}-1$$
    $$z_{\\mathrm{Sohle}}(y)=${WSURF}\\ \\mathrm{m}-h(y)$$</div>
    <p>Vertikal werden ${NZ} terrainfolgende Sigma-Schichten verwendet; jede Zelle hat die lokale Dicke h(y)/${NZ}:</p>
    <div class="frm">$$z(\\sigma)=z_{\\mathrm{Sohle}}(y)+\\sigma\\,h(y),\\qquad \\sigma\\in[0,\\,1]$$</div>
    <h3>11.3 Geschwindigkeitsfeld</h3>
    <div class="frm">$$u(x,y,z)=1{,}15\\;v_m\\;\\big(1-0{,}5\\,\\hat y^{\\,2}\\big)\\;\\sigma^{1/7}$$</div>
    <p>Querprofil parabolisch (Maximum in Flussmitte), Vertikalprofil nach dem 1/7-Potenzgesetz
    turbulenter Gerinneströmungen; der Vorfaktor 1,15 normiert auf die mittlere Geschwindigkeit
    v<sub>m</sub>. An Ein- und Auslass wird eine lokale Senken- bzw. Quellströmung überlagert
    (Punktsenke/-quelle v<sub>r</sub> = Q<sub>Pumpe</sub>/(4π·r²)), die Ansaugung und Einleitung abbildet. Die reale Einzugszone der Ansaugung ist der Stromröhren-Schlauch mit Querschnitt A = Q<sub>Pumpe</sub>/u.</p>
    <h3>11.4 Wärmepumpe (COP-Kopplung)</h3>
    <div class="frm">$$P_{\\mathrm{el}}=\\frac{P_{\\mathrm{heiz}}}{\\mathrm{COP}},\\qquad
    P_{\\mathrm{entzug}}=P_{\\mathrm{heiz}}\\left(1-\\frac{1}{\\mathrm{COP}}\\right)$$
    $$\\Delta T=\\frac{P_{\\mathrm{entzug}}}{\\rho\\,c_p\\,Q_{\\mathrm{Pumpe}}},\\qquad
    Q_{\\mathrm{Pumpe}}=\\pi\\left(\\frac{d}{2}\\right)^{\\!2} v_{\\mathrm{Ansaug}}$$
    $$T_{\\mathrm{aus}}=\\max\\big(T_{\\mathrm{ein,lokal}}-\\Delta T,\\;0\\,{}^{\\circ}\\mathrm{C}\\big)$$</div>
    <h3>11.5 Einleitung als energieerhaltender Quellterm</h3>
    <p>Das Rohr mündet senkrecht von oben; der Strahl prallt auf die Sohle und durchmischt die
    Wassersäule der Einleitzone. Pro Zeitschritt wird exakt das Wärmedefizit
    P<sub>entzug</sub>·Δt konservativ in die Zone eingebracht (Gauß-Gewichte ŵ, Σŵ = 1,
    sohlwärts betont), begrenzt durch die Auslasstemperatur:</p>
    <div class="frm">$$T_{\\mathrm{Zelle}}\\;\\leftarrow\\;\\max\\!\\left(T_{\\mathrm{Zelle}}
    -\\frac{P_{\\mathrm{entzug}}\\,\\Delta t\\;\\hat w}{\\rho\\,c_p\\,V_{\\mathrm{Zelle}}},\\;T_{\\mathrm{aus}}\\right)$$</div>
    <p>Damit ist der Wärmedefizit-Strom durch jeden Querschnitt stromab exakt gleich der
    Entzugsleistung (verifiziert: Φ(x) = Σ ρ·c<sub>p</sub>·u·(T<sub>Fluss</sub>−T)·dA = P<sub>entzug</sub>
    für alle x hinter dem Auslass). Die Verdünnung mit der Strömung ergibt sich von selbst:
    mehr Durchfluss verteilt dasselbe Defizit auf mehr Wasser → wärmere, kürzere Fahne.</p>
    <h3>11.6 Turbulente Austauschgrößen und Kennzahlen</h3>
    <div class="frm">$$\\alpha_{\\mathrm{eff}}=\\alpha_{\\mathrm{mol}}+0{,}012\\;v_m\\,h_{\\max},\\qquad
    Re=\\frac{v_m\\,h_{\\max}}{\\nu},\\qquad Pe=\\frac{v_m\\,h_{\\max}}{\\alpha_{\\mathrm{eff}}}$$</div>
    <p>mit ν = 1,3·10<sup>−6</sup> m²/s (Wasser, ~10 °C) und α<sub>mol</sub> = 1,4·10<sup>−7</sup> m²/s.</p>
    <h3>11.7 Randbedingungen und Durchmischung</h3>
    <p>Einström-Rand (x = 0): T = T<sub>Fluss</sub> (Dirichlet). Ausström-Rand (x = L): Nullgradient.
    Ufer, Sohle und Wasseroberfläche: adiabat (kein Wärmestrom, kein Austausch mit der Atmosphäre).
    Vollständig durchmischte Abkühlung des Gesamtflusses:</p>
    <div class="frm">$$\\Delta T_{\\mathrm{mix}}=\\frac{P_{\\mathrm{entzug}}}{\\rho\\,c_p\\,Q_{\\mathrm{Fluss}}},\\qquad
    Q_{\\mathrm{Fluss}}=\\tfrac{2}{3}\\,B\\,h_{\\max}\\,v_m
    \\;\\;\\text{(Integral des Parabelprofils)}$$</div></div>
  `;
  // LaTeX-Formeln rendern (KaTeX); Fallback: Rohtext bleibt sichtbar, falls CDN fehlt
  try{
    if(typeof renderMathInElement==='function')
      renderMathInElement(document.getElementById('reportBody'),
        {delimiters:[{left:'$$',right:'$$',display:true}],throwOnError:false});
  }catch(e){}
  // Inhaltsverzeichnis: sanft im Report-Container scrollen
  document.querySelectorAll('#reportBody .toc a').forEach(a=>{
    a.addEventListener('click',ev=>{ev.preventDefault();
      const t=document.getElementById(a.dataset.t);
      if(t)t.scrollIntoView({behavior:'smooth',block:'start'});});
  });
  // Drehansicht: Ziehen = Winkel, Mausrad = Zoomstufe
  const ttimg=document.getElementById('ttimg');
  let dragX=null;
  ttimg.addEventListener('pointerdown',e=>{dragX=e.clientX;ttimg.setPointerCapture(e.pointerId);ttimg.style.cursor='grabbing';});
  ttimg.addEventListener('pointermove',e=>{
    if(dragX===null)return;
    const stepPx=14, moved=Math.trunc((e.clientX-dragX)/stepPx);
    if(moved!==0){TT.idx=((TT.idx-moved)%TT.perRing+TT.perRing)%TT.perRing;dragX+=moved*stepPx;ttimg.src=ttSrc();}
  });
  ttimg.addEventListener('pointerup',()=>{dragX=null;ttimg.style.cursor='grab';});
  ttimg.addEventListener('wheel',e=>{e.preventDefault();
    TT.ring=(e.deltaY<0)?Math.min(TT.rings-1,TT.ring+1):Math.max(0,TT.ring-1);
    ttimg.src=ttSrc();},{passive:false});
  showTab('report');
}
document.getElementById('btnReport').addEventListener('click',generateReport);
document.getElementById('btnRegen').addEventListener('click',generateReport);
document.getElementById('reportBody').addEventListener('click',e=>{
  if(e.target.classList.contains('rGen')) generateReport();
});
document.getElementById('btnPdf').addEventListener('click',()=>{
  /* window.print() wird in eingebetteten Vorschauen (Sandbox-iframe) stumm blockiert.
     Daher: im Top-Level-Fenster direkt drucken; eingebettet stattdessen den Bericht als
     eigenständige HTML-Datei herunterladen, die sich beim Öffnen selbst zum Druck anbietet. */
  let embedded=false;
  try{ embedded=(window.self!==window.top); }catch(e){ embedded=true; }
  if(!embedded){ window.print(); return; }
  const body=document.getElementById('reportBody').innerHTML;
  // relevante Report-Stile aus dem Stylesheet extrahieren
  let css='';
  for(const sheet of document.styleSheets){
    try{ for(const r of sheet.cssRules){
      if(r.cssText && r.cssText.indexOf('#reportBody')>=0) css+=r.cssText.replace(/#reportBody/g,'.rb')+'\n';
    }}catch(e){}
  }
  const doc=`<!DOCTYPE html><html lang="de"><head><meta charset="utf-8">
<title>Flusswasserwärmepumpe – Simulationsbericht</title>
<link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/KaTeX/0.16.9/katex.min.css">
<style>
body{margin:0;background:#eef1f4;font-family:Segoe UI,Arial,sans-serif}
.bar{position:sticky;top:0;background:#16202b;color:#dfe9f2;padding:9px 16px;font-size:13px;display:flex;gap:12px;align-items:center}
.bar button{font-size:13px;padding:6px 16px;border-radius:6px;cursor:pointer;background:#1b2734;color:#dfe9f2;border:1px solid #3fc1c9}
.rb{max-width:900px;margin:18px auto;background:#fff;color:#1c2733;padding:42px 52px;border-radius:4px;
box-shadow:0 3px 16px rgba(0,0,0,.25);font-size:13.5px;line-height:1.55}
${css}
#ttview .tthint{display:none}
@media print{.bar{display:none}body{background:#fff}.rb{box-shadow:none;max-width:none;margin:0;padding:10mm 12mm}
.rb h2{page-break-after:avoid}.rb .sec{page-break-inside:avoid}}
</style></head><body>
<div class="bar"><button onclick="window.print()">⬇ Als PDF drucken/speichern</button>
<span>Diese Datei ist der vollständige Bericht – über den Button als PDF speichern.</span></div>
<div class="rb">${body}</div>
<script>document.querySelectorAll('.toc a').forEach(a=>a.addEventListener('click',e=>{e.preventDefault();
const t=document.getElementById(a.dataset.t||a.getAttribute('href').slice(1));if(t)t.scrollIntoView({behavior:'smooth'});}));<\/script>
</body></html>`;
  const blob=new Blob([doc],{type:'text/html;charset=utf-8'});
  const a=document.createElement('a');
  a.href=URL.createObjectURL(blob); a.download='Simulationsbericht_Flusswaermepumpe.html';
  document.body.appendChild(a); a.click();
  setTimeout(()=>{URL.revokeObjectURL(a.href); a.remove();},0);
});

/* =========================================================================
   (7) INITIALISIERUNG & ANIMATIONSSCHLEIFE
   ========================================================================= */
function resize(){
  let w=sceneEl.clientWidth, h=sceneEl.clientHeight;
  if(!w||!h){ w=window.innerWidth-350; h=window.innerHeight-280; } // Fallback
  w=Math.max(w,200); h=Math.max(h,200);
  renderer.setSize(w,h,false);
  renderer.domElement.style.width='100%';
  renderer.domElement.style.height='100%';
  camera.aspect=w/h; camera.updateProjectionMatrix();
}
window.addEventListener('resize',resize);
if(window.ResizeObserver){ try{ new ResizeObserver(()=>resize()).observe(sceneEl); }catch(e){} }

function init(){
 try{
  if(!renderer) return;
  computeDerived();
  T.fill(P.tIn);
  buildEnv(); buildWater(); buildPump(); placePlanes(); buildVectors();
  refreshLabels(); refreshDerived();
  relax(80);              // Anfangsfeld einschwingen
  fieldStats(); buildIso(); buildFlow();
  resize(); setView('oblique'); drawColorbar();
  updateGridInfo(+document.getElementById('gridRes').value);
  Object.entries(TG).forEach(([id,obj])=>obj.visible=document.getElementById(id).checked);
  circParticles.visible=document.getElementById('t_stream').checked;
  jetParticles.visible=document.getElementById('t_stream').checked;
  suctParticles.visible=document.getElementById('t_stream').checked;
  if(isoGroup)  isoGroup.visible =document.getElementById('t_iso').checked;
  if(flowGroup) flowGroup.visible=document.getElementById('t_flow').checked;
  updatePlanes();
  controls.update();
  renderer.render(scene,camera);   // sofortiges Erst-Rendering
  animate();
 }catch(e){ showError('init(): '+(e&&e.stack?e.stack:e)); }
}

let frame=0, animErr=false;
function animate(){
  requestAnimationFrame(animate);
 try{
  // kontinuierliche, quasi-transiente Lösung (hält Feld im Gleichgewicht)
  relax(2);
  fieldStats();
  updatePlanes();
  updateParticles(0.4);
  updateCircuit(0.4);
  updateJet(0.4);
  updateSuction(0.4);
  updateBoat(0.4);
  updateGauges();
  // Rohrmündungen nach Temperatur einfärben (Auslass sichtbar kälter)
  if(frame%6===0){ drawColorbar(); updateKPI(); updateCharts();
    const cv=renderer.domElement, hud=document.getElementById('hud');
    if(hud) hud.textContent=`Canvas ${cv.clientWidth}×${cv.clientHeight}px · Buffer ${cv.width}×${cv.height} · `+
      `Objekte ${scene.children.length} · Kamera (${camera.position.x|0},${camera.position.y|0},${camera.position.z|0})`;
  }
  if(frame%24===0){
    if(document.getElementById('t_iso').checked)  buildIso();
    if(document.getElementById('t_flow').checked) buildFlow();
  }
  controls.update();
  renderer.render(scene,camera);
  frame++;
 }catch(e){ if(!animErr){animErr=true; showError('animate(): '+(e&&e.stack?e.stack:e));} }
}

window.addEventListener('load',init);
