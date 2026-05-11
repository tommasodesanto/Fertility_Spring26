// All interactive widgets for tour.html. Inlined by assemble_tour.py
// just before </body>. Depends on window.TOUR being defined (it is — set
// up in the inline <script> earlier in the document).

// Widget 1 — utility heatmap (pure JS, no JSON dependency) -------------
(function utilityHeatmap() {
  const svg = document.getElementById("u-heatmap");
  if (!svg) return;
  const ids = ["alpha","sigma","cbarn","hjump","parity"];
  const get = id => document.getElementById("u-" + id);
  const update = () => {
    const a  = +get("alpha").value;
    const sg = +get("sigma").value;
    const cn = +get("cbarn").value;
    const hj = +get("hjump").value;
    const n  = +get("parity").value;
    document.getElementById("u-alpha-val").textContent = a.toFixed(2);
    document.getElementById("u-sigma-val").textContent = sg.toFixed(2);
    document.getElementById("u-cbarn-val").textContent = cn.toFixed(3);
    document.getElementById("u-hjump-val").textContent = hj.toFixed(2);
    const cBar0 = 0.10, hBar0 = 0.50, hBarn = 0.10;
    const cFloor = cBar0 + cn * n;
    const hFloor = hBar0 + hBarn * n + (n >= 1 ? hj : 0);
    const cMin = cFloor + 0.10, cMax = 4.0;
    const hMin = hFloor + 0.10, hMax = 8.0;
    const W = 70, H = 35;
    const w = +svg.getAttribute("width"), h = +svg.getAttribute("height");
    const padL = 36, padB = 28, padT = 8, padR = 8;
    const cellW = (w - padL - padR) / W;
    const cellH = (h - padT - padB) / H;
    let parts = [];
    let uMin = +Infinity, uMax = -Infinity;
    const grid = [];
    for (let j = 0; j < H; j++) {
      const hh = hMin + (j / (H - 1)) * (hMax - hMin);
      const row = [];
      for (let i = 0; i < W; i++) {
        const cc = cMin + (i / (W - 1)) * (cMax - cMin);
        const cArg = cc - cFloor, hArg = hh - hFloor;
        let u = NaN;
        if (cArg > 0 && hArg > 0) {
          const inner = Math.pow(cArg, a) * Math.pow(hArg, 1 - a);
          u = (Math.pow(inner, 1 - sg) - 1) / (1 - sg);
          if (Number.isFinite(u)) {
            uMin = Math.min(uMin, u);
            uMax = Math.max(uMax, u);
          }
        }
        row.push(u);
      }
      grid.push(row);
    }
    if (!Number.isFinite(uMin)) { uMin = 0; uMax = 1; }
    for (let j = 0; j < H; j++) {
      for (let i = 0; i < W; i++) {
        const u = grid[j][i];
        const t = Number.isFinite(u) ? Math.max(0, Math.min(1, (u - uMin) / (uMax - uMin || 1))) : -1;
        let fill;
        if (t < 0) fill = "#f0eee8";
        else {
          // diverging-ish: dark blue (low) -> mid -> warm yellow (high)
          const r = Math.round(43 + (240 - 43) * t);
          const g = Math.round(87 + (200 - 87) * t);
          const b = Math.round(128 + (110 - 128) * t);
          fill = `rgb(${r},${g},${b})`;
        }
        parts.push(`<rect x="${(padL + i * cellW).toFixed(2)}" y="${(padT + (H - 1 - j) * cellH).toFixed(2)}" width="${(cellW + 0.5).toFixed(2)}" height="${(cellH + 0.5).toFixed(2)}" fill="${fill}"/>`);
      }
    }
    parts.push(`<text x="${padL}" y="${h - 8}" font-size="11" fill="#444">c → (${cMin.toFixed(2)} … ${cMax.toFixed(1)})</text>`);
    parts.push(`<text x="${w - padR}" y="${h - 8}" font-size="11" text-anchor="end" fill="#444">u ∈ [${uMin.toFixed(2)}, ${uMax.toFixed(2)}]</text>`);
    parts.push(`<text x="6" y="${h / 2 + 4}" font-size="11" fill="#444" transform="rotate(-90 12 ${h / 2})">h ↑ (${hMin.toFixed(2)} … ${hMax.toFixed(1)})</text>`);
    svg.innerHTML = parts.join("");
  };
  ids.forEach(id => get(id).addEventListener("input", update));
  update();
})();

// Helper: populate a <select> ----------------------------------------
function _fillSelect(id, values, formatter, defaultIdx) {
  const el = document.getElementById(id);
  if (!el) return;
  el.innerHTML = values.map((v, i) => {
    const sel = (i === (defaultIdx || 0)) ? " selected" : "";
    return `<option value="${v}"${sel}>${formatter ? formatter(v) : v}</option>`;
  }).join("");
}

// Widget 2 — tenure decision tree (data-driven) ----------------------
(function tenureTree() {
  const svg = document.getElementById("tt-svg");
  if (!svg) return;
  const data = TOUR.data;
  const axes = data.policies && data.policies.axes;
  if (!axes) {
    svg.innerHTML = '<text x="20" y="40" fill="#a84040" font-size="13">no slice data available</text>';
    return;
  }
  _fillSelect("tt-j", axes.j, v => `j=${v} (age ${+v + 18})`, Math.min(3, axes.j.length - 1));
  _fillSelect("tt-n", Array.from({length: axes.npar}, (_, i) => i), v => `n=${v}`);
  _fillSelect("tt-cs", axes.cs, v => `cs=${v}`);
  _fillSelect("tt-to", Array.from({length: axes.nt}, (_, i) => i), v => v == 0 ? "rent" : `own ${v}`);
  _fillSelect("tt-i", Array.from({length: axes.I}, (_, i) => i), v => v == 0 ? "Periphery" : "Center");
  const bGrid = (data.grids && data.grids.b_grid) || Array.from({length: 60}, (_, i) => i);
  const bSlider = document.getElementById("tt-b");
  bSlider.max = bGrid.length - 1;
  bSlider.value = Math.min(20, bGrid.length - 1);
  const update = () => {
    const j = +document.getElementById("tt-j").value;
    const n = +document.getElementById("tt-n").value;
    const cs = +document.getElementById("tt-cs").value;
    const to = +document.getElementById("tt-to").value;
    const i = +document.getElementById("tt-i").value;
    const bIdx = +bSlider.value;
    document.getElementById("tt-b-val").textContent = `b ≈ ${TOUR.fmt(bGrid[bIdx], 2)}`;
    const key = `${j}|${n}|${cs}|${to}|${i}`;
    const slice = data.policies.slices[key];
    if (!slice) {
      svg.innerHTML = `<text x="20" y="40" fill="#a84040" font-size="12">no slice for ${key} — pick a different combination</text>`;
      return;
    }
    const tn = slice.tn[bIdx];
    const branches = ["rent", "own 1", "own 2", "own 3", "own 4"].slice(0, axes.nt);
    const w = +svg.getAttribute("width"), h = +svg.getAttribute("height");
    const branchW = 110;
    const totalW = branches.length * branchW + (branches.length - 1) * 12;
    const startX = (w - totalW) / 2;
    let parts = [
      `<text x="${w / 2}" y="22" text-anchor="middle" font-size="14" font-family="Inter">Tenure choice at b = ${TOUR.fmt(bGrid[bIdx], 2)}: <tspan font-weight="600" fill="#3a8060">${branches[tn]}</tspan></text>`,
      `<circle cx="${w / 2}" cy="56" r="9" fill="#2b5780"/>`,
      `<text x="${w / 2}" y="80" text-anchor="middle" font-size="11" fill="#555" font-family="Inter">origin: ${i == 0 ? "Periphery" : "Center"}, t_o = ${branches[to]}</text>`,
    ];
    branches.forEach((label, k) => {
      const x = startX + k * (branchW + 12);
      const y = 130;
      const isOpt = (k === tn);
      parts.push(`<line x1="${w / 2}" y1="65" x2="${x + branchW / 2}" y2="${y}" stroke="${isOpt ? '#3a8060' : '#bbb'}" stroke-width="${isOpt ? 2 : 1}" opacity="${isOpt ? 1 : 0.5}"/>`);
      parts.push(`<rect x="${x}" y="${y}" width="${branchW}" height="50" fill="${isOpt ? '#e6f5ec' : '#fff'}" stroke="${isOpt ? '#3a8060' : '#bbb'}" stroke-width="${isOpt ? 2 : 1}" rx="4"/>`);
      parts.push(`<text x="${x + branchW / 2}" y="${y + 22}" text-anchor="middle" font-size="13" font-family="Inter" font-weight="${isOpt ? 600 : 400}" fill="${isOpt ? '#3a8060' : '#555'}">${label}</text>`);
      parts.push(`<text x="${x + branchW / 2}" y="${y + 40}" text-anchor="middle" font-size="11" fill="#666" font-family="JetBrains Mono">tn = ${k}</text>`);
    });
    parts.push(`<text x="20" y="${h - 36}" font-size="12" fill="#444" font-family="Inter">savings $b'$ at chosen tenure: <tspan font-family="JetBrains Mono">${TOUR.fmt(slice.bp[bIdx], 3)}</tspan></text>`);
    parts.push(`<text x="20" y="${h - 18}" font-size="12" fill="#444" font-family="Inter">consumption $c$ at chosen tenure: <tspan font-family="JetBrains Mono">${TOUR.fmt(slice.c[bIdx], 3)}</tspan></text>`);
    svg.innerHTML = parts.join("");
  };
  ["tt-j","tt-n","tt-cs","tt-to","tt-i","tt-b"].forEach(id => {
    const el = document.getElementById(id);
    if (el) el.addEventListener("input", update);
  });
  update();
})();

// Widget 3a — location logit temperature ------------------------------
(function locLogit() {
  const svg = document.getElementById("ll-bars");
  if (!svg) return;
  const baseV = [-0.20, 0.30];
  const update = () => {
    const k = +document.getElementById("ll-kappa").value;
    document.getElementById("ll-kappa-val").textContent = k.toFixed(2);
    const m = Math.max.apply(null, baseV);
    const ex = baseV.map(v => Math.exp((v - m) / k));
    const Z = ex.reduce((a, b) => a + b, 0);
    const probs = ex.map(e => e / Z);
    TOUR.bars(svg, ["Periphery", "Center"], probs, { max: 1, color: "#5b87b3" });
  };
  document.getElementById("ll-kappa").addEventListener("input", update);
  update();
})();

// Widget 3b — fertility logit temperature -----------------------------
(function fertLogit() {
  const svg = document.getElementById("fl-bars");
  if (!svg) return;
  const baseV = [-0.05, 0.10, -0.05, -0.30];
  const update = () => {
    const k = +document.getElementById("fl-kappa").value;
    document.getElementById("fl-kappa-val").textContent = k.toFixed(2);
    const m = Math.max.apply(null, baseV);
    const ex = baseV.map(v => Math.exp((v - m) / k));
    const Z = ex.reduce((a, b) => a + b, 0);
    const probs = ex.map(e => e / Z);
    TOUR.bars(svg, ["n'=0", "n'=1", "n'=2", "n'=3"], probs, { max: 1, color: "#3a8060" });
  };
  document.getElementById("fl-kappa").addEventListener("input", update);
  update();
})();

// Widget 4 — backward induction sweep --------------------------------
(function backwardSweep() {
  const svg = document.getElementById("bw-svg");
  if (!svg) return;
  const J = (TOUR.data.grids && TOUR.data.grids.J) || 60;
  let j = J - 1;
  let timer = null;
  let stageTick = 0;
  const stages = [
    {label: "Vd", desc: "savings + tenure-cond", color: "#5b87b3"},
    {label: "VH", desc: "tenure max", color: "#3a8060"},
    {label: "VI", desc: "location logit", color: "#a84040"},
    {label: "V",  desc: "fertility logit", color: "#c8961f"},
  ];
  const render = () => {
    const w = +svg.getAttribute("width"), h = +svg.getAttribute("height");
    document.getElementById("bw-age").textContent = `age j = ${j}  (year ${j + 18})`;
    const cursorX = 50 + (J - 1 - j) / Math.max(J - 1, 1) * (w - 100);
    let parts = [
      `<line x1="50" y1="${h - 36}" x2="${w - 50}" y2="${h - 36}" stroke="#bbb"/>`,
      `<text x="50" y="${h - 14}" font-size="11" fill="#666" font-family="Inter">j=${J - 1} (terminal)</text>`,
      `<text x="${w - 50}" y="${h - 14}" font-size="11" text-anchor="end" fill="#666" font-family="Inter">j=0 (entry)</text>`,
      `<line x1="${cursorX}" y1="${h - 50}" x2="${cursorX}" y2="${h - 28}" stroke="#2b5780" stroke-width="2"/>`,
      `<circle cx="${cursorX}" cy="${h - 36}" r="6" fill="#2b5780"/>`,
      `<text x="${cursorX}" y="${h - 60}" text-anchor="middle" font-size="11" fill="#2b5780" font-family="Inter">j = ${j}</text>`,
    ];
    const totalW = 4 * 110 + 3 * 14;
    const startX = (w - totalW) / 2;
    stages.forEach((s, k) => {
      const x = startX + k * 124, y = 24;
      const active = (stageTick % 4) === k;
      parts.push(`<rect x="${x}" y="${y}" width="110" height="48" fill="${active ? '#e6f5ec' : '#fff'}" stroke="${active ? s.color : '#bbb'}" stroke-width="${active ? 2 : 1}" rx="4"/>`);
      parts.push(`<text x="${x + 55}" y="${y + 21}" text-anchor="middle" font-size="14" font-family="JetBrains Mono" font-weight="600" fill="${active ? s.color : '#444'}">${s.label}</text>`);
      parts.push(`<text x="${x + 55}" y="${y + 38}" text-anchor="middle" font-size="10.5" fill="${active ? s.color : '#666'}" font-family="Inter">${s.desc}</text>`);
      if (k < 3) {
        const ax = x + 110, ay = y + 24;
        parts.push(`<path d="M${ax} ${ay} l8 0 m-3 -3 l3 3 l-3 3" stroke="${active ? s.color : '#bbb'}" stroke-width="${active ? 1.6 : 1}" fill="none"/>`);
      }
    });
    svg.innerHTML = parts.join("");
  };
  document.getElementById("bw-play").onclick = function () {
    if (timer) {
      clearInterval(timer); timer = null;
      this.textContent = "Play";
      return;
    }
    this.textContent = "Pause";
    timer = setInterval(() => {
      stageTick = (stageTick + 1) % 4;
      if (stageTick === 0) j = j > 0 ? j - 1 : J - 1;
      render();
    }, 220);
  };
  document.getElementById("bw-step").onclick = function () {
    stageTick = (stageTick + 1) % 4;
    if (stageTick === 0) j = j > 0 ? j - 1 : J - 1;
    render();
  };
  document.getElementById("bw-reset").onclick = function () {
    j = J - 1; stageTick = 0;
    if (timer) { clearInterval(timer); timer = null; document.getElementById("bw-play").textContent = "Play"; }
    render();
  };
  render();
})();

// Widget 5 — policy explorer ------------------------------------------
(function policyExplorer() {
  const data = TOUR.data;
  const axes = data.policies && data.policies.axes;
  if (!axes) return;
  _fillSelect("exp-j", axes.j, v => `j=${v} (age ${+v + 18})`, Math.min(4, axes.j.length - 1));
  _fillSelect("exp-n", Array.from({length: axes.npar}, (_, i) => i), v => `n=${v}`);
  _fillSelect("exp-cs", axes.cs, v => `cs=${v}`);
  _fillSelect("exp-to", Array.from({length: axes.nt}, (_, i) => i), v => v == 0 ? "rent" : `own ${v}`);
  _fillSelect("exp-i", Array.from({length: axes.I}, (_, i) => i), v => v == 0 ? "Periphery" : "Center");
  const slicesEl = document.getElementById("explorer-axis-summary");
  if (slicesEl) {
    slicesEl.textContent = `j ∈ {${axes.j.join(",")}}, n ∈ [0..${axes.npar - 1}], cs ∈ {${axes.cs.join(",")}}, t_o ∈ [0..${axes.nt - 1}], i ∈ [0..${axes.I - 1}] — ${Object.keys(data.policies.slices).length} slices`;
  }
  const bGrid = (data.grids && data.grids.b_grid) || Array.from({length: 60}, (_, i) => i);
  const update = () => {
    const j = +document.getElementById("exp-j").value;
    const n = +document.getElementById("exp-n").value;
    const cs = +document.getElementById("exp-cs").value;
    const to = +document.getElementById("exp-to").value;
    const i = +document.getElementById("exp-i").value;
    document.getElementById("exp-j-years").textContent = `(${j + 18} years old)`;
    const key = `${j}|${n}|${cs}|${to}|${i}`;
    const s = data.policies.slices[key];
    const els = ["exp-bp","exp-c","exp-tn","exp-loc","exp-fert"].map(id => document.getElementById(id));
    if (!s) {
      els.forEach(el => el.innerHTML = `<text x="10" y="20" font-size="11" fill="#a84040">no slice for ${key}</text>`);
      document.getElementById("exp-summary").innerHTML = `<span class="small" style="color:#a84040">slice ${key} not in dump</span>`;
      return;
    }
    TOUR.lineChart(els[0], bGrid, s.bp, { color: "#2b5780" });
    TOUR.lineChart(els[1], bGrid, s.c,  { color: "#a84040" });
    TOUR.lineChart(els[2], bGrid, s.tn, { color: "#3a8060", yMin: -0.2, yMax: axes.nt - 0.8 });
    if (s.loc_probs) {
      const probsCenter = s.loc_probs.map(row => row[1]);
      TOUR.lineChart(els[3], bGrid, probsCenter, { color: "#5b87b3", yMin: 0, yMax: 1, label: "Pr(Center | b)" });
    } else {
      els[3].innerHTML = '<text x="10" y="20" font-size="11" fill="#888">no loc_probs in slice</text>';
    }
    if (s.fert_probs) {
      const mid = Math.floor(s.fert_probs.length / 2);
      const ps = s.fert_probs[mid];
      TOUR.bars(els[4], ps.map((_, k) => `n'=${k}`), ps, { max: 1, color: "#c8961f" });
    } else {
      els[4].innerHTML = '<text x="10" y="20" font-size="11" fill="#888">fertility logit not active here<tspan x="10" dy="14">(only at j fertile, n=0, cs=0)</tspan></text>';
    }
    const mid = Math.floor(s.bp.length / 2);
    const tnAtMid = s.tn[mid];
    const tnLabel = tnAtMid === 0 ? "rent" : `own ${tnAtMid}`;
    document.getElementById("exp-summary").innerHTML =
      `<div><strong>State:</strong> j=${j} (age ${j + 18}), n=${n}, cs=${cs}, t_o=${to == 0 ? 'rent' : 'own ' + to}, i=${i == 0 ? 'P' : 'C'}</div>` +
      `<div style="margin-top:6px;"><strong>At median b (≈${TOUR.fmt(bGrid[mid], 2)}):</strong></div>` +
      `<div class="mono">b' = ${TOUR.fmt(s.bp[mid], 3)}</div>` +
      `<div class="mono">c  = ${TOUR.fmt(s.c[mid], 3)}</div>` +
      `<div class="mono">tn* = ${tnAtMid} (${tnLabel})</div>` +
      (s.loc_probs ? `<div class="mono">P(Center | b_med) = ${TOUR.fmt(s.loc_probs[mid][1], 3)}</div>` : "");
  };
  ["exp-j","exp-n","exp-cs","exp-to","exp-i"].forEach(id => {
    const el = document.getElementById(id);
    if (el) el.addEventListener("input", update);
  });
  update();
})();

// GE replay console ---------------------------------------------------
(function geReplay() {
  const trace = TOUR.data.ge_trace || [];
  const traj = document.getElementById("ge-trajectory");
  const errSvg = document.getElementById("ge-error");
  const log = document.getElementById("ge-log");
  if (!trace.length || !traj || !errSvg || !log) {
    if (traj) traj.innerHTML = '<text x="20" y="40" font-size="12" fill="#888">no ge_trace in dump</text>';
    return;
  }
  let cur = 0, timer = null;
  const xs = trace.map(t => t.p[0]);
  const ys = trace.map(t => t.p[1]);
  const xMin = Math.min.apply(null, xs), xMax = Math.max.apply(null, xs);
  const yMin = Math.min.apply(null, ys), yMax = Math.max.apply(null, ys);
  const ws = +traj.getAttribute("width"), hs = +traj.getAttribute("height");
  const padL = 36, padB = 30, padT = 14, padR = 14;
  const sx = x => padL + (x - xMin) / (xMax - xMin || 1) * (ws - padL - padR);
  const sy = y => hs - padB - (y - yMin) / (yMax - yMin || 1) * (hs - padT - padB);

  const renderTraj = () => {
    let parts = [
      `<line x1="${padL}" y1="${hs - padB}" x2="${ws - padR}" y2="${hs - padB}" stroke="#bbb"/>`,
      `<line x1="${padL}" y1="${padT}" x2="${padL}" y2="${hs - padB}" stroke="#bbb"/>`,
      `<text x="${padL}" y="${hs - 8}" font-size="10" fill="#666">p_P (Periphery): ${TOUR.fmt(xMin, 2)}..${TOUR.fmt(xMax, 2)}</text>`,
      `<text x="6" y="${padT + 8}" font-size="10" fill="#666">p_C: ${TOUR.fmt(yMax, 2)}</text>`,
      `<text x="6" y="${hs - padB}" font-size="10" fill="#666">${TOUR.fmt(yMin, 2)}</text>`,
    ];
    for (let k = 0; k <= cur; k++) {
      const t = trace[k];
      if (k > 0) {
        const p = trace[k - 1];
        parts.push(`<line x1="${sx(p.p[0]).toFixed(1)}" y1="${sy(p.p[1]).toFixed(1)}" x2="${sx(t.p[0]).toFixed(1)}" y2="${sy(t.p[1]).toFixed(1)}" stroke="#5b87b3" stroke-width="1.5" opacity="0.7"/>`);
      }
      parts.push(`<circle cx="${sx(t.p[0]).toFixed(1)}" cy="${sy(t.p[1]).toFixed(1)}" r="${k === cur ? 5 : 3}" fill="${k === cur ? '#a84040' : (t.mode === 'F' ? '#2b5780' : '#5b87b3')}"/>`);
    }
    parts.push(`<text x="${ws - padR}" y="${padT + 12}" text-anchor="end" font-size="11" font-family="Inter" fill="#444">iter ${cur + 1}/${trace.length}</text>`);
    traj.innerHTML = parts.join("");

    const xsi = trace.slice(0, cur + 1).map((_, i) => i + 1);
    const ysi = trace.slice(0, cur + 1).map(t => Math.log10(Math.max(t.err, 1e-6)));
    TOUR.lineChart(errSvg, xsi, ysi, { color: "#a84040", label: "log₁₀(err)", xMin: 1, xMax: trace.length });

    log.innerHTML = trace.slice(0, cur + 1).map(t =>
      `${t.mode}${String(t.iter).padStart(2)}: p=[${t.p.map(x => x.toFixed(3)).join(",")}]  err_p=${t.err_p.toFixed(4)}  err_e=${t.err_e.toFixed(4)}  err=${t.err.toFixed(4)}  own=${(t.own_rate * 100).toFixed(1)}%  TFR=${t.TFR.toFixed(2)}`
    ).join("\n");
    log.scrollTop = log.scrollHeight;
    document.getElementById("ge-iter").textContent = `iter ${cur + 1} / ${trace.length}`;
  };

  document.getElementById("ge-play").onclick = function () {
    if (timer) {
      clearInterval(timer); timer = null;
      this.textContent = "Play";
      return;
    }
    this.textContent = "Pause";
    timer = setInterval(() => {
      if (cur < trace.length - 1) cur++;
      else { clearInterval(timer); timer = null; document.getElementById("ge-play").textContent = "Play"; }
      renderTraj();
    }, 380);
  };
  document.getElementById("ge-restart").onclick = function () {
    if (timer) { clearInterval(timer); timer = null; document.getElementById("ge-play").textContent = "Play"; }
    cur = 0;
    renderTraj();
  };
  renderTraj();
})();

// Moments table populator ---------------------------------------------
(function momentsTable() {
  const tbody = document.getElementById("moments-body");
  const m = TOUR.data.moments;
  if (!tbody || !m) return;
  const targets = m.targets || {};
  const model = m.model || {};
  const weights = m.weights || {};
  const keys = Object.keys(targets).sort();
  if (!keys.length) {
    const s = TOUR.data.summary || {};
    tbody.innerHTML =
      `<tr><td><code>TFR</code></td><td class="num">—</td><td class="num">${TOUR.fmt(s.TFR, 2)}</td><td class="num">—</td><td class="num">—</td></tr>` +
      `<tr><td><code>own_rate</code></td><td class="num">—</td><td class="num">${TOUR.fmt(s.own_rate, 3)}</td><td class="num">—</td><td class="num">—</td></tr>`;
    return;
  }
  tbody.innerHTML = keys.map(k => {
    const t = targets[k];
    const v = model[k];
    const w = weights[k];
    const miss = (Number.isFinite(t) && Number.isFinite(v) && Math.abs(t) > 1e-6) ?
      ((v - t) / Math.abs(t) * 100).toFixed(1) + "%" : "—";
    return `<tr>
      <td><code>${k}</code></td>
      <td class="num">${TOUR.fmt(t, 3)}</td>
      <td class="num">${TOUR.fmt(v, 3)}</td>
      <td class="num" style="color:${Number.isFinite(v) && Math.abs((v - t) / Math.abs(t || 1)) > 0.20 ? '#a84040' : '#3a8060'};">${miss}</td>
      <td class="num">${TOUR.fmt(w, 1)}</td>
    </tr>`;
  }).join("");
})();

// Cross-highlight on dual-view panels ---------------------------------
// Identifiers in math marked .xh-tag and in code marked .xh-tag link by data-xh.
(function crossHighlight() {
  document.querySelectorAll(".dual [data-xh]").forEach(el => {
    el.style.cursor = "help";
    el.addEventListener("mouseenter", () => {
      const xh = el.getAttribute("data-xh");
      document.querySelectorAll(`.dual [data-xh="${xh}"]`).forEach(other => {
        other.classList.add("xh-active");
      });
    });
    el.addEventListener("mouseleave", () => {
      document.querySelectorAll(".xh-active").forEach(other => other.classList.remove("xh-active"));
    });
  });
})();
