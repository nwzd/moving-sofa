// Gerver's optimal sofa in the rotating-hallway form A(α) = (f(α) cos α, g(α) sin α).
// Everything on the page is computed here from the closed forms and their constants.
(() => {
  'use strict';

  const PI = Math.PI;
  const K = {
    phi: 0.039177364790083641, theta: 0.68130150938272489,
    a: 1.2103224220726888, K1: 0.4034408073575629,
    b1: -0.5276245980267846, b2: 0.9202583851606376, K2: -0.3054160633413417, L2: 0.4724066197508055,
    c: 0.6260455228484659, L3: 0.8896264790032219,
  };
  const ELLIPSE = { r: 0.605514, t: 0.667834 };
  const AREA = '2.2195317';
  const deg = d => d * PI / 180;

  // Inner corner A(α) of the hallway, α in [0, π].
  function corner(al) {
    const m = al > PI / 2, a2 = m ? PI - al : al, s = a2 / 2, cs = Math.cos(s), sn = Math.sin(s);
    let X, Y;
    if (a2 <= 2 * K.phi) {
      X = K.a * Math.cos(a2) - Math.sin(a2) / 4 - cs + sn / 2 + K.K1;
      Y = K.a * Math.sin(a2) + (1 + Math.cos(a2)) / 4 - sn - cs / 2;
    } else if (a2 <= 2 * K.theta) {
      const P = K.b2 + K.b1 * s - s * s / 4, Q = 1 + K.b1 - s / 2;
      X = P * cs + Q * sn + K.K2;
      Y = P * sn - Q * cs + K.L2;
    } else {
      const p = K.c - PI / 4, w = PI / 4 - s;
      X = p * (cs - sn) + w * (cs + sn);
      Y = p * (cs + sn) - w * (cs - sn) + K.L3;
    }
    return { x: m ? -X : X, y: Y };
  }
  const fOf = al => Math.abs(Math.cos(al)) < 1e-9 ? (K.c - PI / 4 + 1) / Math.SQRT2 : corner(al).x / Math.cos(al);
  const gOf = al => Math.sin(al) < 1e-9 ? K.a - 0.5 : corner(al).y / Math.sin(al);
  const BREAKS = [2 * K.phi, 2 * K.theta, PI - 2 * K.theta, PI - 2 * K.phi].map(b => b * 180 / PI);

  const css = name => getComputedStyle(document.documentElement).getPropertyValue(name).trim();
  const SVGNS = 'http://www.w3.org/2000/svg';
  function el(tag, attrs, parent, text) {
    const e = document.createElementNS(SVGNS, tag);
    for (const [k, v] of Object.entries(attrs)) e.setAttribute(k, v);
    if (text != null) e.textContent = text;
    if (parent) parent.appendChild(e);
    return e;
  }

  // ---------- r = f(α), t = g(α) chart ----------
  function fgChart(svg, box, withEllipse) {
    const { left, right, y0, y1, vmin, vmax, ticks, labelY } = box;
    const cx = a => left + a / 180 * (right - left);
    const cy = v => y1 - (v - vmin) / (vmax - vmin) * (y1 - y0);
    const axis = 'stroke: var(--wall); stroke-width: 1';
    const quiet = 'fill: var(--muted)';
    el('line', { x1: left, x2: left, y1: y0, y2: y1, style: axis }, svg);
    el('line', { x1: left, x2: right, y1: y1, y2: y1, style: axis }, svg);
    for (const v of ticks) {
      el('line', { x1: left - 5, x2: left, y1: cy(v), y2: cy(v), style: axis }, svg);
      el('text', { x: left - 9, y: cy(v) + 4, 'text-anchor': 'end', 'font-size': 12, style: quiet }, svg, v.toFixed(2));
    }
    for (const a of [0, 45, 90, 135, 180]) {
      el('line', { x1: cx(a), x2: cx(a), y1: y1, y2: y1 + 5, style: axis }, svg);
      el('text', { x: cx(a), y: labelY, 'text-anchor': 'middle', 'font-size': 12, style: quiet }, svg, a + '°');
    }
    for (const b of BREAKS) {
      el('line', { x1: cx(b), x2: cx(b), y1: y0, y2: y1, style: 'stroke: var(--wall); stroke-dasharray: 2 3; opacity: 0.7' }, svg);
    }
    if (withEllipse) {
      for (const [v, name, dy] of [[ELLIPSE.r, 'best ellipse r', -7], [ELLIPSE.t, 'best ellipse t', 16]]) {
        el('line', { x1: left, x2: right, y1: cy(v), y2: cy(v), style: 'stroke: var(--muted); stroke-dasharray: 6 4' }, svg);
        el('text', { x: cx(40), y: cy(v) + dy, 'text-anchor': 'middle', 'font-size': 12, style: quiet }, svg, `${name} = ${v.toFixed(4)}`);
      }
    }
    const series = [[fOf, 'var(--accent)', 'r = f(α)'], [gOf, 'var(--accent-2)', 't = g(α)']];
    for (const [fn, color, name] of series) {
      let d = '';
      for (let a = 0; a <= 180; a += 0.5) d += (a ? ' L' : 'M') + cx(a).toFixed(1) + ' ' + cy(fn(deg(a))).toFixed(1);
      el('path', { d, style: `fill: none; stroke: ${color}; stroke-width: 2.5` }, svg);
      el('text', { x: right + 13, y: cy(fn(PI)) + 4, 'font-size': 13, 'font-weight': 600, style: `fill: ${color}` }, svg, name);
    }
    return { cx, cy };
  }

  // ---------- Canvas scene: hallway, corner path, carved material ----------
  const VIEW = { x0: -1.85, x1: 1.85, y0: -0.45, y1: 1.45 };
  const M = 6;
  const INNER = [[0, 0], [0, -M], [-M, -M], [-M, 0]];                    // region below both inner walls
  const OUTER = [[1, -M], [M, -M], [M, M], [-M, M], [-M, 1], [1, 1]];     // region beyond the outer walls
  const CORRIDOR = [[1, 1], [-M, 1], [-M, 0], [0, 0], [0, -M], [1, -M]];
  const OUTER_WALL = [[-M, 1], [1, 1], [1, -M]];
  const INNER_WALL = [[-M, 0], [0, 0], [0, -M]];
  const STEP = deg(0.25);

  function pose(al, q) {
    const A = corner(al), t = al / 2, c = Math.cos(t), s = Math.sin(t);
    return [A.x + c * q[0] - s * q[1], A.y + s * q[0] + c * q[1]];
  }

  class Scene {
    constructor(canvas) {
      this.canvas = canvas;
      this.ctx = canvas.getContext('2d');
      this.mat = document.createElement('canvas');
      this.fin = document.createElement('canvas');
      this.width = 0;
      this.resize(true);
    }
    resize(force) {
      const w = this.canvas.clientWidth || 720;
      if (!force && w === this.width) return false;
      this.width = w;
      const h = w * (VIEW.y1 - VIEW.y0) / (VIEW.x1 - VIEW.x0);
      this.canvas.style.height = h + 'px';
      this.dpr = window.devicePixelRatio || 1;
      this.W = Math.round(w * this.dpr);
      this.H = Math.round(h * this.dpr);
      this.S = this.W / (VIEW.x1 - VIEW.x0);
      for (const c of [this.canvas, this.mat, this.fin]) { c.width = this.W; c.height = this.H; }
      this.col = { accent: css('--accent'), ink: css('--ink'), wall: css('--wall'), corridor: css('--corridor'), panel: css('--panel') };
      this.buildFinal();
      this.carved = -1;
      return true;
    }
    px(x) { return (x - VIEW.x0) * this.S; }
    py(y) { return (VIEW.y1 - y) * this.S; }
    path(ctx, pts, close) {
      ctx.beginPath();
      pts.forEach(([x, y], i) => (i ? ctx.lineTo(this.px(x), this.py(y)) : ctx.moveTo(this.px(x), this.py(y))));
      if (close) ctx.closePath();
    }
    slab(ctx) {
      ctx.globalCompositeOperation = 'source-over';
      ctx.clearRect(0, 0, this.W, this.H);
      ctx.fillStyle = this.col.accent;
      ctx.fillRect(this.px(VIEW.x0), this.py(1), (VIEW.x1 - VIEW.x0) * this.S, this.S);
    }
    remove(ctx, al) {
      ctx.globalCompositeOperation = 'destination-out';
      ctx.fillStyle = '#000';
      for (const shape of [INNER, OUTER]) { this.path(ctx, shape.map(q => pose(al, q)), true); ctx.fill(); }
      ctx.globalCompositeOperation = 'source-over';
    }
    buildFinal() {
      const f = this.fin.getContext('2d');
      this.slab(f);
      for (let a = 0; a < PI; a += STEP) this.remove(f, a);
      this.remove(f, PI);
    }
    carveTo(al) {
      const m = this.mat.getContext('2d');
      if (this.carved < 0 || al < this.carved) { this.slab(m); this.remove(m, 0); this.carved = 0; }
      while (this.carved + STEP < al) { this.carved += STEP; this.remove(m, this.carved); }
      this.remove(m, al);
      this.carved = al;
    }
    trace(ctx, upTo, dashed) {
      const pts = [];
      for (let a = 0; a < upTo; a += deg(0.5)) { const A = corner(a); pts.push([A.x, A.y]); }
      const A = corner(upTo); pts.push([A.x, A.y]);
      ctx.save();
      ctx.strokeStyle = this.col.ink;
      ctx.lineWidth = 1.6 * this.dpr;
      if (dashed) ctx.setLineDash([5 * this.dpr, 4 * this.dpr]);
      this.path(ctx, pts, false);
      ctx.stroke();
      ctx.restore();
    }
    hallway(ctx, al, corridor) {
      if (corridor) { ctx.fillStyle = this.col.corridor; this.path(ctx, CORRIDOR.map(q => pose(al, q)), true); ctx.fill(); }
      ctx.strokeStyle = this.col.wall;
      ctx.lineWidth = 2 * this.dpr;
      ctx.lineJoin = 'miter';
      for (const wall of [OUTER_WALL, INNER_WALL]) { this.path(ctx, wall.map(q => pose(al, q)), false); ctx.stroke(); }
      const A = corner(al);
      ctx.fillStyle = this.col.ink;
      ctx.beginPath();
      ctx.arc(this.px(A.x), this.py(A.y), 4 * this.dpr, 0, 2 * PI);
      ctx.fill();
    }
    renderTurn(al) {
      const ctx = this.ctx;
      ctx.clearRect(0, 0, this.W, this.H);
      ctx.globalAlpha = 0.14; ctx.drawImage(this.fin, 0, 0);
      this.carveTo(al);
      ctx.globalAlpha = 0.36; ctx.drawImage(this.mat, 0, 0);
      ctx.globalAlpha = 1;
      this.trace(ctx, al, false);
      this.hallway(ctx, al, true);
    }
    renderShape() {
      const ctx = this.ctx;
      ctx.clearRect(0, 0, this.W, this.H);
      ctx.globalAlpha = 0.32; ctx.drawImage(this.fin, 0, 0);
      ctx.globalAlpha = 1;
      this.hallway(ctx, PI / 2, false);
      this.trace(ctx, PI, true);
      ctx.fillStyle = this.col.ink;
      for (let d = 0; d <= 180; d += 30) {
        const A = corner(deg(d));
        ctx.beginPath(); ctx.arc(this.px(A.x), this.py(A.y), 3 * this.dpr, 0, 2 * PI); ctx.fill();
      }
      ctx.font = `600 ${14 * this.dpr}px ${css('--sans')}`;
      ctx.textAlign = 'center';
      ctx.fillText(`Gerver's sofa · area ${AREA}`, this.px(0), this.py(0.84));
    }
  }

  function watch(scene, draw) {
    new ResizeObserver(() => { if (scene.resize(false)) draw(); }).observe(scene.canvas);
    matchMedia('(prefers-color-scheme: dark)').addEventListener('change', () => { scene.resize(true); draw(); });
  }

  // ---------- Figure: the animated turn ----------
  function setupTurn() {
    const canvas = document.getElementById('turn');
    if (!canvas) return;
    const scene = new Scene(canvas);
    const svg = document.getElementById('turn-chart');
    const box = { left: 56, right: 630, y0: 14, y1: 128, vmin: 0.58, vmax: 0.72, ticks: [0.60, 0.64, 0.68, 0.72], labelY: 148 };
    const { cx, cy } = fgChart(svg, box, false);
    const cursor = el('line', { y1: box.y0 - 4, y2: box.y1, style: 'stroke: var(--ink); stroke-width: 1.25' }, svg);
    const dotF = el('circle', { r: 5, style: 'fill: var(--accent)' }, svg);
    const dotG = el('circle', { r: 5, style: 'fill: var(--accent-2)' }, svg);
    const readout = document.getElementById('turn-readout');
    const button = document.getElementById('turn-play');
    const slider = document.getElementById('turn-slider');

    const reduced = matchMedia('(prefers-reduced-motion: reduce)').matches;
    let alpha = reduced ? PI / 3 : 0, playing = !reduced, visible = false, hold = 0, last = null, raf = null;

    function draw() {
      scene.renderTurn(alpha);
      const d = alpha * 180 / PI, f = fOf(alpha), g = gOf(alpha);
      slider.value = d.toFixed(1);
      for (const [k, v] of [['x1', cx(d)], ['x2', cx(d)]]) cursor.setAttribute(k, v);
      dotF.setAttribute('cx', cx(d)); dotF.setAttribute('cy', cy(f));
      dotG.setAttribute('cx', cx(d)); dotG.setAttribute('cy', cy(g));
      readout.textContent = `α = ${d.toFixed(1)}°  ·  r = f(α) = ${f.toFixed(4)}  ·  t = g(α) = ${g.toFixed(4)}`;
    }
    function tick(now) {
      raf = null;
      if (!playing || !visible) { last = null; return; }
      if (last !== null) {
        const dt = Math.min(now - last, 100);
        if (alpha >= PI) { hold += dt; if (hold >= 1500) { alpha = 0; hold = 0; } }
        else alpha = Math.min(PI, alpha + dt / 9000 * PI);
      }
      last = now;
      draw();
      raf = requestAnimationFrame(tick);
    }
    const kick = () => { if (raf === null && playing && visible) raf = requestAnimationFrame(tick); };
    const label = () => { button.textContent = playing ? 'Pause' : 'Play'; };

    button.addEventListener('click', () => {
      playing = !playing;
      if (playing && alpha >= PI) { alpha = 0; hold = 0; }
      last = null; label(); kick();
    });
    slider.addEventListener('input', () => {
      playing = false; label();
      alpha = deg(Number(slider.value));
      draw();
    });
    new IntersectionObserver(entries => { visible = entries.some(e => e.isIntersecting); kick(); }, { threshold: 0.1 }).observe(canvas);
    watch(scene, draw);
    label();
    draw();
  }

  // ---------- Figure: f and g with the best ellipse ----------
  function setupChart() {
    const svg = document.getElementById('fg-chart');
    if (!svg) return;
    fgChart(svg, { left: 56, right: 630, y0: 14, y1: 312, vmin: 0.58, vmax: 0.72,
      ticks: [0.58, 0.60, 0.62, 0.64, 0.66, 0.68, 0.70, 0.72], labelY: 332 }, true);
    el('text', { x: 343, y: 354, 'text-anchor': 'middle', 'font-size': 12, style: 'fill: var(--muted)' }, svg, 'α (the hallway has turned by α/2)');
  }

  // ---------- Figure: the final shape ----------
  function setupShape() {
    const canvas = document.getElementById('shape-canvas');
    if (!canvas) return;
    const scene = new Scene(canvas);
    const draw = () => scene.renderShape();
    watch(scene, draw);
    draw();
  }

  if (window.renderMathInElement) {
    window.renderMathInElement(document.body, { delimiters: [{ left: '$$', right: '$$', display: true }], throwOnError: false });
  }
  setupTurn();
  setupChart();
  setupShape();
})();
