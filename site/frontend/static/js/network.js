/* Interactive STRING network (P12.4): lightweight force-directed layout, drag, channel/score
   filters, hover highlight, 2-hop. No external dependency. Progressive enhancement over the
   server-rendered static SVG. Data is read from a <script type="application/json" id="net-data">. */
(function () {
  "use strict";
  var dataEl = document.getElementById("net-data");
  var host = document.getElementById("net-interactive");
  var staticEl = document.getElementById("net-static");
  if (!dataEl || !host) return;
  var G;
  try { G = JSON.parse(dataEl.textContent); } catch (e) { return; }
  if (!G || !G.nodes || G.nodes.length < 2) return;

  var W = G.width, H = G.height, SVGNS = "http://www.w3.org/2000/svg";
  var cx = W / 2, cy = H / 2;

  // model
  var nodes = G.nodes.map(function (n) { return Object.assign({}, n, { vx: 0, vy: 0, fixed: n.is_center }); });
  var byId = {}; nodes.forEach(function (n) { n.x = n.x || cx; n.y = n.y || cy; byId[n.rv] = n; });
  var edges = G.edges.map(function (e) { return { s: byId[e.source], t: byId[e.target], color: e.color, w: e.w, op: e.op, kind: e.kind, score: e.score, vis: true }; })
                     .filter(function (e) { return e.s && e.t; });
  nodes.forEach(function (n) { n.deg = 0; }); edges.forEach(function (e) { e.s.deg++; e.t.deg++; });

  // svg scaffold
  var svg = document.createElementNS(SVGNS, "svg");
  svg.setAttribute("viewBox", "0 0 " + W + " " + H);
  svg.setAttribute("width", "100%");
  svg.style.minWidth = "620px"; svg.style.display = "block"; svg.style.touchAction = "none";
  var gEdges = document.createElementNS(SVGNS, "g");
  var gNodes = document.createElementNS(SVGNS, "g");
  svg.appendChild(gEdges); svg.appendChild(gNodes);

  edges.forEach(function (e) {
    var l = document.createElementNS(SVGNS, "line");
    l.setAttribute("stroke", e.color); l.setAttribute("stroke-width", e.w); l.setAttribute("stroke-opacity", e.op);
    e.el = l; gEdges.appendChild(l);
  });
  nodes.forEach(function (n) {
    var a = document.createElementNS(SVGNS, "a"); a.setAttributeNS("http://www.w3.org/1999/xlink", "href", n.href); a.setAttribute("href", n.href);
    var c = document.createElementNS(SVGNS, "circle");
    c.setAttribute("r", n.r); c.setAttribute("fill", n.color);
    c.setAttribute("stroke", n.is_center ? "#111" : "#fff"); c.setAttribute("stroke-width", n.is_center ? 2.5 : 1);
    if (n.hypothetical) c.setAttribute("stroke-dasharray", "2.5,2");
    c.style.cursor = "pointer";
    var ttl = document.createElementNS(SVGNS, "title"); ttl.textContent = n.title; c.appendChild(ttl);
    var t = document.createElementNS(SVGNS, "text");
    t.setAttribute("text-anchor", n.is_center ? "middle" : "start"); t.setAttribute("font-size", n.is_center ? 11 : 10);
    t.setAttribute("fill", n.is_center ? "#fff" : "#333"); t.setAttribute("font-weight", n.is_center ? 700 : 400);
    t.setAttribute("pointer-events", "none"); t.textContent = n.is_center ? n.name.slice(0, 9) : n.name;
    a.appendChild(c); a.appendChild(t); gNodes.appendChild(a);
    n.el = c; n.txt = t; n.link = a;
    bindDrag(n, a, c);
    bindHover(n, c);
  });

  host.appendChild(svg);
  host.hidden = false;
  if (staticEl) staticEl.hidden = true;

  // ---- force simulation ----
  var alpha = 1, EDGE_LEN = Math.min(W, H) / 6, REPULSE = 5200;
  function tick() {
    for (var i = 0; i < nodes.length; i++) {
      var a = nodes[i];
      for (var j = i + 1; j < nodes.length; j++) {
        var b = nodes[j];
        var dx = a.x - b.x, dy = a.y - b.y, d2 = dx * dx + dy * dy || 1, d = Math.sqrt(d2);
        var f = REPULSE / d2, ux = dx / d, uy = dy / d;
        a.vx += ux * f; a.vy += uy * f; b.vx -= ux * f; b.vy -= uy * f;
      }
    }
    edges.forEach(function (e) {
      if (!e.vis) return;
      var dx = e.t.x - e.s.x, dy = e.t.y - e.s.y, d = Math.sqrt(dx * dx + dy * dy) || 1;
      var f = (d - EDGE_LEN) * 0.02, ux = dx / d, uy = dy / d;
      e.s.vx += ux * f; e.s.vy += uy * f; e.t.vx -= ux * f; e.t.vy -= uy * f;
    });
    nodes.forEach(function (n) {
      n.vx += (cx - n.x) * 0.008; n.vy += (cy - n.y) * 0.008;   // centering
      if (n.fixed || n.drag) { n.vx = n.vy = 0; return; }
      n.vx *= 0.82; n.vy *= 0.82;
      n.x += Math.max(-25, Math.min(25, n.vx * alpha));
      n.y += Math.max(-25, Math.min(25, n.vy * alpha));
      n.x = Math.max(n.r, Math.min(W - n.r, n.x)); n.y = Math.max(n.r, Math.min(H - n.r, n.y));
    });
    render();
    alpha *= 0.985;
    if (alpha > 0.03) requestAnimationFrame(tick);
  }
  function render() {
    edges.forEach(function (e) {
      e.el.setAttribute("x1", e.s.x.toFixed(1)); e.el.setAttribute("y1", e.s.y.toFixed(1));
      e.el.setAttribute("x2", e.t.x.toFixed(1)); e.el.setAttribute("y2", e.t.y.toFixed(1));
    });
    nodes.forEach(function (n) {
      n.el.setAttribute("cx", n.x.toFixed(1)); n.el.setAttribute("cy", n.y.toFixed(1));
      n.txt.setAttribute("x", (n.is_center ? n.x : n.x + n.r + 3).toFixed(1));
      n.txt.setAttribute("y", (n.y + 4).toFixed(1));
      if (!n.is_center) n.txt.setAttribute("text-anchor", "start");
    });
  }
  function reheat() { alpha = Math.max(alpha, 0.5); requestAnimationFrame(tick); }

  // ---- drag ----
  function svgPt(evt) {
    var p = svg.createSVGPoint(); p.x = evt.clientX; p.y = evt.clientY;
    return p.matrixTransform(svg.getScreenCTM().inverse());
  }
  function bindDrag(n, a, c) {
    var moved = false;
    c.addEventListener("pointerdown", function (evt) {
      evt.preventDefault(); moved = false; n.drag = true; c.setPointerCapture(evt.pointerId);
      var p0 = svgPt(evt), ox = n.x - p0.x, oy = n.y - p0.y;
      function mv(e2) { var p = svgPt(e2); n.x = p.x + ox; n.y = p.y + oy; moved = true; render(); }
      function up(e2) {
        n.drag = false; c.releasePointerCapture(evt.pointerId);
        c.removeEventListener("pointermove", mv); c.removeEventListener("pointerup", up);
        if (moved) reheat();
      }
      c.addEventListener("pointermove", mv); c.addEventListener("pointerup", up);
    });
    a.addEventListener("click", function (e) { if (moved) { e.preventDefault(); moved = false; } });
  }

  // ---- hover highlight ----
  function bindHover(n, c) {
    c.addEventListener("mouseenter", function () {
      var keep = {}; keep[n.rv] = 1;
      edges.forEach(function (e) { if (e.s === n) keep[e.t.rv] = 1; if (e.t === n) keep[e.s.rv] = 1; });
      nodes.forEach(function (m) { m.el.style.opacity = keep[m.rv] ? 1 : 0.15; m.txt.style.opacity = keep[m.rv] ? 1 : 0.15; });
      edges.forEach(function (e) { e.el.style.opacity = (e.s === n || e.t === n) && e.vis ? 0.95 : 0.05; });
    });
    c.addEventListener("mouseleave", function () {
      nodes.forEach(function (m) { m.el.style.opacity = 1; m.txt.style.opacity = 1; });
      edges.forEach(function (e) { e.el.style.opacity = e.vis ? e.op : 0; });
    });
  }

  // ---- filters (channel + score) ----
  function applyFilters() {
    var chans = {};
    document.querySelectorAll(".net-chan").forEach(function (ch) { chans[ch.value] = ch.checked; });
    var minScore = parseInt((document.getElementById("net-score") || {}).value || "0", 10);
    edges.forEach(function (e) {
      e.vis = (chans[e.kind] !== false) && e.score >= minScore;
      e.el.style.display = e.vis ? "" : "none";
    });
    var connected = {}; nodes.forEach(function (m) { connected[m.rv] = m.is_center; });
    edges.forEach(function (e) { if (e.vis) { connected[e.s.rv] = 1; connected[e.t.rv] = 1; } });
    nodes.forEach(function (m) {
      var on = connected[m.rv] || m.hop === 1;   // toujours montrer le centre et le 1er cercle
      m.link.style.display = on ? "" : "none";
    });
    reheat();
  }
  document.querySelectorAll(".net-chan").forEach(function (ch) { ch.addEventListener("change", applyFilters); });
  var sc = document.getElementById("net-score");
  if (sc) { sc.addEventListener("input", function () { var o = document.getElementById("net-score-val"); if (o) o.textContent = sc.value; applyFilters(); }); }

  requestAnimationFrame(tick);
})();
