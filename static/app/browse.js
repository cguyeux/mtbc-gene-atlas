/* Navigateur génomique, version cliente, pour l'atlas statique.
 *
 * Port fidèle de `_browse_track()` (backend/app.py) : même bornage du span
 * [3000, 400000], même choix de cadre de coordonnées, même géométrie de flèches,
 * mêmes graduations. Le serveur rendait ce SVG à chaque requête ; ici le calcul
 * se fait dans le navigateur à partir de /data/browse_track.json, donc aucune
 * instance à faire vivre et un zoom instantané.
 *
 * Paramètres d'URL identiques à l'ancienne route (`center`, `pos`, `span`), de
 * sorte que les liens des fiches (`/browse?center=Rv0001&span=40000`) marchent tels quels.
 */
(function () {
  "use strict";
  var BASE = (window.ATLAS_BASE || "");
  var W = 1000, DEFAULT_SPAN = 40000;
  var AH = 15, TIP = 7, Y_PLUS = 26, Y_MINUS = 63;
  var COLOR = { requalified: "#2e7d32", family_assigned: "#1565c0", dark: "#8a8a8a" };
  var genes = null;

  function esc(s) {
    return String(s == null ? "" : s).replace(/[&<>"']/g, function (c) {
      return { "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[c];
    });
  }

  /* Reproduit la sélection du serveur : le cadre (MTBC0 vs H37Rv) est celui du gène
     centré. Mélanger les deux dessinerait une fausse synténie (cf. commentaire de
     _browse_track : Rv2512a apparaissait entre Rv2491/Rv2492 alors qu'il en est à ~23 kb). */
  function buildTrack(centerRv, centerPos, span) {
    span = Math.max(3000, Math.min(parseInt(span, 10) || DEFAULT_SPAN, 400000));
    var cg = centerRv ? genes.byRv[centerRv] : null;
    var center = null;
    if (cg && cg.s != null && cg.e != null) center = Math.floor((cg.s + cg.e) / 2);
    if (center === null) center = parseInt(centerPos, 10) || 2200000;
    var half = Math.floor(span / 2);
    var w0 = Math.max(0, center - half), w1 = w0 + span;
    var frame = (cg && cg.f) || "MTBC0";

    var items = genes.rows.filter(function (g) {
      return g.f === frame && g.s <= w1 && g.e >= w0;
    }).sort(function (a, b) { return a.s - b.s; });

    var scale = W / span;
    var drawn = items.map(function (g) {
      var x = Math.max(0, (g.s - w0) * scale);
      var w = Math.max(5, (Math.min(g.e, w1) - Math.max(g.s, w0)) * scale);
      var plus = (g.st || "+") !== "-";
      var y = plus ? Y_PLUS : Y_MINUS;
      var tip = Math.min(TIP, w);
      var pts = plus
        ? [x + "," + y, (x + w - tip).toFixed(1) + "," + y, (x + w).toFixed(1) + "," + (y + AH / 2).toFixed(1),
           (x + w - tip).toFixed(1) + "," + (y + AH), x.toFixed(1) + "," + (y + AH)].join(" ")
        : [(x + w).toFixed(1) + "," + y, (x + tip).toFixed(1) + "," + y, x.toFixed(1) + "," + (y + AH / 2).toFixed(1),
           (x + tip).toFixed(1) + "," + (y + AH), (x + w).toFixed(1) + "," + (y + AH)].join(" ");
      return {
        rv: g.rv, name: g.n || g.rv, color: COLOR[g.v] || "#bdbdbd", pts: pts,
        cx: (x + w / 2).toFixed(1), ty: plus ? y - 4 : y + AH + 11,
        label: w >= 26, product: g.p || "", verdict: g.v || "",
        isCenter: Boolean(centerRv) && g.rv === centerRv,
      };
    });

    var step = Math.max(500, Math.round(span / 6 / 500) * 500);
    var ticks = [];
    for (var t = (Math.floor(w0 / step) + 1) * step; t < w1; t += step) {
      ticks.push({ x: ((t - w0) * scale).toFixed(1), kb: Math.round(t / 1000) });
    }
    return {
      items: drawn, w0: w0, w1: w1, span: span, center: center, frame: frame, ticks: ticks,
      left: Math.max(0, center - span), right: center + span,
      zin: Math.max(3000, Math.floor(span / 2)), zout: Math.min(400000, span * 2),
    };
  }

  /* Port à l'identique de la macro `track_svg` (_track.html) : mêmes hauteurs,
     mêmes contours, mêmes graduations, pour que le rendu statique soit
     indiscernable du rendu serveur. */
  function svgHtml(tr) {
    var H = 120;
    var p = ['<svg viewBox="0 0 ' + W + ' ' + H + '" width="100%" ' +
             'style="min-width:640px; display:block" font-family="system-ui, sans-serif">',
             '<line x1="0" y1="52" x2="' + W + '" y2="52" stroke="#e0e0e0" stroke-width="1"/>',
             '<text x="4" y="12" font-size="9" fill="#aaa">+ strand</text>',
             '<text x="4" y="88" font-size="9" fill="#aaa">− strand</text>'];
    tr.items.forEach(function (g) {
      p.push('<a href="' + BASE + '/gene/' + esc(g.rv) + '">');
      p.push("<title>" + esc(g.name) + " (" + esc(g.rv) + ") — " + esc(g.verdict || "n/a") +
             (g.product ? ": " + esc(g.product) : "") + "</title>");
      p.push('<polygon points="' + g.pts + '" fill="' + g.color + '" stroke="' +
             (g.isCenter ? "#111" : "#ffffff") + '" stroke-width="' + (g.isCenter ? 2 : 0.6) +
             '" opacity="' + (g.isCenter ? 1 : 0.9) + '"/>');
      if (g.label) {
        p.push('<text x="' + g.cx + '" y="' + g.ty + '" text-anchor="middle" font-size="10" fill="' +
               (g.isCenter ? "#111" : "#444") + '" font-weight="' + (g.isCenter ? 600 : 400) + '">' +
               esc(g.name) + "</text>");
      }
      p.push("</a>");
    });
    tr.ticks.forEach(function (t) {
      p.push('<line x1="' + t.x + '" y1="100" x2="' + t.x + '" y2="105" stroke="#bbb"/>');
      p.push('<text x="' + t.x + '" y="116" text-anchor="middle" font-size="9" fill="#999">' +
             String(t.kb).replace(/\B(?=(\d{3})+(?!\d))/g, " ") + " kb</text>");
    });
    p.push("</svg>");
    return p.join("");
  }

  function link(params) {
    var p = new URLSearchParams();
    Object.keys(params).forEach(function (k) { if (params[k] != null && params[k] !== "") p.set(k, params[k]); });
    return BASE + "/browse/?" + p.toString();
  }

  function render() {
    var p = new URLSearchParams(location.search);
    var tr = buildTrack((p.get("center") || "").trim(), p.get("pos"), p.get("span") || DEFAULT_SPAN);

    document.getElementById("bw-svg").innerHTML = svgHtml(tr);
    document.getElementById("bw-meta").innerHTML =
      "window <code>" + tr.w0.toLocaleString("en-US").replace(/,/g, " ") + "</code>–<code>" +
      tr.w1.toLocaleString("en-US").replace(/,/g, " ") + "</code> bp (" +
      Math.round(tr.span / 1000) + " kb) &middot; frame <code>" + esc(tr.frame) + "</code> &middot; " +
      tr.items.length + " gene" + (tr.items.length === 1 ? "" : "s");

    var nav = document.getElementById("bw-nav");
    if (nav) {
      // mêmes libellés et même classe .btn que le gabarit serveur (browse.html)
      nav.innerHTML =
        '<a class="btn" href="' + link({ pos: tr.center - Math.floor(tr.span / 2), span: tr.span }) + '">&#9664; Pan left</a>' +
        '<a class="btn" href="' + link({ pos: tr.center + Math.floor(tr.span / 2), span: tr.span }) + '">Pan right &#9654;</a>' +
        '<a class="btn" href="' + link({ pos: tr.center, span: tr.zin }) + '">&#43; Zoom in</a>' +
        '<a class="btn" href="' + link({ pos: tr.center, span: tr.zout }) + '">&minus; Zoom out</a>';
      // navigation sans rechargement : on garde l'URL partageable
      Array.prototype.forEach.call(nav.querySelectorAll("a"), function (a) {
        a.addEventListener("click", function (ev) {
          ev.preventDefault();
          history.pushState({}, "", a.getAttribute("href"));
          render();
        });
      });
    }
  }

  function boot() {
    fetch(BASE + "/data/browse_track.json")
      .then(function (r) { return r.json(); })
      .then(function (rows) {
        var byRv = {};
        rows.forEach(function (g) { byRv[g.rv] = g; });
        genes = { rows: rows, byRv: byRv };
        render();
        window.addEventListener("popstate", render);
        var f = document.getElementById("bw-jump");
        if (f) {
          f.addEventListener("submit", function (ev) {
            ev.preventDefault();
            var v = f.elements.center ? f.elements.center.value.trim() : "";
            history.pushState({}, "", link(/^\d+$/.test(v) ? { pos: v } : { center: v }));
            render();
          });
        }
      })
      .catch(function (e) {
        var el = document.getElementById("bw-meta");
        if (el) el.textContent = "Could not load the genome track: " + e;
      });
  }

  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", boot);
  else boot();
})();
