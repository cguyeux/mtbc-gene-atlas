/* Réseau d'interactions STRING (ego-network), version cliente, pour l'atlas statique.
 *
 * Port de `_network_graph()` + de la macro `network_svg` (backend/app.py, _network.html) :
 * même sélection BFS des partenaires, même placement circulaire, mêmes couleurs de
 * canal, mêmes épaisseurs d'arête. Les données viennent de /data/network.json
 * (top 24 partenaires par gène avec leurs canaux : ~0,3 Mo gzippé), ce qui est très
 * inférieur aux ~353 Mo qu'aurait coûté le figement d'un graphe par gène et par
 * profondeur — et cela garde le re-centrage instantané.
 *
 * Paramètres d'URL identiques à l'ancienne route (`center`, `hops`, `top`).
 */
(function () {
  "use strict";
  var BASE = (window.ATLAS_BASE || "");
  var NET_W = 900, NET_H = 560, MAX_NODES = 46;
  var VERDICT_COLOR = { requalified: "#2e7d32", family_assigned: "#1565c0", dark: "#8a8a8a" };
  var DB = null; // { partners: {rv: [{rv,s,ch,cd}]}, meta: {rv: {n,v,h,p}} }

  function esc(s) {
    return String(s == null ? "" : s).replace(/[&<>"']/g, function (c) {
      return { "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[c];
    });
  }

  function topPartners(rv, k) {
    var ps = (DB.partners[rv] || []).slice();
    ps.sort(function (a, b) { return (b.s || 0) - (a.s || 0); });
    return ps.slice(0, k);
  }

  /* Port de _edge_kind() : le canal qui domine décide de la couleur. */
  function edgeKind(p) {
    var ch = p.ch || {};
    if ((ch.experimental || 0) >= 400 || (ch.database || 0) >= 400) return ["experimental", "#2e7d32"];
    if (p.cd || (ch.neighborhood || 0) >= 400 || (ch.cooccurence || 0) >= 400 || (ch.fusion || 0) >= 400)
      return ["context", "#ef6c00"];
    return ["other", "#9e9e9e"];
  }

  function buildGraph(centerRv, top, hops) {
    if (!DB.meta[centerRv] || !DB.partners[centerRv]) return null;
    top = Math.max(3, Math.min(top || 12, 24));
    hops = hops >= 2 ? 2 : 1;

    var hop = {}, order = [];
    hop[centerRv] = 0; order.push(centerRv);
    var h1 = topPartners(centerRv, top);
    h1.forEach(function (p) {
      if (!(p.rv in hop)) { hop[p.rv] = 1; order.push(p.rv); }
    });
    if (hops === 2) {
      h1.forEach(function (p) {
        topPartners(p.rv, 6).forEach(function (q) {
          if (!(q.rv in hop) && Object.keys(hop).length < MAX_NODES) {
            hop[q.rv] = 2; order.push(q.rv);
          }
        });
      });
    }
    order = order.filter(function (rv) { return DB.meta[rv]; });
    var nodeset = {};
    order.forEach(function (rv) { nodeset[rv] = true; });

    var cx = NET_W / 2, cy = NET_H / 2;
    var R1 = Math.min(NET_W, NET_H) / 2 - 130, R2 = Math.min(NET_W, NET_H) / 2 - 40;
    var h1list = order.filter(function (rv) { return hop[rv] === 1; });
    var h2list = order.filter(function (rv) { return hop[rv] === 2; });
    var pos = {}; pos[centerRv] = [cx, cy];
    h1list.forEach(function (rv, i) {
      var a = 2 * Math.PI * i / Math.max(1, h1list.length) - Math.PI / 2;
      pos[rv] = [cx + R1 * Math.cos(a), cy + R1 * Math.sin(a)];
    });
    h2list.forEach(function (rv, i) {
      var a = 2 * Math.PI * i / Math.max(1, h2list.length) - Math.PI / 2;
      pos[rv] = [cx + R2 * Math.cos(a), cy + R2 * Math.sin(a)];
    });

    var nodes = order.map(function (rv) {
      var m = DB.meta[rv], x = pos[rv][0], y = pos[rv][1], isC = rv === centerRv;
      return {
        rv: rv, name: m.n || rv, x: +x.toFixed(1), y: +y.toFixed(1),
        r: isC ? 22 : (m.h ? 12 : 13),
        color: VERDICT_COLOR[m.v] || "#bdbdbd", isCenter: isC, hypothetical: Boolean(m.h),
        href: isC ? BASE + "/gene/" + rv : BASE + "/network/?center=" + rv,
        tx: +(x + (x >= cx ? 16 : -16)).toFixed(1), anchor: x >= cx ? "start" : "end",
        title: (m.n || rv) + " (" + (m.v || "n/a") + ") — " + (m.p || ""),
      };
    });

    var edges = [], seen = {};
    order.forEach(function (a) {
      (DB.partners[a] || []).forEach(function (q) {
        var b = q.rv;
        if (!nodeset[b] || a === b) return;
        var key = [a, b].sort().join("|");
        if (seen[key]) return;
        seen[key] = true;
        var kc = edgeKind(q), sc = q.s || 0;
        edges.push({
          color: kc[1], x1: +pos[a][0].toFixed(1), y1: +pos[a][1].toFixed(1),
          x2: +pos[b][0].toFixed(1), y2: +pos[b][1].toFixed(1),
          w: +(0.8 + 3.0 * sc / 1000).toFixed(2), op: +(0.22 + 0.5 * sc / 1000).toFixed(2),
        });
      });
    });
    return { nodes: nodes, edges: edges, centerRv: centerRv, nPartners: h1list.length,
             nNodes: nodes.length, hops: hops, top: top };
  }

  function svgHtml(g) {
    var p = ['<svg viewBox="0 0 ' + NET_W + " " + NET_H + '" width="100%" ' +
             'style="min-width:620px; display:block" font-family="system-ui, sans-serif">'];
    g.edges.forEach(function (e) {
      p.push('<line x1="' + e.x1 + '" y1="' + e.y1 + '" x2="' + e.x2 + '" y2="' + e.y2 +
             '" stroke="' + e.color + '" stroke-width="' + e.w + '" stroke-opacity="' + e.op + '"/>');
    });
    g.nodes.forEach(function (n) {
      p.push('<a href="' + esc(n.href) + '"><title>' + esc(n.title) + "</title>");
      p.push('<circle cx="' + n.x + '" cy="' + n.y + '" r="' + n.r + '" fill="' + n.color +
             '" stroke="' + (n.isCenter ? "#111" : "#ffffff") + '" stroke-width="' +
             (n.isCenter ? 2.5 : 1) + '"' + (n.hypothetical ? ' stroke-dasharray="2.5,2"' : "") + "/>");
      if (n.isCenter) {
        p.push('<text x="' + n.x + '" y="' + (n.y + 4) + '" text-anchor="middle" font-size="11" ' +
               'font-weight="700" fill="#ffffff">' + esc(n.name.slice(0, 9)) + "</text>");
      } else {
        p.push('<text x="' + n.tx + '" y="' + (n.y + 4) + '" text-anchor="' + n.anchor +
               '" font-size="10" fill="#333">' + esc(n.name) + "</text>");
      }
      p.push("</a>");
    });
    p.push("</svg>");
    return p.join("");
  }

  function render() {
    var q = new URLSearchParams(location.search);
    var center = (q.get("center") || "").trim();
    var hops = parseInt(q.get("hops") || "1", 10);
    var top = parseInt(q.get("top") || "12", 10);
    var host = document.getElementById("nw-graph");
    var meta = document.getElementById("nw-meta");
    if (!host) return;

    if (!center) {
      host.innerHTML = "";
      if (meta) meta.innerHTML = '<span class="muted">Enter a gene above to see its STRING ego-network.</span>';
      return;
    }
    var g = buildGraph(center, top, hops);
    if (!g) {
      host.innerHTML = "";
      if (meta) meta.innerHTML = '<span class="muted">No STRING network for <code>' +
        esc(center) + "</code>.</span>";
      return;
    }
    host.innerHTML = svgHtml(g);
    if (meta) {
      var other = g.hops >= 2 ? 1 : 2;
      meta.innerHTML = "<code>" + esc(g.centerRv) + "</code> &middot; " + g.nPartners +
        " direct partner" + (g.nPartners === 1 ? "" : "s") + " &middot; " + g.nNodes +
        " nodes at " + g.hops + " hop" + (g.hops === 1 ? "" : "s") +
        ' &middot; <a href="' + BASE + '/network/?center=' + esc(g.centerRv) + "&hops=" + other + '">' +
        "show " + other + " hop" + (other === 1 ? "" : "s") + "</a>";
    }
  }

  function boot() {
    fetch(BASE + "/data/network.json")
      .then(function (r) { return r.json(); })
      .then(function (db) {
        DB = db;
        render();
        window.addEventListener("popstate", render);
        document.addEventListener("click", function (ev) {
          var a = ev.target.closest ? ev.target.closest('a[href*="/network/?"]') : null;
          if (!a) return;
          ev.preventDefault();
          history.pushState({}, "", a.getAttribute("href"));
          render();
        });
        var f = document.getElementById("nw-jump");
        if (f) {
          f.addEventListener("submit", function (ev) {
            ev.preventDefault();
            var v = (f.elements.center ? f.elements.center.value : "").trim();
            history.pushState({}, "", BASE + "/network/?center=" + encodeURIComponent(v));
            render();
          });
        }
      })
      .catch(function (e) {
        var meta = document.getElementById("nw-meta");
        if (meta) meta.textContent = "Could not load the network data: " + e;
      });
  }

  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", boot);
  else boot();
})();
