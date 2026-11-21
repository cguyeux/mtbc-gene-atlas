/* Vue cliente de /genes pour la version statique de l'atlas.
 *
 * Reproduit ce que faisait la route FastAPI `/genes` : recherche plein texte,
 * filtres (verdict, hypothetical, enriched, to review), pagination de 50.
 * Lit les mêmes paramètres d'URL que l'ancienne route, donc TOUS les liens
 * existants restent valides (`/genes?hypo=1`, `/genes?q=katG&page=2`, ainsi que
 * les formulaires GET de base.html et home.html).
 *
 * Source de données : /data/genes_index.json, produit par bin/build_static.py.
 */
(function () {
  "use strict";
  var BASE = (window.ATLAS_BASE || "");
  var PER_PAGE = 50;
  var VERDICT_LABEL = { requalified: "Resolved", family_assigned: "Family assigned", dark: "Still unknown" };

  function qs() { return new URLSearchParams(location.search); }
  function esc(s) {
    return String(s == null ? "" : s).replace(/[&<>"']/g, function (c) {
      return { "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[c];
    });
  }

  /* Recherche par POSITION génomique, reprise de _parse_coord() côté serveur :
     accepte 2367322, 1.5Mb, 1500kb. Le cadre de coordonnées est ambigu (MTBC0 pour
     les loci ancrés, H37Rv pour les non ancrés) : comme le serveur, on renvoie les
     correspondances des deux cadres plutôt que d'en masquer une partie. */
  function parseCoord(q) {
    var s = String(q).trim().toLowerCase().replace(/,/g, ".").replace(/\s/g, "");
    var m;
    if ((m = s.match(/^([\d.]+)mb$/))) return Math.round(parseFloat(m[1]) * 1e6);
    if ((m = s.match(/^([\d.]+)kb$/))) return Math.round(parseFloat(m[1]) * 1e3);
    if ((m = s.match(/^(\d+)bp$/))) return parseInt(m[1], 10);
    if (/^\d+$/.test(s)) return parseInt(s, 10);
    return null;
  }

  function filterRows(rows, p) {
    var q = (p.get("q") || "").trim().toLowerCase();
    var verdict = p.get("verdict") || "";
    var hypo = p.get("hypo") === "1";
    var enriched = p.get("enriched") === "1";
    var review = p.get("review") === "1";
    var pos = q ? parseCoord(q) : null;

    return rows.filter(function (g) {
      if (verdict && g.v !== verdict) return false;
      if (hypo && !g.h) return false;
      if (enriched && !g.e) return false;
      if (review && !g.w) return false;
      if (!q) return true;
      if (pos !== null && g.s != null && g.e2 != null && g.s <= pos && g.e2 >= pos) return true;
      return (g.rv + " " + g.n + " " + g.m + " " + g.l + " " + g.r).toLowerCase().indexOf(q) !== -1;
    });
  }

  function rowHtml(g) {
    var pfam = (g.d || []).map(function (d) { return "<code>" + esc(d) + "</code>"; }).join(" ");
    var badge = g.v
      ? '<span class="badge badge-' + esc(g.v) + '">' + esc(VERDICT_LABEL[g.v] || g.v) + "</span>"
      : (g.e ? '<span class="badge">enriched</span>' : "");
    return '<tr' + (g.h ? ' class="hypo"' : "") + ">" +
      '<td><a href="' + BASE + '/gene/' + esc(g.rv) + '">' + esc(g.n || g.rv) + "</a><br><small>" + esc(g.rv) + "</small></td>" +
      "<td><code>" + esc(g.m) + "</code></td>" +
      '<td class="muted">' + esc(g.l) + "</td>" +
      "<td>" + esc(g.r || "—") + "</td>" +
      "<td>" + pfam + "</td>" +
      "<td>" + badge + "</td></tr>";
  }

  function pagerHtml(page, nPages, p) {
    if (nPages <= 1) return "";
    function link(n, label) {
      var q2 = new URLSearchParams(p); q2.set("page", n);
      return '<a href="' + BASE + '/genes/?' + q2.toString() + '">' + label + "</a>";
    }
    return (page > 1 ? link(page - 1, "&larr; prev") : "") +
      "<span>page " + page + " / " + nPages + "</span>" +
      (page < nPages ? link(page + 1, "next &rarr;") : "");
  }

  function render(rows) {
    var p = qs();
    var matched = filterRows(rows, p);
    var page = Math.max(1, parseInt(p.get("page") || "1", 10) || 1);
    var nPages = Math.max(1, Math.ceil(matched.length / PER_PAGE));
    if (page > nPages) page = nPages;
    var slice = matched.slice((page - 1) * PER_PAGE, page * PER_PAGE);

    document.getElementById("gl-count").textContent =
      matched.length + " gene" + (matched.length === 1 ? "" : "s") + " matched.";
    document.getElementById("gl-rows").innerHTML = slice.map(rowHtml).join("");
    document.getElementById("gl-pager").innerHTML = pagerHtml(page, nPages, p);

    // refléter l'état courant dans le formulaire, comme le faisait le gabarit serveur
    var f = document.getElementById("gl-filters");
    if (f) {
      f.q.value = p.get("q") || "";
      f.verdict.value = p.get("verdict") || "";
      f.hypo.checked = p.get("hypo") === "1";
      f.enriched.checked = p.get("enriched") === "1";
      f.review.checked = p.get("review") === "1";
    }
  }

  function boot() {
    fetch(BASE + "/data/genes_index.json")
      .then(function (r) { return r.json(); })
      .then(function (rows) {
        render(rows);
        // navigation arrière/avant sans rechargement
        window.addEventListener("popstate", function () { render(rows); });
        var f = document.getElementById("gl-filters");
        if (f) {
          f.addEventListener("submit", function (ev) {
            ev.preventDefault();
            var p = new URLSearchParams();
            if (f.q.value.trim()) p.set("q", f.q.value.trim());
            if (f.verdict.value) p.set("verdict", f.verdict.value);
            if (f.hypo.checked) p.set("hypo", "1");
            if (f.enriched.checked) p.set("enriched", "1");
            if (f.review.checked) p.set("review", "1");
            history.pushState({}, "", BASE + "/genes/" + (p.toString() ? "?" + p.toString() : ""));
            render(rows);
          });
        }
      })
      .catch(function (e) {
        document.getElementById("gl-count").textContent = "Could not load the gene index: " + e;
      });
  }

  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", boot);
  else boot();
})();
