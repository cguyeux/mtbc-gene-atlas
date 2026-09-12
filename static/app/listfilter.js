/* Filtrage client générique des listes secondaires (/rnas, /microproteins) de
 * l'atlas statique.
 *
 * Ces deux vues étaient filtrées côté serveur par des `select` et des cases à
 * cocher adossés à des colonnes SQL. En statique, on ne peut pas reproduire ces
 * filtres à l'identique sans réimplémenter leur sémantique, et un filtre qui a
 * l'air de marcher sans marcher est pire qu'un filtre absent. On les remplace
 * donc par UNE recherche libre, qui filtre les lignes réellement présentes dans
 * la page — la liste étant figée en entier, aucune donnée n'est hors de portée.
 */
(function () {
  "use strict";

  function boot() {
    var input = document.getElementById("lf-q");
    var table = document.querySelector("table.genes, table.table");
    if (!input || !table) return;
    var body = table.tBodies[0];
    if (!body) return;
    var rows = Array.prototype.slice.call(body.rows);
    var counter = document.getElementById("lf-count");
    var total = rows.length;

    // texte de chaque ligne, calculé une fois
    var haystack = rows.map(function (tr) { return tr.textContent.toLowerCase(); });

    function apply() {
      var q = input.value.trim().toLowerCase();
      var shown = 0;
      for (var i = 0; i < rows.length; i++) {
        var ok = !q || haystack[i].indexOf(q) !== -1;
        rows[i].hidden = !ok;
        if (ok) shown++;
      }
      if (counter) {
        counter.textContent = q
          ? shown + " of " + total + " shown"
          : total + " entr" + (total === 1 ? "y" : "ies");
      }
    }

    input.addEventListener("input", apply);
    apply();
  }

  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", boot);
  else boot();
})();
