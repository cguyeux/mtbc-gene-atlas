#!/usr/bin/env bash
# Déploiement du site annotation_mtbc sur Scaleway Serverless Containers.
# Cristallise le cycle manuel en un outil réutilisable, avec les pièges déjà rencontrés intégrés :
#   - sous-commande `redeploy` (PAS `deploy`, inexistante dans le scw actuel) ;
#   - `-o json` pour parser la sortie (tabulée par défaut) ;
#   - ATTENDRE la fin du push avant l'update (latence d'indexation du registry, sinon "image not found") ;
#   - poll du statut jusqu'à `ready`, en couvrant l'état d'erreur (silence != succès).
# L'ingest (garde-fou de cohérence P6.2 + toutes les couches) tourne DANS le build (Dockerfile.serverless).
# Usage:  deploy/deploy.sh v42
set -euo pipefail

TAG="${1:?usage: deploy.sh vN  (ex: deploy.sh v42)}"
REGION=fr-par
REGISTRY=rg.fr-par.scw.cloud/mtbc/gene-atlas
CONTAINER=af2030f6-a098-4549-a292-1b9565dd61ae      # container 'mtbc'
SITE_DIR="$(cd "$(dirname "$0")/.." && pwd)"        # .../site

cd "$SITE_DIR"

echo "[1/5] build image ($TAG)…"
docker build --network=host -f deploy/Dockerfile.serverless -t mtbc-gene-atlas:latest .

echo "[2/5] tag + push $TAG (bloque jusqu'à la fin du push)…"
docker tag mtbc-gene-atlas:latest "$REGISTRY:$TAG"
docker push "$REGISTRY:$TAG"

# attend que le container repasse `ready` (poll). Couvre l'état error. Renvoie 0=ready, 1=error, 2=timeout.
wait_ready() {
  for i in $(seq 1 60); do
    st=$(scw container container get "$CONTAINER" region="$REGION" -o json 2>/dev/null \
          | python3 -c "import sys,json;print(json.load(sys.stdin).get('status',''))" 2>/dev/null || true)
    echo "   poll $i: status=$st"
    [ "$st" = "ready" ] && return 0
    printf '%s' "$st" | grep -qi error && return 1
    sleep 8
  done
  return 2
}

echo "[3/5] update container -> $TAG…"
scw container container update "$CONTAINER" image="$REGISTRY:$TAG" region="$REGION" -o json >/dev/null
# `update image=` fait basculer le container en état transitoire `updating` ET tire déjà la nouvelle
# image ; il faut ATTENDRE `ready` avant tout redeploy, sinon « transient state error » (bug v43).
echo "   attente stabilisation post-update…"
wait_ready || { echo "post-update: état non-ready"; exit 1; }

echo "[4/5] redeploy (redeploy, pas deploy ; sur container ready) …"
scw container container redeploy "$CONTAINER" region="$REGION" -o json >/dev/null || echo "   (redeploy non nécessaire : l'update a déjà déployé)"

echo "[5/5] attente ready finale…"
if wait_ready; then
  echo "DEPLOYED $TAG (ready) — https://mtbc.gclab.fr"; exit 0
else
  rc=$?
  [ "$rc" -eq 1 ] && { echo "DEPLOY ERROR"; exit 1; }
  echo "TIMEOUT waiting for ready"; exit 2
fi
