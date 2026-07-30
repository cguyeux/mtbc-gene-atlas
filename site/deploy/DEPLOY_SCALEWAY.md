# Déploiement sur Scaleway Serverless Containers (SQLite embarqué)

Ce guide déploie le site comme une image Docker autonome (base SQLite construite
dans l'image au moment du build). Pas de base externe, scale-to-zero, coût quasi
nul à l'arrêt. Le formulaire de feedback fonctionne mais n'est pas persistant
(remis à zéro à froid) ; c'est volontaire pour un atlas en lecture.

Le déploiement se fait en deux outils :
- la console web Scaleway (https://console.scaleway.com) pour créer les ressources,
- Docker en ligne de commande pour pousser l'image (déjà construite et testée).

Région utilisée dans les exemples : `fr-par`. Adapter si besoin.

## Prérequis (une seule fois)

1. Un compte Scaleway et un Projet (par défaut « default »).
2. Une clé API : console > IAM > API Keys > Generate API key. On obtient un
   Access Key (commence par `SCW...`) et un Secret Key (un UUID). Garder le
   Secret Key, il n'est affiché qu'une fois.
3. Docker installé localement (déjà le cas ici).

## Étape 1 - L'image (déjà faite)

L'image est construite et testée en local sous le tag `mtbc-gene-atlas:latest`.
Pour la (re)construire :

```
cd site && docker build --network=host -f deploy/Dockerfile.serverless -t mtbc-gene-atlas:latest .
```

Test local rapide :

```
docker run -d --network host -e PORT=8085 --name mtbc-test mtbc-gene-atlas:latest
```

Ouvrir http://127.0.0.1:8085 puis arrêter avec `docker rm -f mtbc-test`.

## Étape 2 - Créer un Container Registry

Console > Container Registry > Create namespace. Nom au choix, par exemple `mtbc`,
région `fr-par`, visibilité Private. L'endpoint sera `rg.fr-par.scw.cloud/mtbc`.

## Étape 3 - Connecter Docker au registry

Le nom d'utilisateur est littéralement `nologin`, le mot de passe est le Secret Key.

```
docker login rg.fr-par.scw.cloud -u nologin --password-stdin
```

La commande attend le mot de passe sur l'entrée standard : coller le Secret Key
puis Entrée (ou `echo "<SECRET_KEY>" | docker login rg.fr-par.scw.cloud -u nologin --password-stdin`).

## Étape 4 - Taguer et pousser l'image

Remplacer `mtbc` par le nom du namespace créé à l'étape 2.

```
docker tag mtbc-gene-atlas:latest rg.fr-par.scw.cloud/mtbc/gene-atlas:latest
```
```
docker push rg.fr-par.scw.cloud/mtbc/gene-atlas:latest
```

## Étape 5 - Créer le Serverless Container

Console > Serverless > Containers > Create namespace (par exemple `mtbc`), puis
Create container :
- Image : choisir `rg.fr-par.scw.cloud/mtbc/gene-atlas:latest` (le registry du même compte est proposé automatiquement).
- Port : `8080`.
- Ressources : 256 Mo de mémoire, 140 mvCPU (largement suffisant ; 128/70 marche aussi).
- Scale : min 0 (scale-to-zero), max 1 (ou 2).
- Déployer.

Scaleway fournit alors une URL HTTPS du type
`https://gene-atlas-xxxx.functions.fnc.fr-par.scw.cloud`. C'est l'adresse publique du site.

## Étape 6 - Mettre à jour le contenu plus tard

Après un nouveau batch d'annotation :

```
cd .. && python analyses/phase4_export_json.py data/pilot_batch_01.txt && cp data/gene_xref.tsv site/content/gene_xref.tsv
```

Puis reconstruire, repousser, redéployer :

```
cd site && docker build --network=host -f deploy/Dockerfile.serverless -t mtbc-gene-atlas:latest . && docker tag mtbc-gene-atlas:latest rg.fr-par.scw.cloud/mtbc/gene-atlas:latest && docker push rg.fr-par.scw.cloud/mtbc/gene-atlas:latest
```

IMPORTANT - piège du tag `latest` : Scaleway met l'image en cache par tag et ne
la re-télécharge PAS au redéploiement si le tag est identique. Utiliser un tag
versionné à chaque mise à jour (`v2`, `v3`, ...) et pointer le container dessus,
sinon le redeploy continue de servir l'ancienne image.

```
docker tag mtbc-gene-atlas:latest rg.fr-par.scw.cloud/mtbc/gene-atlas:v2 && docker push rg.fr-par.scw.cloud/mtbc/gene-atlas:v2
```

Puis dans la console : Settings du container > Image > choisir le tag `v2` > Deploy.

## Notes

- Feedback non persistant : remis à zéro au scale-to-zero. Pour le conserver,
  router vers un email (variables type `*_FEEDBACK_MAIL_*` comme dans atlas_mtbc)
  ou brancher une Managed Database PostgreSQL.
- Coût : à l'arrêt, proche de zéro ; on paie le temps de calcul par requête et le
  stockage de l'image dans le registry.
- Domaine personnalisé : la config du container permet d'ajouter un domaine custom
  (CNAME) si vous en avez un.

## Annexe - Tout en CLI scw (alternative à la console)

```
scw init
```
```
scw registry namespace create name=mtbc region=fr-par
```
```
scw container namespace create name=mtbc region=fr-par
```

Récupérer l'ID du namespace container (`scw container namespace list`), puis :

```
scw container container create namespace-id=<NAMESPACE_ID> name=gene-atlas registry-image=rg.fr-par.scw.cloud/mtbc/gene-atlas:latest port=8080 min-scale=0 max-scale=1 memory-limit=256 cpu-limit=140 region=fr-par
```
```
scw container container redeploy <CONTAINER_ID> region=fr-par
```

`scw container container get <CONTAINER_ID>` affiche l'URL (`DomainName`).

### Mise à jour d'un container existant (cycle de batch)

**Raccourci (recommandé) : tout le cycle est scripté dans `deploy/deploy.sh`.** Une seule commande
fait build → tag → push → update → redeploy → attente `ready`, avec les pièges déjà intégrés :

```
deploy/deploy.sh v42
```

À la main, après `docker push ...:vN`, deux commandes suffisent (container id obtenu via
`scw container container list region=fr-par`) :

```
scw container container update <CONTAINER_ID> image=rg.fr-par.scw.cloud/mtbc/gene-atlas:vN region=fr-par
```
```
scw container container redeploy <CONTAINER_ID> region=fr-par
```

Pièges rencontrés : la sous-commande est `redeploy` (et NON `deploy`, qui n'existe pas dans le
scw actuel → affiche l'aide sans rien faire ; corrigé 2026-07-04 lors du déploiement v41) ;
**`update image=…` bascule le container en état transitoire `updating` ET tire déjà la nouvelle
image (l'update SUFFIT à déployer) : un `redeploy` LANCÉ IMMÉDIATEMENT après échoue en
« transient state error : container is in a transient state 'updating' » (rencontré au v43).
→ ATTENDRE le retour à `ready` après l'update avant tout redeploy ; `deploy.sh` encapsule cette
attente (fonction `wait_ready`) et tolère l'échec du redeploy (redondant après update).** ;
l'argument est `image=` (et non `registry-image=`) ; sortie tabulée par défaut, ajouter `-o json`
pour parser ; `scw`
exige un `default_organization_id` dans `~/.config/scw/config.yaml` (chez Scaleway,
l'org et le projet "default" partagent souvent le même UUID ; à défaut, le
récupérer via `GET https://api.scaleway.com/account/v3/projects/<project_id>`).
