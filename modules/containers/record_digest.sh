#!/bin/bash
# ─── Record what a pulled container actually is ──────────────────────────────
# Usage: record_digest.sh <sif> <image_uri> <apptainer_cachedir>
#        record_digest.sh <sif> <image_uri> --existing
#
# --existing is for a SIF pulled before manifests were introduced: the registry
# is not consulted (a lookup today would describe today's :latest, not the
# build on disk), so IMAGE_DIGEST is recorded as unknown.
#
# Writes <sif>.manifest next to the SIF with the image URI, pull date, the
# registry image digest (the sha256 Docker Hub lists under "Digest" for that
# tag, so a floating tag like :latest is pinned to a known build) and the
# sha256 of the SIF file itself. The manifest is small and is committed to the
# repo, so git history records which build each run used.
#
# The image digest comes from Apptainer's OCI cache: after a pull the image
# lives in <cachedir>/oci-tmp/<digest>/. If that is unavailable (shub://, or a
# cache that was already cleaned) a HEAD request to the registry is tried.

set -uo pipefail

SIF="$1"
IMAGE="$2"
CACHEDIR="${3:-}"
MANIFEST="${SIF}.manifest"

digest=""
note=""

if [ "$CACHEDIR" = "--existing" ]; then
    CACHEDIR=""
    note="pre-existing SIF; digest not recorded at pull time. Delete the SIF and rerun the pipeline to re-pull and pin."
fi

# 1. Apptainer OCI cache directory name
if [ -n "$CACHEDIR" ] && [ -d "${CACHEDIR}/oci-tmp" ]; then
    digest=$(ls -1 "${CACHEDIR}/oci-tmp" 2>/dev/null | grep -E '^[0-9a-f]{64}$' | head -1)
    [ -n "$digest" ] && digest="sha256:${digest}"
fi

# 2. Registry HEAD request (Docker Hub and gcr.io/ghcr.io style registries)
if [ -z "$digest" ] && [ -z "$note" ] && [[ "$IMAGE" == docker://* ]]; then
    ref="${IMAGE#docker://}"
    tag="${ref##*:}"; repo="${ref%:*}"
    [ "$tag" = "$repo" ] && tag="latest"
    case "$repo" in
        */*/*) registry="${repo%%/*}"; path="${repo#*/}" ;;      # us.gcr.io/x/y
        *.*/*) registry="${repo%%/*}"; path="${repo#*/}" ;;      # ghcr.io/x
        */*)   registry="registry-1.docker.io"; path="$repo" ;;  # satijalab/seurat
        *)     registry="registry-1.docker.io"; path="library/$repo" ;;
    esac
    accept="application/vnd.docker.distribution.manifest.list.v2+json,application/vnd.oci.image.index.v1+json,application/vnd.docker.distribution.manifest.v2+json,application/vnd.oci.image.manifest.v1+json"
    auth=""
    if [ "$registry" = "registry-1.docker.io" ]; then
        token=$(curl -fsSL "https://auth.docker.io/token?service=registry.docker.io&scope=repository:${path}:pull" 2>/dev/null | sed -n 's/.*"token":"\([^"]*\)".*/\1/p')
        [ -n "$token" ] && auth="Authorization: Bearer ${token}"
    fi
    digest=$(curl -fsSI -H "Accept: ${accept}" ${auth:+-H "$auth"} "https://${registry}/v2/${path}/manifests/${tag}" 2>/dev/null \
             | tr -d '\r' | awk -F': ' 'tolower($1)=="docker-content-digest"{print $2}')
fi

sif_sha=$(sha256sum "$SIF" | cut -d' ' -f1)
short="${digest#sha256:}"; short="${short:0:12}"

hub_url=""
if [ -n "$digest" ] && [[ "$IMAGE" == docker://* ]] && [[ "$IMAGE" != docker://*.*/* ]]; then
    ref="${IMAGE#docker://}"; tag="${ref##*:}"; repo="${ref%:*}"
    [ "$tag" = "$repo" ] && tag="latest"
    hub_url="https://hub.docker.com/layers/${repo}/${tag}/images/sha256-${digest#sha256:}"
fi

cat > "$MANIFEST" <<EOF
IMAGE=${IMAGE}
PULLED=$(date '+%Y-%m-%d %H:%M:%S')
IMAGE_DIGEST=${digest:-unknown}
IMAGE_DIGEST_SHORT=${short:-unknown}
SIF=$(basename "$SIF")
SIF_SHA256=${sif_sha}
SIF_SIZE=$(du -h "$SIF" | cut -f1)
HUB_URL=${hub_url}
NOTE=${note}
EOF

echo "[INFO] Container manifest written to ${MANIFEST}"
cat "$MANIFEST"
[ -n "$digest" ] || [ -n "$note" ] || echo "[WARNING] Could not determine image digest for ${IMAGE}"
exit 0
