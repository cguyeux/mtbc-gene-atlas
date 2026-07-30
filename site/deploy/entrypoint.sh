#!/bin/sh
# Serverless entrypoint. The catalogue DB is baked into the image at build time
# (read-only at /app/data/db.sqlite). We copy it to the writable, ephemeral
# /tmp so the app can also accept feedback writes within a running instance.
# Feedback is NOT persistent (lost on scale-to-zero) -- by design for a
# read-mostly atlas; route it to email/Managed DB later if persistence matters.
set -e
cp -n /app/data/db.sqlite /tmp/app.sqlite 2>/dev/null || true
exec uvicorn backend.app:app --host 0.0.0.0 --port "${PORT:-8080}"
