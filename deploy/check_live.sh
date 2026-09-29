#!/usr/bin/env bash
# After DNS and TLS: check every public page of https://kinfolk.mdeller.com returns 200.
#   bash deploy/check_live.sh
set -uo pipefail
BASE="${1:-https://kinfolk.mdeller.com}"
fail=0
for p in / /k/ABL1 /resemblance /rooms /black-sheep /dark /almanac /about /api/kinases /api/summary/induced_fit \
         /static/data/tree_sequence.json; do
  code=$(curl -s -o /dev/null -w '%{http_code}' --max-time 20 "$BASE$p")
  printf '%-32s %s\n' "$p" "$code"; [ "$code" = 200 ] || fail=1
done
s=$(curl -s --max-time 20 "$BASE/k/ABL1" | grep -o 'structures/[^"]*cif.gz' | head -1)
[ -n "$s" ] && printf '%-32s %s\n' "structure (gzip)" "$(curl -s -o /dev/null -w '%{http_code}' -H 'Accept-Encoding: gzip' "$BASE/$s")"
exit $fail
