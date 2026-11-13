# shellcheck shell=bash
# fetch-pinned.sh: the one download rule of the asset builders in
# scripts/ci/ (build-callability.sh, build-cn-intervals.sh,
# build-sv-rule-regions.sh). Source it; it defines two functions.
#
#   fetch_pinned DIR NAME URL SHA256
#     Makes DIR/NAME the file at URL with that sha256. A file already there
#     with the right sha256 is kept, so DIR can be a cache that outlives the
#     run. Otherwise it prints "Downloading URL" on stderr (the e2e cases
#     count these lines), downloads to a temporary file of its own beside
#     DIR/NAME (so two runs on one cache never share one) and renames it
#     when the sha256 matches. A 403 or 404 fails at once: a missing or
#     forbidden file stays so. Any other failure is tried again 6 times,
#     30 s apart, for about five minutes: EBI's FTP server refused every
#     connection for minutes at a time on 2026-11-06 (curl exit 7). A
#     connection that does not open in 30 s, or a transfer under 1 KB/s for
#     60 s, counts as a failure, so a stalled server cannot hold the build.
#     Returns 1 with an ERROR line when the file cannot be had.
#
#   pinned_sha256 FILE
#     The sha256 of FILE, from sha256sum or shasum.
#
# og never downloads: the builders fetch here, then hand the folder to
# og dev build.

pinned_sha256() {
  if command -v sha256sum > /dev/null 2>&1; then
    sha256sum "$1" | cut -d' ' -f1
  else
    shasum -a 256 "$1" | cut -d' ' -f1
  fi
}

fetch_pinned() {
  local dir=$1 name=$2 url=$3 want=$4 file part code got
  local attempt=0 started=$SECONDS
  file="${dir}/${name}"
  if [ -f "$file" ] && [ "$(pinned_sha256 "$file")" = "$want" ]; then
    return 0
  fi
  mkdir -p "$dir" || return 1
  part=$(mktemp "${file}.part.XXXXXX") || return 1
  echo "Downloading ${url}" >&2
  while ! code=$(curl -fsSL --connect-timeout 30 --speed-limit 1024 --speed-time 60 \
      -w '%{http_code}' -o "$part" "$url"); do
    case $code in
      403 | 404)
        rm -f "$part"
        echo "ERROR: cannot download ${url}: HTTP ${code}" >&2
        return 1
        ;;
    esac
    attempt=$((attempt + 1))
    # Seven attempts at most. Like curl's retry-max-time, this bounds when
    # another attempt starts; an active transfer keeps its stall limits.
    if [ "$attempt" -ge 7 ] || [ "$((SECONDS - started + 30))" -ge 300 ]; then
      rm -f "$part"
      echo "ERROR: cannot download ${url}" >&2
      return 1
    fi
    sleep 30
  done
  got=$(pinned_sha256 "$part")
  if [ "$got" != "$want" ]; then
    rm -f "$part"
    echo "ERROR: ${name} has sha256 ${got}, want ${want}" >&2
    return 1
  fi
  chmod 644 "$part"
  mv "$part" "$file"
}
