# shellcheck shell=bash
# fetch-pinned.sh: the one download rule of the asset builders in
# scripts/ci/ (build-callability.sh, build-cn-intervals.sh,
# build-sv-rule-regions.sh). Source it; it defines two functions.
#
#   fetch_pinned DIR NAME URL SHA256
#     Makes DIR/NAME the file at URL with that sha256. A file already there
#     with the right sha256 is kept, so DIR can be a cache that outlives the
#     run. Otherwise it prints "Downloading URL" on stderr (the e2e cases
#     count these lines), downloads to DIR/NAME.part and renames it when the
#     sha256 matches. A 403 or 404 fails at once: a missing or forbidden file
#     stays so. Any other failure is tried again 6 times, 30 s apart, for
#     about five minutes: EBI's FTP server refused every connection for
#     minutes at a time on 2026-11-06 (curl exit 7). Returns 1 with an ERROR
#     line when the file cannot be had.
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
  local dir=$1 name=$2 url=$3 want=$4 file code got
  file="${dir}/${name}"
  if [ -f "$file" ] && [ "$(pinned_sha256 "$file")" = "$want" ]; then
    return 0
  fi
  mkdir -p "$dir" || return 1
  echo "Downloading ${url}" >&2
  if ! code=$(curl -fsSL -w '%{http_code}' -o "${file}.part" "$url"); then
    case $code in
      403 | 404)
        rm -f "${file}.part"
        echo "ERROR: cannot download ${url}: HTTP ${code}" >&2
        return 1
        ;;
    esac
    if ! curl -fsSL --retry 6 --retry-all-errors --retry-delay 30 --retry-max-time 300 \
        -o "${file}.part" "$url"; then
      rm -f "${file}.part"
      echo "ERROR: cannot download ${url}" >&2
      return 1
    fi
  fi
  got=$(pinned_sha256 "${file}.part")
  if [ "$got" != "$want" ]; then
    rm -f "${file}.part"
    echo "ERROR: ${name} has sha256 ${got}, want ${want}" >&2
    return 1
  fi
  mv "${file}.part" "$file"
}
