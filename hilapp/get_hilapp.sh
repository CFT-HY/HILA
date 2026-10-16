#!/usr/bin/env bash
#
# get_hilapp.sh - download a prebuilt, statically linked hilapp binary matching this
# hila source tree, and install it to hilapp/bin/hilapp.
#
# Binaries are built by the GitHub workflow .github/workflows/hilapp-binaries.yml for
# Linux x86_64 and aarch64, for each change of hilapp/ on the main branch.  The release tag
# is  hilapp-s<signature>-<tree>, where <signature> is HILA_SIGNATURE_NUMBER and <tree> the
# git tree hash of the hilapp/ directory.  Thus the binary corresponds exactly to the
# checked-out hilapp source.
#
# Usage:  hilapp/get_hilapp.sh [options]
#   -c, --compatible   if there is no binary for this exact source version, use the newest
#                      binary with the same signature number
#   -f, --force        download even if the installed binary is already this version
#   -r, --repo <o/r>   GitHub repository (default CFT-HY/HILA, or $HILAPP_REPO)
#   -h, --help         this help
#
# For testing, the download locations can be changed with $HILAPP_DOWNLOAD_URL (default
# https://github.com/<repo>/releases/download) and $HILAPP_API_URL (default
# https://api.github.com/repos/<repo>/releases?per_page=100).
#
# Needs git, curl or wget, tar and sha256sum (or shasum).  After "git pull" run this again;
# it does nothing if the binary is up to date.

set -euo pipefail

repo="${HILAPP_REPO:-CFT-HY/HILA}"
compatible=0
force=0

usage() {
    sed -n '3,/^$/s/^# \{0,1\}//p' "$0"
    exit "${1:-0}"
}

while [ $# -gt 0 ]; do
    case "$1" in
    -c | --compatible) compatible=1 ;;
    -f | --force) force=1 ;;
    -r | --repo)
        shift
        repo="${1:?--repo needs owner/name}"
        ;;
    -h | --help) usage 0 ;;
    *)
        echo "get_hilapp.sh: unknown option '$1'" >&2
        usage 1
        ;;
    esac
    shift
done

error() {
    echo "get_hilapp.sh: ERROR: $*" >&2
    exit 1
}

hilapp_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
hila_dir="$(dirname "${hilapp_dir}")"
bin_dir="${hilapp_dir}/bin"

# architecture
case "$(uname -s)-$(uname -m)" in
Linux-x86_64 | Linux-amd64) arch=x86_64 ;;
Linux-aarch64 | Linux-arm64) arch=aarch64 ;;
*) error "no prebuilt hilapp for $(uname -s) $(uname -m), build it from source (see hilapp/Makefile)" ;;
esac

# source version: signature number and git tree hash of hilapp/
sig=$(awk '/^#define HILA_SIGNATURE_NUMBER/ {print $3}' "${hila_dir}/libraries/hila_signatures.h")
[ -n "${sig}" ] || error "could not read HILA_SIGNATURE_NUMBER from libraries/hila_signatures.h"

command -v git > /dev/null || error "git is needed to identify the hilapp source version"
tree=$(git -C "${hila_dir}" rev-parse --short=12 HEAD:hilapp 2> /dev/null) ||
    error "could not get the git version of ${hilapp_dir} (not a git checkout?)"

tag="hilapp-s${sig}-${tree}"

if ! git -C "${hila_dir}" diff --quiet HEAD -- hilapp/src hilapp/Makefile libraries/hila_signatures.h; then
    echo "get_hilapp.sh: WARNING: hilapp/src has local changes, which are not in the prebuilt binary" >&2
fi

# downloader
if command -v curl > /dev/null; then
    fetch() { curl -fsSL --retry 3 -o "$2" "$1"; }
    fetch_stdout() { curl -fsSL --retry 3 "$1"; }
elif command -v wget > /dev/null; then
    fetch() { wget -q -O "$2" "$1"; }
    fetch_stdout() { wget -q -O - "$1"; }
else
    error "curl or wget is needed"
fi

if command -v sha256sum > /dev/null; then
    sha_check() { sha256sum -c --status "$1"; }
elif command -v shasum > /dev/null; then
    sha_check() { shasum -a 256 -c --status "$1"; }
else
    error "sha256sum or shasum is needed"
fi

download_url="${HILAPP_DOWNLOAD_URL:-https://github.com/${repo}/releases/download}"
api_url="${HILAPP_API_URL:-https://api.github.com/repos/${repo}/releases?per_page=100}"

release_file="${bin_dir}/.hilapp_release"

asset="hilapp-linux-${arch}.tar.gz"
tmp=$(mktemp -d)
trap 'rm -rf "${tmp}"' EXIT

download() {
    # $1 = tag; returns non-zero if the release or asset does not exist
    local url="${download_url}/$1/${asset}"
    fetch "${url}" "${tmp}/${asset}" 2> /dev/null && fetch "${url}.sha256" "${tmp}/${asset}.sha256"
}

if [ "${force}" = 0 ] && [ -x "${bin_dir}/hilapp" ] && [ -f "${release_file}" ] &&
    [ "$(cat "${release_file}")" = "${tag}" ]; then
    echo "hilapp is up to date (${tag})"
    exit 0
fi

echo "Downloading ${tag} for linux-${arch} from ${download_url}"

if ! download "${tag}"; then
    if [ "${compatible}" = 0 ]; then
        cat >&2 << EOF
get_hilapp.sh: no prebuilt hilapp for this source version (${tag}).
  Binaries are built for the main branch of ${repo}; this checkout may be on another
  branch, have unpushed commits, or the build may not have finished yet.  Options:
    hilapp/get_hilapp.sh --compatible    newest binary with the same signature ${sig}
    make -C hilapp install               build from source (needs llvm/clang development files)
EOF
        exit 1
    fi

    # newest release with the same signature (releases are listed newest first)
    compat=$(fetch_stdout "${api_url}" |
        grep -o "\"tag_name\": *\"hilapp-s${sig}-[0-9a-f]*\"" | head -1 | sed 's/.*"\(hilapp-[^"]*\)"/\1/') || true
    [ -n "${compat}" ] || error "no prebuilt hilapp with signature ${sig} found"
    [ -f "${release_file}" ] && [ "$(cat "${release_file}")" = "${compat}" ] && [ "${force}" = 0 ] && {
        echo "hilapp is up to date with the newest compatible binary (${compat})"
        exit 0
    }
    echo "get_hilapp.sh: WARNING: using ${compat}, which is not built from exactly this hilapp source" >&2
    download "${compat}" || error "download of ${compat} failed"
    tag="${compat}"
fi

(cd "${tmp}" && sha_check "${asset}.sha256") || error "checksum of ${asset} does not match"

mkdir -p "${tmp}/x"
tar -C "${tmp}/x" -xzf "${tmp}/${asset}"
"${tmp}/x/hilapp" --help > /dev/null 2>&1 || error "downloaded hilapp does not run on this system"

mkdir -p "${bin_dir}"
if [ -e "${bin_dir}/hilapp" ]; then
    mv -f "${bin_dir}/hilapp" "${bin_dir}/hilapp.previous"
fi
mv "${tmp}/x/hilapp" "${bin_dir}/hilapp"
echo "${tag}" > "${release_file}"

echo "Installed ${bin_dir}/hilapp (${tag})"
sed 's/^/    /' "${tmp}/x/BUILD_INFO"
