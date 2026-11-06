#!/usr/bin/env bash
# Clone and build AlVin and CoxIter at pinned commits (fresh, independent of any
# local checkouts elsewhere). Sources/builds live in tools/src (gitignored).
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
SRC="$HERE/src"
mkdir -p "$SRC"

ALVIN_COMMIT="${ALVIN_COMMIT:-}"
COXITER_COMMIT="${COXITER_COMMIT:-}"

clone() {  # name url commit
  local name=$1 url=$2 commit=$3
  if [ ! -d "$SRC/$name/.git" ]; then
    git clone --recursive "$url" "$SRC/$name"
  fi
  if [ -n "$commit" ]; then
    git -C "$SRC/$name" checkout -q "$commit"
    git -C "$SRC/$name" submodule update --init --recursive
  fi
  git -C "$SRC/$name" rev-parse HEAD > "$HERE/$name.commit"
}

clone CoxIter https://github.com/rgugliel/CoxIter.git "$COXITER_COMMIT"
clone AlVin   https://github.com/rgugliel/AlVin.git   "$ALVIN_COMMIT"

# Build-only port of AlVin to igraph 1.0 / C++17 / arm64 (see patches/port_alvin.py).
git -C "$SRC/AlVin" checkout -q -- .
python3 "$HERE/patches/port_alvin.py" "$SRC/AlVin"

# Apple CLT clang on this SDK cannot find libc++ headers; use Homebrew LLVM.
LLVM=/opt/homebrew/opt/llvm/bin
CMAKE_ARGS=(-DCMAKE_BUILD_TYPE=Release
            -DCMAKE_CXX_FLAGS=-I/opt/homebrew/include
            -DCMAKE_EXE_LINKER_FLAGS=-L/opt/homebrew/lib)
# AlVin: OpenMP disabled -- openblas (via igraph) loads a second libomp and aborts.
# This only makes AlVin single-threaded; results are unaffected.
for name in CoxIter AlVin; do
  extra=()
  [ "$name" = AlVin ] && extra=(-DCMAKE_DISABLE_FIND_PACKAGE_OpenMP=TRUE)
  rm -rf "$SRC/$name/build"; mkdir -p "$SRC/$name/build/output"
  (cd "$SRC/$name/build" && CXX=$LLVM/clang++ CC=$LLVM/clang \
     cmake "${CMAKE_ARGS[@]}" ${extra[@]+"${extra[@]}"} .. >/dev/null && make -j8 >/dev/null)
done
echo "AlVin   @ $(cat "$HERE/AlVin.commit")"
echo "CoxIter @ $(cat "$HERE/CoxIter.commit")"
