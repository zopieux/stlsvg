#!/bin/sh

set -ex

path=$(nix-build | tee | tail -1)
install -Dm644 -t dist "$path/stlsvg.js" "$path/stlsvg.wasm" "$path/index.html"
