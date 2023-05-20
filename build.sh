#!/bin/sh

set -ex

nix build .
install -Dm644 -t dist result/stlsvg.js result/stlsvg.wasm result/index.html
