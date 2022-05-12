let nixpkgs = import <nixpkgs> {};
in nixpkgs.callPackage ./stlsvg.nix {}
