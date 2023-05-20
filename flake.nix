{
  description = "stlsvg";
  inputs = {
    # More recent nixpkgs makes wasm-ld error out on duplicate symbols in the
    # -D<X>_LIBRARIES=.so overrides. -_-
    nixpkgs.url = "github:nixos/nixpkgs/22.05";
  };
  outputs = { self, nixpkgs }:
    let
      system = "x86_64-linux";
      pkgs = import nixpkgs { inherit system; };
    in
    {
      packages.${system}.stlsvg = pkgs.callPackage ./stlsvg.nix {};
      defaultPackage.${system} = self.packages.${system}.stlsvg;
    };
}
