{
  description = "stlsvg";
  inputs = {
    nixpkgs.url = "github:nixos/nixpkgs?ref=nixos-24.05";
    flake-utils.url = "github:numtide/flake-utils";
  };
  outputs = { self, flake-utils, nixpkgs }: flake-utils.lib.eachDefaultSystem (system:
    let pkgs = import nixpkgs { inherit system; }; in rec
    {
      defaultPackage = pkgs.callPackage ./stlsvg.nix { };
      devShell = pkgs.mkShell {
        inputsFrom = [ defaultPackage ];
      };
    });
}
