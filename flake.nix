{
  description = "stlsvg";
  inputs = {
    nixpkgs.url = "github:nixos/nixpkgs?ref=203f2ddbe3a48ede1b20b3b86bc8664b311b512d";
    flake-utils.url = "github:numtide/flake-utils";
  };
  outputs = { self, flake-utils, nixpkgs }: flake-utils.lib.eachDefaultSystem (system:
    let pkgs = import nixpkgs { inherit system; }; in rec
    {
      defaultPackage = pkgs.callPackage ./stlsvg.nix { };
      devShell = pkgs.mkShell {
        packages = [
          pkgs.entr
          pkgs.cmake
          pkgs.emscripten
          (pkgs.writeShellScriptBin "configure-stlsvg" ''
            set -ex
            ${pkgs.lib.trim defaultPackage.configureScript} \
              "$@"
          '')
        ];
        inherit (defaultPackage) mygmp mympfr myboost mycgal;
      };
    });
}
