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
          pkgs.watchexec
          pkgs.cmake
          pkgs.emscripten
          pkgs.cgal
          pkgs.gmp
          pkgs.mpfr
          pkgs.boost
          (pkgs.writeShellScriptBin "configure-stlsvg" ''
            set -ex
            # Configure Native Build
            cmake -B build-native -S src -DCMAKE_BUILD_TYPE=Debug -DCGAL_DO_NOT_WARN_ABOUT_CMAKE_BUILD_TYPE=TRUE

            # Configure Emscripten Build
            ${pkgs.lib.trim defaultPackage.configureScript} \
              -B build \
              -S src \
              "$@"
          '')
          (pkgs.writeShellScriptBin "watch-build" ''
            set -e
            watchexec -w src -- "cmake --build build-native & ; cmake --build build & ; wait"
          '')
        ];
        inherit (defaultPackage) mygmp mympfr myboost mycgal;
      };
    });
}
