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
      devShell = with pkgs; mkShell {
        nativeBuildInputs = [ cmake ];
        buildInputs = [
          cgal_5
          mpfr
          gmp
          defaultPackage.myboost
        ];
        packages = [ entr ];
        # inputsFrom = [ defaultPackage ];
        # buildInputs = with pkgs; [ cmake ];
        # # propagatedBuildInputs = with defaultPackage; [ mycgal mygmp mympfr ];
        # BOOST_ROOT = "${defaultPackage.myboost}";
        # GMPDIR = "${defaultPackage.mygmp}";
      };
    });
}
