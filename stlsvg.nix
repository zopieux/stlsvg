{ emscriptenStdenv, gmp, mpfr, cgal_5, cmake }:

let
  stdenv = emscriptenStdenv;
  buildPhase = ''
    emmake make -j$(nproc)
  '';
  NIX_CFLAGS_COMPILE = "";
  EM_CACHE = "/tmp/devnull";

  mygmp = (gmp.override { inherit stdenv; }).overrideDerivation (old: {
    doCheck = false;
    dontStrip = true;
    enableParallelBuilding = true;
    inherit NIX_CFLAGS_COMPILE EM_CACHE buildPhase;
    # Jesus fucking Christ.
    postPatch = ''
      substituteInPlace configure \
        --replace 'as_fn_set_status $1' 'echo set_status $1' \
        --replace '  exit $1' '  echo exit $1'
    '';
    configurePhase = ''
      HOME=$TMPDIR
      ABI=x32 emconfigure ./configure --disable-assembly --disable-fat --prefix=$out
    '';
  });

  mympfr = (mpfr.override { inherit stdenv; }).overrideDerivation (old: {
    propagatedBuildInputs = [ mygmp ];
    doCheck = false;
    dontStrip = true;
    enableParallelBuilding = true;
    outputs = [ "out" "dev" ];
    inherit NIX_CFLAGS_COMPILE EM_CACHE buildPhase;
    configurePhase = ''
      HOME=$TMPDIR
      EMCONFIGURE_JS=1 emconfigure ./configure --with-gmp-include=${mygmp.dev}/include --with-gmp-lib=${mygmp}/lib --prefix=$out
    '';
  });

  myboost = builtins.fetchTarball {
    url = "http://downloads.sourceforge.net/project/boost/boost/1.79.0/boost_1_79_0.tar.bz2";
    sha256 = "080fr3y6xyb02k7zdq647rzmsfxzic47yjzqj2kvmqhgkpsj42m1";
  };

  mycgal = (cgal_5.override { inherit stdenv; }).overrideDerivation (old: {
    propagatedBuildInputs = [ mygmp mympfr ];
    nativeBuildInputs = [ cmake ];
    doCheck = false;
    dontStrip = true;
    enableParallelBuilding = true;
    inherit NIX_CFLAGS_COMPILE EM_CACHE;
    prePatch = ''
      substituteInPlace cmake/modules/CGAL_SetupCGALDependencies.cmake \
        --replace 'find_package(Threads REQUIRED)' ' ' \
        --replace 'target_link_libraries(''${target} INTERFACE Threads::Threads)' ' '
    '';
    configurePhase = ''
      HOME=$TMPDIR
      emcmake cmake . $cmakeFlags \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_INSTALL_PREFIX=$out \
        -DGMP_INCLUDE_DIR:STRING=${mygmp.dev}/include -DGMP_LIBRARIES:STRING=${mygmp}/lib/libgmp.so \
        -DMPFR_INCLUDE_DIR:STRING=${mympfr.dev}/include -DMPFR_LIBRARIES:STRING=${mympfr}/lib/libmpfr.so \
        -DBoost_INCLUDE_DIR:STRING=${myboost} \
        -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON \
        -DCGAL_HAS_NO_THREADS=ON \
        -DWITH_CGAL_ImageIO=OFF -DWITH_CGAL_Qt5:BOOL=OFF \
        -DCMAKE_CXX_FLAGS="-v -U__SSE2_MATH__ --ignore-dynamic-linking -DCGAL_HAS_NO_THREADS -U__GNUG__ -DCGAL_NO_ASSERTIONS -DCGAL_FORWARD -DBOOST_MATH_DISABLE_STD_FPCLASSIFY -DBOOST_NO_NATIVE_LONG_DOUBLE_FP_CLASSIFY -DBOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS"
    '';
  });

in
stdenv.mkDerivation rec {
  name = "stlsvg";
  src = ./src;
  enableParallelBuilding = true;
  propagatedBuildInputs = [ mycgal mygmp mympfr ];
  nativeBuildInputs = [ cmake ];
  dontStrip = true;
  EM_CACHE = "/build/.cache";
  configurePhase = ''
    HOME=$TMPDIR
    emcmake cmake . $cmakeFlags \
      -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_INSTALL_PREFIX=$out \
      -DGMP_INCLUDE_DIR:STRING=${mygmp.dev}/include -DGMP_LIBRARIES:STRING=${mygmp}/lib/libgmp.so \
      -DMPFR_INCLUDE_DIR:STRING=${mympfr.dev}/include -DMPFR_LIBRARIES:STRING=${mympfr}/lib/libmpfr.so \
      -DBoost_INCLUDE_DIR:STRING=${myboost}
  '';
  checkPhase = "";
  installPhase = ''
    install -Dm644 -t "$out" stlsvg.js stlsvg.wasm
    install -Dm644 stlsvg.html "$out/index.html"
  '';
}
