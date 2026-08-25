#!/usr/bin/env bash
#
# Build R from source into a user-owned prefix. No root, no package manager.
#
# Probes the machine for each of R's mandatory dependencies and builds only the
# ones that are missing or too old, then builds R against them. Everything is
# linked with an RPATH pointing at the prefix, so the resulting R finds its
# libraries without LD_LIBRARY_PATH -- which is what makes it survive batch
# jobs, cron, and `nohup` shells that never source your profile. That is the
# failure mode this script exists to prevent:
#
#     exec/R: error while loading shared libraries: libicuuc.so.66
#
# An env file is still written for convenience, but R does not depend on it.
#
# Usage:
#     ./install_r_local.sh                          # latest R, ~/.local/r-local
#     ./install_r_local.sh --prefix /nobackup/$USER/r-local
#     ./install_r_local.sh --version 4.4.2 --jobs 16
#     ./install_r_local.sh --with-readline --with-icu
#
# The default --prefix is under $HOME, which HPC clusters commonly cap at a
# few GB -- easy to blow through, since building openssl/curl/R leaves object
# files and intermediate artifacts several times larger than the final
# install. If you hit "Disk quota exceeded", point the whole prefix at
# scratch space instead:
#     ./install_r_local.sh --prefix /nobackup/$USER/r-local
# or keep the small final install under $HOME and only redirect the large
# build/ directory (deleted incrementally as each component finishes, so this
# rarely matters once --prefix is already on scratch):
#     ./install_r_local.sh --build-dir /nobackup/$USER/r-build
#
# Re-running is safe: each component drops a stamp file and is skipped on the
# next pass, and anything already installed under $PREFIX is detected by a
# real compile+link probe regardless of stamps -- so a quota failure partway
# through does not lose the components that finished before it hit. --force
# rebuilds everything.
#
# NOTE ON INTEGRITY: tarballs are fetched over HTTPS from each project's
# canonical host. Checksums are deliberately not hardcoded here -- a wrong
# pinned hash is worse than none, since it fails closed on a legitimate
# upstream re-release. Verify by hand if your threat model needs it.

set -euo pipefail

# --------------------------------------------------------------- defaults ---
PREFIX="${HOME}/.local/r-local"
BUILD_DIR=""                       # empty => "${PREFIX}/build"
R_VERSION=""                       # empty => detect the latest from CRAN
JOBS="$( { nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4; } )"
[ "$JOBS" -gt 8 ] && JOBS=8
WITH_READLINE=0                    # off: Rscript does not need a line editor
WITH_ICU=0                         # off: R runs fine without it, and ICU is
                                   #      the library that broke the last R
FORCE=0
KEEP_BUILD=0
NATIVE=0

# Pinned dependency versions. xz is deliberately held at 5.4.x: 5.6.0 and
# 5.6.1 shipped the CVE-2024-3094 backdoor.
ZLIB_V=1.3.1
BZIP2_V=1.0.8
XZ_V=5.4.6
PCRE2_V=10.44
OPENSSL_V=3.0.14
CURL_V=8.8.0
NCURSES_V=6.4
READLINE_V=8.2
ICU_V=74-2
ICU_UNDER=74_2
R_FALLBACK=4.4.2                   # used only if CRAN cannot be reached

usage() {
    sed -n '2,30p' "$0" | sed 's/^# \{0,1\}//'
    exit 0
}

while [ $# -gt 0 ]; do
    case "$1" in
        --prefix)        PREFIX="$2"; shift 2 ;;
        --build-dir)     BUILD_DIR="$2"; shift 2 ;;
        --version)       R_VERSION="$2"; shift 2 ;;
        --jobs|-j)       JOBS="$2"; shift 2 ;;
        --with-readline) WITH_READLINE=1; shift ;;
        --with-icu)      WITH_ICU=1; shift ;;
        --native)        NATIVE=1; shift ;;
        --force)         FORCE=1; shift ;;
        --keep-build)    KEEP_BUILD=1; shift ;;
        -h|--help)       usage ;;
        *) echo "unknown option: $1 (try --help)" >&2; exit 2 ;;
    esac
done

PREFIX="$(mkdir -p "$PREFIX" && cd "$PREFIX" && pwd)"
# Build artifacts (extracted sources, .o files) are far larger than the final
# installed libraries and, unlike $PREFIX, are safe to point at ephemeral
# scratch: nothing outside this script ever reads $BUILD. Defaults alongside
# $PREFIX for backward compatibility, but --build-dir /nobackup/... /or
# pointing --prefix there directly is the fix for a home-quota failure.
BUILD="${BUILD_DIR:-${PREFIX}/build}"
BUILD="$(mkdir -p "$BUILD" && cd "$BUILD" && pwd)"
STAMPS="${BUILD}/stamps"
LOGS="${BUILD}/logs"
mkdir -p "$STAMPS" "$LOGS"

say()  { printf '\n\033[1m==> %s\033[0m\n' "$*"; }
info() { printf '    %s\n' "$*"; }
die()  { printf '\n\033[1;31mERROR: %s\033[0m\n' "$*" >&2; exit 1; }

# ---------------------------------------------------------------- toolchain --
say "Checking the toolchain"
for tool in make tar; do
    command -v "$tool" >/dev/null || die "$tool not found on PATH"
done
CC="${CC:-$(command -v gcc || command -v cc || true)}"
CXX="${CXX:-$(command -v g++ || command -v c++ || true)}"
FC="${FC:-$(command -v gfortran || true)}"
[ -n "$CC" ]  || die "no C compiler found; try 'module load gcc'"
[ -n "$CXX" ] || die "no C++ compiler found; try 'module load gcc'"
[ -n "$FC" ]  || die "no Fortran compiler (gfortran) found. R cannot be built
       without one and bootstrapping GCC here is out of scope.
       Try 'module avail gcc' and load a toolchain that includes gfortran."
export CC CXX FC F77="$FC"
info "CC=$CC"
info "CXX=$CXX"
info "FC=$FC"
info "jobs=$JOBS  prefix=$PREFIX"

DL=""
command -v curl >/dev/null && DL=curl
[ -z "$DL" ] && command -v wget >/dev/null && DL=wget
[ -z "$DL" ] && die "neither curl nor wget found; cannot download sources"

fetch() {  # fetch <url> <output>
    local url="$1" out="$2"
    [ -s "$out" ] && { info "cached $(basename "$out")"; return 0; }
    info "downloading $(basename "$out")"
    if [ "$DL" = curl ]; then
        curl -fsSL --retry 3 -o "$out.part" "$url"
    else
        wget -q -t 3 -O "$out.part" "$url"
    fi
    mv "$out.part" "$out"
}

# ------------------------------------------------------------- build env ----
# Every dependency lands in $PREFIX; the RPATH is what frees the result from
# LD_LIBRARY_PATH at runtime.
export CPPFLAGS="-I${PREFIX}/include ${CPPFLAGS:-}"
export LDFLAGS="-L${PREFIX}/lib -L${PREFIX}/lib64 -Wl,-rpath,${PREFIX}/lib -Wl,-rpath,${PREFIX}/lib64 ${LDFLAGS:-}"
export PKG_CONFIG_PATH="${PREFIX}/lib/pkgconfig:${PREFIX}/lib64/pkgconfig:${PKG_CONFIG_PATH:-}"
export LD_LIBRARY_PATH="${PREFIX}/lib:${PREFIX}/lib64:${LD_LIBRARY_PATH:-}"
export PATH="${PREFIX}/bin:${PATH}"
OPT="-O2"
[ "$NATIVE" -eq 1 ] && OPT="-O2 -march=native"

stamped() { [ "$FORCE" -eq 0 ] && [ -f "${STAMPS}/$1" ]; }
stamp()   { touch "${STAMPS}/$1"; }

# can_link <header> <link-flags> [extra-cflags] -- true if a TU compiles+links
can_link() {
    local hdr="$1" libs="$2" extra="${3:-}"
    local t; t="$(mktemp -d)"
    printf '#include <%s>\nint main(void){return 0;}\n' "$hdr" > "$t/t.c"
    # shellcheck disable=SC2086
    if $CC $extra ${CPPFLAGS} "$t/t.c" -o "$t/t" ${LDFLAGS} $libs >/dev/null 2>&1
    then rm -rf "$t"; return 0; else rm -rf "$t"; return 1; fi
}

# build_autotools <name> <tarball-url> <srcdir-name> [configure args...]
build_autotools() {
    local name="$1" url="$2" dir="$3"; shift 3
    stamped "$name" && { info "$name already built, skipping"; return 0; }
    say "Building $name"
    local tb="${BUILD}/$(basename "$url")"
    fetch "$url" "$tb"
    rm -rf "${BUILD:?}/${dir}"
    tar -xf "$tb" -C "$BUILD"
    (
        cd "${BUILD}/${dir}"
        ./configure --prefix="$PREFIX" CFLAGS="$OPT -fPIC" "$@"
        make -j"$JOBS"
        make install
    ) > "${LOGS}/${name}.log" 2>&1 || {
        tail -30 "${LOGS}/${name}.log" >&2
        die "$name failed to build; full log at ${LOGS}/${name}.log"
    }
    stamp "$name"
    rm -rf "${BUILD:?}/${dir}"
    info "$name installed"
}

# ---------------------------------------------------------- dependencies ----
say "Probing for R's mandatory dependencies"

if can_link zlib.h -lz; then info "zlib      present"; else
    build_autotools zlib "https://zlib.net/fossils/zlib-${ZLIB_V}.tar.gz" \
        "zlib-${ZLIB_V}"
fi

if can_link bzlib.h -lbz2; then info "bzip2     present"; else
    # bzip2 has no configure script, and its shared-library target is awkward
    # to install cleanly. A static -fPIC build links into R perfectly well and
    # removes a runtime dependency entirely.
    if ! stamped bzip2; then
        say "Building bzip2"
        tb="${BUILD}/bzip2-${BZIP2_V}.tar.gz"
        fetch "https://sourceware.org/pub/bzip2/bzip2-${BZIP2_V}.tar.gz" "$tb"
        rm -rf "${BUILD:?}/bzip2-${BZIP2_V}"
        tar -xf "$tb" -C "$BUILD"
        (
            cd "${BUILD}/bzip2-${BZIP2_V}"
            make -j"$JOBS" CFLAGS="$OPT -fPIC -D_FILE_OFFSET_BITS=64"
            make install PREFIX="$PREFIX"
        ) > "${LOGS}/bzip2.log" 2>&1 || {
            tail -30 "${LOGS}/bzip2.log" >&2; die "bzip2 failed; see ${LOGS}/bzip2.log"; }
        stamp bzip2
        rm -rf "${BUILD:?}/bzip2-${BZIP2_V}"
        info "bzip2 installed"
    else info "bzip2 already built, skipping"; fi
fi

if can_link lzma.h -llzma; then info "xz/lzma   present"; else
    build_autotools xz \
        "https://github.com/tukaani-project/xz/releases/download/v${XZ_V}/xz-${XZ_V}.tar.gz" \
        "xz-${XZ_V}" --disable-xzdec --disable-lzmadec --disable-doc
fi

if can_link pcre2.h -lpcre2-8 "-DPCRE2_CODE_UNIT_WIDTH=8"; then
    info "pcre2     present"
else
    build_autotools pcre2 \
        "https://github.com/PCRE2Project/pcre2/releases/download/pcre2-${PCRE2_V}/pcre2-${PCRE2_V}.tar.gz" \
        "pcre2-${PCRE2_V}" --enable-pcre2-8 --enable-unicode --enable-jit
fi

# curl: R requires >= 7.28 and, in practice, working TLS for install.packages().
NEED_CURL=1
if command -v curl-config >/dev/null && can_link curl/curl.h -lcurl; then
    cv="$(curl-config --version | awk '{print $2}')"
    if [ "$(printf '%s\n7.28.0\n' "$cv" | sort -V | head -1)" = "7.28.0" ]; then
        info "libcurl   present ($cv)"; NEED_CURL=0
    else
        info "libcurl   too old ($cv), building"
    fi
fi
if [ "$NEED_CURL" -eq 1 ]; then
    if can_link openssl/ssl.h "-lssl -lcrypto"; then
        info "openssl   present"
    else
        if ! stamped openssl; then
            say "Building openssl"
            tb="${BUILD}/openssl-${OPENSSL_V}.tar.gz"
            fetch "https://www.openssl.org/source/openssl-${OPENSSL_V}.tar.gz" "$tb"
            rm -rf "${BUILD:?}/openssl-${OPENSSL_V}"
            tar -xf "$tb" -C "$BUILD"
            (
                cd "${BUILD}/openssl-${OPENSSL_V}"
                ./config --prefix="$PREFIX" --openssldir="${PREFIX}/ssl" shared
                make -j"$JOBS"
                make install_sw
            ) > "${LOGS}/openssl.log" 2>&1 || {
                tail -30 "${LOGS}/openssl.log" >&2; die "openssl failed; see ${LOGS}/openssl.log"; }
            stamp openssl
            rm -rf "${BUILD:?}/openssl-${OPENSSL_V}"
            info "openssl installed"
        else info "openssl already built, skipping"; fi
    fi
    build_autotools curl "https://curl.se/download/curl-${CURL_V}.tar.gz" \
        "curl-${CURL_V}" --with-openssl --without-libpsl --disable-ldap
fi

if [ "$WITH_READLINE" -eq 1 ]; then
    can_link curses.h -lncurses || build_autotools ncurses \
        "https://ftp.gnu.org/gnu/ncurses/ncurses-${NCURSES_V}.tar.gz" \
        "ncurses-${NCURSES_V}" --with-shared --enable-widec --without-manpages
    can_link readline/readline.h -lreadline || build_autotools readline \
        "https://ftp.gnu.org/gnu/readline/readline-${READLINE_V}.tar.gz" \
        "readline-${READLINE_V}"
fi

if [ "$WITH_ICU" -eq 1 ]; then
    if ! stamped icu; then
        say "Building ICU"
        tb="${BUILD}/icu4c-${ICU_UNDER}-src.tgz"
        fetch "https://github.com/unicode-org/icu/releases/download/release-${ICU_V}/icu4c-${ICU_UNDER}-src.tgz" "$tb"
        rm -rf "${BUILD:?}/icu"
        tar -xf "$tb" -C "$BUILD"
        (
            cd "${BUILD}/icu/source"
            ./configure --prefix="$PREFIX" --disable-tests --disable-samples
            make -j"$JOBS"
            make install
        ) > "${LOGS}/icu.log" 2>&1 || {
            tail -30 "${LOGS}/icu.log" >&2; die "ICU failed; see ${LOGS}/icu.log"; }
        stamp icu
        rm -rf "${BUILD:?}/icu"
        info "ICU installed"
    else info "ICU already built, skipping"; fi
fi

# ------------------------------------------------------------------- R ------
if [ -z "$R_VERSION" ]; then
    say "Resolving the latest R release from CRAN"
    R_VERSION="$(
        { [ "$DL" = curl ] && curl -fsSL https://cran.r-project.org/src/base/R-4/ \
          || wget -qO- https://cran.r-project.org/src/base/R-4/ ; } 2>/dev/null \
        | grep -oE 'R-4\.[0-9]+\.[0-9]+\.tar\.gz' \
        | sed 's/^R-//; s/\.tar\.gz$//' | sort -V | tail -1 || true)"
    if [ -z "$R_VERSION" ]; then
        R_VERSION="$R_FALLBACK"
        info "CRAN unreachable; falling back to R $R_VERSION"
    else
        info "latest is R $R_VERSION"
    fi
fi

R_PREFIX="${PREFIX}/R-${R_VERSION}"
if stamped "R-${R_VERSION}"; then
    info "R ${R_VERSION} already built, skipping"
else
    say "Building R ${R_VERSION}  (this is the long one)"
    tb="${BUILD}/R-${R_VERSION}.tar.gz"
    fetch "https://cran.r-project.org/src/base/R-4/R-${R_VERSION}.tar.gz" "$tb"
    rm -rf "${BUILD:?}/R-${R_VERSION}"
    tar -xf "$tb" -C "$BUILD"

    RCONF=( --prefix="$R_PREFIX" --with-x=no --with-tcltk=no --disable-java )
    # --enable-R-shlib is deliberately NOT set. It routes every call through
    # libR.so and costs a few percent; this R is used for timing measurements.
    [ "$WITH_READLINE" -eq 1 ] || RCONF+=( --with-readline=no )
    [ "$WITH_ICU" -eq 1 ]      || RCONF+=( --with-ICU=no )
    # No --with-blas/--with-lapack: R's bundled reference BLAS is single
    # threaded and deterministic, which is what the benchmark controls assume.
    # An external threaded BLAS would silently reintroduce parallelism.
    (
        cd "${BUILD}/R-${R_VERSION}"
        ./configure "${RCONF[@]}" \
            CFLAGS="$OPT" CXXFLAGS="$OPT" FFLAGS="$OPT" FCFLAGS="$OPT"
        make -j"$JOBS"
        make install
    ) > "${LOGS}/R.log" 2>&1 || {
        tail -40 "${LOGS}/R.log" >&2
        die "R failed to build; full log at ${LOGS}/R.log"
    }
    stamp "R-${R_VERSION}"
    [ "$KEEP_BUILD" -eq 1 ] || rm -rf "${BUILD:?}/R-${R_VERSION}"
fi

# -------------------------------------------------------------- env file ----
ENV_FILE="${PREFIX}/env.sh"
cat > "$ENV_FILE" <<EOF
# Source this to put the locally built R first on PATH.
#     source ${ENV_FILE}
#
# The R binary carries an RPATH to ${PREFIX}/lib, so it resolves its own
# libraries without LD_LIBRARY_PATH. That is set below anyway for any helper
# tool that needs it, but R does not depend on it -- which is why this R keeps
# working inside batch jobs that never read your profile.
export R_LOCAL_PREFIX="${PREFIX}"
export PATH="${R_PREFIX}/bin:${PREFIX}/bin:\${PATH}"
export LD_LIBRARY_PATH="${PREFIX}/lib:${PREFIX}/lib64:\${LD_LIBRARY_PATH:-}"

# Keep user-installed R packages beside this R, not in a version-shared dir.
export R_LIBS_USER="${PREFIX}/library"
EOF
mkdir -p "${PREFIX}/library"

# --------------------------------------------------------------- verify -----
say "Verifying"
export PATH="${R_PREFIX}/bin:${PATH}"
RS="${R_PREFIX}/bin/Rscript"
[ -x "$RS" ] || die "expected $RS to exist and be executable"

# R does not bake an rpath to its OWN lib dir into exec/R -- it relies on the
# bin/R and bin/Rscript wrapper scripts to set LD_LIBRARY_PATH=$R_HOME/lib
# before exec'ing the real binary. libRblas.so and libRlapack.so live there,
# not under $PREFIX/lib where the global rpath above points, so ldd'ing
# exec/R directly reports them "not found" on every build, good or bad --
# that is not a defect, it is how R has always shipped. Match the wrapper's
# environment so the check reflects how R actually runs.
R_LIB="${R_PREFIX}/lib/R/lib"
if LD_LIBRARY_PATH="${R_LIB}:${LD_LIBRARY_PATH}" \
        ldd "${R_PREFIX}/lib/R/bin/exec/R" 2>/dev/null | grep -q "not found"; then
    LD_LIBRARY_PATH="${R_LIB}:${LD_LIBRARY_PATH}" \
        ldd "${R_PREFIX}/lib/R/bin/exec/R" | grep "not found" >&2
    die "the new R still has unresolved libraries (see above)"
fi
info "no unresolved shared libraries"

"$RS" --vanilla -e 'cat(R.version.string, "\n")' \
    || die "Rscript runs but exited non-zero"

# The benchmark drives kmeans through Rscript, so exercise that exact path.
"$RS" --vanilla -e '
  set.seed(1); X <- matrix(rnorm(2000), ncol = 10)
  km <- kmeans(X, centers = X[1:3, ], iter.max = 50, algorithm = "Hartigan-Wong")
  cat("kmeans smoke test OK:", km$iter, "iters, tot.withinss",
      format(km$tot.withinss, digits = 8), "\n")
  cat("libcurl capability:", capabilities("libcurl"), "\n")' \
    || die "kmeans smoke test failed"

# Each component already cleaned its own extracted tree right after
# installing; this only catches --keep-build leftovers and the cached
# tarballs, which are a few hundred MB total, not the multi-GB build trees.
[ "$KEEP_BUILD" -eq 1 ] || find "$BUILD" -mindepth 1 -maxdepth 1 -type d \
    ! -name stamps ! -name logs -exec rm -rf {} + 2>/dev/null || true

say "Done"
cat <<EOF
    R ${R_VERSION} is installed at:
        ${R_PREFIX}

    Use it in this shell:
        source ${ENV_FILE}

    Make it permanent:
        echo 'source ${ENV_FILE}' >> ~/.bashrc

    For a batch job, source it inside the job script -- not only in .bashrc,
    which non-interactive shells skip. The RPATH means R itself will work
    either way; PATH is what needs setting.

    Then, from benchmarks/:
        .venv/bin/python scripts/bench_twitter.py --dry-run
EOF
