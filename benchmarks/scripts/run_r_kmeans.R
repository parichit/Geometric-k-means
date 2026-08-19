#!/usr/bin/env Rscript
# One timed R `stats::kmeans` fit, driven by bench_twitter.py.
#
# Invoked once per (k, seed) so that R's fits interleave with the Python ones
# and all three implementations see the same thermal state.
#
# The dataset arrives as a raw little-endian float64 blob in COLUMN-MAJOR
# order, written by the Python side after it has already dropped non-numeric
# columns and applied --subsample. Going through a raw blob rather than the
# original CSV means R clusters bit-identical data to scikit-learn and
# geokmeans, and it costs a readBin instead of a multi-hundred-MB read.csv.
#
# Initial centroids come from the same inits/<dataset>_k<K>_seed<S>.csv every
# other implementation reads, at %.18e precision.
#
# Energy is sampled from RAPL *inside* this process, around the kmeans() call
# only, so R's startup and data load are outside the measured region exactly
# as they are for the Python implementations.
#
# Results are written as key=value lines to --out; stdout carries only logs.

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i)) {
    if (is.null(default)) stop("missing required argument: ", flag)
    return(default)
  }
  args[[i + 1L]]
}

data_path <- get_arg("--data")
n         <- as.integer(get_arg("--n"))
d         <- as.integer(get_arg("--d"))
init_path <- get_arg("--init")
k         <- as.integer(get_arg("--k"))
max_iter  <- as.integer(get_arg("--max-iter", "300"))
algorithm <- get_arg("--algorithm", "Hartigan-Wong")
out_path  <- get_arg("--out")

# ---------------------------------------------------------------- RAPL ----
# Mirrors scripts/energy.py: package-* and dram are summed, core/uncore/psys
# are not (they would double-count), and deltas are taken modulo the counter
# range because energy_uj is free-running and wraps.
rapl_domains <- function() {
  root <- Sys.getenv("GEOKMEANS_RAPL_ROOT", "/sys/class/powercap")
  if (!dir.exists(root)) return(list())
  doms <- list()
  pkgs <- list.dirs(root, recursive = FALSE, full.names = TRUE)
  pkgs <- pkgs[grepl("/intel-rapl:[0-9]+$", pkgs)]
  add <- function(doms, path) {
    nm <- try(trimws(readLines(file.path(path, "name"), warn = FALSE)[1]),
              silent = TRUE)
    mx <- try(as.numeric(readLines(file.path(path, "max_energy_range_uj"),
                                   warn = FALSE)[1]), silent = TRUE)
    ok <- try(as.numeric(readLines(file.path(path, "energy_uj"),
                                   warn = FALSE)[1]), silent = TRUE)
    if (inherits(nm, "try-error") || inherits(mx, "try-error") ||
        inherits(ok, "try-error")) return(doms)
    kind <- if (startsWith(nm, "package")) "package"
            else if (nm == "dram") "dram" else "other"
    c(doms, list(list(path = path, name = nm, max = mx, kind = kind)))
  }
  for (p in sort(pkgs)) {
    doms <- add(doms, p)
    subs <- list.dirs(p, recursive = FALSE, full.names = TRUE)
    for (s in sort(subs[grepl("/intel-rapl:[0-9]+:[0-9]+$", subs)])) {
      doms <- add(doms, s)
    }
  }
  doms
}

read_rapl <- function(doms) {
  vapply(doms, function(x)
    as.numeric(readLines(file.path(x$path, "energy_uj"), warn = FALSE)[1]),
    numeric(1))
}

DOMS <- rapl_domains()

# ---------------------------------------------------------------- data ----
con <- file(data_path, "rb")
# column-major on the wire, so matrix() consumes it without a transpose
X <- matrix(readBin(con, "double", n = as.numeric(n) * as.numeric(d),
                    size = 8, endian = "little"), nrow = n, ncol = d)
close(con)

C0 <- as.matrix(read.csv(init_path, header = FALSE))
dimnames(C0) <- NULL
stopifnot(nrow(C0) == k, ncol(C0) == d)

cat(sprintf("R %s | %d x %d | k=%d | %s | iter.max=%d\n",
            getRversion(), n, d, k, algorithm, max_iter))

# ----------------------------------------------------------------- fit ----
warnings_seen <- character(0)
before <- if (length(DOMS)) read_rapl(DOMS) else numeric(0)
t0 <- Sys.time()
km <- withCallingHandlers(
  stats::kmeans(X, centers = C0, iter.max = max_iter, algorithm = algorithm),
  warning = function(w) {
    warnings_seen <<- c(warnings_seen, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
after <- if (length(DOMS)) read_rapl(DOMS) else numeric(0)

# --------------------------------------------------------------- energy ----
energy <- c(total = NA_real_, pkg = NA_real_, dram = NA_real_, watts = NA_real_)
if (length(DOMS)) {
  deltas <- vapply(seq_along(DOMS), function(i) {
    ((after[i] - before[i]) %% (DOMS[[i]]$max + 1)) / 1e6
  }, numeric(1))
  kinds <- vapply(DOMS, function(x) x$kind, character(1))
  pkg <- sum(deltas[kinds == "package"])
  dram <- sum(deltas[kinds == "dram"])
  energy <- c(total = pkg + dram, pkg = pkg, dram = dram,
              watts = if (elapsed > 0) (pkg + dram) / elapsed else NA_real_)
}

# Hartigan-Wong has no closed-form distance count and stats::kmeans exposes no
# counter, so this is the same analytic convention used for scikit-learn --
# a full n*k assignment pass per outer iteration plus the initial one. For
# Hartigan-Wong it is an UPPER BOUND: the optimal-transfer stage only scans
# live clusters and the quick-transfer stage touches two centres per point.
n_dist <- (as.numeric(km$iter) + 1) * as.numeric(n) * as.numeric(k)

converged <- !any(grepl("did not converge", warnings_seen))

lines <- c(
  sprintf("implementation=r_%s", tolower(gsub("[^A-Za-z]", "_", algorithm))),
  sprintf("r_algorithm=%s", algorithm),
  sprintf("n_iterations=%d", as.integer(km$iter)),
  sprintf("wall_clock_seconds=%.9f", elapsed),
  sprintf("sse=%.9e", sum(km$tot.withinss)),
  sprintf("n_distance_calculations=%.0f", n_dist),
  sprintf("distance_count_method=analytic_upper_bound"),
  sprintf("converged=%s", if (converged) "True" else "False"),
  sprintf("energy_total_joules=%s", format(energy[["total"]], digits = 12)),
  sprintf("energy_pkg_joules=%s", format(energy[["pkg"]], digits = 12)),
  sprintf("energy_dram_joules=%s", format(energy[["dram"]], digits = 12)),
  sprintf("power_watts=%s", format(energy[["watts"]], digits = 12)),
  sprintf("r_warnings=%s", paste(warnings_seen, collapse = " | "))
)
writeLines(lines, out_path)
cat(sprintf("  %d iters  %.3fs  SSE=%.6e%s\n", km$iter, elapsed,
            sum(km$tot.withinss),
            if (length(warnings_seen)) "  [warned]" else ""))
