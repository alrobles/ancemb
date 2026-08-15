libs <- file.path(R_PACKAGE_DIR, "libs", R_ARCH)
dir.create(libs, recursive = TRUE, showWarnings = FALSE)
for (file in c("symbols.rds", Sys.glob(paste0("*", SHLIB_EXT)))) {
  if (file.exists(file)) {
    file.copy(file, file.path(libs, file))
  }
}

inst_stan <- file.path("..", "inst", "stan")
if (dir.exists(inst_stan)) {
  warning(
    "Stan models in inst/stan/ are deprecated in {instantiate} ",
    ">= 0.0.4.9001 (2024-01-03). Please put them in src/stan/ instead."
  )
  if (file.exists("stan")) {
    warning("src/stan/ already exists. Not copying models from inst/stan/.")
  } else {
    message("Copying inst/stan/ to src/stan/.")
    fs::dir_copy(path = inst_stan, new_path = "stan")
  }
}

bin <- file.path(R_PACKAGE_DIR, "bin")
if (!file.exists(bin)) {
  dir.create(bin, recursive = TRUE, showWarnings = FALSE)
}
bin_stan <- file.path(bin, "stan")
fs::dir_copy(path = "stan", new_path = bin_stan)

is_r_cmd_check <- function() {
  r_package_dir <- if (exists("R_PACKAGE_DIR", inherits = TRUE)) {
    get("R_PACKAGE_DIR", inherits = TRUE)
  } else {
    Sys.getenv("R_PACKAGE_DIR", "")
  }
  in_check_dir <- grepl("\\.Rcheck", r_package_dir, fixed = FALSE)
  !is.na(Sys.getenv("_R_CHECK_PACKAGE_NAME_", NA)) ||
    tolower(Sys.getenv("_R_CHECK_LICENSE_")) %in% c("true", "1") ||
    "CheckExEnv" %in% search() ||
    in_check_dir
}

if (instantiate::stan_cmdstan_exists() && !is_r_cmd_check()) {
  callr::r(
    func = function(bin_stan) {
      instantiate::stan_package_compile(
        models = instantiate::stan_package_model_files(path = bin_stan)
      )
    },
    args = list(bin_stan = bin_stan),
    show = TRUE,
    stderr = "2>&1"
  )
} else if (is_r_cmd_check()) {
  message(
    "Skipping Stan pre-compilation under R CMD check; ",
    "models will compile on first use."
  )
} else {
  warning(
    "CmdStan not found. The Stan models will not be pre-compiled during ",
    "installation. Install CmdStan and reinstall 'ancemb' to use the models."
  )
}
