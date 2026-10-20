## Source code ##

## Forecast settings (edit here, then re-run the pipeline and re-render the dashboard) ##
forecast_n_reference       <- 3L     # number of prior same-stage cohorts averaged per transition
forecast_regime_aware      <- TRUE   # TRUE: adjust for the transfer-into-FCBpre policy (30 units + GPA gate); FALSE: historical-average model
forecast_regime_post_from  <- "F25"  # first term whose transitions are post-policy (F25>S26 is the first affected step)
forecast_regime_uplift     <- 0      # optional scenario: add to P(FCBpre->FCBud) at stage 3>4 (e.g., 0.03); 0 = none

source(here::here("Rscripts/markov_01_etl.R"))
source(here::here("Rscripts/markov_04_matrixview.R"))
source(here::here("Rscripts/markov_02_decomposition.R"))
source(here::here("Rscripts/markov_03_forecast_S27.R"))
source(here::here("Rscripts/markov_05_waterfall.R"))
source(here::here("Rscripts/markov_06_backcast_v2.R"))

# Optional dashboard render at end of pipeline.
# The dashboard can be rendered manually in RStudio when needed.
render_dashboard <- FALSE
dashboard_qmd <- "enrollment_dashboard.qmd"

find_r_exe <- function() {
	for (candidate in c("Rscript", "R")) {
		cmd <- Sys.which(candidate)
		if (nzchar(cmd[[1]])) {
			return(normalizePath(cmd[[1]], winslash = "/", mustWork = FALSE))
		}
	}

	windows_candidates <- c(
		"C:/Program Files/R/R-4.5.1/bin/x64/Rscript.exe",
		"C:/Program Files/R/R-4.5.1/bin/x64/R.exe",
		"C:/Program Files/R/R-4.5.2/bin/x64/Rscript.exe",
		"C:/Program Files/R/R-4.5.2/bin/x64/R.exe",
		"C:/Program Files/R/R-4.5.3/bin/x64/Rscript.exe",
		"C:/Program Files/R/R-4.5.3/bin/x64/R.exe",
		"C:/Program Files/R/R-4.6.0/bin/x64/Rscript.exe",
		"C:/Program Files/R/R-4.6.0/bin/x64/R.exe",
		"C:/Program Files/R/R-4.6.1/bin/x64/Rscript.exe",
		"C:/Program Files/R/R-4.6.1/bin/x64/R.exe"
	)
	found <- windows_candidates[file.exists(windows_candidates)]
	if (length(found)) {
		return(normalizePath(found[[1]], winslash = "/", mustWork = FALSE))
	}

	""
}

if (isTRUE(render_dashboard)) {
	dashboard_path <- here::here(dashboard_qmd)
	dashboard_html <- here::here("enrollment_dashboard.html")
	prior_html_exists <- file.exists(dashboard_html)
	prior_html_mtime <- if (prior_html_exists) file.info(dashboard_html)$mtime else NA_real_

	resolve_quarto_bin <- function() {
		candidate_cmd <- Sys.which("quarto")
		if (nzchar(candidate_cmd) && file.exists(candidate_cmd)) {
			return(normalizePath(candidate_cmd, winslash = "/", mustWork = FALSE))
		}

		env_qpath <- Sys.getenv("QUARTO_PATH", unset = "")
		if (nzchar(env_qpath)) {
			env_candidates <- c(
				file.path(env_qpath, "quarto.exe"),
				file.path(env_qpath, "quarto.cmd"),
				file.path(env_qpath, "quarto"),
				env_qpath
			)
			env_hit <- env_candidates[file.exists(env_candidates)]
			if (length(env_hit)) {
				return(normalizePath(env_hit[[1]], winslash = "/", mustWork = FALSE))
			}
		}

		windows_candidates <- c(
			"C:/Program Files/Quarto/bin/quarto.exe",
			"C:/Program Files/Quarto/bin/quarto.cmd",
			"C:/PROGRA~1/Quarto/bin/quarto.exe",
			"C:/PROGRA~1/Quarto/bin/quarto.cmd",
			"C:/Program Files/RStudio/resources/app/bin/quarto/bin/quarto.exe",
			"C:/Program Files/RStudio/resources/app/bin/quarto/bin/quarto.cmd",
			"C:/PROGRA~1/RStudio/RESOUR~1/app/bin/quarto/bin/quarto.exe",
			"C:/PROGRA~1/RStudio/RESOUR~1/app/bin/quarto/bin/quarto.cmd"
		)
		hit <- windows_candidates[file.exists(windows_candidates)]
		if (length(hit)) {
			return(normalizePath(hit[[1]], winslash = "/", mustWork = FALSE))
		}

		""
	}

	if (!file.exists(dashboard_path)) {
		warning(sprintf("Dashboard source not found: %s", dashboard_path))
	} else {
		rendered <- FALSE
		quarto_bin <- resolve_quarto_bin()

		if (nzchar(quarto_bin)) {
			Sys.setenv(QUARTO_PATH = dirname(quarto_bin))
			Sys.setenv(QUARTO_R = find_r_exe())
			cli_status <- tryCatch(
				system2(
					quarto_bin,
					args = c("render", dashboard_path, "--log-level", "info"),
					error = function(e) e
				),
				error = function(e) e
			)

			if (is.numeric(cli_status) && identical(cli_status, 0L)) {
				message(sprintf("Rendered dashboard via Quarto CLI: %s", dashboard_qmd))
				rendered <- TRUE
			} else {
				message("Quarto CLI render failed in this runtime. Trying R package fallback.")
			}
		}

		if (!rendered && requireNamespace("quarto", quietly = TRUE)) {
			render_result <- tryCatch(
				{
					quarto::quarto_render(
						dashboard_path,
						output_file = basename(dashboard_html),
						execute_dir = dirname(dashboard_path),
						quiet = TRUE
					)
					TRUE
				},
				error = function(e) {
					warning(sprintf("R quarto fallback failed: %s", conditionMessage(e)))
					FALSE
				}
			)

			if (isTRUE(render_result)) {
				message(sprintf("Rendered dashboard via R quarto package: %s", dashboard_qmd))
				rendered <- TRUE
			}
		}

		if (!rendered) {
			message("Skipped dashboard render: no usable Quarto executable was found for this runtime.")
		} else if (!file.exists(dashboard_html)) {
			warning(sprintf("Dashboard render completed but %s was not written.", dashboard_html))
		} else {
			post_html_mtime <- file.info(dashboard_html)$mtime
			if (is.na(prior_html_mtime) || is.na(post_html_mtime) || !identical(post_html_mtime, prior_html_mtime)) {
				message(sprintf("Verified dashboard updated: %s", dashboard_html))
			} else {
				message(sprintf("Dashboard render completed; %s already exists and was not changed.", dashboard_html))
			}
		}
	}
}
