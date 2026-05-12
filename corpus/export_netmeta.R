# corpus/export_netmeta.R -- exports every netmeta built-in NMA dataset
# to the sheafnma standard CSV schema (study, treat1, treat2, effect, se).
#
# Run from the SheafNMA repo root:
#   Rscript corpus/export_netmeta.R
#
# Output: corpus/data/<dataset_name>.csv (one per dataset).
# Records: corpus/data/_manifest.json (dataset, n_treatments, n_edges,
# n_studies, netmeta_version, R_version, exported_at).

suppressWarnings(suppressMessages({
  library(netmeta)
  library(jsonlite)
}))

out_dir <- file.path("corpus", "data")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

datasets <- as.character(data(package = "netmeta")$results[, "Item"])
manifest <- list()

# Helper: finalise a contrast data.frame and write if valid
finalise <- function(out, name, out_dir) {
  out <- out[is.finite(out$effect) & is.finite(out$se) & out$se > 0, ]
  if (nrow(out) < 3) {
    cat("  SKIP -- <3 usable rows\n")
    return(NULL)
  }
  fn <- file.path(out_dir, paste0(name, ".csv"))
  write.csv(out, fn, row.names = FALSE, fileEncoding = "UTF-8")
  n_treat <- length(unique(c(out$treat1, out$treat2)))
  list(
    dataset      = name,
    file         = paste0(name, ".csv"),
    n_treatments = n_treat,
    n_edges      = nrow(out),
    n_studies    = length(unique(out$study))
  )
}

for (name in datasets) {
  cat("Exporting:", name, "...\n")
  res <- tryCatch({
    e <- new.env()
    data(list = name, envir = e)
    df <- get(name, envir = e)
    nms <- names(df)

    # ------------------------------------------------------------------
    # FORMAT A: contrast-level with TE / seTE
    # Datasets: Senn2013, Linde2016
    # ------------------------------------------------------------------
    if ("TE" %in% nms && "seTE" %in% nms &&
        "treat1" %in% nms && "treat2" %in% nms) {
      studlab <- if ("studlab" %in% nms) df$studlab else seq_len(nrow(df))
      out <- data.frame(
        study  = as.character(studlab),
        treat1 = as.character(df$treat1),
        treat2 = as.character(df$treat2),
        effect = as.numeric(df$TE),
        se     = as.numeric(df$seTE),
        stringsAsFactors = FALSE
      )
      finalise(out, name, out_dir)

    # ------------------------------------------------------------------
    # FORMAT B: contrast-level with lnOR / selnOR
    # Datasets: Linde2016 fallback
    # ------------------------------------------------------------------
    } else if ("lnOR" %in% nms && "selnOR" %in% nms &&
               "treat1" %in% nms && "treat2" %in% nms) {
      studlab <- if ("id" %in% nms) df$id else seq_len(nrow(df))
      out <- data.frame(
        study  = as.character(studlab),
        treat1 = as.character(df$treat1),
        treat2 = as.character(df$treat2),
        effect = as.numeric(df$lnOR),
        se     = as.numeric(df$selnOR),
        stringsAsFactors = FALSE
      )
      finalise(out, name, out_dir)

    # ------------------------------------------------------------------
    # FORMAT C: arm-based long format (treatment, event/r, total/N)
    # Datasets: Baker2009 (treatment/exac/total),
    #           Dogliotti2014 (treatment/stroke/total),
    #           Dong2013 (treatment/death/randomized),
    #           Gurusamy2011 (treatment/death/n),
    #           Woods2010 (treatment/r/N)
    # ------------------------------------------------------------------
    } else if (any(c("treatment", "Treatment") %in% nms)) {
      treat_col <- if ("treatment" %in% nms) "treatment" else "Treatment"
      study_col <- if ("study" %in% nms) "study" else if ("id" %in% nms) "id" else if ("author" %in% nms) "author" else NULL
      if (is.null(study_col)) {
        cat("  SKIP -- no study/id/author column found\n"); next
      }
      # Detect event column
      event_col <- intersect(c("exac", "stroke", "death", "r", "event"), nms)[1]
      n_col     <- intersect(c("total", "n", "randomized", "N"), nms)[1]
      if (is.na(event_col) || is.na(n_col)) {
        cat("  SKIP -- cannot identify event/n columns\n"); next
      }
      tmp <- data.frame(
        treatment = df[[treat_col]],
        event     = as.numeric(df[[event_col]]),
        n         = as.numeric(df[[n_col]]),
        studlab   = df[[study_col]],
        stringsAsFactors = FALSE
      )
      pw <- pairwise(treat = treatment, event = event, n = n,
                     studlab = studlab, data = tmp, sm = "OR")
      out <- data.frame(
        study  = as.character(pw$studlab),
        treat1 = as.character(pw$treat1),
        treat2 = as.character(pw$treat2),
        effect = as.numeric(pw$TE),
        se     = as.numeric(pw$seTE),
        stringsAsFactors = FALSE
      )
      finalise(out, name, out_dir)

    # ------------------------------------------------------------------
    # FORMAT D: wide multi-arm (Treatment1/2/3, y1/sd1/n1, ...)
    # Datasets: Franchini2012, parkinson, Stowe2010
    # ------------------------------------------------------------------
    } else if (any(c("Treatment1", "t1") %in% nms) &&
               any(c("y1", "y.1") %in% nms)) {
      # Determine treat columns prefix
      if ("Treatment1" %in% nms) {
        tcols <- c("Treatment1", "Treatment2", "Treatment3")
      } else {
        tcols <- c("t1", "t2", "t3")
      }
      study_col <- if ("Study" %in% nms) "Study" else if ("study" %in% nms) "study" else NULL
      if (is.null(study_col)) {
        cat("  SKIP -- no Study/study column\n"); next
      }
      arms <- list()
      for (i in 1:3) {
        tc <- if ("Treatment1" %in% nms) paste0("Treatment", i) else paste0("t", i)
        yc <- paste0("y", i); sdc <- paste0("sd", i); nc <- paste0("n", i)
        if (all(c(tc, yc, sdc, nc) %in% nms)) {
          arms[[i]] <- data.frame(
            treatment = as.character(df[[tc]]),
            mean      = as.numeric(df[[yc]]),
            sd        = as.numeric(df[[sdc]]),
            n         = as.numeric(df[[nc]]),
            studlab   = as.character(df[[study_col]]),
            stringsAsFactors = FALSE
          )
        }
      }
      long_df <- do.call(rbind, arms)
      long_df <- long_df[!is.na(long_df$treatment) & nchar(trimws(long_df$treatment)) > 0, ]
      pw <- pairwise(treat = treatment, mean = mean, sd = sd, n = n,
                     studlab = studlab, data = long_df, sm = "MD")
      out <- data.frame(
        study  = as.character(pw$studlab),
        treat1 = as.character(pw$treat1),
        treat2 = as.character(pw$treat2),
        effect = as.numeric(pw$TE),
        se     = as.numeric(pw$seTE),
        stringsAsFactors = FALSE
      )
      finalise(out, name, out_dir)

    # ------------------------------------------------------------------
    # FORMAT E: smokingcessation -- wide binary (event1/n1/event2/n2/...)
    # ------------------------------------------------------------------
    } else if (all(c("event1", "n1", "event2", "n2", "treat1", "treat2") %in% nms)) {
      rows <- list()
      study_ids <- seq_len(nrow(df))
      for (i in 1:3) {
        ec <- paste0("event", i); nc <- paste0("n", i); tc <- paste0("treat", i)
        if (all(c(ec, nc, tc) %in% nms)) {
          rows[[i]] <- data.frame(
            treatment = as.character(df[[tc]]),
            event     = as.numeric(df[[ec]]),
            n         = as.numeric(df[[nc]]),
            studlab   = as.character(study_ids),
            stringsAsFactors = FALSE
          )
        }
      }
      long_df <- do.call(rbind, rows)
      long_df <- long_df[!is.na(long_df$treatment) & long_df$n > 0, ]
      pw <- pairwise(treat = treatment, event = event, n = n,
                     studlab = studlab, data = long_df, sm = "OR")
      out <- data.frame(
        study  = as.character(pw$studlab),
        treat1 = as.character(pw$treat1),
        treat2 = as.character(pw$treat2),
        effect = as.numeric(pw$TE),
        se     = as.numeric(pw$seTE),
        stringsAsFactors = FALSE
      )
      finalise(out, name, out_dir)

    # ------------------------------------------------------------------
    # FORMAT F: Linde2015 -- wide multi-arm with treatment1/2/3 + n/resp
    # ------------------------------------------------------------------
    } else if (all(c("treatment1", "treatment2", "n1", "resp1", "n2", "resp2") %in% nms)) {
      study_col <- if ("id" %in% nms) "id" else if ("author" %in% nms) "author" else NULL
      if (is.null(study_col)) {
        cat("  SKIP -- no study identifier\n"); next
      }
      rows <- list()
      for (i in 1:3) {
        tc <- paste0("treatment", i); nc <- paste0("n", i); ec <- paste0("resp", i)
        if (all(c(tc, nc, ec) %in% nms)) {
          rows[[i]] <- data.frame(
            treatment = as.character(df[[tc]]),
            event     = as.numeric(df[[ec]]),
            n         = as.numeric(df[[nc]]),
            studlab   = as.character(df[[study_col]]),
            stringsAsFactors = FALSE
          )
        }
      }
      long_df <- do.call(rbind, rows)
      long_df <- long_df[!is.na(long_df$treatment) & nchar(trimws(long_df$treatment)) > 0 &
                           is.finite(long_df$event) & is.finite(long_df$n) & long_df$n > 0, ]
      pw <- pairwise(treat = treatment, event = event, n = n,
                     studlab = studlab, data = long_df, sm = "OR")
      out <- data.frame(
        study  = as.character(pw$studlab),
        treat1 = as.character(pw$treat1),
        treat2 = as.character(pw$treat2),
        effect = as.numeric(pw$TE),
        se     = as.numeric(pw$seTE),
        stringsAsFactors = FALSE
      )
      finalise(out, name, out_dir)

    # ------------------------------------------------------------------
    # FORMAT G: dietaryfat -- wide binary with d1/d2/d3 events + years
    # (treat1/treat2/treat3 already named; d=deaths, years=person-time)
    # Use Poisson / rate approach or fallback to arm-based pairwise
    # ------------------------------------------------------------------
    } else if (all(c("treat1", "treat2", "d1", "d2") %in% nms) && "years1" %in% nms) {
      # dietaryfat has no n (patient count) but has years (person-years).
      # We use incidence rates: treat as approximate, derive logRR from rates.
      # pairwise() supports rate data: pairwise(treat=, event=, time=, ...)
      rows <- list()
      study_ids <- if ("ID" %in% nms) df$ID else seq_len(nrow(df))
      for (i in 1:3) {
        tc <- paste0("treat", i); dc <- paste0("d", i); yc <- paste0("years", i)
        if (all(c(tc, dc, yc) %in% nms)) {
          rows[[i]] <- data.frame(
            treatment = as.character(df[[tc]]),
            event     = as.numeric(df[[dc]]),
            time      = as.numeric(df[[yc]]),
            studlab   = as.character(study_ids),
            stringsAsFactors = FALSE
          )
        }
      }
      long_df <- do.call(rbind, rows)
      long_df <- long_df[!is.na(long_df$treatment) & nchar(trimws(long_df$treatment)) > 0 &
                           is.finite(long_df$event) & is.finite(long_df$time) &
                           long_df$time > 0, ]
      pw <- pairwise(treat = treatment, event = event, time = time,
                     studlab = studlab, data = long_df, sm = "IRR")
      out <- data.frame(
        study  = as.character(pw$studlab),
        treat1 = as.character(pw$treat1),
        treat2 = as.character(pw$treat2),
        effect = as.numeric(pw$TE),
        se     = as.numeric(pw$seTE),
        stringsAsFactors = FALSE
      )
      finalise(out, name, out_dir)

    } else {
      cat("  SKIP -- unrecognised format\n  (cols:", paste(nms, collapse=", "), ")\n")
      next
    }
  }, error = function(err) {
    cat("  ERROR:", conditionMessage(err), "\n")
    NULL
  })
  if (!is.null(res)) manifest[[name]] <- res
}

manifest_path <- file.path(out_dir, "_manifest.json")
write_json(list(
  datasets        = manifest,
  netmeta_version = as.character(packageVersion("netmeta")),
  r_version       = R.version.string,
  exported_at     = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
), manifest_path, auto_unbox = TRUE, pretty = TRUE)

cat("Done. Manifest at:", manifest_path, "\n")
cat("Datasets exported:", length(manifest), "/", length(datasets), "\n")
