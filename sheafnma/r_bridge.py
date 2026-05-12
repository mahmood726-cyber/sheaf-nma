"""Bridge to R `netmeta` for dataset loading and node-splitting."""
from __future__ import annotations

import json
import os
import shutil
import subprocess
import tempfile
from pathlib import Path


def rscript_path() -> str | None:
    """Locate Rscript: env var, PATH, then known Windows install."""
    if env := os.environ.get("RSCRIPT"):
        if Path(env).exists():
            return env
    found = shutil.which("Rscript")
    if found:
        return found
    win = r"C:\Program Files\R\R-4.5.2\bin\Rscript.exe"
    if Path(win).exists():
        return win
    return None


def _run_rscript(script: str, timeout: int = 120) -> str:
    rs = rscript_path()
    if rs is None:
        raise RuntimeError("Rscript not on PATH; install R 4.5+ and netmeta package")
    with tempfile.NamedTemporaryFile("w", suffix=".R", delete=False, encoding="utf-8") as f:
        f.write(script)
        path = f.name
    try:
        proc = subprocess.run(
            [rs, "--vanilla", path],
            capture_output=True, text=True, timeout=timeout, check=False,
        )
    finally:
        Path(path).unlink(missing_ok=True)
    if proc.returncode != 0:
        raise RuntimeError(
            f"Rscript failed (exit {proc.returncode})\n"
            f"stdout:\n{proc.stdout}\nstderr:\n{proc.stderr}"
        )
    return proc.stdout


def list_netmeta_datasets() -> list[str]:
    """Names of every dataset in the netmeta R package."""
    script = r"""
suppressWarnings(suppressMessages(library(netmeta)))
suppressWarnings(suppressMessages(library(jsonlite)))
d <- data(package = "netmeta")$results
items <- as.character(d[, "Item"])
cat(toJSON(items, auto_unbox = FALSE))
"""
    out = _run_rscript(script)
    return json.loads(out.strip())


def load_netmeta_dataset(name: str) -> dict:
    """Load a netmeta dataset and return it in sheafnma standard format.

    netmeta datasets vary in format (long contrast vs. arm-based); we use
    pairwise() to normalise to a contrast-level data frame.
    """
    script = f"""
suppressWarnings(suppressMessages({{
  library(netmeta)
  library(jsonlite)
}}))
data("{name}")
df <- get("{name}")
# Try contrast-level first
contrast_cols <- intersect(c("TE", "seTE", "treat1", "treat2", "studlab"), names(df))
if (length(contrast_cols) >= 4) {{
  out <- data.frame(
    study  = as.character(df$studlab),
    treat1 = as.character(df$treat1),
    treat2 = as.character(df$treat2),
    effect = as.numeric(df$TE),
    se     = as.numeric(df$seTE),
    stringsAsFactors = FALSE
  )
}} else {{
  # arm-based: convert via pairwise()
  arm_cols <- c("treatment", "studlab", "event", "n")
  if (all(arm_cols %in% names(df))) {{
    pw <- pairwise(treat = treatment, event = event, n = n,
                   studlab = studlab, data = df, sm = "OR")
    out <- data.frame(
      study  = as.character(pw$studlab),
      treat1 = as.character(pw$treat1),
      treat2 = as.character(pw$treat2),
      effect = as.numeric(pw$TE),
      se     = as.numeric(pw$seTE),
      stringsAsFactors = FALSE
    )
  }} else {{
    stop("dataset format not recognised: ", "{name}")
  }}
}}
out <- out[is.finite(out$effect) & is.finite(out$se) & out$se > 0, ]
cat(toJSON(out, dataframe = "rows", auto_unbox = FALSE))
"""
    out = _run_rscript(script)
    rows = json.loads(out.strip())
    if not rows:
        raise ValueError(f"netmeta dataset {name!r} produced no usable contrasts")
    nodes = sorted({r["treat1"] for r in rows} | {r["treat2"] for r in rows})
    edges = [
        {"study": r["study"], "treat1": r["treat1"], "treat2": r["treat2"],
         "effect": float(r["effect"]), "se": float(r["se"])}
        for r in rows
    ]
    return {"nodes": nodes, "edges": edges, "source": f"netmeta::{name}"}


def run_netsplit(network: dict) -> dict:
    """Shell to netmeta::netsplit for per-edge node-splitting p-values."""
    rows = ",\n".join(
        f'  list(studlab="{e["study"]}", treat1="{e["treat1"]}", '
        f'treat2="{e["treat2"]}", TE={e["effect"]}, seTE={e["se"]})'
        for e in network["edges"]
    )
    script = f"""
suppressWarnings(suppressMessages({{
  library(netmeta); library(jsonlite)
}}))
contrasts <- list(
{rows}
)
df <- do.call(rbind, lapply(contrasts, as.data.frame))
nm <- netmeta(TE = df$TE, seTE = df$seTE,
              treat1 = df$treat1, treat2 = df$treat2,
              studlab = df$studlab, sm = "MD", common = FALSE, random = TRUE)
ns <- netsplit(nm)
# ns$random has columns: comparison, TE, seTE, lower, upper, statistic, p
# 'comparison' is "treat1<sep>treat2"; sep is nm$sep.trts (default ":")
out_df <- ns$random
if (is.null(out_df) || nrow(out_df) == 0) out_df <- ns$compare.random
sep <- nm$sep.trts
if (is.null(sep) || nchar(sep) == 0) sep <- ":"
# Split comparison into treat1 / treat2
split_parts <- strsplit(as.character(out_df$comparison), sep, fixed = TRUE)
out_df$treat1 <- sapply(split_parts, function(x) x[1])
out_df$treat2 <- sapply(split_parts, function(x) if (length(x) >= 2) x[2] else NA_character_)
# Find p column (named "p" in netmeta 3.x)
pcol <- grep("^p$", names(out_df), value = TRUE)
if (length(pcol) == 0) pcol <- grep("^p", names(out_df), value = TRUE)[1]
out_df$p_value <- as.numeric(out_df[[pcol[1]]])
out_df <- out_df[, c("treat1", "treat2", "p_value")]
out_df <- out_df[is.finite(out_df$p_value) & !is.na(out_df$treat1) & !is.na(out_df$treat2), ]
cat(toJSON(out_df, dataframe = "rows", auto_unbox = FALSE))
"""
    out = _run_rscript(script, timeout=300)
    rows_out = json.loads(out.strip()) if out.strip() else []
    return {"edges": [
        {"treat1": r["treat1"], "treat2": r["treat2"],
         "p_value": float(r["p_value"])}
        for r in rows_out
    ]}
