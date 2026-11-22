# Run from the repository root. Set RAMAN_DATA_DIR if the dataset is elsewhere.
data_dir <- Sys.getenv("RAMAN_DATA_DIR", unset = file.path("data", "261023"))
out_dir <- Sys.getenv("RAMAN_OUTPUT_DIR", unset = file.path("results", "analysis_out"))

if (!dir.exists(data_dir)) {
  stop("Update data_dir to the folder containing +SPORE and soil blank.")
}

packages <- c("dplyr", "ggplot2", "mclust")
missing_packages <- packages[
  !vapply(packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop(
    "Install required packages: ",
    paste(missing_packages, collapse = ", ")
  )
}

library(dplyr)
library(ggplot2)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Load and match spectra
get_group <- function(files) {
  paths <- gsub(
    "\\\\",
    "/",
    normalizePath(files, winslash = "/", mustWork = FALSE)
  )
  ifelse(
    grepl("/\\+SPORE/", paths),
    "+SPORE",
    ifelse(grepl("/soil blank/", paths), "soil blank", NA_character_)
  )
}

all_files <- list.files(
  data_dir,
  pattern = "\\.txt$",
  recursive = TRUE,
  full.names = TRUE
)
all_paths <- gsub(
  "\\\\",
  "/",
  normalizePath(all_files, winslash = "/", mustWork = FALSE)
)
file_kind <- ifelse(
  grepl("/background/", all_paths),
  "background",
  ifelse(
    grepl("/processed/", all_paths),
    "processed",
    ifelse(grepl("^original-", basename(all_files)), "original", "other")
  )
)
cat("File inventory\n")
print(addmargins(table(get_group(all_files), file_kind, useNA = "ifany")))

raw_files <- list.files(
  data_dir,
  pattern = "^original-.*\\.txt$",
  recursive = TRUE,
  full.names = TRUE
)
processed_files <- list.files(
  data_dir,
  pattern = "^processed-.*\\.txt$",
  recursive = TRUE,
  full.names = TRUE
)

if (!length(raw_files) || !length(processed_files)) {
  stop("Original and processed spectrum files were not both found.")
}

read_spectrum <- function(file) {
  x <- read.table(
    file,
    header = FALSE,
    sep = "",
    colClasses = "numeric",
    comment.char = "",
    quote = ""
  )
  if (ncol(x) != 2L) stop("Expected two columns in: ", file)
  names(x) <- c("raman_shift", "intensity")
  x
}

spectra <- lapply(raw_files, read_spectrum)
processed_spectra_unordered <- lapply(processed_files, read_spectrum)
reference_axis <- spectra[[1]]$raman_shift

raw_keys <- paste(
  get_group(raw_files),
  sub("^original-", "", basename(raw_files)),
  sep = "|"
)
processed_keys <- paste(
  get_group(processed_files),
  sub("^processed-", "", basename(processed_files)),
  sep = "|"
)

if (anyDuplicated(raw_keys) || anyDuplicated(processed_keys)) {
  stop("Duplicate spectrum filenames were found.")
}

processing_order <- match(raw_keys, processed_keys)
if (anyNA(processing_order)) {
  stop("Some original spectra do not have a matching processed file.")
}
processed_spectra <- processed_spectra_unordered[processing_order]

same_axis <- function(x, axis) {
  length(x$raman_shift) == length(axis) &&
    isTRUE(all.equal(x$raman_shift, axis, tolerance = 1e-8))
}
raw_axis_ok <- vapply(spectra, same_axis, logical(1), axis = reference_axis)
processed_axis_ok <- vapply(
  processed_spectra,
  same_axis,
  logical(1),
  axis = reference_axis
)

if (!all(raw_axis_ok) || !all(processed_axis_ok)) {
  stop("Some spectra do not share the reference Raman-shift axis.")
}

raw_matrix <- do.call(rbind, lapply(spectra, function(x) x[["intensity"]]))
processed_matrix <- do.call(
  rbind,
  lapply(processed_spectra, function(x) x[["intensity"]])
)

sample_id <- sub(
  "\\.txt$",
  "",
  sub("^original-", "", basename(raw_files))
)

qc <- data.frame(
  group = get_group(raw_files),
  sample_id = sample_id,
  date = substr(sample_id, 1, 8),
  hour = substr(sample_id, 10, 11),
  raw_max = apply(raw_matrix, 1, max),
  raw_median = apply(raw_matrix, 1, median),
  processed_max = apply(processed_matrix, 1, max),
  processed_q99 = apply(
    processed_matrix, 1, quantile, probs = 0.99
  ),
  stringsAsFactors = FALSE
) %>%
  mutate(block = paste(group, date, hour, sep = " | "))

cat("\nSpectra matched:", nrow(qc), "\n")
cat("Raw spectra on shared axis:", sum(raw_axis_ok), "of", length(raw_axis_ok), "\n")
cat(
  "Processed spectra on shared axis:",
  sum(processed_axis_ok), "of", length(processed_axis_ok), "\n"
)
cat(
  "Raman-shift range:",
  min(reference_axis), "to", max(reference_axis), "cm^-1\n"
)

qc_summary <- qc %>%
  group_by(group, date, hour) %>%
  summarise(
    n_spectra = n(),
    median_raw_max = median(raw_max),
    median_raw_counts = median(raw_median),
    median_processed_max = median(processed_max),
    median_processed_q99 = median(processed_q99),
    .groups = "drop"
  )
cat("\nIntensity QC by acquisition block\n")
print(qc_summary, n = Inf, digits = 3)

raw_max_index <- max.col(raw_matrix, ties.method = "first")
processed_max_index <- max.col(processed_matrix, ties.method = "first")
extreme_check <- qc %>%
  mutate(
    raw_max_shift = reference_axis[raw_max_index],
    processed_max_shift = reference_axis[processed_max_index]
  ) %>%
  arrange(desc(processed_max)) %>%
  slice_head(n = 5) %>%
  select(
    group, sample_id, raw_max, raw_max_shift,
    processed_max, processed_max_shift, processed_q99
  )
cat("\nFive spectra with the highest processed maximum\n")
print(extreme_check, n = Inf, digits = 3)

# Background and band scores
background_window <- reference_axis >= 1800 & reference_axis <= 2600
background_qc <- data.frame(
  group = qc$group,
  sample_id = qc$sample_id,
  block = qc$block,
  original_background = apply(
    raw_matrix[, background_window, drop = FALSE],
    1,
    median
  ),
  processed_background = apply(
    processed_matrix[, background_window, drop = FALSE],
    1,
    median
  )
)

band_defs <- list(
  b659 = list(center = 659, left = c(610, 635), right = c(690, 715)),
  b824 = list(center = 824, left = c(780, 800), right = c(850, 870)),
  b1017 = list(center = 1017, left = c(955, 985), right = c(1055, 1085)),
  b1395 = list(center = 1395, left = c(1330, 1360), right = c(1420, 1450)),
  b1575 = list(center = 1575, left = c(1520, 1545), right = c(1610, 1640))
)

score_band <- function(M, center, left, right) {
  peak_idx <- abs(reference_axis - center) <= 10
  left_idx <- reference_axis >= left[1] & reference_axis <= left[2]
  right_idx <- reference_axis >= right[1] & reference_axis <= right[2]
  if (!any(peak_idx) || !any(left_idx) || !any(right_idx)) {
    stop("A peak or flank window has no Raman-shift points.")
  }
  peak_height <- apply(M[, peak_idx, drop = FALSE], 1, max)
  flank_baseline <- (
    rowMeans(M[, left_idx, drop = FALSE]) +
      rowMeans(M[, right_idx, drop = FALSE])
  ) / 2
  peak_height - flank_baseline
}

processed_band_scores <- as.data.frame(
  lapply(
    band_defs,
    function(b) score_band(
      processed_matrix, b$center, b$left, b$right
    )
  )
)
raw_band_scores <- as.data.frame(
  lapply(
    band_defs,
    function(b) score_band(raw_matrix, b$center, b$left, b$right)
  )
)

band_features <- bind_cols(
  qc %>% select(group, sample_id, date, hour, block),
  processed_band_scores,
  background_qc %>% select(original_background, processed_background)
)

block_summary <- band_features %>%
  group_by(group, date, hour, block) %>%
  summarise(
    n_spectra = n(),
    across(
      all_of(names(band_defs)),
      median,
      .names = "median_{.col}"
    ),
    median_original_background = median(original_background),
    median_processed_background = median(processed_background),
    .groups = "drop"
  )
cat("\nMedian band scores and background by acquisition block\n")
print(block_summary, n = Inf, digits = 3)

safe_spearman <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 3 || sd(x[keep]) == 0 || sd(y[keep]) == 0) {
    return(NA_real_)
  }
  cor(x[keep], y[keep], method = "spearman")
}

blocks <- unique(band_features$block)
band_names <- names(band_defs)
other_bands <- setdiff(band_names, "b1017")

feature_correlations <- bind_rows(lapply(blocks, function(current_block) {
  d <- band_features[band_features$block == current_block, ]
  data.frame(
    block = current_block,
    n_spectra = nrow(d),
    b659 = safe_spearman(d$b1017, d$b659),
    b824 = safe_spearman(d$b1017, d$b824),
    b1395 = safe_spearman(d$b1017, d$b1395),
    b1575 = safe_spearman(d$b1017, d$b1575)
  )
}))
cat("\nWithin-block correlations with the 1017 score\n")
print(feature_correlations, row.names = FALSE, digits = 3)

processing_correlations <- bind_rows(lapply(blocks, function(current_block) {
  use <- band_features$block == current_block
  values <- vapply(
    band_names,
    function(name) {
      safe_spearman(raw_band_scores[[name]][use],
                    processed_band_scores[[name]][use])
    },
    numeric(1)
  )
  data.frame(
    block = current_block,
    n_spectra = sum(use),
    as.list(values),
    check.names = FALSE
  )
}))
cat("\nRaw-to-processed rank correlations by band\n")
print(processing_correlations, row.names = FALSE, digits = 3)

background_correlations <- bind_rows(lapply(blocks, function(current_block) {
  d <- band_features[band_features$block == current_block, ]
  data.frame(
    block = current_block,
    n_spectra = nrow(d),
    rho_original_background =
      safe_spearman(d$b1017, d$original_background),
    rho_processed_background =
      safe_spearman(d$b1017, d$processed_background)
  )
}))
cat("\n1017 score correlations with background proxies\n")
print(background_correlations, row.names = FALSE, digits = 3)

# CH normalization and time checks
ch_score <- score_band(
  processed_matrix,
  center = 2930,
  left = c(2750, 2800),
  right = c(3050, 3100)
)

ratio_df <- data.frame(
  group = band_features$group,
  sample_id = band_features$sample_id,
  date = band_features$date,
  hour = band_features$hour,
  block = band_features$block,
  score_1017 = band_features$b1017,
  ch_score = ch_score,
  ratio_1017_ch = band_features$b1017 / ch_score,
  b659 = band_features$b659,
  b824 = band_features$b824,
  b1395 = band_features$b1395,
  b1575 = band_features$b1575,
  processed_background = band_features$processed_background
)

ch_summary <- bind_rows(lapply(blocks, function(current_block) {
  x <- ratio_df$ch_score[ratio_df$block == current_block]
  data.frame(
    block = current_block,
    n_spectra = length(x),
    min = min(x),
    q05 = as.numeric(quantile(x, 0.05)),
    median = median(x),
    max = max(x),
    n_nonpositive = sum(x <= 0)
  )
}))
cat("\nCH-stretch score by acquisition block\n")
print(ch_summary, row.names = FALSE, digits = 3)

ratio_summary <- bind_rows(lapply(blocks, function(current_block) {
  d <- ratio_df[ratio_df$block == current_block, ]
  data.frame(
    block = current_block,
    n_spectra = nrow(d),
    median_1017 = median(d$score_1017),
    median_CH = median(d$ch_score),
    median_ratio = median(d$ratio_1017_ch),
    ratio_q25 = as.numeric(quantile(d$ratio_1017_ch, 0.25)),
    ratio_q75 = as.numeric(quantile(d$ratio_1017_ch, 0.75))
  )
}))
cat("\n1017/CH score by acquisition block\n")
print(ratio_summary, row.names = FALSE, digits = 3)

cat("\nHighest-ratio spectra in each group\n")
top_ratio_spectra <- bind_rows(
  ratio_df %>%
    filter(group == "soil blank") %>%
    arrange(desc(ratio_1017_ch)) %>%
    slice_head(n = 10),
  ratio_df %>%
    filter(group == "+SPORE") %>%
    arrange(desc(ratio_1017_ch)) %>%
    slice_head(n = 10)
)
print(
  top_ratio_spectra %>%
    select(
      group, sample_id, block, ratio_1017_ch, score_1017,
      ch_score, b1395, b1575, processed_background
    ),
  n = Inf,
  digits = 3
)

time_df <- ratio_df %>%
  mutate(
    timestamp = as.POSIXct(
      sample_id,
      format = "%Y%m%d-%H%M%OS",
      tz = "UTC"
    )
  )
if (anyNA(time_df$timestamp)) {
  stop("Some timestamps could not be parsed from the filenames.")
}

time_correlations <- bind_rows(lapply(split(time_df, time_df$block), function(d) {
  elapsed_min <- as.numeric(
    difftime(d$timestamp, min(d$timestamp), units = "mins")
  )
  data.frame(
    block = d$block[1],
    n_spectra = nrow(d),
    duration_min = round(max(elapsed_min), 1),
    rho_score_1017 = safe_spearman(elapsed_min, d$score_1017),
    rho_ratio = safe_spearman(elapsed_min, d$ratio_1017_ch),
    rho_CH = safe_spearman(elapsed_min, d$ch_score),
    rho_background =
      safe_spearman(elapsed_min, d$processed_background)
  )
}))
cat("\nWithin-block time correlations\n")
print(time_correlations, row.names = FALSE, digits = 3)

# Exploratory mixture models
spore_blocks <- unique(ratio_df$block[ratio_df$group == "+SPORE"])

run_bic <- function(x, block_name, score_type) {
  bic <- mclust::mclustBIC(
    x,
    G = 1:3,
    modelNames = "V",
    verbose = FALSE
  )
  data.frame(
    block = block_name,
    score_type = score_type,
    components = 1:3,
    BIC = as.numeric(bic)
  )
}

bic_results <- bind_rows(lapply(spore_blocks, function(current_block) {
  use <- ratio_df$block == current_block
  bind_rows(
    run_bic(
      band_features$b1017[use],
      current_block,
      "Unnormalized 1017 score"
    ),
    run_bic(
      ratio_df$ratio_1017_ch[use],
      current_block,
      "CH-normalized ratio"
    )
  )
}))
cat("\nExploratory mixture-model BIC values\n")
print(bic_results, row.names = FALSE, digits = 3)

ratio_fits <- setNames(
  lapply(spore_blocks, function(current_block) {
    use <- ratio_df$block == current_block
    mclust::Mclust(
      ratio_df$ratio_1017_ch[use],
      G = 1:3,
      modelNames = "V",
      verbose = FALSE
    )
  }),
  spore_blocks
)

component_signature <- bind_rows(lapply(spore_blocks, function(current_block) {
  use <- ratio_df$block == current_block
  fit <- ratio_fits[[current_block]]
  ordered_means <- order(as.numeric(fit$parameters$mean))
  component <- match(fit$classification, ordered_means)
  max_posterior <- apply(fit$z, 1, max)
  d <- ratio_df[use, ]
  f <- band_features[use, ]
  raw <- raw_band_scores[use, ]
  bg <- background_qc$processed_background[use]

  bind_rows(lapply(seq_along(ordered_means), function(k) {
    in_component <- component == k
    data.frame(
      block = current_block,
      selected_G = fit$G,
      component_low_to_high = k,
      n_spectra = sum(in_component),
      median_ratio = median(d$ratio_1017_ch[in_component]),
      median_CH = median(d$ch_score[in_component]),
      median_processed_background = median(bg[in_component]),
      median_max_posterior = median(max_posterior[in_component]),
      median_1017_processed = median(f$b1017[in_component]),
      median_1395_processed = median(f$b1395[in_component]),
      median_1575_processed = median(f$b1575[in_component]),
      median_1017_raw = median(raw$b1017[in_component]),
      median_1395_raw = median(raw$b1395[in_component]),
      median_1575_raw = median(raw$b1575[in_component])
    )
  }))
}))
cat("\nExploratory ratio-component signatures\n")
print(component_signature, row.names = FALSE, digits = 3)

assignment_overlap <- bind_rows(lapply(spore_blocks, function(current_block) {
  use <- ratio_df$block == current_block
  fit_raw <- mclust::Mclust(
    band_features$b1017[use],
    G = 1:3,
    modelNames = "V",
    verbose = FALSE
  )
  fit_ratio <- ratio_fits[[current_block]]
  raw_component <- match(
    fit_raw$classification,
    order(as.numeric(fit_raw$parameters$mean))
  )
  ratio_component <- match(
    fit_ratio$classification,
    order(as.numeric(fit_ratio$parameters$mean))
  )
  out <- as.data.frame(table(
    unnormalized_component = raw_component,
    CH_component = ratio_component
  ))
  out$block <- current_block
  out$selected_G_unnormalized <- fit_raw$G
  out$selected_G_ratio <- fit_ratio$G
  out
}))
cat("\nOverlap of unnormalized and CH-normalized assignments\n")
print(assignment_overlap, row.names = FALSE)

denominator_sensitivity <- bind_rows(lapply(spore_blocks, function(current_block) {
  use <- ratio_df$block == current_block
  ch <- ratio_df$ch_score[use]
  ratio <- ratio_df$ratio_1017_ch[use]
  scores <- band_features[use, ]

  bind_rows(lapply(c(0, 0.01, 0.05, 0.10), function(q) {
    cutoff <- if (q == 0) -Inf else as.numeric(quantile(ch, q))
    keep <- ch >= cutoff
    fit <- mclust::Mclust(
      ratio[keep],
      G = 1:3,
      modelNames = "V",
      verbose = FALSE
    )
    high_component <- which.max(as.numeric(fit$parameters$mean))
    in_high <- fit$classification == high_component
    kept_scores <- scores[keep, ]
    data.frame(
      block = current_block,
      CH_cutoff = if (q == 0) "None" else paste0("Below ", 100 * q, "% quantile"),
      n_removed = sum(!keep),
      selected_G = fit$G,
      high_component_n = sum(in_high),
      high_component_pct = 100 * mean(in_high),
      median_CH_high = median(ch[keep][in_high]),
      median_1017_high = median(kept_scores$b1017[in_high]),
      median_1395_high = median(kept_scores$b1395[in_high]),
      median_1575_high = median(kept_scores$b1575[in_high])
    )
  }))
}))
cat("\nSensitivity to low CH denominators\n")
print(denominator_sensitivity, row.names = FALSE, digits = 3)

# Figures
label_group <- function(x) {
  ifelse(x == "+SPORE", "Spore-spiked soil", "Soil blank")
}
group_colors <- c(
  "Spore-spiked soil" = "#2166AC",
  "Soil blank" = "#D55E00"
)
plot_meta <- band_features %>%
  mutate(group_label = label_group(group))

fingerprint_idx <- reference_axis >= 600 & reference_axis <= 1800
fingerprint_axis <- reference_axis[fingerprint_idx]
fingerprint_matrix <- processed_matrix[, fingerprint_idx, drop = FALSE]

spectra_long <- data.frame(
  sample_index = rep(seq_len(nrow(fingerprint_matrix)),
                     each = length(fingerprint_axis)),
  raman_shift = rep(fingerprint_axis,
                    times = nrow(fingerprint_matrix)),
  intensity = as.vector(t(fingerprint_matrix))
) %>%
  left_join(
    plot_meta %>%
      mutate(sample_index = row_number()) %>%
      select(sample_index, group_label),
    by = "sample_index"
  )

spectra_summary <- spectra_long %>%
  group_by(group_label, raman_shift) %>%
  summarise(
    median_intensity = median(intensity, na.rm = TRUE),
    lower_50 = quantile(intensity, 0.25, na.rm = TRUE),
    upper_50 = quantile(intensity, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

p_fingerprint <- ggplot(
  spectra_summary,
  aes(
    x = raman_shift,
    y = median_intensity,
    color = group_label,
    fill = group_label
  )
) +
  geom_ribbon(
    aes(ymin = lower_50, ymax = upper_50),
    alpha = 0.10,
    color = NA
  ) +
  geom_line(linewidth = 1.2) +
  geom_vline(
    xintercept = c(659, 824, 1017, 1395, 1575),
    color = "grey45",
    linetype = "dashed",
    linewidth = 0.45
  ) +
  scale_color_manual(values = group_colors) +
  scale_fill_manual(values = group_colors, guide = "none") +
  coord_cartesian(xlim = c(600, 1800)) +
  labs(
    title = "Fingerprint spectra by sample group",
    subtitle = "Median and middle 50%; groups collected on different dates",
    x = expression("Raman shift (cm"^{-1}*")"),
    y = "Processed intensity (counts)",
    color = NULL
  ) +
  theme_bw(base_size = 15) +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 13, hjust = 0.5),
    axis.title = element_text(size = 17, face = "bold"),
    axis.text = element_text(size = 14, color = "black"),
    legend.position = "bottom",
    legend.text = element_text(size = 14),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.4)
  )

print(p_fingerprint)
ggsave(
  file.path(out_dir, "Figure1_fingerprint_spectra.png"),
  p_fingerprint, width = 11, height = 7, dpi = 600
)
ggsave(
  file.path(out_dir, "Figure1_fingerprint_spectra.pdf"),
  p_fingerprint, width = 11, height = 7
)

ratio_plot_df <- ratio_df %>%
  mutate(
    group_label = factor(
      label_group(group),
      levels = c("Spore-spiked soil", "Soil blank")
    )
  )

p_ratio <- ggplot(
  ratio_plot_df,
  aes(x = group_label, y = ratio_1017_ch, color = group_label)
) +
  geom_hline(
    yintercept = 0,
    color = "grey45",
    linetype = "dashed",
    linewidth = 0.5
  ) +
  geom_boxplot(
    aes(fill = group_label),
    width = 0.34,
    outlier.shape = NA,
    alpha = 0.18,
    linewidth = 0.9,
    show.legend = FALSE
  ) +
  geom_jitter(
    width = 0.13,
    size = 1.8,
    alpha = 0.35,
    show.legend = FALSE
  ) +
  scale_color_manual(values = group_colors, guide = "none") +
  scale_fill_manual(values = group_colors, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.06, 0.10))) +
  labs(
    title = "CH-normalized 1017 score",
    x = NULL,
    y = "1017/CH score"
  ) +
  theme_bw(base_size = 15) +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title.y = element_text(size = 17, face = "bold"),
    axis.text = element_text(size = 14, color = "black"),
    axis.text.x = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.4)
  )

print(p_ratio)
ggsave(
  file.path(out_dir, "Figure2_CH_normalized_1017.png"),
  p_ratio, width = 8.5, height = 7, dpi = 600
)
ggsave(
  file.path(out_dir, "Figure2_CH_normalized_1017.pdf"),
  p_ratio, width = 8.5, height = 7
)

band_plot_df <- band_features %>%
  mutate(
    group_label = label_group(group)
  ) %>%
  filter(
    is.finite(b1017),
    is.finite(b1395),
    is.finite(b1575)
  ) %>%
  group_by(group_label) %>%
  mutate(panel = paste0(group_label, " (n = ", n(), ")")) %>%
  ungroup()

panel_levels <- c(
  paste0(
    "Spore-spiked soil (n = ",
    sum(band_plot_df$group_label == "Spore-spiked soil"),
    ")"
  ),
  paste0(
    "Soil blank (n = ",
    sum(band_plot_df$group_label == "Soil blank"),
    ")"
  )
)
band_plot_df$panel <- factor(band_plot_df$panel, levels = panel_levels)

p_band_cooccurrence <- ggplot(
  band_plot_df,
  aes(x = b1017, y = b1395, color = b1575)
) +
  geom_hline(
    yintercept = 0, color = "grey55",
    linetype = "dashed", linewidth = 0.5
  ) +
  geom_vline(
    xintercept = 0, color = "grey55",
    linetype = "dashed", linewidth = 0.5
  ) +
  geom_point(size = 2.8, alpha = 0.78) +
  facet_wrap(~panel, ncol = 2) +
  scale_color_viridis_c(
    name = expression("1575 cm"^{-1}*" score"),
    option = "C",
    guide = guide_colorbar(
      direction = "horizontal",
      title.position = "top",
      title.hjust = 1,
      label.position = "bottom",
      barwidth = grid::unit(6, "cm"),
      barheight = grid::unit(0.45, "cm")
    )
  ) +
  labs(
    title = "Ca-DPA-associated band relationships",
    x = expression("1017 cm"^{-1}*" score (counts)"),
    y = expression("1395 cm"^{-1}*" score (counts)")
  ) +
  theme_bw(base_size = 16) +
  theme(
    plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 15, color = "black"),
    strip.text = element_text(size = 15, face = "bold"),
    strip.background = element_rect(fill = "grey95", color = "grey35"),
    legend.position = "bottom",
    legend.justification = "right",
    legend.box.just = "right",
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 13),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.4)
  )

print(p_band_cooccurrence)
ggsave(
  file.path(out_dir, "Figure3_band_relationships.png"),
  p_band_cooccurrence, width = 12, height = 8, dpi = 600
)
ggsave(
  file.path(out_dir, "Figure3_band_relationships.pdf"),
  p_band_cooccurrence, width = 12, height = 8
)

ratio_hist_df <- ratio_df %>%
  mutate(
    group_label = factor(
      label_group(group),
      levels = c("Spore-spiked soil", "Soil blank")
    )
  ) %>%
  filter(is.finite(ratio_1017_ch))

x_min <- min(ratio_hist_df$ratio_1017_ch)
x_max <- max(ratio_hist_df$ratio_1017_ch)
x_pad <- 0.05 * (x_max - x_min)

p_ratio_distribution <- ggplot(
  ratio_hist_df,
  aes(x = ratio_1017_ch)
) +
  geom_histogram(
    aes(y = after_stat(density), fill = group_label),
    binwidth = 0.025,
    boundary = 0,
    position = "identity",
    alpha = 0.4,
    color = "white",
    linewidth = 0.25
  ) +
  geom_density(
    aes(color = group_label),
    linewidth = 1.1,
    linetype = "dashed",
    adjust = 1.1,
    show.legend = FALSE
  ) +
  geom_vline(
    xintercept = 0,
    color = "grey40",
    linetype = "dotted",
    linewidth = 0.6
  ) +
  scale_fill_manual(
    name = "Sample group",
    values = group_colors,
    labels = c("Spore-spiked soil", "Soil blank")
  ) +
  scale_color_manual(
    values = group_colors,
    guide = "none"
  ) +
  coord_cartesian(xlim = c(x_min - x_pad, x_max + x_pad)) +
  labs(
    title = "1017/CH score distribution",
    x = "1017/CH score",
    y = "Density"
  ) +
  theme_bw(base_size = 15) +
  theme(
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 15, face = "bold"),
    axis.text = element_text(size = 15, color = "black"),
    legend.position = "bottom",
    legend.title = element_text(size = 15, face = "bold"),
    legend.text = element_text(size = 15),
    legend.key.size = grid::unit(0.7, "cm"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.4)
  )

print(p_ratio_distribution)
ggsave(
  file.path(out_dir, "Figure4_1017_CH_distribution.png"),
  p_ratio_distribution, width = 10, height = 7, dpi = 600
)
ggsave(
  file.path(out_dir, "Figure4_1017_CH_distribution.pdf"),
  p_ratio_distribution, width = 10, height = 7
)
