## Shared plotting helpers for the QC reports (tSNR, framewise displacement).
##
## Contract: the data frame passed in has columns
##   ID   subject factor
##   y    the dependent variable
##   grp  the factor shown on the x axis (day, or run within a day)
##   blk  optional second factor that groups the x axis into blocks (echo).
##        Omit it, or pass block = NULL, for a single-block plot.
## add_cells() builds the `cell` column that the panel actually plots.

p_stars <- function(p) {
  ifelse(p < 0.001, "***",
  ifelse(p < 0.01,  "**",
  ifelse(p < 0.05,  "*", NA_character_)))
}

# Lowest free vertical tier per bracket, so that overlapping spans stack
# instead of being drawn on top of each other.
pack_tiers <- function(xmin, xmax, pad = 0.35) {
  tier <- integer(length(xmin))
  for (i in seq_along(xmin)) {
    j <- seq_len(i - 1L)
    t <- 1L
    while (any(tier[j] == t & xmin[i] <= xmax[j] + pad & xmax[i] >= xmin[j] - pad)) {
      t <- t + 1L
    }
    tier[i] <- t
  }
  tier
}

# One box per grp (x blk) combination. With a block factor the order is
# block-major, so that the groups being compared sit next to each other.
add_cells <- function(d, block = NULL) {
  d <- dplyr::mutate(d, grp = droplevels(grp))
  if (is.null(block)) {
    return(dplyr::mutate(d, cell = factor(as.character(grp), levels = levels(d$grp))))
  }
  d   <- dplyr::mutate(d, blk = droplevels(.data[[block]]))
  lev <- paste0(rep(levels(d$blk), each = nlevels(d$grp)), "_", levels(d$grp))
  dplyr::mutate(d, cell = factor(paste0(as.character(blk), "_", as.character(grp)),
                                 levels = lev))
}

## Significance brackets -------------------------------------------------------
## Drawn only for comparisons that survive correction. `pw` is an rstatix
## pairwise table with group1, group2, p.adj, plus the block column if blocked.
##
## Heights are computed per block and just above that block's own data: with
## blocks on very different scales (echo 1 is ~3x echo 4) a single panel-wide
## height would strand the low block's brackets far above their boxes.
## With log_y = TRUE the stacking is multiplicative, so the gaps stay even.
sig_geom <- function(pw, d, block = NULL, log_y = FALSE,
                     step_frac = 0.055, tick_frac = 0.014) {
  grp <- levels(d$grp)
  tf  <- if (log_y) log10 else identity
  itf <- if (log_y) function(z) 10^z else identity

  span <- diff(range(tf(d$y)))
  step <- step_frac * span
  tick <- tick_frac * span

  if (is.null(block)) {
    d  <- dplyr::mutate(d, .blk = "all")
    pw <- dplyr::mutate(pw, .blk = "all")
  } else {
    d  <- dplyr::mutate(d,  .blk = as.character(.data[[block]]))
    pw <- dplyr::mutate(pw, .blk = as.character(.data[[block]]))
  }
  blk_levels <- unique(d$.blk)

  y_blk <- d |>
    dplyr::group_by(.blk) |>
    dplyr::summarise(y0 = max(tf(y)), .groups = "drop")

  b <- pw |>
    dplyr::ungroup() |>
    dplyr::mutate(label = p_stars(p.adj)) |>
    dplyr::filter(!is.na(label)) |>
    dplyr::mutate(
      e    = match(.blk, blk_levels),
      xmin = (e - 1) * length(grp) + pmin(match(group1, grp), match(group2, grp)),
      xmax = (e - 1) * length(grp) + pmax(match(group1, grp), match(group2, grp))
    ) |>
    dplyr::left_join(y_blk, by = ".blk") |>
    dplyr::arrange(e, xmax - xmin, xmin) |>
    dplyr::group_by(e) |>
    dplyr::mutate(tier = pack_tiers(xmin, xmax)) |>
    dplyr::ungroup() |>
    dplyr::mutate(y_t = y0 + tier * step, y = itf(y_t))

  if (nrow(b) == 0) return(list())

  # Four points per bracket: down-tick, across, down-tick.
  paths <- data.frame(
    id = rep(seq_len(nrow(b)), each = 4),
    x  = as.vector(rbind(b$xmin, b$xmin, b$xmax, b$xmax)),
    y  = itf(as.vector(rbind(b$y_t - tick, b$y_t, b$y_t, b$y_t - tick)))
  )

  list(
    geom_path(data = paths, aes(x = x, y = y, group = id),
              colour = "black", linewidth = 0.4, inherit.aes = FALSE),
    geom_text(data = b, aes(x = (xmin + xmax) / 2, y = y, label = label),
              vjust = -0.15, size = 4.6, colour = "black", inherit.aes = FALSE),
    expand_limits(y = itf(max(b$y_t) + 1.8 * step))
  )
}

## Theme ----------------------------------------------------------------------
qc_theme <- function(base_font = 14) {
  theme_classic(base_size = base_font) +
    theme(
      panel.grid    = element_blank(),
      axis.line     = element_line(colour = "black", linewidth = 0.5),
      axis.ticks    = element_line(colour = "black", linewidth = 0.5),
      axis.text     = element_text(colour = "black"),
      axis.title    = element_text(colour = "black"),
      axis.title.y  = element_text(margin = margin(r = 8)),
      plot.title    = element_text(face = "bold", size = base_font + 1),
      plot.subtitle = element_text(size = base_font - 2),
      plot.caption  = element_text(size = base_font - 3, hjust = 0, colour = "grey30"),
      plot.margin   = margin(6, 8, 14, 6)
    )
}

## Panel ----------------------------------------------------------------------
qc_panel <- function(d, y_lab, x_title, title, subtitle, caption,
                     pw = NULL, block = NULL, log_y = FALSE,
                     short = identity, block_prefix = "", base_font = 14) {
  n_g <- nlevels(d$grp)

  p <- ggplot(d, aes(x = cell, y = y)) +
    geom_boxplot(
      fill = "white", colour = "black",
      linewidth = 0.5, width = 0.62,
      outlier.shape = NA          # outliers are already shown as points below
    ) +
    # Jitter is horizontal only, so no point is drawn at a false value.
    geom_point(
      position = position_jitter(width = 0.16, height = 0, seed = 42),
      shape = 21, fill = "black", colour = "white",
      size = 2.1, stroke = 0.35, alpha = 0.9
    ) +
    scale_x_discrete(labels = function(z) short(sub("^[^_]+_", "", z))) +
    coord_cartesian(clip = "off") +
    labs(title = title, subtitle = subtitle, caption = caption,
         x = x_title, y = y_lab) +
    theme(axis.title.x = element_text(
      margin = margin(t = if (is.null(block)) 8 else 30)))

  if (log_y) {
    p <- p + scale_y_log10()
  } else {
    p <- p + scale_y_continuous(expand = expansion(mult = c(0.04, 0.06)))
  }

  # One label per block, centred under that block's boxes.
  if (!is.null(block)) {
    blk <- data.frame(
      x     = (seq_len(nlevels(d$blk)) - 1) * n_g + (n_g + 1) / 2,
      label = paste0(block_prefix, levels(d$blk))
    )
    p <- p + geom_text(data = blk, aes(x = x, y = -Inf, label = label),
                       vjust = 3.5, size = base_font / .pt, colour = "black",
                       inherit.aes = FALSE)
  }

  if (!is.null(pw)) p <- p + sig_geom(pw, d, block = block, log_y = log_y)
  p
}
