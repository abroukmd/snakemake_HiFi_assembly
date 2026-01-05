#!/usr/bin/env Rscript

library(data.table)
library(dplyr)
library(zoo)
library(ggplot2)
library(ggforce)
library(RColorBrewer)
library(tidyr)

# ---------------------------
# Argument parsing
# ---------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 4) {
  stop("Usage: ploidy_plot.R <sample_name> <fai> <depth> <output_pdf>",
       call. = FALSE)
}

name      <- args[1]
fai_file  <- args[2]
depth_file<- args[3]
out_file  <- args[4]

cat("Processing:", name, "\n")
cat("FAI:", fai_file, "\n")
cat("DEPTH:", depth_file, "\n")
cat("OUTPUT:", out_file, "\n")

# ---------------------------
# Color map
# ---------------------------
cols <- c("grey70", "#1f78b4", "#33a02c", "#ff7f00", "#e31a1c")
names(cols) <- c("Not present","n","2n","3n","4n")

# ---------------------------
# Helper functions
# ---------------------------
mode_fn <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

hysteresis_assign <- function(x, tol = 0.15) {
  out <- x
  last <- x[1]
  for (i in seq_along(x)) {
    if (abs(x[i] - last) <= tol) {
      out[i] <- last
    } else {
      last <- x[i]
      out[i] <- last
    }
  }
  return(out)
}

enforce_min_run <- function(x, min_len = 3) {  # patched: was 7
  r <- rle(x)
  idx <- rep(r$lengths >= min_len, r$lengths)
  x[!idx] <- NA
  fill <- zoo::na.locf(x, na.rm = FALSE)
  fill <- zoo::na.locf(fill, fromLast = TRUE, na.rm = FALSE)
  return(fill)
}

# ---------------------------
# Load FAI
# ---------------------------
chr <- fread(fai_file, header = FALSE)[, .(chromosome = V1, length = V2)]
if (nrow(chr) == 0)
  stop("FAI file contains zero contigs; cannot plot.")

chr_order <- chr$chromosome
genome_size <- sum(chr$length)

# ---------------------------
# Load depth
# ---------------------------
input <- fread(depth_file, sep="\t",
               col.names = c("chromosome","position","depth"))

if (nrow(input) == 0)
  stop("Depth file has zero rows — cannot plot.")

input$depth <- as.integer(input$depth)

# ---------------------------
# Compute bin sizes
# ---------------------------
target_bins <- 10000

chr <- chr %>%
  mutate(
    bin_size = round(length / (target_bins * (length / genome_size))),
    bin_size = pmin(bin_size, length),
    bin_size = pmax(bin_size, 500)     # patched: min 500bp (was 1000)
  )

# ---------------------------
# Binning
# ---------------------------
binned_list <- list()

for (chr_name in chr$chromosome) {
  chr_data <- input[input$chromosome == chr_name]

  if (nrow(chr_data) == 0)
    next

  bin_size <- chr$bin_size[chr$chromosome == chr_name]

  tmp <- chr_data %>%
    mutate(bin = floor(position / bin_size) * bin_size) %>%
    group_by(chromosome, bin) %>%
    summarise(mean_depth = mean(depth, na.rm = TRUE), .groups = "drop")

  binned_list[[chr_name]] <- tmp
}

if (length(binned_list) == 0)
  stop("No depth bins created — depth may not match FAI contigs.")

binned_depth <- bind_rows(binned_list) %>%
  arrange(chromosome, bin)

cat("DEBUG: bins per chromosome:\n\n")
print(table(binned_depth$chromosome))

cat("\nDEBUG: head(binned_depth):\n")
print(head(binned_depth))

# ---------------------------
# Smoothing & ploidy
# ---------------------------
# patched: k=5 instead of 15
binned_depth <- binned_depth %>%
  group_by(chromosome) %>%
  mutate(smoothed_depth = rollmedian(mean_depth, k = 5, fill = NA,
                                     align = "center")) %>%
  ungroup()

# fallback if smoothing erased too much
if (all(is.na(binned_depth$smoothed_depth))) {
  cat("WARNING: rollmedian produced all NA — using raw depth\n")
  binned_depth$smoothed_depth <- binned_depth$mean_depth
}

avg <- mean(binned_depth$mean_depth, na.rm = TRUE) / 2

binned_depth <- binned_depth %>%
  mutate(ploidy = smoothed_depth / avg) %>%
  filter(ploidy <= 6)   # keep reasonable tail

# ---------------------------
# Ploidy state calling
# ---------------------------
binned_depth <- binned_depth %>%
  arrange(chromosome, bin) %>%
  group_by(chromosome) %>%
    mutate(
      ploidy_num  = pmax(0, pmin(4, ploidy)),
      state_raw   = round(ploidy_num),
      state_hyst  = hysteresis_assign(ploidy_num, tol = 0.20),   # patched: tol=0.20
      state_mode  = rollapply(state_hyst, width = 7, FUN = mode_fn,
                              fill = NA, align = "center"),      # patched: width 7
      state_final = enforce_min_run(state_mode, min_len = 3)     # patched: min_len=3
    ) %>%
  ungroup() %>%
  mutate(col_final = factor(state_final, levels = 0:4, labels = names(cols)))

# ---------------------------
# Join chromosome lengths
# ---------------------------
binned_depth <- binned_depth %>%
  left_join(chr, by="chromosome") %>%
  mutate(chromosome = factor(chromosome, levels = chr_order))

cat("DEBUG: rows remaining:", nrow(binned_depth), "\n")
cat("DEBUG: chromosomes in final dataset:\n")
print(levels(binned_depth$chromosome))

cat("DEBUG: First few rows:\n")
print(head(binned_depth))

# ---------------------------
# Fallback: allow empty NA states (old script behavior)
# ---------------------------
if (nrow(binned_depth) == 0) {
  cat("\nERROR: No chromosomes remain after filtering — writing diagnostic PDF\n")

  pdf(out_file, width = 12, height = 5)
  plot.new()
  title(main = paste("No ploidy signal detected for", name))
  text(0.5, 0.5,
       "All bins eliminated after filtering.\nThis genome may be too fragmented or smoothing too strong.",
       cex = 1.2)
  dev.off()
  quit(save = "no", status = 0)
}

# ---------------------------
# Plot
# ---------------------------
p <- ggplot(binned_depth, aes(x = bin, y = ploidy)) +
  geom_blank(aes(x = length)) +
  geom_point(aes(color = col_final), size = 0.6, na.rm = FALSE) +  # keep NA
  scale_color_manual(values = cols, name = "Estimated ploidy", drop = FALSE) +
  coord_cartesian(ylim = c(0, 7)) +
  labs(x = "Position (bp)", y = "Ploidy", title = name) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 14),
    axis.text.x = element_blank(),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    plot.title  = element_text(size = 24, hjust = 0.5),
    legend.text = element_text(size = 18),
    legend.title= element_text(size = 20, face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  ) +
  guides(color = guide_legend(ncol = 4, override.aes = list(size = 6))) +
  ggforce::facet_row(~ chromosome, scales = "free_x", space = "free")

# ---------------------------
# Write PDF
# ---------------------------
ggsave(out_file, p, width = 18, height = 10, dpi = 300)
