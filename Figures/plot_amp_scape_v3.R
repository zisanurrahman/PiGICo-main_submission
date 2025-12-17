#!/usr/bin/env Rscript

# =================== Paths ===================
root <- "/Users/zisan/SGI_Paper/AMP/amp_scape_output"
edges_path    <- file.path(root, "edges.tsv")
nodes_path    <- file.path(root, "nodes.tsv")
clusters_path <- file.path(root, "clusters.tsv")
hits_path_1   <- file.path(root, "ampsphere_best_hits_known.tsv")  # preferred
hits_path_2   <- file.path(root, "ampsphere_best_hits.tsv")        # fallback
hits_path     <- if (file.exists(hits_path_1)) hits_path_1 else hits_path_2
out_prefix    <- file.path(root, "amp_scape_topN_plus_novel_clean")

# =================== Tuning ===================
min_weight_scaled <- 0.20       # prune weak edges (0..1 scale)
top_n_components  <- 35         # keep N largest + all with ≥1 Novel node
point_size_range  <- c(1.4, 5.5) # smaller nodes to reduce clutter
label_novel       <- FALSE      # set TRUE to label novel nodes

# =================== Packages ===================
need <- c("dplyr","readr","igraph","tidygraph","ggraph","graphlayouts","ggplot2","scales","ggrepel","grid")
for (p in need) if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = "https://cloud.r-project.org")
library(dplyr); library(readr); library(igraph); library(tidygraph); library(ggraph); library(graphlayouts)
library(ggplot2); library(scales); library(ggrepel); library(grid)

# =================== Load ===================
stopifnot(file.exists(edges_path), file.exists(nodes_path), file.exists(clusters_path))
edges <- suppressWarnings(read_tsv(edges_path, show_col_types = FALSE))
nodes <- suppressWarnings(read_tsv(nodes_path, show_col_types = FALSE))
clus  <- suppressWarnings(read_tsv(clusters_path,
                                   col_names = c("id","cluster"),
                                   col_types = cols(id = col_character(), cluster = col_character()),
                                   show_col_types = FALSE))
if (!all(c("q","t","weight") %in% names(edges))) stop("edges.tsv must have q,t,weight")
if (!("id" %in% names(nodes))) stop("nodes.tsv must have id")

# Hits (lookup approach)
hits <- if (file.exists(hits_path)) {
  suppressWarnings(read_tsv(hits_path, show_col_types = FALSE,
                            col_types = cols(query = col_character(),
                                             target = col_character(),
                                             pid = col_double(),
                                             status = col_character(),
                                             .default = col_guess())))
} else NULL

# =================== Normalize + attach clusters ===================
nodes$id <- trimws(as.character(nodes$id))
edges$q  <- trimws(as.character(edges$q))
edges$t  <- trimws(as.character(edges$t))
nodes <- nodes %>% left_join(clus, by = "id")

# =================== Build Present/Novel via lookup; drop no-hit ===================
nr <- nrow(nodes)
nodes$kn_status3    <- factor(rep("no_hit", nr), levels = c("known","novel","no_hit"))
nodes$ampsphere_pid <- rep(NA_real_, nr)

if (!is.null(hits)) {
  hits$query  <- trimws(as.character(hits$query))
  hit_map <- hits %>%
    transmute(query,
              status = ifelse(tolower(status) == "known", "known", "novel"),
              pid    = as.numeric(pid))
  idx <- match(nodes$id, hit_map$query)
  in_hits <- !is.na(idx)
  nodes$kn_status3[in_hits]    <- factor(hit_map$status[idx[in_hits]], levels = c("known","novel","no_hit"))
  nodes$ampsphere_pid[in_hits] <- hit_map$pid[idx[in_hits]]
}

# Collapse to two classes & DROP no_hit before building the graph
nodes <- nodes %>%
  mutate(kn_status2 = case_when(
    kn_status3 == "known" ~ "Present",
    kn_status3 == "novel" ~ "Novel",
    TRUE                  ~ NA_character_
  )) %>%
  filter(!is.na(kn_status2))  # <- drop no_hit entirely

# Source -> shape
nodes$source_group <- dplyr::case_when(
  grepl("PiBac", nodes$id) ~ "PiBac",
  grepl("SGI",   nodes$id) ~ "SGI",
  TRUE                     ~ "Other"
)
shape_map <- c(PiBac = 16, SGI = 17, Other = 18)  # circle, triangle, diamond

# =================== Edge scaling & keep edges among kept nodes ===================
wmax <- suppressWarnings(max(edges$weight, na.rm = TRUE)); if (!is.finite(wmax)) wmax <- 1
edges <- edges %>% mutate(weight_scaled = if (wmax <= 1.5) weight else weight/100)
edges <- edges %>%
  filter(weight_scaled >= min_weight_scaled) %>%
  filter(q %in% nodes$id & t %in% nodes$id)

# =================== Build graph (no no-hit nodes) ===================
g <- graph_from_data_frame(
  d = edges %>% select(q, t, weight_scaled) %>% rename(weight = weight_scaled),
  directed = FALSE,
  vertices = nodes %>% distinct(id, .keep_all = TRUE) %>% rename(name = id)
)

tg <- as_tbl_graph(g) %>%
  activate(nodes) %>%
  mutate(
    deg = centrality_degree(),
    pr  = centrality_pagerank(),
    cc  = group_components()
  )

# =================== Keep top N components + any with ≥1 Novel ===================
nd_cc <- as_tibble(tg, active = "nodes") %>%
  mutate(is_novel = kn_status2 == "Novel") %>%
  group_by(cc) %>%
  summarise(n = n(), has_novel = any(is_novel), .groups = "drop")

keep_top   <- nd_cc %>% arrange(desc(n)) %>% slice_head(n = top_n_components) %>% pull(cc)
keep_novel <- nd_cc %>% filter(has_novel) %>% pull(cc)
keep_cc    <- union(keep_top, keep_novel)

tg <- tg %>% activate(nodes) %>% filter(cc %in% keep_cc)

# =================== Node size = %ID (fallback degree) ===================
# =================== Node size = %ID to AMPsphere (fallback degree) ===================
nd <- as_tibble(tg, active = "nodes")
pid_vals <- suppressWarnings(as.numeric(nd$ampsphere_pid))
if (max(pid_vals, na.rm = TRUE) <= 1.5) pid_vals <- pid_vals * 100  # convert fractions to percent
pid_vals[!is.finite(pid_vals)] <- NA_real_
pid_vals <- pmin(pmax(pid_vals, 0), 100)  # clamp to [0,100]

# if missing, fallback to degree (rescaled 0–100 for consistency)
if (all(is.na(pid_vals))) {
  pid_vals <- rescale(nd$deg, to = c(0, 100))
}

tg <- tg %>% activate(nodes) %>% mutate(ampsphere_pid = pid_vals)


# =================== Layout (weighted) ===================
set.seed(42)
lay <- tryCatch(create_layout(tg, layout = "igraph", algorithm = "fr"),
                error = function(e) create_layout(tg, layout = "fr"))

# =================== Plot ===================
pal <- c(Present = "#1F77B4", Novel = "#1B9E77")
#"#1F77B4", RF = "darkslateblue", XGB = "#1B9E77")

title_txt <- sprintf(
  "AMP network — top %d components + all novel-containing components, min edge weight ≥ %.2f",
  top_n_components, min_weight_scaled
)

p <- ggraph(lay) +
  geom_edge_link0(aes(width = weight), colour = "grey", alpha = 1, show.legend = FALSE) +
  scale_edge_width(range = c(0.1, 0.8), guide = "none")+
  geom_node_point(
    aes(color = kn_status2, shape = source_group, size = ampsphere_pid),
    stroke = 0.1,   # border thickness
    fill = NA       # no fill, just outline
  )+
  scale_color_manual(values = pal, name = "AMPSphere") +
  scale_shape_manual(values = shape_map, name = "Source") +
  scale_size_continuous(
    name = "ID to AMPSphere",
    range = point_size_range,
    limits = c(0, 100),
    breaks = c(25, 50, 75, 100),
    labels = function(x) paste0(x, "%")
  )+
  theme_void(base_size = 11) +
  theme(
    legend.position   = "right",
    legend.key.size   = unit(0.9, "cm"),
    plot.title        = element_text(face = "bold")
  )+
  theme(legend.position = "top", legend.direction = "horizontal", legend.title = element_text(size = 12, face = 'bold'),
        legend.text = element_text(size = 10)) +
  guides(
    color = guide_legend(nrow = 1, byrow = TRUE),
    shape = guide_legend(nrow = 1, byrow = TRUE),
    size  = guide_legend(nrow = 1, byrow = TRUE)
  )


  #ggtitle(title_txt)
print(p)
if (label_novel) {
  p <- p + geom_node_text(
    aes(label = name),
    data = function(d) d %>% dplyr::filter(kn_status2 == "Novel"),
    size = 2.7, colour = "black", repel = TRUE, box.padding = 0.3, point.padding = 0.2,
    segment.size = 0.2
  )
}
print(p)
ggsave(paste0(out_prefix, "_plot.png"), p, width = 7.5, height = 6, dpi = 600, bg = "white")
ggsave(paste0(out_prefix, "_plot.pdf"), p, width = 13, height = 9)
message("[saved] ", paste0(out_prefix, "_plot.png"))
message("[saved] ", paste0(out_prefix, "_plot.pdf"))

# =================== Console summary ===================
cat("\n=== Plot summary ===\n")
cat(sprintf("Nodes: %d  Edges: %d\n", gorder(tg), gsize(tg)))
print(as_tibble(tg, active = "nodes") %>% count(kn_status2, sort = TRUE) %>% rename(count = n))
print(as_tibble(tg, active = "nodes") %>% count(source_group, sort = TRUE) %>% rename(count = n))


# =================== Number big networks (>5), label, and export members ===================

library(dplyr); library(tidyr); library(readr); library(ggplot2); library(ggpubr); library(purrr)

# 1) Component sizes and numbering (only components with >5 nodes)
nodes_tbl <- as_tibble(tg, active = "nodes") %>% select(name, cc)
cc_sizes  <- nodes_tbl %>% count(cc, name = "size")
big_cc    <- cc_sizes %>% filter(size > 4) %>% arrange(desc(size)) %>%
  mutate(network_number = row_number()) %>% select(cc, size, network_number)

nodes_tbl <- nodes_tbl %>% left_join(big_cc, by = "cc")

# 2) Annotate numbers on the plot (centroid per big component)
# use the layout already computed: 'lay' has x/y per node in same order as nodes_tbl
# After you build nodes_tbl and lay_df / lab_pos:

lay_df <- as.data.frame(lay)
stopifnot(nrow(lay_df) == nrow(nodes_tbl))
lay_df$network_number <- nodes_tbl$network_number  # attach the numbers to node positions

lab_pos <- lay_df %>%
  dplyr::filter(!is.na(network_number)) %>%
  dplyr::group_by(network_number) %>%
  dplyr::summarise(x = mean(x), y = mean(y), .groups = "drop")

p_labeled <- p +
  ggplot2::geom_text(
    data = lay_df %>% dplyr::filter(!is.na(network_number)),
    aes(x = x, y = y, label = network_number),
    size = 3.2, fontface = "bold", color = "black",
    inherit.aes = FALSE
  ) +
  ggplot2::geom_label(
    data = lab_pos,
    aes(x = x, y = y, label = network_number),
    size = 3.4, fontface = "bold",
    fill = scales::alpha("white", 0.7),
    color = "black",
    label.size = 0.2,
    inherit.aes = FALSE
  )


print(p_labeled)
ggsave(paste0(out_prefix, "_plot_NUMBERED.png"), p_labeled, width = 7.5, height = 6, dpi = 600, bg = "white")
ggsave(paste0(out_prefix, "_plot_NUMBERED.pdf"), p_labeled, width = 13, height = 9)

# 3) Export members for numbered networks
members_csv <- paste0(out_prefix, "_network_members_gt5.csv")
members_df <- nodes_tbl %>%
  filter(!is.na(network_number)) %>%
  transmute(network_number, node_id = name, component_id = cc, component_size = big_cc$size[match(cc, big_cc$cc)]) %>%
  arrange(network_number, node_id)
readr::write_csv(members_df, members_csv)

message("[saved] ", members_csv)

# =================== Join RPKM and build per-network abundance & q0 ===================

# ---- paths to RPKM (adjust if needed) ----
# =================== Join RPKM via SGI_/PiBac_ TAG and build per-network abundance & q0 ===================

suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(readr); library(ggplot2); library(purrr); library(scales) })

# ---- paths to RPKM (adjust if needed) ----
PWD_dir     <- "/Users/zisan/SGI_Paper/results/AlignmentResults_Guelph"  # PWD (unhealthy)
Healthy_dir <- "/Users/zisan/SGI_Paper/results/AlignmentResults_PigCat"  # Healthy

# read all rpkm.tsv files under a directory
read_rpkm_dir <- function(dir_path, group_label) {
  files <- list.files(dir_path, pattern = "rpkm\\.tsv$", full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
  if (!length(files)) return(tibble(sample = character(), group = character(), Contig = character(), RPKM = double()))
  map_dfr(files, function(f) {
    df <- suppressWarnings(readr::read_tsv(f, show_col_types = FALSE))
    cn <- names(df)
    contig_col <- cn[grepl("^contig$", cn, ignore.case = TRUE)][1]
    rpkm_col   <- cn[grepl("^rpkm$",   cn, ignore.case = TRUE)][1]
    if (is.na(contig_col) || is.na(rpkm_col)) return(NULL)
    tibble(
      sample = basename(f) %>% sub("\\.rpkm\\.tsv$", "", ., ignore.case = TRUE),
      group  = group_label,
      Contig = as.character(df[[contig_col]]),
      RPKM   = suppressWarnings(as.numeric(df[[rpkm_col]]))
    )
  })
}

rpkm_pwd <- read_rpkm_dir(PWD_dir, "PWD")
rpkm_hth <- read_rpkm_dir(Healthy_dir, "Healthy")
rpkm_all <- bind_rows(rpkm_pwd, rpkm_hth) %>% filter(is.finite(RPKM))

# ----- KEY FIX: match on SGI_/PiBac_ tag -----
# From node_id (e.g., "smORF_105|SGI_469"), take the token after '|', then canonicalize to SGI_### or PiBac_###
node_tag_map <- members_df %>%
  mutate(tag = sub("^.*\\|", "", node_id)) %>%                              # after last '|'
  mutate(tag = sub("^((SGI|PiBac)_[0-9]+).*$", "\\1", tag, perl = TRUE)) %>%
  distinct(network_number, node_id, tag)

# From Contig (e.g., "PiBac_003_16", "SGI_469_2"), extract first SGI_/PiBac_ token
rpkm_all2 <- rpkm_all %>%
  mutate(tag = sub("^.*?((SGI|PiBac)_[0-9]+).*$", "\\1", Contig, perl = TRUE)) %>%
  filter(tag %in% node_tag_map$tag)

# Join and aggregate to node, then to network per sample
node_abund <- rpkm_all2 %>%
  inner_join(node_tag_map, by = "tag") %>%
  group_by(sample, group, network_number, node_id) %>%
  summarise(node_abundance = sum(RPKM, na.rm = TRUE), .groups = "drop")

net_abund <- node_abund %>%
  group_by(sample, group, network_number) %>%
  summarise(rel_abundance = sum(node_abundance, na.rm = TRUE), .groups = "drop")

# q0 richness within network = number of nodes present (>0) per sample
net_q0 <- node_abund %>%
  group_by(sample, group, network_number) %>%
  summarise(q0 = sum(node_abundance > 0, na.rm = TRUE), .groups = "drop")

# Save long tables
abund_csv <- paste0(out_prefix, "_network_abundance_long.csv")
q0_csv    <- paste0(out_prefix, "_network_q0_long.csv")
readr::write_csv(net_abund, abund_csv)
readr::write_csv(net_q0,    q0_csv)
message("[saved] ", abund_csv)
message("[saved] ", q0_csv)

# =================== Plots: abundance and q0 per numbered network ===================

keep_nets <- sort(unique(members_df$network_number))

abund_plot_df <- net_abund %>%
  filter(network_number %in% keep_nets) %>%
  mutate(
    group = factor(group, levels = c("Healthy","PWD")),
    network_number = factor(network_number, levels = keep_nets)
  )

q0_plot_df <- net_q0 %>%
  filter(network_number %in% keep_nets) %>%
  mutate(
    group = factor(group, levels = c("Healthy","PWD")),
    network_number = factor(network_number, levels = keep_nets)
  )

# Guard + quick diagnostics
message("abund_plot_df rows: ", nrow(abund_plot_df))
message("q0_plot_df rows: ", nrow(q0_plot_df))
if (nrow(abund_plot_df) == 0L || nrow(q0_plot_df) == 0L) {
  message("\n--- DEBUG ---")
  message("Examples mapped tags (RPKM): ", paste(utils::head(unique(rpkm_all2$tag), 10), collapse = ", "))
  message("Examples mapped tags (nodes): ", paste(utils::head(unique(node_tag_map$tag), 10), collapse = ", "))
  stop("No data to facet: no matched RPKM to numbered networks. Check tag extraction or file paths.")
}

# 1) Relative abundance (Healthy vs PWD), faceted by network, free y, with panel border
p_ab <- ggplot(abund_plot_df, aes(x = group, y = rel_abundance, fill = group)) +
  geom_violin(trim = FALSE, width = 0.9, alpha = 0.75, color = NA) +
  geom_boxplot(width = 0.15, outlier.size = 0.7, alpha = 0.95) +
  facet_wrap(~ network_number, scales = "free_y",
             labeller = labeller(network_number = function(v) paste0("Network ", v))) +
  scale_fill_manual(values = c(Healthy = "#4C78A8", PWD = "#F58518")) +
  labs(x = NULL, y = expression("Abundance (Sample"^-1*")"), fill = NULL,
       title = "AMP network abundance per sample (summed RPKM per network)") +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "top",
    strip.text = element_text(size = 9, face = "bold"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  )
print(p_ab)
ggsave(paste0(out_prefix, "_network_abundance_faceted.png"), p_ab, width = 12, height = 9, dpi = 300, bg = "white")
ggsave(paste0(out_prefix, "_network_abundance_faceted.pdf"), p_ab, width = 12, height = 9, bg = "white")

# 2) q0 richness (nodes present) per network, same faceting and style
p_q0 <- ggplot(q0_plot_df, aes(x = group, y = q0, fill = group)) +
  geom_violin(trim = FALSE, width = 0.9, alpha = 0.75, color = NA) +
  geom_boxplot(width = 0.15, outlier.size = 0.7, alpha = 0.95) +
  facet_wrap(~ network_number, scales = "free_y",
             labeller = labeller(network_number = function(v) paste0("Network ", v))) +
  scale_fill_manual(values = c(Healthy = "#4C78A8", PWD = "#F58518")) +
  labs(x = NULL, y = "q0 richness (nodes present)", fill = NULL,
       title = "AMP network richness (q0) per sample") +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "top",
    strip.text = element_text(size = 9, face = "bold"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  )
print(p_q0)
ggsave(paste0(out_prefix, "_network_q0_faceted.png"), p_q0, width = 12, height = 9, dpi = 300, bg = "white")
ggsave(paste0(out_prefix, "_network_q0_faceted.pdf"), p_q0, width = 12, height = 9, bg = "white")

# 3) Also save per-network individual panels
save_one_panel <- function(df, yvar, ylab, stub) {
  nets <- sort(unique(df$network_number))
  for (nn in nets) {
    d <- df %>% filter(network_number == nn)
    p1 <- ggplot(d, aes(x = group, y = .data[[yvar]], fill = group)) +
      geom_violin(trim = FALSE, width = 0.9, alpha = 0.75, color = NA) +
      geom_boxplot(width = 0.15, outlier.size = 0.7, alpha = 0.95) +
      scale_fill_manual(values = c(Healthy = "#4C78A8", PWD = "#F58518")) +
      labs(title = paste0("Network ", nn), x = NULL, y = ylab, fill = NULL) +
      theme_minimal(base_size = 11) +
      theme(
        legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
      )
    ggsave(sprintf("%s_%03d.png", stub, as.integer(nn)), p1, width = 4.0, height = 3.2, dpi = 300, bg = "white")
    ggsave(sprintf("%s_%03d.pdf", stub, as.integer(nn)), p1, width = 4.0, height = 3.2, bg = "white")
  }
}

save_one_panel(abund_plot_df, "rel_abundance", expression("Abundance (Sample"^-1*")"),
               paste0(out_prefix, "_network_abundance_panel"))
save_one_panel(q0_plot_df, "q0", "q0 richness (nodes present)",
               paste0(out_prefix, "_network_q0_panel"))

# =================== Combo abundance + richness (q0) with median slope lines ===================

SWAP_LABELS <- TRUE  # set TRUE to display Healthy as "PWD" and PWD as "Healthy"

# Keep only numbered networks
keep_nets <- sort(unique(members_df$network_number))

# Merge abundance & q0 per sample×group×network
mat <- net_abund %>%
  inner_join(net_q0, by = c("sample","group","network_number")) %>%
  filter(network_number %in% keep_nets) %>%
  mutate(
    group_disp = if (SWAP_LABELS) dplyr::recode(group, Healthy = "PWD", PWD = "Healthy") else group
  )

# Scale q0 to the abundance axis within each network so both can be plotted on the same y
scales_tbl <- mat %>%
  group_by(network_number) %>%
  summarise(
    maxA = max(rel_abundance, na.rm = TRUE),
    maxQ = max(q0,            na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    maxA = ifelse(is.finite(maxA) & maxA > 0, maxA, 1),
    maxQ = ifelse(is.finite(maxQ) & maxQ > 0, maxQ, 1)
  )

mat <- mat %>%
  left_join(scales_tbl, by = "network_number") %>%
  mutate(q0_scaled = q0 / maxQ * maxA)

# Medians per network×group for slope lines
med_tbl <- mat %>%
  group_by(network_number, group_disp) %>%
  summarise(
    med_abund     = median(rel_abundance, na.rm = TRUE),
    med_q0_scaled = median(q0_scaled,     na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(x = ifelse(group_disp == "Healthy", 1, 2))  # x positions for slopes

slope_abund <- med_tbl %>%
  select(network_number, group_disp, x, y = med_abund) %>%
  arrange(network_number, x)

slope_q0 <- med_tbl %>%
  select(network_number, group_disp, x, y = med_q0_scaled) %>%
  arrange(network_number, x)

# Direction annotations (which group higher by median)
dir_labels <- med_tbl %>%
  tidyr::pivot_wider(names_from = group_disp, values_from = c(med_abund, med_q0_scaled)) %>%
  mutate(
    ab_dir = case_when(
      `med_abund_PWD`     > `med_abund_Healthy`     ~ "Abundance: PWD higher",
      `med_abund_PWD`     < `med_abund_Healthy`     ~ "Abundance: Healthy higher",
      TRUE ~ "Abundance: Tie"
    ),
    q0_dir = case_when(
      `med_q0_scaled_PWD` > `med_q0_scaled_Healthy` ~ "Richness: PWD higher",
      `med_q0_scaled_PWD` < `med_q0_scaled_Healthy` ~ "Richness: Healthy higher",
      TRUE ~ "Richness: Tie"
    )
  ) %>%
  select(network_number, ab_dir, q0_dir)

# y positions for annotations per facet
ypos <- mat %>%
  group_by(network_number) %>%
  summarise(y.max = max(rel_abundance, na.rm = TRUE), .groups = "drop")

annot <- dir_labels %>%
  left_join(ypos, by = "network_number") %>%
  mutate(x = 1.5, y1 = y.max * 1.05, y2 = y.max * 1.13)

# Facet label helper
net_lab <- function(v) paste0("Network ", v)

# =================== PLOT (faceted) ===================
p_combo <- ggplot(mat, aes(x = group_disp, y = rel_abundance)) +
  # abundance distribution
  geom_violin(aes(fill = group_disp), width = 0.9, trim = FALSE, color = NA, alpha = 0.70) +
  geom_boxplot(width = 0.15, outlier.size = 0.7, alpha = 0.95) +
  # median slope: abundance (solid)
  geom_line(data = slope_abund, aes(x = x, y = y, group = network_number, color = "Median abundance"),
            linewidth = 0.5, inherit.aes = FALSE, na.rm = TRUE) +
  geom_point(data = slope_abund, aes(x = x, y = y, color = "Median abundance"),
             size = 1.5, inherit.aes = FALSE, shape = 16, na.rm = TRUE) +
  # median slope: richness q0 (scaled) (dashed)
  geom_line(data = slope_q0, aes(x = x, y = y, group = network_number, color = "Median richness (q0, scaled)"),
            linewidth = 0.5, linetype = "22", inherit.aes = FALSE, na.rm = TRUE) +
  geom_point(data = slope_q0, aes(x = x, y = y, color = "Median richness (q0, scaled)"),
             size = 1.5, inherit.aes = FALSE, shape = 17, na.rm = TRUE) +
  # direction annotations
  #geom_text(data = annot, aes(x = x, y = y1, label = ab_dir),
   #         size = 3.0, fontface = "bold", inherit.aes = FALSE) +
  #geom_text(data = annot, aes(x = x, y = y2, label = q0_dir),
  #          size = 3.0, fontface = "italic", inherit.aes = FALSE) +
  facet_wrap(~ network_number, scales = "free_y",
             labeller = labeller(network_number = net_lab)) +
  scale_fill_manual(values = c(Healthy = "darkslateblue", PWD = "#1F77B4")) +
  scale_color_manual(
    name = NULL,
    values = c("Median abundance" = "grey",
               "Median richness (q0, scaled)" = "black")
  ) +
  labs(
    x = NULL,
    y = expression("Abundance (Sample"^-1*")"),
    fill = NULL
  ) +
  theme_minimal(base_size = 08) +
  theme(
    legend.position = "top",
    legend.box = "horizontal",
    strip.text = element_text(size = 9, face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_text(size = 10),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  )

print(p_combo)
ggsave(paste0(out_prefix, "_network_combo_abundance_plus_q0Slope_faceted.png"),
       p_combo, width = 12, height = 9, dpi = 300, bg = "white")
ggsave(paste0(out_prefix, "_network_combo_abundance_plus_q0Slope_faceted.pdf"),
       p_combo, width = 12, height = 9, bg = "white")

# =================== Optional: save each network panel separately ===================
unique_nets <- sort(unique(mat$network_number))
for (nn in unique_nets) {
  dat_g   <- dplyr::filter(mat,          network_number == nn)
  sab_g   <- dplyr::filter(slope_abund,  network_number == nn)
  sq0_g   <- dplyr::filter(slope_q0,     network_number == nn)
  annot_g <- dplyr::filter(annot,        network_number == nn)
  
  p_one <- ggplot(dat_g, aes(x = group_disp, y = rel_abundance)) +
    geom_violin(aes(fill = group_disp), width = 0.9, trim = FALSE, color = NA, alpha = 0.70) +
    geom_boxplot(width = 0.15, outlier.size = 0.7, alpha = 0.95) +
    geom_line(data = sab_g, aes(x = x, y = y, group = network_number, color = "Median abundance"),
              linewidth = 0.5, inherit.aes = FALSE, na.rm = TRUE) +
    geom_point(data = sab_g, aes(x = x, y = y, color = "Median abundance"),
               size = 1.5, inherit.aes = FALSE, shape = 16, na.rm = TRUE) +
    geom_line(data = sq0_g, aes(x = x, y = y, group = network_number, color = "Median richness (q0, scaled)"),
              linewidth = 0.5, linetype = "22", inherit.aes = FALSE, na.rm = TRUE) +
    geom_point(data = sq0_g, aes(x = x, y = y, color = "Median richness (q0, scaled)"),
               size = 1.5, inherit.aes = FALSE, shape = 16, na.rm = TRUE) +
    #geom_text(data = annot_g, aes(x = x, y = y1, label = ab_dir),
     #         size = 3.0, fontface = "bold", inherit.aes = FALSE) +
    #geom_text(data = annot_g, aes(x = x, y = y2, label = q0_dir),
     #         size = 3.0, fontface = "italic", inherit.aes = FALSE) +
    scale_fill_manual(values = c(Healthy = "darkslateblue", PWD = "#1F77B4")) +
    scale_color_manual(name = NULL,
                       values = c("Median abundance" = "grey",
                                  "Median richness (q0, scaled)" = "black")) +
    labs(title = paste0("Network ", nn), x = NULL,
         y = expression("Abundance (Sample"^-1*")"), fill = NULL) +
    theme_minimal(base_size = 08) +
    theme(
      legend.position = "top",
      legend.box = "horizontal",
      strip.text = element_text(size = 9, face = "bold"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      axis.text.y = element_blank(),
      axis.title.y = element_text(size=07),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
    )
  
  ggsave(sprintf("%s_network_combo_abundance_plus_q0Slope_%03d.png", out_prefix, as.integer(nn)),
         p_one, width = 1.5, height = 2, dpi = 600, bg = "white")
  ggsave(sprintf("%s_network_combo_abundance_plus_q0Slope_%03d.pdf", out_prefix, as.integer(nn)),
         p_one, width = 4.2, height = 3.4, bg = "white")
}
