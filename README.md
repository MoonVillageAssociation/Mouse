# Mouse

Create contrast_map.csv to define data and contrasts
```
ID,GLDS,stat_col,flip,label
1,GLDS-606,Stat_(Wild Type & Space Flight)v(Wild Type & Ground Control),FALSE,606_colon_29d
2,GLDS-673,Stat_(Space Flight)v(Ground Control),FALSE,673_soleus_38d
3,GLDS-674,Stat_(Space Flight & 29 week & On ISS)v(Ground Control & 29 week & On Earth),FALSE,674_kidney_58d
4,GLDS-248,Stat_(Space Flight & ~60 day & On ISS & Carcass)v(Ground Control & ~60 day & On Earth & Carcass),FALSE,248_lung_56d
5,GLDS-589,Stat_(Space Flight & 29 week)v(Ground Control & 29 week),FALSE,589_brain_58d
```

## KEGG GSEA

Setup
```
library(tidyverse)
library(clusterProfiler)
```

Specify data to analyze from contrast_map.csv
```
TARGET_ID <- 4
```

Setup KEGG
```
organism_kegg <- "mmu"
entrez_col <- "ENTREZID"
```

Read contrast table
```
contrast_df <- read_csv("contrast_map.csv", show_col_types = FALSE) %>%
  filter(ID == TARGET_ID)

stopifnot(nrow(contrast_df) == 1)

TARGET_GLDS <- contrast_df$GLDS[1]
```

Function to run KEGG GSEA
```
run_kegg_gsea <- function(df, stat_col, flip = FALSE) {

  ranks <- df %>%
    filter(!is.na(.data[[stat_col]]),
           !is.na(.data[[entrez_col]])) %>%
    distinct(.data[[entrez_col]], .keep_all = TRUE) %>%
    mutate(stat_use = if (flip) - .data[[stat_col]] else .data[[stat_col]]) %>%
    select(!!entrez_col, stat_use) %>%
    deframe() %>%
    sort(decreasing = TRUE)

  gseKEGG(
    geneList     = ranks,
    organism     = organism_kegg,
    pvalueCutoff = 0.05,
    verbose      = FALSE
  )
}
```

Run
```
csv_file <- paste0(
  TARGET_GLDS,
  "_rna_seq_differential_expression_GLbulkRNAseq.csv"
)

message("Processing ", TARGET_GLDS)

df <- read_csv(csv_file, show_col_types = FALSE)

contrast <- contrast_df

message("  Contrast: ", contrast$label)

if (!contrast$stat_col %in% colnames(df)) {
  stop("Stat column not found: ", contrast$stat_col)
}

gsea <- run_kegg_gsea(
  df       = df,
  stat_col = contrast$stat_col,
  flip     = contrast$flip
)

if (!is.null(gsea) && nrow(gsea@result) > 0) {

  res <- as.data.frame(gsea@result)
  res$GLDS     <- TARGET_GLDS
  res$Contrast <- contrast$label
  res$StatUsed <- contrast$stat_col

  out_file <- paste0(
    TARGET_GLDS, "_", contrast$label, "_KEGG_GSEA.csv"
  )

  write_csv(res, out_file)
  message("    Saved: ", out_file)
} else {
  warning("No significant KEGG pathways")
}
```

## KEGG Gene Lists 
Expands significant KEGG pathways to full KEGG gene lists and filtered to expressed genes only

Setup
```
library(tidyverse)
library(KEGGREST)
```

Specify data to analyze from contrast_map.csv
```
TARGET_ID <- 4
```

Read contrast table
```
contrast_table <- "contrast_map.csv"

contrast <- read_csv(contrast_table, show_col_types = FALSE) %>%
  filter(ID == TARGET_ID)

stopifnot(nrow(contrast) == 1)

glds  <- contrast$GLDS
label <- contrast$label

kegg_file <- paste0(
  glds, "_", label, "_KEGG_GSEA.csv"
)

background_file <- paste0(
  glds, "_rna_seq_differential_expression_GLbulkRNAseq.csv"
)

message("Processing ", glds, " (", label, ")")
```

Load data
```
gsea <- read_csv(kegg_file, show_col_types = FALSE)

bg_genes <- read_csv(background_file, show_col_types = FALSE) %>%
  filter(!is.na(ENTREZID)) %>%
  pull(ENTREZID) %>%
  unique()
```

Fetch KEGG gene sets
```
kegg_full <- map_dfr(gsea$ID, function(pid) {

  kg <- tryCatch(keggGet(pid), error = function(e) NULL)
  if (is.null(kg)) return(NULL)

  genes <- kg[[1]]$GENE
  if (is.null(genes) || length(genes) < 2) return(NULL)

  entrez <- genes[seq(1, length(genes), by = 2)]

  tibble(
    ID = pid,
    ENTREZID = entrez
  )
}) %>%
  filter(ENTREZID %in% bg_genes) %>%
  distinct()

out_file <- paste0(
  glds, "_", label, "_KEGG_FULL_PATHWAY_GENES_EXPRESSED.csv"
)

write_csv(kegg_full, out_file)
message("Saved: ", out_file)
```


## Enrichment of element-interacting genes in KEGG pathways
Monte-Carlo enrichment of CTD element gene sets in KEGG pathways compared to a random background

Setup
```
library(tidyverse)
library(org.Mm.eg.db)
library(AnnotationDbi)
```

Specify data to analyze from contrast_map.csv
```
TARGET_ID <- 4
```

Specify metals of interest
```
metal_files <- c(
  "Al.csv",
  "B.csv",
  "Ba.csv",
  "Ca.csv",
  "Cd.csv",
  "Cr.csv",
  "Cu.csv",
  "Fe.csv",
  "K.csv",
  "Li.csv",
  "Mo.csv",
  "Na.csv",
  "Ni.csv",
  "P.csv",
  "Pb.csv",
  "Zn.csv"
)
```

Setup permutations
```
set.seed(123)
n_perm <- 1000
```

Read contrast table
```
contrast_table <- "contrast_map.csv"

contrast <- read_csv(contrast_table, show_col_types = FALSE) %>%
  filter(ID == TARGET_ID)

stopifnot(nrow(contrast) == 1)

glds  <- contrast$GLDS
label <- contrast$label
```

Input files
```
kegg_pathway_file <- paste0(
  glds, "_", label, "_KEGG_FULL_PATHWAY_GENES_EXPRESSED.csv"
)

kegg_gsea_file <- paste0(
  glds, "_", label, "_KEGG_GSEA.csv"
)

background_file <- paste0(
  glds, "_rna_seq_differential_expression_GLbulkRNAseq.csv"
)

message("Processing ", glds, " (", label, ")")
```

Background universe
```
bg <- read_csv(background_file, show_col_types = FALSE) %>%
  filter(!is.na(ENTREZID)) %>%
  pull(ENTREZID) %>%
  unique()
```

KEGG gene sets
```
kegg_sets <- read_csv(kegg_pathway_file, show_col_types = FALSE) %>%
  distinct(ID, ENTREZID)

 kegg_desc <- read_csv(kegg_gsea_file, show_col_types = FALSE) %>%
  dplyr::select(ID, Description) %>%
  distinct()
```

Check for enrichment of element-interacting genes
```
all_metal_results <- list()

for (metal_file in metal_files) {

  metal <- tools::file_path_sans_ext(basename(metal_file))

  metal_genes <- read_csv(metal_file, show_col_types = FALSE) %>%
    pull(`Gene ID`) %>%
    intersect(bg)

  res <- kegg_sets %>%
    group_by(ID) %>%
    summarise(
      pathway_genes = list(unique(ENTREZID)),
      n_pathway = length(pathway_genes),
      .groups = "drop"
    ) %>%
    rowwise() %>%
    mutate(
      observed_k = length(intersect(pathway_genes, metal_genes)),
      null_k = list({
        replicate(
          n_perm,
          length(intersect(
            sample(bg, n_pathway),
            metal_genes
          ))
        )
      }),
      p_value = mean(unlist(null_k) >= observed_k),
      GLDS = glds,
      Metal = metal
    ) %>%
    ungroup() %>%
    mutate(padj = p.adjust(p_value, method = "BH")) %>%
    filter(padj < 0.05) %>%
    left_join(kegg_desc, by = "ID") %>%
    dplyr::select(
      GLDS, Metal,
      ID, Description,
      n_pathway,
      observed_k,
      p_value, padj
    )

  all_metal_results[[metal]] <- res
}

final <- bind_rows(all_metal_results)

out_file <- paste0(
  glds, "_", label, "_KEGG_Pathway_Metal_MC_SIGNIFICANT.csv"
)

write_csv(final, out_file)
message("Saved: ", out_file)
```

## Filtering
Select pathways of interest from *_KEGG_Pathway_Metal_MC_SIGNIFICANT.csv
Generate custom lists with *_selected_kegg_pathways.csv

## Plot

Setup
```
library(tidyverse)
library(cowplot)
```

Specify data to analyze from contrast_map.csv
```
TARGET_ID <- 4
```

Read contrast table
```
contrast_table <- "contrast_map.csv"

contrast <- read_csv(contrast_table, show_col_types = FALSE) %>%
  filter(ID == TARGET_ID)

stopifnot(nrow(contrast) == 1)

glds_id <- contrast$GLDS
label   <- contrast$label

message("Processing ", glds_id, " (", label, ")")
```

Input
```
# KEGG GSEA results
kegg_csv <- paste0(
  glds_id, "_", label, "_KEGG_GSEA.csv"
)

# Metal enrichment results
metal_csv <- paste0(
  glds_id, "_", label, "_KEGG_Pathway_Metal_MC_SIGNIFICANT.csv"
)

# Hand-selected KEGG pathways (YOU curate this)
pathway_txt <- paste0(
  glds_id, "_", label, "_selected_kegg_pathways.csv"
)

# Output plot
out_file <- paste0(
  glds_id, "_", label, "_NES_plus_MetalHeatmap.png"
)
```
Load data
```
kegg  <- read_csv(kegg_csv, show_col_types = FALSE)
metal <- read_csv(metal_csv, show_col_types = FALSE)

pathway_csv <- paste0(
  glds_id, "_", label, "_selected_kegg_pathways.csv"
)

pathways_of_interest <- read_csv(pathway_csv, show_col_types = FALSE) %>%
  pull(ID) %>%
  unique()

stopifnot(length(pathways_of_interest) > 0)
```

Prep NES data
```nes_df <- kegg %>%
  filter(ID %in% pathways_of_interest) %>%
  mutate(
    Pathway = str_remove(
      Description,
      " - Mus musculus \\(house mouse\\)"
    )
  ) %>%
  arrange(NES)

# Enforce shared ordering between plots
pathway_levels <- nes_df$Pathway
```
Prep element data
```
metal_df <- metal %>%
  left_join(
    kegg %>% select(ID, Description),
    by = "ID"
  ) %>%
  filter(ID %in% pathways_of_interest) %>%
  mutate(
    Pathway = str_remove(
      Description.y,
      " - Mus musculus \\(house mouse\\)"
    ),
    Metal = factor(Metal, levels = c("Al", "Cu", "Fe", "Ti")),
    padj_plot = pmax(padj, 1e-300)
  ) %>%
  mutate(
    Pathway = factor(Pathway, levels = pathway_levels)
  )
```

NES barplot
```
fixed_margin <- margin(5, 5, 5, 5)

p_nes <- ggplot(
  nes_df,
  aes(x = NES, y = factor(Pathway, levels = pathway_levels))
) +
  geom_col(fill = "darkred") +
  labs(
    x = "NES",
    y = NULL,
    title = expression("Colon " * italic("(OSD-667)"))
  ) +
  theme_minimal(base_size = 7) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.3),
    axis.text.y  = element_text(size = 7),
    plot.margin  = fixed_margin,
    text         = element_text(size = 7, color = "black")
  )
```

Element heatmap
```
p_metal <- ggplot(
  metal_df,
  aes(x = Metal, y = Pathway, fill = -log10(padj_plot))
) +
  geom_tile() +
  scale_fill_gradient(
    low = "#A7C7E7",
    high = "#003366",
    na.value = "transparent",
    name = expression(-log[10]~FDR)
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 7) +
  theme(
    plot.margin      = fixed_margin,
    panel.background = element_blank(),
    panel.grid       = element_blank(),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.3),
    axis.text.y      = element_blank(),
    axis.ticks.y     = element_blank(),
    legend.key.width  = unit(3, "mm"),
    legend.key.height = unit(4, "mm"),
    text             = element_text(size = 7, color = "black")
  )
```

Combine & save
```
p <- plot_grid(
  p_nes + theme(plot.margin = margin(5, 1, 5, 25)),
  p_metal + theme(plot.margin = margin(5, 5, 5, 1)),
  nrow = 1,
  align = "h",
  axis = "tb",
  rel_widths = c(5.3, 2)
)

p <- ggdraw(p)

ggsave(
  out_file,
  p,
  width  = 5,
  height = 1.1,
  dpi    = 600,
  bg     = "transparent"
)

message("Saved: ", out_file)
```
