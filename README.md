# MouseAnalysis

## KEGG

Create contrast_map.csv to define data and contrasts
```
ID,GLDS,stat_col,flip,label
1,GLDS-606,Stat_(Wild Type & Space Flight)v(Wild Type & Ground Control),FALSE,WT_SF_vs_GC
2,GLDS-673,Stat_(Space Flight)v(Ground Control),FALSE,SF_vs_GC
3,GLDS-674,Stat_(Space Flight & 12 week & On ISS)v(Ground Control & 12 week & On Earth),FALSE,SF12w_ISS_vs_GC12w_Earth
4,GLDS-674,Stat_(Space Flight & 29 week & On ISS)v(Ground Control & 29 week & On Earth),FALSE,SF29w_ISS_vs_GC29w_Earth
```

Specify data to analyze from contrast_map.csv
```
TARGET_ID <- 4
```

Setup
```
library(tidyverse)
library(clusterProfiler)

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

## KEGG
...
