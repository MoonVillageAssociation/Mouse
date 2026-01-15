# MouseAnalysis

## KEGG

Setup
```
library(tidyverse)
library(clusterProfiler)

organism_kegg <- "mmu"
entrez_col <- "ENTREZID"
```

Define contrasts from RNA-seq data
```
contrast_map <- list(

  "GLDS-606" = list(
    list(
      stat_col = "Stat_(Wild Type & Space Flight)v(Wild Type & Ground Control)",
      flip     = FALSE,
      label    = "WT_SF_vs_GC"
    )
  ),

  "GLDS-673" = list(
    list(
      stat_col = "Stat_(Space Flight)v(Ground Control)",
      flip     = FALSE,
      label    = "SF_vs_GC"
    )
  ),

  "GLDS-674" = list(
    list(
      stat_col = "Stat_(Space Flight & 12 week & On ISS)v(Ground Control & 12 week & On Earth)",
      flip     = FALSE,
      label    = "SF12w_ISS_vs_GC12w_Earth"
    ),
    list(
      stat_col = "Stat_(Space Flight & 29 week & On ISS)v(Ground Control & 29 week & On Earth)",
      flip     = FALSE,
      label    = "SF29w_ISS_vs_GC29w_Earth"
    )
  )
)
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

Loop
```
for (glds in names(contrast_map)) {

  csv_file <- paste0(
    glds,
    "_rna_seq_differential_expression_GLbulkRNAseq.csv"
  )

  message("Processing ", glds)

  df <- read_csv(csv_file, show_col_types = FALSE)

  for (contrast in contrast_map[[glds]]) {

    message("  Contrast: ", contrast$label)

    if (!contrast$stat_col %in% colnames(df)) {
      warning("    Stat column not found: ", contrast$stat_col)
      next
    }

    gsea <- run_kegg_gsea(
      df       = df,
      stat_col = contrast$stat_col,
      flip     = contrast$flip
    )

    if (is.null(gsea) || nrow(gsea@result) == 0) {
      warning("    No significant KEGG pathways")
      next
    }

    res <- as.data.frame(gsea@result)
    res$GLDS     <- glds
    res$Contrast <- contrast$label
    res$StatUsed <- contrast$stat_col

    out_file <- paste0(
      glds, "_", contrast$label, "_KEGG_GSEA.csv"
    )

    write_csv(res, out_file)
    message("    Saved: ", out_file)
  }
}
```

