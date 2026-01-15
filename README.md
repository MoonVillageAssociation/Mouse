# MouseAnalysis

## KEGG

'''
library(tidyverse)
library(clusterProfiler)

organism_kegg <- "mmu"
entrez_col <- "ENTREZID"
'''

Define contrasts from RNA-seq data
'''
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
'''

