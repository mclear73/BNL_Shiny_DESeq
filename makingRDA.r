library(DESeq2)
library(readr)

data <- readr::read_csv("Chlam_Counts.csv", col_types = readr::cols())


createDF <- function(filesz) {
  data <- readr::read_csv(filesz, col_types = readr::cols())
#  gene.ID <- data %>% pull(1)
#  raw.counts <- data[,-1]
#  rownames(raw.counts) <- gene.ID
  return(data)
}


ret$vst <- createDF("vst.csv")
ret$rownorm <- createDF("rownorm.csv")
ret$raw <- createDF("raw.csv")
ret$cpm <- createDF("cpm.csv")
ret$rlog <- createDF("rlog.csv")

save(ret, file="HomRNA.rda")

load("genavi.rda")

edit(ret$vst)
