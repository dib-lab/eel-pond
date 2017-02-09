install.packages.auto <- function(x) {
  x <- as.character(substitute(x))
  if(isTRUE(x %in% .packages(all.available=TRUE))) {
    eval(parse(text = sprintf("require(\"%s\")", x)))
  } else {
    #update.packages(ask= FALSE) #update installed packages.
    eval(parse(text = sprintf("install.packages(\"%s\", dependencies = TRUE,repos=structure(c(CRAN='https://cran.cnr.berkeley.edu/')))", x)))
  }
  if(isTRUE(x %in% .packages(all.available=TRUE))) {
    eval(parse(text = sprintf("require(\"%s\")", x)))
  } else {
    source("http://bioconductor.org/biocLite.R")
    #biocLite(character(), ask=FALSE) #update installed packages.
    eval(parse(text = sprintf("biocLite(\"%s\")", x)))
    eval(parse(text = sprintf("require(\"%s\")", x)))
  }
}
install.packages.auto("edgeR")

files <- c("0Hour_ATCACG_L002_R1_001.extract.quant.counts",
           "0Hour_ATCACG_L002_R1_002.extract.quant.counts",
           "0Hour_ATCACG_L002_R1_003.extract.quant.counts",
           "0Hour_ATCACG_L002_R1_004.extract.quant.counts",
           "0Hour_ATCACG_L002_R1_005.extract.quant.counts",
           "6Hour_CGATGT_L002_R1_001.extract.quant.counts",
           "6Hour_CGATGT_L002_R1_002.extract.quant.counts",
           "6Hour_CGATGT_L002_R1_003.extract.quant.counts",
           "6Hour_CGATGT_L002_R1_004.extract.quant.counts",
           "6Hour_CGATGT_L002_R1_005.extract.quant.counts"
)


labels=c("0Hour_1", "0Hour_2", "0Hour_3", "0Hour_4", "0Hour_5",
"6Hour_1", "6Hour_2", "6Hour_3",
"6Hour_4", "6Hour_5")

data <- readDGE(files)

#print(data)
#head(data$counts)

###

group <- c(rep("0Hour",5), rep("12Hour",4), rep("18Hour", 8),
          rep("24Hour", 11), rep("6Hour", 5))

dge = DGEList(counts=data, group=group)
dge <- estimateCommonDisp(dge)
dge <- estimateTagwiseDisp(dge)

# make an MA-plot of 0 vs 6 hour

et <- exactTest(dge, pair=c("0Hour", "6Hour"))
etp <- topTags(et, n=100000)
etp$table$logFC = -etp$table$logFC
pdf("nema-edgeR-MA-plot.pdf")
plot(
  etp$table$logCPM,
  etp$table$logFC,
  xlim=c(-3, 20), ylim=c(-12, 12), pch=20, cex=.3,
  col = ifelse( etp$table$FDR < .2, "red", "black" ) )
dev.off()

# plot MDS
pdf("nema-edgeR-MDS.pdf")
plotMDS(dge, labels=labels)
dev.off()

# output CSV for 0-6 hr
write.csv(etp$table, "nema-edgeR-0v6.csv")

