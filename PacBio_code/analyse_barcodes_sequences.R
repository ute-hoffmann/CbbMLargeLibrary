library(ggplot2)
barcodes <- read_csv("output/selected_barcodes_count.csv", col_names = c("mut", "occurrences"))
barcodes <- barcodes[2:nrow(barcodes), ]
p <- ggplot(barcodes, aes(x=occurrences)) + geom_histogram(fill="#69b3a2") + theme_light()
p
ggsave("final_set_barcode_hist.pdf", plot=p, width=15, height=9, units="cm")

length(unique(barcodes$mut))
nrow(subset(barcodes, barcodes$occurrences > 1))