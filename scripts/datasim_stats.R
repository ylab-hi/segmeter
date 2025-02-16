#################### params ################
setwd("/Users/sej9799/Documents/projects/segmeter/simdata") # set wd to datadir
simname <- "sim_001"
intvlnums <- c("10", "100", "1K", "10K", "100K")
intvlmap <- c("10"=10, "100"=100, "1K"=1000, "10K"=10000, "100K"=100000)

library(ggplot2)
library(gtools)

plotsdir <- normalizePath(file.path(getwd(), "..", "plots/", "simdata"), mustWork = FALSE)
dir.create(plotsdir, showWarnings = FALSE, recursive = TRUE)

# takes a long time (only use once)
# stats = data.frame(
#   intvlnum = integer(),
#   intvl_size = double(),
#   gap_size = double()
# )
# 
# current_chr <- ""
# for (i in intvlnums) {
#   intvl <- intvlmap[i]
#   filename <- paste0(getwd(), "/sim/", simname, "/BED/ref/", i, "_sorted.bed")
#   # read in file and iterate through it
#   df <- read.table(filename, header=FALSE, stringsAsFactors=FALSE)
#   colnames(df) <- c("chr", "start", "end", "name")
#   
#   for (row in 1:nrow(df)) {
#     if(current_chr != df$chr[row]) {
#       intvlsize = df$end[row] - df$start[row]
#       gapsize = df$start[row]
#       newrow <- data.frame(intvlnum = intvl, 
#                            intvl_size = intvlsize,
#                            gap_size = gapsize)
#       stats <- rbind(stats, newrow)
#       current_chr = df$chr[row]
#       
#     } else {
#       intvlsize <- df$end[row] - df$start[row]
#       gapsize <- df$start[row] - df$end[row-1]
#       newrow <- data.frame(intvlnum = intvl, 
#                            intvl_size = intvlsize,
#                            gap_size = gapsize)
#       stats <- rbind(stats, newrow)
#       current_chr = df$chr[row]
#     }
#   }
# }

# create boxplot
p_box_gap <- ggplot(stats, aes(x = factor(intvlnum), y = gap_size, fill=factor(intvlnum))) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,
               outlier.size=4) +
  # geom_jitter(width = 0.2, alpha = 0.5, color = "blue") + 
  theme_minimal()
ggsave(file = paste0(plotsdir, "/", simname, "_gap_size_boxplot.pdf"), plot = p_box_gap)

p_box_intvl <- ggplot(stats, aes(x = factor(intvlnum), y = intvl_size, fill=factor(intvlnum))) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,
               outlier.size=4) +
  # geom_jitter(width = 0.2, alpha = 0.5, color = "blue") + 
  theme_minimal()
ggsave(file = paste0(plotsdir, "/", simname, "_intvl_size_boxplot.pdf"), plot = p_box_intvl)

# stats about chromosomes/intervals (within the chrnums, chrlens)
for (i in intvlnums) {
  filename_chrnums <- paste0(getwd(), "/sim/", simname, "/BED/", i, "_chrnums.txt")
  chrnums <- read.csv(filename_chrnums, header=FALSE, stringsAsFactors=FALSE, sep="\t")
  colnames(chrnums) <- c("chr", "num")
  
  chrnums$chr <- factor(chrnums$chr, levels = mixedsort(unique(chrnums$chr)))
  p_bar_chrnums <- ggplot(chrnums, aes(x = chr, y = num)) +
    geom_bar(stat="identity") +
    theme_minimal()
  ggsave(file = paste0(plotsdir, "/", simname, "_", i, "_chrnums_barplot.pdf"), plot = p_bar_chrnums)
  
  filename_chrlens <- paste0(getwd(), "/sim/", simname, "/BED/", i, "_chromlens.txt")
  chrlens <- read.csv(filename_chrlens, header=FALSE, stringsAsFactors=FALSE, sep="\t")
  colnames(chrlens) <- c("chr", "len")
  
  chrlens$chr <- factor(chrlens$chr, levels = mixedsort(unique(chrlens$chr)))
  p_hist_chrlens <- ggplot(chrlens, aes(x=chr, y=len)) + 
    geom_bar(stat="identity") +
    theme_minimal()
  ggsave(file = paste0(plotsdir, "/", simname, "_", i, "_chrlens_barplot.pdf"), plot = p_hist_chrlens)
}

