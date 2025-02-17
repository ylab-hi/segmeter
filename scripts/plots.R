setwd("/Users/sej9799/Documents/projects/segmeter/simdata/")
simname <- "sim_001" # name of simulation data to use

library(ggplot2)
library(dplyr)


datapath <- "simdata/" # remember to add / at the end

perfect <- list()
data_basic <- list()

intvlnums <- c("10","100","1K","10K","100K")
intvlmap <- {}
intvlmap["10"] <- 10
intvlmap["100"] <- 100
intvlmap["1K"] <- 1000
intvlmap["10K"] <- 10000
intvlmap["100K"] <- 100000

# tools <- c("tabix", "bedtools", "giggle", "bedtk", "ucsc", "gia", "granges", "bedops", "igd", "ailist")
tools <- c("bedops", "bedmaps")

# queries
intvlqueries <- c("perfect", "5p-partial", "3p-partial", "enclosed", "contained")
gapqueries <- c("perfect-gap", "left-adjacent-gap", "right-adjacent-gap", "mid-gap1", "mid-gap2")
cmplxqueries <- c("mult")

# percent <- paste0(seq(10,100, by=10), '%')
percent <- seq(10,100,10)
replicates <- c("bench_001", "bench_002", "bench_003")

data_complete <- NULL
for (i in intvlnums) {
  for (j in tools) {
    for (k in percent) { # 
      repldata <- NULL
      for (l in 1:length(replicates)) {
        filename <- paste0(getwd(), "/bench/", replicates[[l]], "/", j, "/", i, "/stats/", i, "_query_stats_", k, ".txt")
        tmpdata <- read.csv(filename, header=TRUE,sep="\t")
        # rename time and memory columns (according to the replicates)
        colnames(tmpdata)[c(4,5)] <- c(paste0("time_", replicates[[l]]), paste0("mem_", replicates[[l]]))
        tmpdata <- cbind(toolname = rep(j, nrow(tmpdata)), tmpdata) # add column for toolname
        # add subset (insert after data_type column)
        tmpdata <- cbind(tmpdata[,1:3], subset = rep(k, nrow(tmpdata)), tmpdata[,4:6]) # rearrange columns
        # rename query type column
        tmpdata$query_type <- c(intvlqueries, gapqueries, cmplxqueries)
        
        basic <- tmpdata[tmpdata$data_type == "basic", c(6,7)]
        complex <- tmpdata[tmpdata$data_type == "complex", c(6,7)]
        
        total_rows <- data.frame(
          toolname = rep(j,2),
          intvlnum = rep(intvlmap[[i]],2),
          data_type = c("basic","complex"),
          subset = rep(k,2),
          query_type = rep("all",2),
          time = c(sum(basic[,1]), sum(complex[,1])),
          mem = c(max(basic[,2]), max(complex[,2]))
        )
        colnames(total_rows)[c(6,7)] <- c(paste0("time_", replicates[[l]]), paste0("mem_", replicates[[l]]))
        
        # add to tmpdata
        tmpdata <- rbind(tmpdata, total_rows)
        
        if(is.null(repldata)) {
          repldata <- tmpdata
        } else {
          repldata <- cbind(repldata, tmpdata[,c(6,7)])
        }
      }
      
      # correct memory for some tools (e.g., bedmaps, bedops, bedtools, bedtools_sorted, bedtools_tabix, giggle, igd)
      # *1000000/1024
      if (j == "bedmaps" || 
          j == "bedops" ||
          j == "bedtools" ||
          j == "bedtools_sorted" ||
          j == "bedtools_tabix" ||
          j == "giggle" ||
          j == "gia" ||
          j == "igd") {
        repldata[,7:ncol(repldata)] <- (repldata[,7:ncol(repldata)]*1000000)/1024
      }
      
      
      # calculate stats
      idx_time <- seq(6, 6+length(replicates)*2-1, by=2)
      repldata$mean_time <- rowMeans(repldata[,idx_time])
      repldata$sd_time <- apply(repldata[,idx_time], 1, sd)
      repldata$se_time <- repldata$sd_time/sqrt(length(replicates))
      repldata$ci_time <- repldata$se_time*1.96
      
      # same for memory
      idx_mem <- seq(7, 7+length(replicates)*2-1, by=2)
      repldata$mean_mem <- rowMeans(repldata[,idx_mem])
      repldata$sd_mem <- apply(repldata[,idx_mem], 1, sd)
      repldata$se_mem <- repldata$sd_mem/sqrt(length(replicates))
      repldata$ci_mem <- repldata$se_mem*1.96
      
      
      if (is.null(data_complete)) {
        data_complete <- repldata
      } else {
        data_complete <- rbind(data_complete, repldata)
      }
    }
  }
}


###### PLOTS #########
plotsdir_runtime <- normalizePath(file.path(getwd(), "..", "plots/", "runtime"), mustWork = FALSE)
dir.create(plotsdir_runtime, showWarnings = FALSE, recursive = TRUE)

plotsdir_memory <- normalizePath(file.path(getwd(), "..", "plots/", "memory"), mustWork = FALSE)
dir.create(plotsdir_memory, showWarnings = FALSE, recursive = TRUE)

############# BASIC ####################

# basic runtime (plot for each query type) - 100%
pd <- position_dodge(0.1)
for (i in c(intvlqueries, gapqueries)) {
  # get subset
  tmpdata <- data_complete[data_complete$query_type == i & data_complete$subset == 100,]
  p <- ggplot(tmpdata, aes(x = intvlnum, y = mean_time, color = toolname, group = toolname)) + 
    geom_errorbar(aes(ymin=mean_time-se_time, ymax=mean_time+se_time), width=.1, position=pd) +
    geom_line(size = .75, position = pd) + geom_point(size = 7, position = pd, aes(shape = toolname)) + 
    scale_x_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    scale_y_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    scale_x_log10() + # Use log scale for better visualization of interval numbers
    scale_y_log10() + # Use log scale for better visualization of time
    scale_shape_manual(values = c(0:14)) +
    labs(
      title = paste0(i, " queries ", "100%"),
      x = "interval Number",
      y = "time in seconds (log scale)",
      color = "Tool"
    ) +
    theme_minimal()
  
  ggsave(paste0(plotsdir_runtime, "/", i, "_time_queries_100.pdf"), plot = p)
}

# do the same but for total time and memory
all_runtime <- data_complete[data_complete$query_type == "all" & data_complete$subset == 100 & data_complete$data_type == "basic",]
p_all <- ggplot(all_runtime, aes(x = intvlnum, y = mean_time, color = toolname, group = toolname)) + 
  geom_errorbar(aes(ymin=mean_time-se_time, ymax=mean_time+se_time), width=.1, position=pd) +
  geom_line(size = .75, position = pd) + geom_point(size = 7, position = pd, aes(shape = toolname)) + 
    scale_x_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    scale_y_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    scale_x_log10() + # Use log scale for better visualization of interval numbers
    scale_y_log10() + # Use log scale for better visualization of time
    scale_shape_manual(values = c(0:14)) +
    labs(
      title = paste0("All time queries ", "100%"),
      x = "interval Number",
      y = "time in seconds (log scale)",
      color = "Tool"
    ) +
    theme_minimal()
ggsave(paste0(plotsdir_runtime, "/all_time_queries_100.pdf"), plot = p_all)

# create similar plot for runtime (barplot) - we consider only 100,000 intervals (since its the maximum anyway)
mostintvls <- intvlmap[[intvlnums[[length(intvlnums)]]]] # get highest number of interval (since thats the maximum mem)
basic_all_highintvl <- data_complete[data_complete$intvlnum == mostintvls &
                             data_complete$query_type == "all" & 
                             data_complete$data_type == "basic" &
                             data_complete$subset == 100,]

basic_all_highintvl <- ggplot(basic_all_highintvl, aes(x = toolname, y = mean_mem, fill = toolname)) + 
  geom_bar(stat="identity", position=pd) + 
  geom_errorbar(aes(ymin=mean_mem-se_mem, ymax=mean_mem+se_mem), width=.1, position=pd)
ggsave(paste0(plotsdir_memory, "/memory_basic_allintervals.pdf"), plot = basic_all_highintvl)

# p_mem_basic_all <- ggplot(all_runtime, aes(x = intvlnum, y = mean_mem, color = toolname, group = toolname)) + 
#   # geom_errorbar(aes(ymin=mean_mem-se_mem, ymax=mean_mem+se_mem), width=.1, position=pd) +
#   geom_line(size = .75, position = pd) + geom_point(size = 7, position = pd, aes(shape = toolname)) + 
#     # scale_x_log10(
#     #   breaks = scales::trans_breaks("log10", function(x) 10^x),
#     #   labels = scales::trans_format("log10", scales::math_format(10^.x))
#     # ) +
#     scale_y_log10(
#       breaks = scales::trans_breaks("log10", function(x) 10^x),
#       labels = scales::trans_format("log10", scales::math_format(10^.x))
#     ) +
#     # scale_x_log10() + # Use log scale for better visualization of interval numbers
#     scale_y_log10() + # Use log scale for better visualization of time
#     scale_shape_manual(values = c(0:14)) +
#     labs(
#       title = paste0(i, " queries ", "100%"),
#       x = "interval Number",
#       y = "memory in MB",
#       color = "Tool"
#     ) +
#     theme_minimal()
  
  


# stacked bar plot for basic queries (to see the fraction for each type to the total runtime)
# this is average over all intvlnums
stackeddata <- data_complete[data_complete$query_type %in% c(intvlqueries,gapqueries),]
stacked_order <- c(intvlqueries, gapqueries)

query_fraction_data <- stackeddata %>%
  group_by(toolname, query_type) %>%  # group by both variables at once
  summarise(mean_time_per_query = mean(mean_time), .groups = "drop") %>%  # calculate means
  group_by(toolname) %>%  # regroup by tool for proportion calculation
  mutate(total = sum(mean_time_per_query),
         proportion = mean_time_per_query/total * 100) %>%
  ungroup()

ps <- ggplot(query_fraction_data, aes(x=toolname, y=proportion, 
                                      fill=factor(query_type, levels=stacked_order))) +  # use factor to set order
  geom_bar(stat="identity", position="fill") +
  scale_y_continuous(labels=scales::percent_format(scale=1)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position="right"
  ) +
  labs(x="Tool", 
       y="Percentage", 
       fill="Query Type")
ggsave(paste0(plotsdir_runtime, "/basic_query_fraction.pdf"), plot = ps)









##### complex runtime
# determine number of interval in each subset (simdata)
cmplx_subset_size = data.frame(intvlnums=rep(intvlnums, each=length(percent)),
                               percent=rep(percent, times=length(intvlnums)),
                               numqueries=rep(0, length(intvlnums)*length(percent)))
# for each intvlnum (determine number of complex queries)
for (i in intvlnums) {
  for (j in percent) {
    cmplx_subset_multfile <- paste0(getwd(), "/sim/", simname, "/BED/complex/query/mult/", i, "_", j, "bin.bed")
    # checks if subset file is empty (can happen for small number of intervals and high fraction)  
    if(file.size(cmplx_subset_multfile) != 0) {
      cmplx_subset_size[cmplx_subset_size$intvlnums == i & cmplx_subset_size$percent == j,3] <- nrow(read.table(cmplx_subset_multfile))
    }
  }
}


## plot per intvlnums
plotsdir_runtime_cmplx <- paste0(plotsdir_runtime, "/cmplx")
dir.create(plotsdir_runtime_cmplx, showWarnings = FALSE, recursive = TRUE)

for (i in intvlnums) {
  
  cmplxdata <- data_complete[data_complete$intvlnum == intvlmap[[i]] & 
                                 data_complete$data_type == "complex" &
                                 data_complete$query_type == "mult",]
  
  # divide the mean time by the number of queries (according to intvlnum and subset)
  for (j in percent) {
    subset_size <- cmplx_subset_size[cmplx_subset_size$intvlnums == i & cmplx_subset_size$percent == j,3]
    if(subset_size != 0) {
      cmplxdata[cmplxdata$intvlnum == intvlmap[[i]] & cmplxdata$subset == j,]$mean_time <- cmplxdata[cmplxdata$intvlnum == intvlmap[[i]] & cmplxdata$subset == j,]$mean_time/subset_size
      cmplxdata[cmplxdata$intvlnum == intvlmap[[i]] & cmplxdata$subset == j,]$se_time <- cmplxdata[cmplxdata$intvlnum == intvlmap[[i]] & cmplxdata$subset == j,]$se_time/subset_size
      cmplxdata[cmplxdata$intvlnum == intvlmap[[i]] & cmplxdata$subset == j,]$sd_time <- cmplxdata[cmplxdata$intvlnum == intvlmap[[i]] & cmplxdata$subset == j,]$sd_time/subset_size
      cmplxdata[cmplxdata$intvlnum == intvlmap[[i]] & cmplxdata$subset == j,]$ci_time <- cmplxdata[cmplxdata$intvlnum == intvlmap[[i]] & cmplxdata$subset == j,]$ci_time/subset_size
    }
  }
  
  # plot
  p_cmplxdata_intvl <- ggplot(cmplxdata, aes(x=subset, y=mean_time,  color = toolname, group = toolname)) + 
    geom_line(size = .75) + geom_point(size = 7, aes(shape = toolname)) + 
    geom_errorbar(aes(ymin=mean_time-se_time, ymax=mean_time+se_time), width=.1, position=pd) +
    scale_y_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    scale_y_log10() + # Use log scale for better visualization of time
    scale_y_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    scale_shape_manual(values = c(0:14)) +
    labs(
      title =  paste0("Complex Queries (", i, ")"),
      x = "Complex Queries Subset",
      y = "time in seconds (log scale)",
      color = "Tool"
    ) +
    theme_minimal()
  
  ggsave(paste0(plotsdir_runtime_cmplx, "/", i, "_multqueries_time.pdf"), plot=p_cmplxdata_intvl, width=300, height=300, unit="mm")
}







# plot all
tmpdata <- data_complete[data_complete$data_type == "complex" & data_complete$query_type == "all",]
# merge data for different subset (in new data.frame)

intvlnums_numeric <- unname(mapply(function(x) intvlmap[[x]], intvlnums))
all_complex <- data.frame(toolname = rep(tools, each=length(intvlnums)),
                                         intvlnums = rep(intvlnums_numeric, times=length(tools)))
for (r in replicates) {
  zero_col <- rep(0.0, length(tools)*length(intvlnums_numeric))
  time_col <- paste0("time_", r)
  mem_col <- paste0("mem_", r)
  all_complex <- cbind(all_complex, setNames(list(zero_col), time_col))
  all_complex <- cbind(all_complex, setNames(list(zero_col), mem_col))
}

for (i in tools) {
  for (j in intvlnums) {
    entry <- data_complete[data_complete$data_type == "complex" & 
                             data_complete$query_type == "all" & 
                             data_complete$toolname == i & 
                             data_complete$intvlnum == intvlmap[[j]],]
    
    for (r in replicates) {
      timecolumn <- paste0("time_", r)
      memcolumn <- paste0("mem_", r)
      coltime_sum <- sum(entry[,timecolumn])
      colmax_mem <- max(entry[,memcolumn])
      
      all_complex[all_complex$toolname == i & all_complex$intvlnums == intvlmap[[j]], timecolumn] <- coltime_sum
      all_complex[all_complex$toolname == i & all_complex$intvlnums == intvlmap[[j]], memcolumn] <- colmax_mem
      
    }
  }
}

# calculate stats
idx_time <- seq(3, 3+length(replicates)*2-1, by=2)
all_complex$mean_time <- rowMeans(all_complex[,idx_time])
all_complex$sd_time <- apply(all_complex[,idx_time], 1, sd)
all_complex$se_time <- all_complex$sd_time/sqrt(length(replicates))
all_complex$ci_time <- all_complex$se_time*1.96

idx_mem <- seq(4, 4+length(replicates)*2-1, by=2)
all_complex$mean_mem <- rowMeans(all_complex[,idx_mem])
all_complex$sd_mem <- apply(all_complex[,idx_mem], 1, sd)
all_complex$se_mem <- all_complex$sd_mem/sqrt(length(replicates))
all_complex$ci_mem <- all_complex$se_mem*1.96

# filter tools (only keep the ones with perfect precision)
# tools_to_keep <- c("tabix", "bedops", "gia", "igd", "giggle", "bedtk", "ailist", "bedtools", "ucsc", "granges")
tools_to_keep <- c("bedops", "bedmaps")
all_complex <- all_complex[all_complex$toolname %in% tools_to_keep,]
                
p_all_complex <- ggplot(all_complex, aes(x = intvlnums, y = mean_time, color = toolname, group = toolname)) + 
  geom_line(size = .75) + geom_point(size = 7, aes(shape = toolname)) + 
  geom_errorbar(aes(ymin=mean_time-se_time, ymax=mean_time+se_time), width=.1, position=pd) +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  # scale_x_log10() + # Use log scale for better visualization of interval numbers
  # scale_y_log10() + # Use log scale for better visualization of time
  scale_shape_manual(values = c(0:14)) +
  labs(
    title =  "Complex Queries (All)",
    x = "interval Number",
    y = "time in seconds (log scale)",
    color = "Tool"
  ) +
  theme_minimal()

ggsave(paste0(plotsdir_runtime, "/complex_all.pdf"), plot=p_all_complex)


# create similar plot for runtime (barplot) - we consider only 100,000 intervals (since its the maximum anyway)
pd <- position_dodge(width = 0.1)  # Ensure consistent positioning

p_mem_complex <- ggplot(all_complex, aes(x = toolname, y = mean_mem, fill = toolname, group=toolname)) + 
  geom_bar(stat = "identity", position = pd) +
  geom_errorbar(aes(ymin = mean_mem - se_mem, ymax = mean_mem + se_mem), 
                width = .1, position = pd) +  # Adjust width for better visibility
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x), 
                labels = scales::trans_format("log10", scales::math_format(10^.x)))
ggsave(paste0(plotsdir_memory, "/memory_complex_allintervals.pdf"), plot = p_mem_complex)



#### INDEX ######
tools_idx <- c("tabix", "bedtools_sorted", "giggle", "igd")

df_idx <- data.frame()
  


for (i in tools_idx) {
  for (j in intvlnums) {
    df_entry <- data.frame(toolname = i, intvlnum = j)
    for (k in replicates) {
      idx_stats_filename <- paste0(getwd(), "/bench/", k, "/", i, "/", j, "/", j, "_idx_stats.txt")
      idx_stats <- read.csv(idx_stats_filename, header=TRUE,sep="\t")
      colnames(idx_stats) <- c("intvlnum", "time", "mem", "size")
      
      df_entry[paste0("time_", k)] <- idx_stats$time
      
      if (i == "giggle") {
        df_entry[paste0("mem_", k)] <- (idx_stats$mem*1000000)/1024
      } else {
        df_entry[paste0("mem_", k)] <- idx_stats$mem
      }
      
      
      
      # df_entry[paste0("size_", k)] <- idx_stats$size
    }
    df_idx <- rbind(df_idx, df_entry)
  }
}

# calculate stats
idx_time <- seq(3, 3+length(replicates)*2-1, by=2)
df_idx$mean_time <- rowMeans(df_idx[,idx_time]) 
df_idx$sd_time <- apply(df_idx[,idx_time], 1, sd)
df_idx$se_time <- df_idx$sd_time/sqrt(length(replicates))
df_idx$ci_time <- df_idx$se_time*1.96

idx_mem <- seq(4, 4+length(replicates)*2-1, by=2)
df_idx$mean_mem <- rowMeans(df_idx[,idx_mem])
df_idx$sd_mem <- apply(df_idx[,idx_mem], 1, sd)
df_idx$se_mem <- df_idx$sd_mem/sqrt(length(replicates))
df_idx$ci_mem <- df_idx$se_mem*1.96

plotsdir_index <- normalizePath(file.path(getwd(), "..", "plots/", "index"), mustWork = FALSE)
dir.create(plotsdir_index, showWarnings = FALSE, recursive = TRUE)

#barplot index time
plot_idx_time <- ggplot(df_idx, aes(x = toolname, y = mean_time, fill = toolname)) + 
  geom_bar(stat="identity", position = pd) + 
  geom_errorbar(aes(ymin=mean_time-se_time, ymax=mean_time+se_time), width=.1)
ggsave(paste0(plotsdir_index, "/index_time.pdf"), plot = plot_idx_time)


plot_idx_mem <- ggplot(df_idx, aes(x = toolname, y = mean_mem, fill = toolname)) + 
  geom_bar(stat="identity", position = pd) + 
  geom_errorbar(aes(ymin=mean_mem-se_time, ymax=mean_mem+se_time), width=.1)
ggsave(paste0(plotsdir_index, "/index_mem.pdf"), plot = plot_idx_mem)




   













  
  
        
        
        
#         
