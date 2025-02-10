setwd("/Users/sej9799/Documents/projects/segmeter/simdata/")

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

tools <- c("tabix", "bedtools", "giggle", "bedtk", "ucsc", "gia", "granges", "bedops", "igd", "ailist")

# queries
intvlqueries <- c("perfect", "5p-partial", "3p-partial", "enclosed", "contained")
gapqueries <- c("perfect-gap", "left-adjacent-gap", "right-adjacent-gap", "mid-gap1", "mid-gap2")
cmplxqueries <- c("mult")

# percent <- paste0(seq(10,100, by=10), '%')
percent <- seq(10,100,10)
replicates <- c("bench_001", "bench_002", "bench_003")


data_complete <- NULL
for (i in intvlnums) {
  data_intvl[[i]] <- list()
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
        
        # print(basic)
        # print(sum(basic[,1]))
        # 
        # 
        
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
        
        # print(tmpdata)
        # print(total_rows)
        
        # add to tmpdata
        tmpdata <- rbind(tmpdata, total_rows)
        
        if(is.null(repldata)) {
          repldata <- tmpdata
        } else {
          repldata <- cbind(repldata, tmpdata[,c(6,7)])
        }
      }
      
      # calculate stats
      idx_time <- seq(6, 6+length(replicates)*2-1, by=2)
      repldata$mean_time <- rowMeans(repldata[,idx_time])
      repldata$sd_time <- apply(repldata[,idx_time], 1, sd)
      repldata$se_time <- repldata$sd_time/sqrt(length(replicates))
      repldata$ci_time <- repldata$se_time*1.96
      
      if (is.null(data_complete)) {
        data_complete <- repldata
      } else {
        data_complete <- rbind(data_complete, repldata)
      }
    }
  }
}


###### PLOTS #########
plotsdir <- normalizePath(file.path(getwd(), "..", "plots/", "runtime"), mustWork = FALSE)
dir.create(plotsdir, showWarnings = FALSE, recursive = TRUE)

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
  
  ggsave(paste0(plotsdir, "/", i, "_time_queries_100.pdf"), plot = p)
}

# do the same but for total time
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
      title = paste0(i, " queries ", "100%"),
      x = "interval Number",
      y = "time in seconds (log scale)",
      color = "Tool"
    ) +
    theme_minimal()
ggsave(paste0(plotsdir, "/all_time_queries_100.pdf"), plot = p_all)




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
ggsave(paste0(plotsdir, "/basic_query_fraction.pdf"), plot = ps)

### complex runtime
# plot all
tmpdata <- data_complete[data_complete$data_type == "complex" & data_complete$query_type == "all",]
# merge data for different subset (in new data.frame)

for (i in tools) {
  for (j in intvlnums) {
    tmpdata <- data_complete[data_complete$data_type == "complex" & data_complete$query_type == "all" & data_complete$toolname == i & data_complete$intvlnum == intvlmap[[j]],]
  }
}







for (i in intvlnums) {
  
}






  
  
        
        
        
#         
#         filename <- paste0(datapath, "bench/",k,"/",i,"/stats/",i,"_query_stats_",l,".txt")
#         tmpdata <- read.csv(filename, header=TRUE,sep="\t")
#         
#         toolname <- rep(j, nrow(tmpdata)) # toolname column
#         tmpdata <- cbind(toolname, tmpdata)
#         
#         data_intvl[[i]] <- rbind(data_intvl[[i]], tmpdata)
#       }
#       
#       
#       
#       filename <- paste0(datapath, "bench/",k,"/",i,"/stats/",i,"_query_stats_")
#       
#       
#       
#       filename <- paste0(getwd(),"/bench/",j,"/",i,"_query_stats_100.txt")
#       tmpdata <- read.csv(paste(getwd(),"/bench/",j,"/",i,"_query_stats_100.txt",sep=""), header=TRUE,sep="\t")
#       
#       toolname <- rep(j, nrow(tmpdata)) # toolname column
#       tmpdata <- cbind(toolname, tmpdata)
#       print(tmpdata) 
#       
#       all[[k]] <- tmpdata
#     }
#     
#     
#     
#     
#     
#     filename <- paste0(getwd(),"/bench/",j,"/",i,"_query_stats_100.txt")
#     tmpdata <- read.csv(paste(getwd(),"/bench/",j,"/",i,"_query_stats_100.txt",sep=""), header=TRUE,sep="\t")
#     
#     
#     toolname <- rep(j, nrow(tmpdata)) # toolname column
#     tmpdata <- cbind(toolname, tmpdata)
#     print(tmpdata) 
#     
#     intvl <- tmpdata[tmpdata$query_type %in% paste0(intvlqueries,"_100%"), ]
#     gap <- tmpdata[tmpdata$query_type %in% paste0(gapqueries,"_100%"), ]
#    
#     print(intvl)
# 
#     data_intvl[[i]] <- rbind(data_intvl[[i]], intvl)
#     data_intvl[[i]] <- rbind(data_intvl[[i]], gap)
#   }
# }
# 
# data_queries <- list()
# data_queries[["total"]] <- data.frame(toolname = rep(tools, each=length(intvlnums)),
#                     intvlnum = rep(intvlnums_numeric, times=length(tools)),
#                     data_type = rep("both", times=length(tools)*length(intvlnums)),
#                     query_type = rep("total", times=length(tools)*length(intvlnums)),
#                     time = rep(0, times=length(tools)*length(intvlnums)),
#                     max_RSS.MB. = rep(0, times=length(tools)*length(intvlnums)))
# 
# for (i in c(intvlqueries,gapqueries)) {
#   data_queries[[i]] <- list()
#   for (j in intvlnums) {
#     tmpdata <- data_intvl[[j]]
#     qtype <- tmpdata[tmpdata$query_type %in% paste0(i,"_100%"), ]
#     data_queries[[i]] <- rbind(data_queries[[i]], qtype)
# 
#     for (k in 1:nrow(qtype)) {
#       total <- data_queries[["total"]]
#       qtype_toolname <- qtype[k,]$toolname
#       qtype_intvlnum <- qtype[k,]$intvlnum
# 
#       time <- qtype[k,]$time
#       mem <- qtype[k,]$max_RSS.MB.
#       tmp <- total[total$intvlnum == qtype_intvlnum & total$toolname == qtype_toolname,]$time
#       tmp <- tmp + time
#       total[total$intvlnum == qtype_intvlnum & total$toolname == qtype_toolname,]$time <- tmp
#       data_queries[["total"]] <- total
#     }
#   }
# }
# # 
# # 
# # create folder plots
# dir.create("plots", showWarnings = FALSE)
# 
# ##### PLOTS #######
# for (i in 1:length(data_queries)) {
#   p <- ggplot(data_queries[[i]], aes(x = intvlnum, y = time, color = toolname)) +
#     geom_line(size = 1) +
#     geom_point(size = 2) +
#     scale_x_log10(
#       breaks = scales::trans_breaks("log10", function(x) 10^x),
#       labels = scales::trans_format("log10", scales::math_format(10^.x))
#     ) +
#     scale_y_log10(
#       breaks = scales::trans_breaks("log10", function(x) 10^x),
#       labels = scales::trans_format("log10", scales::math_format(10^.x))
#     ) +
#     scale_x_log10() + # Use log scale for better visualization of interval numbers
#     scale_y_log10() + # Use log scale for better visualization of time
#     labs(
#       title = paste0(names(data_queries)[[i]], " queries"),
#       x = "Interval Number (log scale)",
#       y = "Time (log seconds)",
#       color = "Tool"
#     ) +
#     theme_minimal()
# 
#   ggsave(paste0(getwd(),"/plots/", names(data_queries)[i], "time_queries", ".pdf"), plot = p)
# }
# 
# 
# 
# 
# # ggplot(data_intvl$`100K`, aes(x=toolname, y=time, fill=query_type)) + geom_bar(stat="identity")
