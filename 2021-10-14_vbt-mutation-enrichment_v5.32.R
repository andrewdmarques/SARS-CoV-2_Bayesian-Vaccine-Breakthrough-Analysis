# Load libraries.
library(readxl)
library(lubridate)
library(dplyr)
library(plyr)
library(ggplot2)
library(reshape)
library(plotly)
library(bayesplot)

# Function to graph mutations of interest#######################################

# Define the function that Smooths the data out.

curate_moi <- function(file_mutation, file_genome, time_bin, gene, title_graph, title_legend, rationale, moi){
  
  # Count how many genomes were sequenced for the time bin
  b_original <- read_excel(file_genome)
  b1 <- b_original[c("sample_date","lineage","rationale","lab_id", "trial")]
  
  # Use only the random samples from the master file
  #rationale <- "vaccine breakthrough"
  b2 <- b1
  if(rationale == "surveillance"){
    pattern <- c("surveillance", "hospitalized","asymptomatic")
    b2 <- dplyr::filter(b2, grepl(paste(pattern, collapse="|"),rationale))
    b2 <- subset(b2, !grepl("breakthrough", b2$rationale))
  }
  if(rationale == "vaccine breakthrough"){
    b2 <- dplyr::filter(b2, grepl(c("vaccine breakthrough"),rationale))
    b2 <- subset(b2, rationale != "vaccine breakthrough (<14 days)")
  }
  if(rationale == "s drop"){
    b2 <- dplyr::filter(b2, grepl("s drop",rationale))
  }
  if(rationale == "other"){
    pattern <- c("outbreak",NA ,"university","special request","vaccine breakthrough (<14 days)")
    b2 <- dplyr::filter(b2, grepl(paste(pattern, collapse="|"),  rationale))
  }
  
  # Count the number of genomes for each bin.
  b3 <- b2
  b3$sample_date <- as.Date(b3$sample_date, "%Y-%m-%d")
  b3$sample_date_bin <- b3$sample_date
  # Summarize the data
  b3$counts <- 'y'
  b3$sample_date_bin <- floor_date(b3$sample_date, time_bin) # Adds a column with date for the week only
  b4 <- ddply(b3, c("sample_date_bin", "counts"), summarise, n = length(counts)) # Makes a new dataframe with the counts of each lineage for each of the week categories
  
  # Count the mutations for the graphing data frame.
  # Open Excel file for mutations. 
  a_original <- readxl::read_excel(file_mutation)
  
  # Select only mutations in spike.
  a1 <- dplyr::filter(a_original, genes == gene)
  
  ### Correct the deletions for the mutations.
  for(i in 1:length(a1$VSP)){
    if(grepl("del", a1$type[i])){
      a1$type[i] <- paste(gsub(" ", "_", a1$type[i]), "_", as.character(a1$POS[i]), sep = "")}
    if(grepl("silent", a1$type[i])){
      a1$type[i] <- paste(gsub(" ", "_", a1$type[i]), "_", as.character(a1$POS[i]), sep = "")}}
  
  # Pair down to only the required columns.
  a2 <- a1[c("date","type", "VSP")]
  
  # Create the labels with only the mutations of interest.
  a2$mutant_label <- 'empty'
  
  for(j in 1:length(moi)){
    match <- moi[j]
    for(i in 1:length(a2$type)){
      if(a2$type[i] == match){
        a2$mutant_label[i] <- match
      }
      if(a2$type[i] != match && a2$mutant_label[i] == 'empty'){
        a2$mutant_label[i] <- 'other'
      }
    }
  }
  
  # Filter the mutation samples to only include those with the specified rationale.
  a3 <- a2
  a3 <- dplyr::filter(a3, mutant_label != 'other')
  a3$filter <- 'n'
  for(i in 1:length(a3$filter)){
    for(j in 1:length(b3$lab_id)){
      if(a3$VSP[i] == b3$lab_id[j]){
        a3$filter[i] <- 'y'
        break
      }
    }
  }
  a4 <- dplyr::filter(a3, filter == "y")
  
  # Summarize the counts of the mutations for each time bin.
  a4$date_bin <- a4$date
  a4$date_bin <- floor_date(a4$date, time_bin) # Adds a column with date for the week only
  
  # Remove the duplicate samples (two identical samples with the same mutations).
  a4 <- dplyr::distinct(a4)
  
  a5 <- ddply(a4, c("date_bin", "mutant_label"), summarise, n = length(mutant_label)) # Makes a new dataframe with the counts of each lineage for each of the week categories
  
  # Make a list of the unique dates and lineages
  mutant_unique <- unique(a5[c("mutant_label")])
  date_bin_unique <- unique(a5[c('date_bin')])
  
  # Add important lineages to unique lineages. This makes sure the color assigned to lineages is consistent.
  for(li in 1:length(moi)){
    mutant_unique <- mutant_unique %>% add_row(mutant_label = moi[li])
  }
  mutant_unique <- unique(mutant_unique[c("mutant_label")])
  row.names(mutant_unique) <- NULL
  
  # Make the columns of the graph with the repeating lineages for each date.
  mutant_repeat <- mutant_unique
  for(d in 2:lengths(date_bin_unique)){
    mutant_repeat <- Map(c, mutant_repeat, mutant_unique)
  }
  
  # Make the column of the graph with the repeating date, where there are enough dates for each lineage.
  date_repeat <- list(date = date_bin_unique[1,1])
  for(d in 1:lengths(date_bin_unique)){
    for(l in 1:lengths(mutant_unique)){
      date_repeat <-Map(c, date_repeat, date_bin_unique[d,1])
    }
  } 
  lengths(date_repeat)
  
  # Remove the initial repeated date from the creation of this new list.
  date_repeat <- date_repeat$date[-2]
  date_repeat <- data.frame(date_repeat)
  lengths(date_repeat)
  
  # Combine the repeated dates and lineages into a single data frame.
  a6 <- c(date_repeat, mutant_repeat)
  a6$date_repeat <- format(as.POSIXct(a6$date_repeat,format='%Y-%m-%d %H:%M:%S'),format='%Y-%m-%d')
  a6 <- data.frame(a6)
  
  # Make a list of zero counts to be temporary placeholders.  
  count_zero <- rep(c(0),each=length(a6$date))
  count_zero <- list(counts = count_zero)
  a7 <- c(date_repeat, mutant_repeat, count_zero)
  a7 <- data.frame(a7)
  a7$date_repeat <- as.Date(a7$date_repeat)
  # a7$date_repeat <- format(as.POSIXct(a7$date_repeat,format='%Y-%m-%d %H:%M:%S'),format='%Y-%m-%d')
  
  # Loop through the data to determine which dates have counts for each lineage and replace zero placeholders with appropriate counts
  
  for(i in 1:length(a7$date_repeat)){
    #i <- 4
    date_search <- a7$date_repeat[i]
    mutant_label_search <- a7$mutant_label[i]
    data <- subset(a5, a5$date_bin == date_search)
    data <- subset(data, data$mutant_label == mutant_label_search)
    if(length(data$n) == 1){
      a7$counts[i] <- data$n
    }
  }
  
  # Add the genome count for the bin to the mutation bin.
  a7$counts_total <- 0
  a <- a7
  for(i in 1:length(a7$date)){
    for(j in 1:length(b4$sample_date_bin)){
      if(a7$date[i] == b4$sample_date_bin[j]){
        a7$counts_total[i] <- b4$n[j]
      }
    }
  }
  
  # Determine the percent of the variant for that time.
  a7$percent <- a7$counts/a7$counts_total * 100
  write.csv(a7, paste(gene,"_Bayesian/",gene, "_", rationale, ".csv", sep = "")) # XXX
  
  # Prepare the data frame to be plotted
  a8 <- b4
  a8 <- subset(a8,select = -c(counts))
  
  colnames(a8)[colnames(a8) == "sample_date_bin"] <- "date"
  for(h in 1:length(moi)){
    a8[moi[h]] <- 0
    for(i in 1:length(a8$date)){
      # h <- 1
      # i <- 2
      date_search <- a8$date[i]
      mutant_label_search <- moi[h]
      data <- subset(a7, a7$date_repeat == date_search)
      data <- subset(data, data$mutant_label == mutant_label_search)
      if(length(data$percent) == 1){
        a8[i,h+2] <- data$percent
      }
    }
  }
  a9 <- subset(a8,select = -c(n))
  return(a9)
}

# Determine the mutations of interest.
determine_moi <- function(file_mutation, gene, rationale){
  
  # Determine the most important mutations to track.
  mut1 <- read_excel(file_mutation)
  
  # Filter to look at only the gene of interest.
  mut1 <- dplyr::filter(mut1, genes == gene)
  
  ### Correct the deletions for the mutations.
  for(i in 1:length(mut1$VSP)){
    if(grepl("del", mut1$type[i])){
      mut1$type[i] <- paste(gsub(" ", "_", mut1$type[i]), "_", as.character(mut1$POS[i]), sep = "")}
    if(grepl("silent", mut1$type[i])){
      mut1$type[i] <- paste(gsub(" ", "_", mut1$type[i]), "_", as.character(mut1$POS[i]), sep = "")}}
  
  # Make a column for the mutation including the gene and type.
  mut1$mutation <- mut1$type
  
  # Get the counts for each mutation.
  mut2 <- data.frame(table(mut1$mutation))
  colnames(mut2) <- c("mutation", "count")
  
  # Order the mutation data frame from highest to lowest.
  mut3 <- mut2[order(-mut2$count),]
  row.names(mut3) <- NULL
  
  # Use only mutations that occur 5% of genomes or more.
  # Determine the number of samples.
  samples <- length(unique(mut1$VSP))
  cutoff <- round(0.05*samples)
  mut4 <- subset(mut3, count >= cutoff)
  moi <- as.character(mut4$mutation)
  return(moi)
}

# Plot the mutations.
dotplot <- function(spike_surv, spike_vbt, gene){
  for(i in 2:length(colnames(spike_surv))){
    
    # Get the data for the Surv
    #i <- 3
    data1 <- spike_surv[1]
    data1 <- cbind(data1,spike_surv[i])
    colnames(data1) <- c("date", "surv")
    data1$surv <- data1$surv/100
    
    # Get the data for the VBT
    mutation <- colnames(spike_surv[i])
    data2 <- spike_vbt[1]
    vbt_mutations <- colnames(spike_vbt)
    k <- 0
    # Determine which mutation should be drawn from the vaccine breakthrough group.
    for(j in 1:length(vbt_mutations)){
      if(vbt_mutations[j]==mutation){
        k <- j
      }}
    
    # Pull the vaccine breakthrough for the associated mutation.
    if(k > 0){
      data2 <- cbind(data2,spike_vbt[k])
      colnames(data2) <- c("date", "vbt")
      data2$vbt <- data2$vbt/100
    }else{
      data2$vbt <- 0
    }
    
    # Combine the data into one data frame.
    date_unique <- c(data1$date,data2$date)
    date_unique <- unique(date_unique)
    
    # Make a new data blank data frame.
    data3 <- data.frame(matrix(NA, nrow = length(date_unique), ncol = 3))
    colnames(data3) <- c("date", "surv", "vbt")
    
    # Fill in the date.
    data3$date <- date_unique
    
    # Fill in the Surveillance proportions.
    for(j in 1:length(data1$date)){
      for(k in 1:length(data3$date)){
        if(data1$date[j] == data3$date[k]){
          data3$surv[k] <- data1$surv[j]}}}
    
    # Fill in the breakthrough proportions.
    for(j in 1:length(data2$date)){
      for(k in 1:length(data3$date)){
        if(data2$date[j] == data3$date[k]){
          data3$vbt[k] <- data2$vbt[j]}}}
    
    # Replace all NA with 0.
    data4 <- data3
    data4[is.na(data4)] <- -1
    
    # Plot the data.
    data <- data4
    
    dir.create(gene)
    title <- paste(gene, "/mutation-proportion_", mutation,".pdf", sep = "")
    pdf(file=title, width = 8.5, height = 11)
    title <- paste("Proportion of Mutation: ", mutation, sep = "")
    dotchart(data$surv, main = title, pch = 21, labels = data$date, bg = "blue",
             pt.cex = 2.0, xlim = c(0,1) + c(0, 0))
    points(data$vbt, 1:nrow(data), col = "red", pch = 19, cex = 1.5)
    par(xpd=TRUE)
    legend(0,-1,  
           legend = c("Surveillance", "Vaccine Breakthrough"), 
           col = c("blue", "red"), 
           pch = c(19,19), 
           bty = "n", 
           pt.cex = 2, 
           cex = 1.2, 
           text.col = "black", 
           horiz = F , 
           inset = c(0.1, 0.1))
    par(xpd=FALSE)
    coord_flip()
    dev.off()
  }
}

# Define function to get the dates for the bayesian graph.
bayesian_date <- function(){
  a1 <- read.csv(paste(gene,"_Bayesian/", gene, "_surveillance.csv", sep = ""))
  b1 <- read.csv(paste(gene,"_Bayesian/", gene, "_vaccine breakthrough.csv", sep = ""))
  c1 <- read.csv(paste(gene,"_Bayesian/", gene, "_s drop.csv", sep = ""))
  d1 <- read.csv(paste(gene,"_Bayesian/", gene, "_other.csv", sep = ""))
  
  # Determine the max and min dates.
  a1$date_repeat <- as.Date(a1$date_repeat)
  b1$date_repeat <- as.Date(b1$date_repeat)
  c1$date_repeat <- as.Date(c1$date_repeat)
  d1$date_repeat <- as.Date(d1$date_repeat)
  date1 <- c(a1$date_repeat, b1$date_repeat, c1$date_repeat, d1$date_repeat)
  date2 <- unique(date1)
  date2 <- as.Date(date2)
  date_max <- max(date2)
  date_min <- min(date2)
  
  # Get all of the dates in the date range. 
  date_length <- as.numeric(abs(as.Date(date_min) - as.Date(date_max))) + 1
  date3 <- seq(as.Date(date_min), by = "day", length.out = date_length)
  date4 <- floor_date(date3, time_bin)
  date4 <- unique(date4)
  return(date4)
}

# Define function to prepare data for bayesian input.
bayesian_prep <- function(){
  a1 <- read.csv(paste(gene,"_Bayesian/", gene, "_surveillance.csv", sep = ""))
  b1 <- read.csv(paste(gene,"_Bayesian/", gene, "_vaccine breakthrough.csv", sep = ""))
  c1 <- read.csv(paste(gene,"_Bayesian/", gene, "_s drop.csv", sep = ""))
  d1 <- read.csv(paste(gene,"_Bayesian/", gene, "_other.csv", sep = ""))
  
  # Determine the max and min dates.
  a1$date_repeat <- as.Date(a1$date_repeat)
  b1$date_repeat <- as.Date(b1$date_repeat)
  c1$date_repeat <- as.Date(c1$date_repeat)
  d1$date_repeat <- as.Date(d1$date_repeat)
  date1 <- c(a1$date_repeat, b1$date_repeat, c1$date_repeat, d1$date_repeat)
  date2 <- unique(date1)
  date2 <- as.Date(date2)
  date_max <- max(date2)
  date_min <- min(date2)
  
  # Get all of the dates in the date range. 
  date_length <- as.numeric(abs(as.Date(date_min) - as.Date(date_max))) + 1
  date3 <- seq(as.Date(date_min), by = "day", length.out = date_length)
  date4 <- floor_date(date3, time_bin)
  date4 <- unique(date4)
  date_length <- length(date4)
  date4 <- rep(date4, each = 2)
  date_real_repeat <- rep(date4, each = length(unique(a1$mutant_label)))
  date <- 1:date_length
  date <- rep(date, each = 2)
  date_repeat <- rep(date, each = length(unique(a1$mutant_label)))
  
  # Adapt the data frame to have all the dates.
  bay <- data.frame(matrix(NA, nrow = length(date_repeat), ncol = 6))
  colnames(bay) <- c("date_real","date", "rationale", "mutation", "category", "count")
  bay$date_real <- date_real_repeat
  bay$date <- date_repeat
  mutations <- rep(unique(a1$mutant_label), date_length)
  mutations <- rep(mutations, each = 2)
  bay$mutation <- mutations
  
  # Populate the data frame with random.
  bay1 <- bay
  for(i in 1:length(bay1$date)){
    #i <- 179
    sub <- subset(a1, as.character(a1$date_repeat) == (bay1$date_real[i]))
    sub <- subset(sub, sub$mutant_label == bay1$mutation[i])
    if((i %% 2) == 0) {
      bay1$category[i] <- "FALSE"
      bay1$count[i] <- sub$counts_total[1] - sub$counts[1]
      bay1$rationale[i] <- "random"
    } else {
      bay1$category[i] <- "TRUE"
      bay1$count[i] <- sub$counts[1]
      bay1$rationale[i] <- "random"
    }
  }
  bay1[is.na(bay1)] <- 0
  
  # Populate the data frame with vaccine.
  bay2 <- bay
  for(i in 1:length(bay2$date)){
    #i <- 1
    sub <- subset(b1, b1$date_repeat == bay2$date_real[i])
    sub <- subset(sub, sub$mutant_label == bay2$mutation[i])
    if((i %% 2) == 0) {
      bay2$category[i] <- "FALSE"
      bay2$count[i] <- sub$counts_total[1] - sub$counts[1]
      bay2$rationale[i] <- "vaccine"
    } else {
      bay2$category[i] <- "TRUE"
      bay2$count[i] <- sub$counts[1]
      bay2$rationale[i] <- "vaccine"
    }
  }
  bay2[is.na(bay2)] <- 0
  
  # Populate the data frame with s drop.
  bay3 <- bay
  for(i in 1:length(bay3$date)){
    #i <- 1
    sub <- subset(c1, c1$date_repeat == bay3$date_real[i])
    sub <- subset(sub, sub$mutant_label == bay3$mutation[i])
    if((i %% 2) == 0) {
      bay3$category[i] <- "FALSE"
      bay3$count[i] <- sub$counts_total[1] - sub$counts[1]
      bay3$rationale[i] <- "s drop"
    } else {
      bay3$category[i] <- "TRUE"
      bay3$count[i] <- sub$counts[1]
      bay3$rationale[i] <- "s drop"
    }
  }
  bay3[is.na(bay3)] <- 0
  
  # Populate the data frame with other.
  bay4 <- bay
  for(i in 1:length(bay4$date)){
    #i <- 1
    sub <- subset(d1, d1$date_repeat == bay4$date_real[i])
    sub <- subset(sub, sub$mutant_label == bay4$mutation[i])
    if((i %% 2) == 0) {
      bay4$category[i] <- "FALSE"
      bay4$count[i] <- sub$counts_total[1] - sub$counts[1]
      bay4$rationale[i] <- "other"
    } else {
      bay4$category[i] <- "TRUE"
      bay4$count[i] <- sub$counts[1]
      bay4$rationale[i] <- "other"
    }
  }
  bay4[is.na(bay4)] <- 0
  
  # Combine the data frames to one data frame.
  bay5 <- rbind(bay1, bay2, bay3, bay4)
  
  return(bay5)
}

# Combine the pdf plots together.
combine_pdf <- function(dir){
  # Determine the pdfs that are present.
  files1 <- list.files(dir)
  files2 <- ""
  for(i in 1:length(files1)){
    if(grepl(".pdf", files1[i])){
      files2 <- c(files2, files1[i])
    }
  }
  files2 <- files2[-1]
  
  for(i in 1:length(files2)){
    files2[i] <- paste(dir, "/", files2[i], sep = "")
  }
  
  # Combine them with the other one
  library(pdftools)
  title <- paste(dir, "_Summary_v5.pdf", sep = "")
  pdf_combine(files2, output = title)
}

# Combine the summary statistics files together.
combine_statistics_each_gene <- function(gene, moi){
  dir <- paste(gene,"_Bayesian", sep = "")
  # Determine the summary csvs that are present.
  files1 <- list.files(dir)
  files2 <- ""
  for(i in 1:length(files1)){
    if(grepl("Odds-Enrichment-Summary.csv", files1[i])){
      files2 <- c(files2, files1[i])
    }
  }
  files2 <- files2[-1]
  
  # Make a data frame that will hold all of the contents of the summary statistics.
  sum1 <- data.frame(matrix(NA, nrow = length(files2), ncol = 6))
  colnames(sum1) <- c("Mutation" , "Mean", "Lower 95% CrI", "Upper 95% CrI", "Gene", "Position")
  
  # Fill the data frame with the mutations.
  for(i in 1:length(moi)){
    for(j in 1:length(files2)){
      if(grepl(moi[i], files2[j])){
        # i <- 6
        sum1$Mutation[i] <- moi[i]
        # print(j)
        temp <- read.csv(paste(dir, "/", files2[j], sep = ""))
        sum1$Mean[i] <- temp$Mean[1]
        sum1$`Lower 95% CrI`[i] <- temp$Lower.95..CrI[1]
        sum1$`Upper 95% CrI`[i] <- temp$Upper.95..CrI[1]
        sum1$Gene[i] <- gene
        sum1$Position[i] <- temp$Position[1]
      }}}
  # Save the summary csv.
  write.csv(sum1, paste(dir, "_Summary_v5_Odds-Enrichment.csv", sep =""), quote=FALSE, row.names = F)

  # Make a plot for each of the distributions.
  files3 <- ""
  for(i in 1:length(files1)){
    if(grepl(".rds", files1[i])){
      temp <- paste(dir, "/", files1[i], sep = "") 
      files3 <- c(files3, temp)
    }
  }
  files3 <- files3[-1]
  
  # Determine the number of rows required
  moi_length <- length(moi)
  nrows <- ceiling(moi_length/3)
  
  # Save the generated figures.
  pdf(paste(gene, "_Bayesian_Summary_v5_Density_Plots.pdf", sep = ""), width=10, height=nrows*2) # Use this line to make the summary for all the mutations.
  par(mfrow = c(nrows,3), mar = c(2,1,2,1))
  
  # Which mutations failed to converge:
  #black_list <- "D614G, silent_3037, P314L"
  black_list <- ""
  
  # Make a plot for each of the mutations and combine them together.
  for(i in 1:length(moi)){
    # i <- 1
    for(j in 1:length(files3)){
      if(grepl(moi[i],files3[j])){
        if(grepl(moi[i], black_list) == F){
        # i <- 1
        # j <- 13
          print(files3[j])
        stan <- readRDS(files3[j])
        temp <- as.matrix(stan)
        dens <- density(exp(temp[,'vaccineChange']))
        plot(dens, xlab='', ylab = "", bty='n', yaxs='i', log = "x", xlim = c(0.01,100), ylim = c(0,3.5), xaxt = 'n', xaxs = 'i', yaxt = 'n', 
             main = moi[i] # Use this line to make the summary for all the mutations.
        ) 
        polygon(dens$x, dens$y, 
                col = "grey") # Toggle between this selection (all grey) and the line above (by predifined colors)
        abline(v=1,lty=2)
        Axis(at = c(0.001, 0.01, 0.1, 1, 10, 100, 1000), side = 1, labels = c("", expression(10^-2), "",  1, "",  expression(10^2), ""), las = 1, lwd.ticks = 2)
        Axis(at = c(0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
             , side = 1, labels = c("","","","","","","","","","","","","","","","","","","","","","","","","",""), las = 1)
        Axis(at = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90, 200, 300, 400, 500, 600, 700, 800, 900), 
             side = 1, labels = c("","","","","","","","","","","","","","","","","","","","","","","","","",""), las = 1)
      }}
    }
  }
  dev.off()
}

# Combine the summary statistics from all of the genes into one file.
combine_statistics_all_genes <- function(){
  # Collect all the data from the data frame.
  df_spike <- read.csv("S_Bayesian_Summary_v5_Odds-Enrichment.csv")
  df_orf1ab <- read.csv("ORF1ab_Bayesian_Summary_v5_Odds-Enrichment.csv")
  df_n <- read.csv("N_Bayesian_Summary_v5_Odds-Enrichment.csv")
  df_m <- read.csv("M_Bayesian_Summary_v5_Odds-Enrichment.csv")
  
  # Combine all the data together.
  df <- rbind(df_spike, df_orf1ab, df_n, df_m)
  
  # Sort the data frame by nucleotide position.
  df <- df[order(df$Position),]
  
  # Mark the significant mutations.
  df$Significant <- ""
  for(i in 1:length(df$Mutation)){
    if(df$Mean[i] > 1){
      if(df$Lower.95..CrI[i] > 1){
        df$Significant[i] <- "significant increase"
      }
    }else
      if(df$Upper.95..CrI[i] < 1){
        df$Significant[i] <- "significant decrease"
      }
  }
  
  # Add the protein affected by the mutation. This data comes from https://genome.ucsc.edu/cgi-bin/hgTracks?db=wuhCor1&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=NC_045512v2%3A1%2D29903&hgsid=1179262955_MHChLGp9MpiVCd4iW2gBAbg7LoZk
  df$Protein <- ""
  for(i in 1:length(df$Mutation)){
    if(df$Position[i] >= 266 && df$Position[i] < 805){df$Protein[i] <- "Nsp1"}
    if(df$Position[i] >= 806 && df$Position[i] < 2719){df$Protein[i] <- "Nsp2"}
    if(df$Position[i] >= 2720 && df$Position[i] < 8554){df$Protein[i] <- "Nsp3"}
    if(df$Position[i] >= 8555 && df$Position[i] < 10054){df$Protein[i] <- "Nsp4"}
    if(df$Position[i] >= 10055 && df$Position[i] < 10972){df$Protein[i] <- "3CL-PRO"}
    if(df$Position[i] >= 10973 && df$Position[i] < 11842){df$Protein[i] <- "Nsp6"}
    if(df$Position[i] >= 11843 && df$Position[i] < 12091){df$Protein[i] <- "Nsp7"}
    if(df$Position[i] >= 12092 && df$Position[i] < 12685){df$Protein[i] <- "Nsp8"}
    if(df$Position[i] >= 12686 && df$Position[i] < 13024){df$Protein[i] <- "Nsp9"}
    if(df$Position[i] >= 13025 && df$Position[i] < 13441){df$Protein[i] <- "Nsp10"}
    if(df$Position[i] >= 13442 && df$Position[i] < 16236){df$Protein[i] <- "Pol"}
    if(df$Position[i] >= 16237 && df$Position[i] < 18039){df$Protein[i] <- "Hel"}
    if(df$Position[i] >= 18040 && df$Position[i] < 19620){df$Protein[i] <- "ExoN"}
    if(df$Position[i] >= 19621 && df$Position[i] < 20658){df$Protein[i] <- "Nsp15"}
    if(df$Position[i] >= 20659 && df$Position[i] < 21552){df$Protein[i] <- "Nsp16"}
    if(df$Position[i] >= 21563 && df$Position[i] < 25381){df$Protein[i] <- "Spike"}
    # if(df$Position[i] >= 21599 && df$Position[i] < 23617){df$Protein[i] <- "S1"}
    # if(df$Position[i] >= 23618 && df$Position[i] < 25381){df$Protein[i] <- "S2"}
    if(df$Position[i] >= 25393 && df$Position[i] < 26217){df$Protein[i] <- "ORF3a"}
    if(df$Position[i] >= 26245 && df$Position[i] < 26469){df$Protein[i] <- "Envelope"}
    if(df$Position[i] >= 26523 && df$Position[i] < 27188){df$Protein[i] <- "Membrane"}
    if(df$Position[i] >= 27202 && df$Position[i] < 27384){df$Protein[i] <- "ORF6"}
    if(df$Position[i] >= 27439 && df$Position[i] < 27756){df$Protein[i] <- "ORF7a"}
    if(df$Position[i] >= 27756 && df$Position[i] < 27884){df$Protein[i] <- "ORF7b"}
    if(df$Position[i] >= 27939 && df$Position[i] < 28256){df$Protein[i] <- "ORF8"}
    if(df$Position[i] >= 28274 && df$Position[i] < 29530){df$Protein[i] <- "Nucleocapsid"}
    if(df$Position[i] >= 29558 && df$Position[i] < 29671){df$Protein[i] <- "ORF10"}
  }
  
  # Reorder the colunms to be more logical.
  df <- df[, c(6, 5, 8, 1, 2, 3, 4, 7)]
  
  colnames(df) <- c("Genomic Position", "Gene", "Protein", "Mutation", "Mean", "Lower 95% CrI", "Upper 95% CrI",  "Significant")
  # df <- df[order(df$Significant,decreasing = TRUE),] # This will make the most significant mutations come to the top. 
  
  # Save the summary csv.
  write.csv(df, paste("Summary_v5_Odds-Enrichment.csv", sep =""), quote=FALSE, row.names = F)
}

# Combine all predicted point mutations plots.
combine_predicted_plots <- function(date){
  date <- date_all
  genes <- c("S", "ORF1ab", "N", "M")
  
  # # Adjust the dates so that they match
  date <- as.POSIXct(date)
  
  files3 <- ""
  for(i in 1:length(genes)){
    gene <- genes[i]
    
    # Determine the proprtion over time files that are present.
    dir <- paste(gene,"_Bayesian", sep = "")
    files1 <- list.files(dir)
    files2 <- ""
    for(i in 1:length(files1)){
      if(grepl("_Proportion-Over-Time.csv", files1[i])){
        files2 <- c(files2, files1[i])}}
    files2 <- files2[-1]
    for(i in 1:length(files2)){
      files2[i] <- paste(dir, "/", files2[i], sep = "")
    }
    
    # Save all files for all genes to one variable.
    files3 <- c(files3,files2)
  }
  files3 <- files3[-1]
  
  # Determine all of the mutations that are present.
  mutations <- gsub("Bayesian/", "", files3)
  mutations <- gsub("_Proportion-Over-Time.csv", "", mutations)
  
  # Make a blank data frame with all of the data. 
  data1 <- data.frame(matrix(NA, nrow = length(date), ncol = length(files3)+1))
  colnames(data1) <- c("date", mutations)
  
  # Fill in the dates for the data frame of all mutations.
  data1$date <- date
  
  # Populate the data frame with the proportion from the proportion over time data file.
  data2 <- data1
  data2$date <- as.POSIXct(data2$date)
  for(i in 1:length(files3)){
    temp <- read.csv(files3[i])
    temp$Date <- as.POSIXct(temp$Date)
    for(j in 1:length(temp$Date)){
      for(k in 1:length(date)){
        if(temp$Date[j] == date[k]){
          data2[k,i+1] <- temp$Proportion[j]
        }
      }
    }
  }
  data2$date <- as.character(as.Date(data2$date))
  
  # Get the data in the correct orientation for analysis.
  library(reshape2)
  data3 <- reshape2::melt(data2)
  data3$date <- as.Date(data3$date)
  colnames(data3) <- c("date", "mutation", "proportion")
  
  # Determine the colors
  # Assign mutations to either alpha, delta, or both. The default will be black.
  alpha <- "ORF1ab_del_9_11288N_R203K,N_G204R,S_N501Y,S_P681H,S_del_3_21991,S_del_6_21765,S_T716I,ORF1ab_silent_14676,ORF1ab_silent_5986,S_S982A,ORF1ab_silent_913,ORF1ab_T1001I,ORF1ab_I2230T,ORF1ab_silent_15279,N_S235F,S_A570D,N_D3L,ORF1ab_A1708D,ORF1ab_silent_16176,ORF8_Y73C,ORF8_Q27*,S_D1118H,ORF8_R52I"
  
  delta <- "ORF1ab_T3255I, intergenic_del_1_28271, N_D377Y,M_I82T,S_L452R,S_D950N,N_R203M,ORF1ab_G662S,
  ORF1ab_P1000L,ORF3a_S26L,S_P681R,N_D63G,S_T478K,S_T19R,ORF7a_T120I,ORF7a_V82A,S_del_6_22029,ORF8_del_6_28248,
  ORF1ab_A1918V,ORF1ab_P2046L,N_G215C,ORF1ab_P2287S,ORF1ab_silent_8986,ORF1ab_T3646A,
  ORF1ab_silent_11332,ORF1ab_A1306S,ORF1ab_V2930L,ORF7b_T40I,ORF1ab_P1640L,ORF1ab_V3718A,
  ORF1ab_P309L,  ORF1ab_silent_1267,ORF1ab_silent_12946,ORF1ab_P1570L,ORF7a_L116F,
  ORF1ab_T3750I,ORF1ab_silent_13019,S_A222V,ORF1ab_S1898F,S_K417N,ORF1ab_silent_5584,ORF1ab_S538A,
  ORF1ab_M2417I,ORF7b_E33V,ORF1ab_D2183Y,S_silent_25279,S_W258L,ORF1ab_I1257V,S_silent_24208,N_silent_29095,
  ORF3a_L53F,ORF7a_silent_27507,ORF1ab_silent_2509,N_silent_29227,S_V1104L,ORF1ab_K669N,ORF1ab_silent_14187,
  ORF1ab_I3731V,ORF3a_silent_25806,ORF1ab_silent_3055,ORF1ab_A1146T,S_V70F,ORF1ab_silent_8062,
  ORF1ab_silent_20436,ORF3a_D222Y,ORF1ab_P3447H,ORF1ab_silent_17406,ORF1ab_silent_18468,ORF3a_E239Q,
  N_silent_29509,N_silent_29050,ORF1ab_silent_18360,N_A208S,ORF1ab_A599V"
  
  both <- "ORF1ab_silent_3037,S_D614G,intergenic_NA,ORF1ab_P314L"
  colors <- c("")
  for(i in 1:length(mutations)){
    temp_color <- "black"
    if(grepl(mutations[i],alpha)){
      temp_color <- "#33a02c"}
    if(grepl(mutations[i],delta)){
      temp_color <- "#f70028"
    }
    if(grepl(mutations[i],both)){
      temp_color <- "#f57e4e"}
    colors <- c(colors, temp_color)
  }
  colors <- colors[-1]
  
  h <- 5
  w <- h*1.618
  pdf("Summary_v5_Point-Mutations.pdf", height = h, width = w)
  ggplot(data = data3, aes(x=date, y=proportion)) + 
    geom_line(aes(colour=mutation)) + 
    ggtitle("Proportion of Point Mutations in Surveillance Genomes") + 
    theme(legend.position = "none", panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
          panel.background = element_rect(fill = "white", colour = "white"),
          plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 45, hjust=1.1)) + 
    scale_x_date(date_labels = "%b %Y", limits = c(as.Date("2020-05-03"), NA), expand = c(0, 0), date_breaks = "1 months") + 
    scale_y_continuous(expand = c(0, 0), limits = c(0,1)) +
    xlab("Date") + ylab("Proportion of Genomes") + 
    scale_color_manual(values=colors)
  dev.off()
  }

# Define function for bayesian calculations.
bayesian_calculations <- function(){
  for(i in 1:length(moi)){
  # i <- 37
    mutation_target <- moi[i]
    # mutation_target <- "N501Y"
    
    # Reformat into the dat structure.
    bay6 <- subset(bay5, bay5$mutation == mutation_target)
    bay7 <- structure(bay6$count,
                      .Dim = c(2L, length(bay6$date)/8, 4L), # 8 comes from 4 different groups (random, vaccine, s drop, other) and 2 categories (True and False); 2x4=8
                      .Dimnames = list(c("TRUE", "FALSE"), unique(bay6$date), 
                                       c("random", "vaccine", "s drop", "other")))
    dat <- bay7
    source('functions.R')
    options(mc.cores = parallel::detectCores())
    options(mc.cores = 17)
    rstan::rstan_options(auto_write = TRUE)
    
    # Uncomment these lines to run the model again.
    mod <- rstan::stan_model("model_mutations.stan")
    fit<-runMutationStan(dat[,,'random'],dat[,,'s drop'],dat[,,'vaccine'],mod, nIter=2000)
    # Save the stan object.
    dir <- paste(gene,"_Bayesian", sep = "")
    dir.create(dir)
    saveRDS(fit$stan, paste(dir, "/bayesian_", mutation_target,".rds", sep = ""))
    fit_stan <- fit$stan
    
    # Plot the graph of the specific mutation with the confidence interval
    # Prepare the stan data object.
    mut_stan <- as.matrix(fit_stan)
    mut_stan <- invLogit(mut_stan)
    # Make the credible interval
    cri<-apply(mut_stan[,grep('means',colnames(mut_stan))],2,quantile,c(.025,.975))
    mut_stan <- apply(mut_stan[,grep('means',colnames(mut_stan))],2,mean)
    
    # Convert the week counts back to the week names.
    mut_stan2 <- data.frame(matrix(NA, nrow = length(mut_stan), ncol = 4))
    colnames(mut_stan2) <- c("Date", "Proportion", "2.5%", "97.5%")
    mut_stan2$Date <- as.POSIXct(date)
    # Add the proportions
    for(i in 1:length(mut_stan)){
      mut_stan2$Proportion[i] <- mut_stan[i]
    }
    # Add the credible interval
    for(i in 1:(length(cri)/2)){
      mut_stan2$`2.5%`[i] <- cri[1,i]
      mut_stan2$`97.5%`[i] <- cri[2,i]
    }
    
    # Save the estimated proportion as .csv file.
    write.csv(mut_stan2, paste(dir, "/", mutation_target, "_Proportion-Over-Time.csv", sep =""), quote=FALSE, row.names = F)
    
    # Plot the data for a smoothing curve.
    library(ggplot2)
    
    plot_occurance <- ggplot() +
      geom_line(data = mut_stan2, aes(x=Date, y=Proportion, color="Proportion"), lwd = 1.01) +
      geom_ribbon(data = mut_stan2, aes(x=Date, ymin = `2.5%`, ymax = `97.5%`), alpha = 0.1, fill = "black")+
      scale_color_manual(name = "", values = c("Proportion" = "black"))+
      scale_x_datetime(date_labels = "%b %Y", date_breaks = "month") +
      ylim(0, 1) + 
      ylab("Proportion") + 
      ggtitle(mutation_target) + 
      theme(axis.text = element_text(size = 12), axis.text.x = element_text(angle = 45, hjust=1.1),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.position = "none",plot.title = element_text(hjust = 0.5))
    plot_occurance   
    
    # Plot the data for the amount that the mutation stnds out from the surveillance data set. 
    # Prepare the data frame
    mut_stan3 <- as.matrix(fit_stan)
    mut_stan4 <- data.frame(mut_stan3)
    mut_stan5 <- mut_stan4$vaccineChange
    mut_stan5 <- exp(mut_stan5)       # convert to fold change
    mut_stan6 <- data.frame(matrix(NA, nrow = length(mut_stan5), ncol = 2))
    colnames(mut_stan6) <- c("count", "Fold Enrichment in Odds of Vaccine Breakthrough")
    mut_stan6$`Fold Enrichment in Odds of Vaccine Breakthrough` <- mut_stan5
    mut_stan6$count <- "Enrichment"
    
    # Save odds enrichment.
    write.csv(mut_stan6, paste(dir, "/", mutation_target, "_Odds-Enrichment.csv", sep =""), quote=FALSE, row.names = F)
    
    # Determine the credible interval and save it as a file.
    mut_stan7 <- mut_stan6[order(mut_stan6$`Fold Enrichment in Odds of Vaccine Breakthrough`),]
    row.names(mut_stan7) <- NULL  
    count_95 <- round(length(mut_stan7$count)*0.975)
    count_05 <- round(length(mut_stan7$count)*0.025)
    mut_stan8 <- data.frame(matrix(NA, nrow = 1, ncol = 6))
    colnames(mut_stan8) <- c("Mutation", "Mean", "Lower 95% CrI", "Upper 95% CrI", "Gene", "Position")
    mut_stan8$Mutation <- mutation_target
    mut_stan8$`Upper 95% CrI` <- mut_stan7$`Fold Enrichment in Odds of Vaccine Breakthrough`[count_95]
    mut_stan8$`Lower 95% CrI` <- mut_stan7$`Fold Enrichment in Odds of Vaccine Breakthrough`[count_05]
    mut_stan8$Mean <- mean(mut_stan7$`Fold Enrichment in Odds of Vaccine Breakthrough`)
    mut_stan8$Gene <- gene
    # Fill in the position for the mutation
    mut1 <- read_excel(file_mutation)
    mut1 <- subset(mut1, genes == gene)
    for(i in 1:length(mut1$VSP)){
      if(grepl("del", mut1$type[i])){
        mut1$type[i] <- paste(gsub(" ", "_", mut1$type[i]), "_", as.character(mut1$POS[i]), sep = "")}
      if(grepl("silent", mut1$type[i])){
        mut1$type[i] <- paste(gsub(" ", "_", mut1$type[i]), "_", as.character(mut1$POS[i]), sep = "")}}
    mut1 <- subset(mut1, type == mutation_target)
    mut_stan8$Position <- mut1$POS[1]
    write.csv(mut_stan8, paste(dir, "/", mutation_target, "_Odds-Enrichment-Summary.csv", sep =""), quote=FALSE, row.names = F)
    
    # Plot enrichment.
    require(scales)
    plot_enrichment <- ggplot(mut_stan6, aes(y=`Fold Enrichment in Odds of Vaccine Breakthrough`)) + 
      ggtitle("Enrichment") + 
      geom_boxplot(fill="black", alpha=0.2, outlier.shape = NA) + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            axis.ticks.x = element_blank(), axis.text.x = element_blank()) +  
      #scale_y_continuous(limits = c(0.01, 100, breaks = seq(-5, 5, by = 0.5))) + 
      scale_y_log10(limits = c(0.01, 100), labels = scales::number_format(accuracy = 0.01))
    plot_enrichment
    
    # Combine the two plots together.
    library(grid)
    library(gridExtra)
    plots <- grid.arrange(plot_occurance, plot_enrichment, nrow = 1, widths = c(5, 1))
    # plots
    
    # save the plot.
    title <- paste(dir, "/bayesian_", mutation_target,".pdf", sep = "")
    pdf(file=title, width = 8, height = 4)
    grid.arrange(plot_occurance, plot_enrichment, nrow = 1, widths = c(4,1.1))
    dev.off()
  }
  
  # Combine all of the bayesian plots into one pdf
  combine_pdf(paste(gene,"_Bayesian", sep = ""))
}

# Defining variables to graph mutations of interest______________________________________________________________________________

# Shared Files
file_mutation   <- "2021-09-23_mutations.xlsx"
file_genome     <- "2021-09-25_genomeMetaData.xlsx"

# Spike
time_bin        <- "week"
gene            <- "S"
rationale <-  "surveillance"
dir.create(paste(gene,"_Bayesian/", sep = ""))
moi <- determine_moi(file_mutation, gene, rationale)
spike_surv <- curate_moi(file_mutation, file_genome, time_bin, gene, title_graph, title_legend, rationale, moi)
rationale <- "vaccine breakthrough"
spike_vbt <- curate_moi(file_mutation, file_genome, time_bin, gene, title_graph, title_legend, rationale, moi)
rationale <- "s drop"
spike_sdrop <- curate_moi(file_mutation, file_genome, time_bin, gene, title_graph, title_legend, rationale, moi)
rationale <- "other"
spike_other <- curate_moi(file_mutation, file_genome, time_bin, gene, title_graph, title_legend, rationale, moi)

date <- bayesian_date()
date_all <- date
bay5 <- bayesian_prep()

# Perform bayesian calculations.
bayesian_calculations()
combine_statistics_each_gene(gene, moi)

# Orf1ab
time_bin        <- "week"
gene            <- "ORF1ab"
rationale <-  "surveillance"
dir.create(paste(gene,"_Bayesian/", sep = ""))
moi <- determine_moi(file_mutation, gene, rationale)
ORF1ab_surv <- curate_moi(file_mutation, file_genome, time_bin, gene, title_graph, title_legend, rationale, moi)
rationale <- "vaccine breakthrough"
ORF1ab_vbt <- curate_moi(file_mutation, file_genome, time_bin, gene, title_graph, title_legend, rationale, moi)
rationale <- "s drop"
ORF1ab_sdrop <- curate_moi(file_mutation, file_genome, time_bin, gene, title_graph, title_legend, rationale, moi)
rationale <- "other"
ORF1ab_other <- curate_moi(file_mutation, file_genome, time_bin, gene, title_graph, title_legend, rationale, moi)

date <- bayesian_date()
bay5 <- bayesian_prep()

# Perform bayesian calculations.
bayesian_calculations()
combine_statistics_each_gene(gene, moi)

# Nucleocapsid
time_bin        <- "week"
gene            <- "N"
rationale <-  "surveillance"
dir.create(paste(gene,"_Bayesian/", sep = ""))
moi <- determine_moi(file_mutation, gene, rationale)
nucleocapsid_surv <- curate_moi(file_mutation, file_genome, time_bin, gene, title_graph, title_legend, rationale, moi)
rationale <- "vaccine breakthrough"
nucleocapsid_vbt <- curate_moi(file_mutation, file_genome, time_bin, gene, title_graph, title_legend, rationale, moi)
rationale <- "s drop"
nucleocapsid_sdrop <- curate_moi(file_mutation, file_genome, time_bin, gene, title_graph, title_legend, rationale, moi)
rationale <- "other"
nucleocapsid_other <- curate_moi(file_mutation, file_genome, time_bin, gene, title_graph, title_legend, rationale, moi)

date <- bayesian_date()
bay5 <- bayesian_prep()

# Perform bayesian calculations.
bayesian_calculations()
combine_statistics_each_gene(gene, moi)

# Membrane
time_bin        <- "week"
gene            <- "M"
rationale <-  "surveillance"
dir.create(paste(gene,"_Bayesian/", sep = ""))
moi <- determine_moi(file_mutation, gene, rationale)
nucleocapsid_surv <- curate_moi(file_mutation, file_genome, time_bin, gene, title_graph, title_legend, rationale, moi)
rationale <- "vaccine breakthrough"
nucleocapsid_vbt <- curate_moi(file_mutation, file_genome, time_bin, gene, title_graph, title_legend, rationale, moi)
rationale <- "s drop"
nucleocapsid_sdrop <- curate_moi(file_mutation, file_genome, time_bin, gene, title_graph, title_legend, rationale, moi)
rationale <- "other"
nucleocapsid_other <- curate_moi(file_mutation, file_genome, time_bin, gene, title_graph, title_legend, rationale, moi)

date <- bayesian_date()
bay5 <- bayesian_prep()

# Perform bayesian calculations.
bayesian_calculations()
combine_statistics_each_gene(gene, moi)

# Combine all of the enrichment data together.
combine_statistics_all_genes()

# Combine all of the plots of point mutations estimated by the bayesian model.
combine_predicted_plots(date_all)
