library("abc")
library(bayestestR)
# setwd("") maybe will be needed in order to run the project from the folder
rm(list = ls())
glob_stat <<- 0 # It will count the number of datasets I have traversed

# Function that calculates the w statistic
W_Statistic <- function() {
  out <- 0
  for (i in 1:49)
  {
    out <- out + 1 / i
  }
  return(out)
}

# Function that calculates the k statistic
K_Statistic <- function(file, wt) {
  k_stat <- 0
  i <- 1
  dataset <- scan(file, what = "", sep = ",", quote = "\"", )

  # For each line in the dataset(ka8e dataset xwrizetai me kenh gramh ara gia ka8e ena apo auta) I split it to each character
  dataset <- scan(file, what = "", sep = ",", quote = "\"", )
  for (i in 1:(length(dataset) - 1)) {
    for (j in (i + 1):length(dataset)) {
      dif_temp <- mapply(function(x, y) sum(x != y), strsplit(dataset[i], ""), strsplit(dataset[j], ""))
      k_stat <- k_stat + dif_temp
    }
  }

  # I write each k I found into the k_stat that I will delete afterwards
  write.table(k_stat / length(dataset), file = wt, row.names = FALSE, col.names = FALSE, append = TRUE, quote = F, eol = "\n")
  return(k_stat / length(dataset))
}

# Function that Combines K file and W file into 1 Final statistic file
combine_stats <- function(filek, filew) {
  k_stat <- rep(NA, 10000)
  w_stat <- rep(NA, 10000)
  i <- 1
  conk <- file(filek, "r")
  conw <- file(filew, "r")

  print("Combining...")

  for (i in 1:10000) {
    linek <- readLines(conk, n = 1)
    linew <- readLines(conw, n = 1)

    k_stat[i] <- linek
    w_stat[i] <- linew
  }

  # The df contains all statistics of k,w so I simply write that into ms.stats
  df <- data.frame(
    s_k = k_stat,
    s_w = w_stat
  )
  write.table(df, "ms.stats.txt", row.names = FALSE, quote = F)

  close(conk)
  close(conw)
}

# Function that splits 10000 datasets into each dataset and then passes it to the appropriate functions
processFile <- function(filepath, sp) {
  con <- file(filepath, "r")
  i <- 1
  while (i < sp) {
    line <- readLines(con, n = 1)
    i <- i + 1
  }

  line <- readLines(con, n = 1)

  # Tmp file will have 1 dataset each time that will be transfered to k and w functions
  write.table(line, file = "tmp.txt", row.names = FALSE, col.names = FALSE, quote = F, eol = "\n")

  reposition_sp <- 0
  j <- 0
  col <- 0
  if (length(line) != 0) {
    while (j < 50) {
      line <- readLines(con, n = 1)
      write.table(line, file = "tmp.txt", row.names = FALSE, col.names = FALSE, quote = F, append = TRUE, eol = "\n")
      reposition_sp <- reposition_sp + 1
      if (j == 0) {
        col <- nchar(line)
      }
      j <- j + 1
    }
  }

  data <- file("tmp.txt")
  k_stat <- K_Statistic(data, "k_stat.txt")
  w_stat <- col / W_Statistic()
  print(paste0("For dataset number: ", glob_stat, " K stat is: ", k_stat, " and W stat is: ", w_stat))
  write.table(w_stat, file = "w_stat.txt", row.names = FALSE, col.names = FALSE, quote = F, append = TRUE, eol = "\n")
  close(data)

  # If I found EOF I must stop searching the file(end of 10000 datasets)
  line <- readLines(con, n = 1)
  if (length(line) == 0) {
    print("EOF")
    glob_stat <<- glob_stat + 1
    return("Zero")
  }

  close(con)

  sp_offset <- sp + reposition_sp + 1
  # glob_stat counts the dataset
  glob_stat <<- glob_stat + 1
  return(sp_offset)
}

# Function that will calculate the s,w statistic of obs file.
Observe <- function(file) {
  data <- file(file)
  dataset <- scan(file, what = "", sep = ",", quote = "\"", )
  k <- K_Statistic(data, "obs.txt")
  w <- nchar(dataset[1]) / W_Statistic()

  # This data frame(df) will cointain all the  statistics from the observation that I have previously calculated.So I pass them into the data frame and then the data frame to the obs.txt
  df <- data.frame(
    s_k = k,
    s_w = w
  )
  write.table(df, "obs.txt", row.names = FALSE, quote = F)

  close(data)
}

# Function that using a while loop finds the statistics and the combines them into the final txt
CalculateStats <- function() {
  output <- processFile("ms_sim_final.out", 1)
  while (glob_stat < 10000) {
    output <- processFile("ms_sim_final.out", output)
  }

  combine_stats("k_stat.txt", "w_stat.txt")
}

# Function that solves the 1,3 question
SumUp <- function(obs, stats, param) {
  observation <- read.table(obs, header = TRUE)
  observation
  sims <- read.table(stats, header = TRUE)
  pars <- read.table(param, header = FALSE)

  simstats <- sims
  names(simstats)

  dim(observation)
  dim(simstats)

  myabc <- abc(target = observation, param = pars, sumstat = simstats, method = "loclinear", tol = 0.01)
  summary(myabc)

  ci(myabc$adj.values, method = "HDI")

  pdf("plots.pdf")
  plot(density(pars[, 1]))

  d.prior <- density(pars[, 1])
  d.pos1 <- density(myabc$unadj.values[, 1])
  d.pos <- density(myabc$adj.values[, 1])

  plot(d.prior$x, d.prior$y, ylim = c(0, max(d.prior$y, d.pos1$y, d.pos$y)), type = "l", col = "black")
  points(d.pos1$x, d.pos1$y, type = "l", col = "blue")
  points(d.pos$x, d.pos$y, type = "l", col = "red")

  dev.off()
}

# Function that solves the 2 question
Ce <- function(obs, stats, param) {
  observation <- read.table(obs, header = TRUE)
  observation
  sims <- read.table(stats, header = TRUE)
  pars <- read.table(param, header = FALSE)

  simstats <- sims
  names(simstats)

  dim(observation)
  dim(simstats)

  myabc <- abc(target = observation, param = pars, sumstat = simstats, method = "loclinear", tol = 0.01)

  ci(myabc$adj.values, method = "HDI")
}

# Function that removes files that I dont need
CleanFiles <- function() {
  file.remove("tmp.txt")
}

# To Calculate obs statistics run
# Observe("ms_obs_final.out.txt")

# To calculate the k,w stats run
# CalculateStats()

# To find the answers for the questions 1,3 and make the pdf with the plots run
# SumUp("obs.txt", "ms.stats.txt", "pars_final.txt")

# To find the answer for the question 2(95% HDI: [78.27, 149.62]) run
# Ce("obs.txt", "ms.stats.txt", "pars_final.txt")

# To clean the directory of tmp files simply run
# CleanFiles()
