library(tidyr)
library(lubridate)
library(stringr)
library(dplyr)
library(data.table)

#
#	file names used below
#	these are defined separately so they can be overwritten as needed for PID code
#

GFNAME_jhu_death_data_padded <<- "usa/data/jhu_death_data_padded.rds"
GFNAME_nyt_death_data_padded <<- "usa/data/nyt_death_data_padded.rds"
GFNAME_states <<- "usa/data/states.csv"
GFNAME_weighted_ifr <<- "usa/data/weighted_ifr_states.RDS"
GFNAME_global_mobility_report <<- 'usa/data/Global_Mobility_Report.csv'
GFNAME_us_population <<- "usa/data/us_population.rds"
GFNAME_us_regions <<- "usa/data/usa-regions.csv"

#
#	functions
#

smooth_fn = function(x, days=3) {
  return(rollmean(x, days, align = "right"))

}

read_death_data <- function(source, smooth = FALSE){
  if (source == "jhu"){
    death_data <- readRDS(GFNAME_jhu_death_data_padded)
  } else if (source ==  "nyt"){
    death_data <- readRDS(GFNAME_nyt_death_data_padded)
  }
  death_data <- death_data[which(death_data$date >= ymd("2020-02-1")),]
  
  if (smooth == TRUE){
    k = 3
    death_data$daily_deaths = c(rep(0, k-1), as.integer(smooth_fn(death_data$daily_deaths, days=k)))
  }
  
  #read regions
  regions <- read.csv('usa/data/usa-regions.csv', stringsAsFactors = FALSE)
  death_data <- merge(death_data, regions, by = 'code')
  return(death_data)
}

read_ifr_data <- function(){
  state_ifr <- readRDS(GFNAME_weighted_ifr)
  return(state_ifr)
  
}


read_pop_count_us = function(path_to_file = GFNAME_us_population){
  pop_count <- readRDS(path_to_file)
  pop_count <- pop_count[which(!is.na(pop_count$code)),] %>%
    reshape2::melt(id.vars = c("Region", "code")) %>%
    rename(age = variable, pop = value, state = Region)
  pop_count <- pop_count %>%
    group_by(state) %>%
    summarise(popt:= sum(pop)) 

  return(pop_count)
}