require(data.table)
require(rstan)
require(EnvStats)
require(covid19AgeModel)
pkg.dir <- system.file(package = "covid19AgeModel" )

## command line parsing if any
args <- list(    
  indir.results= '/rds/general/project/ratmann_covid19/live/age_renewal_usa/base_age_fsq_mobility_201015e8_cmdstanv-4states_Oct21_Levin_K12open/base_age_fsq_mobility_201015e8_cmdstanv-4states_Oct21_Levin_K12open-2501614[1].pbs',
  location.index= 1,
  with.flow=1,
  forecast.period= 90, # days
  school.reopen= 1,
  school_level="K5",
  multiplier_cntct_school_opening = 1,
  counterfactual.scenario = 1
)


args_line <-  as.list(commandArgs(trailingOnly=TRUE))
if(length(args_line) > 0) 
{
  stopifnot(args_line[[1]]=='-indir.results')
  stopifnot(args_line[[3]]=='-location.index')
  stopifnot(args_line[[5]]=='-with.flow')
  stopifnot(args_line[[7]]=='-forecast.period')
  stopifnot(args_line[[9]]=='-school.reopen')
  stopifnot(args_line[[11]]=='-multiplier_cntct_school_opening')
  stopifnot(args_line[[13]]=='-school_level')
  stopifnot(args_line[[15]]=='-counterfactual.scenario')
  
  args <- list()    
  args[['indir.results']] <- args_line[[2]]  
  args[['location.index']] <- args_line[[4]]  
  args[['with.flow']] <- args_line[[6]]  
  args[['forecast.period']] <- args_line[[8]]  
  args[['school.reopen']] <- args_line[[10]] 
  args[['multiplier_cntct_school_opening']] <- args_line[[12]]
  args[['school_level']] <- args_line[[14]]
  args[['counterfactual.scenario']] <- args_line[[16]]
}

#	print args
str(args)

#	read args
indir.results <- args$indir.results
location.index <- as.integer(args$location.index)
with.flow <- args$with.flow
forecast.period <- as.numeric(args$forecast.period)
multiplier_cntct_school_opening <- as.numeric(args$multiplier_cntct_school_opening)
school.reopen <- args$school.reopen == '1'
school_level <- as.character(args$school_level)
counterfactual.scenario = as.numeric(args$counterfactual.scenario)

#	read Stan input data and add location.index to parallelise computations
cat('\nReading Stan input data...')
infile.stanin <- list.files(indir.results, pattern='.*_stanin.RData$', recursive=TRUE)
stopifnot(length(infile.stanin)>0)
stopifnot(length(infile.stanin)<=1)
tmp <- load(file.path(indir.results, infile.stanin))
stopifnot(c('args','stan_data')%in%tmp)
stan_data$LOCATION_PROCESSING_IDX <- location.index

# reset args
forecast <- forecast.period
args[['work_dir']] <- getwd()
pkg.dir <- system.file(package = "covid19AgeModel" )
args$file_contact_intensities_outbreak_China <- file.path(pkg.dir, "data-v110", "estimate_contact_intensities_outbreak_China.rds")
args$file_school_Ox_NPI_data <- file.path(pkg.dir, "data-v110", "OxCGRT_US_subnational_06Oct2020.csv") 
args$file_school_EduWeek_data <- file.path(pkg.dir, "data-v110", "Coronavirus_and_School_Closures_200910.csv") 

# dev for new args
if(is.null(args$with_EduWeek_data)) args$with_EduWeek_data = 0
if(is.null(args$with_ifr_fixed_time_decay)) args$with_ifr_fixed_time_decay = 0
if(is.null(args$effect_weeks)) args$effect_weeks <- 2

# determine which levels have closed/reopened
args$with_elementary_schools_status = 0
args$with_elementary_high_schools_status = 0
if(school_level == "K5"){
  args$with_elementary_schools_status = 1
}
if(school_level == "K12"){
  args$with_elementary_high_schools_status = 1
}

num_days_sim <- (max(death_data$date) - min(death_data$date) + 1 + forecast)[[1]]

# reprocess data
processed_data <- make_stan_data_core(states = args$states, 
                                      death_data = death_data, 
                                      deathByAge_data = deathByAge_data,                                           
                                      ifr.by.age = ifr.by.age,                                         
                                      serial_interval = serial_interval,
                                      pop_info = pop_info, 
                                      dcontact= dcontact,     
                                      num_days_sim = num_days_sim, 
                                      forecast = forecast)

if(args$with_google_mobility && 
   args$decouple_mobility_trends_into_baseline_plus_rest_parts==0 &&
   args$decouple_mobility_trends_into_decline_plus_upswing_parts==0 && 
   args$decouple_mobility_trends_into_baseline_plus_decline_plus_upswing_parts==0)
{
  cat('\nAdding google mobility covariates to stan_data ...')
  processed_data <- stan_data_add_google_mobility(processed_data, mobility_data, args$with_avg_mobility_data)
}
if(args$with_fsq_mobility && 
   args$decouple_mobility_trends_into_baseline_plus_rest_parts==0 &&
   args$decouple_mobility_trends_into_decline_plus_upswing_parts==0 && 
   args$decouple_mobility_trends_into_baseline_plus_decline_plus_upswing_parts==0)
{  
  cat('\nAdding Foursquare mobility covariates to stan_data ...')
  processed_data <- stan_data_add_mobility_trends_fsq(processed_data, mobility_data, log=TRUE)
}
if( args$with_google_mobility && 
    args$decouple_mobility_trends_into_decline_plus_upswing_parts==1 )
{	
  cat('\nAdding 2 part decoupled Google mobility trends to stan_data ... ')
  processed_data <- stan_data_add_google_mobility_covariates_after_breakpoint(processed_data)
}
if( args$with_fsq_mobility==1 && 
    args$decouple_mobility_trends_into_decline_plus_upswing_parts==1 )
{		
  cat('\nAdding decoupled FSQ mobility trends (decline + upswing) to stan_data ... ')
  processed_data <- stan_data_add_decoupled_mobility_trends_decline_plus_upswing_fsq(processed_data, mobility_data, log=TRUE)	
}
if( args$with_fsq_mobility==1 && 
    args$decouple_mobility_trends_into_baseline_plus_rest_parts==1 )
{		
  cat('\nAdding decoupled FSQ mobility trends (baseline + rest) to stan_data ... ')
  processed_data <- stan_data_add_decoupled_mobility_trends_decline_plus_rest_fsq(processed_data, mobility_data, log=TRUE)	
}
if( args$with_fsq_mobility==1 && 
    args$decouple_mobility_trends_into_baseline_plus_decline_plus_upswing_parts==1 )
{		
  cat('\nAdding decoupled FSQ mobility trends (baseline + decline + upswing) to stan_data ... ')
  processed_data <- stan_data_add_decoupled_mobility_trends_3_parts_fsq(processed_data, mobility_data, log=TRUE)	
}
if( args$with_estimated_pairwise_mobility_trends)
{	
  # process_chi_emodo_mobilty(mobility_data,pop_info,infile_emodo='~/contacts_and_mobility/contacts_by_age_20200729.csv')
  cat('\n Add pairwise mobility trend ... ')	
  processed_data <- stan_data_add_pairwise_mobility_trends(processed_data,file.path(pkg.dir,'data-v110','chi_emodo_mobility_trends.rds'))
}
if( args$with_biweekly_upswing_time_effect | 
    args$with_biweekly_upswing_time_effect_for_all_loc |
    args$with_ifr_fixed_time_decay
)
{
  cat('\nAdding upswing time effects ... ')
  processed_data <- stan_data_add_upswing_time_effect(
    processed_data, 
    mobility_data, 
    effect_weeks=args$effect_weeks)		 	
}
if( args$with_contact_intensities_zhang && !args$with_EduWeek_data )
{
  reopening_date <- as.Date("2020-08-24")
  cat('\nAdding school status ... ')
  processed_data <- stan_data_add_school_status(
    processed_data, 
    args$file_school_Ox_NPI_data, 
    school.reopen, 
    reopening_Date= reopening_date
  )    
}
if( args$with_contact_intensities_zhang && args$with_EduWeek_data)
{
  default_school_reopening_date <- as.Date("2020-08-24")
  cat('\nAdding school status with Education Weekly data ... ')
  processed_data <- stan_data_add_school_status_EduWeek(
    processed_data = processed_data, 
    path_to_file_school_Ox_data = args$file_school_Ox_NPI_data, 
    path_to_file_school_EduWeek_data = args$file_school_EduWeek_data,
    default_school_reopening_date = default_school_reopening_date,
    counterfactual_scenario = as.numeric(school.reopen == 0 | counterfactual.scenario == 1),
    forecast_with_schools_reopened = school.reopen
  )    
}
if( args$with_contact_intensities_zhang && args$with_elementary_schools_status ){
  cat('\nAdding contact intensities during elementary school closure from Zhang et al. and mutliplier... ')
  processed_data <- stan_data_add_contact_school_closure(
    processed_data, 
    dcontact, 
    args$file_contact_intensities_outbreak_China, 
    args$multiplier_cntct_school_closure
  )
}
if( args$with_contact_intensities_zhang && args$with_elementary_high_schools_status ){
  cat('\nAdding contact intensities during elementary and high school closure from Zhang et al. and mutliplier... ')
  processed_data <- stan_data_add_contact_elementary_high_school_closure(
    processed_data, 
    dcontact, 
    args$file_contact_intensities_outbreak_China, 
    args$multiplier_cntct_school_closure
  )
}
if( args$with_contact_intensities_zhang && !args$with_EduWeek_data && args$with_elementary_schools_status )
{
  cat('\nAdding contact intensities after school re-opening from Zhang et al. and mutliplier... ')
  processed_data <- stan_data_add_contact_school_opening(
    processed_data, 
    dcontact, 
    args$file_contact_intensities_outbreak_China, 
    args$multiplier_cntct_school_closure, 
    multiplier_cntct_school_opening=multiplier_cntct_school_opening, 
    min_pc_contacts=0.05,
    reopening_date= reopening_date
  )	  
}
if(args$with_contact_intensities_zhang && args$with_EduWeek_data && args$with_elementary_schools_status){
  cat('\nAdding contact intensities after school re-opening from Zhang et al. and mutliplier... ')
  processed_data <- stan_data_add_contact_school_opening_EduWeek(
    processed_data, 
    dcontact, 
    path_to_file_contact_intensities_outbreak_China = args$file_contact_intensities_outbreak_China, 
    multiplier_cntct_school_closure=args$multiplier_cntct_school_closure, 
    multiplier_cntct_school_opening=multiplier_cntct_school_opening, 
    min_pc_contacts=0.05
  )	  
}
if(args$with_contact_intensities_zhang && args$with_EduWeek_data && args$with_elementary_high_schools_status){
  cat('\nAdding contact intensities after elementary and high school re-opening from Zhang et al. and mutliplier... ')
  processed_data <- stan_data_add_contact_elementary_high_school_opening_EduWeek(
    processed_data, 
    dcontact, 
    args$file_contact_intensities_outbreak_China, 
    args$multiplier_cntct_school_closure, 
    multiplier_cntct_school_opening=multiplier_cntct_school_opening, 
    min_pc_contacts=0.05
  )	  
}

# reset stan_data
for (name in names(processed_data$stan_data)) {
  stan_data[[name]] <- processed_data$stan_data[[name]]
}

if(with.flow=='1')
{  
  stan_data <- gqs_add_stan_data_for_flows(stan_data,dates)
}

if(is.null(stan_data$rev_iar_daysSinceInfection))
{
  stan_data <- gqs_add_stan_data_for_iar(stan_data)
}

# if it is a counterfactual scenario, the school effect is the supplied school reopening multiplier 
if(counterfactual.scenario & school.reopen){
  stan_data[['counterfactual_school_effect']] = multiplier_cntct_school_opening
} else{ # else the school effect is elt_school_intv_effect estimated in the analysis
  stan_data[['counterfactual_school_effect']] = -1
}


#	read stanfits
cat('\nReading Stanfit ...')
infile.stanfits <- list.files(indir.results, pattern='.*_stanout.RData$', recursive=TRUE)
stopifnot(length(infile.stanfits)>0)
stopifnot(length(infile.stanfits)<=1)
tmp <- load(file.path(indir.results, infile.stanfits[1]))

#	find stan model and stan model for generating quantities
file_stanModel <- gsub('.*covid19AgeModel/(.*)','\\1',args$file_stanModel)
file_stanModel <- file.path(pkg.dir, file_stanModel)
file_stanModel_gqs <- gsub('cmdstanv','gqs',file_stanModel)
stopifnot( file.exists(file_stanModel_gqs) )

cat('\nCompiling gqs model file ...')
m2 <- rstan::stan_model(file_stanModel_gqs)

cat('\nGenerating quantities ...')
draws <- as.matrix(fit)
fit2 <- rstan::gqs(m2, data=stan_data, draws=draws)
fit.gqs <- rstan::extract(fit2)

#	save stan_data
if(location.index==1)
{
  multiplier_name <- (multiplier_cntct_school_opening)*100
  tmp <- file.path(indir.results, paste0(basename(args$job_dir),'_stan_data', '_sensitivity_school_reopen_', as.integer(school.reopen), '_counterfactual_', counterfactual.scenario,  '_multiplier_', multiplier_name, '_level_', school_level,'.RDS'))
  if(as.integer(school.reopen) == 0)
  {
    tmp <- file.path(indir.results, paste0(basename(args$job_dir),'_stan_data', '_sensitivity_school_reopen_', as.integer(school.reopen), '_level_', school_level, '.RDS'))
  }
  saveRDS(stan_data, file=tmp)	  
}

#	save fit.gqs, try up to 5 times
multiplier_name <- (multiplier_cntct_school_opening)*100
file <- file.path(indir.results, paste0(basename(args$job_dir), '_location', location.index,'_stangqs', '_sensitivity_school_reopen_', as.integer(school.reopen), '_counterfactual_', counterfactual.scenario, '_multiplier_', multiplier_name, '_level_', school_level, '.RDS'))
if(as.integer(school.reopen) == 0)
{
  file <- file.path(indir.results, paste0(basename(args$job_dir), '_location', location.index,'_stangqs', '_sensitivity_school_reopen_', as.integer(school.reopen), '_level_', school_level, '.RDS'))
}
io_saveRDS(fit.gqs, args[['work_dir']], dirname(file), basename(file), check_if_saved_n=10)

cat('\nFinished base-ages-generate-quantities-forecast.r ...')