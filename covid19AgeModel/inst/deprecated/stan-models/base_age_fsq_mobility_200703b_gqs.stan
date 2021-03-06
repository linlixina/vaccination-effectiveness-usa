functions {


   // function equivalent to %in% on R from https://discourse.mc-stan.org/t/stan-equivalent-of-rs-in-function/3849
  int r_in(int pos,int[] pos_var) {
   
    for (p in 1:(size(pos_var))) {
       if (pos_var[p]==pos) {
       // can return immediately, as soon as find a match
          return 1;
       } 
    }
    return 0;
  }

  /* 
   * returns multiplier on the rows of the contact matrix over time for one country
   */
  matrix country_impact(
      vector beta,
      real upswing_rdeff_local,
      int N2,
      int A,
      int COVARIATES_N,
      real avg_cntct_local,
      matrix[] covariates_local
      )
  {

    // scaling of contacts after intervention effect on day t in location m
    matrix[N2,A] impact_intv;
    
    // define multipliers for contacts in each location
    impact_intv = rep_matrix( 0, N2, A);
    //for(i in 1:COVARIATES_N)
    //{
    impact_intv += beta[1] * covariates_local[1] + (beta[2] + upswing_rdeff_local) * covariates_local[2];
    //}
    impact_intv = exp( impact_intv );
        
    return (impact_intv);
  }
  
    matrix country_EcasesByAge(// parameters
                              real R0_local,
                              real e_cases_N0_local,
                              vector beta,
                              real upswing_rdeff_local,
                              row_vector log_relsusceptibility_age,
                              row_vector log_reltransmissibility_age,
                              // data
                              int N0,
                              int N2,
                              int A,
                              int COVARIATES_N,
                              int SI_CUT,
                              // wkend_idx[1:WKEND_IDX_N[m],m]
                              int[] wkend_idx_local,
                              real avg_cntct_local,
                              matrix[] covariates_local,
                              matrix cntct_weekends_mean_local,
                              matrix cntct_weekdays_mean_local,
                              row_vector rev_serial_interval,
                              row_vector popByAge_abs_local,
                              int init_A
                              )
  {
    // dummy variables
    // real<lower=0> saturation_adj;
    real zero = 0.0;
    
    // probability of infection given contact in location m
    real rho0 = R0_local / avg_cntct_local;
    
    // expected new cases by calendar day, age, and location under self-renewal model
    // and a container to store the precomputed cases by age
    matrix[N2,A] E_casesByAge = rep_matrix( zero, N2, A );

    matrix[N2,A] impact_intv = country_impact(beta,
                                               upswing_rdeff_local,
                                               N2,
                                               A,
                                               COVARIATES_N,
                                               avg_cntct_local,
                                               covariates_local);
    
    // init expected cases by age and location in first N0 days
    E_casesByAge[1:N0,init_A] = rep_vector( e_cases_N0_local, N0 );
  
    // calculate expected cases by age and country under self-renewal model after first N0 days
    // and adjusted for saturation
    for (t in (N0+1):N2)
    {
      int start_idx_rev_serial = SI_CUT-t+2;
      int start_idx_E_casesByAge = t-SI_CUT;
      row_vector[A] prop_susceptibleByAge = rep_row_vector(1.0, A) - (rep_row_vector(1.0, t-1) * E_casesByAge[1:(t-1),:] ./ popByAge_abs_local);
      if(start_idx_rev_serial < 1) {
        start_idx_rev_serial = 1;
      }
      if(start_idx_E_casesByAge < 1) {
        start_idx_E_casesByAge = 1;
      }
      // TODO can t we vectorise this?
      for(a in 1:A){
        if(prop_susceptibleByAge[a] < 0) { // account for values of Ecases > pop at initalization
          prop_susceptibleByAge[a] = 0;
        }
      }
      {
        row_vector[A] tmp_row_vector_A = rev_serial_interval[start_idx_rev_serial:SI_CUT] * E_casesByAge[start_idx_E_casesByAge:(t-1)];
        tmp_row_vector_A .*= impact_intv[t,];
        tmp_row_vector_A .*= ( rho0 * exp(log_reltransmissibility_age) );
        if(r_in(t, wkend_idx_local) == 1){
            E_casesByAge[t] = tmp_row_vector_A * cntct_weekends_mean_local;
        }else{
            E_casesByAge[t] = tmp_row_vector_A * cntct_weekdays_mean_local;
        }
        E_casesByAge[t] .*= prop_susceptibleByAge;
        E_casesByAge[t] .*= exp(log_relsusceptibility_age);
        //E_casesByAge[t] .*= impact_intv[t,];
      }
    }
  
    return(E_casesByAge);
  }
  
  
  matrix country_Rta(// parameters
      real rho0_local,
      row_vector log_relsusceptibility_age,
      row_vector log_reltransmissibility_age,
      matrix impact_intv,
      matrix E_casesByAge_local,
      // data
      int N2,
      int A,
      int[] wkend_idx_local,
      matrix cntct_weekends_mean_local,
      matrix cntct_weekdays_mean_local,
      row_vector popByAge_abs_local
      )
  {
        matrix[N2,A] RtByAge;
        matrix[N2,A] prop_susceptibleByAge;
        
        for(a in 1:A)
        {
            prop_susceptibleByAge[1,a] = 0;
            prop_susceptibleByAge[2:N2,a] = cumulative_sum( E_casesByAge_local[1:(N2-1),a] ) / popByAge_abs_local[a];
            for( t in 1:N2)
            {
                if( prop_susceptibleByAge[t,a]>1 )
                {
                    prop_susceptibleByAge[t,a]= 1.;
                }
            }
        }
        prop_susceptibleByAge = rep_matrix(1.0, N2, A) - prop_susceptibleByAge;
        RtByAge = prop_susceptibleByAge;
        RtByAge .*= rep_matrix( exp(log_relsusceptibility_age), N2);
        RtByAge *= rho0_local;
        for(t in 1:N2)
        {
            if(r_in(t, wkend_idx_local) == 1)
            {
                // RtByAge[t,:] = ( (impact_intv[t,:] .* RtByAge[t,:]) * (cntct_weekends_mean[m]') ) .* impact_intv[t,:];
                RtByAge[t,:] = ( RtByAge[t,:] * (cntct_weekends_mean_local') ) .* impact_intv[t,:];
            }
            else
            {
                // RtByAge[t,:] = ( (impact_intv[t,:] .* RtByAge[t,:]) * (cntct_weekdays_mean[m]') ) .* impact_intv[t,:];
                RtByAge[t,:] = ( RtByAge[t,:] * (cntct_weekdays_mean_local') ) .* impact_intv[t,:];
            }
        }
        return(RtByAge);
  }
  
    matrix country_lambdaByAge(// parameters
    real R0_local,
    vector beta,
    real upswing_rdeff_local,
    row_vector log_relsusceptibility_age,
    row_vector log_reltransmissibility_age,
    // data
    int N0,
    int N2,
    int A,
    int COVARIATES_N,
    int SI_CUT,
    // wkend_idx[1:WKEND_IDX_N[m],m]
    int[] wkend_idx_local,
    real avg_cntct_local,
    matrix[] covariates_local,
    matrix cntct_weekends_mean_local,
    matrix cntct_weekdays_mean_local,
    row_vector rev_serial_interval,
    matrix E_casesByAge
    )
    {
      // dummy variables
      // real<lower=0> saturation_adj;
     real zero = 0.0;
     int tflow =0;
     matrix[N2,A] impact_intv = country_impact(beta,
                                               upswing_rdeff_local,
                                               N2,
                                               A,
                                               COVARIATES_N,
                                               avg_cntct_local,
                                               covariates_local);
      
      row_vector[A] tmp_row_vector_A;
      matrix[A,A] tmp_lambda;
      matrix[N2, A] lambdaByAge;
      
      // probability of infection given contact in location m
      real rho0 = R0_local / avg_cntct_local;
    
      // calculate expected cases by age and country under self-renewal model after first N0 days
      // and adjusted for saturation
      lambdaByAge[1:N0, 1:A] = rep_matrix(0.,N0,A);
      
      for (t in (N0+1):N2)
      {
        int start_idx_rev_serial = SI_CUT-t+2;
        int start_idx_E_casesByAge = t-SI_CUT;
        
          if(start_idx_rev_serial < 1) {
            start_idx_rev_serial = 1;
          }
          if(start_idx_E_casesByAge < 1) {
            start_idx_E_casesByAge = 1;
          }

            tmp_row_vector_A = rev_serial_interval[start_idx_rev_serial:SI_CUT] * E_casesByAge[start_idx_E_casesByAge:(t-1)];
            tmp_row_vector_A .*= impact_intv[t,];
            tmp_row_vector_A .*= ( rho0 * exp(log_reltransmissibility_age) );
            
            if(r_in(t, wkend_idx_local) == 1){
              lambdaByAge[t] = tmp_row_vector_A * cntct_weekends_mean_local;
            }else{
              lambdaByAge[t] = tmp_row_vector_A * cntct_weekdays_mean_local;
            }
            
            lambdaByAge[t] .*= exp(log_relsusceptibility_age);
            //tmp_lambda .*= rep_matrix( impact_intv[t,],A);
            
      }
      return(lambdaByAge);
    }
    
    matrix[] country_EflowsByLowDimAge(// parameters
        real R0_local,
        real e_cases_N0_local,
        vector beta,
        real upswing_rdeff_local,
        row_vector log_relsusceptibility_age,
        row_vector log_reltransmissibility_age,
        // data
        int N0,
        int N2,
        int A,
        int COVARIATES_N,
        int SI_CUT,
        // wkend_idx[1:WKEND_IDX_N[m],m]
        int[] wkend_idx_local,
        real avg_cntct_local,
        matrix[] covariates_local,
        matrix cntct_weekends_mean_local,
        matrix cntct_weekdays_mean_local,
        row_vector rev_serial_interval,
        row_vector popByAge_abs_local,
        int init_A,
        int Nflow,
        int[] N_low_dim_age_band,
        int nrow_low_dim_age_band,
        int[,] low_dim_age_band,
        matrix E_casesByAge,
        int[] week_start_index
        )
    {
      // dummy variables
      // real<lower=0> saturation_adj;
     real zero = 0.0;
     int tflow=0;
     matrix[N2,A] impact_intv = country_impact(beta,
                                               upswing_rdeff_local,
                                               N2,
                                               A,
                                               COVARIATES_N,
                                               avg_cntct_local,
                                               covariates_local);
      
      row_vector[A] tmp_row_vector_A;
      matrix[A,A] tmp_flow;
      matrix[nrow_low_dim_age_band,nrow_low_dim_age_band] flow[Nflow];
      
      // probability of infection given contact in location m
      real rho0 = R0_local / avg_cntct_local;
    
      for (i in 1:Nflow){
        flow[i]=rep_matrix( zero,nrow_low_dim_age_band,nrow_low_dim_age_band );
      }
      
      // calculate expected cases by age and country under self-renewal model after first N0 days
      // and adjusted for saturation
      for (t in (N0+1):N2)
      {
        int start_idx_rev_serial = SI_CUT-t+2;
        int start_idx_E_casesByAge = t-SI_CUT;
        row_vector[A] prop_susceptibleByAge = rep_row_vector(1.0, A) - (rep_row_vector(1.0, t-1) * E_casesByAge[1:(t-1),:] ./ popByAge_abs_local);

          if(start_idx_rev_serial < 1) {
            start_idx_rev_serial = 1;
          }
          if(start_idx_E_casesByAge < 1) {
            start_idx_E_casesByAge = 1;
          }
          for(a in 1:A){
            if(prop_susceptibleByAge[a] < 0) { // account for values of Ecases > pop at initalization
                 prop_susceptibleByAge[a] = 0;
            }
          }

          if (t >= week_start_index[1]){
            if(r_in(t,week_start_index)){
              tflow += 1;
            }
            tmp_row_vector_A = rev_serial_interval[start_idx_rev_serial:SI_CUT] * E_casesByAge[start_idx_E_casesByAge:(t-1)];
            tmp_row_vector_A .*= impact_intv[t,];
            tmp_row_vector_A .*= ( rho0 * exp(log_reltransmissibility_age) );
            
            if(r_in(t, wkend_idx_local) == 1){
              tmp_flow = rep_matrix(tmp_row_vector_A', A) .* cntct_weekends_mean_local;
            }else{
              tmp_flow = rep_matrix(tmp_row_vector_A', A) .* cntct_weekdays_mean_local;
            }
            
            tmp_flow .*= rep_matrix(exp(log_relsusceptibility_age),A);
            tmp_flow .*= rep_matrix(prop_susceptibleByAge,A);
            //tmp_flow .*= rep_matrix( impact_intv[t,],A);
  
            for (a in 1:nrow_low_dim_age_band){
              for (aa in 1:nrow_low_dim_age_band){
                flow[tflow][a,aa] += sum( to_array_1d(tmp_flow[ low_dim_age_band[a,1:N_low_dim_age_band[a]],low_dim_age_band[aa,1:N_low_dim_age_band[aa]] ]));
              }
            }
          }
      }
      return(flow);
    }
  
    matrix[] country_EflowsByHighDimAge(// parameters
        real R0_local,
        real e_cases_N0_local,
        vector beta,
        real upswing_rdeff_local,
        row_vector log_relsusceptibility_age,
        row_vector log_reltransmissibility_age,
        // data
        int N0,
        int N2,
        int A,
        int COVARIATES_N,
        int SI_CUT,
        // wkend_idx[1:WKEND_IDX_N[m],m]
        int[] wkend_idx_local,
        real avg_cntct_local,
        matrix[] covariates_local,
        matrix cntct_weekends_mean_local,
        matrix cntct_weekdays_mean_local,
        row_vector rev_serial_interval,
        row_vector popByAge_abs_local,
        int init_A,
        int full_flows_on_week_n,
        int[] full_flows_on_week,
        int[] full_flows_on_day,
        matrix E_casesByAge
        )
    {
      // dummy variables
      // real<lower=0> saturation_adj;
     real zero = 0.0;
     int tflow = 0;
     matrix[N2,A] impact_intv = country_impact(beta,
                                               upswing_rdeff_local,
                                               N2,
                                               A,
                                               COVARIATES_N,
                                               avg_cntct_local,
                                               covariates_local);
      
      row_vector[A] tmp_row_vector_A;
      matrix[A,A] tmp_flow;
      matrix[A,A] full_flow[full_flows_on_week_n];
      
      // probability of infection given contact in location m
      real rho0 = R0_local / avg_cntct_local;
      
      for (i in 1:full_flows_on_week_n){
        full_flow[i]=rep_matrix( zero, A, A);
      }
      
      // calculate expected cases by age and country under self-renewal model after first N0 days
      // and adjusted for saturation
      for (t in (N0+1):N2)
      {
        int start_idx_rev_serial = SI_CUT-t+2;
        int start_idx_E_casesByAge = t-SI_CUT;
        row_vector[A] prop_susceptibleByAge = rep_row_vector(1.0, A) - (rep_row_vector(1.0, t-1) * E_casesByAge[1:(t-1),:] ./ popByAge_abs_local);
        
          if(start_idx_rev_serial < 1) {
            start_idx_rev_serial = 1;
          }
          if(start_idx_E_casesByAge < 1) {
            start_idx_E_casesByAge = 1;
          }
          for(a in 1:A){
            if(prop_susceptibleByAge[a] < 0) { // account for values of Ecases > pop at initalization
                 prop_susceptibleByAge[a] = 0;
            }
          }

          if (r_in(t,full_flows_on_day)){
            if(r_in(t,full_flows_on_week)){
              tflow += 1;
            }
            tmp_row_vector_A = rev_serial_interval[start_idx_rev_serial:SI_CUT] * E_casesByAge[start_idx_E_casesByAge:(t-1)];
            tmp_row_vector_A .*= impact_intv[t,];
            tmp_row_vector_A .*= ( rho0 * exp(log_reltransmissibility_age) );
            
            if(r_in(t, wkend_idx_local) == 1){
              tmp_flow = rep_matrix(tmp_row_vector_A', A) .* cntct_weekends_mean_local;
            }else{
              tmp_flow = rep_matrix(tmp_row_vector_A', A) .* cntct_weekdays_mean_local;
            }
            
            tmp_flow .*= rep_matrix(exp(log_relsusceptibility_age),A);
            tmp_flow .*= rep_matrix(prop_susceptibleByAge,A);
            //tmp_flow .*= rep_matrix(impact_intv[t,],A);
            
            full_flow[tflow] += tmp_flow;
              
          }
      }
      return(full_flow);
    }
  
  matrix check_country_EflowsByLowDimAge(// parameters
                                int nrow_low_dim_age_band,
                                int Nflow,
                                int N0,
                                int N2,
                                matrix E_casesByAge,
                                matrix[] flow,
                                int[] N_low_dim_age_band,
                                int[,] low_dim_age_band,
                                int[] week_start_index
                              )
  {
    real zero = 0.0;
    matrix[Nflow,nrow_low_dim_age_band] E_casesByLowDimAge=rep_matrix(zero,Nflow,nrow_low_dim_age_band );
    matrix[Nflow,nrow_low_dim_age_band] flow_cases;
    matrix[Nflow,nrow_low_dim_age_band] differences;
    int tflow = 0;
    
    for (t in (N0+1):N2)
    {
      if(t >= week_start_index[1]){
        if(r_in(t,week_start_index)){
          tflow += 1;
        }
        for (a in 1:nrow_low_dim_age_band){
          E_casesByLowDimAge[tflow,a] += sum(to_array_1d(E_casesByAge[t,low_dim_age_band[a,1:N_low_dim_age_band[a]]]));
        }
      }
    }
    
    for (t in 1:Nflow){
          for (a in 1:nrow_low_dim_age_band){
            flow_cases[t,a] = sum(to_array_1d((flow[t][,a])));
          }
    }
    
    differences = E_casesByLowDimAge - flow_cases;
    return(differences);
  }
  
  matrix country_EdeathsByAge(// parameters
                              real R0_local,
                              real e_cases_N0_local,
                              vector beta,
                              real upswing_rdeff_local,
                              row_vector log_relsusceptibility_age,
                              row_vector log_reltransmissibility_age,
                              real phi,
                              // data
                              int N0,
                              int N2,
                              int A,
                              int COVARIATES_N,
                              int SI_CUT,
                              int[] wkend_idx_local,
                              real avg_cntct_local,
                              matrix[] covariates_local,
                              matrix cntct_weekends_mean_local,
                              matrix cntct_weekdays_mean_local,
                              row_vector rev_ifr_daysSinceInfection,
                              row_vector log_ifr_age_base,
                              real log_ifr_age_rnde_mid1_local,
                              real log_ifr_age_rnde_mid2_local,
                              real log_ifr_age_rnde_old_local,
                              row_vector rev_serial_interval,
                              row_vector popByAge_abs_local,
                              int init_A
                              )
  {
    // dummy variables
    // real<lower=0> saturation_adj;

    real zero = 0.0;
    
    // expected deaths by calendar day (1st dim) age (2nd dim) and country (3rd dim), under self-renewal model
    // and expected deaths by calendar day (rows) and country (cols), under self-renewal model
    //matrix<lower=0>[N2,A] E_deathsByAge[M];
    matrix[N2,A] E_deathsByAge = rep_matrix( zero, N2, A );
  
    // expected new cases by calendar day, age, and location under self-renewal model
    // and a container to store the precomputed cases by age
    matrix[N2,A] E_casesByAge = country_EcasesByAge(R0_local,
                                    e_cases_N0_local,
                                    beta,
                                    upswing_rdeff_local,
                                    log_relsusceptibility_age,
                                    log_reltransmissibility_age,
                                    N0,
                                    N2,
                                    A,
                                    COVARIATES_N,
                                    SI_CUT,
                                    wkend_idx_local,
                                    avg_cntct_local,
                                    covariates_local,
                                    cntct_weekends_mean_local,
                                    cntct_weekdays_mean_local,
                                    rev_serial_interval,
                                    popByAge_abs_local,
                                    init_A);

    // calculate expected deaths by age and country
    E_deathsByAge[1] = 1e-15 * E_casesByAge[1];
    for (t in 2:N2)
    {
      E_deathsByAge[t] = rev_ifr_daysSinceInfection[(N2-(t-1)+1):N2 ] * E_casesByAge[1:(t-1)];
    }
    E_deathsByAge .*= rep_matrix(exp(   log_ifr_age_base +
            append_col(append_col(append_col(
                rep_row_vector(0., 4),
                rep_row_vector(log_ifr_age_rnde_mid1_local, 6)),
                rep_row_vector(log_ifr_age_rnde_mid2_local, 4)),
                rep_row_vector(log_ifr_age_rnde_old_local, 4))
            ), N2);
  
    return(E_deathsByAge);
  }
}

data {
  int<lower=1> M; // number of countries
  int<lower=1> N0; // number of initial days for which to estimate infections
  int<lower=1> N[M]; // days of observed data for country m. each entry must be <= N2
  int<lower=1> N2; // days of observed data + # of days to forecast
  int<lower=1> A; // number of age bands
  int<lower=1> SI_CUT; // number of days in serial interval to consider
  int<lower=1> COVARIATES_N; // number of days in serial interval to consider
  int WKEND_IDX_N[M]; // number of weekend indices in each location
  //	data
  real pop[M];
  matrix<lower=0, upper=1>[A,M] popByAge; // proportion of age bracket in population in location
  int epidemicStart[M];
  int deaths[N2, M]; // reported deaths -- the rows with i > N contain -1 and should be ignored
  int<lower=0> wkend_idx[N2,M]; //indices of 1:N2 that correspond to weekends in location m
  matrix[N2,A] covariates[M, COVARIATES_N]; // predictors for fsq contacts by age
  // data by age
  int<lower=0> M_AD; // number of countries with deaths by age data
  int<lower=1> dataByAgestart[M_AD]; // start of death by age data
  int deathsByAge[N2, A, M_AD]; // reported deaths by age -- the rows with i < dataByAgestart[M_AD] contain -1 and should be ignored + the column with j > A2[M_AD] contain -1 and should be ignored 
  int<lower=2> A_AD[M_AD]; // number of age groups reported 
  matrix[A, A] map_age[M_AD]; // map the age groups reported with 5 y age group -- the column with j > A2[M_AD] contain -1 and should be ignored
  int map_country[M,2]; // first column indicates if country has death by age date (1 if yes), 2 column map the country to M_AD
  //	priors
  matrix[A,A] cntct_weekdays_mean[M]; // mean of prior contact rates between age groups on weekdays
  matrix[A,A] cntct_weekends_mean[M]; // mean of prior contact rates between age groups on weekends
  real<upper=0> hyperpara_ifr_age_lnmu[A];  // hyper-parameters for probability of death in age band a log normal mean
  real<lower=0> hyperpara_ifr_age_lnsd[A];  // hyper-parameters for probability of death in age band a log normal sd
  row_vector[N2] rev_ifr_daysSinceInfection; // probability of death s days after infection in reverse order
  row_vector[SI_CUT] rev_serial_interval; // fixed pre-calculated serial interval using empirical data from Neil in reverse order
  int<lower=1, upper=A> init_A; // age band in which initial cases occur in the first N0 days
  //
  int<lower=1, upper=M> LOCATION_PROCESSING_IDX;
  int nrow_low_dim_age_band;
  int ncol_low_dim_age_band;
  int N_low_dim_age_band[nrow_low_dim_age_band];
  int low_dim_age_band[nrow_low_dim_age_band,ncol_low_dim_age_band];
  int debug_flow;
  int nrow_week_start_index;
  int week_start_index[nrow_week_start_index,M];
  int Nflow_country[M];
  int nrow_full_flows_on_week_country;
  int full_flows_on_week_country[nrow_full_flows_on_week_country,M];
  int full_flows_on_week_n_country[M];
  int nrow_full_flows_on_day_country;
  int full_flows_on_day_country[nrow_full_flows_on_day_country,M];
}

transformed data {
  vector<lower=0>[M] avg_cntct;
  vector[A] ones_vector_A = rep_vector(1.0, A);
  row_vector[A] ones_row_vector_A = rep_row_vector(1.0, A);
  int trans_deaths[M, N2]; // reported deaths -- the rows with i > N contain -1 and should be ignored
  matrix[M,A] popByAge_abs; 
  
  for( m in 1:M )
  {
    avg_cntct[m] = popByAge[:,m]' * ( cntct_weekdays_mean[m] * ones_vector_A ) * 5./7.;
    avg_cntct[m] += popByAge[:,m]' * ( cntct_weekends_mean[m] * ones_vector_A ) * 2./7.;

    trans_deaths[m,:] = deaths[:,m];
    
    popByAge_abs[m,:] = popByAge[:,m]' * pop[m]; // pop by age is a proportion of pop by age and pop is the absolute number 
  }
}

parameters {
  vector<lower=0>[M] R0; // R0
  real<lower=0> kappa; // variance parameter for country-specific R0  
  real<lower=0> tau; // prior rate of expected number of cases per day in the first N0 days, for each country
  real<lower=0> e_cases_N0[M]; // expected number of cases per day in the first N0 days, for each country
  vector[COVARIATES_N] beta; // regression coefficients for time varying multipliers on contacts
  real upswing_rnde[M];
  real<lower=0> sd_upswing_rnde;
  real<lower=0> phi; // overdispersion parameter for likelihood model
  row_vector<upper=0>[A] log_ifr_age_base; // probability of death for age band a
  row_vector[M] log_ifr_age_rnde_mid1;
  row_vector[M] log_ifr_age_rnde_mid2;
  row_vector[M] log_ifr_age_rnde_old;
  real<lower=0> sd_log_ifr_age_rnde_mid1;
  real<lower=0> sd_log_ifr_age_rnde_mid2;
  real<lower=0> sd_log_ifr_age_rnde_old;
  row_vector[2] log_relsusceptibility_age_reduced;
  row_vector[2] log_reltransmissibility_age_reduced;
  real<lower=0> sd_log_reltransmissibility_age;
}

generated quantities
{
    real rho0;
    matrix<lower=0>[N2,A] E_casesByAge;
    matrix<lower=0>[N2,A] E_deathsByAge;
    vector<lower=0>[N2] E_deaths;
    matrix<lower=0>[N2,A] RtByAge;
    matrix<lower=0>[N2,A] lambdaByAge;
    vector<lower=0>[N2] Rt;
    matrix[Nflow_country[LOCATION_PROCESSING_IDX],nrow_low_dim_age_band] differences; 
    matrix[nrow_low_dim_age_band,nrow_low_dim_age_band] flow[Nflow_country[LOCATION_PROCESSING_IDX]];
    matrix[A,A] full_flow[full_flows_on_week_n_country[LOCATION_PROCESSING_IDX]];
    row_vector[A] log_relsusceptibility_age;
    row_vector[A] log_reltransmissibility_age;

    log_relsusceptibility_age = append_col( append_col( log_relsusceptibility_age_reduced[ { 1, 1, 1 } ],
        rep_row_vector(0., 10) ),
        log_relsusceptibility_age_reduced[ { 2,2,2,2,2 } ]
        );
    log_reltransmissibility_age = append_col( append_col( log_reltransmissibility_age_reduced[ { 1, 1, 1 } ],
        rep_row_vector(0., 10) ),
        log_reltransmissibility_age_reduced[ { 2,2,2,2,2 } ]
        );
    // generate expected cases by age + expected deaths by age
    {
        int m = LOCATION_PROCESSING_IDX;
        matrix[N2,A] tmp;
        int Nflow = Nflow_country[m];
        int full_flows_on_day[nrow_full_flows_on_day_country];
        int full_flows_on_week[nrow_full_flows_on_week_country];
        int full_flows_on_week_n;
        matrix[N2,A] impact_intv =
            country_impact(beta,
                upswing_rnde[m],
                N2,
                A,
                COVARIATES_N,
                avg_cntct[m],
                covariates[m]
                );
                
        full_flows_on_day = full_flows_on_day_country[,m];
        full_flows_on_week = full_flows_on_week_country[,m];
        full_flows_on_week_n=full_flows_on_week_n_country[m];
        
        E_casesByAge =
            country_EcasesByAge(R0[m],
                e_cases_N0[m],
                beta,
                upswing_rnde[m],
                log_relsusceptibility_age,
                log_reltransmissibility_age,
                N0,
                N2,
                A,
                COVARIATES_N,
                SI_CUT,
                wkend_idx[1:WKEND_IDX_N[m],m],
                avg_cntct[m],
                covariates[m],
                cntct_weekends_mean[m],
                cntct_weekdays_mean[m],
                rev_serial_interval,
                popByAge_abs[m,],
                init_A
                );
                
      flow = 
        country_EflowsByLowDimAge(R0[m],
                                  e_cases_N0[m],
                                  beta,
                                  upswing_rnde[m],
                                  log_relsusceptibility_age,
                                  log_reltransmissibility_age,
                                  N0,
                                  N2,
                                  A,
                                  COVARIATES_N,
                                  SI_CUT,
                                  wkend_idx[1:WKEND_IDX_N[m],m],
                                  avg_cntct[m],
                                  covariates[m],
                                  cntct_weekends_mean[m],
                                  cntct_weekdays_mean[m],
                                  rev_serial_interval,
                                  popByAge_abs[m,],
                                  init_A,
                                  Nflow,
                                  N_low_dim_age_band,
                                  nrow_low_dim_age_band,
                                  low_dim_age_band,
                                  E_casesByAge,
                                  week_start_index[,m]
                                  );
   
       lambdaByAge = country_lambdaByAge(R0[m],
                                  beta,
                                  upswing_rnde[m],
                                  log_relsusceptibility_age,
                                  log_reltransmissibility_age,
                                  N0,
                                  N2,
                                  A,
                                  COVARIATES_N,
                                  SI_CUT,
                                  wkend_idx[1:WKEND_IDX_N[m],m],
                                  avg_cntct[m],
                                  covariates[m],
                                  cntct_weekends_mean[m],
                                  cntct_weekdays_mean[m],
                                  rev_serial_interval,
                                  E_casesByAge
                                  );
                                  
                                  
        if(debug_flow){
            differences = check_country_EflowsByLowDimAge(nrow_low_dim_age_band,
                      Nflow,
                      N0,
                      N2,
                      E_casesByAge,
                      flow,
                      N_low_dim_age_band,
                      low_dim_age_band,
                      week_start_index[,m]
                      );
        }
          
        full_flow = country_EflowsByHighDimAge(R0[m],
                    e_cases_N0[m],
                    beta,
                    upswing_rnde[m],
                    log_relsusceptibility_age,
                    log_reltransmissibility_age,
                    N0,
                    N2,
                    A,
                    COVARIATES_N,
                    SI_CUT,
                    wkend_idx[1:WKEND_IDX_N[m],m],
                    avg_cntct[m],
                    covariates[m],
                    cntct_weekends_mean[m],
                    cntct_weekdays_mean[m],
                    rev_serial_interval,
                    popByAge_abs[m,],
                    init_A,
                    full_flows_on_week_n,
                    full_flows_on_week,
                    full_flows_on_day,
                    E_casesByAge);
                                    
                      
        E_deathsByAge =
            country_EdeathsByAge(
                R0[m],
                e_cases_N0[m],
                beta,
                upswing_rnde[m],
                log_relsusceptibility_age,
                log_reltransmissibility_age,
                phi,
                N0,
                N2,
                A,
                COVARIATES_N,
                SI_CUT,
                wkend_idx[1:WKEND_IDX_N[m],m],
                avg_cntct[m],
                covariates[m],
                cntct_weekends_mean[m],
                cntct_weekdays_mean[m],
                rev_ifr_daysSinceInfection,
                log_ifr_age_base,
                log_ifr_age_rnde_mid1[m],
                log_ifr_age_rnde_mid2[m],
                log_ifr_age_rnde_old[m],
                rev_serial_interval,
                popByAge_abs[m,],
                init_A
                );
    
    
        // generate total expected deaths
        E_deaths = E_deathsByAge * ones_vector_A;
    
        // generate other model parameters
        rho0 = R0[m] / avg_cntct[m];
    
        // generate R_ta
        RtByAge =
            country_Rta(
                R0[m] / avg_cntct[m],
                log_relsusceptibility_age,
                log_reltransmissibility_age,
                impact_intv,
                E_casesByAge,
                N2,
                A,
                wkend_idx[1:WKEND_IDX_N[m],m],
                cntct_weekends_mean[m],
                cntct_weekdays_mean[m],
                popByAge_abs[m,]
                );
        
        // generate Rt as weighted avg of R_ta
        Rt = RtByAge * popByAge[:,m];
    }
}



