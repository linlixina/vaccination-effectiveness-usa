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
        real dip_rdeff_local,
        real[] upswing_timeeff_reduced_local,
        int N2,
        int A,
        int A_CHILD,
        int[] AGE_CHILD,
        int COVARIATES_N,
        matrix[] covariates_local,
        row_vector upswing_age_rdeff_local,
        int[] upswing_timeeff_map_local
        )
    {
        // scaling of contacts after intervention effect on day t in location m
        matrix[N2,A] impact_intv;
    
        // define multipliers for contacts in each location
        impact_intv = (beta[2] + rep_matrix(upswing_age_rdeff_local, N2)) .* covariates_local[3];
        
        // expand upswing time effect
        impact_intv .*= rep_matrix( to_vector( upswing_timeeff_reduced_local[ upswing_timeeff_map_local ]), A );
        
        // add other coeff*predictors
        impact_intv += covariates_local[1];
        impact_intv += ( beta[1] + dip_rdeff_local ) * covariates_local[2];
        
        impact_intv = exp( impact_intv );
        
        // impact_intv set to 1 for children
        impact_intv[:,AGE_CHILD] = rep_matrix(1.0, N2, A_CHILD);
        
        return (impact_intv);
    }
  
    row_vector country_lambdaByAge_t(real R0_local_strain,
                              real avg_cntct_local,
                              row_vector log_relsusceptibility_age,
                              real impact_intv_children_effect,
                              real impact_intv_onlychildren_effect,
                              matrix impact_intv,
                              // data
                              int elementary_school_reopening_idx_local,
                              int elt_reopening_index_local,
                              vector SCHOOL_STATUS_local,
                              int A,
                              int A_CHILD,
                              int SI_CUT,
                              int[] wkend_idx_local,
                              matrix cntct_weekends_mean_local,
                              matrix cntct_weekdays_mean_local,
                              matrix cntct_school_closure_weekends_local,
                              matrix cntct_school_closure_weekdays_local,
                              matrix cntct_elementary_school_reopening_weekends_local,
                              matrix cntct_elementary_school_reopening_weekdays_local,
                              row_vector rev_serial_interval,
                              matrix E_casesByAge,
                              int t
                              )
        {
    
        row_vector[A] tmp_row_vector_A;
        row_vector[A] tmp_row_vector_A_no_impact_intv;
        row_vector[A] tmp_row_vector_A_with_children_impact_intv;
        row_vector[A] tmp_row_vector_A_with_children_and_childrenchildren_impact_intv;
        matrix[A,A] tmp_matrix_A;
        matrix[A,A] cntct_counterfactual_school_closure;
        int start_idx_rev_serial = SI_CUT-t+2;
        int start_idx_E_casesByAge = t-SI_CUT;
        
        // probability of infection given contact in location m
        real rho0 = R0_local_strain / avg_cntct_local;
        
        // expected cases at time t
        row_vector[A] E_lambdaByAge_t;
        
        if(start_idx_rev_serial < 1)
            {
                start_idx_rev_serial = 1;
            }
        if(start_idx_E_casesByAge < 1)
        {
            start_idx_E_casesByAge = 1;
        }
        
        tmp_row_vector_A = rev_serial_interval[start_idx_rev_serial:SI_CUT] * E_casesByAge[start_idx_E_casesByAge:(t-1)];
        tmp_row_vector_A *= rho0;
        tmp_row_vector_A_no_impact_intv = tmp_row_vector_A;
        tmp_row_vector_A .*= impact_intv[t,];
            
        if(r_in(t, wkend_idx_local) == 1)
        {
            if(SCHOOL_STATUS_local[t] == 0 && t<elementary_school_reopening_idx_local)
            {
                E_lambdaByAge_t =
                    append_col(
                        tmp_row_vector_A_no_impact_intv * cntct_weekends_mean_local[,1:A_CHILD],
                        tmp_row_vector_A_no_impact_intv[1:A_CHILD] * cntct_weekends_mean_local[1:A_CHILD,(A_CHILD+1):A] +
                            (tmp_row_vector_A[(A_CHILD+1):A] * cntct_weekends_mean_local[(A_CHILD+1):A,(A_CHILD+1):A]) .* impact_intv[t,(A_CHILD+1):A]
                        );
            }
            else if(SCHOOL_STATUS_local[t] == 0 && t>=elementary_school_reopening_idx_local)
            {
                tmp_row_vector_A_with_children_impact_intv = tmp_row_vector_A;
                tmp_row_vector_A_with_children_impact_intv[1:A_CHILD] *= impact_intv_children_effect;
                tmp_row_vector_A_with_children_and_childrenchildren_impact_intv = tmp_row_vector_A_with_children_impact_intv;
                tmp_row_vector_A_with_children_and_childrenchildren_impact_intv[1:A_CHILD] *= impact_intv_onlychildren_effect;
                
                E_lambdaByAge_t =
                    append_col(
                        tmp_row_vector_A_with_children_and_childrenchildren_impact_intv * cntct_elementary_school_reopening_weekends_local[,1:A_CHILD],
                        tmp_row_vector_A_with_children_impact_intv * cntct_elementary_school_reopening_weekends_local[,(A_CHILD+1):A]
                        );
                E_lambdaByAge_t[1:A_CHILD] *= impact_intv_children_effect;
                E_lambdaByAge_t[(A_CHILD+1):A] .*= impact_intv[t,(A_CHILD+1):A];
            }
            // counterfactual scenario schools remain closed
            else if(SCHOOL_STATUS_local[t] == 1 && t>=elt_reopening_index_local)
            {
                // for the counterfactual scenarios, adjust the contact matrix so that it is at most as large as the inferred contact matrix
                tmp_matrix_A = rep_matrix( append_col(
                    impact_intv[t, 1:A_CHILD] * impact_intv_children_effect,
                    impact_intv[t, (A_CHILD+1):A]
                    ), A );
                tmp_matrix_A = tmp_matrix_A .* cntct_elementary_school_reopening_weekends_local .* (tmp_matrix_A');
                tmp_matrix_A[1:A_CHILD,1:A_CHILD] *= impact_intv_onlychildren_effect;
                cntct_counterfactual_school_closure = cntct_school_closure_weekends_local;
                for( a in 1:A_CHILD)
                {
                    // set children -> any
                    for( aprime in 1:A)
                    {
                        cntct_counterfactual_school_closure[a, aprime] =
                            min([
                                cntct_counterfactual_school_closure[a, aprime],
                                tmp_matrix_A[a, aprime]
                                ]);
                    }
                    // set adult -> children
                    for( aprime in (A_CHILD+1):A)
                    {
                        cntct_counterfactual_school_closure[aprime, a] =
                            min([
                                cntct_counterfactual_school_closure[aprime, a],
                                tmp_matrix_A[aprime, a]
                                ]);
                    }
                }
                // calculate cases as in school closure scenario
                E_lambdaByAge_t =
                    append_col(
                        tmp_row_vector_A_no_impact_intv * cntct_counterfactual_school_closure[,1:A_CHILD],
                        tmp_row_vector_A_no_impact_intv[1:A_CHILD] * cntct_counterfactual_school_closure[1:A_CHILD,(A_CHILD+1):A] +
                            (tmp_row_vector_A[(A_CHILD+1):A] * cntct_counterfactual_school_closure[(A_CHILD+1):A,(A_CHILD+1):A]) .* impact_intv[t,(A_CHILD+1):A]
                        );
            }
            else
            {
                E_lambdaByAge_t =
                    append_col(
                        tmp_row_vector_A_no_impact_intv * cntct_school_closure_weekends_local[,1:A_CHILD],
                        tmp_row_vector_A_no_impact_intv[1:A_CHILD] * cntct_school_closure_weekends_local[1:A_CHILD,(A_CHILD+1):A] +
                            (tmp_row_vector_A[(A_CHILD+1):A] * cntct_school_closure_weekends_local[(A_CHILD+1):A,(A_CHILD+1):A]) .* impact_intv[t,(A_CHILD+1):A]
                        );
            }
        }
        else
        {
            if(SCHOOL_STATUS_local[t] == 0 && t<elementary_school_reopening_idx_local)
            {
                E_lambdaByAge_t =
                    append_col(
                        tmp_row_vector_A_no_impact_intv * cntct_weekdays_mean_local[,1:A_CHILD],
                        tmp_row_vector_A_no_impact_intv[1:A_CHILD] * cntct_weekdays_mean_local[1:A_CHILD,(A_CHILD+1):A] +
                            (tmp_row_vector_A[(A_CHILD+1):A] * cntct_weekdays_mean_local[(A_CHILD+1):A,(A_CHILD+1):A]) .* impact_intv[t,(A_CHILD+1):A]
                        );
            }
            else if(SCHOOL_STATUS_local[t] == 0 && t>=elementary_school_reopening_idx_local)
            {
                tmp_row_vector_A_with_children_impact_intv = tmp_row_vector_A;
                tmp_row_vector_A_with_children_impact_intv[1:A_CHILD] *= impact_intv_children_effect;
                tmp_row_vector_A_with_children_and_childrenchildren_impact_intv = tmp_row_vector_A_with_children_impact_intv;
                tmp_row_vector_A_with_children_and_childrenchildren_impact_intv[1:A_CHILD] *= impact_intv_onlychildren_effect;
                
                E_lambdaByAge_t =
                    append_col(
                        tmp_row_vector_A_with_children_and_childrenchildren_impact_intv * cntct_elementary_school_reopening_weekdays_local[,1:A_CHILD],
                        tmp_row_vector_A_with_children_impact_intv * cntct_elementary_school_reopening_weekdays_local[,(A_CHILD+1):A]
                        );
                E_lambdaByAge_t[1:A_CHILD] *= impact_intv_children_effect;
                E_lambdaByAge_t[(A_CHILD+1):A] .*= impact_intv[t,(A_CHILD+1):A];
            }
            // counterfactual scenario schools remain closed
            else if(SCHOOL_STATUS_local[t] == 1 && t>=elt_reopening_index_local)
            {
                // for the counterfactual scenarios, adjust the contact matrix so that it is at most as large as the inferred contact matrix
                tmp_matrix_A = rep_matrix( append_col(
                    impact_intv[t, 1:A_CHILD] * impact_intv_children_effect,
                    impact_intv[t, (A_CHILD+1):A]
                    ), A );
                tmp_matrix_A = tmp_matrix_A .* cntct_elementary_school_reopening_weekdays_local .* (tmp_matrix_A');
                tmp_matrix_A[1:A_CHILD,1:A_CHILD] *= impact_intv_onlychildren_effect;
                cntct_counterfactual_school_closure = cntct_school_closure_weekdays_local;
                for( a in 1:A_CHILD)
                {
                    // set children -> any
                    for( aprime in 1:A)
                    {
                        cntct_counterfactual_school_closure[a, aprime] =
                            min([
                                cntct_counterfactual_school_closure[a, aprime],
                                tmp_matrix_A[a, aprime]
                                ]);
                    }
                    // set adult -> children
                    for( aprime in (A_CHILD+1):A)
                    {
                        cntct_counterfactual_school_closure[aprime, a] =
                            min([
                                cntct_counterfactual_school_closure[aprime, a],
                                tmp_matrix_A[aprime, a]
                                ]);
                    }
                }
                // calculate cases as in school closure scenario
                E_lambdaByAge_t =
                    append_col(
                        tmp_row_vector_A_no_impact_intv * cntct_counterfactual_school_closure[,1:A_CHILD],
                        tmp_row_vector_A_no_impact_intv[1:A_CHILD] * cntct_counterfactual_school_closure[1:A_CHILD,(A_CHILD+1):A] +
                            (tmp_row_vector_A[(A_CHILD+1):A] * cntct_counterfactual_school_closure[(A_CHILD+1):A,(A_CHILD+1):A]) .* impact_intv[t,(A_CHILD+1):A]
                        );
            }
            else
            {
                E_lambdaByAge_t =
                    append_col(
                        tmp_row_vector_A_no_impact_intv * cntct_school_closure_weekdays_local[,1:A_CHILD],
                        tmp_row_vector_A_no_impact_intv[1:A_CHILD] * cntct_school_closure_weekdays_local[1:A_CHILD,(A_CHILD+1):A] +
                            (tmp_row_vector_A[(A_CHILD+1):A] * cntct_school_closure_weekdays_local[(A_CHILD+1):A,(A_CHILD+1):A]) .* impact_intv[t,(A_CHILD+1):A]
                        );
            }
        }
        
        E_lambdaByAge_t .*= exp(log_relsusceptibility_age);
        
        return(E_lambdaByAge_t);
  }
  
    row_vector country_Rta_t(real rho0_local,
        row_vector prop_susceptibleByAge,
        row_vector log_relsusceptibility_age,
        real impact_intv_children_effect,
        real impact_intv_onlychildren_effect,
        matrix impact_intv,
        // data
        int elementary_school_reopening_idx_local,
        int elt_reopening_index_local,
        int N2,
        vector SCHOOL_STATUS_local,
        int A,
        int A_CHILD,
        int[] wkend_idx_local,
        matrix cntct_weekends_mean_local,
        matrix cntct_weekdays_mean_local,
        matrix cntct_school_closure_weekends_local,
        matrix cntct_school_closure_weekdays_local,
        matrix cntct_elementary_school_reopening_weekends_local,
        matrix cntct_elementary_school_reopening_weekdays_local,
        row_vector popByAge_abs_local,
        int t){
          
        row_vector[A] tmp_row_vector_A;
        row_vector[A] tmp_row_vector_A_no_impact_intv;        
        row_vector[A] tmp_row_vector_A_with_children_impact_intv;
        row_vector[A] tmp_row_vector_A_with_children_and_childrenchildren_impact_intv;
        matrix[A,A] tmp_matrix_A;
        matrix[A,A] cntct_counterfactual_school_closure;
        row_vector[A] RtByAge_t;
        
        RtByAge_t = rep_row_vector(1.0, A) - prop_susceptibleByAge ;
        RtByAge_t .*= exp(log_relsusceptibility_age);
        RtByAge_t *= rho0_local;
        
        tmp_row_vector_A_no_impact_intv = RtByAge_t;
        tmp_row_vector_A = impact_intv[t,:] .* RtByAge_t;
        
        if(r_in(t, wkend_idx_local) == 1)
        {
            if(SCHOOL_STATUS_local[t] == 0 && t<elementary_school_reopening_idx_local)
            {
                RtByAge_t =
                    append_col(
                        tmp_row_vector_A_no_impact_intv * (cntct_weekends_mean_local')[:,1:A_CHILD],
                        tmp_row_vector_A_no_impact_intv[1:A_CHILD] * (cntct_weekends_mean_local')[1:A_CHILD,(A_CHILD+1):A] +
                            (tmp_row_vector_A[(A_CHILD+1):A] * (cntct_weekends_mean_local')[(A_CHILD+1):A,(A_CHILD+1):A]) .* impact_intv[t,(A_CHILD+1):A]
                        );
            }
            else if(SCHOOL_STATUS_local[t] == 0 && t>=elementary_school_reopening_idx_local)
            {
                tmp_row_vector_A_with_children_impact_intv = tmp_row_vector_A;
                tmp_row_vector_A_with_children_impact_intv[1:A_CHILD] *= impact_intv_children_effect;
                tmp_row_vector_A_with_children_and_childrenchildren_impact_intv = tmp_row_vector_A_with_children_impact_intv;
                tmp_row_vector_A_with_children_and_childrenchildren_impact_intv[1:A_CHILD] *= impact_intv_onlychildren_effect;
                                    
                RtByAge_t =
                    append_col(
                        tmp_row_vector_A_with_children_and_childrenchildren_impact_intv * cntct_elementary_school_reopening_weekends_local[,1:A_CHILD],
                        tmp_row_vector_A_with_children_impact_intv * cntct_elementary_school_reopening_weekends_local[,(A_CHILD+1):A]
                        );
                RtByAge_t[1:A_CHILD] *= impact_intv_children_effect;
                RtByAge_t[(A_CHILD+1):A] .*= impact_intv[t,(A_CHILD+1):A];
            }
            // counterfactual scenario schools remain closed
            else if(SCHOOL_STATUS_local[t] == 1 && t>=elt_reopening_index_local)
            {
                // for the counterfactual scenarios, adjust the contact matrix so that it is at most as large as the inferred contact matrix
                tmp_matrix_A = rep_matrix( append_col(
                    impact_intv[t, 1:A_CHILD] * impact_intv_children_effect,
                    impact_intv[t, (A_CHILD+1):A]
                    ), A );
                tmp_matrix_A = tmp_matrix_A .* cntct_elementary_school_reopening_weekends_local .* (tmp_matrix_A');
                tmp_matrix_A[1:A_CHILD, 1:A_CHILD] *= impact_intv_onlychildren_effect;
                cntct_counterfactual_school_closure = cntct_school_closure_weekends_local;
                for( a in 1:A_CHILD)
                {
                    // set children -> any
                    for( aprime in 1:A)
                    {
                        cntct_counterfactual_school_closure[a, aprime] =
                            min([
                                cntct_counterfactual_school_closure[a, aprime],
                                tmp_matrix_A[a, aprime]
                                ]);
                    }
                    // set adult -> children
                    for( aprime in (A_CHILD+1):A)
                    {
                        cntct_counterfactual_school_closure[aprime, a] =
                            min([
                                cntct_counterfactual_school_closure[aprime, a],
                                tmp_matrix_A[aprime, a]
                                ]);
                    }
                }
                RtByAge_t =
                    append_col(
                        tmp_row_vector_A_no_impact_intv * (cntct_counterfactual_school_closure')[:,1:A_CHILD],
                        tmp_row_vector_A_no_impact_intv[1:A_CHILD] * (cntct_counterfactual_school_closure')[1:A_CHILD,(A_CHILD+1):A] +
                            (tmp_row_vector_A[(A_CHILD+1):A] * (cntct_counterfactual_school_closure')[(A_CHILD+1):A,(A_CHILD+1):A]) .* impact_intv[t,(A_CHILD+1):A]
                        );
            }
            else
            {
                RtByAge_t =
                    append_col(
                        tmp_row_vector_A_no_impact_intv * (cntct_school_closure_weekends_local')[:,1:A_CHILD],
                        tmp_row_vector_A_no_impact_intv[1:A_CHILD] * (cntct_school_closure_weekends_local')[1:A_CHILD,(A_CHILD+1):A] +
                            (tmp_row_vector_A[(A_CHILD+1):A] * (cntct_school_closure_weekends_local')[(A_CHILD+1):A,(A_CHILD+1):A]) .* impact_intv[t,(A_CHILD+1):A]
                        );
            }
        }
        else
        {
            if(SCHOOL_STATUS_local[t] == 0 && t<elementary_school_reopening_idx_local)
            {
                RtByAge_t =
                    append_col(
                        tmp_row_vector_A_no_impact_intv * (cntct_weekdays_mean_local')[:,1:A_CHILD],
                        tmp_row_vector_A_no_impact_intv[1:A_CHILD] * (cntct_weekdays_mean_local')[1:A_CHILD,(A_CHILD+1):A] +
                            (tmp_row_vector_A[(A_CHILD+1):A] * (cntct_weekdays_mean_local')[(A_CHILD+1):A,(A_CHILD+1):A]) .* impact_intv[t,(A_CHILD+1):A]
                        );
            }
            else if(SCHOOL_STATUS_local[t] == 0 && t>=elementary_school_reopening_idx_local)
            {
                tmp_row_vector_A_with_children_impact_intv = tmp_row_vector_A;
                tmp_row_vector_A_with_children_impact_intv[1:A_CHILD] *= impact_intv_children_effect;
                tmp_row_vector_A_with_children_and_childrenchildren_impact_intv = tmp_row_vector_A_with_children_impact_intv;
                tmp_row_vector_A_with_children_and_childrenchildren_impact_intv[1:A_CHILD] *= impact_intv_onlychildren_effect;
                                    
                RtByAge_t =
                    append_col(
                        tmp_row_vector_A_with_children_and_childrenchildren_impact_intv * cntct_elementary_school_reopening_weekdays_local[,1:A_CHILD],
                        tmp_row_vector_A_with_children_impact_intv * cntct_elementary_school_reopening_weekdays_local[,(A_CHILD+1):A]
                        );
                RtByAge_t[1:A_CHILD] *= impact_intv_children_effect;
                RtByAge_t[(A_CHILD+1):A] .*= impact_intv[t,(A_CHILD+1):A];
            }
            else if(SCHOOL_STATUS_local[t] == 1 && t>=elt_reopening_index_local)
            {
                // for the counterfactual scenarios, adjust the contact matrix so that it is at most as large as the inferred contact matrix
                tmp_matrix_A = rep_matrix( append_col(
                    impact_intv[t, 1:A_CHILD] * impact_intv_children_effect,
                    impact_intv[t, (A_CHILD+1):A]
                    ), A );
                tmp_matrix_A = tmp_matrix_A .* cntct_elementary_school_reopening_weekdays_local .* (tmp_matrix_A');
                tmp_matrix_A[1:A_CHILD, 1:A_CHILD] *= impact_intv_onlychildren_effect;
                cntct_counterfactual_school_closure = cntct_school_closure_weekdays_local;
                for( a in 1:A_CHILD)
                {
                    // set children -> any
                    for( aprime in 1:A)
                    {
                        cntct_counterfactual_school_closure[a, aprime] =
                            min([
                                cntct_counterfactual_school_closure[a, aprime],
                                tmp_matrix_A[a, aprime]
                                ]);
                    }
                    // set adult -> children
                    for( aprime in (A_CHILD+1):A)
                    {
                        cntct_counterfactual_school_closure[aprime, a] =
                            min([
                                cntct_counterfactual_school_closure[aprime, a],
                                tmp_matrix_A[aprime, a]
                                ]);
                    }
                }
                RtByAge_t =
                    append_col(
                        tmp_row_vector_A_no_impact_intv * (cntct_counterfactual_school_closure')[:,1:A_CHILD],
                        tmp_row_vector_A_no_impact_intv[1:A_CHILD] * (cntct_counterfactual_school_closure')[1:A_CHILD,(A_CHILD+1):A] +
                            (tmp_row_vector_A[(A_CHILD+1):A] * (cntct_counterfactual_school_closure')[(A_CHILD+1):A,(A_CHILD+1):A]) .* impact_intv[t,(A_CHILD+1):A]
                        );
            }
            else
            {
                RtByAge_t =
                   append_col(
                        tmp_row_vector_A_no_impact_intv * (cntct_school_closure_weekdays_local')[:,1:A_CHILD],
                        tmp_row_vector_A_no_impact_intv[1:A_CHILD] * (cntct_school_closure_weekdays_local')[1:A_CHILD,(A_CHILD+1):A] +
                            (tmp_row_vector_A[(A_CHILD+1):A] * (cntct_school_closure_weekdays_local')[(A_CHILD+1):A,(A_CHILD+1):A]) .* impact_intv[t,(A_CHILD+1):A]
                        );
            }
        }
        
        return(RtByAge_t);
    }
  
    matrix[] country_EcasesByAge(// parameters
        vector R0_local,
        real e_cases_N0_local,
        row_vector log_relsusceptibility_age,
        real impact_intv_children_effect,
        real impact_intv_onlychildren_effect,
        matrix impact_intv,
        // data
        int N0,
        int elementary_school_reopening_idx_local,
        int elt_reopening_index_local,
        int N2,
        int N_STRAINS,
        int INTRO_NEW_STRAIN_TIME_IDX_local,
        real prop_cases_new_strain_first_day,
        vector SCHOOL_STATUS_local,
        int A,
        int A_CHILD,
        int SI_CUT,
        int[] wkend_idx_local,
        real avg_cntct_local,
        matrix cntct_weekends_mean_local,
        matrix cntct_weekdays_mean_local,
        matrix cntct_school_closure_weekends_local,
        matrix cntct_school_closure_weekdays_local,
        matrix cntct_elementary_school_reopening_weekends_local,
        matrix cntct_elementary_school_reopening_weekdays_local,
        row_vector rev_serial_interval,
        row_vector popByAge_abs_local,
        int N_init_A,
        int[] init_A
        )
    {
        real zero = 0.0;
        real N_init_A_real= N_init_A*1.;
        row_vector[A] prop_susceptibleByAge;
        // expected new cases by calendar day, age, and location under self-renewal model
        // and a container to store the precomputed cases by age
        matrix[N2,A] E_casesByAge[N_STRAINS] = rep_array(rep_matrix( zero, N2, A ), N_STRAINS);
                  
        // initial strain: init expected cases by age and location in first N0 days
        E_casesByAge[1][1:N0, init_A] = rep_matrix( e_cases_N0_local/N_init_A_real, N0, N_init_A);
        
        // calculate expected cases by age and country under self-renewal model after first N0 days
        // and adjusted for saturation
        for (t in (N0+1):N2)
        {
            prop_susceptibleByAge = rep_row_vector(1.0, A) - 
                                    (rep_row_vector(1.0, t-1) * (E_casesByAge[1][1:(t-1),:] + E_casesByAge[N_STRAINS][1:(t-1),:]) ./ popByAge_abs_local);
            
            for(a in 1:A)
            {
                if(prop_susceptibleByAge[a] < 0)
                { // account for values of Ecases > pop at initalization
                    prop_susceptibleByAge[a] = 0;
                } 
            }
            
            E_casesByAge[1][t,:] = country_lambdaByAge_t(R0_local[1],
                                                          avg_cntct_local,
                                                          log_relsusceptibility_age,
                                                          impact_intv_children_effect,
                                                          impact_intv_onlychildren_effect,
                                                          impact_intv,
                                                          elementary_school_reopening_idx_local,
                                                          elt_reopening_index_local,
                                                          SCHOOL_STATUS_local,
                                                          A,
                                                          A_CHILD,
                                                          SI_CUT,
                                                          wkend_idx_local,
                                                          cntct_weekends_mean_local,
                                                          cntct_weekdays_mean_local,
                                                          cntct_school_closure_weekends_local,
                                                          cntct_school_closure_weekdays_local,
                                                          cntct_elementary_school_reopening_weekends_local,
                                                          cntct_elementary_school_reopening_weekdays_local,
                                                          rev_serial_interval,
                                                          E_casesByAge[1],
                                                          t);
            E_casesByAge[1][t,:] .*= prop_susceptibleByAge;
            
            if(t == INTRO_NEW_STRAIN_TIME_IDX_local){
        
                E_casesByAge[N_STRAINS][(t-SI_CUT+1):t,:] = E_casesByAge[1][(t-SI_CUT+1):t,:] * prop_cases_new_strain_first_day;
                E_casesByAge[1][(t-SI_CUT+1):t,:] = E_casesByAge[1][(t-SI_CUT+1):t,:] - E_casesByAge[N_STRAINS][(t-SI_CUT+1):t,:];
            }
            
            if(t > INTRO_NEW_STRAIN_TIME_IDX_local){
              
              prop_susceptibleByAge = rep_row_vector(1.0, A) - 
                                      (((rep_row_vector(1.0, t) * E_casesByAge[1][1:t,:]) + (rep_row_vector(1.0, t-1) * E_casesByAge[N_STRAINS][1:(t-1),:])) ./ popByAge_abs_local);
              
              for(a in 1:A){
                if(prop_susceptibleByAge[a] < 0)
                { // account for values of Ecases > pop at initalization
                    prop_susceptibleByAge[a] = 0;
                } 
              }
            
              E_casesByAge[N_STRAINS][t,:] = country_lambdaByAge_t(R0_local[N_STRAINS],
                                                          avg_cntct_local,
                                                          log_relsusceptibility_age,
                                                          impact_intv_children_effect,
                                                          impact_intv_onlychildren_effect,
                                                          impact_intv,
                                                          elementary_school_reopening_idx_local,
                                                          elt_reopening_index_local,
                                                          SCHOOL_STATUS_local,
                                                          A,
                                                          A_CHILD,
                                                          SI_CUT,
                                                          wkend_idx_local,
                                                          cntct_weekends_mean_local,
                                                          cntct_weekdays_mean_local,
                                                          cntct_school_closure_weekends_local,
                                                          cntct_school_closure_weekdays_local,
                                                          cntct_elementary_school_reopening_weekends_local,
                                                          cntct_elementary_school_reopening_weekdays_local,
                                                          rev_serial_interval,
                                                          E_casesByAge[N_STRAINS],
                                                          t);
                                                          
              E_casesByAge[N_STRAINS][t,:] .*= prop_susceptibleByAge;
                                                  
            }

        }
        return(E_casesByAge);
    }
    
    matrix[] country_lambdaByAge(// parameters
        vector R0_local,
        row_vector log_relsusceptibility_age,
        real impact_intv_children_effect,
        real impact_intv_onlychildren_effect,
        matrix impact_intv,
        matrix[] E_casesByAge,
        // data
        int N0,
        int elementary_school_reopening_idx_local,
        int elt_reopening_index_local,
        int N2,
        int N_STRAINS,
        int INTRO_NEW_STRAIN_TIME_IDX_local,
        real prop_cases_new_strain_first_day,
        vector SCHOOL_STATUS_local,
        int A,
        int A_CHILD,
        int SI_CUT,
        int[] wkend_idx_local,
        real avg_cntct_local,
        matrix cntct_weekends_mean_local,
        matrix cntct_weekdays_mean_local,
        matrix cntct_school_closure_weekends_local,
        matrix cntct_school_closure_weekdays_local,
        matrix cntct_elementary_school_reopening_weekends_local,
        matrix cntct_elementary_school_reopening_weekdays_local,
        row_vector rev_serial_interval
        )
    {
        real zero = 0.0;
        matrix[N2, A] lambdaByAge[N_STRAINS] = rep_array(rep_matrix( zero, N2, A ), N_STRAINS);
                
        
            
        // calculate expected cases by age and country under self-renewal model after first N0 days
        // and adjusted for saturation
  
        for (t in (N0+1):N2)
        {
            lambdaByAge[1][t,:] = country_lambdaByAge_t(R0_local[1],
                                                avg_cntct_local,
                                                log_relsusceptibility_age,
                                                impact_intv_children_effect,
                                                impact_intv_onlychildren_effect,
                                                impact_intv,
                                                elementary_school_reopening_idx_local,
                                                elt_reopening_index_local,
                                                SCHOOL_STATUS_local,
                                                A,
                                                A_CHILD,
                                                SI_CUT,
                                                wkend_idx_local,
                                                cntct_weekends_mean_local,
                                                cntct_weekdays_mean_local,
                                                cntct_school_closure_weekends_local,
                                                cntct_school_closure_weekdays_local,
                                                cntct_elementary_school_reopening_weekends_local,
                                                cntct_elementary_school_reopening_weekdays_local,
                                                rev_serial_interval,
                                                E_casesByAge[1],
                                                t);
                                                
            if(t > INTRO_NEW_STRAIN_TIME_IDX_local){
              lambdaByAge[N_STRAINS][t,:] = country_lambdaByAge_t(R0_local[N_STRAINS],
                                                          avg_cntct_local,
                                                          log_relsusceptibility_age,
                                                          impact_intv_children_effect,
                                                          impact_intv_onlychildren_effect,
                                                          impact_intv,
                                                          elementary_school_reopening_idx_local,
                                                          elt_reopening_index_local,
                                                          SCHOOL_STATUS_local,
                                                          A,
                                                          A_CHILD,
                                                          SI_CUT,
                                                          wkend_idx_local,
                                                          cntct_weekends_mean_local,
                                                          cntct_weekdays_mean_local,
                                                          cntct_school_closure_weekends_local,
                                                          cntct_school_closure_weekdays_local,
                                                          cntct_elementary_school_reopening_weekends_local,
                                                          cntct_elementary_school_reopening_weekdays_local,
                                                          rev_serial_interval,
                                                          E_casesByAge[N_STRAINS],
                                                          t);
                                                  
            }

        }
        
        return(lambdaByAge);
    }
   
    matrix country_contacts(// parameters
 				real impact_intv_children_effect,
				real impact_intv_onlychildren_effect,
        matrix impact_intv,
       // data
        int N0,
        int elementary_school_reopening_idx_local,
        int elt_reopening_index_local,
        int N2,
        vector SCHOOL_STATUS_local,
        int A,
        int A_CHILD,
        int[] wkend_idx_local,
        matrix cntct_weekends_mean_local,
        matrix cntct_weekdays_mean_local,
        matrix cntct_school_closure_weekends_local,
        matrix cntct_school_closure_weekdays_local,
        matrix cntct_elementary_school_reopening_weekends_local,
        matrix cntct_elementary_school_reopening_weekdays_local
        )
    {
        row_vector[A] tmp_row_vector_A;
        row_vector[A] tmp_row_vector_A_no_impact_intv;
        row_vector[A] tmp_row_vector_A_with_children_impact_intv;
	      row_vector[A] tmp_row_vector_A_with_children_and_childrenchildren_impact_intv;
	      matrix[A,A] tmp_matrix_A;
        matrix[A,A] cntct_counterfactual_school_closure;
        matrix[N2, A] contacts_ByAge;

    	contacts_ByAge[1:N0, 1:A] = rep_matrix(0.,N0,A);
    	
    	for (t in (N0+1):N2)
    	{
        
        tmp_row_vector_A_no_impact_intv = rep_row_vector(1., A);
        tmp_row_vector_A = impact_intv[t,];

    		if(r_in(t, wkend_idx_local) == 1)
    		{
    			if(SCHOOL_STATUS_local[t] == 0 && t<elementary_school_reopening_idx_local)
    			{
    				contacts_ByAge[t] =
    					append_col(
    						tmp_row_vector_A_no_impact_intv * (cntct_weekends_mean_local)'[,1:A_CHILD],
    						tmp_row_vector_A_no_impact_intv[1:A_CHILD] * (cntct_weekends_mean_local)'[1:A_CHILD,(A_CHILD+1):A] +
    							(tmp_row_vector_A[(A_CHILD+1):A] * (cntct_weekends_mean_local)'[(A_CHILD+1):A,(A_CHILD+1):A]) .* impact_intv[t,(A_CHILD+1):A]
    					);
    			}
    			else if(SCHOOL_STATUS_local[t] == 0 && t>=elementary_school_reopening_idx_local)
    			{
    				tmp_row_vector_A_with_children_impact_intv = tmp_row_vector_A;
    				tmp_row_vector_A_with_children_impact_intv[1:A_CHILD] *= impact_intv_children_effect;
    				tmp_row_vector_A_with_children_and_childrenchildren_impact_intv = tmp_row_vector_A_with_children_impact_intv;
    				tmp_row_vector_A_with_children_and_childrenchildren_impact_intv[1:A_CHILD] *= impact_intv_onlychildren_effect;
    				contacts_ByAge[t,:] =
    					append_col(
    						tmp_row_vector_A_with_children_and_childrenchildren_impact_intv * (cntct_elementary_school_reopening_weekends_local)'[,1:A_CHILD],
    						tmp_row_vector_A_with_children_impact_intv * (cntct_elementary_school_reopening_weekends_local)'[,(A_CHILD+1):A]
    					);
    				contacts_ByAge[t, 1:A_CHILD] *= impact_intv_children_effect;
    				contacts_ByAge[t, (A_CHILD+1):A] .*= impact_intv[t,(A_CHILD+1):A];
    			}
    			else if(SCHOOL_STATUS_local[t] == 1 && t>=elt_reopening_index_local)
    			{
    				// for the counterfactual scenarios, adjust the contact matrix so that it is at most as large as the inferred contact matrix
    				tmp_matrix_A = rep_matrix( append_col(
    					impact_intv[t, 1:A_CHILD] * impact_intv_children_effect,
    					impact_intv[t, (A_CHILD+1):A]
    				), A );
    				tmp_matrix_A = tmp_matrix_A .* cntct_elementary_school_reopening_weekends_local .* (tmp_matrix_A');
                        tmp_matrix_A[1:A_CHILD, 1:A_CHILD] *= impact_intv_onlychildren_effect;
                        cntct_counterfactual_school_closure = cntct_school_closure_weekends_local;
                        for( a in 1:A_CHILD)
                        {
                            // set children -> any
                            for( aprime in 1:A)
                            {
                                cntct_counterfactual_school_closure[a, aprime] =
                                    min([
                                        cntct_counterfactual_school_closure[a, aprime],
                                        tmp_matrix_A[a, aprime]
                                        ]);
                            }
                            // set adult -> children
                            for( aprime in (A_CHILD+1):A)
                            {
                                cntct_counterfactual_school_closure[aprime, a] =
                                    min([
                                        cntct_counterfactual_school_closure[aprime, a],
                                        tmp_matrix_A[aprime, a]
                                        ]);
                            }
                        }
                        
                        contacts_ByAge[t] =
                            append_col(
                                tmp_row_vector_A_no_impact_intv * (cntct_counterfactual_school_closure)'[,1:A_CHILD],
                                tmp_row_vector_A_no_impact_intv[1:A_CHILD] * (cntct_counterfactual_school_closure)'[1:A_CHILD,(A_CHILD+1):A] +
                                    (tmp_row_vector_A[(A_CHILD+1):A] * (cntct_counterfactual_school_closure)'[(A_CHILD+1):A,(A_CHILD+1):A]) .* impact_intv[t,(A_CHILD+1):A]
                                );
                    }
                    else
                    {
                        contacts_ByAge[t] =
                            append_col(
                                tmp_row_vector_A_no_impact_intv * (cntct_school_closure_weekends_local)'[,1:A_CHILD],
                                tmp_row_vector_A_no_impact_intv[1:A_CHILD] * (cntct_school_closure_weekends_local)'[1:A_CHILD,(A_CHILD+1):A] +
                                    (tmp_row_vector_A[(A_CHILD+1):A] * (cntct_school_closure_weekends_local)'[(A_CHILD+1):A,(A_CHILD+1):A]) .* impact_intv[t,(A_CHILD+1):A]
                                );
                    }
                }
                else
                {
                    if(SCHOOL_STATUS_local[t] == 0 && t<elementary_school_reopening_idx_local)
                    {
                        contacts_ByAge[t] =
                            append_col(
                                tmp_row_vector_A_no_impact_intv * (cntct_weekdays_mean_local)'[,1:A_CHILD],
                                tmp_row_vector_A_no_impact_intv[1:A_CHILD] * (cntct_weekdays_mean_local)'[1:A_CHILD,(A_CHILD+1):A] +
                                    (tmp_row_vector_A[(A_CHILD+1):A] * (cntct_weekdays_mean_local)'[(A_CHILD+1):A,(A_CHILD+1):A]) .* impact_intv[t,(A_CHILD+1):A]
                                );
                    }
                    else if(SCHOOL_STATUS_local[t] == 0 && t>=elementary_school_reopening_idx_local)
                    {
                        tmp_row_vector_A_with_children_impact_intv = tmp_row_vector_A;
                        tmp_row_vector_A_with_children_impact_intv[1:A_CHILD] *= impact_intv_children_effect;
                        tmp_row_vector_A_with_children_and_childrenchildren_impact_intv = tmp_row_vector_A_with_children_impact_intv;
                        tmp_row_vector_A_with_children_and_childrenchildren_impact_intv[1:A_CHILD] *= impact_intv_onlychildren_effect;
                        contacts_ByAge[t,:] =
                            append_col(
                                tmp_row_vector_A_with_children_and_childrenchildren_impact_intv * (cntct_elementary_school_reopening_weekdays_local)'[,1:A_CHILD],
                                tmp_row_vector_A_with_children_impact_intv * (cntct_elementary_school_reopening_weekdays_local)'[,(A_CHILD+1):A]
                                );
                        contacts_ByAge[t, 1:A_CHILD] *= impact_intv_children_effect;
                        contacts_ByAge[t, (A_CHILD+1):A] .*= impact_intv[t,(A_CHILD+1):A];                                        
                    }
                    else if(SCHOOL_STATUS_local[t] == 1 && t>=elt_reopening_index_local)
                    {
                        // for the counterfactual scenarios, adjust the contact matrix so that it is at most as large as the inferred contact matrix
                        tmp_matrix_A = rep_matrix( append_col(
                            impact_intv[t, 1:A_CHILD] * impact_intv_children_effect,
                            impact_intv[t, (A_CHILD+1):A]
                            ), A );
                        tmp_matrix_A = tmp_matrix_A .* cntct_elementary_school_reopening_weekdays_local .* (tmp_matrix_A');
    				tmp_matrix_A[1:A_CHILD, 1:A_CHILD] *= impact_intv_onlychildren_effect;
    				cntct_counterfactual_school_closure = cntct_school_closure_weekdays_local;
    				for( a in 1:A_CHILD)
    				{
    					// set children -> any
    					for( aprime in 1:A)
    					{
    						cntct_counterfactual_school_closure[a, aprime] =
    							min([
    								cntct_counterfactual_school_closure[a, aprime],
    								tmp_matrix_A[a, aprime]
    								]);
    					}
    					// set adult -> children
    					for( aprime in (A_CHILD+1):A)
    					{
    						cntct_counterfactual_school_closure[aprime, a] =
    							min([
    								cntct_counterfactual_school_closure[aprime, a],
    								tmp_matrix_A[aprime, a]
    								]);
    					}
    				}
    				
    				contacts_ByAge[t] =
    					append_col(
    						tmp_row_vector_A_no_impact_intv * (cntct_counterfactual_school_closure)'[,1:A_CHILD],
    						tmp_row_vector_A_no_impact_intv[1:A_CHILD] * (cntct_counterfactual_school_closure)'[1:A_CHILD,(A_CHILD+1):A] +
    							(tmp_row_vector_A[(A_CHILD+1):A] * (cntct_counterfactual_school_closure)'[(A_CHILD+1):A,(A_CHILD+1):A]) .* impact_intv[t,(A_CHILD+1):A]
    					);
    			}
    			else
    			{
    				contacts_ByAge[t] =
    					append_col(
    						tmp_row_vector_A_no_impact_intv * (cntct_school_closure_weekdays_local)'[,1:A_CHILD],
    						tmp_row_vector_A_no_impact_intv[1:A_CHILD] * (cntct_school_closure_weekdays_local)'[1:A_CHILD,(A_CHILD+1):A] +
    							(tmp_row_vector_A[(A_CHILD+1):A] * (cntct_school_closure_weekdays_local)'[(A_CHILD+1):A,(A_CHILD+1):A]) .* impact_intv[t,(A_CHILD+1):A]
    					);
    			}
    		}
    	}
        return(contacts_ByAge);
    }
    
    matrix[] country_Rta(// parameters
        vector rho0_local,
        row_vector log_relsusceptibility_age,
        real impact_intv_children_effect,
        real impact_intv_onlychildren_effect,
        matrix impact_intv,
        matrix[] E_casesByAge_local,
        // data
        int elementary_school_reopening_idx_local,
        int elt_reopening_index_local,
        int N2,
        int N_STRAINS,
        int INTRO_NEW_STRAIN_TIME_IDX_local,
        vector SCHOOL_STATUS_local,
        int A,
        int A_CHILD,
        int[] wkend_idx_local,
        matrix cntct_weekends_mean_local,
        matrix cntct_weekdays_mean_local,
        matrix cntct_school_closure_weekends_local,
        matrix cntct_school_closure_weekdays_local,
        matrix cntct_elementary_school_reopening_weekends_local,
        matrix cntct_elementary_school_reopening_weekdays_local,
        row_vector popByAge_abs_local
        )
    {
      real zero = 0.0;
      matrix[N2,A] RtByAge[N_STRAINS] = rep_array(rep_matrix( zero, N2, A ), N_STRAINS);
      matrix[N2,A] prop_susceptibleByAge[N_STRAINS];
              
      for(a in 1:A)
      {
          prop_susceptibleByAge[1][1,a] = 0;
          prop_susceptibleByAge[1][2:N2,a] = cumulative_sum( E_casesByAge_local[1][1:(N2-1),a] + E_casesByAge_local[N_STRAINS][1:(N2-1),a] ) / popByAge_abs_local[a];
          
          prop_susceptibleByAge[N_STRAINS][1:INTRO_NEW_STRAIN_TIME_IDX_local,a] = prop_susceptibleByAge[1][1:INTRO_NEW_STRAIN_TIME_IDX_local, a];
          prop_susceptibleByAge[N_STRAINS][(INTRO_NEW_STRAIN_TIME_IDX_local+1):N2,a] = (cumulative_sum( E_casesByAge_local[1][1:N2,a])[(INTRO_NEW_STRAIN_TIME_IDX_local+1):N2] + 
                                                                                        cumulative_sum( E_casesByAge_local[N_STRAINS][1:(N2-1),a] )[INTRO_NEW_STRAIN_TIME_IDX_local:(N2-1)]) / popByAge_abs_local[a];  
          for( t in 1:N2)
          {
              if( prop_susceptibleByAge[1][t,a]>1 )
              {
                  prop_susceptibleByAge[1][t,a]= 1.;
              }
              if( prop_susceptibleByAge[N_STRAINS][t,a]>1 )
              {
                  prop_susceptibleByAge[N_STRAINS][t,a]= 1.;
              }
          }
      }
      
      
      for(t in 1:N2)
      {
       RtByAge[1][t,:] = country_Rta_t(rho0_local[1],
                                    prop_susceptibleByAge[1][t,:],
                                    log_relsusceptibility_age,
                                    impact_intv_children_effect,
                                    impact_intv_onlychildren_effect,
                                    impact_intv,
                                    elementary_school_reopening_idx_local,
                                    elt_reopening_index_local,
                                    N2,
                                    SCHOOL_STATUS_local,
                                    A,
                                    A_CHILD,
                                    wkend_idx_local,
                                    cntct_weekends_mean_local,
                                    cntct_weekdays_mean_local,
                                    cntct_school_closure_weekends_local,
                                    cntct_school_closure_weekdays_local,
                                    cntct_elementary_school_reopening_weekends_local,
                                    cntct_elementary_school_reopening_weekdays_local,
                                    popByAge_abs_local,
                                    t);
                       
          if(t > INTRO_NEW_STRAIN_TIME_IDX_local){
             RtByAge[N_STRAINS][t,:] = country_Rta_t(rho0_local[N_STRAINS],
                                    prop_susceptibleByAge[N_STRAINS][t,:],
                                    log_relsusceptibility_age,
                                    impact_intv_children_effect,
                                    impact_intv_onlychildren_effect,
                                    impact_intv,
                                    elementary_school_reopening_idx_local,
                                    elt_reopening_index_local,
                                    N2,
                                    SCHOOL_STATUS_local,
                                    A,
                                    A_CHILD,
                                    wkend_idx_local,
                                    cntct_weekends_mean_local,
                                    cntct_weekdays_mean_local,
                                    cntct_school_closure_weekends_local,
                                    cntct_school_closure_weekdays_local,
                                    cntct_elementary_school_reopening_weekends_local,
                                    cntct_elementary_school_reopening_weekdays_local,
                                    popByAge_abs_local,
                                    t);
          }             
         
      }
      return(RtByAge);
    }
    
    matrix country_EeffcasesByAge(// parameters
        matrix E_casesByAge,
        // data
        int N0,
        int N2,
        int A,
        int SI_CUT,
        row_vector rev_serial_interval
        )
        {
          row_vector[A] tmp_row_vector_A;
          matrix[N2, A] EeffcasesByAge;
              
           // calculate expected cases by age and country under self-renewal model after first N0 days
          EeffcasesByAge[1:N0, 1:A] = rep_matrix(0.,N0,A);
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
            EeffcasesByAge[t] = tmp_row_vector_A;
          }
          return(EeffcasesByAge);
        }
  
    matrix country_EflowsByHighDimAge_AggregateDays(// parameters
        real rho0_local,
        real impact_intv_children_effect,
        real impact_intv_onlychildren_effect,
        matrix impact_intv,
        matrix E_casesByAge,
        row_vector log_relsusceptibility_age,
        // data
        int elementary_school_reopening_idx_local,
        int elt_reopening_index_local,
        vector SCHOOL_STATUS_local,
        int A,
        int A_CHILD,
        int SI_CUT,
        int[] wkend_idx_local,
        matrix cntct_weekends_mean_local,
        matrix cntct_weekdays_mean_local,
        matrix cntct_school_closure_weekends_local,
        matrix cntct_school_closure_weekdays_local,
        matrix cntct_elementary_school_reopening_weekends_local,
        matrix cntct_elementary_school_reopening_weekdays_local,
        row_vector rev_serial_interval,
        row_vector popByAge_abs_local,
        int full_flows_Monday_idx_local,
        int n_days
        )
    {
        real zero;
        int start_idx_rev_serial;
        int start_idx_E_casesByAge;
        row_vector[A] tmp_row_vector_A;
        row_vector[A] tmp_row_vector_A_no_impact_intv;
        row_vector[A] prop_susceptibleByAge;
        row_vector[A] tmp_row_vector_A_with_children_impact_intv;
        matrix[A,A] tmp_matrix_A;
        matrix[A,A] cntct_counterfactual_school_closure;
        matrix[A,A] tmp_flow;
        matrix[A,A] full_flow;
        
        zero = 0.0;
        full_flow = rep_matrix( zero, A, A);
        
        for (t in full_flows_Monday_idx_local:(full_flows_Monday_idx_local+n_days-1))
        {
            start_idx_rev_serial = SI_CUT-t+2;
            start_idx_E_casesByAge = t-SI_CUT;
            prop_susceptibleByAge = rep_row_vector(1.0, A) - (rep_row_vector(1.0, t-1) * E_casesByAge[1:(t-1),:] ./ popByAge_abs_local);
            
            if(start_idx_rev_serial < 1)
            {
                start_idx_rev_serial = 1;
            }
            if(start_idx_E_casesByAge < 1)
            {
                start_idx_E_casesByAge = 1;
            }
            // account for values of Ecases > pop at initalization
            for(a in 1:A)
            {
                if(prop_susceptibleByAge[a] < 0)
                {
                    prop_susceptibleByAge[a] = 0;
                }
            }
          
            tmp_row_vector_A = rev_serial_interval[start_idx_rev_serial:SI_CUT] * E_casesByAge[start_idx_E_casesByAge:(t-1)];
            tmp_row_vector_A *= rho0_local;
            tmp_row_vector_A_no_impact_intv = tmp_row_vector_A;
            tmp_row_vector_A .*= impact_intv[t,];
            tmp_row_vector_A_with_children_impact_intv = tmp_row_vector_A;
            tmp_row_vector_A_with_children_impact_intv[1:A_CHILD] *= impact_intv_children_effect;
            
            if(r_in(t, wkend_idx_local) == 1)
            {
                if(SCHOOL_STATUS_local[t] == 0 && t<elementary_school_reopening_idx_local)
                {
                    tmp_flow = append_col( rep_matrix(tmp_row_vector_A_no_impact_intv', A_CHILD), rep_matrix(tmp_row_vector_A', A-A_CHILD) ) .* cntct_weekends_mean_local;
                    tmp_flow .*= append_row(rep_matrix(rep_row_vector(1.0, A), A_CHILD),rep_matrix(impact_intv[t,],A-A_CHILD));
                }
                else if(SCHOOL_STATUS_local[t] == 0 && t>=elementary_school_reopening_idx_local)
                {
                    tmp_flow = rep_matrix( tmp_row_vector_A_with_children_impact_intv', A );
                    tmp_flow .*= cntct_elementary_school_reopening_weekends_local;
                    tmp_flow[:,1:A_CHILD] *= impact_intv_children_effect;
                    tmp_flow[:,(A_CHILD+1):A] .*= rep_matrix( impact_intv[t, (A_CHILD+1):A], A);
                    tmp_flow[1:A_CHILD,1:A_CHILD] *= impact_intv_onlychildren_effect;
                }
                else if(SCHOOL_STATUS_local[t] == 1 && t>=elt_reopening_index_local)
                {
                    // for the counterfactual scenarios, adjust the contact matrix so that it is at most as large as the inferred contact matrix
                    tmp_matrix_A = rep_matrix( append_col(
                        impact_intv[t, 1:A_CHILD] * impact_intv_children_effect,
                        impact_intv[t, (A_CHILD+1):A]
                        ), A );
                    tmp_matrix_A = tmp_matrix_A .* cntct_elementary_school_reopening_weekends_local .* (tmp_matrix_A');
                    tmp_matrix_A[1:A_CHILD,1:A_CHILD] *= impact_intv_onlychildren_effect;
                    cntct_counterfactual_school_closure = cntct_school_closure_weekends_local;
                    for( a in 1:A_CHILD)
                    {
                        // set children -> any
                        for( aprime in 1:A)
                        {
                            cntct_counterfactual_school_closure[a, aprime] =
                                min([
                                    cntct_counterfactual_school_closure[a, aprime],
                                    tmp_matrix_A[a, aprime]
                                    ]);
                        }
                        // set adult -> children
                        for( aprime in (A_CHILD+1):A)
                        {
                            cntct_counterfactual_school_closure[aprime, a] =
                                min([
                                    cntct_counterfactual_school_closure[aprime, a],
                                    tmp_matrix_A[aprime, a]
                                    ]);
                        }
                    }
                    tmp_flow = append_col( rep_matrix(tmp_row_vector_A_no_impact_intv', A_CHILD), rep_matrix(tmp_row_vector_A', A-A_CHILD) ) .* cntct_counterfactual_school_closure;
                    tmp_flow .*= append_row(rep_matrix(rep_row_vector(1.0, A), A_CHILD), rep_matrix(impact_intv[t,],A-A_CHILD));
                }
                else
                {
                    tmp_flow = append_col( rep_matrix(tmp_row_vector_A_no_impact_intv', A_CHILD), rep_matrix(tmp_row_vector_A', A-A_CHILD) ) .* cntct_school_closure_weekends_local;
                    tmp_flow .*= append_row(rep_matrix(rep_row_vector(1.0, A), A_CHILD),rep_matrix(impact_intv[t,],A-A_CHILD));
                }
            }
            else
            {
                if(SCHOOL_STATUS_local[t] == 0 && t<elementary_school_reopening_idx_local)
                {
                    tmp_flow = append_col( rep_matrix(tmp_row_vector_A_no_impact_intv', A_CHILD), rep_matrix(tmp_row_vector_A', A-A_CHILD) ) .* cntct_weekdays_mean_local;
                    tmp_flow .*= append_row(rep_matrix(rep_row_vector(1.0, A), A_CHILD),rep_matrix(impact_intv[t,],A-A_CHILD));
                }
                else if(SCHOOL_STATUS_local[t] == 0 && t>=elementary_school_reopening_idx_local)
                {
                    tmp_flow = rep_matrix( tmp_row_vector_A_with_children_impact_intv', A );
                    tmp_flow .*= cntct_elementary_school_reopening_weekdays_local;
                    tmp_flow[:,1:A_CHILD] *= impact_intv_children_effect;
                    tmp_flow[:,(A_CHILD+1):A] .*= rep_matrix( impact_intv[t, (A_CHILD+1):A], A);
                    tmp_flow[1:A_CHILD,1:A_CHILD] *= impact_intv_onlychildren_effect;
                }
                else if(SCHOOL_STATUS_local[t] == 1 && t>=elt_reopening_index_local)
                {
                    // for the counterfactual scenarios, adjust the contact matrix so that it is at most as large as the inferred contact matrix
                    tmp_matrix_A = rep_matrix( append_col(
                        impact_intv[t, 1:A_CHILD] * impact_intv_children_effect,
                        impact_intv[t, (A_CHILD+1):A]
                        ), A );
                    tmp_matrix_A = tmp_matrix_A .* cntct_elementary_school_reopening_weekdays_local .* (tmp_matrix_A');
                    tmp_matrix_A[1:A_CHILD,1:A_CHILD] *= impact_intv_onlychildren_effect;
                    cntct_counterfactual_school_closure = cntct_school_closure_weekdays_local;
                    for( a in 1:A_CHILD)
                    {
                        // set children -> any
                        for( aprime in 1:A)
                        {
                            cntct_counterfactual_school_closure[a, aprime] =
                                min([
                                    cntct_counterfactual_school_closure[a, aprime],
                                    tmp_matrix_A[a, aprime]
                                    ]);
                        }
                        // set adult -> children
                        for( aprime in (A_CHILD+1):A)
                        {
                            cntct_counterfactual_school_closure[aprime, a] =
                                min([
                                    cntct_counterfactual_school_closure[aprime, a],
                                    tmp_matrix_A[aprime, a]
                                    ]);
                        }
                    }
                    tmp_flow = append_col( rep_matrix(tmp_row_vector_A_no_impact_intv', A_CHILD), rep_matrix(tmp_row_vector_A', A-A_CHILD) ) .* cntct_counterfactual_school_closure;
                    tmp_flow .*= append_row(rep_matrix(rep_row_vector(1.0, A), A_CHILD),rep_matrix(impact_intv[t,],A-A_CHILD));
                }
                else
                {
                    tmp_flow = append_col( rep_matrix(tmp_row_vector_A_no_impact_intv', A_CHILD), rep_matrix(tmp_row_vector_A', A-A_CHILD) ) .* cntct_school_closure_weekdays_local;
                    tmp_flow .*= append_row(rep_matrix(rep_row_vector(1.0, A), A_CHILD),rep_matrix(impact_intv[t,],A-A_CHILD));
                }
            }
            tmp_flow .*= rep_matrix(exp(log_relsusceptibility_age),A);
            tmp_flow .*= rep_matrix(prop_susceptibleByAge,A);
            full_flow += tmp_flow;
        }
        
        return(full_flow);
    }
    
    
    matrix[] country_EflowsByHighDimAge(// parameters
        real rho0_local,
        real impact_intv_children_effect,
        real impact_intv_onlychildren_effect,
        matrix impact_intv,
        matrix E_casesByAge,
        row_vector log_relsusceptibility_age,
        // data
        int elementary_school_reopening_idx_local,
        int elt_reopening_index_local,
        vector SCHOOL_STATUS_local,
        int A,
        int A_CHILD,
        int SI_CUT,
        int N_WEEKS_FULL_FLOWS_local,
        int[] wkend_idx_local,
        matrix cntct_weekends_mean_local,
        matrix cntct_weekdays_mean_local,
        matrix cntct_school_closure_weekends_local,
        matrix cntct_school_closure_weekdays_local,
        matrix cntct_elementary_school_reopening_weekends_local,
        matrix cntct_elementary_school_reopening_weekdays_local,
        row_vector rev_serial_interval,
        row_vector popByAge_abs_local,
        int[] full_flows_Monday_idx_local,
        int n_days
        )
    {
        matrix[A,A] full_flows[N_WEEKS_FULL_FLOWS_local];
        for (week in 1:N_WEEKS_FULL_FLOWS_local)
        {
            full_flows[week] =
                country_EflowsByHighDimAge_AggregateDays(
                    rho0_local,
                    impact_intv_children_effect,
                    impact_intv_onlychildren_effect,
                    impact_intv,
                    E_casesByAge,
                    log_relsusceptibility_age,
                    elementary_school_reopening_idx_local,
                    elt_reopening_index_local,
                    SCHOOL_STATUS_local,
                    A,
                    A_CHILD,
                    SI_CUT,
                    wkend_idx_local,
                    cntct_weekends_mean_local,
                    cntct_weekdays_mean_local,
                    cntct_school_closure_weekends_local,
                    cntct_school_closure_weekdays_local,
                    cntct_elementary_school_reopening_weekends_local,
                    cntct_elementary_school_reopening_weekdays_local,
                    rev_serial_interval,
                    popByAge_abs_local,
                    full_flows_Monday_idx_local[week],
                    n_days
                    );
        }
        return(full_flows);
    }
  
  
    matrix[] country_EflowsByLowDimAge(// parameters
          real rho0_local,
          real impact_intv_children_effect,
          real impact_intv_onlychildren_effect,
          matrix impact_intv,
          matrix E_casesByAge,
          row_vector log_relsusceptibility_age,
          // data
          int elementary_school_reopening_idx_local,
          int elt_reopening_index_local,
          vector SCHOOL_STATUS_local,
          int A,
          int A_REDUCED,
          int A_CHILD,
          int SI_CUT,
          int N_WEEKS_REDUCED_FLOWS_local,
          int[] wkend_idx_local,
          matrix cntct_weekends_mean_local,
          matrix cntct_weekdays_mean_local,
          matrix cntct_school_closure_weekends_local,
          matrix cntct_school_closure_weekdays_local,
          matrix cntct_elementary_school_reopening_weekends_local,
          matrix cntct_elementary_school_reopening_weekdays_local,
          row_vector rev_serial_interval,
          row_vector popByAge_abs_local,
          int[] reduced_age_bands_map,
          int[] reduced_flows_Monday_idx_local,
          int n_days
          )
    {
        matrix[A,A] full_flow;
        matrix[A_REDUCED,A_REDUCED] tmp_flow;
        matrix[A_REDUCED,A_REDUCED] reduced_flows[N_WEEKS_REDUCED_FLOWS_local];
                
        for (week in 1:N_WEEKS_REDUCED_FLOWS_local)
        {
            full_flow =
                country_EflowsByHighDimAge_AggregateDays(
                    rho0_local,
                    impact_intv_children_effect,
                    impact_intv_onlychildren_effect,
                    impact_intv,
                    E_casesByAge,
                    log_relsusceptibility_age,
                    elementary_school_reopening_idx_local,
                    elt_reopening_index_local,
                    SCHOOL_STATUS_local,
                    A,
                    A_CHILD,
                    SI_CUT,
                    wkend_idx_local,
                    cntct_weekends_mean_local,
                    cntct_weekdays_mean_local,
                    cntct_school_closure_weekends_local,
                    cntct_school_closure_weekdays_local,
                    cntct_elementary_school_reopening_weekends_local,
                    cntct_elementary_school_reopening_weekdays_local,
                    rev_serial_interval,
                    popByAge_abs_local,
                    reduced_flows_Monday_idx_local[week],
                    n_days
                    );
                    
            tmp_flow = rep_matrix(0.0, A_REDUCED, A_REDUCED);
            for (a in 1:A)
            {
                for (b in 1:A)
                {
                    tmp_flow[ reduced_age_bands_map[a], reduced_age_bands_map[b] ] += full_flow[a,b];
                }
            }
            
            reduced_flows[week] = tmp_flow;
        }
        return(reduced_flows);
    }
              
    matrix country_EdeathsByAge(// parameters
        matrix E_casesByAge_local,
        // data
        int N2,
        int A,
        row_vector rev_ifr_daysSinceInfection,
        row_vector log_ifr_age_base,
        real log_ifr_age_rnde_mid1_local,
        real log_ifr_age_rnde_mid2_local,
        real log_ifr_age_rnde_old_local
        )
    {
        real zero = 0.0;
        
        matrix[N2,A] E_deathsByAge = rep_matrix( zero, N2, A );
    
        // calculate expected deaths by age and country
        E_deathsByAge[1] = 1e-15 * E_casesByAge_local[1];
        for (t in 2:N2)
        {
            E_deathsByAge[t] = rev_ifr_daysSinceInfection[(N2-(t-1)+1):N2 ] * E_casesByAge_local[1:(t-1)];
        }
        E_deathsByAge .*= rep_matrix(exp(   log_ifr_age_base +
            append_col(append_col(append_col(
                rep_row_vector(0., 4),
                rep_row_vector(log_ifr_age_rnde_mid1_local, 6)),
                rep_row_vector(log_ifr_age_rnde_mid2_local, 4)),
                rep_row_vector(log_ifr_age_rnde_old_local, 4))
            ), N2);
        E_deathsByAge += 1e-15;
        return(E_deathsByAge);
    }
    
    matrix country_EantibodyByAge(// parameters
        matrix E_casesByAge_local,
        // data
        int N2,
        int A,
        row_vector rev_iar_daysSinceInfection
        )
    {
        real zero = 0.0;
        
        matrix[N2,A] e_antibodyByAge = rep_matrix( zero, N2, A );
    
        // calculate expected deaths by age and country
        e_antibodyByAge[1] = 1e-15 * E_casesByAge_local[1];
        for (t in 2:N2)
        {
            e_antibodyByAge[t] = rev_iar_daysSinceInfection[(N2-(t-1)+1):N2 ] * E_casesByAge_local[1:(t-1)];
        }
        e_antibodyByAge += 1e-15;
        return(e_antibodyByAge);
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
  int<lower=1>  N_init_A; // number of age bands with initial cases
  int WKEND_IDX_N[M]; // number of weekend indices in each location
  int<lower=1, upper=N2> N_IMP; // number of impact invervention time effects
  //	data
  real pop[M];
  matrix<lower=0, upper=1>[A,M] popByAge; // proportion of age bracket in population in location
  int epidemicStart[M];
  int deaths[N2, M]; // reported deaths -- the rows with i > N contain -1 and should be ignored
  int<lower=1> elementary_school_reopening_idx[M]; // time index after which schools reopen for country m IF school status is set to 0. each entry must be <= N2 in forecast
  int<lower=1> elt_reopening_index[M]; // time index after which schools reopen for country m 
  int<lower=0> wkend_idx[N2,M]; //indices of 1:N2 that correspond to weekends in location m
  int<lower=1,upper=N_IMP> upswing_timeeff_map[N2,M]; // map of impact intv time effects to time units in model for each state
  //
  // mobility trends
  //  
  matrix[N2,A] covariates[M, COVARIATES_N]; // predictors for fsq contacts by age
  //
  // death data by age
  //
  int<lower=0> M_AD; // number of countries with deaths by age data
  int<lower=1> dataByAgestart[M_AD]; // start of death by age data
  int deathsByAge[N2, A, M_AD]; // reported deaths by age -- the rows with i < dataByAgestart[M_AD] contain -1 and should be ignored + the column with j > A2[M_AD] contain -1 and should be ignored 
  int<lower=2> A_AD[M_AD]; // number of age groups reported 
  matrix[A, A] map_age[M_AD]; // map the age groups reported with 5 y age group -- the column with j > A2[M_AD] contain -1 and should be ignored
  int map_country[M,2]; // first column indicates if country has death by age date (1 if yes), 2 column map the country to M_AD
  //
  // case data by age
  //
  int smoothed_logcases_weeks_n_max;
  int<lower=1, upper=smoothed_logcases_weeks_n_max> smoothed_logcases_weeks_n[M]; // number of week indices per location
  int smoothed_logcases_week_map[M, smoothed_logcases_weeks_n_max, 7]; // map of week indices to time indices
  real smoothed_logcases_week_pars[M, smoothed_logcases_weeks_n_max, 3]; // likelihood parameters for observed cases
  //
  // school case data
  //
  int school_case_time_idx[M,2];
  real school_case_data[M,4];  
  //
  // school closure status
  //
  matrix[N2,M] SCHOOL_STATUS; // school status, 1 if close, 0 if open
  //
  // school effect in counterfactual scenarions 
  //
  real counterfactual_school_effect;
  //
  // contact matrices
  //
  int A_CHILD; // number of age band for child
  int A_ADULT; // number of age band for adult
  int<lower=1,upper=A> AGE_CHILD[A_CHILD]; // age bands with child
  matrix[A,A] cntct_school_closure_weekdays[M]; // min cntct_weekdays_mean and contact intensities during outbreak estimated in Zhang et al 
  matrix[A,A] cntct_school_closure_weekends[M]; // min cntct_weekends_mean and contact intensities during outbreak estimated in Zhang et al
  matrix[A,A] cntct_elementary_school_reopening_weekdays[M]; // contact matrix for school reopening
  matrix[A,A] cntct_elementary_school_reopening_weekends[M]; // contact matrix for school reopening
  // SARS-CoV-2 strain 
  int N_STRAINS; // total number of SARS-CoV-2 strains
  int INTRO_NEW_STRAIN_TIME_IDX[M]; // time index of the introduction of the new strain
  real rel_transmissibility_new_strain; // transmissibility multiplier new SARS-CoV-2 strains
  real prop_cases_new_strain_first_day; // proportion of new strain among cases on introduction day
  //	priors
  matrix[A,A] cntct_weekdays_mean[M]; // mean of prior contact rates between age groups on weekdays
  matrix[A,A] cntct_weekends_mean[M]; // mean of prior contact rates between age groups on weekends
  real<upper=0> hyperpara_ifr_age_lnmu[A];  // hyper-parameters for probability of death in age band a log normal mean
  real<lower=0> hyperpara_ifr_age_lnsd[A];  // hyper-parameters for probability of death in age band a log normal sd
  row_vector[N2] rev_ifr_daysSinceInfection; // probability of death s days after infection in reverse order
  row_vector[N2] rev_iar_daysSinceInfection; // probability of antibody after infection in reverse order
  row_vector[SI_CUT] rev_serial_interval; // fixed pre-calculated serial interval using empirical data from Neil in reverse order
  int<lower=1, upper=A> init_A[N_init_A]; // age band in which initial cases occur in the first N0 days
  //    for generating quantities
  int<lower=1, upper=M> LOCATION_PROCESSING_IDX;
  int<lower=1> NMAX_WEEKS_FULL_FLOWS;
  int<lower=1, upper=NMAX_WEEKS_FULL_FLOWS> N_WEEKS_FULL_FLOWS[M];
  int<lower=-1, upper=N2-6> full_flows_Monday_idx[NMAX_WEEKS_FULL_FLOWS, M];
  int<lower=1, upper=A> A_REDUCED;
  int<lower=1> NMAX_WEEKS_REDUCED_FLOWS;
  int<lower=1, upper=NMAX_WEEKS_REDUCED_FLOWS> N_WEEKS_REDUCED_FLOWS[M];
  int<lower=1> n_days;
  int<lower=1, upper=A_REDUCED> reduced_age_bands_map[A];       // this is just {1,1,2,2,3,3,3,4,4,4,5,5,5,6,6,6,7,7}
  int<lower=-1, upper=N2-6> reduced_flows_Monday_idx[NMAX_WEEKS_REDUCED_FLOWS, M];
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
  real<lower=0> e_cases_N0[M]; // expected number of cases per day in the first N0 days, for each country
  vector[COVARIATES_N-1] beta; // regression coefficients for time varying multipliers on contacts
  real dip_rnde[M];
  real<lower=0> sd_dip_rnde;  
  real<lower=0> phi; // overdispersion parameter for likelihood model
  row_vector<upper=0>[A] log_ifr_age_base; // probability of death for age band a
  real<lower=0> hyper_log_ifr_age_rnde_mid1;
  real<lower=0> hyper_log_ifr_age_rnde_mid2;
  real<lower=0> hyper_log_ifr_age_rnde_old;
  row_vector<lower=0>[M] log_ifr_age_rnde_mid1;
  row_vector<lower=0>[M] log_ifr_age_rnde_mid2;
  row_vector<lower=0>[M] log_ifr_age_rnde_old;
  row_vector[2] log_relsusceptibility_age_reduced;
  real<lower=0> upswing_timeeff_reduced[N_IMP,M];
  real<lower=0> sd_upswing_timeeff_reduced;
  real<lower=0> hyper_timeeff_shift_mid1;
  vector<lower=0>[M] timeeff_shift_mid1;
  real<lower=0.1,upper=1.0> impact_intv_children_effect;
  real<lower=0> impact_intv_onlychildren_effect;
}

generated quantities
{
    vector[N_STRAINS] rho0;
    real impact_intv_children_effect_gqs;
    real impact_intv_onlychildren_effect_gqs;
    row_vector[A] log_relsusceptibility_age;
    matrix[M,A] timeeff_shift_age;
    matrix[N2,A] impact_intv;
    vector<lower=0>[N2] E_casesByStrain[N_STRAINS];
    vector<lower=0>[N2] E_deathsByStrain[N_STRAINS];
    matrix<lower=0>[N2,A] E_effcasesByAgeByStrain[N_STRAINS];
    matrix<lower=0>[N2,A] RtByAgeByStrain[N_STRAINS];
    matrix<lower=0>[N2,A] E_casesByAge;
    matrix<lower=0>[N2,A] E_deathsByAge;
    vector<lower=0>[N2] E_deaths;
    vector<lower=0>[N2] E_cases;
    matrix<lower=0>[N2,A] E_effcasesByAge;
    matrix<lower=0>[N2,A] E_antibodyByAge;
    vector<lower=0>[N2] E_antibody;
    matrix<lower=0>[N2,A] contacts_ByAge;
    real<lower=0> R0_new_strain_m; // R0 new strain
    vector<lower=0>[N_STRAINS] R0_full_m; // R0 all strains
    // matrix<lower=0>[N2,A] E_casesByAgeByStrain[N_STRAINS];
    // matrix<lower=0>[N2,A] E_deathsByAgeByStrain[N_STRAINS];
    // matrix<lower=0>[N2,A] E_antibodyByAgeByStrain[N_STRAINS];
    // vector<lower=0>[N2] E_antibodyByStrain[N_STRAINS];
    // matrix<lower=0>[N2,A] lambdaByAgeByStrain[N_STRAINS];
    // vector<lower=0>[N2] RtByStrain[N_STRAINS];
    // matrix[A_REDUCED,A_REDUCED] reduced_flowsByStrain[ N_STRAINS, N_WEEKS_REDUCED_FLOWS[LOCATION_PROCESSING_IDX] ];
    // matrix[A,A] full_flowsByStrain[ N_STRAINS, N_WEEKS_FULL_FLOWS[LOCATION_PROCESSING_IDX] ];
    
    // reset "impact_intv_children_effect" and "impact_intv_onlychildren_effect" when we are in the counterfactual re-opening scenario
    impact_intv_children_effect_gqs = impact_intv_children_effect;
    impact_intv_onlychildren_effect_gqs = impact_intv_onlychildren_effect;
    if(counterfactual_school_effect > 0)
    {
        impact_intv_children_effect_gqs = counterfactual_school_effect;
        impact_intv_onlychildren_effect_gqs = 1.;
    }
    
    // susceptibility by age
    log_relsusceptibility_age = append_col( append_col( log_relsusceptibility_age_reduced[ { 1, 1, 1 } ],
        rep_row_vector(0., 10) ),
        log_relsusceptibility_age_reduced[ { 2,2,2,2,2 } ]
        );
        
    timeeff_shift_age =
        append_col(
            append_col(
                rep_matrix(0., M, 4),
                rep_matrix(timeeff_shift_mid1, 6)
                ),
            rep_matrix(0., M, 8)
        );
        
    
    // generate expected cases by age + expected deaths by age
    {
        int m = LOCATION_PROCESSING_IDX;
        matrix[N2,A] E_casesByAgeByStrain[N_STRAINS];
        matrix[N2,A] E_deathsByAgeByStrain[N_STRAINS];
        matrix[N2,A] E_antibodyByAgeByStrain[N_STRAINS];
        vector[N2] E_antibodyByStrain[N_STRAINS];
        matrix[N2,A] lambdaByAgeByStrain[N_STRAINS];
        vector[N2] RtByStrain[N_STRAINS];
        matrix[A_REDUCED,A_REDUCED] reduced_flowsByStrain[ N_STRAINS, N_WEEKS_REDUCED_FLOWS[LOCATION_PROCESSING_IDX] ];
        matrix[A,A] full_flowsByStrain[ N_STRAINS, N_WEEKS_FULL_FLOWS[LOCATION_PROCESSING_IDX] ];
    
        // R0 new strain  
        R0_new_strain_m = R0[m] * rel_transmissibility_new_strain;
        R0_full_m = to_vector([R0[m], R0_new_strain_m]);
        
        // probability of transmission
        rho0 = R0_full_m / avg_cntct[m] ;
                
        impact_intv =
            country_impact(beta,
                dip_rnde[m],
                upswing_timeeff_reduced[,m],
                N2,
                A,
                A_CHILD,
                AGE_CHILD,
                COVARIATES_N,
                covariates[m],
                timeeff_shift_age[m,],
                upswing_timeeff_map[,m]
                );
        
          
          E_casesByAgeByStrain =
            country_EcasesByAge(
                R0_full_m,
                e_cases_N0[m],
                log_relsusceptibility_age,
                impact_intv_children_effect_gqs,
                impact_intv_onlychildren_effect_gqs,
                impact_intv,
                N0,
                elementary_school_reopening_idx[m],
                elt_reopening_index[m],
                N2,
                N_STRAINS,
                INTRO_NEW_STRAIN_TIME_IDX[m],
                prop_cases_new_strain_first_day,
                SCHOOL_STATUS[,m],
                A,
                A_CHILD,
                SI_CUT,
                wkend_idx[1:WKEND_IDX_N[m],m],
                avg_cntct[m],
                cntct_weekends_mean[m],
                cntct_weekdays_mean[m],
                cntct_school_closure_weekends[m],
                cntct_school_closure_weekdays[m],
                cntct_elementary_school_reopening_weekends[m],
                cntct_elementary_school_reopening_weekdays[m],
                rev_serial_interval,
                popByAge_abs[m,],
                N_init_A,
                init_A
                );
        
        E_casesByAge = E_casesByAgeByStrain[1] + E_casesByAgeByStrain[N_STRAINS];
        
        contacts_ByAge =
            country_contacts(
                impact_intv_children_effect_gqs,
                impact_intv_onlychildren_effect_gqs,
                impact_intv,
                N0,
                elementary_school_reopening_idx[m],
                elt_reopening_index[m],
                N2,
                SCHOOL_STATUS[,m],
                A,
                A_CHILD,
                wkend_idx[1:WKEND_IDX_N[m],m],
                cntct_weekends_mean[m],
                cntct_weekdays_mean[m],
                cntct_school_closure_weekends[m],
                cntct_school_closure_weekdays[m],
                cntct_elementary_school_reopening_weekends[m],
                cntct_elementary_school_reopening_weekdays[m]
                );
                
        lambdaByAgeByStrain =
            country_lambdaByAge(
                R0_full_m,
                log_relsusceptibility_age,
                impact_intv_children_effect_gqs,
                impact_intv_onlychildren_effect_gqs,
                impact_intv,
                E_casesByAgeByStrain,
                N0,
                elementary_school_reopening_idx[m],
                elt_reopening_index[m],
                N2,
                N_STRAINS,
                INTRO_NEW_STRAIN_TIME_IDX[m],
                prop_cases_new_strain_first_day,
                SCHOOL_STATUS[,m],
                A,
                A_CHILD,
                SI_CUT,
                wkend_idx[1:WKEND_IDX_N[m],m],
                avg_cntct[m],
                cntct_weekends_mean[m],
                cntct_weekdays_mean[m],
                cntct_school_closure_weekends[m],
                cntct_school_closure_weekdays[m],
                cntct_elementary_school_reopening_weekends[m],
                cntct_elementary_school_reopening_weekdays[m],
                rev_serial_interval
                );
                
 
        // generate R_ta
        RtByAgeByStrain =
            country_Rta(
                rho0,
                log_relsusceptibility_age,
                impact_intv_children_effect_gqs,
                impact_intv_onlychildren_effect_gqs,
                impact_intv,
                E_casesByAgeByStrain,
                elementary_school_reopening_idx[m],
                elt_reopening_index[m],
                N2,
                N_STRAINS,
                INTRO_NEW_STRAIN_TIME_IDX[m],
                SCHOOL_STATUS[,m],
                A,
                A_CHILD,
                wkend_idx[1:WKEND_IDX_N[m],m],
                cntct_weekends_mean[m],
                cntct_weekdays_mean[m],
                cntct_school_closure_weekends[m],
                cntct_school_closure_weekdays[m],
                cntct_elementary_school_reopening_weekends[m],
                cntct_elementary_school_reopening_weekdays[m],
                popByAge_abs[m,]
                );
        
        for(v in 1:N_STRAINS){
          
          reduced_flowsByStrain[v] =
            country_EflowsByLowDimAge(
                rho0[v],
                impact_intv_children_effect_gqs,
                impact_intv_onlychildren_effect_gqs,
                impact_intv,
                E_casesByAgeByStrain[v],
                log_relsusceptibility_age,
                elementary_school_reopening_idx[m],
                elt_reopening_index[m],
                SCHOOL_STATUS[,m],
                A,
                A_REDUCED,
                A_CHILD,
                SI_CUT,
                N_WEEKS_REDUCED_FLOWS[m],
                wkend_idx[1:WKEND_IDX_N[m],m],
                cntct_weekends_mean[m],
                cntct_weekdays_mean[m],
                cntct_school_closure_weekends[m],
                cntct_school_closure_weekdays[m],
                cntct_elementary_school_reopening_weekends[m],
                cntct_elementary_school_reopening_weekdays[m],
                rev_serial_interval,
                popByAge_abs[m,],
                reduced_age_bands_map,
                reduced_flows_Monday_idx[1:N_WEEKS_REDUCED_FLOWS[m], m],
                n_days
                );
   
        full_flowsByStrain[v] =
            country_EflowsByHighDimAge(
                rho0[v],
                impact_intv_children_effect_gqs,
                impact_intv_onlychildren_effect_gqs,
                impact_intv,
                E_casesByAgeByStrain[v],
                log_relsusceptibility_age,
                elementary_school_reopening_idx[m],
                elt_reopening_index[m],
                SCHOOL_STATUS[,m],
                A,
                A_CHILD,
                SI_CUT,
                N_WEEKS_FULL_FLOWS[m],
                wkend_idx[1:WKEND_IDX_N[m],m],
                cntct_weekends_mean[m],
                cntct_weekdays_mean[m],
                cntct_school_closure_weekends[m],
                cntct_school_closure_weekdays[m],
                cntct_elementary_school_reopening_weekends[m],
                cntct_elementary_school_reopening_weekdays[m],
                rev_serial_interval,
                popByAge_abs[m,],
                full_flows_Monday_idx[1:N_WEEKS_FULL_FLOWS[m], m],
                n_days
                );
                
          E_antibodyByAgeByStrain[v] =
            country_EantibodyByAge(
                E_casesByAgeByStrain[v],
                N2,
                A,
                rev_iar_daysSinceInfection);

          E_deathsByAgeByStrain[v] =
             country_EdeathsByAge(
                 E_casesByAgeByStrain[v],
                 N2,
                 A,
                 rev_ifr_daysSinceInfection,
                 log_ifr_age_base,
                 log_ifr_age_rnde_mid1[m],
                 log_ifr_age_rnde_mid2[m],
                 log_ifr_age_rnde_old[m]
                 );

          // effective cases
          E_effcasesByAgeByStrain[v] =
              country_EeffcasesByAge(
                  E_casesByAgeByStrain[v],
                  N0,
                  N2,
                  A,
                  SI_CUT,
                  rev_serial_interval
                  );
                  
          E_antibodyByStrain[v] = E_antibodyByAgeByStrain[v] * ones_vector_A;
        
          // generate Rt as weighted avg of R_ta
          RtByStrain[v] = RtByAgeByStrain[v] * popByAge[:,m];
        
          // generate total expected deaths
          E_deathsByStrain[v] = E_deathsByAgeByStrain[v] * ones_vector_A;
          
          // generate total expected cases
          E_casesByStrain[v] = E_casesByAgeByStrain[v] * ones_vector_A;
          
        }
        
        // sum over the two strains
        E_deathsByAge = E_deathsByAgeByStrain[1] + E_deathsByAgeByStrain[N_STRAINS];
        E_effcasesByAge = E_effcasesByAgeByStrain[1] + E_effcasesByAgeByStrain[N_STRAINS];
        E_antibodyByAge = E_antibodyByAgeByStrain[1] + E_antibodyByAgeByStrain[N_STRAINS];
        
        // sum over the two strains and age
        E_deaths = E_deathsByAge * ones_vector_A;
        E_cases = E_casesByAge * ones_vector_A;
        E_antibody = E_antibodyByAge * ones_vector_A;

    }
}







