#include "get_pbties.h"

#define DEBUG 0


double get_prob(double arg_tot)
{
    double prob=0.0;
    if(isinf(exp(arg_tot))==0){
        prob=exp(arg_tot)/(1.0+exp(arg_tot));
    }
    else{
        prob=1.0;
    }
    return prob;
}


double get_intercept_prob(double mu)
{
    // get the total of the arguments
    double arg_tot = mu;
    double prob = get_prob(arg_tot);

    return prob;
}


double get_I_W(Double_Matrices_cont mat_covars,
                            unsigned int ind,
                            double mu,
                            double lambda_1,
                            double lambda_2,
                            double lambda_3,
                            double lambda_4,
                            double lambda_5,
                            double lambda_6,
                            double lambda_7,
                            double lambda_8,
                            double lambda_9,
                            double lambda_10,
                            double lambda_11,
                            double lambda_12)
{
    // get the time dependent and time independent covariates
    
    double age = mat_covars.matrix[ind*mat_covars.nb_columns+0];
    double sex = mat_covars.matrix[ind*mat_covars.nb_columns+1];
    double ethnic_other = mat_covars.matrix[ind*mat_covars.nb_columns+2];
    double ethnic_black = mat_covars.matrix[ind*mat_covars.nb_columns+3];
    double edu_intermediate = mat_covars.matrix[ind*mat_covars.nb_columns+9];
    double edu_low = mat_covars.matrix[ind*mat_covars.nb_columns+10];
    double flat = mat_covars.matrix[ind*mat_covars.nb_columns+11];
    double bmi = mat_covars.matrix[ind*mat_covars.nb_columns+4];
    double ss_prev = mat_covars.matrix[ind*mat_covars.nb_columns+5];
    double ss_curr = mat_covars.matrix[ind*mat_covars.nb_columns+6];
    double as_prev = mat_covars.matrix[ind*mat_covars.nb_columns+7];
    double as_curr = mat_covars.matrix[ind*mat_covars.nb_columns+8];

    // get the total of the arguments
    double arg_tot = mu + lambda_1*age + lambda_2*sex + lambda_3*ethnic_other + lambda_4*ethnic_black +
                    lambda_5*edu_intermediate + lambda_6*edu_low + lambda_7*flat +
                    lambda_8*bmi + lambda_9*ss_prev + lambda_10*ss_curr + lambda_11*as_prev + lambda_12*as_curr;

    double prob = get_prob(arg_tot);

    // cout << " ind " << ind << " day " << day << " time_cov " << time_cov << " cov " << cov <<
    // " mu " << mu << " lambda_1 " << lambda_1 << " lambda_2 " << lambda_2 << " lambda_3 " << lambda_3 <<
    // " arg_tot " << arg_tot << " prob " << prob << endl;

    return prob;
}


void get_from_I(Double_Matrices_cont mat_P_I_W,
                Double_Matrices_cont mat_P_I_I,
                Double_Matrices_cont mat_P_I_R,
                Double_Matrices_cont mat_covars,
                gsl_vector* theta_temp,
                std::vector<int> start_indicator,
                std::vector<int> final_indicator)
{
    double prob;
    unsigned int n_ind=mat_P_I_W.nb_rows;
    unsigned int n_days=mat_P_I_W.nb_columns;

    double mu_0 = theta_temp->data[0];
    double lambda_1_0 = theta_temp->data[1];
    double lambda_2_0 = theta_temp->data[2];
    double lambda_3_0 = theta_temp->data[3];
    double lambda_4_0 = theta_temp->data[4];
    double lambda_5_0 = theta_temp->data[5];
    double lambda_6_0 = theta_temp->data[6];
    double lambda_7_0 = theta_temp->data[7];
    double lambda_8_0 = theta_temp->data[8];
    double lambda_9_0 = theta_temp->data[9];
    double lambda_10_0 = theta_temp->data[10];
    double lambda_11_0 = theta_temp->data[11];
    double lambda_12_0 = theta_temp->data[12];

    double mu_1 = theta_temp->data[13];

    for(unsigned int ind=0;ind<n_ind;ind++){
        
        double cum_prob = 1;

        // Get the I-W transition
        prob=cum_prob;
        prob=prob*get_I_W(mat_covars,
                            ind,
                            mu_0,
                            lambda_1_0,
                            lambda_2_0,
                            lambda_3_0,
                            lambda_4_0,
                            lambda_5_0,
                            lambda_6_0,
                            lambda_7_0,
                            lambda_8_0,
                            lambda_9_0,
                            lambda_10_0,
                            lambda_11_0,
                            lambda_12_0);
        cum_prob-=prob;

        // cout << "Individual " << ind << " first_day " << first_day << " final_day " << final_day << endl;
        for(unsigned int day=0;day<n_days;day++){mat_P_I_W.matrix[ind*mat_P_I_W.nb_columns+day]=prob;}

        // Get the I-R transition
        prob=cum_prob;
        prob=prob*get_intercept_prob(mu_1);
        cum_prob-=prob;
        for(unsigned int day=0;day<n_days;day++){mat_P_I_R.matrix[ind*mat_P_I_R.nb_columns+day]=prob;}

        // Get the I-I transition
        prob=cum_prob;
        for(unsigned int day=0;day<n_days;day++){mat_P_I_I.matrix[ind*mat_P_I_I.nb_columns+day]=prob;}
        
    }
}


void get_from_W(Double_Matrices_cont mat_P_W_D,
                Double_Matrices_cont mat_P_W_R,
                Double_Matrices_cont mat_P_W_W,
                Double_Matrices_cont mat_covars,
                gsl_vector* theta_temp,
                std::vector<int> start_indicator,
                std::vector<int> final_indicator)
{
    double prob;
    unsigned int n_ind=mat_P_W_D.nb_rows;
    unsigned int n_days=mat_P_W_D.nb_columns;

    double mu_0 = theta_temp->data[0];
    double lambda_1_0 = theta_temp->data[1];
    double lambda_2_0 = theta_temp->data[2];
    double lambda_3_0 = theta_temp->data[3];
    double lambda_4_0 = theta_temp->data[4];
    double lambda_5_0 = theta_temp->data[5];
    double lambda_6_0 = theta_temp->data[6];
    double lambda_7_0 = theta_temp->data[7];
    double lambda_8_0 = theta_temp->data[8];
    double lambda_9_0 = theta_temp->data[9];
    double lambda_10_0 = theta_temp->data[10];
    double lambda_11_0 = theta_temp->data[11];
    double lambda_12_0 = theta_temp->data[12];

    double mu_1 = theta_temp->data[13];

    for(unsigned int ind=0;ind<n_ind;ind++){
        
        double cum_prob = 1;

        // Get the I-W transition
        prob=cum_prob;
        prob=prob*get_I_W(mat_covars,
                            ind,
                            mu_0,
                            lambda_1_0,
                            lambda_2_0,
                            lambda_3_0,
                            lambda_4_0,
                            lambda_5_0,
                            lambda_6_0,
                            lambda_7_0,
                            lambda_8_0,
                            lambda_9_0,
                            lambda_10_0,
                            lambda_11_0,
                            lambda_12_0);
        cum_prob-=prob;

        // cout << "Individual " << ind << " first_day " << first_day << " final_day " << final_day << endl;
        for(unsigned int day=0;day<n_days;day++){mat_P_W_D.matrix[ind*mat_P_W_D.nb_columns+day]=prob;}

        // Get the I-R transition
        prob=cum_prob;
        prob=prob*get_intercept_prob(mu_1);
        cum_prob-=prob;
        for(unsigned int day=0;day<n_days;day++){mat_P_W_R.matrix[ind*mat_P_W_R.nb_columns+day]=prob;}

        // Get the I-I transition
        prob=cum_prob;
        for(unsigned int day=0;day<n_days;day++){mat_P_W_W.matrix[ind*mat_P_W_W.nb_columns+day]=prob;}
        
    }
}


double calculate_likelihood(double likelihood,
                            Double_Matrices_cont mat_P,
                            std::vector<std::vector<std::vector<int> > > tensor_indicator,
                            unsigned int trans)
{
    unsigned int n_ind=mat_P.nb_rows;
    double prob;
    
    // might want to swap this around for speed - will make debugging more difficult though
    for(unsigned int ind=0;ind<n_ind;ind++){
        unsigned int trans_days=tensor_indicator[trans][ind].size();
        if(trans_days>1){
            std::vector<int> ind_indicator = tensor_indicator[trans][ind];
            unsigned int first_day = ind_indicator[1];
            unsigned int final_day = ind_indicator[trans_days-1];
            for(unsigned int day=first_day;day<=final_day;day++){
                prob=mat_P.matrix[ind*mat_P.nb_columns+day];
                double logP=log(prob);
                likelihood+=logP;

                // cout << "Individual: " << ind << " prob: " << prob <<
                // " days: " << day << " logP: " << logP << " cumulative likelihood: " << likelihood << endl;
            }

        }

    }
    return likelihood;
}


double wrapped_likelihood(std::vector<std::vector<std::vector<int> > >  tensor_indicator,
                          std::vector<int> start_indicator,
                          std::vector<int> final_indicator,
                          Double_Matrices_cont mat_covars,
                          Double_Matrices_cont mat_P_I_W,
                          Double_Matrices_cont mat_P_I_I,
                          Double_Matrices_cont mat_P_I_R,
                          Double_Matrices_cont mat_P_W_D,
                          Double_Matrices_cont mat_P_W_R,
                          Double_Matrices_cont mat_P_W_W,
                          gsl_vector* theta)

{
    // Get the from I transitions
    int first_param_i=0;
    int last_param_i=14;
    gsl_vector* theta_i=gsl_vector_alloc(last_param_i-first_param_i);
    int i=0;
    for(int j=first_param_i;j<last_param_i;j++){theta_i->data[i]=theta->data[j];i++;}
    get_from_I(mat_P_I_W,
            mat_P_I_I,
            mat_P_I_R,
            mat_covars,
            theta_i,
            start_indicator,
            final_indicator);
    
    // Get the from W transitions
    int first_param_w=14;
    int last_param_w=28;
    gsl_vector* theta_w=gsl_vector_alloc(last_param_w-first_param_w);
    int w=0;
    for(int j=first_param_w;j<last_param_w;j++){theta_w->data[w]=theta->data[j];w++;}
    get_from_W(mat_P_W_D,
               mat_P_W_R,
               mat_P_W_W,
               mat_covars,
               theta_w,
               start_indicator,
               final_indicator);

    double current_L = 0;
    current_L=calculate_likelihood(current_L, mat_P_I_I, tensor_indicator, 0);
    current_L=calculate_likelihood(current_L, mat_P_I_W, tensor_indicator, 1);
    current_L=calculate_likelihood(current_L, mat_P_I_R, tensor_indicator, 2);
    current_L=calculate_likelihood(current_L, mat_P_W_W, tensor_indicator, 3);
    current_L=calculate_likelihood(current_L, mat_P_W_R, tensor_indicator, 4);
    current_L=calculate_likelihood(current_L, mat_P_W_D, tensor_indicator, 5);;

    return current_L;
}


unsigned int candidate_acceptance(double current_L,
                                  double previous_L)
{
    unsigned int accept=0;
    if(isinf(current_L)==0){
        if(current_L>previous_L){
            if(DEBUG==1){
                cout << "\t\tHigher L, Candidate accepted" << endl;
            }
            accept=1;
        }
        else{
            double acc_ratio=exp(current_L-previous_L);
            double rand_test=myrand();
            if(DEBUG==1){
                cout << "\t\tLower L, likelihood ratio " << acc_ratio << " -- test_sample " << rand_test << endl;
            }
            if(rand_test<=acc_ratio){
                accept=1;
                if(DEBUG==1){
                    cout << "\t\tCandidate accepted" << endl;
                }
            }//end of if (test<ratio)
            else{
                if(DEBUG==1){
                    cout << "\t\tCandidate rejected" << endl;
                }
            }//end of else (test>ratio)
        }//end of else (ratio)
    }
    return accept;
}


unsigned int metropolis_step_multiv(std::vector<std::vector<std::vector<int> > >  tensor_indicator,
                                    std::vector<int> start_indicator,
                                    std::vector<int> final_indicator,
                                    Double_Matrices_cont mat_history,
                                    Double_Matrices_cont mat_covars,
                                    Double_Matrices_cont mat_P_I_W,
                                    Double_Matrices_cont mat_P_I_I,
                                    Double_Matrices_cont mat_P_I_R,
                                    Double_Matrices_cont mat_P_W_D,
                                    Double_Matrices_cont mat_P_W_R,
                                    Double_Matrices_cont mat_P_W_W,
                                    gsl_vector * theta_previous,
                                    gsl_matrix * L,
                                    unsigned int iter,
                                    unsigned int n_params)
{
    // Number of parameters
    
    if(DEBUG==1){
        cout<<"Previous theta:"<<endl;
        display_gsl_vector(theta_previous);
    }
    
    if(DEBUG==1){
        cout<<"Current L (i.e. Cholesky decomposition of current Sigma):"<<endl;
        display_gsl_matrix(L);
    }

    // Step 1.1: Sampling candidate for parameter: random walk

    static gsl_rng *RandomNumberGenerator = gsl_rng_alloc( gsl_rng_taus );
    gsl_vector *theta=gsl_vector_alloc(n_params);
    gsl_ran_multivariate_gaussian(RandomNumberGenerator, theta_previous, L, theta);

    // Set to zero for testing
    // theta->data[0] = 0;
    // theta->data[1] = 0;
    // theta->data[2] = 0;
    // theta->data[3] = 0;
    // theta->data[4] = 0;

    if(DEBUG==1){
        cout<<"Generated sample:"<<endl;
        display_gsl_vector(theta);
    }
    
    // Step 1.2: Computing the likelihood for the candidate mu
    double current_L=wrapped_likelihood(tensor_indicator,
                    start_indicator,
                    final_indicator,
                    mat_covars,
                    mat_P_I_W,
                    mat_P_I_I,
                    mat_P_I_R,
                    mat_P_W_D,
                    mat_P_W_R,
                    mat_P_W_W,
                    theta);

    if(DEBUG==1){
        ostringstream current_L_str;
        current_L_str << current_L;
        string current = "Current L " + current_L_str.str();
        // unsigned int n = 0;
        // for(unsigned int trans=0;trans<n_trans;trans++){
        //     for(unsigned int param=0;param<n_params;param++){
        //         ostringstream trans_str;
        //         ostringstream param_str;
        //         ostringstream theta_str;
        //         trans_str << trans;
        //         param_str << param;
        //         theta_str << theta->data[n];

        //         if (param == 0){
        //             current += " -- Mu_" + trans_str.str() + " " + theta_str.str();
        //         }
        //         else {
        //             current += " -- Lambda" + param_str.str() + "_" + trans_str.str() + " " + theta_str.str();
        //         }
        //         n++;
        //     }
        // }
        cout << current << endl;
    }


    // Step 1.3: Accepting/Rejecting the candidate
    
    double previous_L=mat_history.matrix[(iter-1)*mat_history.nb_columns+1];

    unsigned int accept_candidate=0;
    accept_candidate=candidate_acceptance(current_L,
                                          previous_L);

    if(DEBUG==1){
        cout << "Previous L: " << previous_L << endl;
        cout << "\taccept_candidate " <<  accept_candidate
        << endl;
    }


    // Storing values in history matrix

    if(accept_candidate==1){
        mat_history.matrix[iter*mat_history.nb_columns+1]=current_L;
        for(unsigned int i=0;i<n_params;i++){
            mat_history.matrix[iter*mat_history.nb_columns+2+i]=theta->data[i];
            // cout << theta->data[i] << endl;
        }
    } 
    else {
        mat_history.matrix[iter*mat_history.nb_columns+1]=previous_L;
        for(unsigned int i=0;i<n_params;i++){
            mat_history.matrix[iter*mat_history.nb_columns+2+i]=theta_previous->data[i];
            // cout << theta->data[i] << endl;
        }
    }

    return accept_candidate;
}


double adaptive_step_multiv(Double_Matrices_cont mat_history,
                            gsl_matrix* hat_sigma,
                            unsigned int iter,
                            unsigned int n_iter_adaptive,
                            unsigned int N)
{
    if(DEBUG==1){
        cout << "Updating covariance of the proposal -- iter " << iter << endl;
        cout << "Number of iterations in batch (always the same): " << n_iter_adaptive << endl;
    }

    // Storing values of theta over "n_iter_adaptive" previous iterations in X

    gsl_matrix* X=gsl_matrix_alloc(n_iter_adaptive, N);
    for (unsigned int i=0;i<n_iter_adaptive;i++){
        for(unsigned int k=0;k<N;k++){
            X->data[i * X->tda + k]=mat_history.matrix[(iter-i)*mat_history.nb_columns+2+k];
        }
    }

    if(DEBUG==1){
        cout << "Matrix X:" << endl;
        display_gsl_matrix(X);
    }

    // Computing the covariance

    gsl_ran_multivariate_gaussian_vcov(X, hat_sigma);

    if(DEBUG==1){
        cout << "Empirical covariance:" << endl;
        display_gsl_matrix(hat_sigma);
    }
    
    return 0.0;
}

