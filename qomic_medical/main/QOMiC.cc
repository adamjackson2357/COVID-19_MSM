#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <iostream>
#include "xml_file_write.h"
#include "xml_file_read.h"
#include <vector>
#include <algorithm>
#include <memory>
#include <fstream>
#include <../Routines/struc.h>
#include <../Routines/manip_dates.h>
#include <../Routines/dyn_name.h>
#include <../Routines/matrix_handling.h>
#include <../Routines/rand.h>
#include "../Routines/get_input.h"
#include "../Routines/get_pbties.h"
#include <../Routines/xml_file_read.h>
#include <../Classes/String_Matrices.h>
#include <../Classes/Int_Matrices_cont.h>
#include <../Classes/Double_Matrices_cont.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

#define DEBUG 0

using namespace std;


// File Names
char nameIN[256]; // input file called by switch -i
char nameOUT[256]; // output  file called by switch -o


//******************************************************************
//*main
//******************************************************************


int main(int argc, char *  argv[])
{
    int na=0;
    char filename_in_IDs[1000];
    char filename_in_dates[1000];
    char filename_in_status[1000];
    char filename_in_covars[1000];
    char fic_out_MCMC[1000];
    char filename_in_trans[1000];
    char filename_in_theta[1000];
    char filename_in_sigma[1000];
    
    
    long MY_SEED=0;


    // Initialisation of parameters (default values)

    unsigned int n_iter=0;
    unsigned int burn_in=0;

    unsigned int adaptive=1;
    unsigned int n_iter_adaptive=100;
    unsigned int update_frequency=100;
    double sigma_mu=0.0;

    // Number of parameters and transitions
    unsigned int n_params=5;
    unsigned int n_trans=3;

    na++;
    
    // Reading parameters from command line
    
    while(na < argc){
        if ( 0 == strcmp(argv[na],"-ID") ){
            strcpy(filename_in_IDs,argv[++na]);
            if (na+1==argc) break;
            na++;
        }
        if ( 0 == strcmp(argv[na],"-dates") ){
            strcpy(filename_in_dates,argv[++na]);
            if (na+1==argc) break;
            na++;
        }
        else if ( 0 == strcmp(argv[na],"-status") ){
            strcpy(filename_in_status,argv[++na]);
            if (na+1==argc) break;
            na++;
        }
        else if ( 0 == strcmp(argv[na],"-covars") ){
            strcpy(filename_in_covars,argv[++na]);
            if (na+1==argc) break;
            na++;
        }
        else if ( 0 == strcmp(argv[na],"-trans") ){
            strcpy(filename_in_trans,argv[++na]);
            if (na+1==argc) break;
            na++;
        }
        else if ( 0 == strcmp(argv[na],"-theta") ){
            strcpy(filename_in_theta,argv[++na]);
            if (na+1==argc) break;
            na++;
        }
        else if ( 0 == strcmp(argv[na],"-sigma") ){
            strcpy(filename_in_sigma,argv[++na]);
            if (na+1==argc) break;
            na++;
        }
        else if ( 0 == strcmp(argv[na],"-seed") ){
            MY_SEED=(long)((atoi(argv[++na])));
            if (na+1==argc) break;
            na++;
        }
        else if ( 0 == strcmp(argv[na],"-iter") ){
            n_iter=atoi(argv[++na]);
            if (na+1==argc) break;
            na++;
        }
        else if ( 0 == strcmp(argv[na],"-burn_in") ){
            burn_in=atoi(argv[++na]);
            if (na+1==argc) break;
            na++;
        }
        else if ( 0 == strcmp(argv[na],"-adaptive") ){
            adaptive=atoi(argv[++na]);
            if (na+1==argc) break;
            na++;
        }
        else if ( 0 == strcmp(argv[na],"-iter_adaptive") ){
            n_iter_adaptive=atoi(argv[++na]);
            if (na+1==argc) break;
            na++;
        }
        else if ( 0 == strcmp(argv[na],"-update_frequency") ){
            update_frequency=atoi(argv[++na]);
            if (na+1==argc) break;
            na++;
        }
        else if ( 0 == strcmp(argv[na],"-n_params") ){
            n_params=atoi(argv[++na]);
            if (na+1==argc) break;
            na++;
        }
        else if ( 0 == strcmp(argv[na],"-n_trans") ){
            n_trans=atoi(argv[++na]);
            if (na+1==argc) break;
            na++;
        }
        else if ( 0 == strcmp(argv[na],"-out") ){
            strcpy(fic_out_MCMC,argv[++na]);
            if (na+1==argc) break;
            na++;
        }
        else{
            cout << "Unknown option: " << argv[na] << endl;
            exit(1);
        }
    }

    // Reading transition matrix   
    Int_Matrices_cont mat_col;
    mat_col.Read_from_file(filename_in_trans);
    unsigned int n_from = mat_col.nb_rows;
    unsigned int n_to = mat_col.nb_columns;

    cout << "n_from " << n_from << endl;
    cout << "n_to " << n_to << endl;

    cout << "Col matrix" << endl;
    mat_col.Display_matrix();

    // Get theta vector
    Double_Matrices_cont mat_theta;
    mat_theta.Read_from_file(filename_in_theta);
    gsl_matrix *mat_theta_work=Double_matrices_cont_2_gsl_matrix(mat_theta);
    gsl_vector *theta_init = get_one_vect(0, mat_theta_work);
    mat_theta.Free_double_matrix_cont();
    gsl_matrix_free(mat_theta_work);

    cout<<"Initial parameters theta:"<<endl;
    display_gsl_vector(theta_init);

    // Dynamic writing of output files

    // Dynamic file column outputs
    string history_cols="Iter\tL";
    string posterior_cols="BurnIn\tpost_L";

    std::vector<std::string> theta_string = 
            {"mu_0",
             "lambda_1_0", "lambda_2_0", "lambda_3_0", "lambda_4_0", "lambda_5_0", "lambda_6_0",
             "lambda_7_0", "lambda_8_0", "lambda_9_0", "lambda_10_0", "lambda_11_0", "lambda_12_0",
             "lambda_13_0", "lambda_14_0", "lambda_15_0", "lambda_16_0", "lambda_17_0",
             "lambda_18_0", "lambda_19_0", "lambda_20_0", "lambda_21_0", "lambda_22_0",
            "mu_1",
            "mu_2",
             "lambda_1_2", "lambda_2_2", "lambda_3_2", "lambda_4_2", "lambda_5_2", "lambda_6_2",
             "lambda_7_2", "lambda_8_2", "lambda_9_2", "lambda_10_2", "lambda_11_2", "lambda_12_2",
             "lambda_13_2", "lambda_14_2", "lambda_15_2", "lambda_16_2", "lambda_17_2",
             "lambda_18_2", "lambda_19_2", "lambda_20_2", "lambda_21_2", "lambda_22_2",
            "mu_3"};
    
    std::string param_cols="";
    for(unsigned int i=0;i<theta_string.size();i++){
        std::string temp_string = "\t"+theta_string[i];
        param_cols+=temp_string;
    }

    std::string sigma_cols="";
    for(unsigned int i=0;i<theta_string.size();i++){
        for(unsigned int j=0;j<theta_string.size();j++){
            std::string temp_string = "\tsigma_" + theta_string[i] + "_" + theta_string[j];
            sigma_cols+=temp_string;
        }
    }

    history_cols=history_cols+param_cols+sigma_cols;
    posterior_cols=posterior_cols+param_cols;

    // History file
    string mcmc_output_fname=write_dynamic_fname(fic_out_MCMC, "history.txt");
    ofstream f_out;
    f_out.open(mcmc_output_fname.c_str(),ios::out);
    if(f_out.fail()) {
        cout << "Invalid Path and/or permission rights for " << mcmc_output_fname << " -- run stopped." << endl;
        exit(1);
    }
    else {
        f_out << history_cols <<endl;
    } 

    // Posterior likelihood file
    string posterior_l_output_fname=write_dynamic_fname(fic_out_MCMC, "posterior_l.txt");
    ofstream f_outL;
    f_outL.open(posterior_l_output_fname.c_str(),ios::out);
    if(f_outL.fail()){
        cout << "Invalid Path and/or permission rights for " << posterior_l_output_fname << " -- run stopped." << endl;
        exit(1);
    }
    else {
        f_outL << posterior_cols << endl;
    }

    smyrand((long)(MY_SEED));

    // Summary of the adaptive parameters:
    
    if (adaptive==1){
        cout << "Adaptive algorithm: " << endl;
        cout << "Number of iterations for adaptive: " << n_iter_adaptive << endl;
        cout << "Update frequency: " << update_frequency << endl;
    }
    
    // Reading dates file    
    Double_Matrices_cont mat_dates;
    mat_dates.Read_from_file(filename_in_dates);
    gsl_matrix *mat_dates_temp=Double_matrices_cont_2_gsl_matrix(mat_dates);
    unsigned int n_days = gsl_matrix_max(mat_dates_temp)+1;
    gsl_matrix_free(mat_dates_temp);

    // Reading status file
    Int_Matrices_cont mat_status;
    mat_status.Read_from_file(filename_in_status);

    // Get the number of individuals
    unsigned int n_ind=mat_status.nb_rows;

    cout << "n_days_follow_up= " << n_days << " -- n_ind " << n_ind << endl;

    // Reading covariates file
    Double_Matrices_cont mat_covars;
    mat_covars.Read_from_file(filename_in_covars);
    
    // Calculating state matrices from input
    Int_Matrices_cont mat_from;
    mat_from.Alloc_int_matrix_cont_val(n_ind, n_to, -1);
    get_from_matrix(mat_from, mat_status, mat_dates, n_ind, n_to);

    // cout << "From matrix" << endl;
    // mat_from.Display_matrix();

    Int_Matrices_cont mat_to;
    mat_to.Alloc_int_matrix_cont_val(n_ind, n_to, n_days);
    get_to_matrix(mat_to, mat_status, mat_dates, n_ind, n_to);

    // cout << "To matrix" << endl;
    // mat_to.Display_matrix();

    Int_Matrices_cont mat_pos_from;
    mat_pos_from.Alloc_int_matrix_cont_val(n_ind, n_to, -1);
    get_pos_from_matrix(mat_pos_from, mat_status, n_ind, n_to);

    // cout << "Pos from matrix" << endl;
    // mat_pos_from.Display_matrix();

    Int_Matrices_cont mat_pos_to;
    mat_pos_to.Alloc_int_matrix_cont_val(n_ind, n_to, -1);
    get_pos_to_matrix(mat_pos_to, mat_status, n_ind, n_to);

    // cout << "Pos to matrix" << endl;
    // mat_pos_to.Display_matrix();

    // Create the tensor indicator
    cout << "Initialise indicator" << endl;
    std::vector<std::vector<std::vector<int> > > tensor_indicator;

    cout << "Fill with -1s" << endl;
    initialise_tensor_indicator(tensor_indicator, n_trans, n_ind);

    cout << "Create indicator" << endl;
    create_tensor_indicator(tensor_indicator, mat_from, mat_to, mat_pos_from,
                            mat_pos_to, mat_col, n_ind, n_from);

    // cout << "Indicator tensor" << endl;
    // summary_tensor_indicator(tensor_indicator, n_ind, n_trans);

    // start indicator
    // to save time start at the day of infection
    std::vector<int> start_indicator;
    create_start_indicator(tensor_indicator, start_indicator, n_ind, n_days, n_trans);
    // cout << "Start indicator" << endl;
    // summary_vector_indicator(start_indicator, n_ind);

    // end indicator 
    // to save time end at the final day for the individual
    std::vector<int> final_indicator;
    create_final_indicator(tensor_indicator, final_indicator, n_ind, n_trans);
    // cout << "Final indicator" << endl;
    // summary_vector_indicator(final_indicator, n_ind);

    mat_col.Free_int_matrix_cont();
    mat_from.Free_int_matrix_cont();
    mat_to.Free_int_matrix_cont();
    mat_pos_from.Free_int_matrix_cont();
    mat_pos_to.Free_int_matrix_cont();

    // Initialisation of transition probabilities
    Double_Matrices_cont mat_P_I_W;
    mat_P_I_W.Alloc_double_matrix_cont(n_ind, n_days);

    Double_Matrices_cont mat_P_I_I;
    mat_P_I_I.Alloc_double_matrix_cont(n_ind, n_days);

    Double_Matrices_cont mat_P_I_R;
    mat_P_I_R.Alloc_double_matrix_cont(n_ind, n_days);

    Double_Matrices_cont mat_P_W_D;
    mat_P_W_D.Alloc_double_matrix_cont(n_ind, n_days);

    Double_Matrices_cont mat_P_W_R;
    mat_P_W_R.Alloc_double_matrix_cont(n_ind, n_days);

    Double_Matrices_cont mat_P_W_W;
    mat_P_W_W.Alloc_double_matrix_cont(n_ind, n_days);

    // Number of parameters

    int q=n_params*n_params; // Total covariance of all the parameters
    int ncol=1+1+n_params+q; // iter+L+params
    cout << "Number of columns in history matrix: "<< ncol << endl;

    // Matrix in which all visited parameters values and likelihoods will be stored
    Double_Matrices_cont mat_history;
    mat_history.Alloc_double_matrix_cont(n_iter,ncol);

    // Initialisation of current parameters from configs

    cout<<"Initial parameters theta:"<<endl;
    display_gsl_vector(theta_init);

    // Get Sigma matrix
    Double_Matrices_cont mat_sigma;
    mat_sigma.Read_from_file(filename_in_sigma);
    gsl_matrix* Sigma=Double_matrices_cont_2_gsl_matrix(mat_sigma);

    cout<<"Initial covariance Sigma:"<<endl;
    display_gsl_matrix(Sigma);

    // Cholesky decomposition of the covariance Sigma

    gsl_matrix* L=gsl_matrix_calloc(n_params,n_params);
    for(unsigned int i=0;i<n_params*n_params;i++){
        L->data[i]=Sigma->data[i];
    }

    cout<<"After Cholesky decomposition:"<<endl;
    gsl_linalg_cholesky_decomp(L);
    display_gsl_matrix(L);

    double cd=pow(2.4,2)/n_params;
    cout << "Scaling factor: " << cd << endl;

    // Initialisation of the MCMC
    
    // Computing the likelihood with initial parameters

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
                        theta_init);

    // cout << "P_I_I matrix" << endl;
    // mat_P_I_I.Display_matrix();
    // cout << "P_I_W matrix" << endl;
    // mat_P_I_I.Display_matrix();
    // cout << "P_I_R matrix" << endl;
    // mat_P_I_R.Display_matrix();
    // cout << "P_W_W matrix" << endl;
    // mat_P_W_W.Display_matrix();
    // cout << "P_W_C matrix" << endl;
    // mat_P_W_C.Display_matrix();
    // cout << "P_W_R matrix" << endl;
    // mat_P_W_R.Display_matrix();
    // cout << "P_W_D matrix" << endl;
    // mat_P_W_D.Display_matrix();
    // cout << "P_C_C matrix" << endl;
    // mat_P_C_C.Display_matrix();
    // cout << "P_C_R matrix" << endl;
    // mat_P_C_R.Display_matrix();
    // cout << "P_C_D matrix" << endl;
    // mat_P_C_D.Display_matrix();

    cout << "Starting likelihood " << current_L << endl;

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

    // // Storing parameters values at current iteration

    cout << "Filling iter and current_L:" << endl;   
    mat_history.matrix[0*mat_history.nb_columns]=0;
    mat_history.matrix[0*mat_history.nb_columns+1]=current_L;

    cout << "Filling parameter values:" << endl;
    for(unsigned int i=0;i<n_params;i++){
        mat_history.matrix[0*mat_history.nb_columns+2+i]=theta_init->data[i];
        cout << theta_init->data[i] << endl;
    }

    cout << "Filling covariance values:" << endl;
    for(unsigned int i=0;i<n_params*n_params;i++){
        mat_history.matrix[0*mat_history.nb_columns+2+n_params+i]=Sigma->data[i];
        // cout << Sigma->data[i] << endl;
    }

    // Writing parameters values in dynamic output file
    
    for(unsigned int col=0;col<mat_history.nb_columns;col++){
        f_out << mat_history.matrix[0*mat_history.nb_columns+col]
        << "\t";
    }
    f_out << endl;


    // Initialisation of acceptance counts for each parameter (for computation of acceptance probabilities at the end of the run)
    
    unsigned int nb_accepted_theta=0;
    unsigned int nb_accepted_theta_batch=0;
    unsigned int nb_accepted_theta_after_BI=0;

    unsigned int nb_accepted_mu_batch=0;

    // Starting MCMC algorithm

    for(unsigned int iter=1;iter<n_iter;iter++){
        // Writing iteration ID

        mat_history.matrix[iter*mat_history.nb_columns]=iter;

        if(DEBUG==1){
            cout << "**************" << endl
            << "Iter= " << iter << endl
            << "**************" << endl;
        }
        else{
            if(iter%100==0){
                cout << "**************" << endl
                << "Iter= " << iter << endl
                << "**************" << endl;
            }
        }

        // Storing parameters in the vector theta
        gsl_vector* theta_previous=gsl_vector_alloc(n_params);
        for(unsigned int i=0;i<n_params;i++){
            theta_previous->data[i]=mat_history.matrix[(iter-1)*mat_history.nb_columns+2+i];
        }

        // Metropolis step: sampling new candidate, computing its likelihood and accepting/rejecting it

        unsigned int accept_candidate=metropolis_step_multiv(tensor_indicator,
                                    start_indicator,
                                    final_indicator,
                                    mat_history,
                                    mat_covars,
                                    mat_P_I_W,
                                    mat_P_I_I,
                                    mat_P_I_R,
                                    mat_P_W_D,
                                    mat_P_W_R,
                                    mat_P_W_W,
                                    theta_previous,
                                    L,
                                    iter,
                                    n_params);

        // Updating acceptance count for this parameter
        
        if(accept_candidate==1){
            nb_accepted_theta++;
            nb_accepted_theta_batch++;
            if(iter>burn_in){
                nb_accepted_theta_after_BI++;
            }
        }


        // Storing the covariance values
        
        for(unsigned int i=0;i<n_params*n_params;i++){
            mat_history.matrix[iter*mat_history.nb_columns+2+n_params+i]=Sigma->data[i];
        }

        // Updating scaling of the proposal every n_iter_adaptive iterations

        if (adaptive==1){
            if(iter%update_frequency==0){
                gsl_matrix* hat_sigma=gsl_matrix_alloc(n_params,n_params);
                sigma_mu=adaptive_step_multiv(mat_history,
                                              hat_sigma,
                                              iter,
                                              n_iter_adaptive,
                                              n_params);

                // Scaling factor
                
                for (unsigned int i=0;i<(n_params*n_params);i++){
                    L->data[i]=cd*hat_sigma->data[i];
                    Sigma->data[i]=cd*hat_sigma->data[i];
                }

                // Adding a constant to ensure the matrix is positive definite
                
                for (unsigned int i=0;i<n_params;i++){
                    L->data[i * L->tda + i]+=0.00001;
                    Sigma->data[i * Sigma->tda + i]+=0.00001;
                }

                if (DEBUG==1){
                    cout << "Updated covariance:" << endl;
                    display_gsl_matrix(Sigma);
                }
                 
                gsl_linalg_cholesky_decomp(L);
                
                if (DEBUG==1){
                    cout<<"After Cholesky decomposition:"<<endl;
                    display_gsl_matrix(L);
                }
                
                nb_accepted_mu_batch=0;
            }

        }
        
        for(unsigned int col=0;col<mat_history.nb_columns;col++){
            f_out << mat_history.matrix[iter*mat_history.nb_columns+col]
            << "\t";
        }

        // Close line in dynamic history file
        f_out << endl;
    
    }//end of for iter
    
    // Storing parameters in the vector theta
    gsl_vector* theta_final=gsl_vector_alloc(n_params);
    for(unsigned int i=0;i<n_params;i++){
        theta_final->data[i]=mat_history.matrix[(n_iter-1)*mat_history.nb_columns+2+i];
    }
    
    cout << "//////////////////////////////////////////////////////////////////" << endl
    << "//////////////////////////////////////////////////////////////////" << endl
    << "                      MCMC estimation done" << endl
    << "\tAcceptance rate: " << 100.0*(double)(nb_accepted_theta)/(double)(n_iter-1) << endl
    << "\tAcceptance rate (after BI): " << 100.0*(double)(nb_accepted_theta_after_BI)/(double)(n_iter-burn_in-1) << endl;
    unsigned int n = 0;
    for(unsigned int trans=0;trans<n_trans;trans++){
        for(unsigned int param=0;param<n_params;param++){
            ostringstream trans_str;
            ostringstream param_str;
            ostringstream theta_str;
            trans_str << trans;
            param_str << param;
            theta_str << theta_final->data[n];
            if (param == 0){
                cout << "\tAdapted Mu_" << trans_str.str() << " " << theta_str.str() << endl;
            }
            else {
                cout << "Adapted Lambda_" << param_str.str() << "_" + trans_str.str() << " " << theta_str.str() << endl;
            }
            n++;
        }
        
    }


    
    cout << "//////////////////////////////////////////////////////////////////" << endl
    << "//////////////////////////////////////////////////////////////////" << endl;
    
    // Close dynamic history file
    f_out.close();
    
    cout << "//////////////////////////////////////////////////////////////////" << endl
    << "//////////////////////////////////////////////////////////////////" << endl
    << "                      Post-processing" << endl;
    
    
    // Post-processing: getting the posterior means after different numbers of iterations and computing corresponding likelihoods
    
    bool stop_iter=false;
    unsigned int current_iter=0;
    unsigned int max_BI=100000;
    
    while(!stop_iter){
        // Writing in log file
        cout << "\titer" << current_iter << endl;
        
        
        // Starting at current_BI iterations
        unsigned int current_BI=1000*current_iter;
        
        
        if(current_BI>=n_iter || current_BI>max_BI){
            stop_iter=true;
        } else {
            // Initialisation of the posterior means

            // Storing parameters in the vector theta
            gsl_vector* post_mean_theta=gsl_vector_alloc(n_params);
            for(unsigned int i=0;i<n_params;i++){
                post_mean_theta->data[i]=0;
            }

            // Initialisation of the iterations count
            
            double count=0.0;

            // Computing the posterior means
            
            // Sum of the visited values
            for(unsigned int loop=current_BI;loop<n_iter;loop++){
                count+=1.0;
                for(unsigned int i=0;i<n_params;i++){
                    post_mean_theta->data[i]+=mat_history.matrix[loop*mat_history.nb_columns+2+i];
                }
            }
            
            // Division by number of iterations
            for(unsigned int i=0;i<n_params;i++){
                post_mean_theta->data[i]/=count;
            }

            // Computing the corresponding likelihood
            double post_L=wrapped_likelihood(tensor_indicator,
                    start_indicator,
                    final_indicator,
                    mat_covars,
                    mat_P_I_W,
                    mat_P_I_I,
                    mat_P_I_R,
                    mat_P_W_D,
                    mat_P_W_R,
                    mat_P_W_W,
                    post_mean_theta);
            
            // Writing in the dynamic posterior file
            
            ostringstream current_BI_str;
            ostringstream post_L_str;
            current_BI_str << current_BI;
            post_L_str << post_L;
            string posterior_str = current_BI_str.str() + "\t" + post_L_str.str();
            for(unsigned int i=0; i<n_params; i++){
                ostringstream theta_str;
                theta_str << post_mean_theta->data[i];
                posterior_str += "\t" + theta_str.str();
            }
            f_outL << posterior_str << endl;

        }
        // Incrementing iteration
        current_iter++;
    }
    
    // Writing in the log file
    cout << "//////////////////////////////////////////////////////////////////" << endl
    << "//////////////////////////////////////////////////////////////////" << endl;
    
    
    // Closing dynamic posterior file
    f_outL.close();
    
    
    // Removing all objects
    
    mat_history.Free_double_matrix_cont();
    mat_status.Free_int_matrix_cont();
    mat_covars.Free_double_matrix_cont();

    mat_P_I_W.Free_double_matrix_cont();
    mat_P_I_R.Free_double_matrix_cont();
    mat_P_I_I.Free_double_matrix_cont();
    mat_P_W_D.Free_double_matrix_cont();
    mat_P_W_R.Free_double_matrix_cont();
    mat_P_W_W.Free_double_matrix_cont();
    
    

    cout << "End." << endl;

    return 0;
}

