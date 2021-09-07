/* This file is part of QOMiC.
 *      Copyright (c) Marc Chadeau-Hyam (m.chadeau@imperial.ac.uk)
 *      2013
 *
 * QOMiC is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * QOMiC is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with QOMiC.  If not, see <http://www.gnu.org/licenses/>.
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <memory>
#include <fstream>
#include <../Routines/manip_dates.h>
#include <../Routines/dyn_name.h>
#include <../Routines/matrix_handling.h>
#include <../Routines/rand.h>
#include "../Routines/get_input.h"
#include "../Routines/get_pbties.h"
#include <../Classes/Double_Matrices_cont.h>
#include <../Classes/Int_Matrices_cont.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

#define DEBUG 1

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
    char filename_in_dates[1000];
    char filename_in_status[1000];
    char filename_in_covars[1000];
    char filename_in_trans[1000];
    char filename_in_MCMC[1000];
    char filename_scenario[1000];
    
    //   char path_name_out[1000];
    
    long MY_SEED=0;
    
    // Initialisation of parameters (default values)
    
    unsigned int n_iter=0;
    unsigned int nb_iter=0;
    unsigned int burn_in=0;

    unsigned int n_trans=1;
    unsigned int n_params=0;
    
    na++;
    while(na < argc){
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
        else if ( 0 == strcmp(argv[na],"-seed") ){
            MY_SEED=(long)((atoi(argv[++na])));
            if (na+1==argc) break;
            na++;
        }
        else if ( 0 == strcmp(argv[na],"-iter") ){ // Number of simulations
            n_iter=atoi(argv[++na]);
            if (na+1==argc) break;
            na++;
        }
        else if ( 0 == strcmp(argv[na],"-MCMC") ){
            strcpy(filename_in_MCMC,argv[++na]);
            if (na+1==argc) break;
            na++;
        }
        else if (0 == strcmp(argv[na],"-scenario")){
            strcpy(filename_scenario,argv[++na]);
            if (na+1==argc) break;
            na++;
        }
        else if ( 0 == strcmp(argv[na],"-nb_iter") ){ // Number of iterations in QOMiC output
            nb_iter=atoi(argv[++na]);
            if (na+1==argc) break;
            na++;
        }
        else if ( 0 == strcmp(argv[na],"-burn_in") ){ // Burn-in of QOMiC output
            burn_in=atoi(argv[++na]);
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

    // Get the total number of params
    int q=n_trans*n_trans; // Total covariance of all the parameters
    int ncol=1+1+n_trans+q; // iter+L+params
    cout << "Number of columns in history matrix: "<< ncol << endl;

    // Dynamic writing of output files
    cout << filename_scenario << endl;
    
    // Dynamic writing of output files

    // Summary file (with numbers of transitions)
    string summary_output_fname=write_dynamic_fname(filename_scenario, "transitions.txt");
    ofstream f_out_summary;
    f_out_summary.open(summary_output_fname.c_str(),ios::out);
    if(f_out_summary.fail()) {
        cout << "Invalid Path and/or permission rights for " << summary_output_fname << " -- run stopped." << endl;
        exit(1);
    }
    else {
        cout <<  summary_output_fname << endl;
    }

    // Conditional file (with numbers of transitions)
    string conditional_output_fname=write_dynamic_fname(filename_scenario, "conditional.txt");
    ofstream f_out_conditional;
    f_out_conditional.open(conditional_output_fname.c_str(),ios::out);
    if(f_out_conditional.fail()) {
        cout << "Invalid Path and/or permission rights for " << conditional_output_fname << " -- run stopped." << endl;
        exit(1);
    }
    else {
        cout <<  conditional_output_fname << endl;
    }

    // Probability files
    string I_I_prob_output_fname=write_dynamic_fname(filename_scenario, "P_I_I.txt");
    string I_W_prob_output_fname=write_dynamic_fname(filename_scenario, "P_I_W.txt");
    string I_R_prob_output_fname=write_dynamic_fname(filename_scenario, "P_I_R.txt");
    string W_W_prob_output_fname=write_dynamic_fname(filename_scenario, "P_W_W.txt");
    string W_R_prob_output_fname=write_dynamic_fname(filename_scenario, "P_W_R.txt");
    string W_D_prob_output_fname=write_dynamic_fname(filename_scenario, "P_W_D.txt");

    ofstream f_out_I_I_prob;
    ofstream f_out_I_W_prob;
    ofstream f_out_I_R_prob;
    ofstream f_out_W_W_prob;
    ofstream f_out_W_R_prob;
    ofstream f_out_W_D_prob;
    
    f_out_I_I_prob.open(I_I_prob_output_fname.c_str(),ios::out);
    f_out_I_W_prob.open(I_W_prob_output_fname.c_str(),ios::out);
    f_out_I_R_prob.open(I_R_prob_output_fname.c_str(),ios::out);
    f_out_W_W_prob.open(W_W_prob_output_fname.c_str(),ios::out);
    f_out_W_R_prob.open(W_R_prob_output_fname.c_str(),ios::out);
    f_out_W_D_prob.open(W_D_prob_output_fname.c_str(),ios::out);

    smyrand((long)(MY_SEED));
    
    // Reading dates file
    Double_Matrices_cont mat_dates;
    mat_dates.Read_from_file(filename_in_dates);
    gsl_matrix *mat_dates_temp=Double_matrices_cont_2_gsl_matrix(mat_dates);
    unsigned int n_days = gsl_matrix_max(mat_dates_temp);
    gsl_matrix_free(mat_dates_temp);

    // Reading status file
    Int_Matrices_cont mat_status;
    mat_status.Read_from_file(filename_in_status);

    // Get the number of individuals
    unsigned int n_ind=mat_status.nb_rows;
    mat_status.Free_int_matrix_cont();

    cout << "n_days_follow_up= " << n_days << " -- n_ind " << n_ind << endl;
    
    // Reading covariates file
    Double_Matrices_cont mat_covars;
    mat_covars.Read_from_file(filename_in_covars);

    // Read the MCMC history file
    Double_Matrices_cont mat_history;
    mat_history.Read_from_file(filename_in_MCMC);

    // Select only relevant columns
    Double_Matrices_cont mat_MCMC;
    mat_MCMC.Alloc_double_matrix_cont(nb_iter, n_params);
    for(unsigned int i=0;i<nb_iter;i++){
        for(unsigned int j=0;j<n_params;j++){
            mat_MCMC.matrix[i*mat_MCMC.nb_columns+j]=mat_history.matrix[i*mat_history.nb_columns+(j+2)];
        }
    }
    mat_history.Free_double_matrix_cont();
    // cout << "MCMC matrix" << endl;
    // mat_MCMC.Display_matrix();

    // Calculating state matrices from input

    // Initialisation of the cumulative indicator matrix
    Double_Matrices_cont mat_indicator_cum;
    mat_indicator_cum.Alloc_double_matrix_cont(n_ind, n_trans);

    // Initialisation of the cumulative conditional matrix
    Double_Matrices_cont mat_conditional_cum;
    mat_conditional_cum.Alloc_double_matrix_cont(n_ind, n_trans);

    // Initialisation of cumulative transition probabilities
    Double_Matrices_cont mat_cum_P_I_I;
    Double_Matrices_cont mat_cum_P_I_W;
    Double_Matrices_cont mat_cum_P_I_R;
    Double_Matrices_cont mat_cum_P_W_W;
    Double_Matrices_cont mat_cum_P_W_R;
    Double_Matrices_cont mat_cum_P_W_D;

    mat_cum_P_I_I.Alloc_double_matrix_cont(n_ind, n_days);
    mat_cum_P_I_W.Alloc_double_matrix_cont(n_ind, n_days);
    mat_cum_P_I_R.Alloc_double_matrix_cont(n_ind, n_days);
    mat_cum_P_W_W.Alloc_double_matrix_cont(n_ind, n_days);
    mat_cum_P_W_R.Alloc_double_matrix_cont(n_ind, n_days);
    mat_cum_P_W_D.Alloc_double_matrix_cont(n_ind, n_days);

    // Create all the probability matrices to the end
    std::vector<int> prob_start_indicator;
    initialise_indicator(prob_start_indicator, n_ind, 0);
    // summary_vector_indicator(final_indicator, n_ind);

    // Create all the probability matrices from zero
    std::vector<int> simul_start_indicator;
    for(unsigned int i=0;i<n_ind;i++){
        unsigned int start=n_days;
        for(unsigned int j=0;j<mat_dates.nb_columns;j++){
            unsigned int day=mat_dates.matrix[i*mat_dates.nb_columns+j];
            start=min(start, day);
        }
        simul_start_indicator.push_back(start);
    }
    mat_dates.Free_double_matrix_cont();
    // cout << "Simulation start indicator" << endl;
    // summary_vector_indicator(simul_start_indicator, n_ind);

    // Create all the probability matrices to the end
    std::vector<int> final_indicator;
    initialise_indicator(final_indicator, n_ind, (n_days-1));
    // summary_vector_indicator(final_indicator, n_ind);

    for(unsigned int iter=0;iter<n_iter;iter++){

        // Sample a random simulation ID from estimation file (after burn-in)
        unsigned int current_sim=(unsigned int)(myrandRange(0,nb_iter-1));

        if(DEBUG==1){
            cout << "Iter: " << iter << " Current simulation: " << current_sim << endl;
        }       
        
        // Extract the parameters of the sampled simulation
        gsl_vector* theta=gsl_vector_alloc(n_params);
        for(unsigned int i=0;i<n_params;i++){theta->data[i]=mat_MCMC.matrix[current_sim*mat_MCMC.nb_columns+i];}
        // cout << "Theta" << endl;
        // display_gsl_vector(theta);


        // Initialisation of transition probabilities
        Double_Matrices_cont mat_P_I_I;
        Double_Matrices_cont mat_P_I_W;
        Double_Matrices_cont mat_P_I_R;
        Double_Matrices_cont mat_P_W_W;
        Double_Matrices_cont mat_P_W_R;
        Double_Matrices_cont mat_P_W_D;

        mat_P_I_I.Alloc_double_matrix_cont(n_ind, n_days);
        mat_P_I_W.Alloc_double_matrix_cont(n_ind, n_days);
        mat_P_I_R.Alloc_double_matrix_cont(n_ind, n_days);
        mat_P_W_W.Alloc_double_matrix_cont(n_ind, n_days);
        mat_P_W_R.Alloc_double_matrix_cont(n_ind, n_days);
        mat_P_W_D.Alloc_double_matrix_cont(n_ind, n_days);
    
        // Get the from I transitions
        int first_param_i=0;
        int last_param_i=24;
        gsl_vector* theta_i=gsl_vector_alloc(last_param_i-first_param_i);
        int i=0;
        for(int j=first_param_i;j<last_param_i;j++){theta_i->data[i]=theta->data[j];i++;}
        get_from_I(mat_P_I_W,
                    mat_P_I_I,
                    mat_P_I_R,
                    mat_covars,
                    theta_i,
                    prob_start_indicator,
                    final_indicator);

        // cout << "P_I_I matrix" << endl;
        // mat_P_I_I.Display_matrix();
        // cout << "P_I_W matrix" << endl;
        // mat_P_I_I.Display_matrix();
        // cout << "P_I_R matrix" << endl;
        // mat_P_I_R.Display_matrix();

        // Get the from W transitions
        int first_param_w=24;
        int last_param_w=48;
        gsl_vector* theta_w=gsl_vector_alloc(last_param_w-first_param_w);
        int w=0;
        for(int j=first_param_w;j<last_param_w;j++){theta_w->data[w]=theta->data[j];w++;}
        get_from_W(mat_P_W_D,
                mat_P_W_R,
                mat_P_W_W,
                mat_covars,
                theta_w,
                prob_start_indicator,
                final_indicator);

        // cout << "P_W_W matrix" << endl;
        // mat_P_W_W.Display_matrix();
        // cout << "P_W_R matrix" << endl;
        // mat_P_W_R.Display_matrix();
        // cout << "P_W_D matrix" << endl;
        // mat_P_W_D.Display_matrix();

        // Storing the probabilities for all individuals

        for(unsigned int i=0;i<n_ind;i++){
            for(unsigned int j=0;j<n_days;j++){
                mat_cum_P_I_I.matrix[i*mat_cum_P_I_I.nb_columns+j]+=mat_P_I_I.matrix[i*mat_P_I_I.nb_columns+j];
                mat_cum_P_I_W.matrix[i*mat_cum_P_I_W.nb_columns+j]+=mat_P_I_W.matrix[i*mat_P_I_W.nb_columns+j];
                mat_cum_P_I_R.matrix[i*mat_cum_P_I_R.nb_columns+j]+=mat_P_I_R.matrix[i*mat_P_I_R.nb_columns+j];
                mat_cum_P_W_W.matrix[i*mat_cum_P_W_W.nb_columns+j]+=mat_P_W_W.matrix[i*mat_P_W_W.nb_columns+j];
                mat_cum_P_W_R.matrix[i*mat_cum_P_W_R.nb_columns+j]+=mat_P_W_R.matrix[i*mat_P_W_R.nb_columns+j];
                mat_cum_P_W_D.matrix[i*mat_cum_P_W_D.nb_columns+j]+=mat_P_W_D.matrix[i*mat_P_W_D.nb_columns+j];
            }
        }

        // Simulate transitions for each individual and time point
        for(unsigned int ind=0;ind<n_ind;ind++){

            unsigned int day=simul_start_indicator[ind];
            int sampled_trans=0;
            unsigned int Stop_simul1=0;
            unsigned int Stop_simul2=0;
            vector < double > vect_I_trans;
            vector < double > vect_W_trans;

            unsigned int I_count=0;
            unsigned int W_count=0;

            // cout << "Individual " << ind << endl;

            while(Stop_simul1==0 && day<n_days){

                // Extract the current probabilities of transition
                
                double current_P_I_I=mat_P_I_I.matrix[ind*mat_P_I_I.nb_columns+day];
                double current_P_I_W=mat_P_I_W.matrix[ind*mat_P_I_W.nb_columns+day];
                double current_P_I_R=mat_P_I_R.matrix[ind*mat_P_I_R.nb_columns+day];

                // cout << "Day " << day <<
                // " Probability of transition I-I: " << current_P_I_I <<
                // " I-W: " << current_P_I_W <<
                // " I-R: " << current_P_I_R << endl;

                // Storing all probabilities in vect_I_trans
                vect_I_trans.push_back(current_P_I_I);
                vect_I_trans.push_back(current_P_I_W);
                vect_I_trans.push_back(current_P_I_R);

                // Sampling the state based on transition probabilities
                sampled_trans=SampleFromDiscrete_non_cum(vect_I_trans);
                mat_indicator_cum.matrix[ind*mat_indicator_cum.nb_columns+sampled_trans]++;

                day++;
                vect_I_trans.clear();
                // cout << "Transition " << sampled_trans << endl;
                I_count++;

                if(sampled_trans==0){
                }
                else{
                    mat_conditional_cum.matrix[ind*mat_conditional_cum.nb_columns+sampled_trans]+=I_count;
                    Stop_simul1=1;
                }
            }

            if(sampled_trans==1){
                // cout << "Moved into W" << endl;
                while(Stop_simul2==0 && day<n_days){
                    double current_P_W_W=mat_P_W_W.matrix[ind*mat_P_W_W.nb_columns+day];
                    double current_P_W_R=mat_P_W_R.matrix[ind*mat_P_W_R.nb_columns+day];
                    double current_P_W_D=mat_P_W_D.matrix[ind*mat_P_W_D.nb_columns+day];

                    // cout << "Day " << day <<
                    // " Probability of transition W-W: " << current_P_W_W <<
                    // " W-C: " << current_P_W_C <<
                    // " W-R: " << current_P_W_R  <<
                    // " W-D: " << current_P_W_D << endl;

                    // Storing all three probabilities in vect_I_trans
                    vect_W_trans.push_back(current_P_W_W);
                    vect_W_trans.push_back(current_P_W_R);
                    vect_W_trans.push_back(current_P_W_D);

                    // Sampling the state based on transition probabilities
                    sampled_trans=SampleFromDiscrete_non_cum(vect_W_trans)+3;
                    mat_indicator_cum.matrix[ind*mat_indicator_cum.nb_columns+sampled_trans]++;

                    // cout << "Transition " << sampled_trans << endl;

                    day++;
                    vect_W_trans.clear();
                    W_count++;

                    if(sampled_trans==3){
                    }
                    else{
                        mat_conditional_cum.matrix[ind*mat_conditional_cum.nb_columns+sampled_trans]+=W_count;
                        Stop_simul2=1;
                    }
                }
            }
        }

        mat_P_I_I.Free_double_matrix_cont();
        mat_P_I_W.Free_double_matrix_cont();
        mat_P_I_R.Free_double_matrix_cont();
        mat_P_W_W.Free_double_matrix_cont();
        mat_P_W_R.Free_double_matrix_cont();
        mat_P_W_D.Free_double_matrix_cont();
    }
    
    // Compute the average number of conditional transitions per individual

    for(unsigned int i=0;i<n_ind;i++){ // For loop over individuals
        for(unsigned int j=0;j<n_trans;j++){ // For loop over transition probabilities
            mat_conditional_cum.matrix[i*mat_conditional_cum.nb_columns+j]/=mat_indicator_cum.matrix[i*mat_indicator_cum.nb_columns+j];
            f_out_conditional << mat_conditional_cum.matrix[i*mat_conditional_cum.nb_columns+j]
            << " ";
        }
        f_out_conditional << endl;
    }

    // Compute the average number of transitions per individual
    
    for(unsigned int i=0;i<n_ind;i++){ // For loop over individuals
        for(unsigned int j=0;j<n_trans;j++){ // For loop over transition probabilities
            mat_indicator_cum.matrix[i*mat_indicator_cum.nb_columns+j]/=n_iter;
            f_out_summary << mat_indicator_cum.matrix[i*mat_indicator_cum.nb_columns+j]
            << " ";
        }
        f_out_summary << endl;
    }

    // Compute the average probability of transition for each individual

    for(unsigned int i=0;i<n_ind;i++){ // For loop over individuals
        for(unsigned int day=0;day<n_days;day++){ // For loop over transition probabilities for each day
            mat_cum_P_I_I.matrix[i*mat_cum_P_I_I.nb_columns+day]/=n_iter;    
            f_out_I_I_prob << mat_cum_P_I_I.matrix[i*mat_cum_P_I_I.nb_columns+day] << " ";
        }
        f_out_I_I_prob << endl;

        for(unsigned int day=0;day<n_days;day++){ // For loop over transition probabilities for each day
            mat_cum_P_I_W.matrix[i*mat_cum_P_I_W.nb_columns+day]/=n_iter;    
            f_out_I_W_prob << mat_cum_P_I_W.matrix[i*mat_cum_P_I_W.nb_columns+day] << " ";
        }
        f_out_I_W_prob << endl;

        for(unsigned int day=0;day<n_days;day++){ // For loop over transition probabilities for each day
            mat_cum_P_I_R.matrix[i*mat_cum_P_I_R.nb_columns+day]/=n_iter;    
            f_out_I_R_prob << mat_cum_P_I_R.matrix[i*mat_cum_P_I_R.nb_columns+day] << " ";
        }
        f_out_I_R_prob << endl;

        for(unsigned int day=0;day<n_days;day++){ // For loop over transition probabilities for each day
            mat_cum_P_W_W.matrix[i*mat_cum_P_W_W.nb_columns+day]/=n_iter;    
            f_out_W_W_prob << mat_cum_P_W_W.matrix[i*mat_cum_P_W_W.nb_columns+day] << " ";
        }
        f_out_W_W_prob << endl;

        for(unsigned int day=0;day<n_days;day++){ // For loop over transition probabilities for each day
            mat_cum_P_W_R.matrix[i*mat_cum_P_W_R.nb_columns+day]/=n_iter;    
            f_out_W_R_prob << mat_cum_P_W_R.matrix[i*mat_cum_P_W_R.nb_columns+day] << " ";
        }
        f_out_W_R_prob << endl;

        for(unsigned int day=0;day<n_days;day++){ // For loop over transition probabilities for each day
            mat_cum_P_W_D.matrix[i*mat_cum_P_W_D.nb_columns+day]/=n_iter;    
            f_out_W_D_prob << mat_cum_P_W_D.matrix[i*mat_cum_P_W_D.nb_columns+day] << " ";
        }
        f_out_W_D_prob << endl;
    }


    // Closing the connections
    f_out_summary.close();
    f_out_I_I_prob.close();
    f_out_I_W_prob.close();
    f_out_I_R_prob.close();
    f_out_W_W_prob.close();
    f_out_W_R_prob.close();
    f_out_W_D_prob.close();
    
    // Clearing all objects
    mat_indicator_cum.Free_double_matrix_cont();
    mat_covars.Free_double_matrix_cont();
    mat_MCMC.Free_double_matrix_cont();
    mat_col.Free_int_matrix_cont();

    mat_cum_P_I_I.Free_double_matrix_cont();
    mat_cum_P_I_W.Free_double_matrix_cont();
    mat_cum_P_I_R.Free_double_matrix_cont();
    mat_cum_P_W_W.Free_double_matrix_cont();
    mat_cum_P_W_R.Free_double_matrix_cont();
    mat_cum_P_W_D.Free_double_matrix_cont();
    
    return 0;
}

