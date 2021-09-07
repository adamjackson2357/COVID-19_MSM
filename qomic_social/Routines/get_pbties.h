#ifndef   	GET_PBTIES_H_
# define   	GET_PBTIES_H_
#include <stdio.h>
#include "../Classes/String_Matrices.h"
#include "../Classes/Double_Matrices_cont.h"
#include "../Classes/Int_Matrices_cont.h"
#include <../Routines/matrix_handling.h>
#include "./rand.h"
#include "./manip_dates.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>

#include <vector>
using namespace std;

double get_prob(double arg_tot);

double get_intercept_prob(double mu);

double get_I_W(Double_Matrices_cont mat_covars,
                            unsigned int ind,
                            double mu,
                            double lambda_1,
                            double lambda_2,
                            double lambda_3,
                            double lambda_4,
                            double lambda_5,
                            double lambda_6,
                            double lambda_7);

void get_from_I(Double_Matrices_cont mat_P_I_W,
                Double_Matrices_cont mat_P_I_I,
                Double_Matrices_cont mat_P_I_R,
                Double_Matrices_cont mat_covars,
                gsl_vector* theta_temp,
                std::vector<int> start_indicator,
                std::vector<int> final_indicator);

void get_from_W(Double_Matrices_cont mat_P_W_D,
                Double_Matrices_cont mat_P_W_R,
                Double_Matrices_cont mat_P_W_W,
                Double_Matrices_cont mat_covars,
                gsl_vector* theta_temp,
                std::vector<int> start_indicator,
                std::vector<int> final_indicator);

double calculate_likelihood(double likelihood,
                            Double_Matrices_cont mat_P,
                            std::vector<std::vector<std::vector<int> > > tensor_indicator,
                            unsigned int trans);

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
                          gsl_vector* theta);


unsigned int candidate_acceptance(double current_L,
                                  double previous_L);


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
                                    unsigned int n_params);

double adaptive_step_multiv(Double_Matrices_cont mat_history,
                            gsl_matrix* hat_sigma,
                            unsigned int iter,
                            unsigned int n_iter_adaptive,
                            unsigned int N);

#endif 	    /* !GET_PBTIES_H_ */
