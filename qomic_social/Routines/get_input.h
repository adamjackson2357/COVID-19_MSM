#ifndef   	GET_INPUT_H_
# define   	GET_INPUT_H_
#include <stdio.h>
#include "../Classes/String_Matrices.h"
#include "../Classes/Double_Matrices_cont.h"
#include "../Classes/Int_Matrices_cont.h"
#include "./rand.h"
#include "./manip_dates.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <vector>
using namespace std;

unsigned int get_bi_trans(Int_Matrices_cont mat_trans);

unsigned int get_tot_trans(Int_Matrices_cont mat_trans);

unsigned int get_n_params(gsl_vector* theta,
                          unsigned int n_trans);

void get_col_matrix(Int_Matrices_cont mat_col,
                    Int_Matrices_cont mat_trans,
                    unsigned int n_from,
                    unsigned int n_to);

void change_values(Int_Matrices_cont mat,
                  int value);

void get_from_matrix(Int_Matrices_cont mat_from,
                    Int_Matrices_cont mat_status,
                    Double_Matrices_cont mat_dates,
                    unsigned int n_ind,
                    unsigned int n_to);

void get_to_matrix(Int_Matrices_cont mat_to,
                    Int_Matrices_cont mat_status,
                    Double_Matrices_cont mat_dates,
                    unsigned int n_ind,
                    unsigned int n_to);

void get_diff_matrix(Int_Matrices_cont mat_diff,
                    Int_Matrices_cont mat_from,
                    Int_Matrices_cont mat_to,
                    unsigned int n_ind,
                    unsigned int n_to);

void get_pos_from_matrix(Int_Matrices_cont mat_pos_from,
                    Int_Matrices_cont mat_status,
                    unsigned int n_ind,
                    unsigned int n_to);

void get_pos_to_matrix(Int_Matrices_cont mat_pos_to,
                    Int_Matrices_cont mat_status,
                    unsigned int n_ind,
                    unsigned int n_to);

void initialise_indicator(
  std::vector<int> &start_indicator,
  unsigned int n_ind,
  unsigned int value);

void create_start_indicator(
    std::vector<std::vector<std::vector<int> > > &tensor_indicator,
    std::vector<int> &start_indicator,
    unsigned int n_ind,
    unsigned int n_days,
    unsigned int n_trans);

void create_final_indicator(
    std::vector<std::vector<std::vector<int> > > &tensor_indicator,
    std::vector<int> &final_indicator,
    unsigned int n_ind,
    unsigned int n_trans);
  
void initialise_tensor_indicator(
  vector < vector< vector < int > > > &tensor_indicator,
  unsigned int tot_trans,
  unsigned int n_ind);

void create_tensor_indicator(
    std::vector<std::vector<std::vector<int> > > &tensor_indicator,
    Int_Matrices_cont mat_from,
    Int_Matrices_cont mat_to,
    Int_Matrices_cont mat_pos_from,
    Int_Matrices_cont mat_pos_to,
    Int_Matrices_cont mat_col,
    unsigned int n_ind,
    unsigned int n_from);

void summary_vector_indicator(
  std::vector<int> &vector_indicator,
  unsigned int n_ind);

void summary_tensor_indicator(
  std::vector<std::vector<std::vector<int> > > &tensor_indicator,
  unsigned int n_ind,
  unsigned int tot_trans);

#endif 	    /* !GET_INPUT_H_ */
