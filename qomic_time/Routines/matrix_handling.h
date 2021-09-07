#ifndef MATRIX_HANDLING_H
#define MATRIX_HANDLING_H

#include <iostream>
#include <sstream> 
#include <fstream> 
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <stdlib.h>
#include <stddef.h>
#include <stdarg.h>
#include <string.h>
#include <stdio.h>
#include "../Classes/String_Matrices.h"
#include "../Classes/Double_Matrices_cont.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <vector>

using namespace std; 

void standardize_matrix_cont(Double_Matrices_cont M);
gsl_matrix* Double_matrices_cont_2_gsl_matrix(Double_Matrices_cont source);

unsigned int sum_vector_int(vector <unsigned int> &Myvector);
void get_list_var_in_and_out(vector <unsigned int> &list_columns_var_in,
			     vector <unsigned int> &list_columns_var_out,
			     vector <unsigned int> &is_var_in);

gsl_matrix* get_sub_X(vector <unsigned int> &list_columns_var_in,
		      gsl_matrix *mat_X);
gsl_vector* get_one_vect(unsigned int pos_col,
			 gsl_matrix *mat_X);

gsl_matrix* get_X_reduced(vector <unsigned int> &list_columns_var,
			    gsl_matrix *mat_X);
gsl_matrix* get_X_reduced_and_constant(vector <unsigned int> &list_columns_var,
				       gsl_matrix *mat_X);

gsl_matrix* get_sub_matrix_col(gsl_matrix *mat_X,
			       size_t first_col,
			       size_t last_col);

gsl_matrix* get_sub_matrix_row(gsl_matrix *mat_X,
			       size_t first_row,
			       size_t last_row);
void display_gsl_matrix(gsl_matrix *M);

void display_gsl_vector(gsl_vector *V);
void display_gsl_perm(gsl_permutation *P);

void fill_sub_matrix_row(gsl_matrix *X_red,
			 gsl_matrix *mat_X,
			 size_t first_row,
			 size_t last_row);
void fill_sub_matrix_col(gsl_matrix *X_red,
			 gsl_matrix *mat_X,
			 size_t first_col,
			 size_t last_col);


void display_matrix_var_dim(vector < vector <unsigned int> > &M);
void center_matrix_gsl(gsl_matrix *M);
void display_vector_int(vector < unsigned int> &vector);
void display_vector_double(vector < double > &vector);
unsigned int sum_line_std_mat(vector < vector <unsigned int> > &M,
			      unsigned int line);
void get_list_var_in(vector <unsigned int> &list_columns_var_in,
		     vector <unsigned int> &is_var_in);


#endif
