#include <stdio.h>
#include "../Classes/String_Matrices.h"
#include "../Classes/Double_Matrices_cont.h"
#include "./rand.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <vector>

using namespace std; 

unsigned int get_n_days(gsl_vector *vect_doc,
			  gsl_vector *vect_dor,
              gsl_vector *vect_dod);
