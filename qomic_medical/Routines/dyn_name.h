#ifndef DYN_NAME_H
#define DYN_NAME_H


#include <stdlib.h>
#include <stddef.h>
#include <stdarg.h>
#include <string.h>
#include <stdio.h>
#include <cstring>
#include <iostream>
#include <sstream> 
#include <cstdarg>
#include "struc.h"
using namespace std;

string write_scenario_fname(char * fname,
			     string name,
			     string extension);


string write_param_cols(string param_cols,
                        unsigned int n_trans,
                        unsigned int n_params);


string write_sigma_cols(string sigma_cols,
                        unsigned int n_trans,
                        unsigned int n_params);


string write_dynamic_fname(char * fname,
			                string name);

#endif /* DYN_NAME_H */
