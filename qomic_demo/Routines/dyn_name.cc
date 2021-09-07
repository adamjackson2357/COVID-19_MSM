#include "dyn_name.h"
#include "struc.h"
#include <string>
#include <iostream>
#include <sstream> 



//using namespace std;


string write_scenario_fname(char * fname,
			     string name,
			     string extension)
{
  string separator="/";
  string temp_res1=fname;

  string Result=temp_res1+separator+name+extension;
  return Result;  
}


string write_param_cols(string param_cols,
                        unsigned int n_trans,
                        unsigned int n_params)
{
  for(unsigned int trans=0; trans<n_trans; trans++){
    for(unsigned int param=0; param<n_params; param++){
        ostringstream trans_str;
        ostringstream param_str;
        trans_str << trans;
        param_str << param;

        if (param == 0){
            param_cols += "\tMu_" + trans_str.str();
        }
        else {
            param_cols += "\tLambda_" + param_str.str() + "_" + trans_str.str();
        }
    }
  }
  return param_cols;
}


string write_sigma_cols(string sigma_cols,
                        unsigned int n_trans,
                        unsigned int n_params)
{
  for(unsigned int trans_i=0; trans_i<n_trans; trans_i++){
    for(unsigned int param_i=0; param_i<n_params; param_i++){
      for(unsigned int trans_j=0; trans_j<n_trans; trans_j++){
        for(unsigned int param_j=0; param_j<n_params; param_j++){
          ostringstream trans_i_str, param_i_str, trans_j_str, param_j_str;
          trans_i_str << trans_i;
          param_i_str << param_i;
          trans_j_str << trans_j;
          param_j_str << param_j;

          string sigma_i, sigma_j;
          sigma_i = param_i_str.str() + "_" + trans_i_str.str();
          sigma_j = param_j_str.str() + "_" + trans_j_str.str();

          if(sigma_i == sigma_j) {
              if(param_i==0 & param_j==0){
                  sigma_cols += "\tsigma_mu_" + trans_i_str.str();
              }
              else{
                  sigma_cols += "\tsigma_lambda_" + param_i_str.str() + "_" + trans_i_str.str();
              }
          }
          else{
              sigma_cols += "\tsigma_" + sigma_i + "_" + sigma_j;
          }
        }
      }
    }
  }
  return sigma_cols;
}


string write_dynamic_fname(char * fname,
			                string name)
{
  string temp_res1=fname;
  string Result=temp_res1+name;
  return Result;  
}
