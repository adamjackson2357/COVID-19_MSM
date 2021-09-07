#include "matrix_handling.h"

#define DEBUG 0

using namespace std;


void standardize_matrix_cont(Double_Matrices_cont M)
{
  
  for(unsigned int col=0;col<M.nb_columns;col++){
    double mean=0.0;
    for(unsigned int row=0;row<M.nb_rows;row++){
      mean+=M.matrix[row*M.nb_columns+col];
    }
    mean/=M.nb_rows;
    double temp_var=0.0;
    for(unsigned int row=0;row<M.nb_rows;row++){
      temp_var+=pow((M.matrix[row*M.nb_columns+col]-mean),2.0);
    }
    if(M.nb_rows>1){
      temp_var/=M.nb_rows-1;
    }
    else{
      cout << "USAGE::standardize_matrix :: matrix size <2 to compute variance!!!!! run_stopped " << endl;
      exit(1);
    }
    double Mystd=sqrt(temp_var);
    if(Mystd>0){
      for(unsigned int row=0;row<M.nb_rows;row++){
	M.matrix[row*M.nb_columns+col]/=Mystd;
      }
    }
    else{
      cout << "USAGE::standardize_matrix::  estimated std=0; matrix unchanged. " << endl;
    }
  }
}
void center_matrix_gsl(gsl_matrix *M)
{
  
  for(unsigned int col=0;col<M->size2;col++){
    double mean=0.0;
    for(unsigned int row=0;row<M->size1;row++){
      mean+=M->data[row*M->size2+col];
    }
    mean/=M->size1;
    for(unsigned int row=0;row<M->size1;row++){
      M->data[row*M->size2+col]-=mean;
    }
  }
}

gsl_matrix* Double_matrices_cont_2_gsl_matrix(Double_Matrices_cont source)
{
  gsl_matrix* M=gsl_matrix_calloc(source.nb_rows,source.nb_columns);
  M->data=&source.matrix[0];
  return M;
}



unsigned int sum_vector_int(vector <unsigned int> &Myvector)
{
  unsigned int Mysum=0.0;
  for(unsigned int i=0;i<Myvector.size();i++){
    Mysum+=Myvector[i];
  }
  return Mysum;
}


void get_list_var_in_and_out(vector <unsigned int> &list_columns_var_in,
			     vector <unsigned int> &list_columns_var_out,
			     vector <unsigned int> &is_var_in)
{

  for(unsigned int i=0;i<is_var_in.size();i++){
    if(is_var_in[i]==1){
      list_columns_var_in.push_back(i);
    }
    else{
      list_columns_var_out.push_back(i);
    }
  }
}


void get_list_var_in(vector <unsigned int> &list_columns_var_in,
		     vector <unsigned int> &is_var_in)
{
  
  for(unsigned int i=0;i<is_var_in.size();i++){
    if(is_var_in[i]==1){
      list_columns_var_in.push_back(i);
    }
  }
}



gsl_matrix* get_sub_X(vector <unsigned int> &list_columns_var_in,
		      gsl_matrix *mat_X)
{
  unsigned int n_vars_in=list_columns_var_in.size();
  if (DEBUG){
    cout << "\tIn get_X_gam" << endl
	 << "\t#variables in " << list_columns_var_in.size() << endl;
    for(unsigned int i=0;i<list_columns_var_in.size();i++){
      cout << "\tVariable position (" << i+1 << "/ " << list_columns_var_in.size() << ") : " << list_columns_var_in[i] << endl;
    }
  }
  unsigned int nX=mat_X->size1;
  gsl_matrix *X_gam=gsl_matrix_alloc(nX,n_vars_in);
  for(unsigned int col=0;col<n_vars_in;col++){
    //Copying the the col^th variable (at position list_columns_var_in[col] in mat_X)
    gsl_vector_view current_col = gsl_matrix_column (mat_X,list_columns_var_in[col]);
    //Setting the col^th column of X_gam
    gsl_matrix_set_col (X_gam,col,&current_col.vector);
  }
  return X_gam;
}


gsl_vector* get_one_vect(unsigned int pos_col,
			 gsl_matrix *mat_X)
{
  unsigned int nX=mat_X->size1;
  gsl_vector *selected_col=gsl_vector_alloc(nX);
  gsl_matrix_get_col(selected_col,mat_X,pos_col);
  return selected_col;
}



gsl_matrix* get_X_reduced_and_constant(vector <unsigned int> &list_columns_var,
				       gsl_matrix *mat_X)
{
  unsigned int n_vars=list_columns_var.size();
  if (DEBUG){
    cout << "\tIn get_X_reduced_and_constant" << endl
	 << "\t#variables in " << list_columns_var.size() << endl;
    for(unsigned int i=0;i<list_columns_var.size();i++){
      cout << "\tVariable position (" << i+1 << "/ " << list_columns_var.size() << ") : " << list_columns_var[i] << endl;
    }
  }

  unsigned int nX=mat_X->size1;
  gsl_matrix *X_red=gsl_matrix_alloc(nX,(n_vars)+1);
  
  //The first column is set to 1: the constant term of the regression
  gsl_vector *temp_vector=gsl_vector_alloc(nX);
  gsl_vector_set_all(temp_vector,1.0);
  gsl_matrix_set_col (X_red,0,temp_vector);

  gsl_vector_free(temp_vector);

  //Setting the other columns to the list of vars in
  for(unsigned int col=0;col<n_vars;col++){
    //Copying the the col^th variable (at position list_columns_var_in[col] in mat_X)
    gsl_vector_view current_col = gsl_matrix_column (mat_X,list_columns_var[col]);
    //Setting the col^th column of X_red
    gsl_matrix_set_col (X_red,col+1,&current_col.vector);
  }
  return X_red;
}


gsl_matrix* get_X_reduced(vector <unsigned int> &list_columns_var,
			  gsl_matrix *mat_X)
{
  unsigned int n_vars=list_columns_var.size();
  if (DEBUG){
    cout << "\tIn get_X_reduced" << endl
	 << "\t#variables in " << list_columns_var.size() << endl;
    for(unsigned int i=0;i<list_columns_var.size();i++){
      cout << "\tVariable position (" << i+1 << "/ " << list_columns_var.size() << ") : " << list_columns_var[i] << endl;
    }
  }
  
  unsigned int nX=mat_X->size1;
  gsl_matrix *X_red=gsl_matrix_alloc(nX,(n_vars));
  //Setting the other columns to the list of vars in
  for(unsigned int col=0;col<n_vars;col++){
    //Copying the the col^th variable (at position list_columns_var_in[col] in mat_X)
    gsl_vector_view current_col = gsl_matrix_column (mat_X,list_columns_var[col]);
    //Setting the col^th column of X_red
    gsl_matrix_set_col (X_red,col,&current_col.vector);
  }
  return X_red;
}





gsl_matrix* get_sub_matrix_col(gsl_matrix *mat_X,
			       size_t first_col,
			       size_t last_col)
{
  unsigned int n_vars=last_col-first_col;
  unsigned int nX=mat_X->size1;
  gsl_matrix *X_red=gsl_matrix_alloc(nX,n_vars);
  
  for(unsigned int col=0;col<n_vars;col++){
    //Copying the the col^th variable (at position list_columns_var_in[col] in mat_X)
    gsl_vector_view current_col = gsl_matrix_column (mat_X,first_col+col);
    gsl_matrix_set_col (X_red,col,&current_col.vector);
  }
  return X_red;
}

void fill_sub_matrix_col(gsl_matrix *X_red,
			 gsl_matrix *mat_X,
			 size_t first_col,
			 size_t last_col)
{
  unsigned int n_vars=last_col-first_col;
  
  for(unsigned int col=0;col<n_vars;col++){
    gsl_vector_view current_col = gsl_matrix_column (mat_X,first_col+col);
    gsl_matrix_set_col (X_red,col,&current_col.vector);
  }
}


gsl_matrix* get_sub_matrix_row(gsl_matrix *mat_X,
			       size_t first_row,
			       size_t last_row)
{
  unsigned int n_row=last_row-first_row;
  unsigned int nY=mat_X->size2;
  gsl_matrix *X_red=gsl_matrix_alloc(n_row,nY);
  
  for(unsigned int row=0;row<n_row;row++){
    gsl_vector_view current_row = gsl_matrix_row(mat_X,first_row+row);
    gsl_matrix_set_row(X_red,row,&current_row.vector);
  }
  return X_red;
}



void fill_sub_matrix_row(gsl_matrix *X_red,
			 gsl_matrix *mat_X,
			 size_t first_row,
			 size_t last_row)
{
  unsigned int n_row=last_row-first_row;
  
  for(unsigned int row=0;row<n_row;row++){
    gsl_vector_view current_row = gsl_matrix_row(mat_X,first_row+row);
    gsl_matrix_set_row(X_red,row,&current_row.vector);
  }
}




void display_gsl_matrix(gsl_matrix *M)
{
  cout << "nb_rows " << M->size1 << " -- nb_columns " << M->size2 << endl;
    for(unsigned int i=0;i<M->size1;i++){
      for(unsigned int j=0;j<M->size2;j++){
	cout << M->data[i*(M->size2)+j] << " ";
      }
      cout << endl;
    }
}

void display_gsl_vector(gsl_vector *V)
{
  cout << "vector size " << V->size << endl;
  for(unsigned int i=0;i<V->size;i++){
    cout << V->data[i] << " ";
  }
  cout << endl;
}

void display_gsl_perm(gsl_permutation *P)
{
  cout << "perm size " << P->size << endl;
  for(unsigned int i=0;i<P->size;i++){
    cout << P->data[i] << " ";
  }
  cout << endl;
}

void display_vector_int(vector < unsigned int> &vector)
{
  cout << "Vector size " << vector.size() << endl;
  for(unsigned int i=0;i<vector.size();i++){
    cout << vector[i] << " ";
  }
  cout << endl;
}
void display_vector_double(vector < double> &vector)
{
  cout << "Vector size " << vector.size() << endl;
  for(unsigned int i=0;i<vector.size();i++){
    cout << vector[i] << " ";
  }
  cout << endl;
}


void display_matrix_var_dim(vector < vector <unsigned int> > &M)
{
  for(unsigned int row=0;row<M.size();row++){
    cout << "Row #" << row+1 << " -- nb columns " << M[row].size() << endl;
    for(unsigned int col=0;col<M[row].size();col++){
      cout << M[row][col] << " ";
    }
    if(M[row].size()>0){
      cout << endl;
    }
  }


}

unsigned int sum_line_std_mat(vector < vector <unsigned int> > &M,
			      unsigned int line)
{
  unsigned int tmp_sum=0;

  for(unsigned int col=0;col<M[line].size();col++){
    tmp_sum+=M[line][col];
  }
  return tmp_sum;
}


