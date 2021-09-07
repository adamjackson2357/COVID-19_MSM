#include "get_input.h"

#define DEBUG 0


unsigned int get_bi_trans(Int_Matrices_cont mat_trans)
{
  unsigned int n_trans=0;
  unsigned int n_from=mat_trans.nb_rows;
  unsigned int n_to=mat_trans.nb_columns;

  for(unsigned int i=0;i<n_from;i++){
    for(unsigned int j=0;j<n_to;j++){
      unsigned int trans=mat_trans.matrix[i*mat_trans.nb_columns+j];
      if((trans==1) & (i!=j)){
        n_trans++;
      }  
    }
  }
  return n_trans;
}


unsigned int get_tot_trans(Int_Matrices_cont mat_trans)
{
  unsigned int tot_trans=0;
  unsigned int n_from=mat_trans.nb_rows;
  unsigned int n_to=mat_trans.nb_columns;

  for(unsigned int i=0;i<n_from;i++){
    for(unsigned int j=0;j<n_to;j++){
      unsigned int trans=mat_trans.matrix[i*mat_trans.nb_columns+j];
      if(trans==1){
        tot_trans++;
      }
    }
  }
  return tot_trans;
}


unsigned int get_n_params(gsl_vector* theta,
                          unsigned int n_trans)
{
  int N = theta->size;
  unsigned int n_params = N / n_trans;
  return n_params;
}


void get_col_matrix(Int_Matrices_cont mat_col,
                    Int_Matrices_cont mat_trans,
                    unsigned int n_from,
                    unsigned int n_to)
{
  unsigned int col = 0;
  for(unsigned int i=0;i<n_from;i++){
    for(unsigned int j=0;j<n_to;j++){
      unsigned int trans=mat_trans.matrix[i*mat_trans.nb_columns+j];
      if(trans==1){
        mat_col.matrix[i*mat_col.nb_columns+j]=col;
        col++;
      }
    }
  }
}


void change_values(Int_Matrices_cont mat,
                  int value)
{
  unsigned int n_rows = mat.nb_rows;
  unsigned int n_cols = mat.nb_columns;
  for(unsigned int i=0;i<n_rows;i++){
    for(unsigned int j=0;j<n_cols;j++){
      mat.matrix[i*mat.nb_columns+j]=value;
    }
  }
}

void get_from_matrix(Int_Matrices_cont mat_from,
                    Int_Matrices_cont mat_status,
                    Double_Matrices_cont mat_dates,
                    unsigned int n_ind,
                    unsigned int n_to)
{
  
  for(unsigned int ind=0;ind<n_ind;ind++){
    unsigned int col = 0;
    for(unsigned int j=0;j<n_to;j++){
      unsigned int trans=mat_status.matrix[ind*mat_status.nb_columns+j];
      if(trans==1){
        mat_from.matrix[ind*mat_from.nb_columns+col]=mat_dates.matrix[ind*mat_dates.nb_columns+j];
        col++;
      }
    }
  }
}


void get_to_matrix(Int_Matrices_cont mat_to,
                    Int_Matrices_cont mat_status,
                    Double_Matrices_cont mat_dates,
                    unsigned int n_ind,
                    unsigned int n_to)
{
  
  for(unsigned int ind=0;ind<n_ind;ind++){
    unsigned int col = 0;
    for(unsigned int j=1;j<n_to;j++){
      unsigned int trans=mat_status.matrix[ind*mat_status.nb_columns+j];
      if(trans==1){
        mat_to.matrix[ind*mat_to.nb_columns+col]=mat_dates.matrix[ind*mat_dates.nb_columns+j];
        col++;
      }
    }
  }
}


void get_diff_matrix(Int_Matrices_cont mat_diff,
                    Int_Matrices_cont mat_from,
                    Int_Matrices_cont mat_to,
                    unsigned int n_ind,
                    unsigned int n_to)
{
  
  for(unsigned int ind=0;ind<n_ind;ind++){
    for(unsigned int j=0;j<n_to;j++){
      unsigned int diff=mat_to.matrix[ind*mat_to.nb_columns+j]-mat_from.matrix[ind*mat_from.nb_columns+j];
      mat_diff.matrix[ind*mat_diff.nb_columns+j]=diff;
    }
  }
}


void get_pos_from_matrix(Int_Matrices_cont mat_pos_from,
                    Int_Matrices_cont mat_status,
                    unsigned int n_ind,
                    unsigned int n_to)
{
  
  for(unsigned int ind=0;ind<n_ind;ind++){
    unsigned int col = 0;
    for(unsigned int j=0;j<n_to;j++){
      unsigned int trans=mat_status.matrix[ind*mat_status.nb_columns+j];
      if(trans==1){
        mat_pos_from.matrix[ind*mat_pos_from.nb_columns+col]=j;
        col++;
      }
    }
  }
}


void get_pos_to_matrix(Int_Matrices_cont mat_pos_to,
                    Int_Matrices_cont mat_status,
                    unsigned int n_ind,
                    unsigned int n_to)
{
  
  for(unsigned int ind=0;ind<n_ind;ind++){
    unsigned int col = 0;
    for(unsigned int j=1;j<n_to;j++){
      unsigned int trans=mat_status.matrix[ind*mat_status.nb_columns+j];
      if(trans==1){
        mat_pos_to.matrix[ind*mat_pos_to.nb_columns+col]=j;
        col++;
      }
    }
  }
}


void initialise_indicator(
  std::vector<int> &start_indicator,
  unsigned int n_ind,
  unsigned int value)
{
    for (unsigned int i=0;i<n_ind;i++) {
        start_indicator.push_back(value);
    }
}


void initialise_tensor_indicator(
  std::vector<std::vector<std::vector<int> > > &tensor_indicator,
  unsigned int tot_trans,
  unsigned int n_ind)
{
    for(unsigned int t = 0; t < tot_trans; t++){
        std::vector<std::vector<int> > mat_indicator;

        for (unsigned int i = 0; i < n_ind; i++) {

            std::vector<int> vect_indicator;
            vect_indicator.push_back(-1);
            mat_indicator.push_back(vect_indicator);
            vect_indicator.clear();
        }
        tensor_indicator.push_back(mat_indicator);
        mat_indicator.clear();
    }
}


void create_start_indicator(
    std::vector<std::vector<std::vector<int> > > &tensor_indicator,
    std::vector<int> &start_indicator,
    unsigned int n_ind,
    unsigned int n_days,
    unsigned int n_trans)
{
    for(unsigned int ind=0;ind<n_ind;ind++){
        unsigned int first_day=n_days;
        for (unsigned int t=0;t<n_trans;t++){
            unsigned int trans_days=tensor_indicator[t][ind].size();
            if(trans_days>1){
                unsigned int trans_min=tensor_indicator[t][ind][1];
                first_day=min(first_day, trans_min);
            }
        }
        start_indicator.push_back(first_day);
    }
}


void create_final_indicator(
    std::vector<std::vector<std::vector<int> > > &tensor_indicator,
    std::vector<int> &final_indicator,
    unsigned int n_ind,
    unsigned int n_trans)
{
    for(unsigned int ind=0;ind<n_ind;ind++){
        unsigned int final_day=0;
        for (unsigned int t=0;t<n_trans;t++){
            unsigned int trans_days=tensor_indicator[t][ind].size();
            if(trans_days>1){
                unsigned int trans_max=tensor_indicator[t][ind][trans_days-1];
                final_day=max(final_day, trans_max);
            }
        }
        final_indicator.push_back(final_day);
    }
}


void create_tensor_indicator(
  std::vector<std::vector<std::vector<int> > > &tensor_indicator,
  Int_Matrices_cont mat_from,
  Int_Matrices_cont mat_to,
  Int_Matrices_cont mat_pos_from,
  Int_Matrices_cont mat_pos_to,
  Int_Matrices_cont mat_col,
  unsigned int n_ind,
  unsigned int n_from)
{
  for(unsigned int ind=0;ind<n_ind;ind++){
    for(unsigned int i=0;i<n_from;i++){

      // cout << "ind "<< ind << endl;

      // Get the start end and difference
      int from = mat_from.matrix[ind*mat_from.nb_columns+i];
      int to = mat_to.matrix[ind*mat_to.nb_columns+i];
      int diff = (to-from)-1;

      // Get the positions of the from and to
      int from_pos = mat_pos_from.matrix[ind*mat_pos_from.nb_columns+i];
      int to_pos = mat_pos_to.matrix[ind*mat_pos_to.nb_columns+i];

      // cout << "ind " << ind << " from " << from << " to " << to << " diff "
      // << diff << " from_pos " << from_pos << " to_pos " << to_pos << endl;

      // Get the same transition
      int same = mat_col.matrix[from_pos*mat_col.nb_columns+from_pos];

      // cout << "same " << same << endl;

      // If they have spent time in the state and the state is not absorbing, then add
      if((diff>0 && (from_pos<n_from))){
        // cout << "diff " << diff << endl;
        for(int j=from;j<(to-1);j++){
          // The final A-A transition happened two days before the first day in the new state
          // cout << same << "-" << same << " transition " << " day " << j << endl;
          tensor_indicator[same][ind].push_back(j);
        }
      }

      // If there are no further transitions, break
      if(to_pos==(-1)){break;}

      // Otherwise add the change
      else{
        int change = mat_col.matrix[from_pos*mat_col.nb_columns+to_pos];

        // cout << "change" << endl;
        // cout << same << "-" << change << " transition" << " day " << (to-1) << endl;

        // The A-B transition happened the day before the first day in the new state
        tensor_indicator[change][ind].push_back(to-1);
      }
    }
  }
}


void summary_vector_indicator(
  std::vector<int> &vector_indicator,
  unsigned int n_ind)
{
    for (unsigned int i=0;i<n_ind;i++) {
      cout << " Individual " << i << "  Day " << vector_indicator[i] << endl;
    }
}


void summary_tensor_indicator(
  std::vector<std::vector<std::vector<int> > > &tensor_indicator,
  unsigned int n_ind,
  unsigned int tot_trans)
{
    for (unsigned int i=0;i<n_ind;i++){
        for (unsigned int t=0;t<tot_trans;t++){
          cout << " Individual " << i << " Transition " << t << endl;
            for (unsigned int j=0;j<tensor_indicator[t][i].size();j++){
                 cout << "Day " << tensor_indicator[t][i][j] << " ";
            }   
            cout << endl;
        }
    }
}
