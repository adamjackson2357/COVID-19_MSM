#include "Int_Matrices_cont.h"
#define DEBUG 0

using namespace std;


Int_Matrices_cont::Int_Matrices_cont()
{
  nb_rows=0;
  nb_columns=0;
}

void Int_Matrices_cont::Alloc_int_matrix_cont(unsigned int Rows,
					      unsigned int Columns)
{
  
  nb_rows=Rows;
  nb_columns=Columns;
  
  matrix=new int[Rows*Columns];
  for(unsigned int i=0;i<Rows*Columns;i++){
    matrix[i]=0;
  }
}

void Int_Matrices_cont::Alloc_int_matrix_cont_val(unsigned int Rows,
						    unsigned int Columns,
                unsigned int value)
{
  
  nb_rows=Rows;
  nb_columns=Columns;
  
  matrix=new int[Rows*Columns];
  for(unsigned int i=0;i<Rows*Columns;i++){
    matrix[i]=value;
  }
}


void Int_Matrices_cont::Reset_int_matrix_cont()
{
  for(unsigned int i=0;i<nb_rows*nb_columns;i++){
    matrix[i]=0;
  }
}


void Int_Matrices_cont::Free_int_matrix_cont()
{
  delete[] matrix;
}


void Int_Matrices_cont::Read_from_file(char *filename)
{
  
  ifstream INFILE;
  INFILE.open(filename, ios::in);
  //Checking the file path
  
  if(INFILE.fail()){
    cout << "Invalid Path and/or permission rights for " << filename << " -- run stopped." << endl;
    exit(1);
  }
  
  //  Reading the first element: nb_rows
  INFILE >> nb_rows;
  //  Reading the second element: nb_columns
  INFILE >> nb_columns;
  cout << "Reading file: " << filename
       << " -- #rows: " << nb_rows
       << " -- #cols: " << nb_columns << endl;
  //Allocating memory for the Matrix
  Alloc_int_matrix_cont(nb_rows,
			nb_columns);
  
  //Storing the map file in the matrix.
  for(unsigned int current_line=0;current_line<nb_rows;current_line++){
    for(unsigned int current_column=0;current_column<nb_columns;current_column++){
      if(!INFILE.eof()){
	if(DEBUG){
	  cout << " ligne " << current_line
	       << " -- colonne " << current_column
	       << " --Test " << INFILE.eof() << endl;
	}
	INFILE>>matrix[current_line*nb_columns+current_column];
      }
      else{
	cout << "Missing element at line " << current_line
	     << " and column " << current_column
	     << " in file " << filename << "." << endl
	     << "!!!!Run Stopped!!!!" << endl;
	exit(1);
      }
    }
  } 
}


void Int_Matrices_cont::Display_matrix()
{
  for(unsigned int current_line=0;current_line<nb_rows;current_line++){
    for(unsigned int current_column=0;current_column<nb_columns;current_column++){
      cout << matrix[current_line*nb_columns+current_column] << " ";
    }
    cout << endl;
  }
  cout << endl;
}


void Int_Matrices_cont::Display_matrix_header()
{
  
  cout << nb_rows << endl;
  cout << nb_columns << endl;
  
  for(unsigned int current_line=0;current_line<nb_rows;current_line++){
    for(unsigned int current_column=0;current_column<nb_columns;current_column++){
      cout << matrix[current_line*nb_columns+current_column] << " ";
    }
    cout << endl;
  }
  cout << endl;
  
}


void Int_Matrices_cont::Write_to_file(char *filename)
{
  
  ofstream OUTFILE;
  OUTFILE.open(filename, ios::out);
  //Checking the file path
  
  if(OUTFILE.fail()){
    cout << "Invalid Path and/or permission rights for " << filename << " -- run stopped." << endl;
    exit(1);
  }
  
  //Writing the first element: nb_rows
  OUTFILE << nb_rows << endl;
  //Writing the second element: nb_columns
  OUTFILE << nb_columns << endl;
  
  //Writing the core of the matrix.
  for(unsigned int current_line=0;current_line<nb_rows;current_line++){
    for(unsigned int current_column=0;current_column<nb_columns;current_column++){
      OUTFILE << matrix[current_line*nb_columns+current_column] << " ";
    }
    OUTFILE << endl;
  } 
  OUTFILE.close();
}
void Int_Matrices_cont::Write_to_file_bis(ofstream *OUTFILE)
{
  
  //Writing the first element: nb_rows
  (*OUTFILE) << nb_rows << endl;
  //Writing the second element: nb_columns
  (*OUTFILE) << nb_columns << endl;
  
  //Writing the core of the matrix.
  for(unsigned int current_line=0;current_line<nb_rows;current_line++){
    for(unsigned int current_column=0;current_column<nb_columns;current_column++){
      (*OUTFILE) << matrix[current_line*nb_columns+current_column] << " ";
    }
    (*OUTFILE) << endl;
  } 
}

