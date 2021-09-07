#ifndef DOUBLE_MATRICES_CONT_H
#define DOUBLE_MATRICES_CONT_H 1

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <algorithm>
#include <string.h>


using namespace std;


class Double_Matrices_cont
{
public:
  Double_Matrices_cont();
  ~Double_Matrices_cont(){};

  unsigned int nb_rows;
  unsigned int nb_columns;
  double *matrix;
  

  void Alloc_double_matrix_cont(unsigned int Rows, unsigned int Columns);
  void Alloc_double_matrix_cont_val(unsigned int Rows, unsigned int Columns, double value);
  void Replicate_double_matrix_cont(Double_Matrices_cont Source);
  void Reset_double_matrix_cont();
  void Free_double_matrix_cont();
  void Read_from_file(char *filename);
  void Display_matrix();
  void Display_matrix_header();
  void Write_to_file(char *filename);
  void Write_to_file_bis(ofstream *OUTFILE);


};

#endif /* !defined DOUBLE_MATRICES_CONTH */
