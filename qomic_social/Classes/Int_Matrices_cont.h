#ifndef INT_MATRICES_CONT_H
#define INT_MATRICES_CONT_H 1

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <algorithm>
#include <string.h>


using namespace std;


class Int_Matrices_cont
{
public:
  Int_Matrices_cont();
  ~Int_Matrices_cont(){};

  unsigned int nb_rows;
  unsigned int nb_columns;
  int *matrix;
  

  void Alloc_int_matrix_cont(unsigned int Rows, unsigned int Columns);
  void Alloc_int_matrix_cont_val(unsigned int Rows, unsigned int Columns, unsigned int value);
  void Reset_int_matrix_cont();
  void Free_int_matrix_cont();
  void Read_from_file(char *filename);
  void Display_matrix();
  void Display_matrix_header();
  void Write_to_file(char *filename);
  void Write_to_file_bis(ofstream *OUTFILE);


};

#endif /* !defined INT_MATRICES_CONT_H */
