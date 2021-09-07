#ifndef STRING_MATRICES_H
#define STRING_MATRICES_H 1

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <algorithm>
#include <string.h>


using namespace std;


class String_Matrices
{
public:
  String_Matrices();
  ~String_Matrices(){};

  int nb_rows;
  int nb_columns;
  string **matrix;
  

  void Alloc_string_matrix(int Rows,
			   int Columns);
  void Free_string_matrix();
  void Read_from_file(char *filename);
  void Display_matrix();
  void Display_matrix_header();
  void Write_to_file(char *filename);

};

#endif /* !defined STRING_MATRICES_H */
