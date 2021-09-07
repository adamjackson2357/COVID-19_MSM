#ifndef STRUC_H
#define STRUC_H


struct data
{
  int nb_rows;
  int nb_columns;
  float **matrix;
};
typedef struct data data;


struct data_double
{
  int nb_rows;
  int nb_columns;
  double **matrix;
};
typedef struct data_double data_double;


struct data_3D
{
  int nb_rows;
  int nb_columns;
  int vect_size;
  float ***matrix;
};
typedef struct data_3D data_3D;


struct data_double_3D
{
  int nb_rows;
  int nb_columns;
  int vect_size;
  double ***matrix;
};
typedef struct data_double_3D data_double_3D;

struct data_integer_3D
{
  int nb_rows;
  int nb_columns;
  int vect_size;
  int ***matrix;
};

typedef struct data_integer_3D data_integer_3D;

struct data_integer
{
  int nb_rows;
  int nb_columns;
  int **matrix;
};
typedef struct data_integer data_integer;

struct vect_double
{
  double nb_rows;
  double *matrix;
};
typedef struct vect_double vect_double;

struct data_long
{
  int nb_rows;
  int nb_columns;
  long **matrix;
};
typedef struct data_long data_long;

struct data_long_double
{
  int nb_rows;
  int nb_columns;
  long double **matrix;
};
typedef struct data_long_double data_long_double;
struct data_long_double_3D
{
  int nb_rows;
  int nb_columns;
  int vect_size;
  long double ***matrix;
};
typedef struct data_long_double_3D data_long_double_3D;
#endif /* STRUC_H */
 

