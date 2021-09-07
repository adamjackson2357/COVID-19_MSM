#include "manip_dates.h"

#define DEBUG 0

using namespace std;

unsigned int get_n_days(gsl_vector *vect_doc,
			  gsl_vector *vect_dor,
        gsl_vector *vect_dod)
{
  double min_doc=gsl_vector_min(vect_doc);
  double max_doc=gsl_vector_max(vect_doc);
  double min_dor=gsl_vector_min(vect_dor);
  double max_dor=gsl_vector_max(vect_dor);
  double min_dod=gsl_vector_min(vect_dod);
  double max_dod=gsl_vector_max(vect_dod);

  double max1=max(max_doc,max_dor);
  double max2=max(max1,max_dod);
  unsigned int n_days=max2;

  cout << "Min DOC " << min_doc
       << " -- Max DOC " << max_doc
       << " -- Min DOR " << min_dor
       << " -- Max DOR " << max_dor
       << " -- Min DOD " << min_dod
       << " -- Max DOD " << max_dod
       << " -- n_days " << n_days
        << endl;

  return n_days;
  
}
