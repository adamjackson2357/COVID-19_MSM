#include "rand.h"
#include "struc.h"

extern "C" {
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_sf.h>

}

using namespace std;

static gsl_rng *RandomNumberGenerator = gsl_rng_alloc( gsl_rng_taus );

double myrand()
{
    return( gsl_rng_uniform( RandomNumberGenerator ) );
}

double myrandRange( double Min, double Max )
{
    double Range = Max - Min;
    
    return( Min + Range * gsl_rng_uniform( RandomNumberGenerator ) );
}

void smyrand( long seed )
{
    gsl_rng_set(RandomNumberGenerator, static_cast< uint32_t >( seed ) );
}

double gengam( double bb, double aa )
{
    return( gsl_ran_gamma( RandomNumberGenerator, aa, 1.0 / bb ) );
}

double genbet( double aa, double bb )
{
    return( gsl_ran_beta( RandomNumberGenerator, aa, bb ) );
}

double gennor( double av, double sd )
{
    return( av + gsl_ran_gaussian( RandomNumberGenerator, sd ) );
}
double gengau(double sd )
{
    return(gsl_ran_gaussian( RandomNumberGenerator, sd ) );
}

double gennorLtail( double av, double sd, double a )
{
    return( av - gsl_ran_gaussian_tail( RandomNumberGenerator, -(a-av), sd ) );
}

double gennorRtail( double av, double sd, double a )
{
    return( av + gsl_ran_gaussian_tail( RandomNumberGenerator, a-av, sd ) );
}

int genbinomial( int n, double p )
{
    return( gsl_ran_binomial( RandomNumberGenerator, p, n ) );
}

void genMultinomial(size_t K, unsigned int N, double *pbty, unsigned int *n )
{
    gsl_ran_multinomial(RandomNumberGenerator, K, N, pbty, n );
}



unsigned int genpoi( double mu )
{
    return( gsl_ran_poisson( RandomNumberGenerator, mu ) );
}


uint32_t SampleFromRange(uint32_t MIN, uint32_t MAX)
{	
    return (uint32_t)( myrand()*(MAX-MIN) + MIN );
}


uint32_t SampleFromDiscrete_new( vector<double> &cdf )
{
    unsigned int k = 0;
    double u = myrand();
    //cout << "u = " << u << endl;
    while( u > cdf[k] && k < cdf.size()){
        //cout << "k= " << k << "test " << cdf[k] << endl;
        k++;
    }
    return k;
}


int SampleFromDiscrete_non_cum(vector<double> &pbty)
{
    unsigned int k = 0;
    vector<double> cdf;
    double u = myrand();
    
    cdf.push_back(pbty[0]);
    
//    cout << "Size of the vector: " << pbty.size() << endl;
//    cout << "P[0] " << pbty[0] << " -- cdf [0] " << cdf[0] << endl;
    for(unsigned int i=1;i<pbty.size();i++){
        cdf.push_back(cdf[i-1]+pbty[i]);
//        cout << "i " << i << " -- pbty " << pbty[i] << " -- cdf " << cdf[i] << endl;
    }
    
    
    // Sampling the state (0, 1 or 2) based on the cumulative probabilities in cdf
    
    //cout << "u = " << u << endl;
    while( u > cdf[k] && k < cdf.size()){
        //cout << "k= " << k << "test " << cdf[k] << endl;
        k++;
    }
    return k;
}


double MyMulti(data_integer *Vector,int total)
{
    double res=0.0;
    double temp_log_cum=0.0;
    for(int i=0;i<Vector->nb_rows;i++){
        temp_log_cum+=gsl_sf_lnfact(Vector->matrix[i][0]);
    }
    temp_log_cum-=gsl_sf_lnfact(total);
    
    res=exp(temp_log_cum);
    
    return res;
}


double exponen(double x)
{
    double val;
    
    if(x>301.0){
        val = 1e300;
    }
    else if(x<-301.0){
        val = 0.0;
    }
    else{
        val = exp(x);
    }
    return(val);
}

