#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <omp.h>
#include <unistd.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_errno.h>	// GSL_SUCCESS ...
#include <gsl/gsl_odeiv.h>	// ODE solver

#include <gsl/gsl_types.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_min.h>

#include "nrutil.c"

char *strFileNameDate;
#define FERAL_SETS 22			// observational data
#define EXP_SETS 6		        // number of experimental data sets
#define NUM_PARS 70		        // number of parameters to be passed from main to hood
//#define DIM 15						// number of differential equations

const double h = 0.01;		        // time step

//CK Structure for experimental data!!!
/*
struct dataset{  //  building the structure, to be declared later
	int Date;
	int Tree;
	int Covered;
	double total;
	double fungus;
};

typedef struct
{
	double PARS[NUM_PARS];

	double nuV;
	double nuF;
	double nuR;
	double muF;
	double size_C;

	double R_END[FERAL_SETS+1];
	double INITS[FERAL_SETS+1];
	double POPS[4];
	double EV[100];
	double EF[100];
	double DAY_F[DATA_SETS+1];	    // day where fungal infection first happens
	int MAXT[FERAL_SETS+1];		    // number of days in data set (different for each data set)
	int MAXT2[DATA_SETS+1];		    // number of days in EXPERIMENTAL data set (different for each EXPERIMENTAL data set)
	int MAXT3[DATA_SETS+1];		    // number of days in WEATHER data set (different for each WEATHER data set)
	int ***DATA;				    // 3-dimensional array that holds all the data
	int ***EXPDATA;				    // 3-dimensional array that holds all the EXPERIMENTAL data
	//int ***WDATA;				    // 3-dimensional array that holds all the WEATHER data
	double WDATA[DATA_SETS+1][999][8];
	double test_data[1000][36];

	double AcceptedVect[NUM_PARS];
	double LoopVect[NUM_PARS];

	double parm_low[NUM_PARS];
	double parm_high[NUM_PARS];
	double parm_step[NUM_PARS];

	double MLE[NUM_PARS];
	double MLE_host[FERAL_SETS+1];
	double MLE_initR[FERAL_SETS+1];
	//double MLE_gamma[DATA_SETS+1];

	//double gamma_hood[DATA_SETS+1][100];	// 100 needs to be bigger than parm_inc (hood for each gamma value)
	//double best_gamma[DATA_SETS+1][100];
	double best_initS[FERAL_SETS+1];
	double best_initR[FERAL_SETS+1];

	int parm_inc;

	//size_t dim;
	//int calls;
	double sim_results[55][4];		// 1st entry larger than the number of weeks in any data set
	int th_id;
	int pop;
}STRUCTURE;
*/


typedef struct
{
	double PARS[NUM_PARS];

	int MAXT[FERAL_SETS+1];		    // number of days in data set (different for each data set)

	int MAXT2[EXP_SETS+1];		    // number of days in EXPERIMENTAL data set (different for each EXPERIMENTAL data set)
	int MAXT3[EXP_SETS+1];		    // number of days in WEATHER data set (different for each WEATHER data set)
	int ***DATA;				    // 3-dimensional array that holds all the data
	double ***HostFactors;
	double ***FeralPara;
	double ***ExpPara;
	int ***EXPDATA;				    // 3-dimensional array that holds all the EXPERIMENTAL data
	//double ExpFactors[EXP_SETS+1][999][5];
	double ***ExpFactors;



	double test_data[1000][36];

	double AcceptedVect[NUM_PARS];
	double LoopVect[NUM_PARS];

	double parm_low[NUM_PARS];
	double parm_high[NUM_PARS];
	double parm_step[NUM_PARS];

	double MLE[NUM_PARS];
	double MLE_host[FERAL_SETS+1];
	double MLE_initR[FERAL_SETS+1];
	double best_initS[FERAL_SETS+1];
	double best_initR[FERAL_SETS+1];

	int parm_inc;
	double sim_results[55][4];		// 1st entry larger than the number of weeks in any data set
	int th_id;
	int pop;
	int exp;
}STRUCTURE;

/*
#include "minputdata.h"

#include "random_setup2.h"
#include "filenames3.h"

#include "mbounds.h"

#include "prob_dists.h"

#include "molly_odes.h"

#include "mDDEVF.h"

#include "mlikelihoodtemp.h"
*/

//#include "mDDEVF3_single_test.h"
//#include "mlikelihood3.h"

#include "mboundsb.h"
#include "minputdata.h"
#include "mfilenames.h"

#include "molly_odesH1.h"
#include "mDDEVF3_1.h"
#include "mlikelihood3_m.h"

#include "prob_dists.h"
#include "random_setup2.h"


