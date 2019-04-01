#include "headG1.h"
gsl_rng *r2;
int verbose = 1;
int exrealz = 2;

int crash;
int redo = 0;

int main(int argc, char *argv[])
{
	int test = 0;
	//int test = 66;
	//int test = 99;
int View = 1;//CK// turn to 1 to print line search progress.  turn to 0 to run for real.

STRUCTURE Params;
int pro = 1;//atoi(argv[1]);						// pro and argv[1] are the inputs (argv[i] is the i^th input)

// ------------------------------------- Adustable accuracy vs. speed ------------------------------------------------ //

int num_runs	 = 200;
int calls=2;					// number of stochastic simulations for each parameter and IC set

double parm_inc, host_inc, initR_inc;
parm_inc= 1000.0;		host_inc=9.0;		initR_inc=10.0;
//if (pro==1)	{	parm_inc= 60.0;		host_inc=9.0;		initR_inc=10.0;	} //pro = 0 would be the profile lhood, which we don't do
//else		{	parm_inc=15.0;		host_inc=10.0;	initR_inc=10.0;	}


if (test==66)	{	    num_runs=3;	calls = 3; parm_inc= 10.0; host_inc=19.0;} //MG//  low resolution grid search
if (test==99)	{	    num_runs=5;	calls = 12; parm_inc=20.0; host_inc=19.0;} //MG//  low resolution grid search




// ------------------------------------------------------------------------------------------------------------------ //
int i=0; int j;int ii; int jj; int k;
int run;	            int changer;	    double index, tot_index;

signed int LoopNumber=-1;

int Realizations=3;
double inner_parm;				//double outer_parm;
double inner_parm2;				//double outer_parm;

double nuVholder1;  	//CK// holder for logging nuV on line search

int pop;
double log_pop;
int exp;

// -------------------------------------------- MISER STUFF --------------------------------------------------------- //
inputdata(&Params);				// gets Params.DATA[j][i][0-2] and Params.MAXT[i] from inputdata.h


size_t dim[22]={13, 22, 15, 29, 36, 29, 29, 29, 29, 35, 29, 35, 36, 35, 36, 36, 35, 39, 44, 44, 39, 39};  //MG  length of time observed per site

size_t dim2;  //MG what are you for

// --------------------------------------- Name for Output Files ----------------------------------------------------- //
char strFileName[99];					// from filenames.h

GetString(pro,0,strFileName,98);		fflush(stdout);		//getc(stdin);
FILE *molly_results;

// ---------------------------------------- Random Number Stuff ------------------------------------------------------ //

const gsl_rng_type *T2;
long seed;
seed = time(NULL)*(int)getpid();  //use process id
//seed = -1;
gsl_rng_env_setup ();
gsl_rng_default_seed = seed;

T2 = gsl_rng_default;
r2 = gsl_rng_alloc(T2);
//printf("seed = %d\n", seed);


int num_adj_pars=40;			// number of adjustable parameters


// -------------------- parameter high/low values and increments and fixed parameter values ------------------------- //

//global_fixed_parms(&Params);  // gets Params.PARS[i] for fixed parameters from bounds.h  //MG NO. Get these from random_restart_parms
parm_range_inc(&Params,parm_inc,host_inc,initR_inc,num_adj_pars); // gets Params.parm_set,low,high,R_END from bounds.h
// ------------------------------------ Declare Likelihood Quanitites ----------------------------------------------- //
double pop_lhood, pop_lhood2, pop_err,post_hood;	// population lhood (and posterior lhood) calculated for each initS and initR
double mtotal_lhood;
double pop_best_lhood;					// likelihood and error for best initS and initR
double total_lhood;						// sum of pop_best_lhood over all patches
double best_post_hood;	double best_lhood=0;		// best post_hood and lhood
double prior[num_adj_pars];


double PARS[70];


pop = 22;
Params.PARS[2] = 0.3;		//MG what are you?

double y;
y = (double)gsl_ran_flat(r2,0,1);
//printf("y = %e\n", y);

if(y>0.50){
//printf("random\n");
double z;
for(i = 1; i <= num_adj_pars; i++){
	if(i == 3 || i == 4 || i == 5 || i == 6 || i == 7 || i == 8 || i == 11 || i == 12 || i == 13 || i == 17 || i == 18 || i ==19 || i ==21 || i ==23 || i == 25|| i ==29||i==30||i==31||i==32||i==33||i==35){
		z=(double)gsl_ran_flat(r2,0,1);
		Params.PARS[i] = pow(10,Params.parm_low[i]) + (pow(10,Params.parm_high[i]) - pow(10,Params.parm_low[i]))*z;
		}
	if(i == 9 || i==14 || i ==24 || i==27 || i ==28 || i == 34 || i ==36 ){
		z=(double)gsl_ran_flat(r2,0,1);
		Params.PARS[i] = round(Params.parm_low[i] + (Params.parm_high[i] - Params.parm_low[i])*z);
		}
	//if(i == 14 || i == 24|| i ==34 || i == 36 ){
	//	z=(double)gsl_ran_flat(r2,0,1);
		//Params.PARS[i] = Params.parm_low[i] + (Params.parm_high[i] - Params.parm_low[i])*z;
		//}
	}

	Params.MLE[4] = Params.PARS[4];
	Params.MLE[5] = Params.PARS[5];
	Params.MLE[6] = Params.PARS[6];
	Params.MLE[7] = Params.PARS[7];
	Params.MLE[8] = Params.PARS[8];
	Params.MLE[9] = Params.PARS[9];
	Params.MLE[11] = Params.PARS[11];
	Params.MLE[12] = Params.PARS[12];
	Params.MLE[13] = Params.PARS[13];
	Params.MLE[14] = Params.PARS[14];
	Params.MLE[18] = Params.PARS[18];
	Params.MLE[19] = Params.PARS[19];
	Params.MLE[21] = Params.PARS[21];
	Params.MLE[27] = Params.PARS[27];
	Params.MLE[28] = Params.PARS[28];
	Params.MLE[29] = Params.PARS[29];
	Params.MLE[30] = Params.PARS[30];
	Params.MLE[31] = Params.PARS[31];
	Params.MLE[32] = Params.PARS[32];
	Params.MLE[33] = Params.PARS[33];
	Params.MLE[34] = Params.PARS[34];
	Params.MLE[35] = Params.PARS[35];
	Params.MLE[36] = Params.PARS[36];
	Params.MLE[23] = Params.PARS[23];
	Params.MLE[24] = Params.PARS[24];

}else{
//printf("good\n");

//double CHEAT[25] = {0.673300, 0.007626, 2.680351, 2.0, 0.034365, 2.000000, 0.000000, 0.000000, 5.271877, 7.000000, 0.095217, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000};
///double CHEAT[25] = {4.818167e+00,	3.959668e-02,	3.122251e+00,	9.102355e+00,	8.413903e-01,	7.000000e+00,	1.873870e+00,	4.838239e-01,	1.576408e+00,	4.000000e+00,	2.512002e-02,	1.700965e+00,	1.197538e+00,	5.000000e+00,	1.800000e+01,	5.437387e+00,	3.975619e-01,	1.029817e+00,	2.846302e-02,	3.758114e-02,	7.000000e+00,	2.306498e+01,	6.000000e+00,	1.157221e+01,	6.000000e+00};
//double CHEAT[25] = {4.818,	3.959668e-02,	3.122251e+00,	9.102355e+00,	8.413903e-01,	7.000000e+00,	1.873870e+00,	4.838239e-01,	1.576408e+00,	4.000000e+00,	2.512002e-02,	1.700965e+00,	1.197538e+00,	5.000000e+00,	1.800000e+01,	5.437387e+00,	3.975619e-01,	1.029817e+00,	2.846302e-02,	3.758114e-02,	7.000000e+00,	2.306498e+01,	6.000000e+00,	1.157221e+01,	6.000000e+00};
//double CHEAT[25] = {1.083828e-01,	7.695086e-02,	2.490195e+00,	8.037019e+00,	4.981523e-01,	1.000000e+00,	1.365336e+00,	5.879724e-01,	5.591135e-02,	4.000000e+00,	3.981347e-02,	1.429923e+00,	1.528908e+00,	1.300000e+01,	2.000000e+01,	6.102408e+00,	9.479141e-01,	1.247445e-01,	1.180933e-02,	5.012103e-01,	5.000000e+00,	1.437615e+01,	4.000000e+00,	4.250698e+00,	4.000000e+01};

//double CHEAT[25] = {65.46929,	389.0,	0.559245,	10.76345,	0.0099977,	1.0,	3.316519,	0.00899547,	0.05855177,	9.0,	0.1770363,	6578.014,	9467.699,	16.0,	12.0,	10002.3,	2175.761,	0.00564921,	0.365023,	0.00074776,	9.0,	0.2687516,	14.0,	0.4926673,	13.0};
double CHEAT[25] = {1.584966e+01,	4.693935e-02,	3.643355e-01,	1.585112e+03,	6.308702e-02,	1.000000e+00,
1.083443e+00,	2.453031e-01,	9.331905e-01,	1.100000e+01,	1.585112e-01,	3.665803e+03,	6.434202e+03,	1.000000e+01,	1.600000e+01,
3.453047e+03,	8.071239e+03,	6.743800e-01,	1.797741e-01,	1.000000e+00,	2.000000e+00,
1.479649e-01,	9.000000e+00,	3.076281e-01,	1.700000e+01};


Params.PARS[4] = CHEAT[0];
Params.PARS[5] = CHEAT[1];
Params.PARS[6] = CHEAT[2];
Params.PARS[7] = CHEAT[3];
Params.PARS[8] = CHEAT[4];
Params.PARS[9] = CHEAT[5];
Params.PARS[11] = CHEAT[6];
Params.PARS[12] = CHEAT[7];
Params.PARS[13] = CHEAT[8];
Params.PARS[14] = CHEAT[9];
Params.PARS[18] = CHEAT[10];
Params.PARS[19] = CHEAT[11];
Params.PARS[21] = CHEAT[12];
Params.PARS[27] = CHEAT[13];
Params.PARS[28] = CHEAT[14];
Params.PARS[29] = CHEAT[15];
Params.PARS[30] = CHEAT[16];
Params.PARS[31] = CHEAT[17];
Params.PARS[32] = CHEAT[18];
Params.PARS[33] = CHEAT[19];
Params.PARS[34] = CHEAT[20];
Params.PARS[35] = CHEAT[21];
Params.PARS[36] = CHEAT[22];
Params.PARS[23] = CHEAT[23];
Params.PARS[24] = CHEAT[24];

Params.MLE[4] = CHEAT[0];
Params.MLE[5] = CHEAT[1];
Params.MLE[6] = CHEAT[2];
Params.MLE[7] = CHEAT[3];
Params.MLE[8] = CHEAT[4];
Params.MLE[9] = CHEAT[5];
Params.MLE[11] = CHEAT[6];
Params.MLE[12] = CHEAT[7];
Params.MLE[13] = CHEAT[8];
Params.MLE[14] = CHEAT[9];
Params.MLE[18] = CHEAT[10];
Params.MLE[19] = CHEAT[11];
Params.MLE[21] = CHEAT[12];
Params.MLE[27] = CHEAT[13];
Params.MLE[28] = CHEAT[14];
Params.MLE[29] = CHEAT[15];
Params.MLE[30] = CHEAT[16];
Params.MLE[31] = CHEAT[17];
Params.MLE[32] = CHEAT[18];
Params.MLE[33] = CHEAT[19];
Params.MLE[34] = CHEAT[20];
Params.MLE[35] = CHEAT[21];
Params.MLE[36] = CHEAT[22];
Params.MLE[23] = CHEAT[23];
Params.MLE[24] = CHEAT[24];
}




double VALUE[7] = {87.06402, 4788.836, 0.19767, 1, 0.191644, 1.0, 4.0};


double RandNumsPass[50];
int q;
double error = 0.0;

if(verbose==1){
	for(q = 0; q<1;q++){
		printf("q = %d\n", q);
		FILE *fp;


				///need par value sets here
		fp = fopen("G1_bestparms_output.dat","a");
		fclose(fp);

		for(Params.pop=0;Params.pop<22;Params.pop++){
			for(i=0;i<dim[Params.pop];i++){
				RandNumsPass[i] = gsl_ran_flat(r2,0.0,1.0);
				//printf("RandNumsPass[i] = %e\n", RandNumsPass[i]);
				}
			Hood_Pops(RandNumsPass,dim[Params.pop],Params.PARS);
			}
		}
	exit(1);
	}



int r;
int rlzns;
rlzns = 1;
/*if(verbose==1){
	calls = 2;
	double liketest[rlzns][10];
	for(q = 1; q<=rlzns;q++){
		FILE *fp;

		for(r = 0; r < 5; r++){
			printf("r = %d \n", r);
			double lstorepop;
			double store2pop;
				Params.PARS[4] = VALUE[0];
				//Params.PARS[5] = VALUE[1];
				//Params.PARS[6] = VALUE[1];
				Params.PARS[7] = VALUE[1];
				Params.PARS[8] = VALUE[2];
				Params.PARS[9] = VALUE[3];
				//Params.PARS[11] = VALUE[6];
				//Params.PARS[12] = VALUE[7];
				//Params.PARS[13] = VALUE[4];
				//Params.PARS[14] = VALUE[7];
				Params.PARS[18] = VALUE[4];

				Params.PARS[33] = VALUE[5];
				Params.PARS[34] = VALUE[6];
		for(Params.pop=0;Params.pop<22;Params.pop++){
			pop = Params.pop;
			printf("pop = %d \n", pop);
			for(i=0;i<dim[Params.pop];i++){
				RandNumsPass[i] = gsl_ran_flat(r2,0.0,1.0);
				}//close random loop

			//storepop = Hood_Pops(RandNumsPass,dim[Params.pop],Params.PARS);
				gsl_monte_function G = {&Hood_Pops, dim[pop], &Params };	// declares function calling Hood_Pops.h
				double xl[dim[pop]];	double xu[dim[pop]];	// need to redeclare xl and xu since the size changes
				for (jj=0;jj<=dim[pop];jj++)	{
							xl[jj]=0;
							xu[jj]=1;
							}
			// ----------- Use MISER to call function pop_lhood --------------------------------- //
					gsl_monte_miser_state *s = gsl_monte_miser_alloc(dim[pop]);
					gsl_monte_miser_integrate (&G, xl, xu, dim[pop], calls, r2, s, &pop_lhood, &pop_err);
					gsl_monte_miser_free(s);

					lstorepop = log(pop_lhood);

					if(crash==1){
								lstorepop = -40000.0;
								printf("crashed in lhood file\n");
								//getc(stdin);
					}

				store2pop+=lstorepop;
				printf("store2pop = %e\n", store2pop);
				error += pop_err;
				//printf("calls = %d, lstorepop = %e, store2pop = %e\n", calls, lstorepop, store2pop);
				//getc(stdin);
			lstorepop = 0.0;
			}//close loop over populations

			fp = fopen("G1.dat","a");
			fprintf(fp, "call = %d, parset = %d, likelihood = %e, error = %e\n", q, r, store2pop, error);
			fclose(fp);
			store2pop = 0.0;

		}//close loop over parameter sets

	}//close loop over calls
		exit(1);
}//close verbose*/



srand((unsigned)time(NULL));

run=1;	changer=1;	best_post_hood=-50000;
// -------------------------------- INNER MAIN LOOP (After 'pro' is Fixed) ------------------------------------ //
while (changer>0 && run<=num_runs)	{		//open inner loop

	run ++;		changer=0;
	//printf("run = %d\n", run);
	j = gsl_rng_uniform_int(r2, num_adj_pars);
	for	(ii=j;ii<num_adj_pars+j;ii++)			{ //Looping over parameters
			i=1+ii%num_adj_pars;       //CK// converts from parameter number in the loop to an actual parameter in your list of parameters
			if (i==pro)		{	//printf("PROFILE PARAMETER,MOVE TO NEXT ONE\n") //old and dead
		}
		else if (i==4||i==7||i==8||i==9||i==18||i==33||i==34)	{

					for (inner_parm=Params.parm_low[i];inner_parm<=Params.parm_high[i];inner_parm+=Params.parm_step[i])	{
						if (i==4||i==7||i==8||i==18||i==33||i==34)	{	Params.PARS[i] = pow(10,inner_parm);	}

						else								{	Params.PARS[i] = round(inner_parm);	}

				total_lhood=0;
Params.pop = 0;
crash = 0;
redo = 0;
error = 0.0;
//printf("i = %d\n", i);

	while (Params.pop<22 && crash<1 )	{
		error = 0.0;

		pop=Params.pop;
		//printf("pop = %d\n", pop);
		pop_best_lhood = -1e9;
		int MMAX=Params.MAXT[pop];   //MG//
		MMAX=(MMAX/7);
		int MAXT3=(Params.DATA[pop][MMAX][3])*7;
		gsl_monte_function G = {&Hood_Pops, dim[pop], &Params };	// declares function calling Hood_Pops.h
		//////

		double xl[dim[pop]];	double xu[dim[pop]];	// need to redeclare xl and xu since the size changes

		for (jj=0;jj<=dim[pop];jj++)	{
					xl[jj]=0;
					xu[jj]=1;
				}


			// ----------- Use MISER to call function pop_lhood --------------------------------- //
		gsl_monte_miser_state *s = gsl_monte_miser_alloc(dim[pop]);


		gsl_monte_miser_integrate (&G, xl, xu, dim[pop], calls, r2, s, &pop_lhood, &pop_err);
		gsl_monte_miser_free(s);
		pop_lhood2 = log(pop_lhood);
		if(redo > 0){
			crash = 1;
			}
		if(crash==1){
			pop_lhood2 = -80000.0;
			//printf("crashed in lhood file\n");
		}
		error += pop_err;
		//printf("error = %e\n", error);
		if(error>0.001){
		//	printf("fucking liars\n");
			crash = 1;
			pop_lhood2 = -20000.0;
			//getc(stdin);
		}

		if(pop_lhood2 == 0 || pop_lhood2 > 0 || isinf(pop_lhood2)== 1 || isnan(pop_lhood2)== 1 || isinf(-pop_lhood2)== 1 ){
								//printf("badness\n");
								crash = 1;
								pop_lhood2 = -30000.0;
						}

		pop_best_lhood = pop_lhood2;

		total_lhood += pop_best_lhood;

		Params.pop++;
		//printf("popnext = %d, total_lhood = %e\n", Params.pop, total_lhood);

				} //CK// end of going through patches
				post_hood=total_lhood;
				//printf("post_hood = %e\n", post_hood);

				if (post_hood > best_post_hood)	{


								best_post_hood = post_hood;

								best_lhood= total_lhood;				// max lhood associated with this posterior
								Params.MLE[i] = Params.PARS[i];
								Params.PARS[3] = error;


					changer++;
					//printf("best_post_hood = %d \n", best_post_hood);
				}

			}		//CK// END OF LINE SEARCH FOR SELECTED PARAMETER

			//printf("parameter: %d\t has MLE value=%e\n",i,Params.MLE[i]);
			Params.PARS[i]=Params.MLE[i];     //////CK///////!!!!!!!!  Found the spot where it puts the best value back into PARS for continued checking

		}		//CK/  END OF IF STATEMENT FOR SPECIFIC PARAM NUMBERS
	}		//CK// END OF LOOP FOR TICKING THROUGH PARAM NUMBERS FOR LINE SEARCHES
}		//CK// END OF INNER LOOP!!
gsl_rng_free(r2);


// ------------------------------------------ output results to file  --------------------------------------- //
molly_results = fopen(strFileName,"a");
output_file(&Params,molly_results,best_post_hood,num_adj_pars,num_runs, parm_inc); // prints to output file (filenames.h)
fclose(molly_results);

// ------------------------------------------------------------------------------------------------------- //
free_i3tensor(Params.DATA,0,23,0,37,0,8);
free_i3tensor(Params.EXPDATA,0,7,0,37,0,12);
free_d3tensor(Params.HostFactors,0,23,0,37,0,2);
free_d3tensor(Params.ExpFactors,0,7,0,37,0,2);
free_d3tensor(Params.FeralPara,0,23,0,37,0,4);
free_d3tensor(Params.ExpPara,0,7,0,37,0,4);
}  // closes main
