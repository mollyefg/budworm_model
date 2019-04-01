#include "headN1.h"
gsl_rng *r2;
int verbose = 2;
int exrealz = 2;
int redo = 0;
int crash = 0;

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


if (test==66)	{	    num_runs=3;	calls = 5; parm_inc= 50.0; host_inc=19.0;} //MG//  low resolution grid search
if (test==99)	{	    num_runs=5;	calls = 12; parm_inc=50.0; host_inc=19.0;} //MG//  low resolution grid search




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

double CHEAT[25] = {9.121368e-01,	9.997698e-04,	3.054675e-01,	1.585112e+00,	2.511771e-01,	21.0,	3.483373e-04,	4.655861e-01,	5.013026e+00,	7.000000e+00,	3.802069e-02,	4.570250e-03,	1.995722e+00,	1.000000e+00,	21.0,	2.512002e-01,	8.879259e-01,	5.755194e-01,	9.997698e-05,	3.145405e+00,	1.000000e+00,	1.047273e-01,	1.000000e+00,	1.513701e-02,	41.0};

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

double VALUE[5] = {36.467599, 837.455331, 0.817960, 1.000000, 0.453066};
//double VALUE[5] = {5.470245,	107.6469,	0.7870364,	1.0,	0.6428215};
//double VALUE[5] =  {210.385255, 735.123690, 0.778470, 7.000000, 0.000003};

double RandNumsPass[50];
int q;
int r;
int rlzns;
rlzns = 1;
double error = 0.0;

if(verbose==1){
	for(q = 1; q<2;q++){
			FILE *fp;
			Params.PARS[4] = VALUE[0];
			Params.PARS[7] = VALUE[1];
			Params.PARS[8] = VALUE[2];
			Params.PARS[9] = VALUE[3];
			Params.PARS[18] = VALUE[4];

			for(Params.pop=0;Params.pop<22;Params.pop++){
				for(i=0;i<dim[Params.pop];i++){
					RandNumsPass[i] = gsl_ran_flat(r2,0.0,1.0);
					}
				Hood_Pops(RandNumsPass,dim[Params.pop],Params.PARS);


				}
			}
		exit(1);
	}

if(exrealz==1){
	double liketest[rlzns][10];
	for(q = 1; q<=rlzns;q++){
		FILE *fp;

		for(r = 0; r < 1; r++){
			printf("r = %d \n", r);
			num_runs = 1;
			double lstorepop;
			double store2pop;

			Params.PARS[4] = VALUE[0];
			Params.PARS[7] = VALUE[3];
			Params.PARS[8] = VALUE[4];
			Params.PARS[9] = VALUE[5];
			Params.PARS[18] = VALUE[6];



		for(Params.pop=0;Params.pop<22;Params.pop++){
			pop = Params.pop;
			//printf("pop = %d \n", pop);
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
				store2pop+=lstorepop;
				//printf("calls = %d, lstorepop = %e, store2pop = %e\n", calls, lstorepop, store2pop);
				//getc(stdin);
			lstorepop = 0.0;
			}//close loop over populations

			//fp = fopen("modelN1_2Feb.dat","a");
			//fprintf(fp, "m1 = %e, calls = %d, likelihood = %e\n", Params.PARS[9], calls, store2pop);
			//fclose(fp);
			store2pop = 0.0;

		}//close loop over parameter sets

	}//close loop over calls
		exit(1);
}//close verbose
srand((unsigned)time(NULL));

run=1;	changer=1;	best_post_hood=-50000;
// -------------------------------- INNER MAIN LOOP (After 'pro' is Fixed) ------------------------------------ //
while (changer>0 && run<=num_runs)	{		//open inner loop

	run ++;		changer=0;
	j = gsl_rng_uniform_int(r2, num_adj_pars);
	for	(ii=j;ii<num_adj_pars+j;ii++)			{ //Looping over parameters
			i=1+ii%num_adj_pars;       //CK// converts from parameter number in the loop to an actual parameter in your list of parameters
			if (i==pro)		{	//printf("PROFILE PARAMETER,MOVE TO NEXT ONE\n") //old and dead
		}
		else if (i==4||i==7||i==8||i==9||i==18)	{

					for (inner_parm=Params.parm_low[i];inner_parm<=Params.parm_high[i];inner_parm+=Params.parm_step[i])	{
						if (i==4||i==7||i==8||i==18)	{	Params.PARS[i] = pow(10,inner_parm);	}

						else								{	Params.PARS[i] = round(inner_parm);	}

total_lhood=0;
Params.pop = 0;
crash = 0;
error = 0.0;
redo = 0;
	while (Params.pop<22 && crash<1 )	{
		pop=Params.pop;
		//printf("pop = %d\n", pop);
		//printf("i = %d, Params.PARS[i] = %e\n", i, Params.PARS[i]);
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
			pop_lhood2 = -40000.0;
		}

		error += pop_err;

		if(error>0.001){
			crash = 1;
			pop_lhood2 = -20000.0;
		}

		if(pop_lhood2 == 0 || pop_lhood2 > 0 || isinf(pop_lhood2)== 1 || isnan(pop_lhood2)== 1 || isinf(-pop_lhood2)==1 ){
								crash = 1;
								pop_lhood2 = -30000.0;
						}

		pop_best_lhood = pop_lhood2;

		total_lhood += pop_best_lhood;

		Params.pop++;
	}
				post_hood=total_lhood;

				if (post_hood > best_post_hood)				{//besthood

								best_post_hood = post_hood;
								best_lhood= total_lhood;				// max lhood associated with this posterior
								Params.MLE[i] = Params.PARS[i];
								Params.PARS[3] = error;

					changer++;
					//printf("changer = %d \n", changer);
				}//bestthood

			}		//CK// END OF LINE SEARCH FOR SELECTED PARAMETER

			//printf("parameter: %d\t has MLE value=%e\n",i,Params.MLE[i]);
			Params.PARS[i]=Params.MLE[i];     //////CK///////!!!!!!!!  Found the spot where it puts the best value back into PARS for continued checking

		}		//CK/  END OF IF STATEMENT FOR SPECIFIC PARAM NUMBERS
	}		//CK// END OF LOOP FOR TICKING THROUGH PARAM NUMBERS FOR LINE SEARCHES
}		//CK// END OF INNER LOOP!!
gsl_rng_free(r2);


// ------------------------------------------ output results to file  --------------------------------------- //
molly_results = fopen(strFileName,"a");
output_file(&Params,molly_results,best_post_hood,num_adj_pars,num_runs,parm_inc); // prints to output file (filenames.h)
fclose(molly_results);
//printf("closed\n");
// ------------------------------------------------------------------------------------------------------- //
free_i3tensor(Params.DATA,0,23,0,37,0,8);
free_i3tensor(Params.EXPDATA,0,7,0,37,0,12);
free_d3tensor(Params.HostFactors,0,23,0,37,0,2);
free_d3tensor(Params.ExpFactors,0,7,0,37,0,2);
free_d3tensor(Params.FeralPara,0,23,0,37,0,4);
free_d3tensor(Params.ExpPara,0,7,0,37,0,4);
}  // closes main
