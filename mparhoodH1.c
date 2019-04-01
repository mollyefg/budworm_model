#include "headH1.h"
gsl_rng *r2;
int verbose = 2;
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


if (test==66)	{	    num_runs=3;	calls = 4; parm_inc= 10.0; host_inc=19.0;} //MG//  low resolution grid search
if (test==99)	{	    num_runs=10;	calls = 10; parm_inc=40.0; host_inc=19.0;} //MG//  low resolution grid search




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

if(y>0.40){
//printf("random\n");
double z;
for(i = 1; i <= num_adj_pars; i++){
	if(i == 3 || i == 4 || i == 5 || i == 6 || i == 7 || i == 8 || i == 11 || i == 12 || i == 13 || i == 18 || i ==19 || i ==21 || i ==23 || i == 25|| i ==29||i==30||i==31||i==32||i==33||i==35){
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


//double CHEAT[25] = {5.011988e+00,	9.997698e-04,	1.023293e+00,	6.457880e+00,	4.216844e-01,	2.100000e+01,	0.000000e+00,	0.000000e+00,	1.328159e-01,	7.000000e+00,	6.309428e-03,	0.000000e+00,	0.000000e+00,	0.000000e+00,	0.000000e+00,	0.000000e+00,	0.000000e+00,	0.000000e+00,	0.000000e+00,	1.000000e+00,	0.000000e+00,	0.000000e+00,	0.000000e+00,	0.000000e+00,	0.000000e+00};
//-1.743623e+02,

double CHEAT[25] = {106.1802,	0.1984822,	0.07063172,	5650.509,	0.5011699,	1.0,	1.748988,	0.2483177,	1.0,	13.0,	0.1644602,2929.215,	6047.183,	16.0,	19.0,	6724.44,	375.7131,	0.4819336,	0.1208635,	0.6415618,	8.0,	0.6482538,	4.0,	0.473217,	17.0};

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
Params.PARS[19] = CHEAT[0];
Params.PARS[21] = CHEAT[0];
Params.PARS[27] = CHEAT[5];
Params.PARS[28] = CHEAT[5];
Params.PARS[29] = CHEAT[3];
Params.PARS[30] = CHEAT[3];
Params.PARS[31] = CHEAT[10];
Params.PARS[32] = CHEAT[10];
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
Params.MLE[19] = CHEAT[0];
Params.MLE[21] = CHEAT[0];
Params.MLE[27] = CHEAT[5];
Params.MLE[28] = CHEAT[5];
Params.MLE[29] = CHEAT[3];
Params.MLE[30] = CHEAT[3];
Params.MLE[31] = CHEAT[10];
Params.MLE[32] = CHEAT[10];
Params.MLE[33] = CHEAT[19];
Params.MLE[34] = CHEAT[20];
Params.MLE[35] = CHEAT[21];
Params.MLE[36] = CHEAT[22];
Params.MLE[23] = CHEAT[23];
Params.MLE[24] = CHEAT[24];
}


double VALUE[6] = {(1.862181e+01)/2,	(1.516263e+03)/2,	2.630092e-01,	1.000000e+00,	1.117120e-01,	5.754903e-02};
//-1.685299e+02,

double RandNumsPass[50];
int q;

if(exrealz==1){
	for(q = 0; q<10;q++){
		printf("q = %d\n", q);
		FILE *fp;
				Params.PARS[4] = VALUE[0];
				Params.PARS[5] = VALUE[1];
				Params.PARS[6] = VALUE[2];
				Params.PARS[7] = VALUE[3];
				Params.PARS[8] = VALUE[4];
				Params.PARS[9] = VALUE[5];
				//Params.PARS[11] = VALUE[6];
				//Params.PARS[12] = VALUE[7];
				Params.PARS[13] = VALUE[6];
				//Params.PARS[14] = VALUE[7];
				Params.PARS[18] = VALUE[7];
				Params.PARS[19] = VALUE[0];
				Params.PARS[27] = VALUE[5];
				Params.PARS[29] = VALUE[3];
				Params.PARS[31] = VALUE[7];
				//Params.PARS[33] = VALUE[8];
				//Params.PARS[34] = VALUE[7];
				//Params.PARS[35] = VALUE[6];
				//Params.PARS[36] = VALUE[7];

				Params.PARS[21] = Params.PARS[4];
				Params.PARS[28] = Params.PARS[9];
				Params.PARS[30] = Params.PARS[7];
				Params.PARS[32] = Params.PARS[18];
				//Params.PARS[23] = Params.PARS[13];
				//Params.PARS[24] = Params.PARS[14];

		//fp = fopen("modelD_single.dat","a");
		//fclose(fp);

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
double error = 0.0;
if(verbose==1){
	double liketest[rlzns][10];
	for(q = 1; q<=rlzns;q++){
		FILE *fp;
		//fp = fopen("modelD_single_nosumm.dat","a");
		//fprintf(fp, "call =, parset, likelihood\n");
		//fclose(fp);
		for(r = 0; r < 5; r++){
			printf("r = %d \n", r);
			double lstorepop;
			double store2pop;
			calls = 50;
			error = 0.0;
				Params.PARS[4] = VALUE[0];
				//Params.PARS[5] = VALUE[1];
				//Params.PARS[6] = VALUE[2];
				Params.PARS[7] = VALUE[1];
				Params.PARS[8] = VALUE[2];
				Params.PARS[9] = VALUE[3];
				//Params.PARS[11] = VALUE[6];
				//Params.PARS[12] = VALUE[7];
				Params.PARS[13] = VALUE[4];
				//Params.PARS[14] = VALUE[7];
				Params.PARS[18] = VALUE[5];
				//Params.PARS[19] = VALUE[0];
				//Params.PARS[27] = VALUE[5];
				//Params.PARS[29] = VALUE[3];
				//Params.PARS[31] = VALUE[7];
				//Params.PARS[33] = VALUE[8];
				//Params.PARS[34] = VALUE[7];
				//Params.PARS[35] = VALUE[6];
				//Params.PARS[36] = VALUE[7];

				//Params.PARS[21] = Params.PARS[4];
				//Params.PARS[28] = Params.PARS[9];
				//Params.PARS[30] = Params.PARS[7];
				//Params.PARS[32] = Params.PARS[18];
				//Params.PARS[23] = Params.PARS[13];
				//Params.PARS[24] = Params.PARS[14];

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
					//if(pop==15){
					//printf("lhood_check %e, pop_lhood = %e, error = %e\n", lhood_check/calls, pop_lhood, pop_err);
					if(crash==1){
						//printf("after miser, crash is one\n");
						//getc(stdin);
					}

					//getc(stdin);
					//}
					//lhood_check = 0.0;
					//getc(stdin);
					//printf("error= %e\n", pop_err);
					error+=pop_err;

					//getc(stdin);
					gsl_monte_miser_free(s);
					lstorepop = log(pop_lhood);
				store2pop+=lstorepop;
				//printf("calls = %d, lstorepop = %e, store2pop = %e\n", calls, lstorepop, store2pop);
									//fp = fopen("fuckinglikelihood.dat","a");
									//fprintf(fp, "pop = %d, error = %e, likelihood = %e, loglikelihood = %e\n", pop, pop_err, pop_lhood, lstorepop);
					//fclose(fp);
			lstorepop = 0.0;
			//getc(stdin);
			}//close loop over populations

			fp = fopen("23jan_D_single_verbose.dat","a");
			fprintf(fp, "calls = %d, error = %e, likelihood = %e\n", calls, error, store2pop);
			fclose(fp);
			//printf("error = %e, parset = 'good', total lhood = %e\n", error, store2pop);
			//getc(stdin);

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
	//printf("run = %d\n", run);
	j = gsl_rng_uniform_int(r2, num_adj_pars);
	for	(ii=j;ii<num_adj_pars+j;ii++)			{ //Looping over parameters
			i=1+ii%num_adj_pars;       //CK// converts from parameter number in the loop to an actual parameter in your list of parameters
			if (i==pro)		{	//printf("PROFILE PARAMETER,MOVE TO NEXT ONE\n") //old and dead
		}
		else if (i==4||i==7||i==8||i==9||i==13||i==18||i==34)	{

					for (inner_parm=Params.parm_low[i];inner_parm<=Params.parm_high[i];inner_parm+=Params.parm_step[i])	{
						if (i==4||i==7||i==8||i==13||i==18||i==34)	{	Params.PARS[i] = pow(10,inner_parm);	}

						else								{	Params.PARS[i] = round(inner_parm);	}

total_lhood=0;
Params.pop = 0;
crash = 0;
redo = 0;
error = 0.0;

	while (Params.pop<22 && crash<1 )	{
		pop=Params.pop;
		//printf("pop = %d, crash = %d, error = %e, best_post_hood = %e\n", pop, crash, error, best_post_hood);
		//printf("pop = %d, run%d\n", pop, run);


		pop_best_lhood = -1e9;
		int MMAX=Params.MAXT[pop];   //MG//
		MMAX=(MMAX/7);
		int MAXT3=(Params.DATA[pop][MMAX][3])*7;
		gsl_monte_function G = {&Hood_Pops, dim[pop], &Params };	// declares function calling Hood_Pops.h
		//////
		//printf("before miser: crash = %d\n", crash);
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
			//printf("crashed in lhood file\n");
		}

		error += pop_err;
		//printf("error = %e\n", error);

		//printf("after miser: pop = %d, crash = %d, error = %e\n", pop, crash, error);
		//getc(stdin);
		if(error>0.001){
			//printf("fucking liars\n");
			//printf("pop_lhood2 = %e\n", pop_lhood2);
			crash = 1;
			pop_lhood2 = -20000.0;
			//getc(stdin);
		}

		//if(pop_lhood2 == 0 || pop_lhood2 > 0 || isinf(pop_lhood2)== 1 || isnan(pop_lhood2)== 1 ){
		if(pop_lhood2 == 0 || pop_lhood2 > 0 || isinf(pop_lhood2)== 1 || isnan(pop_lhood2)== 1 || isinf(-pop_lhood2)==1 ){
								//printf("pop_lhood = %e, pop_lhood2 = %e\n", pop_lhood, pop_lhood2);
								//printf("badness\n");
								crash = 1;
								pop_lhood2 = -30000.0;
								//getc(stdin);
						}

		pop_best_lhood = pop_lhood2;

		total_lhood += pop_best_lhood;

		//printf("pop = %d, pop_lhood = %e, poplhood2 = %e, total_lhood = %e, error = %e\n",pop, pop_lhood, pop_lhood2, total_lhood, pop_err);
		//getc(stdin);

		Params.pop++;
		//printf("popnext = %d, crash = %d, total_lhood = %e, best_post_hood = %e\n", Params.pop, crash, total_lhood, best_post_hood);

		//getc(stdin);
				} //CK// end of going through patches
			//	printf("outside loop: error = %e, crash = %d\n", error, crash);

				post_hood=total_lhood;

				if (post_hood > best_post_hood)				{

								//printf("change: post_hood = %e, best_post_hood = %e, error = %e\n", post_hood, best_post_hood, error);
								//getc(stdin);
								best_post_hood = post_hood;
								//molly_results = fopen("R3_single_likescores","a");
								//output_file(&Params,molly_results,best_post_hood,num_adj_pars,num_runs, calls); // prints to output file (filenames.h)
								//fclose(molly_results);

								best_lhood= total_lhood;				// max lhood associated with this posterior
								Params.MLE[i] = Params.PARS[i];
								Params.PARS[3] = error;

					changer++;
					//printf("changer = %d \n", changer);
				}

			}		//CK// END OF LINE SEARCH FOR SELECTED PARAMETER

			//printf("parameter: %d\t has MLE value=%e\n",i,Params.MLE[i]);
			Params.PARS[i]=Params.MLE[i];     //////CK///////!!!!!!!!  Found the spot where it puts the best value back into PARS for continued checking

		}		//CK/  END OF IF STATEMENT FOR SPECIFIC PARAM NUMBERS
	}		//CK// END OF LOOP FOR TICKING THROUGH PARAM NUMBERS FOR LINE SEARCHES
}		//CK// END OF INNER LOOP!!
gsl_rng_free(r2);

/*if(test==66){
printf("model = U3_single\n");
printf("alpha1 = %e, beta = %e, variance = %e\n",Params.PARS[4], Params.PARS[5], Params.PARS[6]);
//printf("alpha2 = %e, alpha3 = %e\n", Params.PARS[19], Params.PARS[21]);
printf("rho1 = %e, phi = %e, m1 = %e\n", Params.PARS[7], Params.PARS[8], Params.PARS[9]);
printf("delta1 = %e\n", Params.PARS[18]);
printf("gamma1 = %e, gamma2 = %e\n", Params.PARS[13], Params.PARS[14]);
//printf("rho2 = %e, m2 = %e, delta2 = %e\n", Params.PARS[29], Params.PARS[27], Params.PARS[31]);
//printf("rho3 = %e, m3 = %e, delta3 = %e\n", Params.PARS[30], Params.PARS[28], Params.PARS[32]);
//printf("gamma3 = %e, gamma4 = %e \n", Params.PARS[33], Params.PARS[34]);
//gamma5 = %e, gamma6 = %e, gamma7 = %e, gamma8 = %e\n",Params.PARS[35], Params.PARS[36], Params.PARS[23], Params.PARS[24]);

printf("LHOOD:%e\n",best_post_hood);
}*/

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
