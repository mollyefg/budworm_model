#include "headH1.h"
int verbose = 2;
int exrealz = 2;
gsl_rng *r;
int crash = 0;
int WAIC = 1;
int redo = 0;

//////////////////Begin Dot Product////////////////////////////////////////

float DotProduct (int Length, double *Holder, double *PCA)
{

  double answer = 0;

  int i;
  for(i=0;i<Length;i++)
    answer += PCA[i]*Holder[i];
  return(answer);

}

//////////////////End dot product/////////////////////////////////////////



int main(int argc, char *argv[])
{
//int test = 66;
//int test = 99;
int test = 0;
STRUCTURE Params;
int pro = 1;//atoi(argv[1]);						// pro and argv[1] are the inputs (argv[i] is the i^th input)

// ------------------------------------- Adustable accuracy vs. speed ------------------------------------------------ //
int num_runs	= 1;
double parm_inc;
int calls=2;					// number of stochastic simulations for each parameter and IC set; miser averages over this
int Realizations=1000000;

parm_inc=500.0;	//put in some rounding for fuck's sake

if (test==66)	{calls = 5;	Realizations = 150; parm_inc=10.0;}
if (test==99)   {calls = 15; Realizations = 100000000; parm_inc=50.0;}

//printf("runs=%d, calls = %d, rlzns= %d\n",num_runs,calls, Realizations);

// ------------------------------------------------------------------------------------------------------------------ //
int i=0; int j;int ii; int jj; int k;
int run;	            int changer;	    double index, tot_index;

int num_adj_pars=40;			// number of adjustable parameters
double inner_parm;				//double outer_parm;
double inner_parm2;				//double outer_parm;
int pop;
double log_pop;

// -------------------------------------------- MISER STUFF --------------------------------------------------------- //

inputdata(&Params);				// gets Params.DATA[j][i][0-2] and Params.MAXT[i] from inputdata.h

size_t dim[22]={13, 22, 15, 29, 36, 29, 29, 29, 29, 35, 29, 35, 36, 35, 36, 36, 35, 39, 44, 44, 39, 39};  //MG  length of time observed per site

// --------------------------------------- Name for Output Files ----------------------------------------------------- //
char strFileName[99];					// from filenames.h
GetString(pro,0,strFileName,98);
fflush(stdout);		//getc(stdin);
FILE *molly_results;

// ---------------------------------------- Random Number Stuff ------------------------------------------------------ //
//MG  try replacing with better random stuff?
gsl_rng *r_seed;
r_seed=random_setup();


const gsl_rng_type *T;
long seed;
seed = time(NULL)*(int)getpid();  //use process id
//seed = -1;
gsl_rng_env_setup ();
gsl_rng_default_seed = seed;

T = gsl_rng_default;
r = gsl_rng_alloc(T);

//printf("seed = %d\n", seed);

double host_inc = 19.0; double initR_inc = 10.0; //MG  what are these even
parm_range_inc(&Params,parm_inc,host_inc,initR_inc,num_adj_pars); // gets Params.parm_set,low,high,R_END from bounds.h

// ------------------------------------ Declare Likelihood Quanitites ----------------------------------------------- //
double pop_lhood, pop_lhood2, pop_err,post_hood;	// population lhood (and posterior lhood) calculated for each initS and initR
double pop_best_lhood;					// likelihood and error for best initS and initR
double total_lhood;						// sum of pop_best_lhood over all patches
double best_post_hood;	double best_lhood=0;		// best post_hood and lhood
double prior[num_adj_pars];

double error;


//******************* Just for WAIC calculation **********************////

  if(WAIC==1){ //Turn on WAIC calculation


int maxweeks = FERAL_SETS;

int pop;
    FILE *fp4;

    int NumParamsIn = 25; //MG
    int NumSets = 115;

    double **PosteriorsIn;
    PosteriorsIn = dmatrix(0, NumSets, 0, NumParamsIn);


    // *************** Change file name here *****************

    i = 0;
    int MaxSEIRParamNum = 0;
    int SEIRParamNum = 0;

    fp4 = fopen("H1 WAIC.in", "r");
    //printf("finished opening file\n");

    //printf("begin write file\n");
    while(fscanf(fp4,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
               &(PosteriorsIn[SEIRParamNum][0]),&(PosteriorsIn[SEIRParamNum][1]),&(PosteriorsIn[SEIRParamNum][2]),&(PosteriorsIn[SEIRParamNum][3]),&(PosteriorsIn[SEIRParamNum][4]),
			   				                 &(PosteriorsIn[SEIRParamNum][5]),&(PosteriorsIn[SEIRParamNum][6]),&(PosteriorsIn[SEIRParamNum][7]),&(PosteriorsIn[SEIRParamNum][8]),
			   				                 &(PosteriorsIn[SEIRParamNum][9]),&(PosteriorsIn[SEIRParamNum][10]),&(PosteriorsIn[SEIRParamNum][11]),&(PosteriorsIn[SEIRParamNum][12]),
			   				                 &(PosteriorsIn[SEIRParamNum][13]),&(PosteriorsIn[SEIRParamNum][14]),&(PosteriorsIn[SEIRParamNum][15]),&(PosteriorsIn[SEIRParamNum][16]),
			   				                 &(PosteriorsIn[SEIRParamNum][17]),&(PosteriorsIn[SEIRParamNum][18]),&(PosteriorsIn[SEIRParamNum][19]),&(PosteriorsIn[SEIRParamNum][20]),
			                   				 &(PosteriorsIn[SEIRParamNum][21]),&(PosteriorsIn[SEIRParamNum][22]),&(PosteriorsIn[SEIRParamNum][23]),
                				 &(PosteriorsIn[SEIRParamNum][24]))!=EOF){
                 SEIRParamNum++;
   		 }
    fclose(fp4);


   // MaxSEIRParamNum = 115;
	  MaxSEIRParamNum = 115;

    //Randomize the number of calls:
   	//calls = (size_t) (gsl_rng_uniform(r) * 300 + 1);
	//printf("calls = %d\n", calls);

    calls=3;

    double PopWiseMean[maxweeks], PopWiseM2Log[maxweeks], PopWiseVarLog[maxweeks],PopWiseMeanLog[maxweeks];
    double tempwisemean[maxweeks], tempwisem2log[maxweeks], tempwisevarlog[maxweeks], tempwisemeanlog[maxweeks];

    for(i=0;i<maxweeks;i++){
      PopWiseMean[i] = 0.0;
      PopWiseM2Log[i] = 0.0;
      PopWiseMeanLog[i] = 0.0;
      PopWiseVarLog[i] = 0.0; //unnecessary?
    }

	int ticker2 = 0;
	int ticker3 = 0;
    for(SEIRParamNum=0;SEIRParamNum<MaxSEIRParamNum;SEIRParamNum++){
	printf("parm = %d\n", SEIRParamNum);
	redo = 0;
	//getc(stdin);
   		 for(i=0;i<maxweeks;i++){
     		 tempwisemean[i] = 0.0;
      		 tempwisem2log[i] = 0.0;
     		 tempwisevarlog[i] = 0.0;
     		 tempwisemeanlog[i] = 0.0; //unnecessary?
   		  }

   	  Params.PARS[4] = PosteriorsIn[SEIRParamNum][0];	//alpha1
      //Params.PARS[5] = PosteriorsIn[SEIRParamNum][1];	//beta
      //Params.PARS[6] = PosteriorsIn[SEIRParamNum][2];	//sigma
      Params.PARS[7] = PosteriorsIn[SEIRParamNum][3];	//rho1
      Params.PARS[8] = PosteriorsIn[SEIRParamNum][4];	//phi

      ////FIX m1 at 1
      //Params.PARS[9] = PosteriorsIn[SEIRParamNum][5];	//m1
		Params.PARS[9] = 1.0;

      Params.PARS[13] = PosteriorsIn[SEIRParamNum][8];	//gamma1
      //Params.PARS[14] = PosteriorsIn[SEIRParamNum][9];	//gamma2

      Params.PARS[18] = PosteriorsIn[SEIRParamNum][10];	//delta1
     // Params.PARS[18] = 1.00023;	//delta1

     // Params.PARS[33] = PosteriorsIn[SEIRParamNum][19];	//delta1
     Params.PARS[34] = PosteriorsIn[SEIRParamNum][20];	//delta1

	//printf("alpha1 = %e, beta = %e, var = %e, rho1 = %e, phi = %e, m1 = %e, gamma1 = %e, gamma2 = %e, delta1 = %e\n", Params.PARS[4],Params.PARS[5],Params.PARS[6],Params.PARS[7],Params.PARS[8],Params.PARS[9],Params.PARS[13],Params.PARS[14],Params.PARS[18]);
		//double LHoodPerPop = 0.0;
		//double AvgLHood = 0.0;

		int ticker = 0;
		error = 0.0;
		crash = 0;
      for(pop=0;pop<FERAL_SETS;pop++){
		Params.pop = pop; //MG  this is so important
				//printf("pop = %d\n", pop);

		//printf("Params.pop = %d\n", Params.pop);
					gsl_monte_function G = { &Hood_Pops, dim[pop], &Params };	// declares function calling Hood_Pops.h
					double xl[dim[pop]];	double xu[dim[pop]];	// need to redeclare xl and xu since the size changes
					for (jj=0;jj<=dim[pop];jj++)	{
							xl[jj]=0;
							xu[jj]=1;
							}
					// ----------- Use MISER to call function pop_lhood --------------------------------- //
			//printf("204\n");
					gsl_monte_miser_state *s = gsl_monte_miser_alloc(dim[pop]);
					gsl_monte_miser_integrate (&G,xl,xu,dim[pop],calls,r,s,&pop_lhood,&pop_err);
					gsl_monte_miser_free(s);
		//printf("208\n");
					error += pop_err;
			 		double res = pop_lhood;
					        			//printf("pop = %d, loop = %d, res = %e\n", pop, SEIRParamNum, res);

						if(crash==1){
						redo = 1;
						//printf("crashed in lhood file\n");
							}

						if(error>0.0025){
						redo = 1;
						crash = 1;
						//printf("error failure\n");
							}

						if(log(pop_lhood) == 0 || log(pop_lhood) > 0 || log(pop_lhood)== 1 || log(pop_lhood)== 1 || log(pop_lhood)== 1 ){
								redo = 1;
								crash = 1;
								//printf("log(pop_lhood) fail = %e\n", log(pop_lhood));
								//getc(stdin);
							}

		if(redo > 0 && ticker > 10){
						ticker = 0;
						//redo = 0;
						//LHoodPerPop = 0.0;
						error = 0.0;
						crash = 0;
		 				ticker2 = ticker2 ++;
		 				//printf("ticker2 = %d\n", ticker2);
		 				MaxSEIRParamNum = MaxSEIRParamNum + 1;
		 				//printf("MaxSEIRParamNum = %d\n", MaxSEIRParamNum);
						//printf("about to break\n");
		 				break;
				}

		if(redo > 0 && ticker <=10 ){
			//printf("fail\n");
			redo = 0;
			crash = 0;
			//LHoodPerPop = 0.0;
			pop = -1;
			Params.pop = 0;
			error = 0.0;
			ticker = ticker ++;
			//printf("ticker = %d\n", ticker);

		 }//if statement
				if(redo == 0){
					//printf("redo = 0\n");

			         //LHoodPerPop += log(res);

					tempwisemean[pop] = res;
					tempwisemeanlog[pop] = log(res);
					tempwisem2log[pop] = log(res)*log(res);
					//printf("redo = %d, crash = %d\n", redo, crash);
					//printf("pop = %d, twml[pop] = %e, twm2l[pop] = %e \n", pop, tempwisemeanlog[pop], tempwisem2log[pop]);
  				  //FILE *fp6;

    			//fp6 = fopen("WAIC_rho1_components.dat","a");

    			//fclose(fp6);
					//printf("adding to temps\n");
			        					}

       // printf("pop = %d, LHoodPerPop: %f\n", pop, LHoodPerPop);
       // printf("PopWiseMean %e, PopeWiseMeanLog = %e, PopWiseM2Log = %e\n", tempwisemean[pop], tempwisemeanlog[pop], tempwisem2log[pop]);
      } //Pop loop
		//printf("redo = %d\n", redo);
		if(redo == 0){
		ticker3 ++;
		//printf("ticker3 = %d\n", ticker3);
		//printf("end of SEIRparamnum %d\n", SEIRParamNum);
		for(i = 0; i < FERAL_SETS; i++){
			pop = i;
			PopWiseMean[pop] += tempwisemean[pop];
			PopWiseMeanLog[pop] += tempwisemeanlog[pop];
       		PopWiseM2Log[pop] += tempwisem2log[pop];
       		//FILE *fp6;
    			//fp6 = fopen("WAIC_E1_components2.dat","a");
        			//fprintf(fp6,"pop = %d, loop = %d, popwisemean = %e, popwisemeanlog = %e, popwisem2log = %e\n", pop, SEIRParamNum, PopWiseMean[pop], PopWiseMeanLog[pop], PopWiseM2Log[pop]);
    			//fclose(fp6);
       		//printf("loop = %d, pop = %d, PopWiseMeanLog[pop] = %e, PopWiseM2Log[pop] = %e\n", SEIRParamNum, pop, PopWiseMeanLog[pop],PopWiseM2Log[pop]);
		}//end of additive pop loop
		}//end of if statement

    } //SEIRParamNum

	SEIRParamNum = MaxSEIRParamNum - ticker2;
	//printf("final SEIRParamNum = %d\n", SEIRParamNum);


    double pWAIC2B = 0.0;
    double TotPopWiseMeanLog = 0.0;

    for(pop=0;pop<FERAL_SETS;pop++){
      pWAIC2B += (PopWiseM2Log[pop]/SEIRParamNum) - (PopWiseMeanLog[pop]/SEIRParamNum)*(PopWiseMeanLog[pop]/SEIRParamNum);
     // printf("PWM2/SEIR = %e, PWML[pop]/SEIR = %e, SEIRParamNum = %d \n", PopWiseM2Log[pop]/SEIRParamNum, PopWiseMeanLog[pop]/SEIRParamNum, SEIRParamNum);
      TotPopWiseMeanLog += log(PopWiseMean[pop]/SEIRParamNum);
       printf("pop = %d, pWAIC2B = %e, TotPopWiseMeanLog = %e \n", pop, pWAIC2B, TotPopWiseMeanLog);
      //fprintf(fp5,"%d %d %d %lf %lf %lf %lf %lf \n", calls, ticker2, pop, pWAIC2B, PopWiseM2Log[pop]/SEIRParamNum,  - (PopWiseMeanLog[pop]/SEIRParamNum)*(PopWiseMeanLog[pop]/SEIRParamNum), TotPopWiseMeanLog, log(PopWiseMean[pop]/SEIRParamNum));
    }


    double AltWAIC = -2*TotPopWiseMeanLog + 2*pWAIC2B;
    /*FILE *fp5;
	 fp5 = fopen("WAIC_rho1_19jan.dat","a");
       fprintf(fp5,"%lf %lf %lf\n", AltWAIC, -2*TotPopWiseMeanLog, 2*pWAIC2B);
    fclose(fp5);*/

	//printf("AltWAIC = %e \n", AltWAIC);
	//getc(stdin);


    FILE *fp5;

    fp5 = fopen("WAIC_H1_31May.dat","a");
        fprintf(fp5,"%d %d %d %d %lf %lf %lf\n", calls, MaxSEIRParamNum, ticker2, SEIRParamNum, TotPopWiseMeanLog, AltWAIC, pWAIC2B);
    fclose(fp5);

    //printf("%d %d %d %lf %lf \n", calls, MaxSEIRParamNum, ticker2, TotPopWiseMeanLog, AltWAIC);

    free_dmatrix(PosteriorsIn, 0, NumSets, 0, NumParamsIn);

    //printf("exit\n");
    exit(1);
  } //if(0)

//******************* Just for WAIC calculation **********************////




/////////////////////////////Begin reading in principal component analysis results//////////////////////

	int NumberOfParams=6;			// 32 parameters - alpha1,beta,sigma,rho1,phi,m1,delta1
									// there's a list taped to your monitor
	double RandomNumber;

	double LogJumpToNew;
	double LogJumpToOld;
	double ProbOfAcceptance;

	double LogOldPosterior;
	double LogNewPosterior;
	double LogNewPrior;

	//const gsl_rng_type * T;

	//gsl_rng * r;
	//gsl_rng_env_setup();
	//T = gsl_rng_default;
	//r = gsl_rng_alloc (T);
//	gsl_rng_set(r, time(NULL));
	//srand(time(NULL));

	double SDpca[NumberOfParams];
	double Coefficients[NumberOfParams*NumberOfParams];
	double Center[NumberOfParams];
	double Scale[NumberOfParams];
	double PC[NumberOfParams];
	double Old_PC[NumberOfParams];
	double Holder[NumberOfParams];
	double PCAparams[NumberOfParams];
	double Old_Params[NumberOfParams];

	//MG  should these be on?
	//double AcceptedVect[NumberOfParams]={0.0};
	//double LoopVect[NumberOfParams]={0.0};

	int a;
	int b;
	int ticker;
	int ticker2;

	int Case;

	int Accepted=0;
	signed int LoopNumber=-1;

	int ParCnt2 = NumberOfParams;

	double SigmaInflation=2.0;  //CK// Dave had 4 as his sigma inflation factor

	double sigma[NumberOfParams];

	run=1;	changer=1;	best_post_hood=-10000000000;

/////////////////////////////Begin reading in principal component analysis results//////////////////////

	FILE *file;

	file=fopen("H1_sd.txt", "r");
	for (a=0; a<(NumberOfParams); a++)
	{
		fscanf(file, "%lf\n", &SDpca[a]);
		//printf("%lf\n", SDpca[a]);
	}
	fclose(file);

	file=fopen("H1_Rotations.txt", "r");
	for (a=0; a<(NumberOfParams*NumberOfParams); a++)
	{
		fscanf(file, "%lf\n", &Coefficients[a]);
		//printf("%lf\n", Coefficients[a]);
	}
	fclose(file);
	file=fopen("H1_Scale.txt", "r");
	for (a=0; a<NumberOfParams; a++)
	{
		fscanf(file, "%lf\n", &Scale[a]);
		//printf("%lf\n", Scale[a]);
	}
	fclose(file);

	file=fopen("H1_Center.txt", "r");

	for (a=0; a<NumberOfParams; a++)
	{
		fscanf(file, "%lf\n", &Center[a]);
		//printf("%lf\n", Center[a]);
	}
	fclose(file);

/////////////////////////////DONE reading in principal component analysis results//////////////////////


//////////////////////////////Start MCMC/////////////////////////////////////////////////////////////////////////////

///Generate first set of PC values////

	for(a=0; a<NumberOfParams; a++){
		sigma[a]=SigmaInflation*SDpca[a];
	}

	for (a=0; a<NumberOfParams; a++){
		PC[a]=gsl_ran_gaussian (r, sigma[a]);    //CK//  PC contains current set of PC values
		//printf("a = %d, %lf\n", a, PC[a]);
	}

///Begin Loop that goes through and tweaks PCs one at a time////

//while (1==1) {     //CK// INFINITE LOOP START!!!!
while (LoopNumber<=Realizations) {     //CK// BOUND LOOP START!!!!
		//printf("LoopNumber = %d\n", LoopNumber);
		LoopNumber=LoopNumber+1;
//MG  what is this for?  if (LoopNumber % NumberOfParams == 1)
//		{
		Case=LoopNumber%NumberOfParams;					//Determines which PC to change

		for (a=0; a<NumberOfParams; a++)
		{
			Old_PC[a]=PC[a];							//Store old PC values:
		}

		PC[Case]=gsl_ran_gaussian (r, sigma[Case]);		//Draw 1 new PC

		for (a=0;a<NumberOfParams; a++)								//Back transform PC's into model parameters
		{
			for (b=0; b<NumberOfParams; b++){
				Holder[b]=Coefficients[a*NumberOfParams+b];
					}


			PCAparams[a]=exp(DotProduct(ParCnt2, Holder, PC)*Scale[a]+Center[a]);
			Old_Params[a]=exp(DotProduct(ParCnt2, Holder, Old_PC)*Scale[a]+Center[a]);
		}

		//N1 parameters: 0    1    2    3    4    5     6
						//alpha1, 0

						//rho1, 1
						//phi, 2
						//m1, 3
						// gamma1, 4
						//delta1, 5
							//gamma4, 6

			///fixed m1:  alpha1 0, rho1 1, phi 2, gamma1 3, delta1 4, gamma4 5

		//PCAparams[3] = round(PCAparams[3]);
		//Old_Params[3] = round(Old_Params[3]);
		//if(PCAparams[3] < 1){
		//	PCAparams[3] = 1;
		//}
		//if(Old_Params[3] < 1){
		//	Old_Params[3] = 1;
		//}


		PCAparams[5] = round(PCAparams[5]);
		Old_Params[5] = round(Old_Params[5]);
		if(PCAparams[5] < 1){
			PCAparams[5] = 1;
		}
		if(Old_Params[5] < 1){
			Old_Params[5] = 1;
		}

		//printf("PCAparams=%e, OldParams=%e\n", PCAparams[5], Old_Params[5]);  //getc(stdin);

		LogJumpToNew= -log(gsl_ran_gaussian_pdf(PC[Case], sigma[Case]));		//Metropolis sampling step: LATER YOU WILL USE THIS TO CORRECT FOR PROPOSAL
		LogJumpToOld= -log(gsl_ran_gaussian_pdf(Old_PC[Case], sigma[Case]));

		//printf("New: %f Old: %f\n", LogJumpToNew, LogJumpToOld);

		//MG can I just do this?

	Params.PARS[4] = PCAparams[0];//alpha1
	//Params.PARS[5] = PCAparams[1];//beta
	//Params.PARS[6] = PCAparams[2];//variance
	Params.PARS[7] = PCAparams[1];//rho1
	Params.PARS[8] = PCAparams[2];//phi
	//Params.PARS[9] = PCAparams[3];//m1
	Params.PARS[13] = PCAparams[3];//gamma1
	Params.PARS[18] = PCAparams[4];//delta1
	Params.PARS[34] = PCAparams[5];//delta1

	///m1:
	Params.PARS[9] = 1.0;

for(i = 0; i < 37; i++){
	if(i == 4 || i == 7 || i == 8 || i == 13 || i == 18 || i == 19 || i == 21 || i == 29 || i == 30 || i == 31 || i == 32 || i == 33 || i == 35 || i == 23){
		if(Params.PARS[i] < pow(10,Params.parm_low[i])){
			Params.PARS[i] = pow(10,Params.parm_low[i]);
				}
		if(Params.PARS[i] > pow(10,Params.parm_high[i])){
			Params.PARS[i] = pow(10,Params.parm_high[i]);
				}
		}
	if(i == 9 || i == 14 || i == 27 || i == 28 || i == 34 || i == 36 || i == 24){
				if(Params.PARS[i] < Params.parm_low[i]){
					Params.PARS[i] = Params.parm_low[i];
				}
				if(Params.PARS[i] > Params.parm_high[i]){
					Params.PARS[i] = Params.parm_high[i];
				}
		}
}
				total_lhood=0;
				Params.pop = 0;
				crash = 0;
				error = 0.0;

				// ----------------------- loop over patch numbers -------------------------------------------- //

				//MG  changed this to index from ZERO -- check for problems
				//printf("yes I'm running\n");

redo = 0;
				while( Params.pop<FERAL_SETS && crash<1 )	{

					pop=Params.pop;
					//printf("LoopNumber = %d, first loop pop = %d\n", LoopNumber, pop);
					pop_best_lhood = -1e9;
					gsl_monte_function G = { &Hood_Pops, dim[pop], &Params };	// declares function calling Hood_Pops.h
					double xl[dim[pop]];	double xu[dim[pop]];	// need to redeclare xl and xu since the size changes
					for (jj=0;jj<=dim[pop];jj++)	{
							xl[jj]=0;
							xu[jj]=1;
							}
					// ----------- Use MISER to call function pop_lhood --------------------------------- //

					gsl_monte_miser_state *s = gsl_monte_miser_alloc(dim[pop]);

					gsl_monte_miser_integrate (&G,xl,xu,dim[pop],calls,r,s,&pop_lhood,&pop_err);

					gsl_monte_miser_free(s);
					pop_lhood2=log(pop_lhood);  //CK// Converting back to log likelihoods for MCMC
					error += pop_err;
					//printf("error = %e\n", error);

					if(crash==1){
					pop_lhood2 = -80000.0;
					//printf("crashed in lhood file\n");
						}

					if(error>0.0025){
					//printf("fucking liars\n");
					//printf("line 481: error/calls = %e\n", error/calls);
					crash = 1;
					pop_lhood2 = -20000.0;
					//getc(stdin);
						}

					if(pop_lhood2 == 0 || pop_lhood2 > 0 || isinf(pop_lhood2)== 1 || isinf(-pop_lhood2)== 1 ||isnan(pop_lhood2)== 1 ){
							crash = 1;
							///printf("crashed\n");
							pop_lhood2 = -30000.0;
							}

					total_lhood += pop_lhood2;
					Params.pop++;
					//printf("end of first loop: Params.pop = %d\n", Params.pop);

				} //CK// end of going through patches

				//printf("total_lhood = %e\n", total_lhood);
				//if(total_lhood > -170){
					if( redo > 0){	//do-over if the likelihood score doesn't contain all experiments
					//printf("first loop: redo = 1; LoopNumber = %d\n", LoopNumber);
					redo = 0;
					total_lhood=0;
					Params.pop = 0;
					crash = 0;
					error = 0.0;
					//calls = 150;

				while( Params.pop<FERAL_SETS && crash<1 )	{

					pop=Params.pop;
					pop_best_lhood = -1e9;
					gsl_monte_function G = { &Hood_Pops, dim[pop], &Params };	// declares function calling Hood_Pops.h
					double xl[dim[pop]];	double xu[dim[pop]];	// need to redeclare xl and xu since the size changes
					for (jj=0;jj<=dim[pop];jj++)	{
							xl[jj]=0;
							xu[jj]=1;
							}
					// ----------- Use MISER to call function pop_lhood --------------------------------- //

					gsl_monte_miser_state *s = gsl_monte_miser_alloc(dim[pop]);

					gsl_monte_miser_integrate (&G,xl,xu,dim[pop],calls,r,s,&pop_lhood,&pop_err);

					gsl_monte_miser_free(s);
					pop_lhood2=log(pop_lhood);  //CK// Converting back to log likelihoods for MCMC
					error += pop_err;
					if(redo > 0){
						//printf("redo = 1 again in first try\n");
						crash = 1;
					}

					if(crash==1){
					pop_lhood2 = -80000.0;
						}

					if(error>0.0025){
					crash = 1;
					pop_lhood2 = -20000.0;
					//getc(stdin);
						}

					if(pop_lhood2 == 0 || pop_lhood2 > 0 || isinf(pop_lhood2)== 1 || isinf(-pop_lhood2)== 1 ||isnan(pop_lhood2)== 1 ){
							crash = 1;
							///printf("crashed\n");
							pop_lhood2 = -30000.0;
							}

					total_lhood += pop_lhood2;
					Params.pop++;
					//printf("after double check: total_lhood = %e\n", total_lhood);
				} //CK// end of going through patches
		}//end of if statement


				LogNewPosterior = 0.0;
				LogNewPosterior=total_lhood;



	Params.PARS[4] =  Old_Params[0];//alpha1
	//Params.PARS[5] =  Old_Params[1];//beta
	//Params.PARS[6] =  Old_Params[2];//variance
	Params.PARS[7] =  Old_Params[1];//rho1
	Params.PARS[8] =  Old_Params[2];//phi
	//Params.PARS[9] =  Old_Params[3];//m1
	Params.PARS[13] = Old_Params[3];//gamma1
	Params.PARS[18] =  Old_Params[4];//delta1
	Params.PARS[34] =  Old_Params[5];//delta1

	Params.PARS[9] = 1.0;

for(i = 0; i < 37; i++){
	if(i == 4 || i == 7 || i == 8 || i == 13 || i == 18 || i == 19 || i == 21 || i == 29 || i == 30 || i == 31 || i == 32 || i == 33 || i == 35 || i == 23){
		if(Params.PARS[i] < pow(10,Params.parm_low[i])){
			Params.PARS[i] = pow(10,Params.parm_low[i]);
				}
		if(Params.PARS[i] > pow(10,Params.parm_high[i])){
			Params.PARS[i] = pow(10,Params.parm_high[i]);
				}
		}
	if(i == 9 || i == 14 || i == 27 || i == 28 || i == 34 || i == 36 || i == 24){
				if(Params.PARS[i] < Params.parm_low[i]){
					Params.PARS[i] = Params.parm_low[i];
				}
				if(Params.PARS[i] > Params.parm_high[i]){
					Params.PARS[i] = Params.parm_high[i];
				}
		}
}
			//  MG prepare for second loop over patches
				total_lhood=0;
				Params.pop = 0;
				crash = 0;
				error = 0.0;
				redo = 0;

				// ----------------------- loop over patch numbers -------------------------------------------- //
				while (Params.pop<FERAL_SETS && crash<1 )	{
					pop=Params.pop;
					//printf("second loop: pop = %d\n", pop);
					pop_best_lhood = -1e9;

					gsl_monte_function G = { &Hood_Pops, dim[pop], &Params };	// declares function calling Hood_Pops.h
					double xl[dim[pop]];	double xu[dim[pop]];	// need to redeclare xl and xu since the size changes
					for (jj=0;jj<=dim[pop];jj++)	{
						xl[jj]=0;
						xu[jj]=1;
					}

					// ----------- Use MISER to call function pop_lhood --------------------------------- //
					gsl_monte_miser_state *s = gsl_monte_miser_alloc(dim[pop]);
					gsl_monte_miser_integrate (&G,xl,xu,dim[pop],calls,r,s,&pop_lhood,&pop_err);
					gsl_monte_miser_free(s);

					pop_lhood2=log(pop_lhood);  //CK// Converting back to log likelihoods for MCMC

					error += pop_err;

					if(crash==1){
					pop_lhood2 = -80000.0;
					}

					if(error>0.0025){
					crash = 1;
					pop_lhood2 = -20000.0;
					//getc(stdin);
						}

					if(pop_lhood2 == 0 || pop_lhood2 > 0 || isinf(pop_lhood2)== 1 || isinf(-pop_lhood2)== 1 ||isnan(pop_lhood2)== 1 ){
							crash = 1;
							pop_lhood2 = -30000.0;
							}

						total_lhood += pop_lhood2;
						Params.pop++;
					//printf("end of second loop: Params.pop = %d\n", Params.pop);
				} //CK// end of going through patches

				//if(total_lhood > -170){
				//  MG prepare for second loop over patches
				if(redo > 0){
				//printf("redo = 1 for second pass; LoopNumber = %d\n", LoopNumber);
				redo = 0;
				total_lhood=0;
				Params.pop = 0;
				crash = 0;
				error = 0.0;
				//calls = 150;

				// ----------------------- loop over patch numbers -------------------------------------------- //
				while (Params.pop<FERAL_SETS && crash<1 )	{
					pop=Params.pop;
					//printf("second loop: pop = %d\n", pop);
					pop_best_lhood = -1e9;

					gsl_monte_function G = { &Hood_Pops, dim[pop], &Params };	// declares function calling Hood_Pops.h
					double xl[dim[pop]];	double xu[dim[pop]];	// need to redeclare xl and xu since the size changes
					for (jj=0;jj<=dim[pop];jj++)	{
						xl[jj]=0;
						xu[jj]=1;
					}

					// ----------- Use MISER to call function pop_lhood --------------------------------- //
					gsl_monte_miser_state *s = gsl_monte_miser_alloc(dim[pop]);
					gsl_monte_miser_integrate (&G,xl,xu,dim[pop],calls,r,s,&pop_lhood,&pop_err);
					gsl_monte_miser_free(s);

					pop_lhood2=log(pop_lhood);  //CK// Converting back to log likelihoods for MCMC

					error += pop_err;

					if(redo > 0){
						//printf("failing at end of second pass & try \n");
						crash = 1;
					}

					if(crash==1){
					pop_lhood2 = -80000.0;
					}

					if(error>0.0025){
					//printf("fucking liars\n");
					//printf("line 656: error/calls = %e\n", error/calls);
					crash = 1;
					pop_lhood2 = -20000.0;
					//getc(stdin);
						}


					if(pop_lhood2 == 0 || pop_lhood2 > 0 || isinf(pop_lhood2)== 1 || isinf(-pop_lhood2)== 1 ||isnan(pop_lhood2)== 1 ){
							crash = 1;
							pop_lhood2 = -30000.0;
							}

						total_lhood += pop_lhood2;
						Params.pop++;
					//printf("end of second check: total_lhood = %e\n", total_lhood);
				} //CK// end of going through patches
		}//end of if statement

				LogOldPosterior=0.0;
				LogOldPosterior=total_lhood;


		//---------------------Done calculating likelihood of OLD PARAM set --------------------------//

		// ----------------------Compare New Param set with Old Param Set ----------------------------------------- //


		ProbOfAcceptance=exp(LogNewPosterior+LogJumpToOld - LogOldPosterior-LogJumpToNew);    //Probability of accepting the new PC

		//printf("Reasons to stay=%.0f\n", LogNewPosterior+LogJumpToOld);
		//printf("Reasons to leave=%.0f\n",LogOldPosterior+LogJumpToNew);
		//printf("ProbOfAcceptance=%f\n", ProbOfAcceptance);
//The larger the value, the more likely to accept:  -LogNewPosterior	-LogJumpToOld 	+ LogOldPosterior	+LogJumpToNew
		//printf("loop: %d LogOldPosterior: %f\t LogNewPosterior: %f\t ProbOfAcceptance: %f\n", LoopNumber, LogOldPosterior, LogNewPosterior, ProbOfAcceptance);
		//getc(stdin);

		Params.LoopVect[Case] = Params.LoopVect[Case] + 1;

		if (ProbOfAcceptance>1 || gsl_rng_uniform_pos (r) < ProbOfAcceptance)   //MH-MCMC algorithm
		{
			LogOldPosterior=LogNewPosterior;
			//Accepted=Accepted+1;
			//printf("Accepted\n");
			Params.AcceptedVect[Case] = Params.AcceptedVect[Case] + 1;
		}

		else
		{
			for (a=0; a<NumberOfParams; a++)
			{
				PC[a]=Old_PC[a];
			}

			//printf("Rejected\n");
		}

		//printf("%f\t%f\t%f\n", LogOldPosterior, LogNewPosterior, ProbOfAcceptance);


char testbuff[128];  //Holds the file names

     FILE *fp, *fp2;   //pointers to files

        char *test[] =
       {
         "H1_30May_",
         "YourSecondFileNameHere", //in case you need two
         "YourThirdFileNameHere", //or 3
         "YourFourthFileNameHere", //or even 4
       };

     char testbuff2[128];
     char buffer0[128],buffer1[128],buffer2[128],buffer3[128];
     char bufferB[128];
     char fname;
     int pid;
     pid = getpid();  //getting the process id

     char *strFileType = ".dat";
     strcpy(buffer0, test[0]);
     sprintf(bufferB, "%d", pid);
     strcat(buffer0, bufferB);
     strcat(buffer0,strFileType);
     test[0] =  buffer0;

     strcpy(testbuff,test[0]);

		//printf("LoopNumber = %d\n", LoopNumber);
		if (LoopNumber % 100 == 0)   	//CK// output results every 5000 loops (probably not best plan but we'll see)
		{
		// ------------------------------------------ output results to file  --------------------------------------- //
			//Params.PARS[0] = Accepted;
			//Params.PARS[1] = LoopNumber;

			//printf("%f\t%f\t%f\n", LogOldPosterior, LogNewPosterior, ProbOfAcceptance);
			//printf("loopnumber = %d\n", LoopNumber);

 			fp2 = fopen(testbuff,"a"); //appending to the end of the file

  	 		 fprintf(fp2,"%f, %d, %d, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f  \n",parm_inc, calls, LoopNumber, error, Params.PARS[4], Params.PARS[5], Params.PARS[6],  Params.PARS[7],
			 Params.PARS[8],
			 Params.PARS[9],
			 Params.PARS[11],
			 Params.PARS[12],
			 Params.PARS[13],
			 Params.PARS[14],
			 Params.PARS[18],
			 Params.PARS[19],
			 Params.PARS[21],
			 Params.PARS[27],
			 Params.PARS[28],
			 Params.PARS[29],
			 Params.PARS[30],
			 Params.PARS[31],
			 Params.PARS[32],
			 Params.PARS[33],
			 Params.PARS[34],
			 Params.PARS[35],
			 Params.PARS[36],
			 Params.PARS[23],
			 Params.PARS[24],

  	 		 LogOldPosterior);   //actual printing to the file

   			 fclose(fp2);
			//molly_results = fopen(strFileName,"a");
			//output_file(&Params,molly_results,LogOldPosterior,num_adj_pars,num_runs,calls); // prints to output file (filenames.h)
			//fclose(molly_results);

			//fflush(stdout);
			//printf("Loopprint = %d\n", LoopNumber);
		// ------------------------------------------------------------------------------------------------------- //
		}


}   //end of the infinite while loop

gsl_rng_free(r);

// ------------------------------------------ output results to file  --------------------------------------- //
molly_results = fopen(strFileName,"a");
output_file(&Params,molly_results,LogOldPosterior,num_adj_pars,num_runs, calls); // prints to output file (filenames.h)
fclose(molly_results);
// ------------------------------------------------------------------------------------------------------- //
free_i3tensor(Params.DATA,0,23,0,37,0,8);
free_i3tensor(Params.EXPDATA,0,7,0,37,0,12);
free_d3tensor(Params.HostFactors,0,23,0,37,0,2);
free_d3tensor(Params.ExpFactors,0,7,0,37,0,2);
free_d3tensor(Params.FeralPara,0,23,0,37,0,4);
free_d3tensor(Params.ExpPara,0,7,0,37,0,4);
return 0;
}
