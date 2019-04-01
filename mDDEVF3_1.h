extern int crash;

double DDEVF(void *Paramstuff,double *RandNumsPass,size_t dim,int pop,int maxy_t,double (*sim_results)[4], double (*exp_results)[3][4])
{
//printf("DDEVF!\n");
STRUCTURE* Params;
Params = (STRUCTURE*) Paramstuff;
int DIM = Params->PARS[9] + Params->PARS[9] + Params->PARS[9] + 5;  //MG ?
int m1 = Params->PARS[9];
int m2 = m1;
int m3 = m1;

double t=h;		double t_next=h;	double t_0=h;
double y_ode[DIM];
double S,I1,E1,I2,E2,I3,E3;
int day;
int week;
double alpha1;
double phi;
double alpha2;
double alpha3;

int ticker;
double RandNums1[dim];
double RandNums2[dim];
double RandNums3[dim];

double sd=sqrt(Params->PARS[6]);
//printf("sd = %e\n", sd);
double Q;

double exp_QD_mortality;
double obs_QD_mortality;

int i,j, k;
for(i = 0; i<9; i++){
	for(j=0;j<4;j++){
		sim_results[i][j]=0.0;
	}
}

for(i=0;i<DIM;i++){
	y_ode[i] = 0.0;
	}

double initialS[22] = {0.3092784,0.32,0.5019608,0.064,0.296,0.656,0.908,1.51,4.48,9.767442,0.04,0.3,1.326316,1.306818,1.283784,2.557143,1.648649,0.12,0.9344262,0.7764706,1.931818,1.577778};

double initialP1[22] = {0.115463917525773,0.018,0.0591836734693878,0.016,0.08,0.056,0.096,0.08,0.21,0.348837209302326,0.015,0.07,0.505263157894737,0.693181818181818,0.689189189189189,0.628571428571429,0.662162162162162,0.05,0.671641791044776,0.235294117647059,0.772727272727273,1.26666666666667};
double initialP2[22] = {0.0144329896907217,0.002,0.0224489795918367,0.00001,0.00001,0.004,0.004,0.02,0.1,0.162790697674419,0.00001,0.005,0.0736842105263158,0.0454545454545455,0.108108108108108,0.214285714285714,0.0405405405405405,0.00001,0.149253731343284,0.0941176470588235,0.159090909090909,0.244444444444444};
double initialP3[22] = {0.00001,0.004,0.0163265306122449,0.00001,0.00001,0.008,0.00001,0.015,0.00001,0.00001,0.00001,0.00001,0.0105263157894737,0.00001,0.00001,0.00001,0.00001,0.00001,0.00001,0.00001,0.00001,0.00001};


// -------------------------- integrate until next stoppage event ---------------------------------- //

    double maxt = dim;

    //setting up RandNumsPass thing
alpha1 = Params->PARS[4];
alpha2 = Params->PARS[4];
alpha3 = Params->PARS[4];

phi = Params->PARS[8];


/*		for(i=0;i<dim;i++){
						RandNums1[i] = alpha1*exp(gsl_cdf_gaussian_Pinv(RandNumsPass[i],sd));  //Greg's version
						RandNums2[i] = alpha2*exp(gsl_cdf_gaussian_Pinv(RandNumsPass[i],sd));  //Greg's version
						RandNums3[i] = alpha3*exp(gsl_cdf_gaussian_Pinv(RandNumsPass[i],sd));  //Greg's version
						//RandNums1[i] = alpha1;
						//RandNums2[i] = alpha2;
						//RandNums3[i] = alpha3;
							//FILE *fp;
							//fp = fopen("modelC_alpha1values.dat","a");
							//fprintf(fp, "%e\n", RandNums1[i]);
							//fclose(fp);
						//printf("sigma = %e, alphas = %e %e %e\n", Params->PARS[6], RandNums1[i], RandNums2[i], RandNums3[i]);
						//getc(stdin);


				}*/
    t = 0;

double initS = initialS[pop];
double initE1 = initialP1[pop];
double initE2 = initialP2[pop];
double initE3 = initialP3[pop];

y_ode[0]=initS;
y_ode[1]=initE1;
y_ode[m1+1] = initE2;
y_ode[m1 + m2 + 1] = initE3;

double sumE1;
double sumE2;
double sumE3;

   // week = 1;
	if(pop==4||pop==9||pop==18||pop==19){
	 week = 2;
 	 }
 	 //else{
    if(pop==0||pop==1||pop==2||pop==3||pop==5||pop==6||pop==7||pop==8||pop==10||pop==11||pop==12||pop==13||pop==14||pop==15||pop==16||pop==17||pop==20||pop==21){
   	 week = 3;
 	 }

 	ticker = 0;
    while(t<maxt){
		//alpha1 = RandNums1[ticker];			//MG   the attack rate of the parasitoids
		//alpha2 = RandNums2[ticker];			//MG   the attack rate of the parasitoids
		//alpha3 = RandNums3[ticker];			//MG   the attack rate of the parasitoids
		Params->PARS[15] = alpha1;		//try this dummy holder
		Params->PARS[20] = alpha1;		//try this dummy holder
		Params->PARS[22] = alpha1;		//try this dummy holder
		phi = 1.0;
		Params->PARS[16] = phi;


		//y_ode[0]=S;
		//y_ode[1]=E1;
		//y_ode[m1+1] = E2;
		//y_ode[m1 + m2 + 1] = E3;
		S = 0.0;
		sumE1 = 0.0;
		sumE2 = 0.0;
		sumE3 = 0.0;

	    t_next = t + 1;
		Q = Params->HostFactors[pop+1][week-1][0];
			if(Q==0){
					Q = 0.00005;
		}
		Params->PARS[10] = Q;

		//printf("before: t = %e,  y_ode[0] = %e, y_ode[1] = %e, y_ode[m1+1] = %e, y_ode[m1+m2+1] = %e\n", t, y_ode[0], y_ode[1], y_ode[m1+1], y_ode[m1+m2+1]);
		t=ODE_Solver(t,t_next,Params,y_ode);		//MG// y_ode[0] is the first equation, 1 is the second, etc.
		//printf("after the solver: t = %e, y_ode[0] = %e, y_ode[1] = %e, y_ode[m1+1] = %e, y_ode[m1+m2+1] = %e\n", t, y_ode[0], y_ode[1], y_ode[m1+1], y_ode[m1+m2+1]);

				S = y_ode[0];

				for(i=1;i<m1+1;i++){
						sumE1 += y_ode[i];
						//E1 += y_ode[i];
					}
				for(i=m1+1;i<m1+m2+1;i++){
						sumE2 += y_ode[i];
						//E2 += y_ode[i];
					}
				for(i=m1+m2+1;i<m1+m2+m3+1;i++){
						sumE3 += y_ode[i];
						//E3 += y_ode[i];
					}
				day = t;
						//printf("after: S = %e, E1 = %e, E2 = %e, E3 = %e\n", S, E1 , E2, E3);
						//getc(stdin);

		if (day%7==0)	{
					week=week+1;
					sim_results[week][0]=S;
					sim_results[week][1]=sumE1;
					sim_results[week][2]=sumE2;
					sim_results[week][3]=sumE3;
					//sim_results[week][1]=E1;
					//sim_results[week][2]=E2;
					//sim_results[week][3]=E3;
					//printf("Params: alpha %e beta %e sigma %e rho %e phi %e m1 %e delta1 %e alpha2 %e m2 %e rho2 %e delta2 %e \n",
											//Params->PARS[4], Params->PARS[5],Params->PARS[6],Params->PARS[7],Params->PARS[8],Params->PARS[9],Params->PARS[18],Params->PARS[19],Params->PARS[27],Params->PARS[29],Params->PARS[31]);
					//printf("week = %d, S = %e, sumE1 = %e, sumE2 = %e, sumE3 = %e \n", week, S, sumE1, sumE2, sumE3);
					//getc(stdin);
					obs_QD_mortality = y_ode[4];
					FILE *fp;
					fp = fopen("G1_mortality_16sept.dat","a");
								fprintf(fp, "pop = %d, initS = %e, initE = %e, week = %d, survivors = %e, para = %e, QD = %e\n", pop, initS, initE1+initE2+initE3, week, S, sim_results[week][1] + sim_results[week][2] + sim_results[week][3], obs_QD_mortality);
					fclose(fp);
				}

			ticker++;

		} //MG close the while loop

 //now let's try adding in experimental data ... !
if(pop==3||pop==4||pop==9||pop==15||pop==16){

	int ticker2 = 0;
	double t2 = 0;
	double t2_next;
	double t2_0;
	double maxt2 = 0;
	int week2;

	int day2;
	int expcount;
	int exp;
 	t2_next=h; t2_0=h;

	int expset;

	//size_t runlength[50] = {7, 7, 14, 15, 36, 7, 7, 14, 15, 36, 7, 7, 14, 14, 35, 6, 8, 6, 9, 6, 14, 14, 15, 15, 35, 29, 21, 6, 8, 6, 7, 7, 13, 13, 13, 14, 33, 27, 19,6, 8, 6, 9, 6, 14, 14, 15, 35, 29, 21};
	int startweek[38] = {2, 3, 3, 5, 2,
						 2, 3, 3, 5, 2,
						 2, 3, 3, 5, 2,

						 1, 2, 3, 4, 5, 1, 2, 3, 4, 1, 2, 3,
						 1, 2, 3, 4, 5, 1, 2, 3, 1, 2, 3};

	size_t startday[38] = {8, 15, 15, 29, 8,
						   8, 15, 15, 29, 8,
						   8, 15, 15, 29, 8,
						   0, 6, 13, 19, 26, 0, 6, 13, 19, 0, 6, 14,
						   0, 6, 14, 20, 29, 0, 6, 14, 0, 6, 14};
	size_t endday[38] = {15, 22, 29, 44, 44,
				         15, 22, 29, 44, 44,
				         15, 22, 29, 43, 43,
				         6, 14, 19, 26, 33, 13, 19, 26, 33, 33, 33, 33,
				         6, 14, 20, 29, 35, 14, 20, 29, 35, 35, 35};

	//double initial_exp_S[48] = {0.068, 0.032, 0.064, 0.056, 0.08, 0.068, 0.024, 0.0, 0.056,	0.296,	0.272,	0.24,	0.136,	0.132,	0.096,	0.0, 1.483333333,	9.76744186,	2.158333333,	1.224,	0.346666667,	0.304,	0.272,	0.0, 0.0,	0.105,	0.3,	0.41,	0.395,	0.255,	0.16,	0.015, 0.01,	0.884615385,	2.557142857,	0.739837398,	0.245,	0.31,	0.105,	0.016666667, 0.015,	0.3875,	1.648648649,	0.18134715,	0.35,	0.345,	0.185,	0.066666667};
	//starting with the site density at the maximum (week 2 or 3) for susceptible and parasitized
	double initial_exp_S[38] = {0.064, 0.064, 0.064, 0.08, 0.064,
								0.296, 0.272, 0.272, 0.136, 0.296,
								9.767, 2.1583, 2.1583, 0.3467, 9.767,
								2.557143, 2.557143, 2.557143, 0.7398374, 0.245, 2.557143, 2.557143, 2.557143, 0.7398374, 2.557143, 2.557143, 2.557143,
								1.648649, 1.648649, 1.648649, 0.1813472, 0.35, 1.648649, 1.648649, 1.648649, 1.648649, 1.648649, 1.648649};

	double initial_exp_I1[38] = {0.016, 0.016, 0.016, 0.00001, 0.016,
								0.08, 0.028, 0.028, 0.020, 0.08,
								0.3488, 0.29167, 0.29167, 0.11333, 0.3488,
								0.62857, 0.62857, 0.62857, 0.5447154, 0.21, 0.62857, 0.62857, 0.62857, 0.5447154, 0.62857, 0.62857, 0.62857,
								0.6621622, 0.6621622, 0.6621622, 0.4818653, 0.060,0.6621622, 0.6621622, 0.6621622, 0.6621622, 0.6621622, 0.6621622};

	double initial_exp_I2[38] = { 0.00001,  0.00001,  0.00001,  0.00001, 0.00001,
								 0.00001, 0.0004, 0.0004, 0.00001, 0.00001,
							   	 0.16279, 0.03333,0.033333, 0.023333, 0.16279,
								 0.21428,  0.21428, 0.21428, 0.08130, 0.020,0.21428,  0.21428, 0.21428, 0.08130, 0.21428,  0.21428, 0.21428,
								 0.04375, 0.04375, 0.04054, 0.04145, 0.01, 0.04375, 0.04375, 0.04054, 0.04375, 0.04375, 0.04054};

	double initial_exp_I3[38] = {0.00001, 0.00001, 0.00001, 0.00001, 0.00001,
								  0.00001,  0.00001, 0.00001, 0.00001, 0.00001,
								 0.00001, 0.083333, 0.083333, 0.0233333, 0.00001,
								 0.00001,  0.00001,  0.00001, 0.02439, 0.005, 0.00001,  0.00001,  0.00001, 0.02439,  0.00001,  0.00001,  0.00001,
								 0.005181347, 0.005181347, 0.005181347, 0.005181347, 0.005,0.005181347, 0.005181347, 0.005181347, 0.005181347, 0.005181347, 0.005181347};


	//starting with the actual site density at that week
	/*double initial_exp_S[38] = {0.032, 0.064, 0.064, 0.08, 0.032,
								0.056, 0.296, 0.296, 0.136, 0.056,
								9.767, 2.1583, 2.1583, 0.3467, 9.767,
								0.045, 0.8846154, 2.557143, 0.7398374, 0.245, 0.045, 0.8846154, 2.557143, 0.7398374, 0.045, 0.8846154, 2.557143,
								0.015, 0.3875, 1.648649, 0.1813472, 0.35,0.015, 0.3875, 1.648649, 0.015, 0.3875, 1.648649};

	double initial_exp_I1[38] = {0.0, 0.016, 0.016, 0.0, 0.0,
								0.08, 0.028, 0.028, 0.020, 0.08,
								0.3488, 0.29167, 0.29167, 0.11333, 0.3488,
								0.0, 0.353846, 0.62857, 0.5447154, 0.21, 0.0, 0.353846, 0.62857, 0.5447154,0.0, 0.353846, 0.62857,
								0.0, 0.08125, 0.6621622, 0.4818653, 0.060,0.0, 0.08125, 0.6621622, 0.0, 0.08125, 0.6621622};

	double initial_exp_I2[38] = {0.0, 0.0, 0.0, 0.0, 0.0,
								 0.0, 0.0004, 0.0004, 0.0, 0.0,
							   	 0.16279, 0.03333,0.033333, 0.023333, 0.16279,
								 0.0, 0.061538, 0.2142, 0.08130, 0.020, 0.0, 0.061538, 0.2142, 0.08130, 0.0, 0.061538, 0.2142,
								 0.005, 0.04375, 0.04054, 0.04145, 0.01, 0.005, 0.04375, 0.04054, 0.005, 0.04375, 0.04054};

	double initial_exp_I3[38] = {0.0, 0.0, 0.0, 0.0, 0.0,
								 0.0, 0.0, 0.0, 0.0, 0.0,
								 0.0, 0.083333, 0.083333, 0.233333, 0.0,
								 0.0, 0.0, 0.0, 0.02439, 0.005, 0.0, 0.0, 0.0, 0.02439, 0.0, 0.0, 0.0,
								 0.0, 0.0, 0.0, 0.005181347, 0.005,0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	*/
	//double initial_exp_I1[48] = {0.0,0.0,0.0160,0.0080,0.0,	0.0,0.0,0.0,0.0040,	0.080,0.0280,0.0160,0.020,0.00,0.0040,0.0,0.10,0.3488372,0.2916667,	0.4120,	0.1133333,0.020,	0.0,	0.0,	0.0050,	0.0750,	0.070,	0.1450,	0.1050,	0.070,	0.0,0.0,	0.0,	0.3538462,	0.6285714,	0.5447154,	0.21,	0.1850,	0.0,	0.0,	0.0,	0.0812500,	0.6621622,	0.4818653,	0.060,	0.045,	0.0,	0.0};
	//double initial_exp_I2[48] = {0.00000,	0.00000,	0.00000,	0.00000,	0.00000,	0.00000,	0.00000,	0.00000,	0.00000,	0.00000,	0.00400,	0.00800,	0.00000,	0.01600,	0.04400,	0.00000,	0.06111,	0.16279,	0.03333,	0.04400,	0.02333,	0.01200,	0.04400,	0.00000,	0.00000,	0.01000,	0.00500,	0.00000,	0.00000,	0.01000,	0.04000,	0.05000,	0.00000,	0.06154,	0.21429,	0.08130,	0.02000,	0.02500,0.05500,	0.02500,	0.00500,	0.04375,	0.04054,	0.04145,	0.01000,	0.01000,	0.05000,	0.02500};
	//double initial_exp_I3[48] = {0.00000,	0.00000,	0.00000,	0.00000,	0.00000,	0.00000,	0.00400,	0.00000,	0.00000,	0.00000,	0.00000,	0.00000,	0.00000,	0.00400,	0.00400,	0.00000,	0.00000,	0.00000,	0.08333,	0.00400,	0.02333,	0.01600,	0.01600,	0.00000,	0.00000,	0.00000,	0.00000,	0.01500,	0.00500,	0.02500,	0.02500,	0.01667,	0.00000,	0.00000,	0.00000,	0.02439,	0.00500,	0.01500,	0.00500,	0.00000,	0.00000,	0.00000,	0.00000,	0.00518,	0.00500,	0.03000,	0.02500,	0.00000};

	//try this quick fix; runlength has to be at least 7 days for the ODES to work
	//for(i = 0; i < 50; i ++){
	//	if(runlength[i] < 7){
	//		runlength[i] = runlength[i] + 1;
	//	}
	//}
for(i = 0; i < 38; i ++){
	if((endday[i] - startday[i])<7){
		endday[i] = endday[i] + 1;
	}
}

	for(i=1;i<m1+m2+m3+1;i++){	//MG ? maybe
		y_ode[i] = 0.0;
			}

	if(pop == 3){
		expset = 1;
		expcount = 5;
		}

	if(pop == 4){
		expset = 2;
		expcount = 5;
		}

	if(pop == 9){
		expset = 3;
		expcount = 5;
		}

	if(pop == 15){
		expset = 4;
		expcount = 12;
		}

	if(pop == 16){
		expset = 5;
		expcount = 11;
		}

int expcounters[6] = {0, 5, 10, 15, 27};
for(i = 0; i<12; i++){				// maximum of 12 experiments per population
	for(j = 0; j<3; j++){			// treatments: 0 = exposed; 1 = partial; 2 = covered
		for(k=0;k<4;k++){			// classes: 0 = susceptible; 1 = Af; 2 = Gf; 3 = Others;
			exp_results[i][j][k]=0.0;
			}
		}
	}
for(exp = 0; exp < expcount; exp++){

int treat;

for(treat = 1; treat < 4; treat ++){

		int checker = expcounters[expset - 1] + exp;
		//printf("checker = %d\n", checker);
		t2 = startday[checker];
		ticker2 = 0.0;
		week2 = startweek[checker];
		maxt2 = endday[checker];


		/*double Sexp = initexpS;
		double Iexp1 = initexpI1;
		double Eexp1 = initexpI1;
		double sumEexp1 = 0.0;
		double Iexp2 = initexpI2;
		double Eexp2 = initexpI2;
		double sumEexp2 = 0.0;
		double Iexp3 = initexpI3;
		double Eexp3 = initexpI3;
		double sumEexp3 = 0.0;
		*/

				double initexpS = initial_exp_S[checker]; //substract one to account for weeks indexed from
				double initexpI1 = initial_exp_I1[checker];
				double initexpI2 = initial_exp_I2[checker];
				double initexpI3 = initial_exp_I3[checker];

		y_ode[0] = initexpS;
		y_ode[1] = initexpI1;
		y_ode[m1+1] = initexpI2;
		y_ode[m1+m2+1] = initexpI3;

		double Sexp;
		double sumEexp1;
		double sumEexp2;
		double sumEexp3;


		if(treat == 1){
				//alpha1 = RandNums1[ticker2];
				//alpha2 = RandNums2[ticker2];
				//alpha3 = RandNums3[ticker2];
				Params->PARS[15] = alpha1;
				Params->PARS[20] = alpha1;
				Params->PARS[22] = alpha1;

				phi = 1.0;
				Params->PARS[16] = phi;
				//printf("alpha1 = %e\n", alpha1);
					}
			if(treat == 2){
				//alpha1 = RandNums1[ticker2];
				//alpha2 = RandNums2[ticker2];
				//alpha3 = RandNums3[ticker2];
				Params->PARS[15] = alpha1;
				Params->PARS[20] = alpha1;
				Params->PARS[22] = alpha1;

				phi = Params->PARS[8];		//actually using the real value for phi here
				Params->PARS[16] = phi;
				//printf("phi = %e\n", phi);

					}
			if(treat == 3){
				//alpha1 = RandNums1[ticker2];
				//alpha2 = RandNums2[ticker2];
				//alpha3 = RandNums3[ticker2];
				Params->PARS[15] = alpha1;
				Params->PARS[20] = alpha1;
				Params->PARS[22] = alpha1;

				phi = 0.0;// phi2 = 0.0; phi3 = 0.0;
				//Params->PARS[16] = phi;

					}
//printf("y_ode[0] = %e\n", y_ode[0]);
    while(t2<maxt2){

	    t2_next = t2 + 1;
		Q = Params->ExpFactors[expset][exp*3][0];
		if(Q==0){
					Q = 0.00005;
			}
		Params->PARS[10] = Q;
		sumEexp1 = 0.0;
		sumEexp2 = 0.0;
		sumEexp3 = 0.0;

		t2=ODE_Solver(t2,t2_next,Params,y_ode);		//MG// y_ode[0] is the first equation, 1 is the second, etc.

				Sexp = y_ode[0];
				//printf("y_ode[0] = %e\n", y_ode[0]);

				for(i=1;i<m1+1;i++){
						sumEexp1 += y_ode[i];
						//Eexp1 += y_ode[i];
					}
				for(i=m1+1;i<m1+m2+1;i++){
							sumEexp2 += y_ode[i];
							//Eexp2 += y_ode[i];
					}
				for(i=m1+m2+1;i<m1+m2+m3+1;i++){
							sumEexp3 += y_ode[i];
							//Eexp3 += y_ode[i];
					}

				day2 = t2;
					if (day2%7==0)	{
						week2=week2+1;

						exp_results[exp][treat-1][0]=Sexp;
						exp_results[exp][treat-1][1]=sumEexp1;
						exp_results[exp][treat-1][2]=sumEexp2;
						exp_results[exp][treat-1][3]=sumEexp3;
					//	printf("week2 = %d, Sexp = %e, P1 = %e, P2 = %e, P3 = %e\n", week2, Sexp, sumEexp1, sumEexp2, sumEexp3);
						//getc(stdin);

						exp_QD_mortality = y_ode[4];

					FILE *fp;
					fp = fopen("G1_exp_mortality.dat","a");
								fprintf(fp, "pop = %d, initS = %e, initE = %e, exp = %d, treat = %d, survivors = %e, para = %e, QD = %e\n", pop, initexpS, initexpI1+initexpI2+initexpI3, exp, treat-1, Sexp, exp_results[exp][treat-1][1] +
								exp_results[exp][treat-1][2] +
								exp_results[exp][treat-1][3], exp_QD_mortality);
					fclose(fp);
						}
				ticker2++;
				} //MG close the while loop that runs the ODES
			}  //MG  close the "for" loop over treatments
		} //MG close the "for" loop over experiments
	}//MG close the "if" statement selecting populations
	return 0;
}	//MG close double DDEVF