// ---------------------------------------------------------------------------------------------------------------- //
double bound(int i,int j)				// bounds on parameters for parhood line search
{
double low;
double high;


if (i==1)	{	low = -5.00;		high = 1.0001;	}	//
//parameter 2 is declared in the parhood file for some reason?
/*
if (i==3)	{	low = -3.500;		high = 1.50001;	}	//MG  eta1 - fitted parameter for d-d in the experimental treatments
else if (i==4)	{	low = -3.0001;			high = 0.70001;	}		//k       //MG   alpha
else if (i==5)	{	low = -3.0001;	high = 0.7001;	}		//MG// beta  ///??
//temporarily restrict sigma to be very small
else if (i==6)	{	low = -6.0001;		high = 0.5001;}		//MG  variance parameter!
else if (i==7)	{	low = -3.0001;			high = 1.0001;}		//MG  rho: the half-saturation constant  //CK//  beta = the size of accumulated rainfall window
else if (i==8)	{	low = -3.0001;			high = 0.0;}		//MG  phi: adjusts attack rate for partial experimental treatment
else if (i==9)	{	low = 1;			high = 21;	}		//MG  m, number of exposed classes // make sure the parm_inc and step size only offers whole numbers!

//parameter 10 holds Q in the diffeq's

else if (i==11)	{	low = -4.5;		high =0.71;	}		//MG   mu1		shape the quality function
else if (i==12)	{	low = -4.5;		high = 0.71;	}		//MG   mu2		shape the quality function

else if (i==13)	{	low = -3.5;			high = 1.5001;	}		//MG gamma1  density-death term
else if (i==14)	{	low = 1;			high = 41;	}		//MG gamma2	 density-death term

//parameters 15 and 16 hold alpha and phi in the odes
else if (i==17)	{	low = -3.5;		high = 1.5001;	}		//MG   eta2		// term for density-dependence in the experimental treatments ONLY (model Q_expDD)

else if (i==18)	{	low = -4.0001;		high = 0.3001;	}		//MG// delta!  rate of progression through Em

//add in: alpha2 and alpha3
else if (i==19)	{	low = -3.0001;			high = 0.3001;	}		//k       //MG   alpha2 [attack rate of parasitoid 2]
//parameter 20 holds alpha2 in odes
else if (i==21)	{	low = -3.0001;			high = 0.3001;	}		//k       //MG   alpha3 [attack rate of parasitoid 3]
//parameter 22 holds alpha3 in odes
*/
if (i==3)	{	low = -3.0001;		high = 0.70001;	}	//MG  eta1 - fitted parameter for d-d in the experimental treatments
else if (i==4)	{	low = -3.0001;	high = 4.0001;	}		//k       //MG   alpha
else if (i==5)	{	low = -3.0001;	high = -0.7001;	}		//MG// beta  ///??
//temporarily restrict sigma to be very small
else if (i==6)	{	low = -4.0001;		high = 0.010;}		//MG  variance parameter!
//else if (i==6)	{	low = -1.50001;		high = 0.5001;}		//MG  variance parameter!

else if (i==7)	{	low = -0.0001;			high = 4.0001;}		//MG  rho: the half-saturation constant  //CK//  beta = the size of accumulated rainfall window
else if (i==8)	{	low = -2.0001;			high = 0.0;}		//MG  phi: adjusts attack rate for partial experimental treatment
else if (i==9)	{	low = 1;			high = 21;	}		//MG  m, number of exposed classes // make sure the parm_inc and step size only offers whole numbers!

//parameter 10 holds Q in the diffeq's

else if (i==11)	{	low = -3.5;		high =0.61;	}		//MG   mu1		shape the quality function
else if (i==12)	{	low = -2.6;		high = -0.22;	}		//MG   mu2		shape the quality function

//MG need to find the previous values for this term before A-D models
//else if (i==13)	{	low = -5.0;			high = 0.1001;	}		//MG gamma1  density-death term
//else if (i==13)	{	low = 0.001;			high = 0.7001;	}		//MG gamma1  density-death term
else if (i==13)	{	low = -5.01;			high = 0.0;	}		//MG gamma1  density-death term

else if (i==14)	{	low = 1;			high = 21;	}		//MG gamma2	 density-death term

//parameters 15 and 16 hold alpha and phi in the odes
else if (i==17)	{	low = -3.5;		high = 1.5001;	}		//MG   eta2		// term for density-dependence in the experimental treatments ONLY (model Q_expDD)

else if (i==18)	{	low = -4.0001;		high = 0.0001;	}		//MG// delta!  rate of progression through Em

//add in: alpha2 and alpha3
else if (i==19)	{	low = -3.0001;	high = 4.0001;	}		//k       //MG   alpha2 [attack rate of parasitoid 2]
//parameter 20 holds alpha2 in odes
else if (i==21)	{	low = -3.0001;	high = 4.0001; }		//k       //MG   alpha3 [attack rate of parasitoid 3]
//parameter 22 holds alpha3 in odes

//add in: phi2 and phi3
//MG  formerly for phi2 and 3 values and holders -- now use 23 and 24 for gamma7 and gamma8
else if (i==23)	{	low = -5.01;			high = 0.0;	}		//MG gamma7  density-death term
else if (i==24)	{	low = 1;			high = 21;	}		//MG gamma8  density-death term

//parameters 25 is the holder for eta1 and eta2

//add in: m2 and m3
else if (i==27)	{	low = 1;				high = 21;	}		//MG  m2, number of exposed classes // make sure the parm_inc and step size only offers whole numbers!
else if (i==28)	{	low = 1;				high = 21;	}		//MG  m3, number of exposed classes // make sure the parm_inc and step size only offers whole numbers!

else if (i==29)	{	low = -0.0001;			high = 4.0001;}	//MG  rho2: the half-saturation constant  //CK//  beta = the size of accumulated rainfall window
else if (i==30)	{	low = -0.0001;			high = 4.0001;}		//MG  rho3: the half-saturation constant  //CK//  beta = the size of accumulated rainfall window

else if (i==31)	{	low = -4.0001;		high = 0.0001;	}		//MG// delta2!  rate of progression through Em2
else if (i==32)	{	low = -4.0001;		high = 0.0001;	}		//MG// delta3!  rate of progression through Em3
else if (i==33)	{	low = -5.01;			high = 0.0;	}		//MG gamma3  density-death term

else if (i==34)	{	low = 1;			high = 21;	}		//MG gamma4	 density-death term

else if (i==35)	{	low = -5.01;			high = 0.0;}		//MG gamma5  density-death term
else if (i==36)	{	low = 1;			high = 21;	}		//MG gamma6	 density-death term
//else if (i==37)	{	low = -2.0;			high = -1.0;	}		//MG gamma7  density-death term
//else if (i==38)	{	low = 1;			high = 51;	}		//MG gamma8	 density-death term

else			{	low = 1;			high = 1;	}

if		(j==1)	return low;
else if (j==2)	return high;
else { printf("PROBLEM WITH BOUNDS ON PARAMETERS!!\n");	return 0;	}
}

// ---------------------------------------------------------------------------------------------------------------- //
void parm_range_inc(void *Paramstuff,double parm_inc,double host_inc,double initR_inc,int num_adj_pars)   // MG  all the parameter ranges. All of them!
{
STRUCTURE* Params;
Params = (STRUCTURE*) Paramstuff;
int i;
Params->parm_step[1]=1;
Params->parm_step[7]=1;
if (num_adj_pars>0)	{
	for (i=1;i<=num_adj_pars;i++)	{
		Params->parm_low[i]  = bound(i,1);	Params->parm_high[i] = bound(i,2);
		if (i==2||i==3||i==4||i==5||i==6||i==7||i==8||i==9||i==10||i==11||i==12||i==13||i==14||i==15||i==16||i==17||i==18||i==19||i==20||i==21||i==22||i==23||i==24||i==25||i==26||i==27||i==28||i==29||i==30||i==31||i==32||i==33||i==34||i==35||i==36)	{
			Params->parm_step[i]=(Params->parm_high[i]-Params->parm_low[i])/parm_inc;
			//printf("i = %d, steps size = %e\n", i, Params->parm_step[i]);
		}
	}
}
}
