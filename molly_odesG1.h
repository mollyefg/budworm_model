extern int crash;

// --------------------------------Begin ODE system of model T3--------------------------------------------//


int fast_odes(double t, const double y[], double dydt[],void *Paramstuff)
{
STRUCTURE* Params;
Params = (STRUCTURE*) Paramstuff;

int i, j;


double alpha1 = Params->PARS[15];   			//MG   the attack rate of the parasitoids
double alpha2 = Params->PARS[15];   			//MG   the attack rate of the parasitoids
double alpha3 = Params->PARS[15];   			//MG   the attack rate of the parasitoids

double phi = Params->PARS[16];

double rho1 = Params->PARS[7];
double rho2 = Params->PARS[7];
double rho3 = Params->PARS[7];

int m1	= Params->PARS[9];  //MG// the number of exposed classes
int m2	= Params->PARS[9];  //MG// the number of exposed classes
int m3	= Params->PARS[9];  //MG// the number of exposed classes

double gamma3 = Params->PARS[33];
double gamma4 = Params->PARS[34];

double delta1 = Params->PARS[18];
double delta2 = Params->PARS[18];
double delta3 = Params->PARS[18];

double summation = 0.0;

double Qmax = 43.2562;
double Q = (Params->PARS[10])/Qmax;

int P1sub, P2sub, P3sub;
P1sub = m1+m2+m3+1;
P2sub = m1+m2+m3+2;
P3sub = m1+m2+m3+3;

// ------------------------------------------ ODEs -------------------------------------------- //
if		(y[0]<.000001)	dydt[0]=0;
else

summation = 0.0;
for(j = 0; j<m1+m2+m3+1; j++){
	summation += y[j];
	}
//getc(stdin);

dydt[0] = (-alpha1*phi*y[0]*y[P1sub])/(1+rho1*y[P1sub]) -(alpha2*phi*y[0]*y[P2sub])/(1 + rho2*y[P2sub]) -(alpha3*phi*y[0]*y[P3sub])/(1 + rho3*y[P3sub]) - (gamma3*Q*summation*pow(y[0],gamma4));
dydt[1] = (alpha1*phi*y[0]*y[P1sub])/(1+rho1*y[P1sub]) - m1*delta1*y[1]  - (gamma3*Q*summation*pow(y[1],gamma4));

for(i=2; i <= m1; i++){
	dydt[i] = m1*delta1*y[i-1] - m1*delta1*y[i] - (gamma3*Q*summation*pow(y[i],gamma4));
	if(isnan(dydt[i])== 1 ){
		crash = 1;
		}
	}

dydt[m1+1] = (alpha2*phi*y[0]*y[P2sub])/(1 + rho2*y[P2sub]) - m2*delta2*y[m1+1]  - (gamma3*Q*summation*pow(y[m1+1],gamma4));

for(i=m1 + 2; i <= m1 + m2; i++){
	dydt[i] = m2*delta2*y[i - 1] - m2*delta2*y[i] - (gamma3*Q*summation*pow(y[i],gamma4));
	if(isnan(dydt[i])== 1 ){
			crash = 1;
		}
	}

dydt[m1+m2+1] = (alpha3*phi*y[0]*y[P3sub])/(1 + rho3*y[P3sub]) -m3*delta3*y[m1+m2+1]  - (gamma3*Q*summation*(pow(y[m1+m2+1],gamma4)));
for(i=m1 + m2 + 2; i <= m1 + m2 + m3; i++){
	dydt[i] = m3*delta3*y[i-1] - m3*delta3*y[i]  - (gamma3*Q*summation*pow(y[i],gamma4));
	if(isnan(dydt[i])== 1 ){
			crash = 1;
		}
	}

dydt[m1+m2+m3+1] = m1*delta1*y[m1];		//Parasite 1
dydt[m1+m2+m3+2] = m2*delta2*y[m1 + m2];		//Parasite 2
dydt[m1+m2+m3+3] = m3*delta3*y[m1 + m2 + m3]; //Parasite 3

//mortality due to direct d-d and quality:
dydt[m1+m2+m3+4] = (gamma3*Q*summation*pow(y[0],gamma4)) +
				   (gamma3*Q*summation*pow(y[1],gamma4)) +
				   (gamma3*Q*summation*pow(y[m1+1],gamma4)) +
				   (gamma3*Q*summation*pow(y[3],gamma4));

return GSL_SUCCESS;
}

// ------------------------------------------  ODE Solver  ----------------------------------------------- //
double ODE_Solver(double t_ode,double t_end,void *Paramstuff,double *y_ode)
{
int i;
int status_ode;
double h_init=1.0e-5;

STRUCTURE* Params;
Params = (STRUCTURE*) Paramstuff;

//int DIM = Params->PARS[9]+2;		//MG  updated
int DIM = Params->PARS[9] + Params->PARS[9] + Params->PARS[9] + 5;

const gsl_odeiv_step_type *solver_ode	= gsl_odeiv_step_rkf45; // Runge-Kutta Felberg (4, 5)

// returns pointer to a newly allocated instance of a stepping function of type 'solver_ode' for a system of DIM dimensions //
gsl_odeiv_step *step_ode	= gsl_odeiv_step_alloc(solver_ode, DIM);

gsl_odeiv_control *tol_ode	= gsl_odeiv_control_standard_new(1.0e-10, 1.0e-5, 1.0, 0.2);
gsl_odeiv_evolve *evol_ode	= gsl_odeiv_evolve_alloc(DIM);
gsl_odeiv_system sys_ode;
sys_ode.function  = fast_odes;
sys_ode.dimension = (size_t)(DIM);
sys_ode.params	  = Params;
sys_ode.jacobian = NULL;


// ----------------------------------- Integrate Over Time ------------------------------------ //
while (t_ode<t_end)	{
	status_ode = gsl_odeiv_evolve_apply(evol_ode, tol_ode, step_ode, &sys_ode, &t_ode, t_end, &h_init, y_ode);

	for (i=0;i<DIM;i++)	{
		if (y_ode[i]<0)		{
		y_ode[i]=0;
		}
	}
}

// -------------------------------------- Clear Memory ----------------------------------------- //
gsl_odeiv_evolve_free(evol_ode);
gsl_odeiv_control_free(tol_ode);
gsl_odeiv_step_free(step_ode);
return (t_end);
}