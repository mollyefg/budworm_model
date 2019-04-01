
void GetString(int j,int i,char* strOut, unsigned int strSize)
{
// ---------------------------------------- Name for Output Files ------------------------------------------- //
char *Prefix=" ";
char *max_lhood_prefix="G1_noexp_";
char *mcmc_prefix="this is a mistake";
char *profile_prefix="intrinsic_";

char *FileName	= " ";
char *Path		= " ";
char *Name		= " ";
char *Type		= ".dat";						// type of output file
char profile[33];
// ---------------------------------------- create time suffixes -------------------------------------------- //
struct tm *ptr;									// struct defined in time.h;  used to hold the time and date
time_t lt;										// type defined in time.h; for storing the calendar time
lt = time(NULL);
ptr = localtime(&lt);							// converts calendar time to local time
char Date[30];									// string to hold date and time
strftime(Date, 30, "%m_%d_%H_%M_%S", ptr);		// adds to Date month, day, year, hour, minute, second from ptr
// ---------------------------------- create paths and allocate memory ------------------------------------- //
if	(j==0)	{										// MCMC
	Path="/home/mgallagher/mcmc_results/";
	Prefix=mcmc_prefix;
	if		(i==1)	{	Name = "parms_";											}
	else if	(i==2)	{	Name = "pc_";												}
	else if (i==3)	{	Name = "acc_";												}
	else if (i==4)	{	Name = "lastpars_";											}
	else			{	printf("bad i value in filenames!!!\n");	getc(stdin);	}
	FileName = (char*)calloc((strlen(Path)+strlen(Prefix)+strlen(Name)+strlen(Type)+strlen(Date)+1), sizeof(char));
}
else if (j==1)	{									// Line Search MLE
	Prefix=max_lhood_prefix;
	if		(i==0)	{	Path="/home/mgallagher/moutput/";			Name="max_lhood";			}
	else if (i==4)	{	Path="/home/mgallagher/moutput/";			Name="L_";				}
	else			{	printf("bad i value in filenames!!!\n");		getc(stdin);	}
	FileName = (char*)calloc((strlen(Path)+strlen(Prefix)+strlen(Name) + strlen(Type)+1), sizeof(char));
}
else if (j>1)	{									// Profile Likelihoods
	Path="/home/mgallagher/profile_data/";
	Prefix=profile_prefix;
	sprintf(profile,"%d",j);
	Name="pro_";
	FileName = (char*)calloc((strlen(Path)+strlen(Prefix)+strlen(Name) + strlen(Type)+strlen(profile) +1), sizeof(char));
}
// ---------------------------------- add strings onto the file name ------------------------------------------ //
strcat(FileName,Path);
strcat(FileName,Prefix);
strcat(FileName,Name);

if		(j==0||(j==1 && i==4))	strcat(FileName,Date);		// add Date to output file name		(MCMC & max lhood)
else if (j>1)					strcat(FileName,profile);	// add profile number to file name	(profile likelihood)

strcat(FileName,Type);

if (j==0)			strncpy(strOut,FileName,strSize);
else if (j>=1)		strncpy(strOut,FileName,strSize);
free(FileName);
return;
}

void output_file(void *Paramstuff,FILE *fp_results,double best_post_hood, int num_adj_pars,int num_runs, double parm_inc)
//void output_file(double best_post_hood)

{
STRUCTURE* Params;
Params = (STRUCTURE*) Paramstuff;
int i,pop;


fprintf(fp_results, "%d,\t", num_runs);
fprintf(fp_results, "%e,\t", parm_inc);
fprintf(fp_results, "%d,\t", 1);	//MG model number

fprintf(fp_results, "%e,\t", Params->PARS[3]);  // ERROR  // Eta 1: model L only
fprintf(fp_results, "%e,\t", Params->PARS[4]);  // alpha1
fprintf(fp_results, "%e,\t", Params->PARS[5]);  // beta
fprintf(fp_results, "%e,\t", Params->PARS[6]);  // variance
fprintf(fp_results, "%e,\t", Params->PARS[7]);  // rho1
fprintf(fp_results, "%e,\t", Params->PARS[8]);  // phi
fprintf(fp_results, "%e,\t", Params->PARS[9]);  // m1
fprintf(fp_results, "%e,\t", Params->PARS[11]);  // mu1
fprintf(fp_results, "%e,\t", Params->PARS[12]);  // mu2
fprintf(fp_results, "%e,\t", Params->PARS[13]);  // gamma1
fprintf(fp_results, "%e,\t", Params->PARS[14]);  // gamma2
fprintf(fp_results, "%e,\t", Params->PARS[17]);  // // Eta 2: model L only


fprintf(fp_results, "%e,\t", Params->PARS[18]);  // delta1

fprintf(fp_results, "%e,\t", Params->PARS[19]);  // alpha2
fprintf(fp_results, "%e,\t", Params->PARS[21]);  // alpha3
fprintf(fp_results, "%e,\t", Params->PARS[27]);  // m2
fprintf(fp_results, "%e,\t", Params->PARS[28]);  // m3
fprintf(fp_results, "%e,\t", Params->PARS[29]);  // rho2
fprintf(fp_results, "%e,\t", Params->PARS[30]);  // rho3
fprintf(fp_results, "%e,\t", Params->PARS[31]);  // delta2
fprintf(fp_results, "%e,\t", Params->PARS[32]);  // delta3

fprintf(fp_results, "%e,\t", Params->PARS[33]);  // gamma3
fprintf(fp_results, "%e,\t", Params->PARS[34]);  // gamma4
fprintf(fp_results, "%e,\t", Params->PARS[35]);  // gamma5
fprintf(fp_results, "%e,\t", Params->PARS[36]);  // gamma6
fprintf(fp_results, "%e,\t", Params->PARS[23]);  // gamma7
fprintf(fp_results, "%e,\t", Params->PARS[24]);  // gamma8

fprintf(fp_results, "%e,\t", best_post_hood);


fprintf(fp_results,"\n");

}
