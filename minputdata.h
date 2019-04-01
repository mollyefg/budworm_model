int inputdata(void *Paramstuff)
{

#define MAX_WEEKS 37	// larger than the number of weeks in any data set
//#define MAX_WEEKS2 100	// larger than the number of weeks in any data set
STRUCTURE* Params;
Params = (STRUCTURE*) Paramstuff;

// loads the data into a matrix and finds the number of weeks in the data set
int Sdata[MAX_WEEKS]; int Vdata[MAX_WEEKS]; int Fdata[MAX_WEEKS]; int Ddata[MAX_WEEKS]; int Cdata[MAX_WEEKS]; int D2data[MAX_WEEKS]; double DBHdata[MAX_WEEKS]; double Hdata[MAX_WEEKS]; int DAYdata[MAX_WEEKS];
int treat[MAX_WEEKS]; int Pdata[MAX_WEEKS]; int Bdata[MAX_WEEKS]; int Tdata[MAX_WEEKS]; int Ldata[MAX_WEEKS]; int stday[MAX_WEEKS]; int endday[MAX_WEEKS]; int sw[MAX_WEEKS]; int ew[MAX_WEEKS];
int Healthdata[MAX_WEEKS];
double dbh[MAX_WEEKS]; double hscore[MAX_WEEKS];
int sitenum[MAX_WEEKS]; int year[MAX_WEEKS];
//int Ddata2[MAX_WEEKS2]; double Rain[MAX_WEEKS2]; double MaxT[MAX_WEEKS2]; double MinT[MAX_WEEKS2]; double AveT[MAX_WEEKS2]; double MaxRH[MAX_WEEKS2]; double MinRH[MAX_WEEKS2]; double AveRH[MAX_WEEKS2];

int Gfdata[MAX_WEEKS];
int Afdata[MAX_WEEKS];
int Tachdata[MAX_WEEKS];
int unkdata[MAX_WEEKS];
int Gfexpdata[MAX_WEEKS];
int Afexpdata[MAX_WEEKS];
int Tachexpdata[MAX_WEEKS];
int unkexpdata[MAX_WEEKS];

int weeks;	int i; int j;
int num_weeks[FERAL_SETS+1];			// used for output to file
int num_weeks2[FERAL_SETS+1];			//CK// used for output to file
int num_weeks3[FERAL_SETS+1];			//CK// used for output to file

int total_days=0;					// the number of days summed over all data sets (for MISER)
Params->DATA = i3tensor(0,FERAL_SETS+1,0,MAX_WEEKS,0,8);		//MG  allocates memory at execution time -- allows as much memory as OS has available
Params->EXPDATA = i3tensor(0,EXP_SETS+1,0,MAX_WEEKS,0,12);  //CK// May need to check this.  Not sure what i3tensor does...
Params->HostFactors = d3tensor(0,FERAL_SETS+1,0,MAX_WEEKS,0,2);
Params->ExpFactors = d3tensor(0,EXP_SETS+1,0,MAX_WEEKS,0,2);  //CK// May need to check this.  Not sure what i3tensor does...

Params->FeralPara = d3tensor(0,FERAL_SETS+1,0,MAX_WEEKS,0,4);
Params->ExpPara = d3tensor(0,EXP_SETS+1,0,MAX_WEEKS,0,4);


int FlagF;

char *file;
char *file_name="m2data";
char *file_name2="m2expdata";  //CK// name for inputing the experimental data
char *file_type=".txt";

char *code;
char *code_name="ftp";

char numbs[6];

/*------------------------------- Data Sets ---------------------------------*/
for (j=1;j<=FERAL_SETS;j++)	{
	weeks=0;	i=0;	FlagF=0;
	FILE *ftp_data;

	sprintf(numbs,"%d",j);

	file = (char*)calloc((strlen(file_name)+strlen(file_type)+strlen(numbs)+1),sizeof(char));
	code = (char*)calloc((strlen(code_name)+strlen(numbs)+1),sizeof(char));

	strcat(file,file_name);
	strcat(file,numbs);
	strcat(file,file_type);

	strcat(code,code_name);
	strcat(code,numbs);

	ftp_data=fopen(file,"r");
	if (ftp_data==0)	{printf("file %d open error \n",j);			}
		//while (fscanf(ftp_data,"%d %d %d %d %d %lf %lf %d %d %d\n",&Sdata[i],&Vdata[i],&Fdata[i], &Ddata[i],&D2data[i],&DBHdata[i],&Hdata[i],&DAYdata[i], &sitenum[i], &year[i])!= EOF)			{
		while (fscanf(ftp_data,"%d %d %d %d %d %lf %lf %d %d %d %d %d %d %d\n",&Sdata[i],&Vdata[i],&Fdata[i], &Ddata[i],&D2data[i],&DBHdata[i],&Hdata[i],&DAYdata[i], &sitenum[i], &year[i], &Afdata[i], &Gfdata[i], &Tachdata[i], &unkdata[i])!= EOF)			{
		//printf("i = %d, %d %d %d %d %d %lf %lf %d %d %d %d %d %d %d\n",i, Sdata[i],Vdata[i],Fdata[i], Ddata[i],D2data[i],DBHdata[i],Hdata[i],DAYdata[i], sitenum[i],year[i],Afdata[i],Gfdata[i],Tachdata[i],unkdata[i]);

		Params->DATA[j][i][0]=Sdata[i];
		Params->DATA[j][i][1]=Vdata[i];
		Params->DATA[j][i][2]=Fdata[i];
		Params->DATA[j][i][3]=Ddata[i];
		Params->DATA[j][i][4]=D2data[i];
		Params->DATA[j][i][5]=DAYdata[i];
		Params->DATA[j][i][6]=sitenum[i];
		Params->DATA[j][i][7]=year[i];
		Params->HostFactors[j][i][0]=DBHdata[i];
		Params->HostFactors[j][i][1]=Hdata[i];
		//put in the multi-species parasitoid data:
		Params->FeralPara[j][i][0] = Afdata[i];
		Params->FeralPara[j][i][1] = Gfdata[i];
		Params->FeralPara[j][i][2] = Tachdata[i];
		Params->FeralPara[j][i][3] = unkdata[i];
		//printf("Params->DATA[j][i][10] = %e\n", Params->DATA[j][i][10]);
		//printf("Params->DATA[j][i][13] = %e\n", Params->DATA[j][i][13]);

		weeks++; i++;
	}

	fclose(ftp_data);

	Params->MAXT[j]=7*(weeks-1);				// number of days
	//printf("j = %d, Params->MAXT[j] = %d\n", j, Params->MAXT[j]);
	total_days += Params->MAXT[j];
	num_weeks[j]=i;
	free(file);
	free(code);

}

/*-------------------------------Experimental Data Sets ---------------------------------*/


for (j=1;j<=EXP_SETS;j++)	{
	weeks=0;	i=0;	FlagF=0;
	FILE *ftp_data;

	sprintf(numbs,"%d",j);

	file = (char*)calloc((strlen(file_name2)+strlen(file_type)+strlen(numbs)+1),sizeof(char));
	code = (char*)calloc((strlen(code_name)+strlen(numbs)+1),sizeof(char));

	strcat(file,file_name2);
	strcat(file,numbs);
	strcat(file,file_type);

	strcat(code,code_name);
	strcat(code,numbs);

	ftp_data=fopen(file,"r");
	if (ftp_data==0)	{printf("file %d open error \n",j);			}

	while (fscanf(ftp_data,"%d %d %d %d %d %d %d %d %d %d %lf %lf %d %d %d %d %d %d\n",&treat[i],&Healthdata[i],&Pdata[i],&Bdata[i],&Tdata[i],&Ldata[i],&stday[i],&endday[i],&sw[i],&ew[i],&dbh[i],&hscore[i],&sitenum[i],&year[i], &Afexpdata[i], &Gfexpdata[i], &Tachexpdata[i], &unkexpdata[i])!= EOF)			{
		Params->EXPDATA[j][i][0]=treat[i];
		Params->EXPDATA[j][i][1]=Healthdata[i]; Params->EXPDATA[j][i][2]=Pdata[i]; Params->EXPDATA[j][i][3]=Bdata[i]; Params->EXPDATA[j][i][4]=Tdata[i];
		Params->EXPDATA[j][i][5]=Ldata[i];Params->EXPDATA[j][i][6]=stday[i];Params->EXPDATA[j][i][7]=endday[i];Params->EXPDATA[j][i][8]=sw[i];Params->EXPDATA[j][i][9]=ew[i];
		Params->ExpFactors[j][i][0]=dbh[i]; Params->ExpFactors[j][i][1]=hscore[i];
		Params->EXPDATA[j][i][10]=sitenum[i];Params->EXPDATA[j][i][11]=year[i];

		//MG adding the multi-species data:
		Params->ExpPara[j][i][0] = Afexpdata[i];
		Params->ExpPara[j][i][1] = Gfexpdata[i];
		Params->ExpPara[j][i][2] = Tachexpdata[i];
		Params->ExpPara[j][i][3] = unkexpdata[i];
		//printf("j = %d, treat*exp %d\n", j, i);
		//printf("S = %d, P = %d\n", Healthdata[i], Pdata[i]);
		//printf("Af = %d, Gf = %d, Tach = %d, unk = %d\n", Afexpdata[i], Gfexpdata[i], Tachexpdata[i], unkexpdata[i]);
		//getc(stdin);
		weeks++; i++;
	}

	fclose(ftp_data);
	Params->MAXT2[j]=weeks-1;				// number of weeks
	total_days += Params->MAXT2[j];
	num_weeks2[j]=i;
	free(file);
	free(code);
}

/*----------------------------end of data sets------------------------------*/

FILE *fp_weeks;
fp_weeks=fopen("weeks.dat","w");

for (j=1;j<=FERAL_SETS;j++)	{

	fprintf(fp_weeks,"%d\t",num_weeks[j]);
}
fclose(fp_weeks);
//free(file);
//free(code);
return 0;
}

