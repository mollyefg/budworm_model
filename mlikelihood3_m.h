extern int verbose;
extern int exrealz;

extern int crash;
extern int redo;

//extern double lhood_check;
//double lpass;

double Like(double *RandNumsPass, size_t dim, void *Paramstuff,  double *ap, double *np, double *gp, double *op, double *ae, double *ne, double *ge, double *oe)

{
		STRUCTURE* Params;
		Params = (STRUCTURE*) Paramstuff;
		int pop = Params->pop;
		double liktemp, lik, likesite, likesitetemp;
		double liktemp2, likesitetemp2;
		int siteyear;
		int MAXT4=(50);
		int w;
		double sim_results[9][4];
		double exp_results[12][3][4];
		int numrows = 9;
		double lamb0[numrows]; double lamb1[numrows]; double lamb2[numrows]; double lamb3[numrows];

		double lamb_S[3]; double lamb_E1[3];  double lamb_E2[3]; double lamb_E3[3]; //3 treatments and 3 parasitoids to consider
		int j, start;
		double p, q, x, y;
		double pe[12][3], qe[12][3], xe[12][3], ye[12][3];
		int experiment, expcount, ticker;
		double pesum, qesum, xesum, yesum;

		int divisor;
		divisor = pow(10,3);
		//printf("divisor = %d\n", divisor);
		double liketest, liket;
		int i;

		int count;
		count = 0;

		//printf("before ddevf\n");
		DDEVF(Params,RandNumsPass,dim,pop,MAXT4,sim_results,exp_results);  //
		//printf("like 34\n");
		liktemp=0.0;
		lik=0.0;
		likesitetemp = 0.0;
		liktemp2 = 0.0;
		likesitetemp2 = 0.0;
		int sw;

///// Not sure if this next part lines up properly; should it be sim_results[w+1] ??
/*if(pop==4||pop==9||pop==18||pop==19){
	 sw = 2;
 	 }
 	 //else{
    if(pop==0||pop==1||pop==2||pop==3||pop==5||pop==6||pop==7||pop==8||pop==10||pop==11||pop==12||pop==13||pop==14||pop==15||pop==16||pop==17||pop==20||pop==21){
   	 sw = 3;
 	 }*/
 	 sw =3;
		 for(w=sw; w<numrows; w++){

	    lamb0[w] = sim_results[w][0];
	    //printf("lamb0 = %e \n", lamb0[w]);
		lamb1[w] = sim_results[w][1];
		lamb2[w] = sim_results[w][2];
		lamb3[w] = sim_results[w][3];
		//printf("w = %d, lamb0 = %e, lamb1 = %e, lamb2 = %e, lamb3 = %e\n", w, lamb0[w], lamb1[w], lamb2[w], lamb3[w]);

	  		}  //close "for" loop


p = 0.0;
q = 0.0;
x = 0.0;
y = 0.0;

	for(j = sw; j < 8; j++){		// calculates for weeks 3 through 7
			if(op[j]>0.0 & lamb3[j]>0.0){
							y = op[j]*log(lamb3[j])-lamb3[j];
							/*if(verbose==1){
														FILE *fp;
														fp = fopen("F2 op.dat","a");
														fprintf(fp, "%d, %d, %e, %e\n", pop, j, op[j], y);
														fclose(fp);
							}*/
						}
			if(op[j]>0.0 & lamb3[j] <= 0.0){
							//redo = 1;
							//printf("redo line 96\n");
							lamb3[j] = op[j]/divisor;
							//printf("divisor=%d\n", divisor);
							y = op[j]*log(lamb3[j])-lamb3[j];
							//gp[j] = gp[j] + op[j];
							count = count++;
							/*if(verbose==1){
							FILE *fp;
							fp = fopen("F2_penalties op.dat","a");
							fprintf(fp, "%d, %d, %e, %e\n", pop, j, op[j], y);
							fclose(fp);
							}*/
						}
			if(op[j]<=0.0){
							y = -lamb3[j];
							//y = 0.0;
							//gp[j] = gp[j] + op[j]
							/*if(verbose==1){
														FILE *fp;
														fp = fopen("F2 op.dat","a");
														fprintf(fp, "%d, %d, %e, %e\n", pop, j, op[j], y);
														fclose(fp);
							}*/
								}
			//printf("op[j] = %e, lamb3[j] = %e, divisor = %d\n", op[j], lamb3[j], divisor);
			//getc(stdin);
			////Glypta parasitoids
			if(gp[j]>0.0 & lamb2[j]>0.0){
							x = gp[j]*log(lamb2[j])-lamb2[j];
							/*if(verbose==1){
							FILE *fp;
							fp = fopen("F2 gp.dat","a");
							fprintf(fp, "%d, %d, %e, %e\n", pop, j, gp[j], x);
							fclose(fp);
							}*/
							}
			if(gp[j]>0.0 & lamb2[j]<=0.0){
							//redo = 1;
							//printf("redo line 111\n");
							lamb2[j] = gp[j]/divisor;
							x = gp[j]*log(lamb2[j])-lamb2[j];
							//ap[j] = ap[j] + gp[j];
							count = count++;
						/*	if(verbose==1){
							FILE *fp;
							fp = fopen("F2_penalties gp.dat","a");
							fprintf(fp, "%d, %d, %e, %e\n", pop, j, gp[j], x);
							fclose(fp);
							}*/
							}
			if(gp[j]<=0.0){
							x = -lamb2[j];
							/*if(verbose==1){
							FILE *fp;
							fp = fopen("F2 gp.dat","a");
							fprintf(fp, "%d, %d, %e, %e\n", pop, j, gp[j], x);
							fclose(fp);
							}	*/						//x = 0.0;
							//ap[j] = ap[j] + gp[j];
							}


			/////Apanteles parasitoids
			if(ap[j]>0.0 & lamb1[j]> 0.0){
							q = ap[j]*log(lamb1[j])-lamb1[j];
							/*if(verbose==1){
							FILE *fp;
							fp = fopen("F2 ap.dat","a");
							fprintf(fp, "%d, %d, %e, %e\n", pop, j, ap[j], q);
							fclose(fp);
							}*/
						}
			if(ap[j]>0.0 & lamb1[j] <=0.0){
							//redo = 1;
							//printf("redo line 128\n");
							lamb1[j] = ap[j]/divisor;
							q = ap[j]*log(lamb1[j])-lamb1[j];

							count = count++;
							/*if(verbose==1){
							FILE *fp;
							fp = fopen("F2_penalties ap.dat","a");
							fprintf(fp, "%d, %d, %e, %e\n", pop, j, ap[j], q);
							fclose(fp);
							}*/

							}
			if(ap[j]<=0.0){
							//q = 0.0;
							q = -lamb1[j];
							/*	if(verbose==1){
							FILE *fp;
							fp = fopen("F2 ap.dat","a");
							fprintf(fp, "%d, %d, %e, %e\n", pop, j, ap[j], q);
							fclose(fp);
							}*/
							}

			///// Hosts
			if(np[j]>0 & lamb0[j]> 0.0){
							p = np[j]*log(lamb0[j])-lamb0[j];

						/*	if(verbose==1){
							FILE *fp;
							fp = fopen("F2 np.dat","a");
							fprintf(fp, "%d, %d, %e, %e\n", pop, j, np[j], p);
							fclose(fp);
							}*/
							}
			if(np[j]>0 & lamb0[j] <= 0.0){
							//redo = 1;
							//printf("redo line 138\n");
							lamb0[j] = np[j]/divisor;
							p = np[j]*log(lamb0[j])-lamb0[j];

							count = count++;
							/*if(verbose==1){
							FILE *fp;
							fp = fopen("F2_penalties np.dat","a");
							fprintf(fp, "%d, %d, %e, %e\n", pop, j, np[j], p);
							fclose(fp);
							}*/
							}
		 	if(np[j]<=0.0){
							p = -lamb0[j];
							/*if(verbose==1){
							FILE *fp;
							fp = fopen("F2 np.dat","a");
							fprintf(fp, "%d, %d, %e, %e\n", pop, j, np[j], p);
							fclose(fp);
							}*/
							//p = 0.0;
							}
//printf("pop = %d, count = %d\n",  pop, count);


			liktemp = p + q + x + y;  //MG
			likesitetemp+= liktemp;  //MG
									/*if(verbose==1){
										FILE *fp;
										fp = fopen("best_realdat_30jan.dat","a");
													fprintf(fp, "%d, %d, %e, %e, %e, %e\n", pop, j, np[j], ap[j], gp[j], op[j]);

										  	fclose(fp);
										}*/
									/*if(verbose==1){
										FILE *fp;
										fp = fopen("E1_lambdas_4Feb.dat","a");
											  			fprintf(fp, "%d, %d, %e, %e, %e, %e\n", pop, j, lamb0[j], lamb1[j], lamb2[j], lamb3[j]);

										  	fclose(fp);
										}*/
										/*if(verbose==1){
											FILE *fp;
											fp = fopen("rho2_pq_30jan.dat","a");
															fprintf(fp, "%d, %d, %e, %e, %e, %e\n", pop, j, p, q, x, y);


											  	fclose(fp);
										}*/

	} // close "for" loop
if(verbose==1){
	FILE *fp;
	fp = fopen("out_L3_lowHlowLlowQ.dat","a");
	//fprintf(fp, "population, lambda1, lambda2, rcc\n");
	for(i = 0; i <8; i ++){
		  			fprintf(fp, "%d, %e, %e, %e, %e\n", pop, lamb0[i], lamb1[i], lamb2[i], lamb3[i]);
		  			//printf("pop = %d, i = %d, lamb0 = %e, lamb1 = %e, lamb2 = %e, lamb3 = %e\n", pop, i, lamb0[i], lamb1[i], lamb2[i], lamb3[i]);
		  			//getc(stdin);
	  			}
	  	fclose(fp);
	}
			//start trying to include the experimental data
				if( pop==3 || pop==4 || pop==9 || pop==15 || pop==16){
				pesum = 0.0;
				qesum = 0.0;
				xesum = 0.0;
				yesum = 0.0;

					if(pop == 3){
						expcount = 5;
						}//pop = 3
					if(pop == 4){
						expcount = 5;
						}
					if(pop == 9){
						expcount = 5;
						}

					if(pop == 15){
						expcount = 12;
						}
					if(pop == 16){
						expcount = 11;
						}
						for(i=0; i<12; i++){
							for(j = 0; j<3; j++){
							pe[i][j] = 0.0;
							qe[i][j] = 0.0;
							xe[i][j] = 0.0;
							ye[i][j] = 0.0;
							}
						}
				ticker = 0;
				for(experiment = 0; experiment < expcount; experiment++){	//loop over experiments for that population
					for(w=0;w<3;w++){	//loop over treatments for that population (
						lamb_S[w] = exp_results[experiment][w][0];
						lamb_E1[w] = exp_results[experiment][w][1];
						lamb_E2[w] = exp_results[experiment][w][2];
						lamb_E3[w] = exp_results[experiment][w][3];
						//printf("pop = %d ,w = %d, lambs = %e, lambe1 = %e, lambe2 = %e, lambe3 = %e\n", pop, w, lamb_S[w], lamb_E1[w], lamb_E2[w], lamb_E3[w]);
						//getc(stdin);
						///catching all bad values:
						if(isnan(lamb_S[w]) == 1 || isnan(lamb_E1[w]) == 1){
							//printf("likelihood redo: pop = %d\n", Params->pop);
							redo=1;
						}



							///// hosts

								if(ne[ticker]>0.0 & lamb_S[w] > 0.0){
									pe[experiment][w] = ne[ticker]*log(lamb_S[w]) - lamb_S[w];
										} //close if ne

								if(ne[ticker]>0.0 & lamb_S[w] <= 0.0){
									//redo = 1;
									//printf("redo 231\n");
									lamb_S[w] = ne[ticker]/divisor;
									pe[experiment][w] = ne[ticker]*log(lamb_S[w]) - lamb_S[w];
										}

								if(ne[ticker]<=0.0){
									//pe[experiment][w] = 0.0;
									pe[experiment][w] = -lamb_S[w];
										}


								/// Apanteles parasitoids

								if(ae[ticker]>0.0 & lamb_E1[w] >0.0){
									qe[experiment][w] = ae[ticker]*log(lamb_E1[w]) - lamb_E1[w];
										}

								if(ae[ticker]>0.0 & lamb_E1[w] <=0.0){
									//redo = 1;
									//printf("redo 249\n");
									lamb_E1[w] = ae[ticker]/divisor;
									qe[experiment][w] = ae[ticker]*log(lamb_E1[w]) - lamb_E1[w];
									}

								if(ae[ticker] <=0.0){
									//qe[experiment][w] = 0.0;
									qe[experiment][w] = -lamb_E1[w];
									}

								//// Glypta parasitoids

								if(ge[ticker]>0.0 & lamb_E2[w] > 0.0){
									xe[experiment][w] = ge[ticker]*log(lamb_E2[w]) - lamb_E2[w];
									}
								if(ge[ticker]>0.0 & lamb_E2[w] <= 0.0){
									//redo = 1;
									//printf("redo 265\n");
									lamb_E2[w] = ge[ticker]/divisor;
									xe[experiment][w] = ge[ticker]*log(lamb_E2[w]) - lamb_E2[w];
									}
								if(ge[ticker] <= 0.0){
									xe[experiment][w] = -lamb_E2[w];
									}


								///// Other parasitoids

								if(oe[ticker]>0.0 & lamb_E3[w] > 0.0){
									ye[experiment][w] = oe[ticker]*log(lamb_E3[w]) - lamb_E3[w];
										}
								if(oe[ticker]>0.0 & lamb_E3[w] <= 0.0){
									//redo = 1;
									//printf("redo 279\n");
									lamb_E3[w] = oe[ticker]/divisor;

									ye[experiment][w] = oe[ticker]*log(lamb_E3[w]) - lamb_E3[w];
									}

								if(oe[ticker]<=0.0){
									ye[experiment][w] = -lamb_E3[w];
									}
								/*if(verbose==1){
										FILE *fp;
										fp = fopen("best_expcomponents.dat","a");
				fprintf(fp, "%d, %d, %d, %e, %e, %e, %e\n", pop, experiment, w, pe[experiment][w], qe[experiment][w],xe[experiment][w],ye[experiment][w]);
										fclose(fp);
										}*/

							//printf("pop = %d, experiment = %d, w = %d, pe[experiment][w] = %e, qe[experiment][w] = %e\n", pop, experiment, w, pe[experiment][w], qe[experiment][w]);
							//printf("experiment = %d, treatment = %d, pe = %e, qe = %e, xe = %e, ye = %e\n", experiment, w, pesum, qesum, xesum, yesum);
								pesum += pe[experiment][w];
								qesum += qe[experiment][w];
								xesum += xe[experiment][w];
								yesum += ye[experiment][w];
								ticker++;

															/*if(verbose==1){
																	FILE *fp;
																	fp = fopen("null_expoutput.dat","a");
																	fprintf(fp, "%d, %d, %d, %e, %e, %e, %e\n", pop, experiment, w, lamb_S[w], lamb_E1[w],lamb_E2[w],lamb_E3[w]);
																	fclose(fp);
																		}*/

							}//close loop over treatments



						//catching worst bad values:
						if(pe[experiment][0]==0&&pe[experiment][1]==0&&pe[experiment][2]==0&&qe[experiment][0]==0&&qe[experiment][1]==0&&qe[experiment][2]==0){
							crash=1;
						}
						//printf("next loop: crash = %d\n", crash);
					}//close loop over experiments
					//printf("line 250, crash = %d\n", crash);
			liktemp2 = pesum + qesum + xesum + yesum;
			likesitetemp2 += liktemp2;

			} // close for loop over pops with experiments
			//check for problems with the parameter set just for experiments:
			if(pop == 3 || pop == 4 || pop == 9 || pop == 15 || pop == 16 ){
				if(likesitetemp2 == 0){
					//likesitetemp2 = -502.0;
					crash = 1;
					//printf("bad experiment value 258\n");
					//printf("crash = %d\n", crash);
				}
			}
			//printf("likelihood line 246 \n");
			if(likesitetemp == 0.0){
							//printf("pop = %d, bad obs value, likesitetemp = 0.0 \n", pop);
							//likesitetemp = -501.0;
							crash = 1;
			}



		//printf("pop = %d, likesitetemp = %e, likesitetemp2 = %e\n", pop, likesitetemp, likesitetemp2);

		likesite = likesitetemp + likesitetemp2;
			 /* if(verbose==1){
							FILE *fp;
							fp = fopen("null_components.dat","a");
							//fprintf(fp, "population, lambda1, lambda2, rcc\n");
							//for(i = 0; i <8; i ++){
								  			fprintf(fp, "%d, %e, %e, %e\n", pop, likesitetemp, likesitetemp2, likesite);
							  			//}
							  			//printf("printed in l\n");
							  	fclose(fp);
							}*/
		//printf("likesite = %e\n", likesite);
		//getc(stdin);
	return(likesite);
}
//// 10, 100, up 2, 4, 6, 8, 10

double Likelihood(double *RandNumsPass,size_t dim,void *Paramstuff)
{
	STRUCTURE* Params;
	Params = (STRUCTURE*) Paramstuff;
	int pop = Params->pop;
	double likesum=0.0;
	int w;

	double N[10], A[10], G[10], O[10];

	for(w = 0; w<8; w++){

		double ntop = (Params->DATA[pop+1][w][0]);
		double nbottom = (Params->DATA[pop+1][w][4]);

		double atop =(Params->FeralPara[pop+1][w][0]);
		double abottom = (Params->DATA[pop+1][w][4]);

		double gtop =(Params->FeralPara[pop+1][w][1]);
		double gbottom = (Params->DATA[pop+1][w][4]);

		double otop =(Params->FeralPara[pop+1][w][2])+(Params->FeralPara[pop+1][w][3]);
		double obottom = (Params->DATA[pop+1][w][4]);
		//printf("atop = %e, gtop = %e, otop = %e, bottom = %e \n", atop, gtop, otop, nbottom);
		//getc(stdin);

		if(ntop<0){
		ntop = 0.0;
		}
		if(atop<0){
		atop = 0.0;
		}
		if(gtop<0){
		gtop = 0.0;
		}
		if(otop<0){
		otop = 0.0;
		}


		if(nbottom>0){
		N[w]=ntop/nbottom;
		}else{
		N[w] = 0.0;

		}
		if(abottom>0){
		A[w]=atop/abottom;
		}else{
		A[w] = 0.0;

		}
		if(gbottom>0){
		G[w]=gtop/gbottom;
		}else{
		G[w] = 0.0;

		}
		if(obottom>0){
		O[w]=otop/obottom;
		}else{
		O[w] = 0.0;

		}

		ntop = 0.0;
		nbottom = 0.0;
		atop = 0.0;
		abottom = 0.0;
		gtop = 0.0;
		gbottom = 0.0;
		otop = 0.0;
		obottom = 0.0;
		}

		double NE[12][3], AE[12][3], GE[12][3], OE[12][3];
		if(pop == 3 || pop == 4 || pop == 9 || pop == 15 || pop == 16 ){
		int i;
		int j; // populations/experiments
		int k; // treatment
		int expset, expcount, experiment;

				for(i = 0; i < 12; i++){
					for(j = 0; j < 3; j++){
						NE[i][j] = 0.0;
						AE[i][j] = 0.0;
						GE[i][j] = 0.0;
						OE[i][j] = 0.0;
									}
							}
					if(pop == 3){
						expset = 1;
						expcount = 5;
						}//pop = 3
					if(pop == 4){
						expset = 2;
						expcount = 5;
						}
					if(pop == 9){
						expset = 3;
						expcount = 5;
						}

					if(pop == 15){
						expset = 5;
						expcount = 12;
						}
					if(pop == 16){
						expset = 6;
						expcount = 11;
						}
			//int expcounters[6] = {0, 5, 10, 15, 27, 39};
			for(experiment = 0; experiment < expcount; experiment ++){
				for(k = 0; k < 3; k++){

				int checker = experiment + k + 2*experiment;
					double netop = (Params->EXPDATA[expset][checker][1]);
					double nebottom = (Params->EXPDATA[expset][checker][3]);

					double aetop =(Params->ExpPara[expset][checker][0]);
					double aebottom = (Params->EXPDATA[expset][checker][3]);

					double getop = (Params->ExpPara[expset][checker][1]);
					double gebottom = (Params->EXPDATA[expset][checker][3]);

					double oetop =(Params->ExpPara[expset][checker][2])+(Params->ExpPara[expset][checker][3]);
					double oebottom = (Params->EXPDATA[expset][checker][3]);

				if(netop<0){
				netop = 0.0;
					}
				if(aetop<0){
				aetop = 0.0;
					}
				if(getop<0){
				getop = 0.0;
					}
				if(oetop<0){
				oetop = 0.0;
					}
				if(nebottom>0){
				NE[experiment][k]=netop/nebottom;
					}else{
				NE[experiment][k] = 0.0;
					}
				if(aebottom>0){
				AE[experiment][k]=aetop/aebottom;
					}else{
				AE[experiment][k] = 0.0;
					}

				if(gebottom>0){
				GE[experiment][k]=getop/gebottom;
					}else{
				GE[experiment][k] = 0.0;
					}
				if(oebottom>0){
				OE[experiment][k]=oetop/oebottom;
					}else{
				OE[experiment][k] = 0.0;
					}
				/*
				if(nebottom>0){
				NE[experiment][k]=netop/nebottom;
					}else{
				NE[experiment][k] = -100.0;
					}
				if(aebottom>0){
				AE[experiment][k]=aetop/aebottom;
					}else{
				AE[experiment][k] = -100.0;
					}

				if(gebottom>0){
				GE[experiment][k]=getop/gebottom;
					}else{
				GE[experiment][k] = -100.0;
					}
				if(oebottom>0){
				OE[experiment][k]=oetop/oebottom;
					}else{
				OE[experiment][k] = -100.0;
					}
					*/

				netop = 0.0;
				nebottom = 0.0;
				aetop = 0.0;
				aebottom = 0.0;
				getop = 0.0;
				gebottom = 0.0;
				oetop = 0.0;
				oebottom = 0.0;
				}// treatment (k) loop

				}//for loop
		}//if statement
/*
/////  Putting the "real" data in a file:
int i;
char string1[]="realdata820_";
char filename1[30];
sprintf(filename1,"%s%d",string1,pop);
//printf("%s\n",filename1);
	FILE *f = fopen(filename1, "w");
	if (f == NULL)
	{
  	  printf("Error opening file!\n");
   	 exit(1);
	}
	fprintf(f, "population, np, ap, fix\n");
		for(i = 0; i <8; i ++){
			fprintf(f, "%d, %e, %e\n", pop, N[i], A[i]);
			}
	fclose(f);
/////  done filing real data
*/
	likesum=Like(RandNumsPass, dim, Paramstuff, A, N, G, O, *AE, *NE, *GE, *OE);  //likelihood for one site, all weeks
	return(likesum);
}

double Hood_Pops(double *RandNumsPass,size_t dim,void *Paramstuff)
{

	double likelihoodtemp;
	likelihoodtemp = 0.0;
	STRUCTURE* Params;
	Params = (STRUCTURE*) Paramstuff;
	//printf("before Likelihood\n");
	likelihoodtemp+=Likelihood(RandNumsPass, dim, Paramstuff);
	//printf("after Likelihood\n");
	double lhood;
	lhood = likelihoodtemp;
	//printf("like file: lhood = %e, exp(lhood) = %e\n", lhood, exp(lhood));
	lhood = exp(lhood);

	if( lhood==1 ){
		//lhood = 250000.0;
		crash = 1;
		//printf("bad final lhood 465\n");
		}
	if( lhood==0 ){
		//lhood = 3500000.0;
		crash = 1;
		//printf("what now, fucker?\n");
		}
		//printf("end of lhood: crash = %d\n", crash);
		return(lhood);
}
