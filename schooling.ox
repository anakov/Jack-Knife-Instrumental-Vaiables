//	Jack-knife Instrumental Variables Estimation of returns to schooling
//	Program replicating the results of Angrist, Imbens, and Krueger (JAE 1999)
//	and estimating from smaller sub-samples

#import <oxstd.h>									 
#include "ols.ox"
#include "tsls.ox"
#include "liml.ox"
#include "jive1.ox"
#include "jive2.ox"
#include "setseed.ox"

main()
{
decl samples, paraMAT, DATA, data, N, n, yob, sob, qob, educ, lwage, instr, inclexog, exclexog, yob_dummy, j, i, X, b_seb;

//	loading data
	DATA    = loadmat("jive.dat");
	data    = DATA;
	samples = 1;	 // set to 100 to replicate new Table 2
	paraMAT = zeros(samples,5);
	N       = rows(DATA);
   
	for (j = 0; j < samples; j++)
	{
	setseed(j+1);
	if (samples>1) data = DATA[floor(N*ranu(50000,1))][];
	n     = rows(data);	
	yob   = data[][0];  // year of birth 
	sob   = data[][1];  // state of birth
	qob   = data[][2];  // quarter of birth
	educ  = data[][3];  // years of education
	lwage = data[][4];  // log wage

//	instruments and exogenous covariates generation
	exclexog = zeros(n,30);
	yob_dummy = zeros(n,10);
	for (i = 0; i < n; i++)
		{
		if (qob[i] != 4)
			exclexog[i][(yob[i]-30)*3 + (qob[i]-1)] = 1;  // 30 yob & qob dummy interactions
		yob_dummy[i][yob[i]-30] = 1;	// 10 year of birth dummies
		}
	instr = yob_dummy~exclexog;
 	X = educ~yob_dummy; // education and year of birth dummies are the included regressors

//	estimation
	b_seb = zeros(11,10);											   
	b_seb = ols(lwage,X)~tsls(lwage,X,instr)~liml(lwage,educ,yob_dummy,exclexog)~jive1(lwage,X,instr)~jive2(lwage,X,instr);
	paraMAT[j][] = b_seb[0][0]~b_seb[0][2]~b_seb[0][4]~b_seb[0][6]~b_seb[0][8];
	}
 	if (samples>1) 
		savemat("table 2.xls", ( (quantilec(paraMAT))'~(sqrt(varc(paraMAT)))' )' ); 
	else
		savemat("table II.xls",b_seb[0][0:1]'~b_seb[0][2:3]'~b_seb[0][4:5]'~b_seb[0][6:7]'~b_seb[0][8:9]');
}