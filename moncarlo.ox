// Jack-Knife Instrumental Variables Estimation (Angrist, Imbens, and Krueger, JAE 1999)
// Program for Replication and Extension of the Monte Carlo Study
// (C) Anton Nakov, 2010

 #import <oxstd.h>
 #include "ols.ox"
 #include "tsls.ox"
 #include "liml.ox"
 #include "jive1.ox"
 #include "jive2.ox"
 #include "stats.ox"
 #include "setseed.ox"

 main()
 {
   decl T, N, L, y, x, Z, e_nu, var, rho, A, r, eo, b_seb, table, i;
   T      = 100;            // sample size: set to 25 to replicate Table 4
   N      = 5000;           // number of simulations
   b_seb  = zeros(N,10);    // matrix with point estimates and standard errors
   table  = zeros(29,7);    // matrix with quantiles, MAE and coverage rate
   rho    = 0;              // error autocorrelation: set to 0.95 to replicate Table 3

// matrices used to build autocorrelated errors and to implement GLS  
   A      = lower(toeplitz(rho^(range(0,T-1))')); 
   r      = rho^(range(1,T))';
   eo     = rann(1,1)*sqrt(1/(1-rho^2));
   L      = unit(T);
   L[1][0]= -rho;
   L      = lower(toeplitz(L[][0]));
   L[0][0]= sqrt(1-rho^2);

   setseed(1);
										 
// Model 1: L = 2, K = 3
   var = <0.25,0.2;0.2,0.25>;
   Z = ones(T,1)~rann(T,2); 
   for (i = 0; i < N; i++)
   {
     e_nu = ((r*eo + A*rann(T,1))~rann(T,1))*(choleski(var))';
     x = 0.3*Z[][1] + e_nu[][1];
     y = x + e_nu[][0];
	 x = L*x;
	 y = L*y;
     b_seb[i][] = ols(y,x)~tsls(y,x,Z)~liml(y,x,<>,Z)~jive1(y,x,Z)~jive2(y,x,Z);
   }
   table[0:4][] = stats(b_seb, N);
		
// Model 2: L = 2, K = 21
   var = <0.25,0.2;0.2,0.25>;
   Z = ones(T,1)~rann(T,20);
   for (i = 0; i < N; i++)
   {
     e_nu = ((r*eo + A*rann(T,1))~rann(T,1))*(choleski(var))';
     x = 0.3*Z[][1] + e_nu[][1];
     y = x + e_nu[][0];
 	 x = L*x;
	 y = L*y;
     b_seb[i][] = ols(y,x)~tsls(y,x,Z)~liml(y,x,<>,Z)~jive1(y,x,Z)~jive2(y,x,Z);
   }
   table[6:10][] = stats(b_seb, N);
	
// Model 3: L = 2, K = 21
   var = <1.0,0.8;0.8,1.0>;
// var = 10*var;              // comment in to replicate new Table 1
   Z = ones(T,1)~rann(T,20);
   for (i = 0; i < N; i++)
   {
     e_nu = ((r*eo + A*rann(T,1))~rann(T,1))*(choleski(var))';
     x = 0.3*(Z[][1]) + 0.3*sumsqrr(Z[][2:20]) + e_nu[][1].*(sumsqrr(Z[][2:20]))/19;
     y = x + e_nu[][0];
 	 x = L*x;
	 y = L*y;
     b_seb[i][] = ols(y,x)~tsls(y,x,Z)~liml(y,x,<>,Z)~jive1(y,x,Z)~jive2(y,x,Z);
   }
   table[12:16][] = stats(b_seb, N);
	
// Model 4: L = 2, K = 21
   var = <0.25,0.2;0.2,0.25>;
   Z = ones(T,1)~rann(T,20);
   for (i = 0; i < N; i++)
   {
     e_nu = ((r*eo + A*rann(T,1))~rann(T,1))*(choleski(var))';
     x = e_nu[][1];
     y = x + e_nu[][0];
	 x = L*x;
	 y = L*y;
     b_seb[i][] = ols(y,x)~tsls(y,x,Z)~liml(y,x,<>,Z)~jive1(y,x,Z)~jive2(y,x,Z);
   }
   table[18:22][] = stats(b_seb, N);
		
// Model 5: L = 2, K = 21
   var = <0.25,0.2;0.2,0.25>;
   Z = ones(T,1)~rann(T,20);
   for (i = 0; i < N; i++)
   {
     e_nu = ((r*eo + A*rann(T,1))~rann(T,1))*(choleski(var))';
     x = 0.3*Z[][1] + e_nu[][1];
     y = x + 0.2*Z[][2] + e_nu[][0];
 	 x = L*x;
	 y = L*y;
     b_seb[i][] = ols(y,x)~tsls(y,x,Z)~liml(y,x,<>,Z)~jive1(y,x,Z)~jive2(y,x,Z);
   }
   table[24:28][] = stats(b_seb, N);

 savemat("table I.xls", table);
 }