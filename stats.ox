//  Calculates the quantiles, Median Absolute Error, and 95% confidence interval coverage
stats (const b_seb, const N)

{
	decl b, quants, MAE, CovRate, seb;

	b		= zeros(N,5);
	quants	= zeros(5,5);
	MAE		= zeros(5,1);
	CovRate	= zeros(5,1);

	b = (b_seb[][0]~b_seb[][2]~b_seb[][4]~b_seb[][6]~b_seb[][8])';
	seb = (b_seb[][1]~b_seb[][3]~b_seb[][5]~b_seb[][7]~b_seb[][9])';
	quants = (quantiler(b-1, <0.10,0.25,0.50,0.75,0.90>));
	MAE = (quantiler(fabs(b-1), <0.50>));
	CovRate = countr((b-1)./seb,<-1.96,1.96>)[][1]/N;

	return quants~MAE~CovRate;
}
