//  OLS estimator
ols (const y, const X)
{
 decl T, l, XX_inv, beta, e, se_beta;

 T = rows(X);
 l = columns(X);

 XX_inv = invertsym(X'X);

 beta = XX_inv*(X'y);
 e = y - X*beta;
 se_beta = sqrt(diagonal(XX_inv*(e'e/(T-l))))';

 return beta~se_beta;
}
