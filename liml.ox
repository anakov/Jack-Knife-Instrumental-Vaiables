//  LIML estimator
//  Ox code by Anton Nakov, 2010
liml(const y1, const Y2, const X, const Z)
{
 decl T, l, r, Mzx, Mx, ZX, Y, Y2X, Y_ZX, YY, YX, Y2X_ZX, inv_ZXZX, eigpart, part, beta, e, sebeta, kappa, eigval;
 
 T = rows(y1);
 l = columns(Y2);
 r = columns(X);
 eigval = ones(T,1);

 ZX       = Z~X;
 Y        = y1~Y2;
 Y2X      = Y2~X;
 Y_ZX     = Y'(ZX);
 YY       = Y'Y;
 YX       = Y'X;
 Y2X_ZX   = (Y2X)'(ZX);
 inv_ZXZX = invertsym((ZX)'(ZX));			  
 
 if (X==<>) eigpart  =  YY;
 else       eigpart  = (YY - (YX)*invertsym(X'X)*(YX'));
	 
 eigen(invert(YY - Y_ZX*inv_ZXZX*Y_ZX') * eigpart, &eigval);
 kappa = min(eigval); // println(kappa);
 
 part  = invert( (Y2X)'(Y2X) - kappa*( (Y2X)'(Y2X) - Y2X_ZX*inv_ZXZX*Y2X_ZX' ) );
 beta  = part *( (Y2X)'(y1)  - kappa*( (Y2X)'(y1)  - Y2X_ZX*inv_ZXZX*((ZX)'(y1)) ) );
 e = y1 - (Y2X)*beta;
 sebeta = sqrt(diagonal(part*(e'e/(T-l-r))))';
 
 return beta~sebeta;
}