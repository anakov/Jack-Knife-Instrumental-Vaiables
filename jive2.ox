//  Jack-Knife Instrumental Variables Estimator JIVE2 (Angrist, Imbens, and Krueger, JAE 1999)
//  Ox code by Anton Nakov, 2010

jive2 (const y, const X, const Z)
{
 decl T, l, Xjive2, h, ZZ_inv, ZX, pi, i, Xj2X, bjive2, e, Xj2Xj2, sebjive2;
 T = rows(X); 
 l = columns(X);	
 Xjive2 = zeros(T,columns(X));
 h = zeros(T,1);
 ZZ_inv = invertsym(Z'Z);
 ZX = Z'X;
 pi = ZZ_inv*ZX;
 for (i = 0; i < T; i++)
 {
  h[i] = Z[i][]*ZZ_inv*Z[i][]';
  Xjive2[i][] = (Z[i][]*pi - h[i]*X[i][])/(1-1/T);
 }
 Xj2X = Xjive2'X;
	
 bjive2 = invert(Xj2X)*(Xjive2'y);

 e = y - X*bjive2;
 Xj2Xj2 = invertsym(Xjive2'Xjive2);
 sebjive2 = sqrt(diagonal(invertsym(Xj2X'*Xj2Xj2*Xj2X)*(e'e/(T-l))))';

 return bjive2~sebjive2;
}