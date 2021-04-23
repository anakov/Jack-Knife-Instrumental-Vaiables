//  Jack-Knife Instrumental Variables Estimator JIVE1 (Angrist, Imbens, and Krueger, JAE 1999)
//  Ox code by Anton Nakov, 2010

jive1 (const y, const X, const Z)
{
 decl T, l, Xjive1, h, ZZ_inv, ZX, pi, i, Xj1X, bjive1, e, Xj1Xj1, sebjive1;
 T = rows(X);
 l = columns(X);
 Xjive1 = zeros(T,columns(X));
 h = zeros(T,1);
 ZZ_inv = invertsym(Z'Z);
 ZX = Z'X;
 pi = ZZ_inv*ZX;
 for (i = 0; i < T; i++)
 {
  h[i] = Z[i][]*ZZ_inv*Z[i][]';
  Xjive1[i][] = (Z[i][]*pi - h[i]*X[i][])/(1-h[i]);
 }
 Xj1X = Xjive1'X;

 bjive1 = invert(Xj1X)*(Xjive1'y);

 e = y - X*bjive1;
 Xj1Xj1 = invertsym(Xjive1'Xjive1);
 sebjive1 = sqrt(diagonal(invertsym(Xj1X'*Xj1Xj1*Xj1X)*(e'e/(T-l))))';
	
 return bjive1~sebjive1;
}