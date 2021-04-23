// 2SLS estimator
tsls(const y, const X, const Z)
{
 decl T, l, ZZ_inv, ZX, XZ, XZ_ZZinv_ZX_inv, b2sls, e, seb2sls;

 T = rows(X);
 l = columns(X);

 ZZ_inv = invertsym(Z'Z);
 ZX = Z'X;
 XZ = (ZX)';
 
 XZ_ZZinv_ZX_inv = invertsym(XZ*ZZ_inv*ZX);

 b2sls = XZ_ZZinv_ZX_inv*(XZ*ZZ_inv*Z'y);
 e = y - X*b2sls;
 seb2sls = sqrt(diagonal(XZ_ZZinv_ZX_inv*(e'e/(T-l))))';

 return b2sls~seb2sls;
}
 
