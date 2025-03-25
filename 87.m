// the following Magma code proves that the only quadratic points on X_0(87) are CM points and we find their multiplicities

X,_,c,_,_:=eqs_quos(87,[[29]]);
C:=c[1,1];
J:=Jacobian(C);
assert RankBound(J) eq 0;
pts:=Chabauty0(J);
assert #pts eq 2; // there are 2 rational cusps, so there are no non-cuspidal points

//There are no CM points on X_0(87) by the results of. Clark, Genao, Pollack and Saia. 
//We need to check the exceptional points on X_0(29) 


d:=-1;
K<w>:=QuadraticField(d);
X0:=ChangeRing(SmallModularCurve(29),K);
P:=X0![w-1,2*w+4];
jinv:=jInvariant(P,29);
E:=EllipticCurveFromjInvariant(jinv);
fac:=Factorization(DivisionPolynomial(E,3));
assert Degree(fac[1,1]) eq 4;
// The 3-division polynomial of E is irreducible, so there are no 3-isogenies, hence no 87-isogenies

d:=-7;
K<w>:=QuadraticField(d);
X0:=ChangeRing(SmallModularCurve(29),K);
P:=X0![1/4*(w+1),1/16*(-11*w-7)];
jinv:=jInvariant(P,29);
E:=EllipticCurveFromjInvariant(jinv);
fac:=Factorization(DivisionPolynomial(E,3));
assert Degree(fac[1,1]) eq 4;
// The 3-division polynomial of E is irreducible, so there are no 3-isogenies, hence no 87-isogenies


//This completes the proof