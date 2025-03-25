// the following Magma code proves that the only quadratic points on X_0(78) are CM points and we find th

load "models_and_maps.m";
load "pullbacks.m";

X,_,c,_,_:=eqs_quos(78,[[39]]);
C:=c[1,1];
J:=Jacobian(C);
assert RankBound(J) eq 0;
pts:=Chabauty0(J);
assert #pts eq 6; // there are 4 rational cusps, so there are 2 non-cuspidal points

j := jmap(X,78); //this takes a while



//We now check the multiplicities of the CM points with CM by orders of discriminatn -3, -12 and -27
K:=QuadraticField(-3);
jinv:=0;
CMpts:=pullback_j(X,K,j,jinv);
assert #CMpts eq 2; 
jinv:=54000;
CMpts:=pullback_j(X,K,j,jinv);
assert #CMpts eq 2; 
jinv:=-12288000;
CMpts:=pullback_j(X,K,j,jinv);
assert #CMpts eq 0; 

//Finally, although from the fact that there are only 4 points in the isogeny class with 39-isogenies for the exceptional points on X_0(39), we do a sanity check and double check that there is no 2-torsion in this isogeny class. 

d:=-7;
K<w>:=QuadraticField(d);
X0:=ChangeRing(SmallModularCurve(39),K);
P:=X0![1/4*(w+3),1/8*(3*w+1)];
jinv:=jInvariant(P,35);
E:=EllipticCurveFromjInvariant(jinv);
assert #TorsionSubgroup(E) eq 1;

//This completes the proof.

