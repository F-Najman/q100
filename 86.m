load "models_and_maps.m";
load "pullbacks.m";

////////////////////////////////////////////////////////////////////////////////////////////
N:=86;
X, _, pairs := eqs_quos(N, [[43]]);
j := jmap(X,86); //this takes a while






\\We determine the multilpicites of the CM points now

K:=QuadraticField(-3);
jinv:=-0;
CMpts:=pullback_j(X,K,j,jinv);
assert #CMpts eq 2;
jinv:=54000;
CMpts:=pullback_j(X,K,j,jinv);
assert #CMpts eq 2;

K:=QuadraticField(-7);
jinv:=-3375;
CMpts:=pullback_j(X,K,j,jinv); 
assert #CMpts eq 6;
jinv:=16581375;
CMpts:=pullback_j(X,K,j,jinv); 
assert #CMpts eq 2;

K:=QuadraticField(-2);
jinv:=8000;
CMpts:=pullback_j(X,K,j,jinv);
assert #CMpts eq 2;


//We now check the exceptional points on $X_0(43). These can be found in the Box's paper. 

K<w>:=QuadraticField(-71);
jinv:=(-49*w-977)/4;
pts:=pullback_j(X,K,j,jinv);
E:=EllipticCurveFromjInvariant(jinv);
assert #TorsionSubgroup(E) eq 1;


K<w>:=QuadraticField(-131);
jinv:=(245508467396487686583118*w-245508467396487686583118)/245508467396487686583118;
E:=EllipticCurveFromjInvariant(jinv);
assert #TorsionSubgroup(E) eq 1;
