// the following Magma code proves that the only quadratic points on X_0(88) are CM points and we find their multiplicities
// X_0(44) has finitely many quadratic points, so we just need to check whether any of the pullbacks are quadratic.


load "models_and_maps.m";
load "pullbacks.m";

X,_,c,_,_:=eqs_quos(88,[[88]]);
j := jmap(X,88); //this takes a while


//It remains to check the multiplicities of the CM points.
K:=QuadraticField(-7);
jinv:=-3375;
CMpts:=pullback_j(X,K,j,jinv); //takes a while
assert #CMpts eq 8;
jinv:=16581375;
CMpts:=pullback_j(X,K,j,jinv); //this takes a long while
assert #CMpts eq 8;