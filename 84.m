// the following Magma code proves that the only quadratic points on X_0(84) are CM points and we find their multiplicities
// X_0(42) has finitely many quadratic points, so we just need to check whether any of the pullbacks are quadratic.


load "models_and_maps.m";
load "pullbacks.m";

X,_,c,_,_:=eqs_quos(84,[[84]]);
j := jmap(X,84); //this takes a while

K:=QuadraticField(-3);
jinv:=0;
CMpts:=pullback_j(X,K,j,jinv); //takes a while
assert #CMpts eq 0; 
jinv:=54000;
CMpts:=pullback_j(X,K,j,jinv); // takes a while
assert #CMpts eq 4; 