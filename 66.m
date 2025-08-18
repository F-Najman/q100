// the following Magma code proves that there are only two (pairs of Galois conjugated) quadratic points on X_0(66), both CM. 

 load "models_and_maps.m";
 X,_,c,_,_:=eqs_quos(66,[[11]]);
 C:=c[1,1]; //C is X_0(66)/w_11
 J:=Jacobian(C);
 RankBounds(J); // this proves that J(C) has rank 0 over Q
 pts:=Chabauty0(J);
 assert #pts eq 4; //there are 4 cusps, so we are done, as X_0(66)/w_11 has 4 cusps. Thus X_0(66)/w_11 has no non-cuspidal rational points
 
 
 //We now need to check whether the quadratic CM points of X_0(33) lift to quadratic points on X_0(66) and check their multiplicities. 
 load "pullbacks.m";
 
j := jmap(X,66); //this takes a while
K:=QuadraticField(-2);
jinv:=8000;
CMpts:=pullback_j(X,K,j,jinv);
assert #CMpts eq 4; 
// this confirms that there are 2 conjugacy classes of quadratic points with j=8000


K:=QuadraticField(-11);
jinv:=-32768;
CMpts:=pullback_j(X,K,j,jinv);
assert #CMpts eq 0; 
// We know from results of Clark, Genao, Pollack and Saia that there are no quadratic points with j=-32768, but we do a sanity check;


//This completes the P^1 parametrized points on X_0(33). It remains to chek that the exceptional points quadratic on X_0(33) do not lift to quadratic points on X_0(66).

d:=-2;
K<w>:=QuadraticField(d);
X0:=ChangeRing(SmallModularCurve(33),K);
P:=X0![w-1,7*w-2];
jinv:=jInvariant(P,33);
E:=EllipticCurveFromjInvariant(jinv);
assert #TorsionSubgroup(E) eq 1;
// This shows that an elliptic curve has no 2-torsion, hence no 2-isogenies, hence no 66-isognies. We need to check only one curve per isogeny class.

d:=-7;
K<w>:=QuadraticField(d);
X0:=ChangeRing(SmallModularCurve(33),K);
P:=X0![1/2*(w+1),-1];
jinv:=jInvariant(P,33);
E:=EllipticCurveFromjInvariant(jinv);
assert #TorsionSubgroup(E) eq 1;

// Again, there is no 2-torsion. This completes the proof.



