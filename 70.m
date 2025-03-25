// the following Magma code proves that there are no quadratic points on X_0(70)

load "models_and_maps.m";
load "pullbacks.m";

X,_,c,_,_:=eqs_quos(70,[[35]]);
C:=c[1,1];
J:=Jacobian(C);
assert RankBound(J) eq 0;
pts:=Chabauty0(J);
assert #pts eq 4; // there are 4 rational cusps, so all the points are cusps

// We know from results of Clark, Genao, Pollack and Saia that there are no quadratic CM points on this curve. As a sanity check, we will check that the exceptional CM point on X_0(35) does not lift to X_0(70).


d:=5;
K<w>:=QuadraticField(d);
X0:=ChangeRing(SmallModularCurve(35),K);
P:=X0![1/2*(-w-1),w+3];
jinv:=jInvariant(P,35);
E:=EllipticCurveFromjInvariant(jinv);
assert #TorsionSubgroup(E) eq 1;
// This shows that an elliptic curve has no 2-torsion, hence no 2-isogenies, hence no 70-isognies. 

