// the following Magma code proves that the X_0(90) has no quadratic points 


load "models_and_maps.m";
load "pullbacks.m";

X,_,c,_,_:=eqs_quos(90,[[45]]);
//j := jmap(X,90); //this takes a while

//We know from results of Clark, Genao, Pollack and Saia that there are no quadratic CM points on this curve. We check anyway that the quadratci CM points on X_0(45) do not lift to quadratic points on X_0(45)

K:=QuadraticField(-11);
jinv:=-32768;
E:=EllipticCurveFromjInvariant(jinv);
assert #TorsionSubgroup(E) eq 1;

//We can do the same for the non-CM points on X_0(45), although this is also not necessary, as we can conclude from the sizes of the isogeny classes that there are no 90-isogenies. 

d:=13;
K<w>:=QuadraticField(d);
X0:=ChangeRing(SmallModularCurve(45),K);
P:=X0![-2,(w+5)/2];
jinv:=jInvariant(P,45);
E:=EllipticCurveFromjInvariant(jinv);
assert #TorsionSubgroup(E) eq 1;
d:=-39;
K<w>:=QuadraticField(d);
X0:=ChangeRing(SmallModularCurve(45),K);
P:=X0![-(w+5)/4,(w+13)/8];
jinv:=jInvariant(P,45);
E:=EllipticCurveFromjInvariant(jinv);
assert #TorsionSubgroup(E) eq 1;

// We are done.