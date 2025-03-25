load "models_and_maps.m";
load "pullbacks.m";

////////////////////////////////////////////////////////////////////////////////////////////
N:=99;
X, _, pairs := eqs_quos(N, [[11]]);
bound := 10000;
pullbacks := pullback_points(X,pairs,N, bound);
print "Number of pairs of quadratic points found as pullbacks =", #pullbacks;
;
//Number of pairs of quadratic points found as pullbacks = 1

print "++++++++++++++++++";
j := jmap(X,N); // Runtime: ~ 10 minutes
for tup in pullbacks do
    K<T> := tup[2];
    P := X(K) ! Eltseq(tup[1][1]);
    T2 := -Coefficient(DefiningPolynomial(Ring(Parent(P))),0); // T^2 = this
    jP := j(P)[1];
    if jP eq 1 and j(P)[2] eq 0 then 
        print "P coordinates:", P, "where T^2 =", T2, "and the point is a cusp";
        continue;
    end if;
    tf, D := HasComplexMultiplication(EllipticCurveWithjInvariant(jP));
    assert tf;
    print "P coordinates:", P, "where T^2 =", T2, "and j-invariant =", jP, "and CM by", D;    
end for;

/*
P coordinates: (0 : -4/33*T : 2/33*T : -1/11*T : -1/33*T : -1 : 1 : 0 : 0) where
T^2 = 33 and j-invariant = -3274057859072*T - 18808030478336
and CM by -99
*/

//////////////////////////////////////////////////////////////////////////////////////////////

//Now we find all the rational points on X_0(99)/w_11

C:=pairs[1,1];// this is X_0(99)/w_11
pts:=PointSearch(C,100);
assert #pts eq 3;

load "coleman.m";

K<v>:=PolynomialRing(Integers());
_<u>:=PolynomialRing(K);
z:=1;
f:=u^4 + 2*u^2*v^2 + 8*u^2*v*z + 14*u^2*z^2 - 3*v^4 -
    6*v^2*z^2 + 24*v*z^3 + 33*z^4; // this is the equation of X_0(99)/w_11 again
data:=coleman_data(f,5,10);
Qpoints:=Q_points(data,10^4);
assert #Qpoints eq 3;

//this proves that X_0(99)/w_11 has 3 rational points 

// We can now conclude that the only possible quadratic points are the CM  points, pullbacks of rational points (which we've seen are CM), and possibly pullbacks of exceptional points. 

//We now check whether the exceptional points on X_0(33) have any 99-isogenies. If there existed a 9-isogeny then the primitive 9-division polynomial (i.e 9-division polynomial quotiented out by the 3-division polynomial would have to have a factor of degree dividing 3. Since it's possible for some curves to have 9-isogenies, but not all, we need to check all the curves in the isogeny class.


d:=-2;
K<w>:=QuadraticField(d);
X0:=ChangeRing(SmallModularCurve(33),K);
P:=X0![-1/2*w,1/4*(-4*w-5)];
jinv:=jInvariant(P,33);
E:=EllipticCurveFromjInvariant(jinv);
fac:=Factorization(DivisionPolynomial(E,9) div DivisionPolynomial(E,3)); 
assert Degree(fac[1,1]) eq 9;

P:=X0![w-1,-5*w-5];
jinv:=jInvariant(P,33);
E:=EllipticCurveFromjInvariant(jinv);
fac:=Factorization(DivisionPolynomial(E,9) div DivisionPolynomial(E,3)); 
assert Degree(fac[1,1]) eq 9;

P:=X0![-1/2*w,w+2];
jinv:=jInvariant(P,33);
E:=EllipticCurveFromjInvariant(jinv);
fac:=Factorization(DivisionPolynomial(E,9) div DivisionPolynomial(E,3)); 
assert Degree(fac[1,1]) eq 9;

P:=X0![w-1,7*w-2];
jinv:=jInvariant(P,33);
E:=EllipticCurveFromjInvariant(jinv);
fac:=Factorization(DivisionPolynomial(E,9) div DivisionPolynomial(E,3)); 
assert Degree(fac[1,1]) eq 9;

d:=-7;
K<w>:=QuadraticField(d);
X0:=ChangeRing(SmallModularCurve(33),K);
P:=X0![1/4*(-3*w+1),1/32*(9*w+93)];
jinv:=jInvariant(P,33);
E:=EllipticCurveFromjInvariant(jinv);
fac:=Factorization(DivisionPolynomial(E,9) div DivisionPolynomial(E,3)); 
assert Degree(fac[1,1]) eq 9;

P:=X0![1/2*(w+1),-1];
jinv:=jInvariant(P,33);
E:=EllipticCurveFromjInvariant(jinv);
fac:=Factorization(DivisionPolynomial(E,9) div DivisionPolynomial(E,3)); 
assert Degree(fac[1,1]) eq 9;

P:=X0![1/4*(-3*w+1),1/4*(9*w+33)];
jinv:=jInvariant(P,33);
E:=EllipticCurveFromjInvariant(jinv);
fac:=Factorization(DivisionPolynomial(E,9) div DivisionPolynomial(E,3)); 
assert Degree(fac[1,1]) eq 9;

P:=X0![1/2*(w+1),-w+1];
jinv:=jInvariant(P,33);
E:=EllipticCurveFromjInvariant(jinv);
fac:=Factorization(DivisionPolynomial(E,9) div DivisionPolynomial(E,3)); 
assert Degree(fac[1,1]) eq 9;





//Finally, we check the multiplicites of the CM points. 


K:=QuadraticField(-2);
jinv:=8000;
CMpts:=pullback_j(X,K,j,jinv); //takes a while
assert #CMpts eq 4;

K:=QuadraticField(-11);
jinv:=-32768;
CMpts:=pullback_j(X,K,j,jinv); //takes a while
assert #CMpts eq 4;

K<T>:=QuadraticField(33);
jinv:=-3274057859072*T - 18808030478336; //takes a while
CMpts:=pullback_j(X,K,j,jinv);
assert #CMpts eq 1; // here note that jinv conjugated also gives this number of points

//this completes the proof







*/ THIS IS redundant!

N:=99;
X, _, pairs := eqs_quos(N, [[9]]);
bound := 10000;
pullbacks := pullback_points(X,pairs,N, bound);
print "Number of pairs of quadratic points found as pullbacks =", #pullbacks
;
//Number of pairs of quadratic points found as pullbacks = 4

print "++++++++++++++++++";
j := jmap(X,N); // Runtime: ~ 10 minutes
for tup in pullbacks do
    K<T> := tup[2];
    P := X(K) ! Eltseq(tup[1][1]);
    T2 := -Coefficient(DefiningPolynomial(Ring(Parent(P))),0); // T^2 = this
    jP := j(P)[1];
    if jP eq 1 and j(P)[2] eq 0 then 
        print "P coordinates:", P, "where T^2 =", T2, "and the point is a cusp";
        continue;
    end if;
    tf, D := HasComplexMultiplication(EllipticCurveWithjInvariant(jP));
    assert tf;
    print "P coordinates:", P, "where T^2 =", T2, "and j-invariant =", jP, "and CM by", D;    
end for;



*/P coordinates: (-1/3*T : 2/9*T : -1/9*T : 0 : -1/9*T : -1/3 : 0 : 1 : 0) where
T^2 = -3 and the point is a cusp
P coordinates: (0 : -4/33*T : 2/33*T : -1/11*T : -1/33*T : -1 : 1 : 0 : 0) where
T^2 = 33 and j-invariant = -3274057859072*T - 18808030478336
and CM by -99
P coordinates: (0 : -3/11*T : 1/11*T : 1/11*T : 1/11*T : 0 : 0 : 1 : 1) where
T^2 = -11 and j-invariant = -32768
and CM by -11
P coordinates: (-1/3*T : -2/9*T : 1/9*T : 0 : 1/9*T : 1/3 : 0 : 1 : 0) where T^2
= -3 and the point is a cusp
*/

*/








































//probably only nedd X_0(99)/w_11!!
//again rank 1

d:=11;
m:=9;

load "new_models.m";
a,b,c,h,e:=eqs_quos(d*m,[[d]]);
C:=c[1,1];
IsHyperelliptic(C);
Genus(C);
pts:=Points(C:Bound:=10);
assert #pts eq 3;
//aut:=Automorphisms(C);
// genus 3 rank>0
// ne želi izračunati quotient
// aut:=Automorphisms(C);
// CurveQuotient(AutomorphismGroup(C,[aut[2]]));




























%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%STARO

X, _, pairs := eqs_quos(99, [[99]]);
bound := 10000;
pullbacks := pullback_points(X,pairs,N, bound);
print "Number of pairs of quadratic points found as pullbacks =", #pullbacks\
;
//Number of pairs of quadratic points found as pullbacks = 3

print "++++++++++++++++++";
j := jmap(X,N); // Runtime: ~ 10 minutes
for tup in pullbacks do
    K<T> := tup[2];
    P := X(K) ! Eltseq(tup[1][1]);
    T2 := -Coefficient(DefiningPolynomial(Ring(Parent(P))),0); // T^2 = this
    jP := j(P)[1];
    if jP eq 1 and j(P)[2] eq 0 then 
        print "P coordinates:", P, "where T^2 =", T2, "and the point is a cusp";
        continue;
    end if;
    tf, D := HasComplexMultiplication(EllipticCurveWithjInvariant(jP));
    assert tf;
    print "P coordinates:", P, "where T^2 =", T2, "and j-invariant =", jP, "and CM by", D;    
end for;



*/

P coordinates: (-1 : -1/2*T : 0 : 1/4*T : 1/4*T : 1/4*T : -1/4*T : 0 : 1) where
T^2 = -2 and j-invariant = 8000
and CM by -8
P coordinates: (1 : -1/2*T : 0 : 1/4*T : 1/4*T : -1/4*T : 1/4*T : 0 : 1) where
T^2 = -2 and j-invariant = 8000
and CM by -8
P coordinates: (0 : -3/11*T : 1/11*T : 1/11*T : 1/11*T : 0 : 0 : 1 : 1) where
T^2 = -11 and j-invariant = -32768
and CM by -11

*/


