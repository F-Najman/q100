//This case (n=93) is missing from the paper. 
//It can be easily solved using the methods of the paper, and it turns out X_0(93) has 12 quadratic points, all of them CM.


//We apply Corollary 3.3. from the paper. The values in question are: n=31, d=31, m_2=1, m_1=3. 
//We just need to find the rational points on the genus 5 curve X_0(93)/w_{31}. The curve X_0(93)* is a genus 2 curve, whose Jacobian has rank 2. 
// The rational points on this curve were determined in:
// Adžaga, Chidambaram, Keller, Padurariu, Rational points on hyperelliptic Atkin-Lehner quotients of modular curves and their coverings, 2022.
// To get a quadratic point on X_0(93), it is necessary (but not sufficient) that the j-invariant of the points listed in the paper are defined over Q or a quadratic field.
// This immediately rules out the non-CM points. Conclusion: the P^1-parametrized non-CM quadratic points on X_0(31) do not lift to quadratic points on X_0(93).

// We now check whether the P^1-isolated points lift. We look at the exceptional points found in:
// P. Bruin, F. Najman, Hyperelliptic modular curves X_0(n) and isogenies of elliptic curves over quadratic fields, LMS J. Comput. Math. 18 (2015) 578-602.

C:=SmallModularCurve(31);
K<w>:=QuadraticField(-3);
C1:=ChangeRing(C,K);
P1:=C1![1/2*(w-1),-2];
P2:=C1![1/2*(w-1), 1/2*(w+7)];
_<x>:=PolynomialRing(K);
j1:=jInvariant(P1,31);
j2:=jInvariant(P2,31);
//These points will lift to a point on X_0(93) if and only if the corresponding curve has a 3-isogeny. This is true if and only if the modular polynomial evaluated with one variable specialized to the j-invariant has a root over K.
p:=ChangeRing(ClassicalModularPolynomial(3), K);
assert IsIrreducible(Evaluate(p,[x,j1]));
assert IsIrreducible(Evaluate(p,[x,j2]));
//This proves that P^1-isolated non-CM quadratic points on X_0(31) do not lift to quadratic points on X_0(93).

//Now from: https://github.com/fsaia/least-cm-degree/blob/master/Least%20Degrees/X0/dcm_list_all_min_orders_X0_10k.m we see that the j-invariants that give quadratic CM points on X_0(93) are:
//[ 0, 54000, -12288000, -32768 ]


load "models_and_maps.m";
load "pullbacks.m";
N:=93;
X,u,v,c,d:=eqs_quos(93,[[31]]);
j := jmap(X,N);//~1 min

CM_j:= [ 0, 54000, -12288000, -32768 ];
for jinv in CM_j do
    if jinv eq -32768 then d:=-11; else d:=-3; end if;
    K:=QuadraticField(d);
    CMpts:=pullback_j(X,K,j,jinv);
    jinv, #CMpts;
end for;
/*takes some time and returns:
0 4
54000 2
-12288000 2
-32768 4
*/


