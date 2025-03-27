
// the following Magma code proves that the only quadratic points on X_0(82) are CM points and we find their multiplicities

load "models_and_maps.m";
load "pullbacks.m";

X,_,c,_,_:=eqs_quos(82,[[41]]);
C:=c[1,1];
assert Genus(C) eq 3;
pts:=PointSearch(C,10);
assert #pts eq 6; // there are 4 rational cusps, so there are 2 non-cuspidal points
j := jmap(X,82); //this takes a while

bound := 10000;
pullbacks := pullback_points(X,c,82, bound); 

print "Number of pairs of quadratic points found as pullbacks =", #pullbacks;
print "++++++++++++++++++";

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
Output:  P coordinates: (-2*T : 0 : 0 : 1 : 0 : -T : 0 : 0 : 1) where T^2 = -1 and
j-invariant = 287496
and CM by -16
P coordinates: (0 : 0 : 0 : 0 : -T : -1/2*T : 4 : 2 : 1) where T^2 = -1 and
j-invariant = 1728
and CM by -4
P coordinates: (-2*T : 0 : 0 : -1 : 0 : T : 0 : 0 : 1) where T^2 = -1 and
j-invariant = 1728
and CM by -4
P coordinates: (0 : 0 : 0 : 0 : 0 : -1/2*T : 1 : 1 : 1) where T^2 = -2 and
j-invariant = 8000
and CM by -8
*/





// the following proves that the quotient has only 6 rational points

//We load the files to do classical Chabauty

load "coleman.m";

K<u>:=PolynomialRing(Integers());
_<v>:=PolynomialRing(K);
z:=1;
f:=u^4 - 2*u^2*v^2 - 4*u^2*v*z - 4*u^2*z^2 + v^4 -
    12*v^3*z + 4*v^2*z^2 + 32*v*z^3;
data:=coleman_data(f,3,10);
Qpoints:=Q_points(data,10^4);
assert #Qpoints eq 6; 
L, v := effective_chabauty(data : Qpoints := Qpoints, e := 40);

//printf "L = %o\nQ-points = %o\nv = %o\n", L, Qpoints, v;
//uncomment the previous line to get a bit more info
//about annihilating differentials and their zeroes (candidate points)

if #L eq #Qpoints then
  printf "found all %o Q-points!\n", #Qpoints;
else
  printf "one has to exclude additional %o points\n", #L - #Qpoints;
end if;


//Assuming that the Q-rational points found generate a finite index subgroup of the Jacobian of X_0(82)/w41,
// we have found all the points. Let's check that. 


X3 := ChangeRing(C, GF(3));
assert IsNonsingular(X3);
PicX3, phi3, psi3 := ClassGroup(X3);

X5 := ChangeRing(C, GF(5));
assert IsNonsingular(X5);
PicX5, phi5, psi5 := ClassGroup(X5);

pseq1 := [2,0,1];
pseq2 := [0,0,1];

d1_mod3 := psi3(Place(X3 ! pseq1) - Place(X3 ! pseq2));
d1_mod5 := psi5(Place(X5 ! pseq1) - Place(X5 ! pseq2));

// If d1 were a torsion point then it would have the same order mod 3 and mod 5.
assert Order(d1_mod3) ne Order(d1_mod5); // So d1 has infinite order

// This shows that we have found all the rational points on X_0(84)/w42




// We conclude that all the pullbacks are CM. We need to check the multiplicities of all CM curves that appear, in theory there could be more points than the ones that we already found. 

K:=QuadraticField(-1);
jinv:=1728;
CMpts:=pullback_j(X,K,j,jinv);
assert #CMpts eq 4; 
jinv:=287496;
CMpts:=pullback_j(X,K,j,jinv);
assert #CMpts eq 2; 
K:=QuadraticField(-2);
jinv:=8000;
CMpts:=pullback_j(X,K,j,jinv);
assert #CMpts eq 2; 

//Now we need to check that the exceptional points on X_0(41) do not have any 2-isogenies. 

d:=-1;
K<w>:=QuadraticField(d);
X0:=ChangeRing(SmallModularCurve(41),K);
P:=X0![1/2*(-w-1),1/4*(w+1)];
jinv:=jInvariant(P,41);
E:=EllipticCurveFromjInvariant(jinv);
assert #TorsionSubgroup(E) eq 1;

//This completes the proof. The curves have no 2-torsion, hence no 2-isogenies, hence no 82-isognies.













