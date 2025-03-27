load "models_and_maps.m";
load "pullbacks.m";
load "coleman.m";
load "model_equation_finder.m";

////////////////////////////////////////////////////////////////////////////////////////////
N:=86;
X0, _, pairs := eqs_quos(N, [[43]]);
C:=pairs[1,1]; // this is X_0(86)/w43
//We need to compute the rational points on X_0(86)/w_43
pts:=PointSearch(C,1000);
assert #pts eq 7;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//We now compute the rational points on X_0(86)/w43
//We first need a nice model for Chabauty
 P3<W,X,Y,Z> := ProjectiveSpace(Rationals(),3);
X86_w43 := X0NQuotient(86, [43]); // to compute the rational points, this model turns out to be much better.
//IsIsomorphic(C, X86_w43); //This takes a couple of minutes to compute, returns true.
Q := find_and_test_model_Q_only(W, W+2*X-Y, Y, W+2*X-Y,X86_w43, x, y); 
Q;
// this is one monic model for X0(86)_w43


p := 13; 
N := 12;
data := coleman_data(Q,p,N);


Qpoints := Q_points(data,1000);
//Assuming that the Q-rational points found generate a finite index subgroup of the Jacobian of X_0(74)/w37,
// the following code proves that these are all Q-points.
// This assumption is checked at the end of this file.

L, v := effective_chabauty(data : Qpoints := Qpoints, e := 40);

if #L eq #Qpoints then
  printf "found all %o Q-points!\n", #Qpoints;
else
  printf "one has to exclude additional %o points\n", #L - #Qpoints;
end if;

//found all 7 Q-points!


// We aim to show that the divisors formed by d1 = pseq1-pseq2 and d2 = pseq1-pseq3 generate a rank 2 subgroup of the Jacobian over Q
// We do this in three steps 

// Write J for the Jacobian of the curve X_0(86)/w_43

// Step 1: We prove that J(Q)_tors is a subgroup of Z/11Z
// Step 2: We verify that d_1 and d_2 both have infinite order by considering their reductions modulo 3 and 5 
// Step 3: We verify that d_1 and d_2 are linearly independent points of the Jacobian




X5 := ChangeRing(C, GF(5));
assert IsNonsingular(X5);
PicX5, phi5, psi5 := ClassGroup(X5);
JF5 := TorsionSubgroup(PicX5);

X7 := ChangeRing(C, GF(7));
assert IsNonsingular(X7);
PicX7, phi7, psi7 := ClassGroup(X7);
JF7 := TorsionSubgroup(PicX7);

X17 := ChangeRing(C, GF(17));
assert IsNonsingular(X17);
PicX17, phi17, psi17 := ClassGroup(X17);
JF17 := TorsionSubgroup(PicX17);

assert GCD([#JF5,#JF7,#JF17]) eq 11;

pts; 
// [ (-1 : 1 : 0 : 0), (1 : 1 : -1 : 1), (-1 : 1 : -1 : 1), (1 : 1 : 0 : 0), (0 : 0: 1 : 0), (2 : 0 : -1 : 1), (-2 : 0 : -1 : 1) ]

pseq:= [[-1,1,0,0], [1,1,-1,1], [-1,1,-1,1], [1 , 1 , 0 , 0], [0 , 0, 1 , 0], [2 , 0 , -1 , 1], [-2 , 0 , -1 , 1]];

n1:=1;
n2:=6;
n3:=3;
n4:=7;


d1_mod5 := psi5(Place(X5 ! pseq[n1]) - Place(X5 ! pseq[n2]));
d2_mod5 := psi5(Place(X5 ! pseq[n3]) - Place(X5 ! pseq[n4]));

d1_mod7 := psi7(Place(X7 ! pseq[n1]) - Place(X7 ! pseq[n2]));
d2_mod7 := psi7(Place(X7 ! pseq[n3]) - Place(X7 ! pseq[n4]));

d1_mod17 := psi17(Place(X17 ! pseq[n1]) - Place(X17 ! pseq[n2]));
d2_mod17 := psi17(Place(X17 ! pseq[n3]) - Place(X17 ! pseq[n4]));

// If d1 were a torsion point then it would have the same order mod 7 and mod 5.
assert Order(d1_mod7) ne Order(d1_mod5); // So d1 has infinite order
// Similarly for d2:
assert Order(d2_mod7) ne Order(d2_mod5); // So d2 has infinite order

// Step 3:

// If d1 and d2 were linearly dependent, 
// then they would generate a group isomorphic to Z x G. 
// Here, Z means the integers and
// G is a subgroup of J(Q)_tors, which in turn is a subgroup of Z/11Z

// It follows that the image of <d_1, d_2> in J(F_5) must be of the form:
// Z/aZ or Z/aZ x Z/11Z for some integer a.

// We compute the image of <d_1, d_2> in J(F_5):

A5 := sub<JF5 | [d1_mod5, d2_mod5] >; 
A7 := sub<JF7 | [d1_mod7, d2_mod7] >; 
A17 := sub<JF17 | [d1_mod17, d2_mod17] >;
//Invariants(A5);Invariants(A7);Invariants(A17);
assert IsIsomorphic(A5,AbelianGroup([10,110]));
// This group is not of the form Z/aZ or Z/aZ x Z/11Z
// So the points d_1 and d_2 are linearly independent.



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//We determine the multilpicites of the CM points now
X:=X0;
j := jmap(X,86); //this takes a while
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

