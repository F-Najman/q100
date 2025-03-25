//load "rank_0.m";
//rank_0_quad_pts(96);
load "models_and_maps.m";
X,u,v,c,d:=eqs_quos(96,[[96]]);
Y:=v[1,1];
assert Genus(Y) eq 3;
assert IsHyperelliptic(Y) eq false;
ClassGroup(ChangeRing(Y,GF(11)));
ClassGroup(ChangeRing(Y,GF(5)));
//we conclude J(Y)<=4x4x8


pts:=PointSearch(Y,100);
load "quadPts.m";
deg2,pls1,pls2:=searchDiv2(Y, 1, true);
deg1:=[Divisor(Place(p)): p in pts];
divs:= deg1 cat deg2;
findGenerators(Y,divs,Divisor(Place(pts[1])),5);
