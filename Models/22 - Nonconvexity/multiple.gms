$title	Optimal Route Choice: Time and Reliability


set		pt	Points on the route /0*20/,
		xy	Planar dimensions /x,y/;

parameter	lamda	Weight on speed (1=shortest path) /1/;

set	h		Hazards /h1,h2/;

table	hloc(h,xy)	Hazard locations

	x	y
h1	0.3	0.6
h2	0.7	0.4;

parameter	var(h)	Variance /h1 0.12, h2 0.08/;

set	ig /0*10/; 
alias (ig,jg);

parameter	xg(ig)	Grid value;
xg(ig) = ig.val/10;
display xg;

parameter	dist(h,ig,jg)	Distance from hazard;
dist(h,ig,jg) = sqrt( sqr(xg(ig)-hloc(h,"x"))/2 + 
		      sqr(xg(jg)-hloc(h,"y"))/2);
option dist:2:1:1;
display dist;

parameter	hz(ig,jg)	Hazard rate;
hz(ig,jg) = 1 - prod(h, 1-exp(-0.5*sqr(dist(h,ig,jg)/var(h))));
option hz:2:1:1;
display hz;

variable	OBJ		Objective function;

nonnegative
variables	X(pt)		Route choice -- X dimension
		Y(pt)		Route choice -- Y dimension
		D		Distance between points,
		P(pt,h)		Proximity from hazard h,
		R(pt,h)		Risk from hazard h on passage from pt-1 to pt,
		Z		Overal hazard rate,
		L		Overall length of route;

equations	objdef, length, proximity, risk, hazard, distance;

*	Define the objective as a weighted sum of travel distance and probability 
*	of reading the target:

objdef..	OBJ =e= lamda * L + (1-lamda) * Z;

*	Total length of the path:

length..		L =e= sum(pt, sqrt(sqr(X(pt+1)-X(pt))+sqr(Y(pt+1)-Y(pt)))) - 1;

*	Hazard of the path = 1 - probability (no failure):

hazard..		Z =e= 1 - prod(pt, prod(h,1-R(pt,h)));

distance(pt+1)..	sqr(D) =e= sqr(X(pt+1)-X(pt)) + sqr(Y(pt+1)-Y(pt));


*	Proximity of point pt to hazard h:

proximity(pt+1,h)..	P(pt+1,h) =e= sqrt(sqr((X(pt+1)+X(pt))/2-hloc(h,"x")) + 
				           sqr((Y(pt+1)+Y(pt))/2-hloc(h,"y")));

*	Risk imposed on passage from point pt to pt+1:

risk(pt+1,h)..		R(pt+1,h) =e= exp(-0.5*sqr(P(pt+1,h)/var(h)))*D;

model routechoice /objdef, length, proximity, risk, hazard, distance/;

*	No risk on the arc leading to point 0:

R.FX("0",h) = 0;

*	Upper bound on total length:

L.UP = 2.2;

*	Points much lie in the unit square:

X.LO(pt) = 0;		X.UP(pt) = 1;
Y.LO(pt) = 0;		Y.UP(pt) = 1;

D.LO     = 0.001;	D.UP = 1;	

*	Assign some initial values:

P.L(pt,h) = 0.1;  D.L = 1/20; X.L(pt) = 1/2;	X.L(pt) = 1/2;

*	Start in the southeast (1,0) and traverse to the northwest (0,1):

X.FX("0")  = 1;		Y.FX("0")  = 0;
X.FX("20") = 0;		Y.FX("20") = 1;

option nlp=snopt;

lamda = 0.5;
solve routechoice using nlp minimizing obj;


set	scn	Angle which determines lamda (degrees) /0,30,60,90,minRisk/

set	theta(scn)	Angle which determines lamda (degrees) /0,30,60,90/
	

parameter	rc(scn,pt,xy)	Route choice,
		step(scn)	Step length;

loop(theta,

	lamda = cos(pi*theta.val/180)/(cos(pi*theta.val/180)+sin(pi*theta.val/180));

	solve routechoice using nlp minimizing obj;

*	Report the location of each step on this route:

	step(theta) = D.L;
	rc(theta,pt,"x") = X.L(pt) + eps;
	rc(theta,pt,"y") = Y.L(pt) + eps;
);
option rc:3:2:1;
display rc;

*	Find a route which runs through the northeast corner:

Y.FX("10") = 1;
X.FX("10") = 1;
lamda = 0.01;
solve routechoice using nlp minimizing obj;
*.rc("origin",pt,"x") = X.L(pt) + eps;
*.rc("origin",pt,"y") = Y.L(pt) + eps;

*	Verify that this is a locally optimal route:

D.LO = 0;
D.UP = 1;
Y.UP("10") = 1;
X.UP("10") = 1;
lamda = 0.01;
solve routechoice using nlp minimizing obj;
rc("minRisk",pt,"x") = X.L(pt) + eps;
rc("minRisk",pt,"y") = Y.L(pt) + eps;
step("minRisk") = D.L;

file kplt /'multiple.plt'/; put kplt; kplt.lw=0; file.nw=10; file.nd=3;

$onput
reset

$routes <<EOD
$offput

loop(scn,
	put '# scn=',scn.tl/;
	loop(pt,
	  put rc(scn,pt,"x"), rc(scn,pt,"y")/;
	);
	put //;
);
put 'EOD'//;


put '$hazards <<EOD'/ 
loop(ig,
	loop(jg,
	  put (ig.val/10),(jg.val/10), hz(ig,jg)/;
        );
	put /;
);
put 'EOD'//;

file.nw=0;
put 'array D[5]'/;
PUT 'D[1] = ',step("0")/;
PUT 'D[2] = ',step("30")/;
PUT 'D[3] = ',step("60")/;
PUT 'D[4] = ',step("90")/;
PUT 'D[5] = ',step("minRisk")/;

$onput
set size square
unset title
# unset key
unset colorbox
unset clabel
unset xtics
unset ytics
unset raxis
unset rtics
unset border
# set object 1 rect from 0,0 to 1,1
# set contour base
# set cntrparam levels discrete 0
# set pm3d
set xrange [-0.1:1.1]
set yrange [-0.1:1.1]
set view map
# set dgrid3d
set palette rgbformulae -34,-35,-36
set pm3d interpolate 10,10
array theta[5]
theta[1] = "distance"
theta[2] = "distance/risk"
theta[3] = "risk/distance"
theta[4] = "risk"
theta[5] = "minRisk"
set key outside right

splot   '$hazards' using 1:2:3 with pm3d notitle, \
  for [i=1:4] '$routes' index i-1 using 1:2:(0.0) with lines lw 2 lc rgb "black" dt i  title theta[i],\
  for [i=5:5] '$routes' index i-1 using 1:2:(0.0) with lines lw 2 lc rgb "dark-orange" dt 1  title theta[i]

pause -1

set terminal pngcairo size 1209,764 enhanced font 'Verdana,8'
set output 'diagonal.png'
splot   '$hazards' using 1:2:3 with pm3d notitle, \
  for [i=1:1] '$routes' index i-1 using 1:2:(0.0) with lines lw 2 lc rgb "black" dt i  title theta[i],\
  for [i=1:1] '$routes' index i-1 using 1:2:(0.0):(D[i]) with circles lw 1 lc rgb "black" notitle

set output 'safest.png'
splot   '$hazards' using 1:2:3 with pm3d notitle, \
  for [i=5:5] '$routes' index i-1 using 1:2:(0.0) with lines lw 2 lc rgb "black" dt i  title theta[i],\
  for [i=5:5] '$routes' index i-1 using 1:2:(0.0):(D[i]) with circles lw 1 lc rgb "black" notitle


set output 'sequence.png'
splot   '$hazards' using 1:2:3 with pm3d notitle, \
  for [i=1:4] '$routes' index i-1 using 1:2:(0.0) with lines lw 2 lc rgb "black" dt i  title theta[i]

set output 'localglobal.png'
splot   '$hazards' using 1:2:3 with pm3d notitle, \
  for [i=3:4] '$routes' index i-1 using 1:2:(0.0) with lines lw 2 lc rgb "black" dt i  title theta[i],\
  for [i=5:5] '$routes' index i-1 using 1:2:(0.0) with lines lw 3 lc rgb "red" dt 1  title theta[i]

$offput
putclose;

*	Generate the figures:

execute 'gnuplot multiple.plt';
