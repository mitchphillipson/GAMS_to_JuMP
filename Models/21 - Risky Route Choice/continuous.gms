$title	Optimal Route Choice: Time and Reliability

*	Define the number of rows/columns in the square grid:

parameter	w	Width of the grid/20/;

$if not set npt $set npt 20

set		pt	Points on the route /0*%npt%/;

parameter	lamda	Weight on speed (1=shortest path) /1/;

*	Start in the southeast and traverse to the northwest:
  
variable	OBJ		Objective function;

nonnegative
variables	X(pt)		Route choice -- X dimension
		Y(pt)		Route choice -- Y dimension
		D		Distance between points;

equations	objdef, distance;

*	Define the objective as a weighted sum of travel distance and probability 
*	of reading the target:

objdef..	OBJ =e= sum(pt, lamda * sqrt(sqr(X(pt+1)-X(pt))+sqr(Y(pt+1)-Y(pt))) 

		- (1-lamda) * log( sqrt(sqr((X(pt+1)+X(pt))/2) + sqr((Y(pt+1)+Y(pt))/2))));

distance(pt+1)..  sqr(X(pt+1)-X(pt)) + sqr(Y(pt+1)-Y(pt)) =e= sqr(D);

model routechoice /objdef, distance/;

X.LO(pt) = 0;
X.UP(pt) = w;
Y.LO(pt) = 0;
Y.UP(pt) = w;
D.UP = w;
D.L = w/card(pt);
D.LO = 0.1/card(pt);
D.L = 0;
X.L(pt) = w/2;
X.L(pt) = w/2;
X.FX("0") = w;
Y.FX("0") = 0;
X.FX("%npt%") = 0;
Y.FX("%npt%") = w;

lamda = 0.5;
solve routechoice using nlp minimizing obj;

set	theta	Angle which determines lamda (degrees) /0,10,20,30,40,50,60,70,80,90/,
	xy	Planar dimensions /x,y/;

parameter	rc(theta,pt,xy)	Route choice;

loop(theta,

	lamda = cos(pi*theta.val/180)/(cos(pi*theta.val/180)+sin(pi*theta.val/180));

	solve routechoice using nlp minimizing obj;

	rc(theta,pt,"x") = X.L(pt) + eps;
	rc(theta,pt,"y") = Y.L(pt) + eps;
);
option rc:3:2:1;
display rc;

file kplt /'continuous.plt'/; put kplt; kplt.lw=0; kplt.nw=0; kplt.nd=0;

$onput
reset
set size square
unset xtics
unset ytics
unset raxis
unset rtics
unset border
$offput

put	'set xrange [0:',w,']'/
	'set yrange [0:',w,']'/
	'set object 1 rect from 0,0 to ',w,',',w//
	'$routes <<EOD'/;

loop(theta,
	kplt.nw=0; kplt.nd=1;
	put '# theta=',(theta.val)/;
	kplt.nw=10; kplt.nd=3;
	loop(pt,
		put rc(theta,pt,"x"),rc(theta,pt,"y")/;
	);
	put //;
);
$onput
EOD
set style line 1 lc rgb 'black' lw 2 dt 1
set style line 2 lc rgb 'red'   lw 2 dt 1
set style line 3 lc rgb 'blue'  lw 2 dt 1
set style line 4 lc rgb 'forest-green' lw 2 dt 1
set style line 5 lc rgb 'black' lw 2 dt 2      
set style line 6 lc rgb 'red'   lw 2 dt 2      
set style line 7 lc rgb 'blue'  lw 2 dt 2      
set style line 8 lc rgb 'forest-green' lw 2 dt 2
set style line 9 lc rgb 'black' lw 2 dt 3      
set style line 10 lc rgb 'red'   lw 2 dt 3      		  
array lstyle[10] = [1,2,3,4,5,6,7,8,9,10]

set terminal pngcairo size 1209,764 enhanced font 'Verdana,8'
set output 'continuous.png'
plot for [i=1:10] '$routes' index i using 1:2 with lines ls lstyle[i] notitle

$offput

putclose;

*	Generate the plot:

execute 'gnuplot continuous.plt -';
