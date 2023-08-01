$title	Optimal Route Choice: Time and Reliability

*	Find a route from starting point to ending point
*	which minimizes risk and distance.

*	Here we discretize space and solve for the 
*	optimal path using linear programming.

*	Define the number of rows/columns in a square grid:

$set	n		20

*	Nodes are in the grid:

$eval	nnode (%n%+1)*(%n%+1)

set	r		Rows and columns /0*%n%/
	n		Nodes		 /1*%nnode%/;

alias (r,ir,ic,jr,jc), (n,i,j,k);

*	Start at the lower right point and end at the upper left:

set	sp(ir,ic)	Starting point /%n%.0/,
	ep(ir,ic)	Ending point   /0.%n%/;

parameter	nc		Node count /0/
set		loc(n,ir,ic)	Node locations on the grid;
loop(n$(n.val=1),
	loop((ir,ic),
	  loc(n+nc,ir,ic) = yes;
	  nc = nc + 1;
));



singleton sets
	start(i)	Starting point
	end(j)		Ending point;

loop(sp, start(i) = yes$loc(i,sp););
loop(ep, end(j)   = yes$loc(j,ep););
display start, end;


set	a(i,j)	Arcs include adjacent nodes;

*	Add arcs which run north-south, east-west or diagonal:

$if not set diag $set diag 4

alias (loc,iloc,jloc);
loop((iloc(i,ir,ic),jloc(j,jr,jc)),
	a(i,j)$((max(abs(ir.val-jr.val),abs(ic.val-jc.val)) <= %diag%) and
		(min(abs(ir.val-jr.val),abs(ic.val-jc.val)) = 1)) = yes;
	a(i,j)$(max(abs(ir.val-jr.val),abs(ic.val-jc.val)) <= 1) = yes;
);
a(i,i) = no;

option a:0:0:1;
display a;

parameter	dist(i,j)	Distance between two nodes,
		w		Width of the grid /%n%/;

w = card(r)-1;
loop((a(i,j),iloc(i,ir,ic),jloc(j,jr,jc)),
	dist(i,j) = sqrt( sqr(ir.val-jr.val)+sqr(ic.val-jc.val))/w;
);

set	h(*)		Hazards
	hloc(h<,ir,ic)	Hazard locations /h0.0.0/;

parameter	pn	"Probability of success at each node -- highest at %n%,%n%";
loop(iloc(i,ir,ic),
	pn(i) = prod(hloc(h,jr,jc), sqrt(sqr(ir.val-jr.val)+sqr(ic.val-jc.val)))/sqrt(2*w*w););
display pn;

parameter	p(i,j)	"Probability of success on each arc"; 
p(a(i,j)) = 0.5  * (pn(i)+pn(j));

parameter	lamda	Weight on speed (1=shortest path) /1/;

*	Start in the southeast and traverse to the northwest:
  
variable	OBJ		Objective function;

display dist;
$exit
nonnegative
variables	X(i,j)		Route choice;

equations	conservation, objdef;

*	Define the objective as a weighted sum of travel distance and probability 
*	of reading the target:

objdef..		OBJ =e= sum(a, X(a) * dist(a)*(lamda - (1-lamda)*log(p(a))));

conservation(k)..	sum(a(i,k), X(a)) + 1$start(k) =e= sum(a(k,j),X(a)) + 1$end(k);

model routechoice /conservation, objdef/;

lamda = 0.5;
solve routechoice using lp minimizing obj;

$exit

parameter	prob	Probability of completion;
prob = prod(a$X.L(a), p(a)**dist(a));
display prob;


variable	OBJD		Dual Objective,
		T(i)		Generalize time cost to reach finish line from node k;

equations	dualobj		Defines OBJDUAL
		dualcon		Defines dual constraints;

dualobj..		OBJD =e= sum(end, T(end)) - sum(start,T(start));

dualcon(a(i,j))..	T(j) + dist(a)*(lamda - (1-lamda)*log(p(a))) =g= T(i);

model dual /dualobj, dualcon/;

T.L(k) = conservation.M(k);
OBJD.L = OBJ.L;
solve dual using lp maximizing OBJD;

set	theta			Angle which determines lamda (degrees) /0,10,20,30,40,50,60,70,80,90/,
	xy			Directions /x,y/
	pt			Points in sequence /1*100/,
	pts(theta,pt)	Points required for this path
	point(i)		Location for unwinding route;

parameter
	npt			Point counter,
	rc(theta,pt,xy)	Route choices;

*	Always include the starting point:
rc(theta,"1","x") = %n%;
rc(theta,"1","y") = 0;
pts(theta,"1") = yes;
loop((theta,pt)$sameas(pt,"1"),
	lamda = cos(pi*theta.val/180)/(cos(pi*theta.val/180)+sin(pi*theta.val/180));
	solve routechoice using lp minimizing obj;

*	Move the arcs to a list of points:

	npt = 1;
	point(i) = yes$start(i);
	while(not point(end),
	  loop((a(i,j),loc(j,ir,ic))$(round(X.L(i,j),5) and point(i)),
	    pts(theta,pt+npt) = yes;	    
	    rc(theta,pt+npt,"x") = ir.val;
	    rc(theta,pt+npt,"y") = ic.val;
	    npt = npt + 1;
	    point(i) = no; 
	    point(j) = yes;
	  );
	);
);
option rc:0:0:1;
display rc;

file kplt /'discrete.plt'/; put kplt; kplt.lw=0;

$onput
reset
set size square
unset xtics
unset ytics
unset raxis
unset rtics
unset border
$offput

put	'set xrange [0:%n%]'/
	'set yrange [0:%n%]'/
	'set object 1 rect from 0,0 to  %n%,%n%'//
	'$routes <<EOD'/;

loop(theta,
	kplt.nw=0; kplt.nd=1;
	put '# theta=',(theta.val)/;
	kplt.nw=10; kplt.nd=3;
	loop(pts(theta,pt),
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
array lstyle[11] = [1,2,3,4,5,6,7,8,9,10]

set terminal pngcairo size 1209,764 enhanced font 'Verdana,8'
set output 'discrete.png'

plot for [i=1:10] '$routes' index i using 1:2 with lines ls lstyle[i] notitle

$offput

putclose;

*	Generate the plot:

execute 'gnuplot discrete.plt -';
