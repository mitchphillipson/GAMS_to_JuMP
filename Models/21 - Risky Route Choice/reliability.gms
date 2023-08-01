$title	Optimal Route Trades off Travel Time and Reliability

$setglobal n 20

set	r	Rows and columns /0*%n%/;

alias (r, ir,ic,jr,jc);

$setglobal i ir,ic
$setglobal j jr,jc

*	Introduce a short-hand reference to nodes in the network:

set	i(%i%)	Nodes in the network;
i(%i%) = yes;

*	Use i, j and k to reference nodes:

alias (i,j,k);

parameter	id(%i%)	Node identifier (integer),
		nid	ID counter /0/;
loop(i,
	nid = nid + 1;
	id(i) = nid + 1;);

set	a(%i%,%j%)	Arcs include adjacent nodes;

*	Add arcs which run north-south, east-west or diagonal:

a(%i%,%j%) = yes$(max(abs(ir.val-jr.val),abs(ic.val-jc.val)) = 1);

*	Add arcs on the adjacent row or column which are 
*	two blocks away from node %i%:

$if not set diag $set diag 2
$setglobal diag %diag%

if (%diag%>1,
  a(%i%,%j%)$(not a(%i%,%j%))
	= yes$( (max(abs(ir.val-jr.val),abs(ic.val-jc.val)) = 2) and
		(min(abs(ir.val-jr.val),abs(ic.val-jc.val)) = 1) ););
if (%diag%>2,
  a(%i%,%j%)$(not a(%i%,%j%))
	= yes$( (max(abs(ir.val-jr.val),abs(ic.val-jc.val)) = 3) and
		(min(abs(ir.val-jr.val),abs(ic.val-jc.val)) = 1)););
if (%diag%>3,
  a(%i%,%j%)$(not a(%i%,%j%))
	= yes$( (max(abs(ir.val-jr.val),abs(ic.val-jc.val)) = 4) and
		(min(abs(ir.val-jr.val),abs(ic.val-jc.val)) = 1)););

a(%i%,%i%) = no;
option a:0:0:1;
display a;

parameter	dist	Distance between two nodes;
dist(a(%i%,%j%)) = sqrt( sqr(ir.val-jr.val)+sqr(ic.val-jc.val))/%n%;

*	Start in the southeast and traverse to the northwest:
  
set	start(%i%)	Starting point /%n%.0/
	end(%i%)	Ending point   /0.%n%/;

set	shortest(%i%,%j%)	Arcs on the diagonal path from s to t;
shortest(%i%,ir-1,ic+1) = yes$(ir.val+ic.val=%n%);
option shortest:0:0:1;
display shortest;

set	dp		Danger points /dp1*dp3/
	dploc(dp,%i%)	Danger points /dp1.5.6, dp2.6.18, dp3.18.6/;

parameter	nrisk(%i%)	"Risk index at each node";
nrisk(%i%) = sum(dploc(dp,%j%), 1/(0.01+sqr(ir.val-jr.val)+sqr(ic.val-jc.val)));
display nrisk;

parameter	p	"Probability of success on each arc",
		pref	Reference probability of failure /0.2/,
		rref	Radius of reference risk level
		pscale	Scale factor;

p(a(i,j)) = 0.5 * (nrisk(i)+nrisk(j));

*	prod (pscale p)**d = pref

*	sum(d * (log(pscale) + log(p)) = log(pref)

*	log(pscale) = (log(pref)-sum(d*log(p)))/sum(d)

pscale = exp((log(pref) - 
	sum(shortest(i,j), dist(i,j)*log(p(i,j))))/sum(shortest(i,j),dist(i,j)));
p(a) = p(a) * pscale;
pref = prod(shortest(i,j), p(i,j)**dist(i,j));

*	pref = sum(d) * pscale/(0.01+sqr(rref))

*	rref = sqrt(pref/(pscale*sum(d)) - 0.01)

rref = sqrt(pref/(pscale*sum(shortest,dist(shortest)))-0.01);
display rref;

display pref;

parameter	rmax	Maximum tolerable risk;

variable	OBJ	Objective function;

nonnegative
variables	X		Route choice;

equations	conservation, risktol, objdef;

*	Define the objective as a weighted sum of travel distance and probability 
*	of reaching the target:

objdef..		OBJ =e= sum(a, X(a) * dist(a));

risktol..		log(rmax) =g= sum(a, X(a) * dist(a) * log(p(a)));

conservation(k)..	sum(a(i,k), X(a)) + 1$start(k) =e= sum(a(k,j),X(a)) + 1$end(k);

model routechoice /conservation, risktol, objdef/;

rmax = pref;
solve routechoice using lp minimizing obj;

variable	OBJD	Dual Objective,
		PRISK	Shadow price on the risk constraint,
		T(%i%)	Generalize time cost to reach finish line from node k;

equations	dualobj	Defines OBJDUAL
		dualcon	Defines dual constraints;

dualobj..		OBJD =e= sum(end, T(end)) - sum(start,T(start)) + log(rmax)*PRISK;

dualcon(a(i,j))..	T(j) + dist(a) - PRISK * dist(a)*log(p(a)) =g= T(i);

model dual /dualobj, dualcon/;

PRISK.L = risktol.M;
T.L(k) = conservation.M(k);
OBJD.L = OBJ.L;
solve dual using lp maximizing OBJD;

set	riskval		Set of values of rmax (x10) /1*10/;

set		rc(riskval,%i%,%j%)	Route choice;
loop(riskval,
	rmax = riskval.val/100;
	solve routechoice using lp minimizing obj;
	rc(riskval,a(i,j))$round(X.L(a),4) = yes;
);
option rc:0:0:1;
display rc;

$set i ir,ic
$set j jr,jc

*	Generate the .png files using GNUPLOT:

file kplt /'figure.plt'/; put kplt; kplt.lw=0;

put 'reset'/;
put "#set terminal pngcairo size 1209,764 enhanced font 'Verdana,8'"/;
put "#set output 'figure_%diag%.png'"/;
put 'set xrange [0:%n%]'/;
put 'set yrange [0:%n%]'/;
put 'set object 1 rect from 0,0 to  %n%,%n%'/;
loop(dploc(dp,%i%),
	put 'set object circle at ',ir.tl,',',ic.tl,'  size ',(rref*%n%),' fc rgb "red"'/;);

*	Show the underlying network:

loop(a(i(%i%),j(%j%))$(id(i)<id(j)),
  put 'set arrow from ',ir.tl,',',ic.tl,' to ',jr.tl,',',jc.tl,
	' nohead ls 0'/;);
put /;

PUT '$network <<EOD'/;
loop(a(i(%i%),j(%j%))$(id(i)<id(j)),
  put ir.tl,',',ic.tl,',',jr.tl,',',jc.tl/;);
put 'EOD'//;

loop(riskval,
	put '$risk',riskval.tl,' <<EOD'/;
	loop(rc(riskval,%i%,%j%),
	  put	ir.tl,',',ic.tl/
		jr.tl,',',jc.tl//;);
	put 'EOD'//; );

*	Draw arrows along the route:

loop(riskval,
  loop(rc(riskval,%i%,%j%),
    put 'set arrow from ',ir.tl,',',ic.tl,' to ',jr.tl,',',jc.tl,
	  ' head ls ',riskval.tl/;));
*	This plot includes only arrows, so we enter a vacuous function and remove the title:

put 'plot NaN notitle'/;

putclose;

*	Generate the plot:

execute '=wgnuplot -persist figure.plt -';
