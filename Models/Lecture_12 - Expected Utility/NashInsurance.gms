$title	SubOptimal Insurance

*	Declare some parameters with assigned values.

parameter	pi	True probability of bad outcome /0.01/
		L	Loss with a bad outcome /0.5/,
		sigma	Elasticity /2/
		rho	Risk exponent,
		kother		Insurance provided by other firms /0/,
		dev		Deviation from Nash equilibrium
		iterlog		Iteration log for diagonalization;

rho = 1 - 1/sigma;

*	Declare variables to be determined in the optimization:

variables	EU	Expected utility,
		C_G	Consumption on a good day,
		C_B	Consumption on a bad day,
		GAMMA	Premium for coverage
		K	Coverage
		PROFIT		Profit of one firm;

GAMMA.LO = 0;
GAMMA.UP = +inf;
C_G.LO = 1e-5;
C_B.LO = 1e-5;

C_G.L = 1;
C_B.L = 1;
K.L = 1;

*	Place a lower bound on expected utility equal to the value
*	were there no insurance.


* Bug in Tom's code here. Macros are a straight string replace,
* this means you need (C), not just C. Otherwise, U(1-L) gets
* replaced with 1-L^rho/rho, which is not what we want. 
$macro U(C)	((C)**rho/rho) 


EU.LO = (1-pi) * U(1) + pi * U(1-L);


equation	profitdef;

profitdef..	PROFIT =e= (GAMMA - PI) * (K - kother);


equations	eudef, budget_g, budget_b;


eudef..		EU =e= (1-pi) * U(C_G) + pi * U(C_B);

*	Consumption in the good state:

budget_G..	C_G =e= 1 - GAMMA * K;

*	Consumption in the bad state:

budget_B..	C_B =e= 1 - L + (1-GAMMA) * K;


$macro MU(c)	(C**(rho-1))

equation coverage;

coverage..	GAMMA * ((1-pi) * MU(C_G) + pi*MU(C_B)) =g= pi * MU(C_B);


model nash /all/;



set	n	Number of symmetric insurance companies /1*5/,
	iter	Nash iterations for diagonalization /iter1*iter25/;


nash.solvelink = 5;
loop(n,
	dev = 1;
	kother = K.L * (n.val-1)/n.val;
	loop(iter$round(dev,4),
	  solve nash using nlp maximizing PROFIT;
	  iterlog(n,iter,"dev") = dev;
	  iterlog(n,iter,"K") = K.L;
	  iterlog(n,iter,"GAMMA") = GAMMA.L;
	  iterlog(n,iter,"PROFIT") = PROFIT.L;
	  iterlog(n,iter,"C_G") = C_G.L;
	  iterlog(n,iter,"C_B") = C_B.L;
	  dev = abs(kother - K.L * (n.val-1)/n.val);
	  kother = K.L * (n.val-1)/n.val;
	);
);

option iterlog:3:2:1;
display iterlog;
