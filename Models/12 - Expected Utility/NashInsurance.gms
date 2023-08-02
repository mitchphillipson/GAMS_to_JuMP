$title	SubOptimal Insurance

*	Declare some parameters with assigned values.

parameter	pi	True probability of bad outcome /0.01/
		L	Loss with a bad outcome /0.5/,
		sigma	Elasticity /2/;

*	GAMS is not case sensitivity, but we following the 
*	convention that parameters (exogenous inputs) are
*	in lower case and variables (endogenous outputs) are
*	written in upper case (except for "L").

*	Declare a parameter whose value will be assigned:

parameter	rho	Risk exponent;

rho = 1 - 1/sigma;

*	Declare variables to be determined in the optimization:

variables	EU	Expected utility,
		C_G	Consumption on a good day,
		C_B	Consumption on a bad day,
		GAMMA	Premium for coverage
		K	Coverage;

*	Declare some equations:

equations	eudef, budget_g, budget_b;

*	Represent the utility function as a macro:

$macro U(C)	((C)**rho/rho)

*	Expected utility:

eudef..		EU =e= (1-pi) * U(C_G) + pi * U(C_B);

*	Consumption in the good state:

budget_G..	C_G =e= 1 - GAMMA * K;

*	Consumption in the bad state:

budget_B..	C_B =e= 1 - L + (1-GAMMA) * K;

model insurance /all/;

C_G.L = 1;
C_B.L = 1;
K.L = 1;


*	Declare a macro to compute marginal utility:

$macro MU(c)	((C)**(rho-1))

*	This is a complementarity constraint -- if the marginal cost
*	exceeds the marginal benefit, then K must be zero.

*		Marginal cost of the insurance	        =g= Marginal benefit

equation coverage;

*	Cost of insurance is the premium (GAMMA) which must be paid
*	in each state of the world.  Benefit is the expected value of the
*	payment made in the bad state.

coverage..	GAMMA * ((1-pi) * MU(C_G) + pi*MU(C_B)) =g= pi * MU(C_B);

*	Declare the model as an equilibrium problem corresponding to the first
*	order conditions of the nonlinear programming model:

model equilibrium /eudef.EU, budget_g.C_G, budget_B.C_B, coverage.K/;

*GAMMA.FX = pi;

*solve equilibrium using mcp;


*$exit

variable	PROFIT		Profit of one firm;

parameter	kother		Insurance provided by other firms /0/,
			dev			Deviation from Nash equilibrium
			iterlog		Iteration log for diagonalization;

equation	profitdef;

profitdef..	PROFIT =e= (GAMMA - PI) * (K - kother);

model nash /all/;

set	n	Number of symmetric insurance companies /1*5/,
	iter	Nash iterations for diagonalization /iter1*iter25/;

GAMMA.LO = 0;
GAMMA.UP = +inf;
C_G.LO = 1e-5;
C_B.LO = 1e-5;

*	Place a lower bound on expected utility equal to the value
*	were there no insurance.

EU.LO = (1-pi) * U(1) + pi * U(1-L);

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

*option iterlog:3:2:1;
display iterlog;

execute_unload "iterlog.gdx", iterlog;