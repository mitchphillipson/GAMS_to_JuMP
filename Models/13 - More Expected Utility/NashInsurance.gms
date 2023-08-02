$title	Nash Equilibrium Insurance

*	Declare some parameters with assigned values (inputs to the analysis)

parameter	pi		Consumer subjective estimate of accident /0.01/,
		L		Loss with a bad outcome /0.5/,
		sigma		Degree of relative risk aversion /0.9/,
		kvalue		Lagged value of insurance provision /0/,
		nfirm		Number of firms /1/;

*	GAMS is not case sensitivity, but we following the 
*	convention that parameters (exogenous inputs) are
*	in lower case and variables (endogenous outputs) are
*	written in upper case (except for "L").

*	-------------------------------------------------------------------------------------------------
*	Declare two parameters whose values will be assigned:

parameter	rho		Primal exponent corresponding to SIGMA ,
		gammamax	Maximum value for GAMMA;

*	The primal-form risk exponent is related 
rho = 1 - 1/sigma;

*	The no-insurance outcome determines the maximum amount which can be
*	charged for coverage:

gammamax = pi * (1-L)**(-1/sigma) / ( pi * (1-L)**(-1/sigma) + 1 - pi )
display gammamax;

*	-------------------------------------------------------------------------------------------------
*	Declare variables and equations to be determined in the optimization:

variables	P_G	Price index for consumption on a good day,
		P_B	Price index for consumption on a bad day,
		C_G	Consumption in the good day,
		C_B	Consumption in the bad day,
		P_C	Consumption price index,
		M	Income,
		K	Coverage
		GAMMA	Premium for coverage
		PROFIT	Firm profit;

nonnegative variables P_G, P_B, P_C, GAMMA, K;

*	Declare some equations:

equations	income, P_Cdef, market_g,market_b, demand_k, profitdef, C_Gdef, C_Bdef;

*	-------------------------------------------------------------------------------------------------
*	Define equations:

income..	M =E= P_G + P_B*(1-L);

P_Cdef..	P_C =e= ( (1-pi) * (P_G/(1-pi))**(1-sigma) + pi * (P_B/pi)**(1-sigma) )**(1/(1-sigma));

*	Consumption in the good state:

C_Gdef..	C_G =e= (M/P_C * (P_C*(1-pi)/P_G)**sigma);

C_Bdef..	C_B =e=	(M/P_C * (P_C*pi/P_B)**sigma);

market_G..	1 =e=  C_G + GAMMA*K;

market_B..	1 - L + K =e= C_B + GAMMA*K ;

demand_K..	GAMMA*(P_G+P_B) =e= P_B;

*	Profit of an individual firm:

profitdef..	PROFIT =e= (GAMMA - pi) * (K - kvalue*(nfirm-1)/nfirm);

*	Declare the model as an equilbrium problem -- this can only be solved with 
*	a fixed value of GAMMA:

model insurance /income.M, P_Cdef.P_C, 
		market_g.P_G, market_b.P_B, 
		C_Gdef, C_Bdef, 
		demand_K.K, profitdef.PROFIT/;

*	Declare the model without associating equations and variables so that it can be solved
*	as an optimization problem with an endogenous value of GAMMA:

model nash /all/;

*	-------------------------------------------------------------------------------------------------
*	Assign initial values:

P_C.L = 1;
P_G.L = 1-pi;
P_B.L = pi;
M.L = 1-pi + pi*(1-L);

P_C.LO = 1e-5;
P_G.LO = 1e-5;
P_B.LO = 1e-5;

P_G.FX = 1-pi; 

C_G.L = (1-pi) + pi*(1-L);
C_B.L = (1-pi) + pi*(1-L);

GAMMA.FX = pi;

K.L = (C_B.L-(1-L))/(1-pi);
GAMMA.L = pi;

*	-------------------------------------------------------------------------------------------------
*	Replicate the competitive equilibrium:

insurance.iterlim = 0;
solve insurance using mcp;
abort$round(insurance.objval,6) "Benchmark replication fails.";


*	-------------------------------------------------------------------------------------------------
*	Declare some parameters for calculation of the Nash equilibria:


parameter	dev		Deviation from Nash equilibrium,
		iterlog		Iteration log for diagonalization,
		equil		Equilibrium values;

*	Save values from the competitive equilibrium:

equil("Competitive","K") = K.L;
equil("Competitive","K/L") = K.L/L;
equil("Competitive","GAMMA") = GAMMA.L;
equil("Competitive","PROFIT") = (GAMMA.L - pi) * K.L;
equil("Competitive","P_B/pi") = P_B.L/pi;
$ondotl
equil("Competitive","C_G") = C_G;
equil("Competitive","C_B") = C_B;


*	-------------------------------------------------------------------------------------------------
*	We solve the Nash equilibrium using projection

set	n	Number of symmetric insurance companies /1*15/,
	iter	Nash iterations for diagonalization /1*25/;

GAMMA.LO = 0;
GAMMA.UP = gammamax;
nash.solvelink = 5;
loop(n,
	dev = 1;
	nfirm = n.val;
	kvalue = K.L;
	loop(iter$round(dev,5),
	  GAMMA.L = gammamax/2;
	  solve nash using nlp maximizing PROFIT;

	  iterlog(n,iter,"dev") = dev;
	  iterlog(n,iter,"K") = K.L;
	  iterlog(n,iter,"GAMMA") = GAMMA.L;
	  iterlog(n,iter,"PROFIT") = (GAMMA.L - pi) * K.L;
	  iterlog(n,iter,"P_B/pi") = P_B.L /pi;
	  dev = abs(kvalue - K.L);
	  kvalue = K.L;

*	Save the equilrium values to report for comparison:
	  equil(n,"dev") = dev;
	  equil(n,"K") = K.L;
	  equil(n,"K/L") = K.L/L;
	  equil(n,"GAMMA") = GAMMA.L;
	  equil(n,"PROFIT") = (GAMMA.L - pi) * K.L;
	  equil(n,"P_B/pi") = P_B.L/pi;
$ondotl
	  equil(n,"C_G") = C_G;
	  equil(n,"C_B") = C_B;
	);
);

option iterlog:3:2:1;
display iterlog,equil;

$exit

*	On Windows the data can be written directly to Excel:
execute_unload 'NashInsurance.gdx',iterlog, equil;
execute 'gdxxrw i=NashInsurance.gdx o=NashInsurance.xlsx par=iterlog rng=PivotData!a2 cdim=0 intastext=n par=equil rng=Equil!a2 cdim=0 intastext=n';
