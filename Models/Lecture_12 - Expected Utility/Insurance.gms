$title	Optimal Insurance

*	Declare some parameters with assigned values.

parameter	pi	True probability of bad outcome /0.01/
		L	Loss with a bad outcome /0.5/,
		gamma	Premium for coverage /0.02/,
		sigma	Elasticity /0.5/;

*	GAMS is not case sensitivity, but we following the 
*	convention that parameters (exogenous inputs) are
*	in lower case and variables (endogenous outputs) are
*	written in upper case (except for "L").

*	Declare a parameter whose value will be assigned:

parameter	rho	Risk exponent;
rho = 1 - 1/sigma;

*	Declare and solve the model:

variables	EU	Expected utility,
		C_G	Consumption on a good day,
		C_B	Consumption on a bad day,
		K	Coverage;

*	Declare some equations:

equations	eudef, budget_g, budget_b;

*	Represent the utility function as a macro:

$macro U(C)	(C**rho/rho)

eudef..		EU =e= (1-pi) * U(C_G) + pi * U(C_B);

budget_G..	C_G =e= 1 - gamma * K;

budget_B..	C_B =e= 1 - L + (1-gamma) * K;

model insurance /all/;

C_G.L = 1; C_B.L = 1; K.L = 1;

solve insurance using nlp maximizing EU;



parameter	solution	Report of model solution for comparison across models;

solution("C_G","Max_EU") = C_G.L;
solution("C_B","Max_EU") = C_B.L;
solution("K","Max_EU") = K.L;

variable	EV	Equivalent variation;

equation	evdef;

evdef..		EV =e= 100 * (( (1-pi) * C_G**rho + pi * C_B**rho )**(1/rho) - 1);

model insurance_ev /budget_G, budget_B, evdef/;

solve insurance_ev using nlp maximizing EV;



solution("C_G","Max_EV") = C_G.L;
solution("C_B","Max_EV") = C_B.L;
solution("K","Max_EV") = K.L;

*	Declare a macro to compute marginal utility:

$macro MU(c)	(C**(rho-1))

*	This is a complementarity constraint -- if the marginal cost
*	exceeds the marginal benefit, then K must be zero.

*		Marginal cost of the insurance	        =g= Marginal benefit

equation coverage;

*	Cost of insurance is the premium (gamma) which must be paid
*	in each state of the world.  Benefit is the expected value of the
*	payment made in the bad state.

coverage..	gamma * ((1-PI) * MU(c_g) + PI*MU(C_B)) =g= PI * MU(C_B);

*	Declare the model as an equilibrium problem corresponding to the first
*	order conditions of the nonlinear programming model:

model equilibrium /eudef.EU, evdef.EV, budget_g.C_G, budget_B.C_B, coverage.K/;

solve equilibrium using mcp;

solution("C_G","Equilibrium") = C_G.L;
solution("C_B","Equilibrium") = C_B.L;
solution("K","Equilibrium") = K.L;

display solution;

$exit


----    100 PARAMETER solution  Report of model solution for comparison across models

         Max_EU      Max_EV  Equilibri~

C_G       0.996       0.996       0.996
C_B       0.701       0.701       0.701
K         0.205       0.205       0.205


