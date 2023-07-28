$title    Gig Economy Model

parameter
	vref            Marginal value of profit at reference profit /2/
	piref           Reference profit calibrating sigma /0.25/
	esub            Elasticity of substitution E vs L /4/
	thetak          Value share capital /0.3/,
	thetaE          Value share of employees /0.3/,
	thetaL          Value share of hired labor /0.2/
	w		Wealth (relative to reference profit) /1/,
	sigma_p		Variance of the output price p /0.5/;

set     s               States of the world /s1*s50/;

parameter
	wage            Wage of gig workers
	p(s)            Output price in state s
	sigma		Dual degree of relative risk aversion
	rho		Primal degree of relative risk aversion
	pi0		Reference owner profit
	gamma		Primal form of esub,
	stochastic	Flag for the stochastic model /1/;

variables
	EU		Expected utility,
	K(s)		Capital stock,
	E(s)		Full time employees,
	KK		Capital stock (stochastic model)
	EE		Employment (stochastic model)
	L(s)		Gig employees

	X(s)		Ending wealth by state,
	PI(s)		Profit (or social surplus) by state,
	Y(s)		Supply = demand by state,
	EL(s)		Labor-employee nest;

equations       eudef, eldef, ydef, pidef, naK, naE, xdef;

eudef..         EU =e= sum(s, (1/card(s)) * X(s)**rho)**(1/rho);

eldef(s)..      EL(s) =e= (thetaE/(thetaE+thetaL) * (E(s)/thetaE)**gamma +
			   thetaL/(thetaE+thetaL) * (L(s)/thetaL)**gamma)**(1/gamma);

ydef(s)..       Y(s) =e= (K(s)/thetak)**thetak * EL(s)**(thetaE+thetaL);

pidef(s)..      PI(s) =e= p(s) * Y(s)  - (K(s) + E(s) + wage*L(s));

xdef(s)..       X(s) =e= (w+PI(s))/(w+pi0);

*    Include non-anticipativity constraints in the stochastic model:

naK(s)$stochastic..    K(s) =e= KK;

naE(s)$stochastic..    E(s) =e= EE;

model gigeconomy /all/;

gamma = 1 - 1/esub;
pi0 = 1 - thetak - thetaE - thetaL;
sigma = (log(1+pi0)-log(1+pi0*piref))/log(vref);
rho = 1 - 1/sigma;

*    Assign initial values close to the reference equilibrium:

PI.L(s)  = pi0;
L.L(s)   = thetaL;
K.L(s)   = thetaK;
E.L(s)   = thetaE;
X.L(s) = 1;
Y.L(s) = 1; 
EL.L(s) = 1;

*    Avoid bad function calls:

X.LO(s)  = 0.001; L.LO(s)=0.01; K.LO(s) = 0.01; E.LO(s) = 0.01;

*	Prices are log-normally distributed:

p(s) = exp(sigma_p*normal(0,1));

display p;



*	Benchmark wage rate:

wage = 1;

stochastic = 1;
solve gigeconomy using nlp maximizing EU;




parameter    results        Summary of results;

*	Generate a reporting "subroutine" which takes a single
*	argument (the scenario identifier):

$onechov    >%gams.scrdir%report.gms
results("%1",s,"Y") = Y.L(s);
results("%1",s,"EL") = EL.L(s);
results("%1",s,"L") = L.L(s)/thetaL;
results("%1",s,"K") = K.L(s)/thetaK;
results("%1",s,"E") = E.L(s)/thetaE;
results("%1",s,"P") = p(s);
results("%1",s,"PI") = PI.L(s);
results("%1",s,"wage") = wage;
$offecho

*$batinclude %gams.scrdir%report ref_stochastic

parameter	rp	Risk premium (%),
		EU_s    Expected utility in stochastic model,
		EU_d    Expected utility in the deterministic model,
		evpi    Expected value of perfect information;


rp = 100 * (sum(s, X.L(s))/card(s) / EU.L - 1);
EU_s = EU.L;


stochastic = 0;
solve gigeconomy using nlp maximizing EU;

*$batinclude %gams.scrdir%report ref_deterministic


EU_d = sum(s,X.L(s))/card(s);
display EU_d, EU_s;

evpi = 100 * (EU_d / EU_s - 1);
display rp, evpi;

display rho;
*option results:3:1:1;
*display results;
