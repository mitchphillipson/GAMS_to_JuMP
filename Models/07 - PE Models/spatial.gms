$title	Calibration and Ad-Valorem Taxes in a Spatial Price Equilibrium Model

*	1. Generate a random dataset.

set     i       Supply regions /s1*s5/
	j	Demand regions /d1*d10/;

parameter       dref(j)   Reference demand
                slim(i)   Maximum supply
                mu(i)     Marginal cost of production,
                c(i,j)    Transport cost;

*	Randomly distributed demand and supply:

c(i,j) = uniform(0,1);
dref(j)  = round(uniform(1,100));
slim(i)  = 2*card(j)/card(i)*round(uniform(1,100));
mu(i)  = uniform(0.5,1.5);

*	2. Formulate a linear program to calibrate benchmark trade flows
*	and market prices.

*	Begin by setting up a linear programming model which determines 
*	reference

nonnegative variables   
		X(i,j)	Shipments from i to j,
		S(i)	Supply in region i,
		D(j)	Demand in region j;

free variable   TOTCOST         Objective function;

equations       objdef, supply, demand;

objdef..        TOTCOST =e= sum((i,j), c(i,j) * X(i,j)) + sum(i, mu(i)*S(i));

*       Orient both equations as >= so that the Lagrange multipliers
*       are non-negative:

supply(i)..     S(i) =g= sum(j, X(i,j));

demand(j)..     sum(i, X(i,j)) =g= D(j);

model bmk /all/;

*       Fix demand and place an upper bound on supply in order
*       that the marginal cost of supply is included in the 
*       shadow prices at the equilibrium point:

S.UP(i) = slim(i);  
D.FX(j) = dref(j);

*	Solve a linear program to evaluate supplies and reference prices:

solve bmk using LP MINIMIZING TOTCOST;

abort$(bmk.solvestat<>1) "bmk model error.";
abort$(bmk.modelstat>2) "bmk model is not solved.";


*	3. Define a parameter for storing results.

parameter	qtylog	Log report of values;

qtylog(i,"bmk","S") = S.L(i);
qtylog(j,"bmk","D") = D.L(j);
qtylog(i,"bmk","pS") = supply.M(i);
qtylog(j,"bmk","pD") = demand.M(j);

*	4. Set up the demand function with the reference price and an
*	elasticity of demand.

parameter	pref(j)         Reference demand price
                epsilon(j)      Demand elasticity at the reference point;

pref(j) = demand.m(j); 
epsilon(j) = uniform(0.5, 2);

*	Demand is no longer fixed:

D.LO(j) = 0;  
D.UP(j) = +inf;

*	5. Formulate a model with price elastic demand by including a
*	quadratic form in the objective function corresponding to 
*	consumer surplus.  We set this up as a minimization problem
*	so that Lagrange multipliers for supply(i) and demand(j) are
*	non-negative.

equation        csurplus        Social surplus with horizontal supply curves (Cs);

csurplus..      TOTCOST =e= sum((i,j), c(i,j) * X(i,j)) + sum(i, mu(i)*S(i))

        - sum(j,  pref(j) * D(j) * (1 + (1-0.5*D(j)/dref(j)) / epsilon(j)));

model elast /supply, demand, csurplus/;

*       Now that we have included the consumer surplus in the objective
*	we no longer need to have fixed demand.  To install the 
*	solution, we need to assign the objective value.  $ondotl
*	means that we don't need to include .L values on the RHS
*	of this assignment.

$ondotl
TOTCOST.L = sum((i,j), c(i,j) * X(i,j)) + sum(i, mu(i)*S(i))
        - sum(j,  pref(j) * D(j) * (1 + (1-0.5*D(j)/dref(j)) / epsilon(j)));

*	Assign shadow prices for the objective function and for the 
*	bounds on D.

csurplus.m = 1;
D.M(j) = 0;

*	6. First use Examiner to verify that we have a calibrated model,
*	then solve the model using CONOPT:

elast.optfile = 1;
option qcp=examiner;
$echo returnGamsPoint 1   > examiner.opt
solve elast using QCP minimizing TOTCOST;

*	N.B. If we use CPLEX here, the calibration fails.  This has been
*	reported.

option qcp=conopt;
solve elast using QCP minimizing TOTCOST;

*	Abort if the model does not calibrate.

parameter
	dlog	Log of demand values,
	dev	Deviation in demand values;

dlog(j,"dref") = dref(j);
dlog(j,"D.L") = D.L(j);
dlog(j,"dev") = abs(dref(j)-D.L(j));
display dlog;

dev = sum(j, abs(D.L(j)-dref(j)));
if (round(dev,5),

*	If the model fails to solve, double check the
*	solution using examiner.

	option qcp=examiner
	solve elast using QCP minimizing TOTCOST;
);
abort$round(dev,3) "Consumer surplus calibration fails!", dev;

qtylog(i,"elast","S") = S.L(i);
qtylog(j,"elast","D") = D.L(j);
qtylog(i,"elast","pS") = supply.M(i);
qtylog(j,"elast","pD") = demand.M(j);


*	7. Incorporate producer surplus in the 
*	objective function with random price elasticities of 
*	supply and calibrated marginal cost:

parameter	eta(i)		Price elasticity of supply,
                sref(i)         Reference supply;

eta(i) = uniform(0.5, 2);
sref(i) = S.L(i); 

equation        ssurplus        Social surplus with price elastic supply;

ssurplus..      TOTCOST =e= sum((i,j), c(i,j) * X(i,j)) 
        - sum(j,  pref(j)  * D(j) * (1 + (1-0.5*D(j)/dref(j)) / epsilon(j)))
        + sum(i,  mu(i)    * S(i) * (1 + (0.5*S(i)/sref(i)-1)/eta(i)));

model mkt /supply, demand, ssurplus/;

*	Same steps involved here -- assign marginals and then use 
*	examiner to verify that we have a solution:

SSURPLUS.M = 1;
TOTCOST.L = sum((i,j), c(i,j) * X(i,j)) 
        - sum(j,  pref(j) * D(j) * (1 + (1-0.5*D(j)/dref(j)) / epsilon(j)))
        + sum(i,  mu(i) * S(i) * (1 + (0.5*S(i)/sref(i)-1)/eta(i)));

solve mkt using QCP minimizing TOTCOST;

dev = sum(j, abs(D.L(j)-dref(j)));
abort$round(dev,3) "Market equilibrium model calibration fails!", dev;


qtylog(i,"mkt","S") = S.L(i);
qtylog(j,"mkt","D") = D.L(j);
qtylog(i,"mkt","pS") = supply.M(i);
qtylog(j,"mkt","pD") = demand.M(j);

option qtylog:3:1:2;
display qtylog;

*	8. Set up a model in which we can apply ad-valorem taxes on 
*	trade.  We do this by adding a specific tax term to the 
*	transport cost.  The specific cost equals the ad-valorem tax
*	rate times the market price.

parameter	p(i)		Price used to convert ad-valorem tax to a specific tax,
		tau(i,j)	Tax rate on exports from i to j;

p(i) = supply.m(i);
tau(i,j) = 0;

equation	taxsurplus	Include taxes on bilateral trade in the objective;

taxsurplus..      TOTCOST =e= sum((i,j), (c(i,j) + p(i)*tau(i,j)) * X(i,j)) 
        - sum(j,  pref(j)  * D(j) * (1 + (1-0.5*D(j)/dref(j)) / epsilon(j)))
        + sum(i,  mu(i)    * S(i) * (1 + (0.5*S(i)/sref(i)-1)/eta(i)));

model mkttax /supply, demand, taxsurplus/;

*	Apply a 25% ad-valorem tax on trade between regions:

tau(i,j) = 0.25;

set	iter	iterations /iter1*iter8/;

parameter	plog(*,iter)	Iteration log;

dev = 1;

*	Run iterations until the price deviation is less
*	than 0.001:

loop(iter$round(dev,3),

	p(i) = supply.m(i);

	plog(i,iter) = p(i);

	solve mkttax using QCP minimizing TOTCOST;

*	Calculate price changes:

	dev = sum(i,abs(p(i)-supply.m(i)));

	plog("dev",iter) = dev;
);

display plog;
