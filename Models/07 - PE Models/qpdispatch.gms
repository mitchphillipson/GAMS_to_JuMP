$title	Simple Model of Electricity Market Dispatch 

*	1. Define sets.

set	s	Load segments (peak to base load) /s1*s9/,
	j	Generating units  /Nuclear, Coal, Gas, Diesel/
	i	Demand categories /
		rsd	Residential, 
		com	Commercial, 
		ind	Industrial /;

*	2. Read data.

table	load(s,*)	Electricity load

*			Demand Shares
*			-------------------		
      	qref	hours	rsd	com	ind
s1    	1.00      310	0.7	0.2	0.1
s2    	0.79      750	0.5	0.2	0.3
s3    	0.59     1020	0.4	0.3	0.3
s4    	0.50     2460	0.4	0.3	0.3
s5    	0.44     1910	0.3	0.3	0.4
s6    	0.40      580	0.2	0.3	0.5
s7    	0.35      580	0.1	0.3	0.6
s8    	0.28      600	0.1	0.3	0.6
s9    	0.23      540	0.1	0.3	0.6;

table	supply(j,*)	Electricity supply technology
                   mc   cap
Nuclear             4   0.30
Coal                6   0.30
Gas                 7   0.20
Diesel              9   0.30;

parameter	mc(j)	Marginal cost of dispatch,
		cap(j)	Capacity constraint,
		qref(s)	Reference demand
		h(s)	Hours;

h(s) = load(s,"hours");
mc(j) = supply(j,"mc");
cap(j) = supply(j,"cap");
qref(s) = load(s,"qref");

*	3. Compute reference prices for each of the segments using 
*	a linear program.

variables	TOTCOST	Objective function (dispatch cost);

nonnegative
variables	Y(j,s)	Dispatch;

equations	costdef, demand;

costdef..	TOTCOST =e= sum((j,s), mc(j)*Y(j,s));

demand(s)..	sum(j, Y(j,s)) =e= qref(s);

model mincost /costdef, demand/;

Y.UP(j,s) = supply(j,"cap");

solve mincost minimizing TOTCOST using LP;

parameter	pref(s)		Reference price,
		dref(i,s)	Reference demand;

pref(s) = demand.m(s);
dref(i,s) = qref(s) * load(s,i);

*	4. Set up a calibrated equilibrium model:

parameter	epsilon(i)	Elasticity of demand /rsd 0.1,  com 0.2,  ind 0.5/;

nonnegative
variables	D(i,s)	Aggregate demand,
		PI(j,s)	Shadow price on capacity,
		P(s)	Market price;

equations aggdemand, supplydemand, profit, capacity;

aggdemand(i,s)..	D(i,s) =e= dref(i,s) * (1 - epsilon(i)*(P(s)/pref(s)-1));

supplydemand(s)..	sum(j,Y(j,s)) =e= sum(i,D(i,s));

profit(j,s)..		mc(j) + PI(j,s) =g= P(s);

capacity(j,s)..		cap(j) =g= Y(j,s);

model equil /aggdemand.D, supplydemand.P, profit.Y, capacity.PI/;

Y.UP(j,s) = inf;
P.L(s) = pref(s);
D.L(i,s) = dref(i,s);
PI.L(j,s) = -Y.M(j,s);

equil.iterlim = 0;
SOLVE equil USING mcp;

*	5. Set up a nonlinear programming model which represents the
*	equilibrium allocation:

variable	SURPLUS		Sum of consumer and producer surplus
		K(j)		Capacity of technology j;

equation	surplusdef	Defines the surplus;

surplusdef..	SURPLUS =e= sum((i,s), pref(s)*D(i,s) *
				(1 + 1/epsilon(i) * (1 - D(i,s)/(2*dref(i,s)))))
			 -  TOTCOST;

MODEL SAMUELSON /surplusdef, supplydemand, costdef, capacity/;

K.FX(j) = supply(j,"cap");
Y.UP(j,s) = +inf;

SOLVE samuelson USING nlp MAXIMIZING surplus;

parameter	ldc(s,*)	Load-dispatch curve;
ldc(s,"Benchmark") = sum(i,D.L(i,s));

*	6. Solve a counterfactual simulation:

mc("Coal") = 2 * mc("Coal");

SOLVE samuelson USING nlp MAXIMIZING surplus;

ldc(s,"Counterfactual") = sum(i,D.L(i,s));
display ldc;

parameters	x	Cumulative load segments;
x(s) = 0;
loop(s, x(s) = x(s-1) + load(s,"hours"););
display x;


$exit

*	7. Display the results:

$setglobal domain s
$setglobal xvalue x
$libinclude plot ldc

