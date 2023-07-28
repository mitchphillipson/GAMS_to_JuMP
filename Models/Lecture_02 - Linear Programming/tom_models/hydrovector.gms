$title	Hydro Power Planning Model -- Vector Format

set
  r	Reserviors /A, B/,
  t	All time periods /march, april, may/
  tf(t)	First period /march/,
  tp(t)	Time periods in planning horizon /march, april/,
  i	Types of power /high, low/;

parameter price(i) prices by power types   /high 5.0, low 3.5 /;

Table data summary of relevant data
                A       B
res_cap         2000    1500
minimum         1200    800
march           200     40
april           130     15
level           1900    850
convrate        400     200
pow_cap         60000   35000;

parameter	htr(r)	Heat rate;
htr(r) = data("convrate",r);

FREE
VARIABLE
        REVENUE Aggregate earnings on sales in March and April;

NONNEGATIVE
VARIABLES
        P(r,t)  Power production
        S(r,t)  Water spilled
        L(r,t)  Reservior level (start of month)
	E(i,tp)	Electricity generation;

*	Specify Simple Bounds on Decision Variables

*	Power can be sold at CHF 5 per MWH up to 50,000 MWH per month: 

E.UP("high",tp) = 50000;

*	Storage capacity:

L.UP(r,t) = data("res_cap",r);

*	Minimum reservior levels:

L.LO(r,t) = data("minimum",r);

*	Levels at the beginning of March are specified exogenously:

L.FX(r,tf) = data("level",r);

*	Power plant capacities:

P.UP(r,tp) = data("pow_cap",r);

equations	revenuedef, sales, level;

*	Unlimited quantities of excess power can be sold for
*	CHF 3.50 per MWH.

*	Revenue therefore equals the sum of high and low price
*	revenue:

revenuedef..	REVENUE =E= sum((i,tp), price(i)*E(i,tp));

*	Total sales in March and April equal generation:

sales(tp)..	sum(i, E(i,tp)) =e= sum(r, P(r,tp));

*	Level in period t+1 equals level in t plus net inflows:

level(r,t+1)..	L(r,t+1) =e= L(r,t) + data(t,r) + S(r-1,t) + (P(r-1,t)/htr(r-1))$htr(r-1) - S(r,t) - P(r,t)/htr(r);

model hydro /all/;

solve hydro using lp maximizing revenue;
