$TITLE  Brewery Profit Maximization

set	j	Products /lager, ale/
	i	Ingredients /
			malt	Malt, 
			yeast	Yeast,
			dehops	German hops, 
			wihops	Wisconsin hops/;

parameter	
	p(j)	Profit by product /lager 120, ale 90/
	s(i)	Supply by ingredient /malt 4800, yeast 1750,
				    dehops 1000, wihops 1750/;

table		a(i,j)	Requirements
			lager	ale
	malt		4	2
	yeast		1	1
	dehops		1	0
	wihops		0	1;

variables       Y(j)	Production levels,
                Z       Profit (maximand);

nonnegative variable   Y;

equations       supply(i)	Ingredient supply
		objprimal	Defines Z;

supply(i)..	sum(j, a(i,j)*Y(j)) =L= s(i);

objprimal..     Z =e= sum(j, p(j)*Y(j));

MODEL   PRIMAL /supply, objprimal/;

solve PRIMAL using LP maximizing Z;

variable	W	Dual objective;

nonnegative
variables	PI(i)	Shadow price of ingredient i;

equations	profit(j)	Dual constraint,
		objdual		Dual objectives;

objdual..	W =e= sum(i, PI(i)*s(i));

profit(j)..	sum(i, PI(i)*a(i,j)) =g= p(j);

model DUAL /objdual, profit/;
solve dual using LP minimizing W;



