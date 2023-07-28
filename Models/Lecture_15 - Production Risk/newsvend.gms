$title	Soren Nielsen's Newsvendor Model

set	i	Papers to be ordered /1*100/;   

alias( i, j);

parameter p(i) "Probility that i papers will be demanded";

parameter lambda	Parameter of the truncated Poisson distribution / 20.5 /;

p("1") = exp(-lambda) * lambda;

*	Make an assignment to the next value.  When i reaches 100, the
*	reference off the end of the set is ignored:

loop(i,	p(i+1) = p(i) * lambda / i.val;);

*	Normalize so the p(i)'s sum to 1:

p(i) = p(i) / sum(j, p(j));
display p;

*	Done calculating p(i).

*	Problem data are here.
parameters
    c_s   "Price at which papers are sold"             / 35 /,
    c_p   "Price at which papers are purchased"        / 20 /,
    c_f   "Price for fish-wrapping (unsold) papers"    /  2 /,
    N     "Max number of papers that can be purchased" / 40 /;

parameter f(i,j) "Profit if i papers are purchased and j are demanded";

F(i,j) = c_s * min(i.val, j.val) + c_f * max(0, i.val-j.val) - c_p * i.val;
display F; 

*	Brute force calculation of the optimal value:

parameter g(i) "Expected profit if i papers are purchased";
g(i)$(i.val<=N) = sum(j,p(j) * F(i,j));
display g;

*	Find the highest possible expected profit, and the associated
*	number:
set		i_hat(i)	Papers orders for which profit is maximized;
parameter	max_prof	The maximum profit;

max_prof = smax(i, g(i)); i_hat(i) = yes$(g(i)=max_prof);
display max_prof, i_hat;

*	Set up the newsvendor problem as a two-stage, stochastic program.                      *

*	We now formulate the problem as a two-stage, linear stochastic program (SP).

*	The first-stage decision is: How many papers to purchase, represented by x.
*	The second-stage decision is: How many papers to sell, represented by z(s).
*	The index s is the "scenario", that is, how many papers are demanded.
*
*	Notice the two-stage decision process (two decision points):
*	1. Purchases in the morning, then wait and see 
*	how many are demanded/sold,
*	2. Decide how many papers to sell as fish wrapping, 
*	in this case a trivial decision.

*	We shift from i to s to capture the idea of scenarios:

alias  (i,s);

nonnegative
variables	X      "How many papers to purchase",
		Z(s)   "How many papers are sold";

*	First stage constraint is a simple upper bound -
*	the maximum number we can purchase:

X.UP = n;

*	Constraints on the second-stage variables relates

EQUATIONS purchased(s);

*	Can't sell more than we have purchased

purchased(s)..	Z(s) =L= X;

*	In scenario s, the number sold cannot exceed 
*	the number demanded.  This constraint can be
*	represented by a simple upper bound on Z(s)

Z.UP(s) = s.val;

*	Formulate the objective: expected profit.  Note that
*	we declare profit(s) as a free variable so that we may
*	have scenarios in which profits are negative.  

*	The NLP objective (EXP_PROFIT) must be a free variable.

variables
	PROFIT(s)     "Profit under each scenario",
	EXP_PROFIT    "Expected profit of the whole day";

equations profit_def(s), exp_profit_def;

profit_def(s)..

	PROFIT(s) =e= c_s * Z(s) + c_f * (X-Z(s)) - c_p * X;

exp_profit_def..

	EXP_PROFIT =E= sum(s, p(s) * PROFIT(s) );

model two_stage /all/;

solve two_stage maximizing EXP_PROFIT using lp;


