$title	Hydro Power Planning Model -- Scalar Format

FREE
VARIABLE
        REVENUE Aggregate earnings on sales in March and April;

NONNEGATIVE
VARIABLES
        PA_MAR  Power production -- plant A in March
        PA_APR  Power production -- plant A in April
        PB_MAR  Power production -- plant B in March
        PB_APR  Power production -- plant B in April
	      
        SA_MAR  Water spilled -- reservior A in March
        SA_APR  Water spilled -- reservior A in April
        SB_MAR  Water spilled -- reservior B in March
        SB_APR  Water spilled -- reservior B in April
	      
        LA_MAR  Level -- reservior A start of March (exogenous)
        LA_APR  Level -- reservior A start of April
        LA_MAY  Level -- reservior A start of May
        LB_MAR  Level -- reservior A start of March (exogenous)
        LB_APR  Level -- reservior A start of April
        LB_MAY  Level -- reservior A start of May
	      
        EH_MAR  High value electricity generation -- March
        EL_MAR  Low value electricity generation -- March
        EH_APR  High value electricity generation -- April
        EL_APR  Low value electricity generation -- April;

*	Specify Simple Bounds on Decision Variables

*	Levels at the beginning of March are specified exogenously:

LA_MAR.FX = 1900;
LB_MAR.FX =  850;

*	Power can be sold at CHF 5 per MWH up to 50,000 MWH per month: 

EH_MAR.UP = 50000;
EH_APR.UP = 50000;

*	Storage capacity:

LA_APR.UP = 2000;
LA_MAY.UP = 2000;
LB_APR.UP = 1500;
LB_MAY.UP = 1500;

*	Minimum reservior levels:

LA_APR.LO = 1200;
LA_MAY.LO = 1200;
LB_APR.LO = 800;
LB_MAY.LO = 800;

*	Power plant capacities:

PA_MAR.UP = 60000;
PA_APR.UP = 60000;
PB_MAR.UP = 35000;
PB_APR.UP = 35000;

equations	revenuedef, sales_mar, sales_apr, 
		alevel_apr, alevel_may, 
		blevel_apr, blevel_may;

*	Unlimited quantities of excess power can be sold for
*	CHF 3.50 per MWH.

*	Revenue therefore equals the sum of high and low price
*	revenue:

revenuedef..	REVENUE =E= 5*(EH_MAR+EH_APR) + 3.5*(EL_MAR+EL_APR);

*	Total sales in March and April equal generation:

sales_mar..	EH_MAR + EL_MAR =e= PA_MAR + PB_MAR;
sales_apr..	EH_APR + EL_APR =e= PA_APR + PB_APR;

*	Level in period t equals level in t-1 plus net inflows:

alevel_apr..	LA_APR =E= LA_MAR + 200 - SA_MAR - PA_MAR / 400;
alevel_may..	LA_MAY =E= LA_APR + 130 - SA_APR - PA_APR / 400;

blevel_apr..	LB_APR =E= LB_MAR + 40 + SA_MAR + PA_MAR/400 - SB_MAR - PB_MAR / 200;
blevel_may..	LB_MAY =E= LB_APR + 15 + SA_APR + PA_APR/400 - SB_APR - PB_APR / 200;

model hydro /all/;

solve hydro using lp maximizing revenue;
