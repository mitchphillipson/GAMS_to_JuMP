$title	Nash Equilibrium Maintenance

set	r	Regions			/ercot, ferc/,
	s	States of nature	/summer, winter, vortex/;

alias (r,rr), (s,ss);

$if not set w0    $set w0    0.25
$if not set lamda $set lamda 1.5
$if not set m0    $set m0    0.4
$if not set sigma $set sigma 0.25
$if not set esub  $set esub  0.25

parameter
	sigma		Degree of relative risk aversion /%sigma%/,
	esub		Elasticity of substitution between electricity and other goods /%esub%/,
	pi(s)		State probabilities (read as days per year) /summer 300, winter 64,  vortex 1/,
	c0		Reference aggregate expenditure (electricity plus other goods) /31/,
	lamda(r,s)	Electricity demand shock (used in the model),
	w(r,s)		Weather shock (used in the model);

*	Scale from days per year to probabilities:

pi(s) = pi(s) / 365;

table	w0(r,s)		Weather shocks (fraction loss)
				summer	winter	vortex
		ferc		0	0.2	0.25
		ercot		0	0	%w0%;

table	lamda0(r,s)	Electricity demand shock
				summer	winter	vortex
		ferc		1	1.2	1.2
		ercot		1	1	%lamda%

parameter	m0(r)	Benchmark level of maintenance /ferc 0.8, ercot %m0%/;

parameter	thetae	Energy value share; thetae = 1/c0;

parameter	alpha		Maintenance cost scale factor (calibrated),
		gamma		Maintenance cost exponent (calibrate);
alpha = 0;
gamma = 0;

nonnegative variables

	K(r)		Capacity
	M(r)		Maintenance

	PE(r,s)		Wholesale price
	PR(r)		Electricity generation resource

	X(r,s)		Sales to other region
	PX(s)		Traded electricity price;


equations	profit_G, profit_K, profit_M, profit_X, 
		market_PE, market_PR, market_PX;


variables
	PEU(r)		Shadow price of expected utility
	PC(r,s)		State-contingent price
	EU(r)		Expected utility
	C(r,s)		Consumption;

equations	PEU_def, PC_def, EU_def, C_def;

*	-----------------------------------------------------------------------------

PC_def(r,s)..		PC(r,s) =e= (thetae*(PE(r,s)*lamda(r,s))**(1-esub) + 1 - thetae)**(1/(1-esub));

PEU_def(r)..		PEU(r) =e= sum(ss, pi(ss) * PC(r,ss)**(1-sigma))**(1/(1-sigma));

*	At income level 1, region r agent can afford 1/PEU(r) units of expected utility:

EU_def(r)..		EU(r) =e= 1/PEU(r);

*	Demand for consumption in region r, state s is a function of the level of aggregate utility
*	and the real price of consumption in that state:

C_def(r,s)..		C(r,s) =e= EU(r) * (PEU(r)/PC(r,s))**sigma;
*	-----------------------------------------------------------------------------

profit_K(r)..		sqrt(PR(r)) =e= sum(s, pi(s)*PE(r,s)*(1 - w(r,s)*(1-M(r))));

profit_M(r)..		alpha /(1 - M(r))**gamma =e= sum(s, pi(s)*PE(r,S)*w(r,s));

profit_X(r,s)..		PE(r,s) =e= PX(s);

*	-----------------------------------------------------------------------------

market_PE(r,s)..	K(r)*(1-w(r,s)*(1-M(r))) =e= C(r,s) * lamda(r,s) * (PC(r,s)/(PE(r,s)*lamda(r,s)))**esub + X(r,s);

market_PR(r)..		0.5 =e= 0.5 * K(r) * 1/sqrt(PR(r));

market_PX(s)..		sum(r, X(r,s)) =e= 0;

*	-----------------------------------------------------------------------------

model network /	profit_K.K, profit_M.M, profit_X.X,
		market_PE.PE, market_PR.PR, market_PX.PX,
		PEU_def.PEU, PC_def.PC, EU_def.EU, C_def.C /;

PEU.L(r) = 1;
PC.L(r,s) = 1;
EU.L(r) = 1;
C.L(r,s) = 1;
PE.L(r,s) = 1;
K.L(r) = 1;
M.FX(r) = 0;
PE.L(r,s) = 1;
PR.L(r) = 1;

X.FX(r,s) = 0;
PX.FX(s) = 1;

*	First run ignores weather disturbances

lamda(r,s) = 1;
w(r,s) = 0;

network.iterlim = 0;
solve network using mcp;

*	Impose the weather shocks:

w(r,s) = w0(r,s);
lamda(r,s) = lamda0(r,s);

*	First solve with exogenous maintenance is used to impute
*	the shadow cost of maintenance:

M.FX(r) = m0(r);

network.iterlim = 10000;
solve network using mcp;

*	Calibrate parameters of the maitenance cost function which are consistent
*	with the assumed maintenance levels in ercot and ferc:

gamma = ( log(sum(s, pi(s)*PE.L("ercot",S)*W("ercot",s))) - log(sum(s, pi(s)*PE.L("ferc",S)*W("ferc",s)) ) ) /
	( log(1-M.L("ferc")) - log(1-M.L("ercot")) );

alpha = sum(s, pi(s)*PE.L("ercot",S)*W("ercot",s)) * (1-M.L("ercot"))**gamma;
display gamma, alpha;

*	Double check that this calibration is correct -- we should be at the solution of the model
*	with endogenous maintenance.  Unfix maintenance levels:

M.UP(r) = 1;
M.LO(r) = 0;

network.iterlim = 0;
solve network using mcp;
abort$round(network.objval,3) "BaU replication error.";
network.iterlim = 10000;

parameter	k0(r)	Reference capacity;
k0(r) = K.L(r);

set	items	Report items /
	EU		Expected utility (index)
	K		Installed capacity (reference = 1)
	CM		Maintenance cost 
	M		Maintenance share
	PR		Generation resource price
	PE		Electricity price
	C		End-use consumption
	PHI		Weather impact net maintenance savings (w*(1-M))
	PC		Consumption price
	X		Electricity net exports
	"surplus$"	Producer and consumer surplus (value)
	"surplus%"	Producer and consumer surplus (% relative to notrade)
	"M%"		Maintenace (%)
	"G"		Net generation
	"E"		Electricity demand (net generation less exports)
	"X%G"		Electricity exports as % of generation /;

set	scn	Scenarios /
		ShortRun	Impacts with fixed generating capacity
		LongRun		Impacts with price-responsive generating capacity,
		Climate		Impacts with climate change and price-responsive capacity /

	trd	Trade structure /
		notrade		ERCOT and FERC remain independent networks
		trade		ERCOT and FERC are interconnected /;

set	scen(scn,trd) /(ShortRun,Longrun,Climate).(trade,notrade)/;

parameter	report(scn,trd,items,*,*)	Summary report;

$onechoV >%gams.scrdir%report.gms
$ondotl
report(%scn%,"EU",r,".") = EU(r);
report(%scn%,"K",r,".")  = K.L(r);
report(%scn%,"CM",r,".") = K.L(r) * alpha/(1-gamma) * (1-(1-M.L(r))**(1-gamma));
report(%scn%,"M",r,".")  = M.L(r);
report(%scn%,"PR",r,".") = PR.L(r);
report(%scn%,"PE",".",s) = PX.L(s);
report(%scn%,"c",r,s)    = C(r,s)/c0;
report(%scn%,"PHI",r,s)  = w(r,s)*(1-M.L(r));
report(%scn%,"PE",r,s)   = PE.L(r,s);
report(%scn%,"PC",r,s)   = PC(r,s);
report(%scn%,"X",r,s)    = X.L(r,s);
$offecho

$set scn '"LongRun","notrade"'
$batinclude %gams.scrdir%report 

X.LO(r,s) = -inf;
X.UP(r,s) = +inf;
PX.LO(s) = 0; PX.UP(s) = +inf;

solve network using mcp;

$set scn '"LongRun","trade"'
$batinclude %gams.scrdir%report 

X.FX(r,s) = 0;
PX.FX(s) = 1;
K.FX(r) = k0(r);

SOLVE network using mcp;


$set scn '"ShortRun","notrade"'
$batinclude %gams.scrdir%report 

X.LO(r,s) = -inf;
X.UP(r,s) = +inf;
PX.LO(s) = 0; PX.UP(s) = +inf;

solve network using mcp;

$set scn '"ShortRun","trade"'
$batinclude %gams.scrdir%report 

*	Increase the number of days with winter vortex:

pi("summer") = 300;
pi("vortex") = 10;
pi("winter") = 65 - 10;
pi(s) = pi(s) / 365;

*	Climate scenarios are long-run:

K.UP(r) = inf; K.LO(r) = 0;

*	No-trade:

X.FX(r,s) = 0;
PX.FX(s) = 1;

solve network using mcp;

$set scn '"Climate","notrade"'
$batinclude %gams.scrdir%report 

*	Trade:

X.LO(r,s) = -inf;
X.UP(r,s) = +inf;
PX.LO(s) = 0; PX.UP(s) = +inf;

solve network using mcp;

$set scn '"Climate","trade"'
$batinclude %gams.scrdir%report 

*	Producer profit (resource revenue less maintenance cost):

parameter	profit(*,*,r)	Producer profit;
profit(scen,r) = 0.5*report(scen,"PR",r,".") - report(scen,"CM",r,".");

report(scen,"surplus$",r,"Consumer") = (report(scen,"EU",r,".") - report("LongRun","notrade","EU",r,"."))*c0;
report(scen,"surplus$",r,"Producer") = profit(scen,r);
report(scen,"surplus$",r,"Total") = (report(scen,"EU",r,".") - report("LongRun","notrade","EU",r,"."))*c0+profit(scen,r);

report(scen,"surplus%",r,"Consumer") = 100 * (report(scen,"EU",r,".")/report("LongRun","notrade","EU",r,".")-1);
report(scen,"surplus%",r,"Producer") = 100 * (profit(scen,r)/profit("Longrun","notrade",r)-1);
report(scen,"surplus%",r,"Total") = 00 * ((report(scen,"EU",r,".")*c0+profit(scen,r))/
			    (report("Longrun","notrade","EU",r,".")*c0+profit("Longrun","notrade",r))-1);

report(scen,"M%",r,".") = 100 * report(scen,"M",r,".");
report(scen,"G",r,s) = report(scen,"K",r,".")*(1-report(scen,"PHI",r,s));
report(scen,"E",r,s) = report(scen,"K",r,".")*(1-report(scen,"PHI",r,s))-report(scen,"X",r,s);
report(scen,"X%G",r,s) = 100*report(scen,"X",r,s)/(report(scen,"K",r,".")*(1-report(scen,"PHI",r,s)));

