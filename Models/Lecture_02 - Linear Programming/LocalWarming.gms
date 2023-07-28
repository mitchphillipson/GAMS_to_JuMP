$title	Vanderbeis Local Warming Model

*	The data source for this model is described here:

set	days /1*31000/

*		https://www.ncei.noaa.gov/pub/data/ghcn/daily/readme.txt


*	Choose a weather station to model (n.b. The dataset name is case-sentivity!)

*	USW00014837	Dane county regional airport (1939-present)

*	You can find dataset names here:

*		https://www.ncei.noaa.gov/pub/data/ghcn/daily/ghcnd-stations.txt	

$if not set ds $set ds USW00014837

*	Use curl to download the dataset.

*	All versions of Microsoft Windows 10 and Windows 11 get curl installed
*	by default. The initial curl version Microsoft shipped was 7.55. 
*	It was upgraded to 7.79. 1 in January 2022.

$if not exist %ds%.dly $call curl --remote-name https://www.ncei.noaa.gov/pub/data/ghcn/daily/all/%ds%.dly

*	Execute a Fortran program to parse the data:

$call translate %ds% "%gams.scrdir% "

*       The readdata.gms program provides input data in a GAMS data
*       exchange format input file:

set	yr(*)	Year
	m	Month /01*12/
	dm	Day of the month /01*31/;

parameter	tavg(yr<,m,dm)	Average temperature (tenths of degrees C) /
$offlisting

*	Note that TAVG corresponds to an average for the period ending at 2400 UTC 
*	rather than local midnight

$include %gams.scrdir%%ds%.gms
$onlisting
/
;




set	d(yr,m,dm)	Set of days for which we have data;
option d<tavg;

parameter	tave	Overall average temperature;
tave = sum(d,tavg(d)/10)/card(d);

*       It also provides numeric indices for those dates and the associated
*       average temperatures:

parameter       day(yr,m,dm)  "Day (integer) associated with dates (Jan 1, 1900 = 1)";
day(d(yr,m,dm)) = jdate(yr.val,m.val,dm.val);

parameter	dmin	Minimum day index
		mindays	Required elements of set days;

dmin = smin(d,day(d));
mindays = smax(d,day(d)) - dmin + 1;
abort$( mindays > card(days)) "Need to increase size of set days:",mindays;

*	Convert to a day index which begins from 1 on the first day for which we have data:

day(d) = day(d) - dmin + 1;


NONNEGATIVE
VARIABLE        DEV(yr,m,dm)  Deviation;

VARIABLES       T(yr,m,dm)	Estmated average temperature,

                X0              Constant coefficient,
                X1              Linear  (warming) coefficient,
                X2              Cosine coefficient -- annual cycle,
                X3              Sine coefficient -- annual cycle,
                X4              Cosine coefficient -- solar cycle,
                X5              Sine coefficient -- solar cycle,

                OBJ             Objective;


equations       tdef, devlb, devub, objdef;

tdef(d)..       T(d) =e= X0 + X1*day(d) + X2*cos(2*pi*day(d)/365.25) + X3*sin(2*pi*day(d)/365.25)
                                        + X4*cos(2*pi*day(d)/(10.7*365.25)) + X5*sin(2*pi*day(d)/(10.7*365.25));

devlb(d)..      DEV(d) =g= tavg(d)/10 - T(d);

devub(d)..      DEV(d) =g= T(d) - tavg(d)/10;

objdef..        OBJ =E= sum(d, DEV(d));

model calib /all/;

*	Turn off the listing file to reduce run time:

option limrow=0, limcol=0, solprint=off;

SOLVE CALIB USING LP MINIMIZING OBJ;






*	V = a*sin(theta) + b*cos(theta)
*	m = sqrt(a**2 + b**2)
*	alpha = (a/m)**2
*	V = m (sqrt(alpha) sin(theta) + sqrt(1-alpha) cos(theta))
*	Given m and alpha we have:
*	a = sqrt(alpha) * m	and	b = sqrt(1-alpha) * m

set	mdl	Models we are computing /L1, L2/;

parameter	solutions(*,mdl)	Solution values;

$set mdl L1

solutions("warming","%mdl%") = X1.L*365.25*100;

solutions("X0","%mdl%") = X0.L;
solutions("X1","%mdl%") = X1.L;
solutions("X2","%mdl%") = X2.L;
solutions("M","%mdl%") = sqrt(sqr(X2.L) + sqr(X3.L));
solutions("alpha","%mdl%") = sqr(X2.L/solutions("M","%mdl%"));
solutions("X3","%mdl%") = X3.L;

solutions("X4","%mdl%") = X4.L;
solutions("X5","%mdl%") = X5.L;
solutions("MS","%mdl%") = sqrt(sqr(X4.L) + sqr(X5.L));
solutions("alphas","%mdl%") = sqr(X4.L/solutions("MS","%mdl%"));

parameter	pivotdata(days,yr,m,dm,*)	Pivot report of data;

*	Store data for the day and average temperature data:

loop((d(yr,m,dm),days)$(days.val=day(d)),
	pivotdata(days,d,"Actual") = tavg(d)/10;
	pivotdata(days,d,"T(%mdl%)") = X0.L + X1.L*day(d) 
				+ X2.L*cos(2*pi*day(d)/365.25)
				+ X3.L*sin(2*pi*day(d)/365.25)
				+ X4.L*cos(2*pi*day(d)/(10.7*365.25))
				+ X5.L*sin(2*pi*day(d)/(10.7*365.25));
	pivotdata(days,d,"GW(%mdl%)") = X1.L*(day(d)-1) + tave;
	pivotdata(days,d,"Annual(%mdl%)") = 
				+ X2.L*cos(2*pi*day(d)/365.25)
				+ X3.L*sin(2*pi*day(d)/365.25);
	pivotdata(days,d,"Solar(%mdl%)") = 
				+ X4.L*cos(2*pi*day(d)/(10.7*365.25))
				+ X5.L*sin(2*pi*day(d)/(10.7*365.25));
);


*       Declare the optimization model: L2 norm minimization.

equations	objlsq	Least squares objective;

objlsq..        OBJ =E= sum(d, sqr(
			X0 + X1*day(d) + X2*cos(2*pi*day(d)/365.25)
                                        + X3*sin(2*pi*day(d)/365.25)
                                        + X4*cos(2*pi*day(d)/(10.7*365.25))
                                        + X5*sin(2*pi*day(d)/(10.7*365.25)) - tavg(d)/10));
model calibL2 /objlsq/;

solve calibL2 using qcp minimizing OBJ;

$set mdl L2

solutions("warming","%mdl%") = X1.L*365.25*100;
solutions("X0","%mdl%") = X0.L;
solutions("X1","%mdl%") = X1.L;
solutions("X2","%mdl%") = X2.L;
solutions("M","%mdl%") = sqrt(sqr(X2.L) + sqr(X3.L));
solutions("alpha","%mdl%") = sqr(X2.L/solutions("M","%mdl%"));
solutions("X3","%mdl%") = X3.L;
solutions("X4","%mdl%") = X4.L;
solutions("X5","%mdl%") = X5.L;
solutions("MS","%mdl%") = sqrt(sqr(X4.L) + sqr(X5.L));
solutions("alphas","%mdl%") = sqr(X4.L/solutions("MS","%mdl%"));
display solutions;

loop((d(yr,m,dm),days)$(days.val=day(d)),
	pivotdata(days,d,"T(%mdl%)") = X0.L + X1.L*day(d) 
				+ X2.L*cos(2*pi*day(d)/365.25)
				+ X3.L*sin(2*pi*day(d)/365.25)
				+ X4.L*cos(2*pi*day(d)/(10.7*365.25))
				+ X5.L*sin(2*pi*day(d)/(10.7*365.25));
	pivotdata(days,d,"GW(%mdl%)") = X1.L*(day(d)-1) + tave;
	pivotdata(days,d,"Annual(%mdl%)") = 
				+ X2.L*cos(2*pi*day(d)/365.25)
				+ X3.L*sin(2*pi*day(d)/365.25);
	pivotdata(days,d,"Solar(%mdl%)") = 
				+ X4.L*cos(2*pi*day(d)/(10.7*365.25))
				+ X5.L*sin(2*pi*day(d)/(10.7*365.25));
);


execute_unload '%ds%.gdx', solutions, pivotdata;
execute 'gdxxrw i=%ds%.gdx o=%ds%.xlsx par=solutions rng=Solutions!a1..c12 par=pivotdata rng=PivotData!a2 cdim=0 intastext=y';
