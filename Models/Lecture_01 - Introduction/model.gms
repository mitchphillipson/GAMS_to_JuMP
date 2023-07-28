$title	Vanderbeis Local Warming Model

*       The readdata.gms program provides input data in a GAMS data
*       exchange format input file:

$gdxin data.gdx

*       The input data file contains labels of the days for which average
*       temperature data are available:

set     d(*)    Dates;
$load d

*       It also provides numeric indices for those dates and the associated
*       average temperatures:

parameter       day(d)  Day (integer) associated with dates
                ave(d)  Average temperatures;

$load day ave

*       Declare the optimization model: L1 norm minimization.

NONNEGATIVE
VARIABLE        DEV(d)  Deviation;

VARIABLES       T(d)            Estmated temperature,

                X0              Constant coefficient,
                X1              Linear  (warming) coefficient,
                X2              Cosine coefficient -- annual cycle,
                X3              Sine coefficient -- annual cycle,
                X4              Cosine coefficient -- solar cycle ,
                X5              Sine coefficient -- solar cycle,

                OBJ             Objective;


equations       tdef, devlb, devub, objdef;

tdef(d)..       T(d) =e= X0 + X1*day(d) + X2*cos(2*pi*day(d)/365.25)
                                        + X3*sin(2*pi*day(d)/365.25)
                                        + X4*cos(2*pi*day(d)/(10.7*365.25))
                                         + X5*sin(2*pi*day(d)/(10.7*365.25));

devlb(d)..      DEV(d) =g= ave(d) - T(d);

devub(d)..      DEV(d) =g= T(d) - ave(d);

objdef..        OBJ =E= sum(d, DEV(d));

model calib /all/;

option limrow=0, limcol=0, solprint=off;

SOLVE CALIB USING LP MINIMIZING OBJ;

display X0.L, X1.L, X2.L, X3.L, X4.L, X5.L;

parameter       warming         Rate of warming (degrees F per century);

