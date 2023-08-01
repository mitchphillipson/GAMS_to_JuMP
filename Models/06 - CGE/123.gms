$title  Static 123 Model Ala Devarjan

$include 123data.gms

parameter thetal  Labor share in cost function,
          thetah  Consumption share in expenditure function,
          thetam  Share parameter in Armington function
          thetaz  Share parameter in transformation function  ;

thetal = ld0*pl0 /(kd0*rk0+ld0*pl0);
thetaz = x0 *px0 /(d0+x0*px0);
thetam = m0 *pm0 /(d0+m0*pm0);
thetah = c0/(c0+l0);

Nonnegative Variables
*$SECTORS:
        Y        Production
        A        Armington composite
        M        Imports
        X        Exports

*$COMMODITIES:
        PD       Domestic price index
        PX       Export price index
        PM       Import price index
        PA       Armington price index
        PL       Wage rate index
        RK       Rental price index
        PFX      Foreign exchange

*$CONSUMERS:
        HH       Private households
        GOVT     Government

*$AUXILIARY:
        TAU      Replacement tax
;


$macro PKL   [({PL*(1+tl)/pl0}**thetal * {RK*(1+tk)/rk0}**(1-thetal))$(esubkl eq 1) + ({thetal * {PL*(1+tl)/pl0}**(1-esubkl) + (1-thetal) * {RK*(1+tk)/rk0}**(1-esubkl)}**{1/(1-esubkl)})$(esubkl ne 1)]
$macro LD   (ld0     * (PKL*pl0/{PL*(1+tl)})**esubkl)
$macro KD   (kd0     * (PKL*rk0/{RK*(1+tk)})**esubkl)

$macro PY   ({thetaz * {PX*(1-tx)/px0}**(1+etadx) + (1-thetaz) * PD**(1+etadx)} **{1/(1+etadx)})
$macro DY   (d0      * (PD/PY)**etadx)
$macro XY   (x0      * (PX*(1-tx)/{px0*PY})**etadx)

$macro PDM  ((thetam *(PM*(1+tm)/pm0)**(1-sigmadm) + (1-thetam)*PD**(1-sigmadm))**(1/(1-sigmadm)))
$macro DA   (d0      *(PDM/PD)**sigmadm)
$macro MA   (m0      *(PDM*pm0/{PM*(1+tm)})**sigmadm)

$macro PH   ({thetah *PA**(1-sigma) + (1-thetah)*PL**(1-sigma)}**{1/(1-sigma)})
$macro C    (c0      *(PH/PA)**sigma *1/PH *HH/{c0+l0})
$macro L    (l0      *(PH/PL)**sigma *1/PH *HH/{c0+l0})

Equations

* Zero profit condition
         profity         domestic production,
         profita         Armington supply,
         profitm         imported goods production
         profitx         exported goods production

* Market clearing condition
         marketd         domestic  goods market,
         marketa         Armington goods market
         marketm         imported  goods market
         marketx         exported  goods market
         marketfx        balance of payment
         marketk         capital market
         marketl         labor   market

* Income balance
         incomeg         government budget
         incomeh         household budget

* Additional constraints
         taudef          equal yield constraint;

profity..  PKL*(ld0*pl0 + kd0*rk0) =g= PY*(d0+x0*px0);

profita..  PDM*(m0*pm0 + d0) =e= PA*a0*(1-ta);

profitm..  PFX*pwm =e= PM;

profitx..  PX =e= PFX*pwx;

marketd..  Y*DY =e= A*DA;

marketa..  A*a0 =g= GOVT/PA + C + i0;

marketm..  M*m0 =e= A*MA;

marketx..  Y*XY =e= X*x0;

marketfx.. X*pwx*x0 - M*pwm*m0 =e= -bopdef ;

marketk..  kd0    =e= Y*KD;

marketl..  ld0+l0 =e= Y*LD + L;

incomeg..  GOVT =e= PFX*bopdef + PA*dtax +  PA*g0*TAU
                  + tx*PX*XY*Y + tk*RK*KD*Y + tl*PL*LD*Y +  tm*PM*MA*A + ta*PA*a0*A;

taudef..   GOVT =e= PA * g0;

incomeh..  HH =e= PL*(ld0+l0) - PA*dtax - PA*g0*TAU + RK*kd0 - PA*i0  ;

model mcp123 /marketd.PD, marketa.PA, marketm.PM, marketx.PX, marketfx.PFX, marketk.RK, marketl.PL,
              profity.Y,  profita.A,  profitm.M,  profitx.X,  incomeg.GOVT, incomeh.HH, taudef.TAU/;


TAU.LO = -INF;
TAU.L  = 0;

Y.L = 1;
A.L = 1;
M.L = 1;
X.L = 1;

PD.L = 1;
PX.L = 1;
PM.L = 1;
PA.L = 1;
PL.L = 1;
RK.L = 1;
PFX.L= 1;

GOVT.L = g0;
HH.FX  = c0+l0;

mcp123.iterlim = 0;
solve mcp123 using mcp;

tm = 0;
mcp123.iterlim = 10000;
solve mcp123 using mcp;
