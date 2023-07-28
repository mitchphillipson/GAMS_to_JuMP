$title    GE Model in Mathiesens Framework -- One Activity and Three Markets


parameter    
        ld0    Labor demand /0.6/
        ls0    Leisure demand /0.4/
        le0    Labor endowment

        kd0    Capital demand /0.4/
        ke0    Capital endowment

        sigma    Elasticity of substitution in production /0.7/
        theta    Labor value share
        alpha    Consumption value share;

le0 = ld0 + ls0;
ke0 = kd0;

theta = ld0/(ld0+kd0);
alpha = 1/(1+ls0);

nonnegative
variables    
        PL    Wage rate
        PK    Return to capital
        PY    Output price
        Y     Production;

* This is not the C from slide 29 of the lecture notes. What is this C? C is unit cost function from slide 22.
* Why is PY not present in this value? Isn't that because PY is output price, not the factor prices?
$macro    C    ((theta*PL**(1-sigma)+(1-theta)*PK**(1-sigma))**(1/(1-sigma)))

* These are from slide 23 of the lecture notes. But why is AY 1? The constraint y = f(k,l) = 1
* Also, where is y/ \bar{y} in AL and AK? Seems like missing
$macro    AY    (1)
$macro    AL    (-ld0*(C/PL)**sigma)
$macro    AK    (-kd0*(C/PK)**sigma)

* Where do these D's come from? And why is DK 0? It seems like Market output & inputs and they are used in the equations below. But I don't know either why DK is 0.

$macro    DY    (alpha*(PL*le0+PK*ke0)/PY)
$macro    DL    ((1-alpha)*(PL*le0+PK*ke0)/PL)
$macro    DK    (0)

* BY is 0 because an "output endowment" doesn't make sense (is this correct?) Yes, only labor and capital endowment.
$macro    BY    (0)
$macro    BL    (le0)
$macro    BK    (ke0)

equations    profit, markety, marketk, marketL;

profit..    -PY*AY - PL*AL - PK*AK =g= 0;

markety..    BY + AY*Y =e= DY;

marketk..    BK + AK*Y =g= DK;

marketL..    BL + AL*Y =g= DL;

model    equilibrium /profit.Y, markety.PY, marketk.PK, marketL.PL/;

PL.L = 1;
PK.L = 1;
PY.L = 1;
Y.L = 1;

equilibrium.iterlim = 0;
solve equilibrium using mcp;

PL.LO = 1e-4;
PK.LO = 1e-4;
PY.FX = 1;

le0 = 1.1*le0;
equilibrium.iterlim = 10000;
solve equilibrium using mcp;
