$TITLE  Brewery Profit Maximization

variables       X       Production of lager,
                Y       Production of ale,
                Z       Profit (maximand);

nonnegative variables   X, Y;

equations       malt    Malt budget,
                yeast   Yeast budget,
                profit  Defines Z;

malt..          4 * X + 2 * Y =L= 4800;

yeast..         X + Y =L= 1750;

profit..        Z =e= 120 * X + 90 * Y;

*       Include hops constraints as upper bounds:

X.UP = 1000; Y.UP = 1750;

MODEL   PRIMAL /malt, yeast, profit/;

solve PRIMAL using LP maximizing Z;


VARIABLE	W	Dual ojective;

NONNEGATIVE
VARIABLES	LAMDA1, LAMDA2, LAMDA3, LAMDA4;

EQUATIONS	objdual, dualX, dualY;

objdual..	W =E= 4800*LAMDA1 + 1750*LAMDA2 + 1000*LAMDA3 + 1500*LAMDA4;

dualX..	4*LAMDA1 + 1*LAMDA2 + 1*LAMDA3 + 0*LAMDA4 =G= 120;

dualY..	2*LAMDA1 + 1*LAMDA2 + 0*LAMDA3 + 1*LAMDA4 =g= 90;

model DUAL /objdual, dualX, dualY/;

solve DUAL using LP minimizing W;
