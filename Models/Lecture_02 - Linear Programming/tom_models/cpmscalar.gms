$title  Critical Path Model -- Scalar Format

FREE
VARIABLE
        T       Completion time;

NONNEGATIVE
VARIABLES
        A       "Start time: find site",
        B       "Start time: find engineers",
        C       "Start time: hire opening act",
        D       "Start time: set radio and TV ads",
        E       "Start time: set up ticket agents",
        F       "Start time: prepare electronics",
        G       "Start time: print advertising",
        H       "Start time: set up transportation",
        I       "Start time: rehearsals",
        J       "Start time: last-minute details";

equations

*       Completion times no sooner than any of the tasks:

t_A, t_B, t_C, t_D, t_E, t_F, t_G, t_H, t_I, t_J 

*       Task which must be completed before subsequent task:

s_AB, s_AC, s_AE, s_BF, s_CD, s_CG, s_CH, s_FI, s_HI, s_IJ;

t_A.. T =G= A+3; t_B.. T =G= B+2; t_C.. T =G= C+6;
t_D.. T =G= D+2; t_E.. T =G= E+3; t_F.. T =G= F+3;
t_G.. T =G= G+5; t_H.. T =G= H+1; t_I.. T =G= I+1.5;
t_J.. T =G= J+2;

s_AB.. A+3 =L= B; s_AC.. A+3 =L= C; s_AE.. A+3 =L= E;
s_BF.. B+2 =L= F; s_CD.. C+6 =L= D; s_CG.. C+6 =L= G;
s_CH.. C+6 =L= H; s_FI.. F+3 =L= I; s_HI.. H+1 =L= I;
s_IJ.. I+1.5 =L= J;

model cpm /all/;
solve cpm using lp minimizing T;
