$title  Optimal Thesis Planning

set     i       Activities /
        A       "Formulate an idea"
        B       "Survey policy literature"
        C       "Code up a small illustrative example"
        D       "Write a note on the policy issues and context"
        E       "Meet with advisor for idea approval"
        F       "Search academic literature"
        G       "Collect empirical data"
        H       "Formulate a more complete illustrative model"
        I       "Meet with advisor for model approval"
        J       "Formulate a calibrated pilot model"
        K       "Estimate key elasticities"
        L       "Meet with advisor for simulation approval"
        M       "Perform policy simulations"
        N       "Complete a thesis draft"
        O       "Obtain editorial feedback and advice"
        P       "Redraft the thesis"
        Q       "Meet with advisor for thesis acceptance"/;

alias (i,j);

set     prec(i,j)       Task precedence /
        A.(B*D), (B*D).E, E.(F*H), G.(K,I), H.I, I.J, J.L, 
        (K,L).M, (M,F).N, N.O, O.P, P.Q/;

parameter       delta(i)        Time to complete /
        A 15, B 4, C 5, D 10, E 4, F 15, G 15, H 5, 
        I 5 , J 10, K 4, L 3, M 4, N 20, O 5,P 2, Q 10/;

variables       S(i)    Start time of activity i,
                T       Total time to completion;

equations       start, finish;

start(prec(i,j))..      S(j) =g= S(i) + delta(i);

finish(i)..             T =g= S(i) + delta(i);

S.FX("A") = 0;

model cpm /all/;
solve cpm using LP minimizing T;

set	critpath(*,*)	Activities on the critical path;

loop(prec(i,j)$round(start.m(i,j)),
	critpath(i,"model1") = yes;
	critpath(j,"model1") = yes;
);

display critpath;

parameter       q(i)    Probability of successfully completing activity i,
                pi(i,j) Probability fail task i and return to task j 
                        /E.A 0.25, L.I 0.25, Q.P 0.15/;

q(i) = 1 - sum(j, pi(i,j));

variable        ED(i)   Expected duration of activity i;

equations       eddef, edstart, edfinish;

eddef(i)..      ED(i) =e= q(i)*delta(i) + sum(j, pi(i,j) * (S(i) - S(j) + ED(i)));

edstart(prec(i,j))..    S(j) =g= S(i) + ED(i);

edfinish(i)..           T =g= S(i) + ED(i);

model edcpm /eddef, edstart, edfinish/;

solve edcpm using LP minimizing T;

loop(prec(i,j)$round(edstart.m(i,j)),
	critpath(i,"model2") = yes;
	critpath(j,"model2") = yes;
);

display critpath;
