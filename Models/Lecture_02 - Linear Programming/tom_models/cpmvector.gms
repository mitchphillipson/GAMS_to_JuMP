$title Critical Path Method -- Concert Planning: An Indexed Formulation

set activity	Activities involved in putting on a concert /
        A       "Find Site",
        B       "Find Engineers",
        C       "Hire Opening Act",
        D       "Set Radio and TV Ads",
        E       "Set Up Ticket Agents",
        F       "Prepare Electronics",
        G       "Print Advertising",
        H       "Set up Transportation",
        I       "Rehearsals",
        J       "Last-Minute Details"/;

*	Create two alternative symbols representing this set:

alias (activity,i,j);

parameter duration(activity) In days /
  A 3, B 2, C 6, D 2, E 3, F 3, G 5, H 1, I 1.5, J 2 /;

set prec(i,j) "Precedence order" /
  A.(B,C,E), B.F, C.(D,G,H), (F,H).I, I.J /;

FREE VARIABLE        T       Time to completion;

NONNEGATIVE VARIABLE S(i)    Starting time for activity i;

EQUATIONS   ctime(i) Completion time,  ptime(i,j) Sequence;

ctime(i)..         T =g= S(i) + duration(i);

ptime(prec(i,j)).. S(i) + duration(i) =L= S(j);

model cpm /all/;   solve cpm using lp minimizing T;
