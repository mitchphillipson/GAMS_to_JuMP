$title  surplus maximization and market equilibrium

set     t /a,b,c,d/;

table   tech    Technology
        cost    cap
a       2       2
b       5       2
c       7       4
d       10      inf;

parameter       c(t)    Cost by technology;
c(t) = tech(t,"cost");

nonnegative variables   P,PS,CS,s,Q(t);
free variable           obj;
equations               price, supply, psurplus, csurplus, objective;

price..         P =e= 10 - S*10/6;

supply..        S =e= sum(t, Q(t));

psurplus..      PS =e= sum(t, (P-c(t))* Q(t));

csurplus..      CS =e= (10 - P)*S/2;

objective..     OBJ =e= CS + PS;

Q.UP(t) = tech(t,"cap");

model equil /all/;
solve equil using nlp maximizing OBJ;
