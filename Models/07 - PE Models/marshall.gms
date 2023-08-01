$title  surplus maximization and market equilibrium

set     t /a,b,c,d/;

table   tech    Technology
        cost    cap
a       2       2
b       5       2
c       7       4
d       10      inf;

parameter       c(t)    Cost by technology,
		k(t)	Capacity;

c(t) = tech(t,"cost");
k(t) = tech(t,"cap");


nonnegative variables   P,PS,CS,s,Q(t);
free variable           obj;
equations               price, supply, psurplus, csurplus, objective;

price..         P =e= 10 - S*10/6;

supply..        S =e= sum(t, Q(t));

psurplus..      PS =e= sum(t, (P-c(t))* Q(t));

csurplus..      CS =e= (10 - P)*S/2;

objective..     OBJ =e= CS + PS;

Q.UP(t) = k(t);

model equil /all/;
solve equil using nlp maximizing OBJ;

nonnegative 
variables	Q(t)	Quantity produced,
		RK(t)	Rental rate,
		D	Market demand,
		P	Market price;

equations	capacity(t), profit(t), demand, market;

*	Capital rent plus unit cost of operation must
*	be no less than the market price:

profit(t)..	RK(t) + c(t) =G= P;

capacity(t)..	k(t) =g= Q(t);

demand..	D =e= (6 * (1-P/10));

market..	sum(t, Q(t)) =g= D;

model marketeq /profit.Q, capacity.RK, demand.D, market.P/;

solve marketeq using mcp;
